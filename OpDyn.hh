#ifndef UTOPIA_MODELS_OPDYN_HH
#define UTOPIA_MODELS_OPDYN_HH

#include <cmath>
#include <functional>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/random.hpp>

#include <utopia/core/apply.hh>
#include <utopia/core/graph.hh>
#include <utopia/core/model.hh>
#include <utopia/core/types.hh>
#include <utopia/data_io/graph_utils.hh>

#include "ageing.hh"
#include "modes.hh"
#include "revision.hh"
#include "utils.hh"


namespace Utopia::Models::OpDyn {

/*! The different modes of the model: */
using modes::Mode;
using modes::Mode::None;
using modes::Mode::Ageing;
using modes::Mode::Media;
using modes::Mode::Ageing_and_Media;

/*!Each user-network node accomodates one user. Each user holds an opinion, is susceptible to others'
 opinions, and has a certain tolerance towards other opinions, which is the
 radius of interaction.*/
struct User {
    double opinion;
    double tolerance;
    double susceptibility;
    unsigned int age;
    size_t used_media;
};

/*!Each media-network node accomodates one medium; the principle is similar
to that of the users. Each medium holds an opinion, is able to convince users
of its stance (persuasiveness), is editorially flexible (i.e. able to change
its stance) (susceptibility), and has a certain tolerance for other media's
views. As in the user case, all parameters are dynamic. */

struct Medium {
    double opinion;
    double tolerance;
    double susceptibility;
    double persuasiveness;
    double ads;
    double ads_normalized;
    unsigned int users;
    unsigned int users_previous;
};

/// Each network edge has a certain weight, which can be negative or positive.
struct Weight {
    double attr;
};

/// The directed network type for the OpDyn Model:
using Network_u =   boost::adjacency_list<
                    boost::setS,        // edges
                    boost::vecS,        // vertices
                    boost::bidirectionalS,
                    User,               // vertex property
                    Weight>;            // edge property

/// The undirected network type for the OpDyn Model:
using Network_m =   boost::adjacency_list<
                    boost::setS,        // edges
                    boost::vecS,        // vertices
                    boost::undirectedS,
                    Medium,             // vertex property
                    Weight>;            // edge property

using OpDynTypes = ModelTypes<>;
using pair_double = std::pair<double,double>;
using pair_int = std::pair<int, int>;



/// The OpDyn Model

template<Mode model_mode=None>
class OpDyn:
    public Model<OpDyn<model_mode>, OpDynTypes>
{
public:
    /// The base model type
    using Base = Model<OpDyn<model_mode>, OpDynTypes>;

    /// Data type that holds the configuration
    using Config = typename Base::Config;

    /// Data type of the group to write model data to, holding datasets
    using DataGroup = typename Base::DataGroup;

    /// Data type for a dataset
    using DataSet = typename Base::DataSet;

    /// Data type of the shared RNG
    using RNG = typename Base::RNG;

private:
    // Base members: _time, _name, _cfg, _hdfgrp, _rng, _monitor

    std::uniform_real_distribution<double> _uniform_distr_prob_val;

    // User properties
    const Config _cfg_u;
    Network_u _nw_u;
    const double _radicalisation_parameter;
    const double _rewiring;
    unsigned int _rewiring_count;
    const double _weighting;

    const unsigned int _life_cycle;
    const double _replacement_rate;
    const pair_int _child_ages;
    const pair_int _parent_ages;
    const pair_int _senior_ages;

    // Media properties
    const Config _cfg_m;
    Network_m _nw_m;
    const int _media_time_constant;
    const unsigned int _num_media;
    pair_double _ads;
    pair_double _attr;

    // datasets and groups

    std::shared_ptr<DataGroup> _grp_nw_u;

    std::shared_ptr<DataGroup> _grp_nw_m;

    std::shared_ptr<DataSet> _dset_vertices_u;

    std::shared_ptr<DataGroup> _grp_edges_u;

    std::shared_ptr<DataSet> _dset_edges_u_initial;

    std::shared_ptr<DataSet> _dset_edges_u_final;

    std::shared_ptr<DataSet> _dset_opinion_u;

    std::shared_ptr<DataSet> _dset_tolerance_u;

    std::shared_ptr<DataSet> _dset_susceptibility_u;

    std::shared_ptr<DataSet> _dset_age_u;

    std::shared_ptr<DataSet> _dset_avg_nb_opinion_u;

    std::shared_ptr<DataSet> _dset_opinion_m;

    std::shared_ptr<DataSet> _dset_users;

    std::shared_ptr<DataSet> _dset_ads;

    std::shared_ptr<DataSet> _dset_rewiring_count;

    std::shared_ptr<DataSet> _dset_out_degree;

    std::shared_ptr<DataSet> _dset_in_degree;

    std::shared_ptr<DataSet> _dset_final_opinion_clusters;

    std::shared_ptr<DataSet> _dset_num_opinion_clusters;

    std::shared_ptr<DataSet> _dset_num_weighted_opinion_clusters;

    std::shared_ptr<DataSet> _dset_rel_bc;

    std::shared_ptr<DataSet> _dset_weights;


public:
    // Constructs the OpDyn model

    template<class ParentModel>
    OpDyn (const std::string name,
                    ParentModel &parent)
    :
        // Initialize first via base model
        Base(name, parent),
        _cfg_u(this->_cfg["nw_u"]),
        _cfg_m(this->_cfg["nw_m"]),

        // initialize network
        _nw_u(this->init_nw_u()),
        _nw_m(this->init_nw_m()),

        _uniform_distr_prob_val(std::uniform_real_distribution<double>(0., 1.)),

        // model parameters
        _life_cycle(get_as<int>("life_cycle", this->_cfg)),
        _replacement_rate(get_as<double>("replacement_rate", this->_cfg)),
        _child_ages(get_as<pair_int>("children", this->_cfg["age_groups"])),
        _parent_ages(get_as<pair_int>("parents", this->_cfg["age_groups"])),
        _senior_ages(get_as<pair_int>("seniors", this->_cfg["age_groups"])),
        _radicalisation_parameter(
                    get_as<double>("radicalisation_parameter", this->_cfg)),
        _rewiring(get_as<double>("rewiring", this->_cfg)),
        _rewiring_count(0),
        _weighting(get_as<double>("weighting", this->_cfg)),

        _media_time_constant(get_as<int>("media_time_constant", this->_cfg)),
        _num_media(num_vertices(_nw_m)),
        _ads(get_as<pair_double>("init_ads", this->_cfg)),
        _attr(get_as<pair_double>("attr", this->_cfg)),

        // create datagroups and datasets
        _grp_nw_u(Utopia::DataIO::create_graph_group(_nw_u, this->_hdfgrp,
                                                    "nw_users")),
        _grp_nw_m(Utopia::DataIO::create_graph_group(_nw_m, this->_hdfgrp,
                                                    "nw_media")),

        _dset_vertices_u(this->create_dset("_vertices", _grp_nw_u,
                        {boost::num_vertices(_nw_u)}, 5)),
        _grp_edges_u(_grp_nw_u->open_group("_edges")),
        _dset_edges_u_initial(_grp_edges_u->open_dataset("0",
                        {2, boost::num_edges(_nw_u)})),
        _dset_edges_u_final(_grp_edges_u->open_dataset("1",
                        {2, boost::num_edges(_nw_u)})),
        _dset_opinion_u(this->create_dset("opinion_u", _grp_nw_u,
                        {boost::num_vertices(_nw_u)}, 5)),
        _dset_tolerance_u(this->create_dset("tolerance_u", _grp_nw_u,
                                          {boost::num_vertices(_nw_u)}, 5)),
        _dset_susceptibility_u(this->create_dset("susceptibility_u", _grp_nw_u,
                                        {boost::num_vertices(_nw_u)}, 5)),
        _dset_age_u(this->create_dset("age_u", _grp_nw_u,
                                          {boost::num_vertices(_nw_u)}, 5)),
        _dset_opinion_m(this->create_dset("opinion_m", _grp_nw_m,
                        {boost::num_vertices(_nw_m)}, 5)),
        _dset_avg_nb_opinion_u(this->create_dset("avg_nb_opinion_u", _grp_nw_u,
                        {boost::num_vertices(_nw_u)}, 5)),
        _dset_users(this->create_dset("user_count", _grp_nw_m,
                        {boost::num_vertices(_nw_m)}, 5)),
        _dset_ads(this->create_dset("ads", _grp_nw_m,
                        {boost::num_vertices(_nw_m)}, 5)),
        _dset_rewiring_count(this->create_dset("rewiring_count", _grp_nw_u,
                        {}, 5)),
        _dset_out_degree(_grp_nw_u->open_dataset("out_degree",
                        {boost::num_vertices(_nw_u)})),
        _dset_in_degree(this->create_dset("in_degree", _grp_nw_u,
                        {boost::num_vertices(_nw_u)}, 5)),
        _dset_num_opinion_clusters(this->create_dset("num_opinion_clusters",
                        _grp_nw_u, {}, 5)),
        _dset_num_weighted_opinion_clusters(this->create_dset(
                        "num_weighted_opinion_clusters", _grp_nw_u, {}, 5)),
        _dset_rel_bc(this->create_dset("rel_bc", _grp_nw_u,
                        {boost::num_vertices(_nw_u)}, 5)),
        _dset_weights(this->create_dset("weights", _grp_nw_u,
                        {boost::num_edges(_nw_u)}, 5))
    {
        this->_log->debug("Constructing the OpDyn Model ...");

        this->initialize_properties();

        this->_log->info("Initialized user network with {} vertices and {} edges",
                         num_vertices(_nw_u), num_edges(_nw_u));
        this->_log->info("Initialized media network with {} vertices and {} edges",
                         num_vertices(_nw_m), num_edges(_nw_m));

        // Write the vertex data once as it does not change
        auto [v, v_end] = boost::vertices(_nw_u);
        auto [e, e_end] = boost::edges(_nw_u);

        _dset_vertices_u->write(v, v_end, [&](auto vd){
                       return boost::get(boost::vertex_index_t(), _nw_u, vd);
        });

        _dset_out_degree->write(v, v_end, [&](auto vd){
                                return out_degree(vd, _nw_u);
        });

        _dset_edges_u_initial->write(e, e_end, [&](auto ed){
                return boost::get(boost::vertex_index_t(), _nw_u,
                                  boost::source(ed, _nw_u));
        });

        _dset_edges_u_initial->write(e, e_end, [&](auto ed){
                return boost::get(boost::vertex_index_t(), _nw_u,
                                  boost::target(ed, _nw_u));
        });

        Utopia::DataIO::save_graph(_nw_m, _grp_nw_m);

        _dset_opinion_u->add_attribute("is_vertex_property", true);
        _dset_opinion_m->add_attribute("is_vertex_property", true);
        _dset_users->add_attribute("is_vertex_property", true);
        _dset_ads->add_attribute("is_vertex_property", true);
        _dset_weights->add_attribute("is_edge_property", true);
        _dset_opinion_u->add_attribute("dim_name__1", "vertex");
        _dset_opinion_u->add_attribute("coords_mode__vertex", "start_and_step");
        _dset_opinion_u->add_attribute("coords__vertex", std::vector<std::size_t>{0, 1});
    }

private:

    // Setup functions .........................................................

    void initialize_properties() {
        this->_log->debug("Initializing network properties ...");

        /// Initialize the media network properties if the media network is turned on; this is done first
        /// to be able to set the media user count to 0.
        if constexpr (model_mode == Media or Ageing_and_Media) {
            for (auto v : range<IterateOver::vertices>(_nw_m)) {

                _nw_m[v].opinion=utils::initialize(
                                   this->_cfg["opinion"]["media"],
                                   *this->_rng);
                _nw_m[v].tolerance=utils::initialize(
                                   this->_cfg["tolerance"]["media"],
                                   *this->_rng);
                _nw_m[v].susceptibility=utils::initialize(
                                   this->_cfg["susceptibility"]["media"],
                                   *this->_rng);
                _nw_m[v].persuasiveness=utils::initialize(
                                   this->_cfg["persuasiveness"]["media"],
                                   *this->_rng);
                _nw_m[v].users = 0.;
                _nw_m[v].ads = 0.;

                // set inter-media attractions
                for (auto e : range<IterateOver::out_edges>(v, _nw_m)) {
                    _nw_m[e].attr=utils::set_init_uniform(_attr, *this->_rng);
                }
            }
        }

        /// Second, initialize the user network properties:
        for (auto v : range<IterateOver::vertices>(_nw_u)) {

            _nw_u[v].age = utils::get_rand_int<RNG>(1, 85, *this->_rng);

            _nw_u[v].opinion=utils::initialize(
                                this->_cfg["opinion"]["users"],
                                *this->_rng);

            _nw_u[v].tolerance=utils::initialize(
                                _nw_u[v].age,
                                this->_cfg["tolerance"]["users"],
                                *this->_rng);

            _nw_u[v].susceptibility=utils::initialize(
                                _nw_u[v].age,
                                this->_cfg["susceptibility"]["users"],
                                *this->_rng);

            if constexpr (model_mode == Media or Ageing_and_Media) {
                // choose random medium
                _nw_u[v].used_media = utils::get_rand_int<RNG>(0,
                                                              _num_media,
                                                              *this->_rng);
                _nw_m[_nw_u[v].used_media].users++;
                _nw_m[_nw_u[v].used_media].ads++;
            }
            // set initial edge weight to 1/out-degree
            for (auto e : range<IterateOver::out_edges>(v, _nw_u)) {
                _nw_u[e].attr = 1. / double(out_degree(v, _nw_u));
            }
        }

        if constexpr (model_mode == Media or Ageing_and_Media) {
            revision::normalize_ads(_nw_m);
        }
    }

    Network_u init_nw_u() {
        this->_log->debug("Creating and initializing the user network ...");
        Network_u nw = Graph::create_graph<Network_u>(_cfg_u, *this->_rng);
        return nw;
    }

    Network_m init_nw_m() {
        if constexpr (model_mode == Media or Ageing_and_Media) {
            this->_log->debug("Creating and initializing the media network ...");
            Network_m nw = Graph::create_graph<Network_m>(_cfg_m, *this->_rng);
            return nw;
        }
    }

public:

    // Runtime functions ......................................................

    /** @brief Iterate a single step
     *  @detail Each step consists of (i) user revision: interaction between
     *  two users, (ii) information revision: interaction between a user and
     *  the media, (iii) media revision: interaction between two media.
     *  Media opinion revision can happen on timescales different to that of the user opinion revision.
     */
    void perform_step () {

        if constexpr (model_mode == None or Media) {
            revision::user_revision<Mode::None> (_nw_u,
                                     _weighting,
                                     _rewiring,
                                     _rewiring_count,
                                     _uniform_distr_prob_val,
                                     _radicalisation_parameter,
                                     *this->_rng);
        }

        if constexpr (model_mode == Media or Ageing_and_Media) {
            revision::information_revision (_nw_u,
                                            _nw_m,
                                            _uniform_distr_prob_val,
                                            _radicalisation_parameter,
                                            *this->_rng);

            if (this->get_time()%_media_time_constant==0) {
                      revision::media_revision(_nw_m, *this->_rng);
            }
        }

        //Perform user ageing once a year (= life_cycle numerical steps)
        if (model_mode == Ageing or Ageing_and_Media) {

            revision::user_revision<Mode::Ageing> (_nw_u,
                                         _weighting,
                                         _rewiring,
                                         _rewiring_count,
                                         _uniform_distr_prob_val,
                                         _radicalisation_parameter,
                                         *this->_rng);

            if (this->get_time()%_life_cycle==1) {
                ageing::ageing (_replacement_rate,
                                _num_media,
                                _child_ages,
                                _parent_ages,
                                _senior_ages,
                                _nw_u,
                                this->_log,
                                *this->_rng,
                                this->_cfg["susceptibility"]["users"]["custom"]);
            }
        }
    }

    /// Monitor model information
    /** @detail Here, functions and values can be supplied to the monitor that
     *          are then available to the frontend. The monitor() function is
     *          _only_ called if a certain emit interval has passed; thus, the
     *          performance hit is small.
     */
    void monitor ()
    {
        // Supply some number -- for illustration -- directly by value
        //this->_monitor.set_entry("some_value", 42);
    }


    /// Write data
    void write_data ()
    {

        /*
        auto get_edges_u = std::make_tuple(
            std::make_tuple("_sources",
                [](auto& _nw_u, auto& ed){
                    return boost::get(boost::vertex_index_t(),
                                                 _nw_u,
                                                 boost::source(ed, _nw_u));}
            ),
            std::make_tuple("_targets",
                [](auto& _nw_u, auto& ed){
                                      boost::get(boost::vertex_index_t(),
                                                 _nw_u,
                                                 boost::target(ed, _nw_u));}
            )
        );

        DataIO::save_graph_properties<Network_u::edge_descriptor>(_nw_u,
                        _grp_nw_u, std::to_string(get_time()), get_edges_u);
        */

        // Get iterators
        auto [v, v_end] = boost::vertices(_nw_u);
        auto [w, w_end] = boost::vertices(_nw_m);
        auto [e, e_end] = boost::edges(_nw_u);

        // opinion_u
        _dset_opinion_u->write( v, v_end,
                                [this](auto vd) {
                                return (float)_nw_u[vd].opinion;
                                });


        //tolerance
        _dset_tolerance_u->write(v, v_end,
                                      [this](auto vd) {
                                          return (float)_nw_u[vd].tolerance;
                                      });

        //susceptibility
        _dset_susceptibility_u->write(v, v_end,
                                      [this](auto vd) {
                                          return (float)_nw_u[vd].susceptibility;
                                      });


        if constexpr (model_mode == Ageing or Ageing_and_Media){
            //user age
            _dset_age_u->write(v, v_end,
                                          [this](auto vd) {
                                              return (unsigned int)_nw_u[vd].age;
                                      });
        }

        if constexpr (model_mode == Media or Ageing_and_Media) {
            // opinion_m
            _dset_opinion_m->write( w, w_end,
                                     [this](auto vd){
                                     return (float)_nw_m[vd].opinion;
                                     });

            // user numbers
            _dset_users->write( w, w_end,
                                     [this](auto vd) -> int {
                                     return (int)_nw_m[vd].users;
                                     });
        }

        // average neighbored opinion
        // _dset_avg_nb_opinion_u->write(v, v_end,
        //     [this](auto vd) -> double {
        //     double avg_nb_o = 0.;
        //     for (auto [nb, nb_end] = boost::adjacent_vertices(vd, _nw_u);
        //          nb!=nb_end; ++nb) {
        //         if (fabs(_nw_u[*nb].opinion - _nw_u[vd].opinion)
        //             <= _tolerance_u) {
        //             avg_nb_o += _nw_u[*nb].opinion
        //                         * _nw_u[boost::edge(vd, *nb, _nw_u).first].attr;
        //         }
        //         else {
        //             avg_nb_o += _nw_u[vd].opinion
        //                         * _nw_u[boost::edge(vd, *nb, _nw_u).first].attr;
        //         }
        //     }
        //     return (float)avg_nb_o;
        //     });

        // final edges of user network

        this->_log->debug("Writing {} edges ....", num_edges(_nw_u));

        if (this->get_time() + this->get_write_every() > this->get_time_max()) {
            _dset_edges_u_final->write(e, e_end,
                [&](auto ed){
                    return boost::get(boost::vertex_index_t(), _nw_u,
                                      boost::source(ed, _nw_u));
                }
            );

            _dset_edges_u_final->write(e, e_end,
                [&](auto ed){
                    return boost::get(boost::vertex_index_t(), _nw_u,
                                      boost::target(ed, _nw_u));
                }
            );

            this->_log->debug("All datasets have been written!");

            // std::vector<std::vector<size_t>> oc =
            //                 Graph_Analysis::opinion_clusters( _nw_u,
            //                                                  _tolerance_u);

            // _dset_final_opinion_clusters =
            //             _grp_nw_u->open_dataset("final_clusters", {oc.size()});

            // _dset_final_opinion_clusters->write(oc.begin(), oc.end(),
            //                     [&](auto c) -> double {
            //                     return (float)c.size();
            //                     });

            /*auto oc = Graph_Analysis::opinion_clusters( _nw_u,
                                                        _tolerance_u);
            _dset_num_opinion_clusters->write(oc.size());

            auto woc = Graph_Analysis::weighted_opinion_clusters(
                                                        _nw_u,
                                                        _tolerance_u);
            _dset_num_weighted_opinion_clusters->write(woc.size()); */ //this graph analysis uses old tolerance_u value, need to look at this in detail some time!!!! 24/10/19

            // auto rel_bc = Graph_Analysis::relative_betweenness_centrality(
            //                                                             _nw_u);

            // _dset_rel_bc->write(v, v_end,
            //                     [&](auto vd) -> double {
            //                     return rel_bc.at(vd);
            //                     });

            // _dset_weights->write(   e, e_end,
            //                         [this](auto ed) -> double {
            //                         return (float)_nw_u[ed].attr;
            //                         });
        }


        // weights
        // _dset_weights->write(   e, e_end,
        //                         [this](auto ed) -> double {
        //                         return (float)_nw_u[ed].attr;
        //                         });

        // ads
        //_dset_ads->write(  w, w_end,
        //                    [this](auto vd) -> double {
        //                    return (float)_nw_m[vd].ads;
        //                    });

        // rewiring count
        // _dset_rewiring_count->write(_rewiring_count);

        // _dset_in_degree->write(v, v_end,
        //                        [this](auto vd) -> double {
        //                        return in_degree(vd, _nw_u);
        //                        });

        // auto oc = Graph_Analysis::opinion_clusters( _nw_u,
        //                                             _tolerance_u);
        // _dset_num_opinion_clusters->write(oc.size());

        // auto woc = Graph_Analysis::weighted_opinion_clusters(
        //                                             _nw_u,
        //                                             _tolerance_u);
        // _dset_num_weighted_opinion_clusters->write(woc.size());
    }

    // Getters and setters ....................................................
    // Add getters and setters here to interface with other model
};

} //namespace

#endif // UTOPIA_MODELS_OPDYN_HH
