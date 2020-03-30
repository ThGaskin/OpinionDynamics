#ifndef UTOPIA_MODELS_OPDYN_REVISION
#define UTOPIA_MODELS_OPDYN_REVISION

#include <cmath>
#include <algorithm>
#include <iostream>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <stdlib.h>
#include <spdlog/spdlog.h>

#include "modes.hh"
#include "update.hh"
#include "utils.hh"

namespace Utopia::Models::OpDyn::revision {

using modes::Mode;
using modes::Mode::Ageing;
using modes::Mode::Ageing_and_Media;

// STEP FUNCTIONS ..............................................................

template <typename VertexDescType, typename NWType, typename RNGType>
void pairwise_weighted_update(
                VertexDescType& v,
                NWType& nw,
                std::uniform_real_distribution<double>& prob_distr,
                RNGType& rng,
                const double radicalisation_parameter)
{
    auto nb = v;

// Choose interaction partner. The probability for choosing neighbour w
// is given by the weight on the edge (v, w).

    double nb_prob_frac = prob_distr(rng);
    double cumulative_weights = 0.;
    for (auto [w, w_end] = adjacent_vertices(v, nw); w!=w_end; ++w) {
        if (cumulative_weights < nb_prob_frac) {
            cumulative_weights += nw[edge(v, *w, nw).first].attr;
        }
        if (cumulative_weights >= nb_prob_frac) {
            nb = *w;
            break;
        }
    }

// The opinion update takes the tolerance and the susceptibility into account.
// The tolerance is updated by considering the change in opinion distance
// from the middle 0.5; the susceptibility decreases with age (if set)

    double old_opinion = nw[v].opinion;

    //Opinion update
    update::opinion(v, nb, nw);

    //Tolerance update
    update::tolerance(v, nw, old_opinion, radicalisation_parameter);

}

// The weights are updated proportionally to the distance from a neighbour's
// opinion to the user's current opinion.
template<Mode model_mode, typename NWType, typename VertexDescType, typename RNGType>
void update_weights(VertexDescType v,
                    NWType& nw,
                    const double weighting,
                    const double rewiring,
                    unsigned int rewiring_count,
                    std::uniform_real_distribution<double> prob_distr,
                    RNGType& rng)
{

    if (out_degree(v, nw) != 0) {
        // Reduce weight dependent on the opinion distance and identify
        // edges that will be cut.
        std::vector<VertexDescType> to_drop;
        double sum_of_reduced_weights = 0.;

        auto [e, e_end] = out_edges(v, nw);
        for (auto next = e; e!=e_end; e = next) {
            ++next;

            // If opinion distance is larger than tolerance, rewire with
            // probability 'rewiring'.
            if (fabs(nw[target(*e, nw)].opinion - nw[source(*e, nw)].opinion)
                                > nw[source(*e, nw)].tolerance) {
                if (prob_distr(rng) < rewiring) {
                    to_drop.push_back(target(*e, nw));
                }
            }

            // Change weight proportionally to the opinion distance and the age difference.
            // Both factors are weighted by 50%.
            // NOTE that for weighting > 1, weights can reach Zero.
            if constexpr (model_mode == Mode::Ageing or Mode::Ageing_and_Media) {
              nw[*e].attr *= (1. - weighting
                              * fabs(nw[target(*e, nw)].opinion - nw[v].opinion))
                              + std::exp(std::log(0.5)/0.5
                              *fabs(nw[target(*e, nw)].age - nw[v].age)/nw[v].age);
            }
            else {
              nw[*e].attr *= (1. - weighting
                              * fabs(nw[target(*e, nw)].opinion - nw[v].opinion));
            }

            if (nw[*e].attr < 0.) {
                nw[*e].attr = 0.;
            }

            sum_of_reduced_weights += nw[*e].attr;
        }

        // Try to find suitable new neighbors for rewiring.
        std::vector<VertexDescType> to_add;
        for (size_t i=0; i!=to_drop.size(); i++) {

            auto w = utils::get_rand_nb(nw, v, rng);
            if (std::find(  to_drop.begin(),
                            to_drop.end(), w) == to_drop.end()) {
                if (out_degree(w, nw) != 0) {
                    w = utils::get_rand_nb(nw, w, rng);
                }
                if (edge(v, w, nw).second or (v==w)) {
                    w = random_vertex(nw, rng);
                }
            }
            else {
                w = random_vertex(nw, rng);
            }

            if ((not edge(v, w, nw).second)
                and (std::find( to_add.begin(),
                                to_add.end(), w) == to_add.end())
                and (v!=w)) {

                to_add.push_back(w);
                sum_of_reduced_weights -=
                                nw[edge(v, to_drop[i], nw).first].attr;
                remove_edge(v, to_drop[i], nw);
            }
        }

        // Determine the initial weight which is given to the added edges.
        double init_weight = 0.;
        if (out_degree(v, nw) != 0) {
            // precision threshold due to possible rounding errors
            if (sum_of_reduced_weights < 10e-5) {
                init_weight = 1. / double(to_add.size());
            }
            else {
                init_weight = sum_of_reduced_weights
                                / double(out_degree(v, nw));
            }
        }
        else {
            init_weight = 1. / double(to_add.size());
        }

        for (size_t i=0; i!=to_add.size(); i++) {
            add_edge(v, to_add[i], {init_weight}, nw);
            rewiring_count ++;
        }
    }
}


// This function normalizes the weights of vertex v to 1.
template<typename VertexDescType, typename NWType>
void normalize_weights(VertexDescType v, NWType& nw) {
    const auto log = spdlog::get("root.OpDyn");

    if (out_degree(v, nw) != 0) {
        double weight_norm = 0.;
        for (auto [e, e_end] = out_edges(v, nw); e!=e_end; ++e) {
            weight_norm += nw[*e].attr;
        }

        if (weight_norm != 0.) {
            for (auto [e, e_end] = out_edges(v, nw); e!=e_end; ++e) {
                  nw[*e].attr /= weight_norm;
                // This should be part of a test case...
                if (std::isnan(nw[*e].attr)) {
                  log->error("NAN weight!");
                }
            }
        }
        else {
            log->warn("All weights are Zero! This node's age: {}", nw[v].age );
        }
    }
}

double make_periodic(double val) {
    if (val < -0.5) {
        return -1.;
    }
    else if (val > 0.5) {
        return 1.;
    }
    else {
        return 0.;
    }
}

// returns absolute distance |x-y| for periodic boundaries
double distance_periodic(double x, double y) {
    return fabs(x - y - make_periodic(x - y));
}

// normalize ad values so that the ad fractions represent
// interaction probabilities
template<typename NWType_m>
void normalize_ads(NWType_m& nw_m) {
    double sum = 0.;
    for (auto [v, v_end] = vertices(nw_m); v!=v_end; ++v) {
        sum += nw_m[*v].ads;
    }
    for (auto [v, v_end] = vertices(nw_m); v!=v_end; ++v) {
        nw_m[*v].ads_normalized = nw_m[*v].ads/sum;
    }
}

// This is the bounded confidence interaction.
std::pair<double,bool> user_char_BC(    double own_opinion,
                                        double new_opinion,
                                        double tolerance,
                                        bool periodic = false)
{
    if (periodic) {
        if (distance_periodic(new_opinion, own_opinion) <= tolerance) {
            return std::make_pair(1., true);
        }
        else {
            return std::make_pair(0., false);
        }
    }
    else {
        if (fabs(new_opinion - own_opinion) <= tolerance) {
            return std::make_pair(1., true);
        }
        else {
            return std::make_pair(0., false);
        }
    }
}



// Main revision functions .....................................................
// In a single time step the update consists of three separate revision
// processes: user-revision, information-revision (between users and media)
// and media-revision.

template<Mode model_mode, typename NWType, typename RNGType>
void user_revision( NWType& nw_u,
                    double weighting,
                    double rewiring,
                    unsigned int& rewiring_count,
                    std::uniform_real_distribution<double> prob_distr,
                    double radicalisation_parameter,
                    RNGType& rng) {

    // choose random vertex that gets a revision opportunity
    auto v = random_vertex(nw_u, rng);

    if (out_degree(v, nw_u) != 0) {

        // pairwise opinion update with bounded confidence
        pairwise_weighted_update(   v,
                                    nw_u,
                                    prob_distr,
                                    rng,
                                    radicalisation_parameter);

//update the weights depending on the opinion distance
        update_weights<model_mode>( v,
                        nw_u,
                        weighting,
                        rewiring,
                        rewiring_count,
                        prob_distr,
                        rng);


        normalize_weights(v, nw_u);

    }
}

template<typename NWType_m, typename RNGType>
void media_revision(NWType_m& nw_m, RNGType& rng) {

// choose random vertex for revision
    auto v = random_vertex(nw_m, rng);

// exponential decay of advertisement impact
    nw_m[v].ads *= 0.9;

// Find the neighbour that is both ideologically closest and has the largest
//user count and adjust. If medium v is already the optimum, there is no change.

    if (out_degree(v, nw_m) != 0) {

        auto fittest_nb = v;
        bool found = false;
        for (auto [w, w_end] = adjacent_vertices(v, nw_m); w!=w_end; ++w) {
            if (fabs(nw_m[v].opinion - nw_m[*w].opinion) <= nw_m[v].tolerance
                && nw_m[*w].users > nw_m[fittest_nb].users) {
                  fittest_nb = *w;
                  found = true;
            }
        }
// a medium shifts its stance towards the more popular competitor.
        if (found) {
            if(fabs(nw_m[v].opinion-nw_m[fittest_nb].opinion)>nw_m[v].tolerance/3) {
              nw_m[v].opinion += nw_m[v].susceptibility
                              * nw_m[edge(v, fittest_nb, nw_m).first].attr
                              * (nw_m[fittest_nb].opinion - nw_m[v].opinion);
            }
            else {
              int sgn = 1;
              if (nw_m[v].opinion != nw_m[fittest_nb].opinion
                  && nw_m[v].opinion < nw_m[fittest_nb].opinion) {
                  sgn = -1;
              }
              if (nw_m[v].opinion == nw_m[fittest_nb].opinion) {
                sgn = rand() % 2;
                sgn = 2*sgn-1;
              }
              nw_m[v].opinion = nw_m[fittest_nb].opinion +
                                sgn * nw_m[v].tolerance/3;
            }
        }
    }

    if (nw_m[v].opinion < 0.) {
        nw_m[v].opinion = 0.;
    }
    if (nw_m[v].opinion > 1.) {
        nw_m[v].opinion = 1.;
    }

//Advertising:
// The more users a medium has, the more it can spend on ads,
//which increases its presence
    nw_m[v].ads = nw_m[v].users;
// normalize advertisment values to probabilities
    normalize_ads(nw_m);
// store current user number for comparison at the next revision opportunity
    nw_m[v].users_previous = nw_m[v].users;
}

    enum class UserChartype {
        bc,
        bc_extended,
        gaussian
    };

template<typename NWType_u, typename NWType_m, typename RNGType>
void information_revision(  NWType_u& nw_u,
                            NWType_m& nw_m,
                            std::uniform_real_distribution<double> prob_distr,
                            const double radicalisation_parameter,
                            RNGType& rng) {

    auto v = random_vertex(nw_u, rng);
    size_t new_medium = 0;
    bool found = false;

    // choose new medium. The probability for choosing medium i is given by the
    // ad-fraction of medium i.
    double new_medium_ad_fraction = prob_distr(rng);
    double sum_ad_fraction = 0.;
    for (auto [m, m_end] = vertices(nw_m); m!=m_end; ++m) {
        if (sum_ad_fraction < new_medium_ad_fraction) {
            sum_ad_fraction += nw_m[*m].ads_normalized;
        }
        if (sum_ad_fraction >= new_medium_ad_fraction) {
            new_medium = *m;
            found = true;
            break;
        }
    }

    std::pair<double, bool> user_char = user_char_BC(   nw_u[v].opinion,
                                    nw_m[new_medium].opinion,
                                    nw_u[v].tolerance);

// At first, check if the new medium is affordable for the user.
// If so, switch to the new medium. The user is also influenced in her
// opinion if the medium's opinion is close enough.
// Users changing their opinions also leads to a change in user tolerance
    double opinion_old = nw_u[v].opinion;
    if (prob_distr(rng) <= user_char.first) {
        if (user_char.second) {
            update::opinion(v, new_medium, nw_u, nw_m);
        }

        update::tolerance(v, nw_u, opinion_old, radicalisation_parameter);

        nw_m[nw_u[v].used_media].users -= 1;
        nw_m[new_medium].users += 1;
        nw_u[v].used_media = new_medium;
    }
}

} // namespace

#endif // UTOPIA_MODELS_OPDYN_REVISION
