#define BOOST_TEST_MODULE test utils

#include <random>
#include <type_traits>

#include <boost/test/unit_test.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/random.hpp>

#include <utopia/core/model.hh>
#include <utopia/core/types.hh>
#include <utopia/core/graph.hh>
#include <utopia/data_io/graph_utils.hh>

#include "../update.hh"
#include "../utils.hh"
#include "../revision.hh"
#include "../OpDyn.hh"

namespace Utopia::Models::OpDyn {

// -- Type definitions --------------------------------------------------------

using vertex = boost::graph_traits<Network_u>::vertex_descriptor;

std::mt19937 rng{};

// -- Fixtures ----------------------------------------------------------------

// Test random network
struct TestNetwork {
    Network_u nw;

    TestNetwork(){
        const unsigned num_vertices = 200;
        const unsigned num_edges = 1000;
        constexpr bool allow_parallel = false;
        constexpr bool allow_self_edges = false;

        boost::generate_random_graph(nw,
                                 num_vertices,
                                 num_edges,
                                 rng,
                                 allow_parallel,
                                 allow_self_edges);
    }
};


// -- Actual test -------------------------------------------------------------

BOOST_FIXTURE_TEST_CASE( test_inits,
                         TestNetwork,
                         * boost::unit_test::tolerance(1e-12))
{

    // 1. Test initialisation functions .........................................
    {
        for (auto v : range<IterateOver::vertices>(nw)) {

            std::pair<double, double> interval=std::make_pair(.8, .8);
            nw[v].tolerance=utils::set_init_uniform(interval, rng);
            BOOST_TEST(nw[v].tolerance == .8);

            std::pair<double, double> opinion_interval=std::make_pair(0., 1.);
            nw[v].opinion=utils::set_init_uniform(opinion_interval, rng);
            BOOST_TEST(0<=nw[v].opinion);
            BOOST_TEST(nw[v].opinion<=1);

            std::pair<double, double> distribution_info=std::make_pair(0.5, 0.2);
            nw[v].susceptibility=utils::set_init_Gauss(distribution_info, rng);
            BOOST_TEST(0<=nw[v].susceptibility);
            BOOST_TEST(nw[v].susceptibility<=1);

            std::pair<int, int> age_interval=std::make_pair(45, 45);
            nw[v].age=utils::set_init_uniform(age_interval, rng);
            BOOST_TEST(nw[v].age == 45);

            age_interval=std::make_pair(1, 100);
            nw[v].age=utils::set_init_uniform(age_interval, rng);
            BOOST_TEST(nw[v].age >=1);
            BOOST_TEST(nw[v].age <=100);

            double weight_sum = 0;
            for (auto e : range<IterateOver::out_edges>(v, nw)) {
                nw[e].attr = 1. / (double) boost::out_degree(v, nw);
                weight_sum+=nw[e].attr;
            }
            if (boost::out_degree(v, nw)!=0) {
                BOOST_TEST(weight_sum == 1.);
            }
        }
    }
}

//2. Test opinion update........................................................

BOOST_FIXTURE_TEST_CASE( test_opinion_update,
                         TestNetwork,
                         * boost::unit_test::tolerance(1e-12)) {

    {
        vertex v = 2;
        vertex nb = 3;

        nw[v].opinion = 0.5;
        nw[nb].opinion = 1.;
        nw[v].susceptibility=1.;
        nw[v].tolerance=1.;
        nw[nb].tolerance=.25;
        update::opinion(v, nb, nw);
        BOOST_TEST(nw[v].opinion == 1.);

        nw[nb].opinion = 0.;
        nw[v].susceptibility = .5;
        update::opinion(v, nb, nw);
        BOOST_TEST(nw[v].opinion == .5);

        update::opinion(v, nb, nw);
        BOOST_TEST(nw[v].opinion == 0.25);

        nw[nb].susceptibility=0.1;
        update::opinion(nb, v, nw);
        BOOST_TEST(nw[nb].opinion == 0.025);
    }
}

//3. Test tolerance update.....................................................
BOOST_FIXTURE_TEST_CASE( test_tolerance_update,
                         TestNetwork,
                         * boost::unit_test::tolerance(1e-12)) {


    {
        double k=2; //the radicalisation parameter
        vertex v = 4;
        vertex nb = 5;

        nw[v].opinion = 0.5;
        double previous_opinion = nw[v].opinion;
        nw[nb].opinion = 1.;
        nw[v].tolerance = 0.5;
        nw[v].susceptibility=1.;

        update::opinion(v, nb, nw);
        update::tolerance(v, nw, previous_opinion, k);
        BOOST_TEST(nw[v].tolerance==pow(0.5, 1.5));

        nw[nb].opinion=0.5+1/pow(12, 0.5);
        previous_opinion = nw[v].opinion;
        update::opinion(v, nb, nw);
        update::tolerance(v, nw, previous_opinion, k);
        BOOST_TEST(nw[v].tolerance==0.5);

        k=0;
        update::opinion(v, nb, nw);
        update::tolerance(v, nw, previous_opinion, k);
        BOOST_TEST(nw[v].tolerance==0.5);

        k=4;
        v=6;
        nb=7;
        nw[v].opinion=0;
        nw[v].tolerance=0.5;
        nw[v].susceptibility=1;
        nw[nb].opinion=0.5;
        update::opinion(v, nb, nw);
        BOOST_TEST(nw[v].opinion==0.5);
        previous_opinion=0;
        update::tolerance(v, nw, previous_opinion, k);
        BOOST_TEST(nw[v].tolerance==1);
    }

}

} // namespace Utopia::Models::OpDyn
