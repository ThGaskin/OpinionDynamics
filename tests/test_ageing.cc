#define BOOST_TEST_MODULE test ageing

#include <random>
#include <type_traits>

#include <boost/test/unit_test.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/random.hpp>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include <utopia/core/model.hh>
#include <utopia/core/apply.hh>
#include <utopia/core/types.hh>
#include <utopia/core/graph.hh>
#include <utopia/data_io/graph_utils.hh>
#include <utopia/data_io/cfg_utils.hh>


#include "../utils.hh"
#include "../revision.hh"
#include "../OpDyn.hh"


namespace Utopia::Models::OpDyn {

// -- Type definitions --------------------------------------------------------



// -- Fixtures ----------------------------------------------------------------


/// A fixture that provides the information necessary to setup consumers
/** \detail All members of this are available in the test cases.
  */


/// Creates the root.OpDyn logger, which is used internally.
/** While this logger is usually set up by the model itself, this is not the
  * case within this test scope. Thus, it needs to be set up manually,
  * otherwise it might lead to segfaults when attempting to use the logger.
  */
auto create_OpDyn_logger(const std::string name="root.OpDyn") {
    // Try to get it from the logger registry
    auto logger = spdlog::get(name);

    if (not logger) {
        // Does not exist yet. Create it.
        logger = spdlog::stdout_color_mt(name);
    }
    return logger;
}


/// Test random network
struct TestNetwork {
    using vertex = boost::graph_traits<Network_u>::vertex_descriptor;
    using RNG = std::mt19937;
    using Config = Utopia::DataIO::Config;

    std::shared_ptr<spdlog::logger> log;
    std::shared_ptr<RNG> rng;
    Network_u nw;
    Config cfg;

    TestNetwork()
    :
        // Set up a dedicated test logger
        log([](){
            auto logger = spdlog::get("test_ageing");

            // Create it only if it does not already exist
            if (not logger) {
                logger = spdlog::stdout_color_mt("test_ageing");
            }

            // Set level and global logging pattern
            logger->set_level(spdlog::level::debug);
            spdlog::set_pattern("[%T.%e] [%^%l%$] [%n]  %v");
            // "[HH:MM:SS.mmm] [level(colored)] [logger]  <message>"

            return logger;
        }()),
        rng(std::make_shared<RNG>(std::random_device()())),
        nw{},
        cfg(YAML::LoadFile("test_config.yml"))
    {
        // Initialize the logger that is used within OpDyn.
        // It need not be stored here because within OpDyn it is retrieved
        // directly from the spdlog logger registry.
        create_OpDyn_logger();

        // Create a test network
        const unsigned num_vertices = 2000;
        const unsigned num_edges = 40000;
        constexpr bool allow_parallel = false;
        constexpr bool allow_self_edges = false;

        BOOST_TEST_MESSAGE("Generating random graph ...");
        boost::generate_random_graph(nw,
                                     num_vertices,
                                     num_edges,
                                     *rng,
                                     allow_parallel,
                                     allow_self_edges);

        BOOST_TEST_MESSAGE("Initializing ages ...");
        for (auto v : range<IterateOver::vertices>(nw)){
            nw[v].age=utils::set_init_uniform(std::make_pair(1, 91), *rng);
        }
    }
};


// -- Tests -------------------------------------------------------------------
// ... using the TestNetwork fixture throughout all test cases
BOOST_FIXTURE_TEST_SUITE(network_suite, TestNetwork)


BOOST_AUTO_TEST_CASE(test_user_collection)
{
    using ageing::user_selection_and_ageing;

    std::pair<int, int> child_ages(0, 10);
    std::pair<int, int> parent_ages(20, 40);
    std::pair<int, int> senior_ages(70, 1000);
    std::vector<vertex> children, parents, peers;
    double replacement_rate=0.01;
    std::cout<<cfg["susceptibility"]["users"]["custom"]["peak"]<<std::endl;
    user_selection_and_ageing(children, parents, peers, child_ages,
                              parent_ages, senior_ages, replacement_rate,
                              nw, *rng, cfg["susceptibility"]["users"]["custom"]);

    BOOST_TEST(children.size()==20);
    BOOST_TEST(parents.size()==20);

    for (int i=0; i<children.size(); ++i) {
        BOOST_TEST(nw[children.at(i)].age>=70);
        BOOST_TEST(nw[parents.at(i)].age>=20);
        BOOST_TEST(nw[parents.at(i)].age<=40);
    }
    for (int i=0; i<peers.size(); ++i){
        BOOST_TEST(nw[peers.at(i)].age<=10);
        BOOST_TEST(nw[peers.at(i)].age>0);
    }
}

BOOST_AUTO_TEST_CASE(test_vertex_removal)
{
    using ageing::user_selection_and_ageing;
    using ageing::remove_edges;

    std::pair<int, int> child_ages(0, 10);
    std::pair<int, int> parent_ages(20, 40);
    std::pair<int, int> senior_ages(70, 1000);
    std::vector<vertex> children, parents, peers;
    double replacement_rate = 0.1;

    BOOST_TEST_PASSPOINT();
    user_selection_and_ageing(children, parents, peers, child_ages,
                              parent_ages, senior_ages, replacement_rate,
                              nw, *rng, cfg["susceptibility"]["users"]["custom"]);

    BOOST_TEST(children.size() == 200);
    BOOST_TEST(parents.size() == 200);

    for (int i=0; i<children.size(); ++i){
        BOOST_TEST_CONTEXT("Child: " << i) {
            BOOST_TEST_CHECKPOINT("Removing all edges of vertex ...");
            remove_edges(children.at(i), nw);

            BOOST_TEST_CHECKPOINT("Asserting degree is now zero ...");
            BOOST_TEST(boost::out_degree(children.at(i), nw) == 0);
            BOOST_TEST(boost::in_degree(children.at(i), nw) == 0);
        }
    }
}


BOOST_AUTO_TEST_SUITE_END()

} // namespace Utopia::Models::OpDyn
