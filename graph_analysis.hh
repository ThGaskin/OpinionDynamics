#ifndef UTOPIA_MODELS_OPINIONET_GRAPH_ANALYSIS
#define UTOPIA_MODELS_OPINIONET_GRAPH_ANALYSIS

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/betweenness_centrality.hpp>
#include <boost/property_map/property_map.hpp>
#include <iostream>
#include <vector>
#include <algorithm>

namespace Utopia {
namespace Models {
namespace Opinionet {
namespace Graph_Analysis {

// HELPER FUNCTIONS ............................................................


// Starting from a given vertex, iteratively collect all vertices in tolerance
// range that are connected through an in edge or out edge.
template<typename NWType>
void fill_opinion_cluster(  size_t v,
                            std::vector<size_t>& c,
                            double tolerance,
                            NWType& nw) {

    if (std::find(c.begin(), c.end(), v) == c.end()) {
        c.push_back(v);
        for (auto [w, w_end]=adjacent_vertices(v, nw); w!=w_end; w++) {
            if (fabs(nw[v].opinion - nw[*w].opinion) <= nw[v].susceptibility) {
                fill_opinion_cluster(*w, c, nw[v].susceptibility, nw);
            }
        }
        for (auto [e, e_end]=in_edges(v, nw); e!=e_end; e++) {
            if (fabs(nw[v].opinion - nw[source(*e, nw)].opinion) <= nw[v].susceptibility) {
                fill_opinion_cluster(source(*e, nw), c, nw[v].susceptibility, nw);
            }
        }
    }
}


// Starting from a given vertex, iteratively collect all vertices in tolerance
// rangethat are connected through an in edge or out edge that have a weight
// larger than a certain threshold.
template<typename NWType>
void fill_weighted_opinion_cluster( size_t v,
                                    std::vector<size_t>& c,
                                    double min_weight,
                                    NWType& nw) {

    if (std::find(c.begin(), c.end(), v) == c.end()) {
        c.push_back(v);
        for (auto [w, w_end]=adjacent_vertices(v, nw); w!=w_end; w++) {
            if ((fabs(nw[v].opinion - nw[*w].opinion) <= nw[v].susceptibility)
                and (nw[edge(v, *w, nw).first].attr
                     * out_degree(v, nw) >= min_weight)) {
                fill_weighted_opinion_cluster(  *w,
                                                c,
                                                min_weight,
                                                nw);
            }
        }
        for (auto [e, e_end]=in_edges(v, nw); e!=e_end; e++) {
            if ((fabs(nw[v].opinion - nw[source(*e, nw)].opinion) <= nw[v].susceptibility)
                and (nw[*e].attr
                     * out_degree(source(*e, nw), nw) >= min_weight)) {
                fill_weighted_opinion_cluster(  source(*e, nw),
                                                c,
                                                min_weight,
                                                nw);
            }
        }
    }
}


// Starting from a given vertex, iteratively collect all adjacent vertices.
template<typename NWType>
void fill_community(size_t v,
                    std::vector<size_t>& c,
                    NWType& nw) {

    if (std::find(c.begin(), c.end(), v) == c.end()) {
        c.push_back(v);
        for (auto [w, w_end]=adjacent_vertices(v, nw); w!=w_end; w++) {
            fill_community(*w, c, nw);
        }
    }
}


// STRUCTURE ANALYSIS FUNCTIONS ................................................


// Calculate the reciprocity for a single node (= fraction of outgoing
// links for which the mutual link exists as well).
template<typename NWType, typename VertexDescType>
double reciprocity(NWType& nw, VertexDescType v) {
    double r = 0.;
    for (auto [w, w_end]=adjacent_vertices(v, nw); w!=w_end; ++w) {
        if (edge(*w, v, nw).second) {
            r += 1.;
        }
    }

    return r / double(out_degree(v, nw));
}


// Calculate the reciprocity of the whole graph (= fraction of mutual links).
template<typename NWType>
double reciprocity(NWType& nw) {
    double r = 0.;
    for (auto [e, e_end]=edges(nw); e!=e_end; ++e) {
        if (edge(target(*e, nw), source(*e, nw), nw).second) {
            r += 1.;
        }
    }

    return r / double(num_edges(nw));
}


// Calculate the betweenness centrality of each vertex.
template<typename NWType>
std::vector<double> betweenness_centrality(NWType& nw) {

    std::vector<double> centrality(num_vertices(nw));
    
    brandes_betweenness_centrality(
        nw,
        boost::make_iterator_property_map(  centrality.begin(),
                                            get(boost::vertex_index, nw),
                                            double())
        );

    return centrality;
}


// Calculate the relative betweenness centrality for each vertex
// (normalized with the highest possible value which would be reached
//  if a node is crossed by every single shortest path).
template<typename NWType>
std::vector<double> relative_betweenness_centrality(NWType& nw) {

    std::vector<double> centrality(num_vertices(nw));
    
    brandes_betweenness_centrality(
        nw,
        boost::make_iterator_property_map(  centrality.begin(),
                                            get(boost::vertex_index, nw),
                                            double())
        );

    relative_betweenness_centrality(
        nw,
        boost::make_iterator_property_map(  centrality.begin(),
                                            get(boost::vertex_index, nw),
                                            double())
        );

    // Division by 2 is needed for directed graphs.
    for (auto& val: centrality) {
        val /= 2.;
    }

    return centrality;
}


// Identify groups of agents with similar (within tolerance range) opinions.
template<typename NWType>
std::vector<std::vector<size_t>> opinion_groups(NWType nw,
                                                double tolerance) {

    // First, get pairs of opinion values and vertices
    std::vector<std::pair<double, size_t>> op_v;
    std::vector<std::vector<size_t>> groups;

    for (auto [v, v_end]=vertices(nw); v!=v_end; v++) {
        op_v.push_back(std::make_pair(nw[*v].opinion, *v));
    }

    // sort along the opinion values
    std::sort(op_v.begin(), op_v.end());

    // loop over opinions and make a cut wherever the opinion distance
    // is larger than the tolerance range
    size_t start = 0;
    for (size_t i=0; i!=op_v.size()-1; i++) {

        if (fabs(op_v[i].first - op_v[i+1].first) >= tolerance) {
            std::vector<size_t> group;

            for (size_t j=start; j!=i+1; j++) {
                group.push_back(op_v[j].second);
            }

            start = i+1;
            groups.push_back(group);
        }
    }

    // Add last group
    std::vector<size_t> group;
    for (size_t j=start; j!=op_v.size(); j++) {
        group.push_back(op_v[j].second);
    }
    groups.push_back(group);

    return groups;
}


// Identify groups of agents with similar (within tolerance range) opinions
// that are connected on the network.
template<typename NWType>
std::vector<std::vector<size_t>> opinion_clusters(NWType nw,
                                                  double tolerance) {
    
    std::vector<std::vector<size_t>> opinion_clusters;
    std::vector<size_t> temp_c;
    bool next;

    // Find all opinion clusters through a loop over all vertices as the
    // source of the cluster.
    for (auto [v, v_end]=vertices(nw); v!=v_end; v++) {

        next = false;

        // If vertex is part of an already discovered opinion cluster
        // its cluster is the same (by definition).
        for (auto& c: opinion_clusters) {
            if (std::find(c.begin(), c.end(), *v) != c.end()) {
                next = true;
            }
        }

        if (next) {
            continue;
        }

        // Else get the opinion cluster around this vertex.
        else {

            temp_c.clear();
            fill_opinion_cluster(*v, temp_c, tolerance, nw);
            opinion_clusters.push_back(temp_c);
        }
    }

    return opinion_clusters;
}


// Identify groups of agents with similar (within tolerance range) opinions
// that are connected on the network (with in or out edges that have a weight
// larger than a certain threshold).
template<typename NWType>
std::vector<std::vector<size_t>> weighted_opinion_clusters( 
                                    NWType nw,
                                    double tolerance,
                                    double min_weight = -1.) {
    
    std::vector<std::vector<size_t>> opinion_clusters;
    std::vector<size_t> temp_c;
    bool next;

    if (min_weight < 0.) {
        min_weight = 0.1;
    }

    // Find all opinion clusters through a loop over all vertices as the
    // source of the cluster.
    for (auto [v, v_end]=vertices(nw); v!=v_end; v++) {

        next = false;

        // If vertex is part of an already discovered opinion cluster
        // its cluster is the same (by definition).
        for (auto& c: opinion_clusters) {
            if (std::find(c.begin(), c.end(), *v) != c.end()) {
                next = true;
            }
        }

        if (next) {
            continue;
        }

        // Else get the opinion cluster around this vertex.
        else {

            temp_c.clear();
            fill_weighted_opinion_cluster(  *v,
                                            temp_c,
                                            min_weight,
                                            nw);
            opinion_clusters.push_back(temp_c);
        }
    }

    return opinion_clusters;
}


// Identify groups of agents that are connected via out-edges.
// NOTE that completely isolated vertices are also identified
//      as closed community.
template<typename NWType>
std::vector<std::vector<size_t>> closed_communities(NWType nw) {
    
    std::vector<std::vector<size_t>> cc;
    std::vector<size_t> temp_c;
    bool next;

    // Find all communities through a loop over all vertices as the
    // source of the community.
    for (auto [v, v_end]=vertices(nw); v!=v_end; v++) {

        next = false;

        // If vertex is part of an already discovered community
        // its community has to be the same (if out-degree > 0).
        for (auto& c: cc) {
            if (std::find(c.begin(), c.end(), *v) != c.end()) {
                next = true;
            }
        }

        // This is the case of a 'loner'.
        if (not next) {
            if (in_degree(*v, nw) < 2) {
                for (auto& c: cc) {
                    for (auto& w: c) {
                        if (edge(*v, w, nw).second) {
                            c.push_back(*v);
                            next = true;
                            break;
                        }
                    }
                }
            }
        }

        if (next) {
            continue;
        }

        // Else get the community originating from the vertex.
        else {

            temp_c.clear();
            fill_community(*v, temp_c, nw);
            cc.push_back(temp_c);
        }
    }

    return cc;
}


} // namespace Graph_Analysis
} // namespace Opinionet
} // namespace Models
} // namespace Utopia

#endif // UTOPIA_MODELS_OPINIONET_GRAPH_ANALYSIS
