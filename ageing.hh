#ifndef UTOPIA_MODELS_OPDYN_AGEING
#define UTOPIA_MODELS_OPDYN_AGEING

#include <cmath>
#include <algorithm>
#include <iostream>

#include <spdlog/spdlog.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/assert.hpp>

#include "utils.hh"
#include "revision.hh"

namespace Utopia::Models::OpDyn::ageing {

using pair_int = std::pair<int, int>;

// USER AGEING .................................................................
template <typename VertexDescType, typename NWType, typename RNGType, typename Config>
void user_selection_and_ageing( std::vector<VertexDescType> &children,
                                std::vector<VertexDescType> &parents,
                                std::vector<VertexDescType> &peers,
                                const pair_int child_ages,
                                const pair_int parent_ages,
                                const pair_int senior_ages,
                                const double replacement_rate,
                                NWType& nw,
                                RNGType& rng,
                                const Config& cfg) {


    /*!
    This function collects the nodes that are to be reinitialised as nodes of age 1, and finds an adequate number
     of peers to rewire to. It concurrently increases the age of each user by 1.
    */

    const int vertices_to_remove = (int)(boost::num_vertices(nw)*replacement_rate);
    int peers_to_add = 0;

    /*find old nodes to reinitialise as children, find an equal number of
    parents, and collect a sufficient number of young peers to reconnect to the
    children */
    apply_rule<IterateOver::vertices, Update::async, Shuffle::on>(
    [&children, &parents, &peers, &child_ages, &parent_ages, &senior_ages, &vertices_to_remove, &peers_to_add, &cfg]
    (const VertexDescType v, NWType& nw){

        /* increase the age of every user
        and adjust the susceptibility accordingly */
        ++nw[v].age;
        nw[v].susceptibility=utils::susceptibility(cfg, nw[v].age);

        if(peers.size() < peers_to_add
           or children.size() < vertices_to_remove
           or parents.size() < vertices_to_remove) {

              //nothing happens for users who are neither children nor parents nor seniors
               if ((child_ages.second<nw[v].age && nw[v].age<parent_ages.first)
                   or (parent_ages.second<nw[v].age && nw[v].age<senior_ages.first)) {
               }

               //seniors are reinitialised as children
            if (senior_ages.first<=nw[v].age && children.size()<vertices_to_remove) {
                   children.push_back(v);
                   /* because we also have a parent with an in- and out-edge,
                   we need two fewer peers per child than before */
                   peers_to_add+=degree(v,nw)-2;
               }

               //select parents
               if (parent_ages.first<=nw[v].age && nw[v].age<=parent_ages.second &&
                   parents.size()<children.size()) {
                   parents.push_back(v);
               }

               //find peers to reconnect new children to
               if (nw[v].age<=child_ages.second && peers.size()<peers_to_add) {
                   peers.push_back(v);
               }
        }
    },
    nw,
    rng
    );
}

/// Removes in- and out-edges to old vertex and normalise the previous peers' weights
template <typename VertexDescType, typename NWType>
void remove_edges (VertexDescType v, NWType& nw) {
    for (auto e : range<IterateOver::in_edges>(v, nw)) {
        VertexDescType w = source(e, nw);
        nw[e].attr = 0.;

        if (out_degree(w, nw) > 1) {
            revision::normalize_weights(w, nw);
        }
    }
    clear_vertex(v, nw);
}


template <typename VertexDescType, typename NWType, typename LoggerType, typename RNGType>
void add_edges( const VertexDescType child,
                const VertexDescType parent,
                const std::vector<VertexDescType> peers,
                int out_deg, int in_deg, int deg,
                int iter_number,
                NWType& nw,
                LoggerType& log,
                RNGType& rng) {

    /// Add edges to new peers in such a way that degree is preserved

    /*add child-parent edge and set weight to 1 if no other out-edges will be
    added, 0.5 else */
    double init_weight = 0.5;
    if (deg <= 2 or out_deg<=1) {init_weight = 1.;}
    add_edge(child, parent, {init_weight}, nw);

    //keep track of how many peers still need to be added
    if(out_deg>0) {--out_deg;}
    else {--in_deg;}

    if(in_deg>0 or out_deg>0) {
      add_edge(parent, child, {0.1}, nw);
      revision::normalize_weights(parent, nw);
    }

    // else, degree=1, and we are done
    else {
      revision::normalize_weights(parent, nw);
      return;
    }

    if(in_deg>0) {--in_deg;}
    else {--out_deg;}

    //if we have used up all available degrees, we are done
    if(in_deg<=0 && out_deg <= 0) {return;}

    /* If we have got to this step, there are spare edge degrees that can be
    rewired to peers. The child opinion will in this case be set to
    0.5*parent.opinion + 0.5 average peer opinion, so we need to collect the
    opinions of the peers we are rewiring to */
    double peer_opinions = 0.;
    bool rewire_fail=false;

    //add edge to each peer, collect peer opinion and set weight
    //out-edges
    for (int j=iter_number; j<out_deg+iter_number; ++j) {
          auto peer = peers.at(j%peers.size());
          while(edge(child, peer, nw).second) {
            peer=random_vertex(nw, rng);
            rewire_fail = true;
          }
          add_edge(child, peer, {0.5/out_deg}, nw);
          revision::normalize_weights(peer, nw);
          peer_opinions += nw[peer].opinion;
    }

    //in-edges
    for (int j=out_deg+iter_number; j<in_deg+out_deg+iter_number; ++j) {
          auto peer = peers.at(j%peers.size());
          while(edge(peer, child, nw).second) {
            peer=random_vertex(nw, rng );
            rewire_fail = true;
          }
          add_edge(peer, child, {0.5/in_deg}, nw);
          revision::normalize_weights(peer, nw);
          peer_opinions += nw[peer].opinion;
    }

    //child opinion = 50% parent opinion + 50% peer average
    nw[child].opinion = 0.5*(nw[parent].opinion + peer_opinions/(in_deg+out_deg));

    /* For low vertex numbers, there may not be enough different peers to rewire
    to. In this case, a random vertex must be picked from the remaining age
    groups to preserve the edge count */
    if(rewire_fail)  {
      log->warn("Failed connecting to peer: edge already exists. Rewiring to a \
random vertex. If warning persists, consider increasing vertex \
count or decreasing replacement rate!");
    }
}

template <typename VertexDescType, typename RNGType, typename NWType, typename Config>
void reinitialize( VertexDescType child,
                   VertexDescType parent,
                   int num_media,
                   NWType& nw,
                   RNGType& rng,
                   Config& cfg){

    nw[child].age=1;

    /* Set the child's opinion to the parent opinion,
    which will be changed if the child has peers */
    nw[child].opinion = nw[parent].opinion;

    nw[child].tolerance = nw[parent].tolerance;

    /* Set the child susceptibility to the value for the susceptibility
    function at age 1 */
    nw[child].susceptibility=utils::susceptibility(cfg, 1);
    nw[child].used_media=utils::get_rand_int(0, num_media, rng);
}

//++++DEBUGGING TESTS++++.......................................................
// none of these should appear, and will be moved to testing
template <typename VertexDescType, typename NWType>
void check_and_test( const VertexDescType child,
                     const VertexDescType parent,
                     const int deg,
                     const int in_deg,
                     const int out_deg,
                     NWType& nw) {

    const auto log = spdlog::get("root");

    const int deg_after = degree(child, nw);
    const int in_deg_after = in_degree(child, nw);
    const int out_deg_after = out_degree(child, nw);

    BOOST_ASSERT((deg==deg_after));
    BOOST_ASSERT((std::abs(in_deg-in_deg_after)<=1));
    BOOST_ASSERT((std::abs(out_deg-out_deg_after)<=1));

    if (deg!=0) {assert((out_deg_after>=1));}

    double weight_sum=0.;
    for(auto [e, e_end]=out_edges(child, nw); e!=e_end; ++e){
      weight_sum+=nw[*e].attr;
    }
    if(fabs(weight_sum-1)>1e-12) {
      log->info("Child weight sum is {}!", weight_sum);
    }

    if(out_degree(parent, nw)!=0) {
        weight_sum=0.;
        for(auto [e, e_end]=out_edges(parent, nw); e!=e_end; ++e){
          weight_sum+=nw[*e].attr;
        }
        if(fabs(weight_sum-1)>1e-12) {
          log->info("Parent weight sum is {}!", weight_sum);
        }
    }
}
//..............................................................................

template <typename NWType, typename LoggerType, typename RNGType, typename Config>
void ageing ( const double replacement_rate,
              int num_media,
              pair_int child_ages,
              pair_int parent_ages,
              pair_int senior_ages,
              NWType& nw,
              LoggerType& log,
              RNGType& rng,
              Config& cfg) {

      using vertex = typename boost::graph_traits<NWType>::vertex_descriptor;

      /* the containers in which the nodes we are rewiring and rewiring to are
      stored */
      std::vector<vertex> children, parents, peers;

      user_selection_and_ageing(children,
                                parents,
                                peers,
                                child_ages,
                                parent_ages,
                                senior_ages,
                                replacement_rate,
                                nw,
                                rng,
                                cfg);

      //check user ageing is possible in this step
      if(parents.size()==0) {
          log->info("There are no parent nodes: no user ageing possible in this step.");
          return;
      }
      if(children.size() > parents.size()){
          log->debug("Discrepancy between children and parent node numbers: \
have {} more children than parents.", children.size()-parents.size());
      }
      log->debug("Reinitialising {} vertices as children... ", children.size());
      log->debug("Available parents: {}", parents.size());
      log->debug("Available peers: {}", peers.size());

      /* since there are more peers than children, we move through the peer
      container at a different speed than through the child container. To ensure
      we do not rewire all children to the same peers, we need to remember where
      in the peers container we have got to */
      int at_peer = 0;

      //iterate over children container, and add and remove edges for each child
      for (int i=0; i<children.size(); ++i){

          /* loop over the children and rewire to parents. Since there may be
          fewer parents than children, some parents may get two or more children */
          const vertex child = children.at(i);
          const vertex parent = parents.at(i%parents.size());
          const int deg = degree(child, nw);
          const int in_deg = in_degree(child, nw);
          const int out_deg = out_degree(child, nw);

          reinitialize(child, parent, num_media, nw, rng, cfg);

          /*if a child has no social connections, cannot rewire, since we must
          preserve the edge count*/
          if(deg==0) {continue; }

          remove_edges(child, nw);

          add_edges(child, parent, peers, out_deg, in_deg, deg, at_peer, nw, log, rng);

          check_and_test(child, parent, deg, in_deg, out_deg, nw);

          /*if the degree is greater than two, we have added peers and need to
          move along in the peer container */
          if (deg>2) {at_peer+=deg-2;}

    } // rewiring

    log->debug("Ageing complete.");
}
} // namespace

#endif // UTOPIA_MODELS_OPDYN_REVISION
