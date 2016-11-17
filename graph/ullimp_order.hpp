#ifndef ULLIMP_ORDER_HPP
#define ULLIMP_ORDER_HPP

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/iteration_macros.hpp>

#include "inv_adj_iteration_macro.hpp"

#include <vector>

#include <iostream>

namespace ullimp {

namespace detail {

template <typename Graph>
int clustering1(typename boost::graph_traits<Graph>::vertex_descriptor u, Graph const & g) {
    using namespace boost;

    int c = 0;
    BGL_FORALL_ADJ_T(u, w, g, Graph) {
        BGL_FORALL_ADJ_T(w, r, g, Graph) {
            if(edge(r, u, g).second || edge(u, r, g).second) 
                ++c;
        }
        BGL_FORALL_INV_ADJ_T(w, r, g, Graph) {
            if(edge(r, u, g).second || edge(u, r, g).second) 
                ++c;
        }
    }
    BGL_FORALL_INV_ADJ_T(u, w, g, Graph) {
        BGL_FORALL_ADJ_T(w, r, g, Graph) {
            if(edge(r, u, g).second || edge(u, r, g).second) 
                ++c;
        }
        BGL_FORALL_INV_ADJ_T(w, r, g, Graph) {
            if(edge(r, u, g).second || edge(u, r, g).second) 
                ++c;
        }
    }
    return c;
}

} // namespace detail

template <typename Graph, typename IndexMap>
std::vector<typename boost::graph_traits<Graph>::vertex_descriptor>
vertex_order_RDEG_CNC(Graph const &g, IndexMap index_map) {
    using namespace boost;

    using vertex_type = typename graph_traits<Graph>::vertex_descriptor;

    auto n = num_vertices(g);

    std::vector<vertex_type> vertex_order(n);
    
    std::vector<int> clustdeg_vec(n);
    auto clustdeg = make_iterator_property_map(std::begin(clustdeg_vec), index_map);

    BGL_FORALL_VERTICES_T(u, g, Graph) {
        clustdeg[u] = detail::clustering1(u, g) + degree(u, g);
    }

    std::vector<bool> avail_vec(n, true);
    auto avail = make_iterator_property_map(std::begin(avail_vec), index_map);

    for(int idx=0; idx<n; ++idx) {
        auto bestn = graph_traits<Graph>::null_vertex();
        int bestv = -1;
        BGL_FORALL_VERTICES_T(i, g, Graph) {
            if(avail[i]) {
                int deg = 0;
                BGL_FORALL_ADJ_T(i, p, g, Graph) {
                    if(!avail[p])
                        deg++;
                }
                BGL_FORALL_INV_ADJ_T(i, p, g, Graph) {
                    if(!avail[p])
                        deg++;
                }
                if(bestn == graph_traits<Graph>::null_vertex() || (deg > bestv || (deg == bestv && clustdeg[i] > clustdeg[bestn]))) {
                    bestn = i;
                    bestv = deg;
                }
            }
        }
        avail[bestn] = false;
        vertex_order[idx] = bestn;
    }

    return vertex_order;
}

template <typename Graph>
std::vector<typename boost::graph_traits<Graph>::vertex_descriptor>
vertex_order_RDEG_CNC(Graph const &g) {
    return vertex_order_RDEG_CNC(g, boost::get(boost::vertex_index, g));
}

} // namespace ullimp

#endif // ULLIMP_ORDER_HPP
