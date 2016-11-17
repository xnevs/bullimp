/**
 * Create two graphs (the pattern and the target graph)
 * and use the UllImp algorithm to count the
 * induced isomorphisms pattern -> graph.
 */

#include <iostream>

#include <boost/graph/adjacency_list.hpp> 
#include <boost/graph/iteration_macros.hpp>

#include "../graph/ullimp_sub_graph_iso.hpp"

using graph_type =  boost::adjacency_list<boost::setS, boost::vecS, boost::bidirectionalS>; 

int main() { 
    using namespace boost;
    using namespace ullimp;

    graph_type pattern;
    auto u0 = add_vertex(pattern);
    auto u1 = add_vertex(pattern);
    auto u2 = add_vertex(pattern);
    add_edge(u0, u1, pattern);
    add_edge(u1, u2, pattern);
    add_edge(u2, u0, pattern);

    graph_type graph;
    auto v0 = add_vertex(graph);
    auto v1 = add_vertex(graph);
    auto v2 = add_vertex(graph);
    auto v3 = add_vertex(graph);
    add_edge(v0, v1, graph);
    add_edge(v1, v2, graph);
    add_edge(v2, v0, graph);
    add_edge(v1, v3, graph);
    add_edge(v3, v0, graph);

    auto index_map_p = get(vertex_index, pattern);
    auto index_map_g = get(vertex_index, graph);

    int count = 0;
    ullimp_subgraph_iso(pattern,
                        graph,
                        [&count, &pattern, &index_map_p, &index_map_g](auto const & iso, auto const & isoInv) {
                            std::cout << "found isomorphism:" << std::endl;
                            BGL_FORALL_VERTICES(u, pattern, graph_type) {
                                std::cout << "  " << index_map_p[u] << " -> " << index_map_g[iso[u]] << std::endl;
                            }
                            ++count;
                            return true;
                        });

    std::cout << "isomorphism count = " << count << std::endl;
}
