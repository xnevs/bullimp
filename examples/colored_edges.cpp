#include <iostream>
#include <string>

#include <boost/graph/adjacency_list.hpp> 
#include <boost/graph/iteration_macros.hpp>

#include "../graph/ullimp_sub_graph_iso.hpp"
#include "../graph/ullimp_order.hpp"

struct Color {
    std::string color;
};

using graph_type =  boost::adjacency_list<boost::setS, boost::vecS, boost::bidirectionalS, boost::no_property, Color>; 

int main() { 
    using namespace boost;
    using namespace ullimp;

    graph_type pattern;
    auto u0 = add_vertex(pattern);
    auto u1 = add_vertex(pattern);
    auto u2 = add_vertex(pattern);

    auto e01 = add_edge(u0, u1, pattern).first;
    auto e12 = add_edge(u1, u2, pattern).first;
    auto e20 = add_edge(u2, u0, pattern).first;

    pattern[e01].color = "red";
    pattern[e12].color = "black";
    pattern[e20].color = "black";

    graph_type graph;
    auto v0 = add_vertex(graph);
    auto v1 = add_vertex(graph);
    auto v2 = add_vertex(graph);
    auto v3 = add_vertex(graph);

    auto f01 = add_edge(v0, v1, graph).first;
    auto f12 = add_edge(v1, v2, graph).first;
    auto f20 = add_edge(v2, v0, graph).first;
    auto f13 = add_edge(v1, v3, graph).first;
    auto f30 = add_edge(v3, v0, graph).first;

    graph[f01].color = "red";
    graph[f12].color = "black";
    graph[f20].color = "black";
    graph[f13].color = "black";
    graph[f30].color = "black";

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
                        },
                        get(vertex_index, pattern),
                        get(vertex_index, graph),
                        vertex_order_RDEG_CNC(pattern),
                        [&pattern, &graph](auto const & e, auto const & f) {
                            return pattern[e].color == graph[f].color;
                        },
                        always_equivalent());

    std::cout << "isomorphism count = " << count << std::endl;
}
