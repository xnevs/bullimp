#ifndef ULLIMP_SUB_GRAPH_ISO_HPP
#define ULLIMP_SUB_GRAPH_ISO_HPP

#include <boost/graph/graph_utility.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/iteration_macros.hpp>

#include "inv_adj_iteration_macro.hpp"

#include <boost/graph/mcgregor_common_subgraphs.hpp> // for always_equivalent

#include "ullimp_order.hpp"

namespace ullimp {

namespace detail {
    using namespace boost;

    template <typename Graph,
              typename IndexMap>
    class Digraph {
        using vertex_type = typename graph_traits<Graph>::vertex_descriptor;
        using vertices_size_type = typename graph_traits<Graph>::vertices_size_type;

        vertices_size_type nodeCnt;

        IndexMap index_map;

        std::vector<char> adjacency_vec;

        using degs_type = iterator_property_map<typename std::vector<vertices_size_type>::iterator,
                                                IndexMap, vertex_type, vertices_size_type &>;
        std::vector<vertices_size_type> indegs_vec;
        std::vector<vertices_size_type> outdegs_vec;

        degs_type indegs;
        degs_type outdegs;

        using neighs_type = iterator_property_map<typename std::vector<std::vector<vertex_type>>::iterator,
                                                  IndexMap, vertex_type, std::vector<vertex_type> &>;
        std::vector<std::vector<vertex_type>> inneighs_vec;
        std::vector<std::vector<vertex_type>> outneighs_vec;
        std::vector<std::vector<vertex_type>> notinneighs_vec;
        std::vector<std::vector<vertex_type>> notoutneighs_vec;

        neighs_type inneighs;
        neighs_type outneighs;
        neighs_type notinneighs;
        neighs_type notoutneighs;

        char & adjacency(vertex_type u, vertex_type v) {
            return adjacency_vec[index_map[u]*nodeCnt + index_map[v]];
        }

      public:
        Digraph(Graph const & g, IndexMap index_map)
                    : nodeCnt(num_vertices(g)),
                      index_map(index_map),
                      adjacency_vec(nodeCnt*nodeCnt, false),
                      indegs_vec(num_vertices(g)), indegs(std::begin(indegs_vec), index_map),
                      outdegs_vec(num_vertices(g)), outdegs(std::begin(outdegs_vec), index_map),
                      inneighs_vec(num_vertices(g)), inneighs(std::begin(inneighs_vec), index_map),
                      outneighs_vec(num_vertices(g)), outneighs(std::begin(outneighs_vec), index_map),
                      notinneighs_vec(num_vertices(g)), notinneighs(std::begin(notinneighs_vec), index_map),
                      notoutneighs_vec(num_vertices(g)), notoutneighs(std::begin(notoutneighs_vec), index_map) {
            BGL_FORALL_VERTICES_T(u, g, Graph) {
                BGL_FORALL_ADJ_T(u, v, g, Graph) {
                    adjacency(u,v) = true;
                    ++outdegs[u];
                    ++indegs[v];
                }
            }
            BGL_FORALL_VERTICES_T(u, g, Graph) {
                inneighs[u].reserve(indegs[u]);
                outneighs[u].reserve(outdegs[u]);
                notinneighs[u].reserve(nodeCnt - indegs[u] - 1);
                notoutneighs[u].reserve(nodeCnt - outdegs[u] - 1);
                BGL_FORALL_VERTICES_T(v, g, Graph) {
                    if(v != u) {
                        if(adjacency(v,u)) {
                            inneighs[u].push_back(v);
                        } else {
                            notinneighs[u].push_back(v);
                        }
                        if(adjacency(u,v)) {
                            outneighs[u].push_back(v);
                        } else {
                            notoutneighs[u].push_back(v);
                        }
                    }
                }
            }
        }

        bool edge(vertex_type u, vertex_type v) {
            return adjacency(u,v);
        }

        std::vector<vertex_type> const & inneigh(vertex_type u) {
            return inneighs[u];
        }
        std::vector<vertex_type> const & outneigh(vertex_type u) {
            return outneighs[u];
        }
        std::vector<vertex_type> const & notinneigh(vertex_type u) {
            return notinneighs[u];
        }
        std::vector<vertex_type> const & notoutneigh(vertex_type u) {
            return notoutneighs[u];
        }
    };

    template <typename GraphSmall,
              typename GraphLarge,
              typename IndexMapSmall,
              typename IndexMapLarge,
              typename VertexOrderSmall,
              typename EdgeEquivalencePredicate,
              typename VertexEquivalencePredicate,
              typename SubGraphIsoMapCallback>
    class UllImp {
        using vertex_small_type = typename graph_traits<GraphSmall>::vertex_descriptor;
        using vertex_large_type = typename graph_traits<GraphLarge>::vertex_descriptor;

        using vertex_small_iter_type = typename graph_traits<GraphSmall>::vertex_iterator;
        using vertex_large_iter_type = typename graph_traits<GraphLarge>::vertex_iterator;

        using vertices_small_size_type = typename graph_traits<GraphSmall>::vertices_size_type;
        using vertices_large_size_type = typename graph_traits<GraphLarge>::vertices_size_type;

        class Map {
            vertices_small_size_type rows;        // number of rows
            vertices_large_size_type cols;        // number of columns

            IndexMapSmall index_map_small;
            IndexMapLarge index_map_large;

            bool* data;                 // current matrix
            
            // history - stack of all changes
            int* history;               // stack of indices of changes
            int index;                  // top of stack
            int lastindex;              // last index
            // shots - snapshosts in history
            int* shots;                 // stack of snapshots
            int shotidx;                // top of stack

          public:
            Map(vertices_small_size_type rows, vertices_large_size_type cols,
                IndexMapSmall const & index_map_small, IndexMapLarge const & index_map_large)
                : rows(rows), cols(cols),
                  index_map_small(index_map_small), index_map_large(index_map_large) {

                auto size = rows * cols;
                data = new bool[size];
                history = new int[size];

                index = -1;
                lastindex = -1;
                shots = new int[size];
                shotidx = -1;
            }
            
            ~Map() {
                delete[] data;
                delete[] history;
                delete[] shots;
            }
            
            bool & operator()(typename graph_traits<GraphSmall>::vertex_descriptor row, typename graph_traits<GraphLarge>::vertex_descriptor col) {
                return data[index_map_small[row] * cols + index_map_large[col]];
            }

            bool operator()(typename graph_traits<GraphSmall>::vertex_descriptor row, typename graph_traits<GraphLarge>::vertex_descriptor col) const {
                return data[index_map_small[row] * cols + index_map_large[col]];
            }
            
            void snapshot() {
                shots[++shotidx] = index;
                reset();
            }
            
            void undo() {
                int stop = shots[shotidx--];
                while (index > stop)
                    data[history[index--]] = true;
            }
            
            void annul(vertex_small_type row, vertex_large_type col) {
                auto idx = index_map_small[row] * cols + index_map_large[col];
                if (data[idx]) {
                    data[idx] = false;
                    history[++index] = idx;
                }
            }
            
            void reset() {
                lastindex = index;
            }
        };

        GraphSmall const & g;         // pattern graph
        GraphLarge const & h;         // target graph

        Digraph<GraphSmall,IndexMapSmall> gg;
        Digraph<GraphLarge,IndexMapLarge> hh;

        IndexMapSmall index_map_small;
        IndexMapLarge index_map_large;

        vertices_small_size_type n;   // n=size of pattern
        vertices_large_size_type m;   // m=size of target

        VertexOrderSmall vertex_order_small;

        EdgeEquivalencePredicate edge_comp;
        VertexEquivalencePredicate vertex_comp;

        SubGraphIsoMapCallback user_callback;

        std::vector<vertex_large_type> iso_vec;
        using iso_type = iterator_property_map<typename std::vector<vertex_large_type>::iterator,
                                               IndexMapSmall, vertex_large_type, vertex_large_type &>;
        iso_type iso;       // isomorphism: iso[i] -> j

        std::vector<vertex_small_type> isoinv_vec;
        using isoinv_type = iterator_property_map<typename std::vector<vertex_small_type>::iterator,
                                                  IndexMapLarge, vertex_small_type, vertex_small_type &>;
        isoinv_type isoinv; // isomorphism inverse

        Map map;       // possible mappings (matrix with snapshots)

        std::vector<std::vector<vertex_small_type>> succInNeigh_vec;
        std::vector<std::vector<vertex_small_type>> succOutNeigh_vec;
        std::vector<std::vector<vertex_small_type>> succNotInNeigh_vec;
        std::vector<std::vector<vertex_small_type>> succNotOutNeigh_vec;
        using succNeigh_type = iterator_property_map<typename std::vector<std::vector<vertex_small_type>>::iterator,
                                                     IndexMapSmall, vertex_small_type, std::vector<vertex_small_type> &>;
        succNeigh_type succInNeigh;
        succNeigh_type succOutNeigh;
        succNeigh_type succNotInNeigh;
        succNeigh_type succNotOutNeigh;

        void init() {
            BGL_FORALL_VERTICES_T(i, g, GraphSmall) {
                auto a = in_degree(i, g);
                auto b = out_degree(i, g);
                BGL_FORALL_VERTICES_T(j, h, GraphLarge) {
                    map(i, j) = in_degree(j, h) >= a && out_degree(j, h) >= b;
                }
            }

            refineMapFull();
        }
        
        bool canMapFull(vertex_small_type row, vertex_large_type col) {
            BGL_FORALL_ADJ_T(row, u, g, GraphSmall) {
                bool all0 = true;
                BGL_FORALL_ADJ_T(col, v, h, GraphLarge) {
                    if(map(u, v)) {
                        all0 = false;
                        break;
                    }
                }
                if(all0)
                    return false;
            }
            BGL_FORALL_INV_ADJ_T(row, u, g, GraphSmall) {
                bool all0 = true;
                BGL_FORALL_INV_ADJ_T(col, v, h, GraphLarge) {
                    if(map(u, v)) {
                        all0 = false;
                        break;
                    }
                }
                if(all0)
                    return false;
            }
            return true;
        }
        
        void refineMapFull() {
            bool changed;
            do {
                changed = false;
                BGL_FORALL_VERTICES_T(i, g, GraphSmall) {
                    BGL_FORALL_VERTICES_T(j, h, GraphLarge) {
                        if(map(i, j) && !canMapFull(i, j)) {
                            map(i, j) = false;
                            changed = true;
                        }
                    }
                }
            } while(changed);
        }
        
        bool canMap(vertex_small_type row, vertex_large_type col) {
            for(auto u : succInNeigh[row]) {
                bool all0 = true;
                for(auto v : hh.inneigh(col)) {
                    if(isoinv[v] == graph_traits<GraphSmall>::null_vertex() && map(u, v)) {
                        all0 = false;
                        break;
                    }
                }
                if(all0) {
                    return false;
                }
            }
            for(auto u : succOutNeigh[row]) {
                bool all0 = true;
                for(auto v : hh.outneigh(col)) {
                    if(isoinv[v] == graph_traits<GraphSmall>::null_vertex() && map(u, v)) {
                        all0 = false;
                        break;
                    }
                }
                if(all0) {
                    return false;
                }
            }
            return true;
        }
        
        void refineMapPartial(vertex_small_type row, vertex_large_type col) {
            for(auto u : succInNeigh[row]) {
                for(auto v : hh.inneigh(col)) {
                    if(isoinv[v] == graph_traits<GraphSmall>::null_vertex() && map(u, v) && !canMap(u, v)) {
                        map.annul(u, v);
                    }
                }
            }
            for(auto u : succOutNeigh[row]) {
                for(auto v : hh.outneigh(col)) {
                    if(isoinv[v] == graph_traits<GraphSmall>::null_vertex() && map(u, v) && !canMap(u, v)) {
                        map.annul(u, v);
                    }
                }
            }
/*
            BGL_FORALL_INV_ADJ_T(row, u, g, GraphSmall) {
                BGL_FORALL_INV_ADJ_T(col, v, h, GraphLarge) {
                    if(isoinv[v] == graph_traits<GraphSmall>::null_vertex() && map(u, v) && !canMap(u, v)) {
                        map.annul(u, v);
                    }
                }
            }
            BGL_FORALL_ADJ_T(row, u, g, GraphSmall) {
                BGL_FORALL_ADJ_T(col, v, h, GraphLarge) {
                    if(isoinv[v] == graph_traits<GraphSmall>::null_vertex() && map(u, v) && !canMap(u, v)) {
                        map.annul(u, v);
                    }
                }
            }
*/
        }
        
        void filter(vertex_small_type row, vertex_large_type col) {
            for(auto u : succInNeigh[row]) {
                for(auto v : hh.notinneigh(col)) {
                    if(isoinv[v] == graph_traits<GraphSmall>::null_vertex()) {
                        map.annul(u, v);
                    }
                }
            }
            for(auto u : succOutNeigh[row]) {
                for(auto v : hh.notoutneigh(col)) {
                    if(isoinv[v] == graph_traits<GraphSmall>::null_vertex()) {
                        map.annul(u, v);
                    }
                }
            }
            for(auto v : hh.inneigh(col)) {
                for(auto u : succNotInNeigh[row]) {
                    map.annul(u, v);
                }
            }
            for(auto v : hh.outneigh(col)) {
                for(auto u : succNotOutNeigh[row]) {
                    map.annul(u, v);
                }
            }
/*
            BGL_FORALL_VERTICES_T(u, g, GraphSmall) {
                if(u != row) {
                    bool row_u = edge(row, u, g).second;
                    bool u_row = edge(u, row, g).second;
                    BGL_FORALL_VERTICES_T(v, h, GraphLarge) {
                        if(v != col && isoinv[v] == graph_traits<GraphSmall>::null_vertex()) {
                            bool col_v = edge(col, v, h).second;
                            bool v_col = edge(v, col, h).second;
                            if( row_u && !col_v) map.annul(u, v);
                            if( u_row && !v_col) map.annul(u, v);
                            if(!row_u &&  col_v) map.annul(u, v);
                            if(!u_row &&  v_col) map.annul(u, v);
                        }
                    }
                }
            }
*/
        }
        
        vertex_large_type nextCol(vertex_small_type row) { //TODO glej nextCol1 v mici ullimp
            if(iso[row] != graph_traits<GraphLarge>::null_vertex())
                isoinv[iso[row]] = graph_traits<GraphSmall>::null_vertex();
            
            /*
            auto col = vertex(index_map_large[iso[row]]+1, h);
            while (col < m &&    // TODO pazi col je vertex descriptor
                    (!map(row, col) ||
                     isoinv[col] != graph_traits<GraphSmall>::null_vertex() ||
                     !vertex_comp(row, col) ||
                     !checkEdges(row, col)))
                col += 1;
            return (col < m ? col : graph_traits<GraphLarge>::null_vertex());
            */
            
            vertex_small_iter_type col_it, last_it;
            std::tie(col_it, last_it) = vertices(h);
            std::advance(col_it, index_map_large[iso[row]]+1);
            while(col_it != last_it &&
                   (!map(row, *col_it) ||
                    isoinv[*col_it] != graph_traits<GraphSmall>::null_vertex() ||
                    !vertex_comp(row, *col_it) ||
                    !checkEdges(row, *col_it)))
                ++col_it;

            return (col_it != last_it ? *col_it : graph_traits<GraphLarge>::null_vertex());
        }
        
        bool run() {
            bool iso_found = false;
            decltype(n) depth = 0;
            auto row_iter = std::begin(vertex_order_small);
            while(true) {
                while(depth < n - 1) {
                    auto row = *row_iter;
                    auto col = nextCol(row);
                    if(col != graph_traits<GraphLarge>::null_vertex()) {
                        iso[row] = col;									// map node -> col
                        isoinv[col] = row;								// mark column as used
                        map.snapshot();
                        filter(row, col);
                        if (depth <= n/2) refineMapPartial(row, col);
                        if (canMap(row, col)) {
                            ++depth;
                            ++row_iter;
                        } else {
                            map.undo();
                        }
                    } else {
                        if(row != graph_traits<GraphSmall>::null_vertex()) {   // backtrack iso
                            iso[row] = graph_traits<GraphLarge>::null_vertex();
                        }			

                        if(depth-- == 0) // all paths tried?
                            return iso_found;

                        --row_iter;
                        map.undo();
                    }
                }

                // all rows mapped (depth == num_vertices(g) - 1)
                
                auto row = *row_iter;

                BGL_FORALL_VERTICES_T(col, h, GraphLarge) {
                    if(isoinv[col] == graph_traits<GraphSmall>::null_vertex() &&
                            map(row, col) && vertex_comp(row, col) && checkEdges(row, col)) {
                        iso[row] = col;
                        isoinv[col] = row;
                        if(!user_callback(iso, isoinv))
                            return true;
                        iso_found = true;
                        iso[row] = graph_traits<GraphLarge>::null_vertex();
                        isoinv[col] = graph_traits<GraphSmall>::null_vertex();
                    }
                }

                map.undo();

                --depth;
                --row_iter;
            }

            return iso_found;
        }

        bool checkEdges(vertex_small_type row, vertex_large_type col) {
            BGL_FORALL_OUTEDGES_T(row, e, g, GraphSmall) {
                auto u = target(e, g);
                auto v = iso[u];
                if(v != graph_traits<GraphLarge>::null_vertex()) {
                    auto pe = edge(col, v, h);
                    if(!pe.second || !edge_comp(e, pe.first))
                        return false;
                }
            }
            BGL_FORALL_INEDGES_T(row, e, g, GraphSmall) {
                auto u = source(e, g);
                auto v = iso[u];
                if(v != graph_traits<GraphLarge>::null_vertex()) {
                    auto pe = edge(v, col, h);
                    if(!pe.second || !edge_comp(e, pe.first))
                        return false;
                }
            }
            return true;
        }

    public:
        UllImp(GraphSmall const & g, GraphLarge const & h, SubGraphIsoMapCallback user_callback,
               IndexMapSmall const & index_map_small, IndexMapLarge const & index_map_large,
               VertexOrderSmall const & vertex_order_small,
               EdgeEquivalencePredicate edge_comp, VertexEquivalencePredicate vertex_comp)
                    : g(g), h(h),
                      gg(g, index_map_small), hh(h, index_map_large),
                      index_map_small(index_map_small), index_map_large(index_map_large),
                      n(num_vertices(g)), m(num_vertices(h)),
                      vertex_order_small(vertex_order_small),
                      edge_comp(edge_comp), vertex_comp(vertex_comp),
                      map(n,m, index_map_small, index_map_large),
                      user_callback(user_callback),
                      iso_vec(n, graph_traits<GraphLarge>::null_vertex()),
                      iso(std::begin(iso_vec), index_map_small),
                      isoinv_vec(m, graph_traits<GraphSmall>::null_vertex()),
                      isoinv(std::begin(isoinv_vec), index_map_large),

                      succInNeigh_vec(n), succInNeigh(std::begin(succInNeigh_vec), index_map_small),
                      succOutNeigh_vec(n), succOutNeigh(std::begin(succOutNeigh_vec), index_map_small),
                      succNotInNeigh_vec(n), succNotInNeigh(std::begin(succNotInNeigh_vec), index_map_small),
                      succNotOutNeigh_vec(n), succNotOutNeigh(std::begin(succNotOutNeigh_vec), index_map_small) {

            std::vector<vertices_small_size_type> orderinv_vec(n);
            using orderinv_type = iterator_property_map<typename std::vector<vertices_small_size_type>::iterator,
                                                        IndexMapSmall, vertex_small_type, vertices_small_size_type &>;
            orderinv_type orderinv(std::begin(orderinv_vec), index_map_small);

            for(int i=0; i<n; ++i) {
                orderinv[vertex_order_small[i]] = i;
            }
            BGL_FORALL_VERTICES_T(u, g, GraphSmall) {
                vertices_small_size_type ordu = orderinv[u];
                for(auto v : gg.inneigh(u)) {
                    if(orderinv[v] > ordu) succInNeigh[u].push_back(v);
                }
                for(auto v : gg.outneigh(u)) {
                    if(orderinv[v] > ordu) succOutNeigh[u].push_back(v);
                }
                for(auto v : gg.notinneigh(u)) {
                    if(orderinv[v] > ordu) succNotInNeigh[u].push_back(v);
                }
                for(auto v : gg.notoutneigh(u)) {
                    if(orderinv[v] > ordu) succNotOutNeigh[u].push_back(v);
                }
            }
        }
        
        bool find() {
            init();
            return run();
        }
    };
} // namespace detail

template<typename GraphSmall,
         typename GraphLarge,
         typename IndexMapSmall,
         typename IndexMapLarge,
         typename VertexOrderSmall,
         typename EdgeEquivalencePredicate,
         typename VertexEquivalencePredicate,
         typename SubGraphIsoMapCallback>
bool ullimp_subgraph_iso(GraphSmall const & graph_small, GraphLarge const & graph_large,
                         SubGraphIsoMapCallback user_callback,
                         IndexMapSmall const & index_map_small, IndexMapLarge const & index_map_large, 
                         VertexOrderSmall const &vertex_order_small,
                         EdgeEquivalencePredicate edge_comp,
                         VertexEquivalencePredicate vertex_comp) {
    detail::UllImp<GraphSmall, GraphLarge, IndexMapSmall, IndexMapLarge, VertexOrderSmall,
                   EdgeEquivalencePredicate, VertexEquivalencePredicate, SubGraphIsoMapCallback> 
        ullimp{graph_small, graph_large, user_callback,
               index_map_small, index_map_large,
               vertex_order_small, edge_comp, vertex_comp};
    return ullimp.find();
}

// All default interface for ullimp_subgraph_iso
template <typename GraphSmall,
          typename GraphLarge,
          typename SubGraphIsoMapCallback>
bool ullimp_subgraph_iso(const GraphSmall & graph_small, const GraphLarge & graph_large, 
                         SubGraphIsoMapCallback user_callback) {
    return ullimp_subgraph_iso(graph_small, graph_large, user_callback,
                               boost::get(boost::vertex_index, graph_small), boost::get(boost::vertex_index, graph_large),
                               vertex_order_RDEG_CNC(graph_small),
                               boost::always_equivalent(), boost::always_equivalent());
}

} // namespace ullimp

#endif // ULLIMP_SUB_GRAPH_ISO_HPP
