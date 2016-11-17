#ifndef INV_ADJ_ITERATION_MACRO_HPP
#define INV_ADJ_ITERATION_MACRO_HPP

#include <boost/graph/iteration_macros.hpp>

#define BGL_FORALL_INV_ADJ_T(UNAME, VNAME, GNAME, GraphType) \
for (std::pair<typename GraphType::inv_adjacency_iterator, \
               typename GraphType::inv_adjacency_iterator> BGL_RANGE(__LINE__) = inv_adjacent_vertices(UNAME, GNAME); \
  BGL_FIRST(__LINE__) != BGL_LAST(__LINE__); BGL_FIRST(__LINE__) = BGL_LAST(__LINE__)) \
for (typename boost::graph_traits<GraphType>::vertex_descriptor VNAME; \
  BGL_FIRST(__LINE__) != BGL_LAST(__LINE__) ? (VNAME = *BGL_FIRST(__LINE__), true) : false; \
   ++BGL_FIRST(__LINE__))

#define BGL_FORALL_INV_ADJ(UNAME, VNAME, GNAME, GraphType) \
for (std::pair<GraphType::inv_adjacency_iterator, \
               GraphType::inv_adjacency_iterator> BGL_RANGE(__LINE__) = inv_adjacent_vertices(UNAME, GNAME); \
  BGL_FIRST(__LINE__) != BGL_LAST(__LINE__); BGL_FIRST(__LINE__) = BGL_LAST(__LINE__)) \
for (boost::graph_traits<GraphType>::vertex_descriptor VNAME; \
  BGL_FIRST(__LINE__) != BGL_LAST(__LINE__) ? (VNAME = *BGL_FIRST(__LINE__), true) : false; \
   ++BGL_FIRST(__LINE__))


#endif // INV_ADJ_ITERATION_MACRO_HPP
