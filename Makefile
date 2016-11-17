CXX=g++
CXX_FLAGS=-std=c++14 -O3 -march=native -mtune=native

all: simple colored_vertices colored_edges compare

simple: examples/simple.cpp graph/ullimp_sub_graph_iso.hpp
	$(CXX) $(CXX_FLAGS) -o simple examples/simple.cpp

colored_vertices: examples/colored_vertices.cpp graph/ullimp_sub_graph_iso.hpp graph/ullimp_order.hpp
	$(CXX) $(CXX_FLAGS) -o colored_vertices examples/colored_vertices.cpp

colored_edges: examples/colored_edges.cpp graph/ullimp_sub_graph_iso.hpp graph/ullimp_order.hpp
	$(CXX) $(CXX_FLAGS) -o colored_edges examples/colored_edges.cpp

compare: examples/compare.cpp graph/ullimp_sub_graph_iso.hpp
	$(CXX) $(CXX_FLAGS) -o compare examples/compare.cpp

@PHONY: clean
clean:
	rm simple
	rm colored_vertices
	rm colored_edges
	rm compare
