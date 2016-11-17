# UllImp

An implementation of an improved Ullman's algorithm for the subgraph isomorphism problem.

## Dependencies

The implementation uses the [Boost Graph Library](http://www.boost.org/doc/libs/1_62_0/libs/graph/doc/)

## Examples

Examples of use are located in the `examples` folder.
You can build them by running `make` in the root directory.

The `compare` example is intended to compare this algorithm with the one present in the Boost Graph Library (VF2).
The inputs to the `compare` executable are intended to be graphs from the [VFLib graph database](http://mivia.unisa.it/datasets/graph-database/vflib/).
