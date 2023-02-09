## Introduction
This is the source code of the SIGMOD 2024 paper "*A Universal Scheme for Partitioned Dynamic Shortest Path Index*" (submitted). Please refer to the paper for the algorithm details.

## Algorithms

The following directory contains the implementation code of index construction, query processing, and index update of corresponding algorithms.

1. P-TD: P-TD and P-TD(Post).
1. P-CH: P-CH and P-CH(Post).
1. P-PLL: P-PLL and P-PLL(Post).
1. CT: CT-DS and CT-TD.
1. HP-DS: HP-DS.

All above methods can be directly run on the example graph *NY*, by using source path `../data`.

## Data
An example graph *NY* and corresponding partition results of `PUNCH` with 64 partitions is provided in directory *data* for your reference.

1. To implement P-TD, P-CH, and P-PLL, you need to generate the partition results (`subgraph_edge`, `subgraph_vertex`, `cut_edges`) and corresponding vertex order file `vertex_order`.
1. To implement CT-DS, CT-TD, you need to generate the vertex order of the whole graph, i.e., `NY.order`.
2. Query OD pair `NY.query` and update OD pair `NY.update` are also needed.


## Dependency

1. `g++` and `boost`
1. `METIS` (for HP-DS).

All the codes are runnable after cmake and make: go to corresponding directory, `cmake -DCMAKE_BUILD_TYPE=Release ./` and `make -j`.
