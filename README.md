## Introduction
This is the source code of the ICDE 2024 paper "*A Universal Scheme for Partitioned Dynamic Shortest Path Index*" (submitted). Please refer to the paper for the algorithm details.

## Algorithms

The following directory contains the implementation code of index construction, query processing, and index update of corresponding algorithms.

1. P-TD: P-TD and P-TD(Post), where P-TD(Post) is equal to Q-Planar while P-TD is equal to FHL.
1. P-CH: P-CH, which is equal to U-Planar.
1. CT-TD: which is equal to Q-Core.
1. CT-CH: which is equal to U-Core.
1. HP-TD: which is equal to UQ-Hierarchy.

P-TD, P-CH and HP-TD can be directly run on the example graph *NY* while CT-TD, CT-CH can be run on *GO*, by using source path `../data`.

## Data
The datasets of this paper source from [http://www.dis.uniroma1.it/challenge9/download.shtml](http://www.dis.uniroma1.it/challenge9/download.shtml), [http://snap.stanford.edu/data](http://snap.stanford.edu/data), and [http://konect.cc/](http://konect.cc/). Please refer to the paper for details.

An example graph *NY* and corresponding partition results of `PUNCH` with 64 partitions is provided in the directory *data* for your reference.

1. To implement P-TD and P-CH, you need to generate the partition results (`subgraph_edge`, `subgraph_vertex`, `cut_edges`) and corresponding vertex order file `vertex_order`.
2. Query OD pair `NY.query` and update OD pair `NY.update` are also needed.


## Dependency

1. `g++` and `boost`
1. `METIS` (for HP-TD).

All the codes are runnable after `cmake` and `make`: go to the corresponding directory, `cmake -DCMAKE_BUILD_TYPE=Release ./` and `make -j`.
