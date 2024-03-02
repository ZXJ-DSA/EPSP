## Introduction
This is the source code of the VLDB 2024 paper "*A Universal Scheme for Partitioned Dynamic Shortest Path Index*" (submitted). Please refer to the paper for the algorithm details.

## Algorithms

The following directory contains the implementation code of index construction, query processing, and index update of corresponding algorithms.

1. PlanarPSP: PCH, PTD, PPLL, P-TD-P, and N-CH-P.
1. CorePeripheryPSP->CoreTreePSP: CT, P-PT-CP, and N-PC-CP.
1. CorePeripheryPSP->SketchPSP: Sketch.
1. HierarchicalPSP: G-tree and N-TS-HP.

## Data
The datasets of this paper source from [http://www.dis.uniroma1.it/challenge9/download.shtml](http://www.dis.uniroma1.it/challenge9/download.shtml), [http://snap.stanford.edu/data](http://snap.stanford.edu/data), and [http://konect.cc/](http://konect.cc/). The processed data (including partitioned data) is available in [OneDrive Link](https://hkustgz-my.sharepoint.com/:f:/g/personal/xzhouby_connect_hkust-gz_edu_cn/EkEOQqUbSMZKioVFPdUvJisBSvhvzn0dR-ubJtpt7pmX5A?e=UWolbO). Please refer to the paper for details.


## Dependency

1. `g++` 
2. `boost`
3. `omp`
1. `METIS` (for Hierarchical PSP).

All the codes are runnable after `cmake` and `make`: go to the corresponding directory, `cmake -DCMAKE_BUILD_TYPE=Release ./` and `make -j`.
