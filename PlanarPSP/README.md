## Introduction
This is the source code of Planar PSP indexes and non-partitioned SP indexes. Please refer to the paper for the algorithm details.

## Algorithms
After `cmake -DCMAKE_BUILD_TYPE=Release ./` and `make -j`, please run `./DPSP` and `./DSP` to obtain the usage instruction of planar PSP index and non-partitioned SP index, respectively.

**For the planar PSP index**,

1. arg 3 specifies the shortest path index: 1 for PCH, 2 for PTD; 3 for PPLL.
1. arg 4 specifies the PSP strategy: 1 for pre-boundary, 2 for no-boundary; 3 for post-boundary.
1. arg 5 specifies the partition method: NC for PUNCH; MT for METIS; SC for SCOTCH; kahypar for KaHyPar; geometric for RCB; Bubble for Bubble; HEP for HEP; CLUGP for CLUGP.
2. Note that before testing the performance of PSP index, you may need to specify arg 10 for preprocessing tasks, such as 1 for partition-based ordering, 2 for same-partition query generation.


Examples are shown as bellow.

* Test PTD with the pre-boundary strategy on NY with the PUNCH partition method and partition number 32:
`./DPSP ../Datasets NY 1 1 NC 32 1 100 150`
* Test N-CH-P on NY with PUNCH partition method and partition number 32:
`./DPSP ../Datasets NY 1 2 NC 32 1 100 150`
* Test P-TD-P on NY with PUNCH partition method and partition number 16:
`./DPSP ../Datasets NY 2 3 NC 16 1 100 150`

**For the non-partitioned SP index**,

1. arg 3 specifies the shortest path index: 1 for CH, 2 for TD (H2H); 3 for PPLL.
1. arg 5 specifies whether to use batch update. Although CH and TD support batch updates. You need to set this parameter to 0 to ensure fairness.
2. Note that you should use the same order file with partitioned counterparts by specifying the order file (arg 9) to ensure fair comparison. For example, use `../Datasets/NY/NY_NC_32/vertex_orderMDE2`.

Examples are shown as bellow.

* Test CH on NY with the corresponding PSP index using the PUNCH partition method and partition number 32:
`./DSP ../Datasets NY 1 1 0 100 1 150 ../Datasets/NY/NY_NC_32/vertex_orderMDE2`
* Test H2H on NY with the corresponding PSP index using the PUNCH partition method and partition number 16:
`./DSP ../Datasets NY 2 1 0 100 1 150 ../Datasets/NY/NY_NC_16/vertex_orderMDE2`


## Data
The datasets of this paper source from [http://www.dis.uniroma1.it/challenge9/download.shtml](http://www.dis.uniroma1.it/challenge9/download.shtml), [http://snap.stanford.edu/data](http://snap.stanford.edu/data), and [http://konect.cc/](http://konect.cc/). 
Additionally, the processed data, including partitioned data, is available in our [OneDrive Link](https://hkustgz-my.sharepoint.com/:f:/g/personal/xzhouby_connect_hkust-gz_edu_cn/EkEOQqUbSMZKioVFPdUvJisBSvhvzn0dR-ubJtpt7pmX5A?e=UWolbO). Please refer to the paper for details.

You can also use the example dataset in the directory `../Datasets` for testing.


## Dependency

1. `g++` 
2. `boost`
3. `omp`

All the codes are runnable after `cmake` and `make`: go to the corresponding directory, `cmake -DCMAKE_BUILD_TYPE=Release ./` and `make -j`.
