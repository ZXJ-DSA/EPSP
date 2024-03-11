## Introduction
This is the source code of Core-tree PSP indexes. Please refer to the paper for the algorithm details.

## Algorithms
After `cmake -DCMAKE_BUILD_TYPE=Release ./` and `make -j`, please run `./DCT` to obtain the usage instruction.

**For the core-tree PSP index**,

1. arg 4 specifies the shortest path index for the tree index: 0 for CH, 1 for TD.
1. arg 8 specifies the PSP strategy: 1 for pre-boundary, 2 for no-boundary; 3 for post-boundary.
2. Note that before testing the performance of PSP index, you may need to specify arg 9 for preprocessing tasks, such as 1 for same-tree query generation.


Examples are shown as bellow.

* Test P-PT-CP with bandwidth 20 on NY:
`./DCT ../../Datasets NY 20 1 1 100 150 3`
* Test N-PC-CP with bandwidth 40 on NY:
`./DCT ../../Datasets NY 20 0 1 100 150 2`
* Test CT with bandwidth 40 on GOOGLE:
`./DCT ../../Datasets GOOGLE 40 1 1 100 150 2`



## Data
The datasets of this paper source from [http://www.dis.uniroma1.it/challenge9/download.shtml](http://www.dis.uniroma1.it/challenge9/download.shtml), [http://snap.stanford.edu/data](http://snap.stanford.edu/data), and [http://konect.cc/](http://konect.cc/). 
Additionally, the processed data, including partitioned data, is available in our [OneDrive Link](https://hkustgz-my.sharepoint.com/:f:/g/personal/xzhouby_connect_hkust-gz_edu_cn/EkEOQqUbSMZKioVFPdUvJisBSvhvzn0dR-ubJtpt7pmX5A?e=UWolbO). Please refer to the paper for details.

You can also use the example dataset in the directory `../../data` for testing.


## Dependency

1. `g++` 
2. `boost`
3. `omp`

All the codes are runnable after `cmake` and `make`: go to the corresponding directory, `cmake -DCMAKE_BUILD_TYPE=Release ./` and `make -j`.
