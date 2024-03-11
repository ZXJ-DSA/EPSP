## Introduction
This is the source code of hierarchical PSP indexes. Please refer to the paper for the algorithm details.

## Algorithms
After `cmake -DCMAKE_BUILD_TYPE=Release ./` and `make -j`, please run `./HPSP` to obtain the usage instruction.

**For the hierarchical PSP index**,

1. arg 3 specifies the leaf node size.
1. arg 4 specifies the PSP index: 1 for G-Tree, 4 for N-TS-HP.
1. arg 5 specifies the task: 1 for index construction; 2 for query processing; 3 for index update.



Examples are shown below.

* Build G-Tree index with leaf node size 128 on NY:
`./HPSP ../../Datasets NY 128 1 1 1 10 150`
* Test G-Tree index query processing with leaf node size 128 on NY:
`./HPSP ../../Datasets NY 128 1 2 1 10 150`
* Test G-Tree index update with leaf node size 128 on NY:
`./HPSP ../../Datasets NY 128 1 3 1 10 150`



## Data
The datasets of this paper source from [http://www.dis.uniroma1.it/challenge9/download.shtml](http://www.dis.uniroma1.it/challenge9/download.shtml), [http://snap.stanford.edu/data](http://snap.stanford.edu/data), and [http://konect.cc/](http://konect.cc/). 
Additionally, the processed data, including partitioned data, is available in our [OneDrive Link](https://hkustgz-my.sharepoint.com/:f:/g/personal/xzhouby_connect_hkust-gz_edu_cn/EkEOQqUbSMZKioVFPdUvJisBSvhvzn0dR-ubJtpt7pmX5A?e=UWolbO). Please refer to the paper for details.

You can also use the example dataset in the directory `../../data` for testing.


## Dependency

1. `g++` 
2. `boost`
3. `omp`
4. `metis`

All the codes are runnable after `cmake` and `make`: go to the corresponding directory, `cmake -DCMAKE_BUILD_TYPE=Release ./` and `make -j`.
