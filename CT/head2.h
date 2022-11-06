/*
 * head.h
 *
 *  Created on: 13 Oct 2022
 *      Author: zhangmengxuan
 */

#ifndef HEAD_H_
#define HEAD_H_

#include <stdio.h>
#include <math.h>
#include <vector>
#include <map>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <iostream>
#include <fstream>
#include <math.h>
#include <unordered_set>
#include <unordered_map>
#include <boost/thread/thread.hpp>
#include <chrono>
#include <string>
#include "Semaphore.h"
#define INF 999999999

using namespace std;
using namespace boost;

/*struct Node{//tree node
	vector<pair<int,pair<int,int>>> vert;//neighID/weight/count
	vector<int> pos, pos2;
	vector<int> dis;
	vector<int> ch;
	int height;
	int pa;//parent
	int uniqueVertex;
	Node(){
		vert.clear();
		pos.clear();
		dis.clear();
		ch.clear();
		pa = -1;
		uniqueVertex = -1;
		height = 0;
	}
};*/

struct Nei{
	int nid;
	int w;
	int c;
};

struct Node{//tree node
	vector<pair<int,pair<int,int>>> vert;//neighID/weight/count
	vector<pair<int,Nei>> neighInf;//posID,neighID,weight,count(for shortcut information maintenance)
	vector<int> pos, pos2;
	vector<int> dis, cnt;//the distance value and corresponding count number
    vector<int> vAncestor;//the ancestors, which is corresponding to dis
	//vector<set<int>> FromNode;
	set<int> changedPos;
	vector<bool> FN;//another succint way of FromNode
	set<int> DisRe;
	vector<int> ch;
	int height, hdepth;//hdepty is the deepest node that a vertex still exists
	int pa;//parent
	int uniqueVertex;//?vertex id of this tree node?
	vector<int> piv;//pivot vetex, used in path retrieval
	Node(){
		vert.clear();
		neighInf.clear();
		pos.clear();
		dis.clear();
		cnt.clear();
		ch.clear();
		pa = -1;
		uniqueVertex = -1;
		height = 0;
		hdepth = 0;
		changedPos.clear();
		FN.clear();
		DisRe.clear();
		piv.clear();
	}
};

struct NodeCore{//tree node in core-periphery
	vector<pair<int,int>> vert;//neighID/weight
	vector<int> pos, pos2;
	vector<int> dis;
	vector<int> disInterface;
	vector<int> ch;
	int height;
	int pa;//parent
	int uniqueVertex;//?key vertex id of this node
	int treeroot;//the subtree root
	NodeCore(){
		vert.clear();
		pos.clear();
		dis.clear();
		disInterface.clear();
		ch.clear();
		pa = -1;
		uniqueVertex = -1;
		height = 0;
		treeroot=-1;
	}
};

class Graph{
public:
	int nodenum;
	int edgenum;
	vector<vector<pair<int,int>>> Neighbor;
	vector<int> DD,DD2;//intermediate variable in Contraction
	void ReadGraph(string filename);
	int threadnum;

	//core-periphery decomposition
	int Width;
	vector<map<int,int>> Emap;
	vector<int> NodeOrder;//nodeID order
	vector<int> vNodeOrder;//order nodeID
	vector<int> rankCore;
	vector<vector<pair<int,int>>> NeighborConCore;
	vector<bool> existCore;
	int HighestOrder;
	vector<NodeCore> TreeCore;
	vector<int> EulerSeqCore;
	vector<int> toRMQCore;
	vector<vector<int>> RMQIndexCore;
	void H2HconCore();
	void CHconsCore();
	void deleteECore(int u,int v);
	void insertECore(int u,int v,int w);
	void makeTreeCore();
	int matchCore(int x,vector<pair<int,int>> &vert);

	//graph partition write & read
	int partiNum;
	void PartitionPreProcess();
	void WritePartition(string filename);
	void ReadPartition(string filename);
	void PartitionPostProcess();
	void OverlayGraph();
	//****************boundary vertex in each partition****************//
	vector<vector<int>> BoundVertex;
	vector<set<int>> BoundVertexSet;
	//****************vertex tag for query answering & vertex order****************//
	vector<int> CoreTag;
	vector<int> BoundTag;
	//****************adjacent list for partitions and core****************//
	vector<vector<vector<pair<int,int>>>> AdjaParti;
    vector<unordered_map<int,vector<pair<int,int>>>> AdjaPartiM;
	vector<vector<pair<int,int>>> AdjaCore;
	vector<map<int,int>> AdjaCoreMap;
	//****************Supported partition for boundary pair****************//
	vector<map<int,set<int>>> SuppPartiID;
	vector<map<int,set<int>>> SuppPartiIDReal;

	//Dijkstra in partition & core
	int Dijkstra(int ID1, int ID2,vector<vector<pair<int,int>>> &Neighbor);
	int DijkstraPath(int ID1, int ID2);
	int DijkstraParti(int ID1, int ID2, int pid);
	int DijkstraCore(int ID1, int ID2);
	int DijkstraCorePath(int ID1, int ID2);


	//Core Index
	vector<unordered_map<int,int>> IndexCore;
	vector<unordered_map<int,vector<int>>> PruningPointCore;//v {c,{u}}
	set<pair<int,int>> NoSupportedPairCore;

	//PLL index construction
	vector<vector<pair<int,int>>> Neighbors;
    unordered_map<int,vector<pair<int,int>>> NeighborsM;
	vector<unordered_map<int,int>> Label;
	vector<unordered_map<int,vector<int>>> PruningPointNew;//v {c,{u}}
	set<pair<int,int>> NoSupportedPair;
	void IndexConst1();
	void DijksPrune1(int nodeID, vector<pair<int,int>>& vp);
	void IncreasePLL(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label,vector<unordered_map<int,vector<int>>> &PruningPointNew,set<pair<int,int>> &NoSupportedPair);
	void DecreasePLL(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label);

	//WPSL index construction
	void IndexConstructMThread2New();
	bool DhopLableRefreshMulti2New(int step);
	void labelMultiThread2New(vector<unordered_map<int,int>>& newDhop, vector<int>& p,int step);
	void threadDistribute(vector<vector<int>>& processID);
	int ShortestDisQuery1(int ID1,int ID2,vector<int>& SupNode, int& d);
	void IncreasePSL(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label,vector<unordered_map<int,vector<int>>> &PruningPointNew,set<pair<int,int>> &NoSupportedPair);
	void DecreasePSL(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label);
	int ShortestDisQuery(int ID1,int ID2,vector<unordered_map<int,int>> &Label);
	int DisQueryVally(int ID1, int ID2, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label);
	int DisQueryPeak(int ID1, int ID2,vector<unordered_map<int,int>> &Label);
	int DisQueryLower1(int ID1, int ID2, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label);
	Semaphore* sm = new Semaphore(threadnum);
	Semaphore* sm1 = new Semaphore(1);
	unordered_map<int, Semaphore*> mSm;
    vector<Semaphore*> vSm;
	vector<vector<unordered_map<int,int>>> LabelStep;
	vector<unordered_map<int,int>> Dhop;
	vector<int> Dvectex;
	vector<bool> DvertexNew;
	vector<unordered_map<int,unordered_map<int,int>>> PruningPointStepNew;//v {c,{u, step}}
	struct tri{int s;
			   int c;
			   int t;
	};

	//H2H index construction
    void H2HindexParallel(bool ifParallel);
    void H2HindexP(pair<int,int> p, bool ifParallel);
	void H2Hindex(bool ifParallel);
    void H2Hindex();
	void CHindex();
	//Semaphore* sm = new Semaphore(threadnum);
	void CHindexMT();//multiple thread CH index construction
	void makeTree();
	void makeIndex();
    //intermediate variable and function used in the H2H index construction
    vector<vector<pair<int,pair<int,int>>>> NeighborCon;
    map<pair<int,int>,vector<int>> SCconNodes;
    vector<map<int, vector<int>>> SCconNodesMT;//multiple thread of SCconNodes
    vector<map<int,pair<int,int>>> E;
    vector<vector<int>> VidtoTNid;//one vertex exist in those tree nodes (nodeID--->tree node rank)
    vector<int> rank;
    int heightMax;
    vector<Node> Tree;
    vector<int> EulerSeq;//prepare for the LCA calculation
    vector<int> toRMQ;
    vector<vector<int>> RMQIndex;//?
    void insertEorder(int u,int v,int w);
    void deleteEorder(int u,int v);
    int match(int x,vector<pair<int,pair<int,int>>> &vert);
    void makeRMQ();
    void makeRMQDFS(int p, int height);
    void makeIndexDFS(int p, vector<int> &list);
    void NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x);
    void insertEMTorder(int u,int v,int w);
    //on the whole partitioned graph
    vector<vector<Node>> Trees;
    vector<vector<int>> toRMQs;
    vector<vector<vector<int>>> RMQIndexs;
    vector<vector<int>> ranks;
    vector<int> heightMaxs;
    vector<map<pair<int,int>,vector<int>>> SCconNodess;
    vector<vector<map<int, vector<int>>>> SCconNodesMTs;///
    vector<vector<vector<int>>> VidtoTNids;

	int QueryH2H(int ID1,int ID2);//shortest distance query with no partition
	int LCAQuery(int _p, int _q);
	int QueryH2HPartition(int ID1, int ID2, int PID);//query within partition
	int LCAQueryPartition(int _p, int _q, int PID);
    void DecreaseH2H(int a,int b, int newW, unordered_map<int,vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax);
	void EachNodeProBDis5(int child,vector<int>& line,set<int>& vertexIDChL, map<int,int>& checkedDis, vector<Node> &Tree, vector<int> &rank);
    void IncreaseH2H(int a, int b, int oldW, int newW, unordered_map<int,vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<int>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid);
	void eachNodeProcessIncrease1(int children, vector<int>& line, int& changelabel, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid);


	//index size
	void indexsizeCTDijk();
	void indexsizeCTH2H();
	void EffiCheck(string filename,int runtimes);

	//Query processing
	int Query(int ID1, int ID2);
	int QueryPartiCore(int ID1, int ID2);
	int QueryPartiParti(int ID1, int ID2);
	int QueryCore(int ID1, int ID2);

    /// Extension query
    bool ifParallel = true;
    bool extUpdate = false;
    vector<unordered_set<int>> PartiVertex;
    vector<unordered_map<int,int>> IndexExt;//extended 2-hop label of vertex in periphery
    void ExtensionIndex(pair<int,int> pidRange);
    void ExtensionIndexConstruct(bool ifParallel);
    int QueryPartiCoreExt(int ID1, int ID2);
    int QueryPartiPartiExt(int ID1, int ID2);

	//Correctness Check
	void CorrectnessCheck();
	void CorrectnessCheckCore();

	//Index update
	void Decrease(int a, int b, int oldW, int newW);
	void Increase(int a, int b, int oldW, int newW);

	void TestDataRead(string filename, double Ratio, vector<pair<pair<int,int>,pair<int,int>>>& Data);
    void ReadUpdate(string filename,vector<pair<pair<int,int>,int>>& TestData);

	//Big Graph Processing
	void StainingMethod(int ID);
	void ODGene(int num, string filename);
	void UpdateGene(int num, string filename);

	/*void WriteCoreGraph(string graphfile);
	void ReadCoreGraph(string filename);
	vector<set<int>> E;
	vector<pair<pair<int,int>,int>> EdgeEnum;*/
};

namespace benchmark {

#define NULLINDEX 0xFFFFFFFF

template<int log_k, typename k_t, typename id_t>
class heap {

public:

	// Expose types.
	typedef k_t key_t;
	typedef id_t node_t;

	// Some constants regarding the elements.
	//static const node_t NULLINDEX = 0xFFFFFFFF;
	static const node_t k = 1 << log_k;

	// A struct defining a heap element.
	struct element_t {
		key_t key;
		node_t element;
		element_t() : key(0), element(0) {}
		element_t(const key_t k, const node_t e) : key(k), element(e) {}
	};


public:

    //Enlarge
    /*inline void enlarge(node_t n) {
        n = 2 * n;
    }*/
	// Constructor of the heap.
	heap(node_t n) : n(0), max_n(n), elements(n), position(n, NULLINDEX) {
	}

	/*heap() {

	}*/

	// Size of the heap.
	inline node_t size() const {
		return n;
	}

	//elements in heap
	inline void elementsInHeap(vector<node_t> &eleinheap){
		eleinheap.assign(n, -1);
		for(int i=0;i<n;i++)
			eleinheap[i] = elements[i].element;
	}

	// Heap empty?
	inline bool empty() const {
		return size() == 0;
	}

	// Extract min element.
	inline void extract_min(node_t &element, key_t &key) {
		assert(!empty());

		element_t &front = elements[0];

		// Assign element and key.
		element = front.element;
		key = front.key;

		// Replace elements[0] by last element.
		position[element] = NULLINDEX;
		--n;
		if (!empty()) {
			front = elements[n];
			position[front.element] = 0;
			sift_down(0);
		}
	}

	inline key_t top() {
		assert(!empty());

		element_t &front = elements[0];

		return front.key;

	}

	inline node_t top_value() {

		assert(!empty());

		element_t &front = elements[0];

		return front.element;
	}

	// Update an element of the heap.
	inline void update(const node_t element, const key_t key) {
	   /* if(key)
        {
	        n = 2 * n;
        }*/
//		cout << "Update" << endl;
			if(element >= max_n)
			{
				max_n = 2 * max_n;

				element_t e(0,0);
				vector<element_t> elementsNull(elements.size(), e);
				elements.insert(elements.end(), elementsNull.begin(), elementsNull.end());
		//		elements.resize(max_n);a
				cout << "n:" << n << "\telement:" << element << "\tmax_n:" << max_n << endl;
				vector<node_t> positionNull(position.size(), NULLINDEX);
				position.insert(position.end(), positionNull.begin(), positionNull.end());

			}
		if (position[element] == NULLINDEX) {
			element_t &back = elements[n];
			back.key = key;
			back.element = element;
			position[element] = n;
			sift_up(n++);
		} else {
			node_t el_pos = position[element];
			element_t &el = elements[el_pos];
			if (key > el.key) {
				el.key = key;
				sift_down(el_pos);
			} else {
				el.key = key;
				sift_up(el_pos);
			}
		}
//		for(int i =0;i<max_n; i++)
//			cout << elements[i].element << "\t" << position[i] << endl;
	}


	// Clear the heap.
	inline void clear() {
		for (node_t i = 0; i < n; ++i) {
			position[elements[i].element] = NULLINDEX;
		}
		n = 0;
	}

	// Cheaper clear.
	inline void clear(node_t v) {
		position[v] = NULLINDEX;
	}

	inline void clear_n() {
		n = 0;
	}


	// Test whether an element is contained in the heap.
	inline bool contains(const node_t element) const {
		return position[element] != NULLINDEX;
	}

	inline key_t eleValue(const node_t element) const{//for the element in heap, return the key value
		return elements[position[element]].key;
	}


protected:

	// Sift up an element.
	inline void sift_up(node_t i) {
		assert(i < n);
		node_t cur_i = i;
		while (cur_i > 0) {
			node_t parent_i = (cur_i-1) >> log_k;
			if (elements[parent_i].key > elements[cur_i].key)
				swap(cur_i, parent_i);
			else
				break;
			cur_i = parent_i;
		}
	}

	// Sift down an element.
	inline void sift_down(node_t i) {
		assert(i < n);

		while (true) {
			node_t min_ind = i;
			key_t min_key = elements[i].key;

			node_t child_ind_l = (i << log_k) + 1;
			node_t child_ind_u = std::min(child_ind_l + k, n);

			for (node_t j = child_ind_l; j < child_ind_u; ++j) {
				if (elements[j].key < min_key) {
					min_ind = j;
					min_key = elements[j].key;
				}
			}

			// Exchange?
			if (min_ind != i) {
				swap(i, min_ind);
				i = min_ind;
			} else {
				break;
			}
		}
	}

	// Swap two elements in the heap.
	inline void swap(const node_t i, const node_t j) {
		element_t &el_i = elements[i];
		element_t &el_j = elements[j];

		// Exchange positions
		position[el_i.element] = j;
		position[el_j.element] = i;

		// Exchange elements
		element_t temp = el_i;
		el_i = el_j;
		el_j = temp;
	}



private:

	// Number of elements in the heap.
	node_t n;

	// Number of maximal elements.
	node_t max_n;

	// Array of length heap_elements.
	vector<element_t> elements;

	// An array of positions for all elements.
	vector<node_t> position;
};




template<int log_k, typename k_t, typename id_t, typename w_t>
class pHeap {

    public:

        // Expose types.
        typedef k_t key_t;
        typedef id_t node_t;
        typedef w_t weight_t;

        // Some constants regarding the elements.
        //static const node_t NULLINDEX = 0xFFFFFFFF;
        static const node_t k = 1 << log_k;

        // A struct defining a heap element.
        struct element_t {
            key_t key;
            node_t element;
            weight_t weight;
            element_t() : key(0), element(0), weight(0) {}
            element_t(const key_t k, const node_t e, const weight_t w) : key(k), element(e), weight(w) {}
        };


    public:

        // Constructor of the heap.
        pHeap(node_t n) : n(0), max_n(n), elements(n), position(n, NULLINDEX) {
        }

        // Size of the heap.
        inline node_t size() const {
            return n;
        }

        //elements in heap
        inline void elementsInHeap(vector<node_t> &eleinheap){
            eleinheap.assign(n, -1);
            for(int i=0;i<n;i++)
                eleinheap[i] = elements[i].element;
        }

        // Heap empty?
        inline bool empty() const {
            return size() == 0;
        }

        // Extract min element.
        inline void extract_min(node_t &element, key_t &key, weight_t &weight) {
            assert(!empty());

            element_t &front = elements[0];

            // Assign element, key and weight.
            element = front.element;
            key = front.key;
            weight = front.weight;

            // Replace elements[0] by last element.
            position[element] = NULLINDEX;
            --n;
            if (!empty()) {
                front = elements[n];
                position[front.element] = 0;
                sift_down(0);
            }
        }

        inline key_t top() {
            assert(!empty());

            element_t &front = elements[0];

            return front.key;

        }

        inline node_t top_value() {

            assert(!empty());

            element_t &front = elements[0];

            return front.element;
        }

        // Update an element of the heap.
        inline void update(const node_t element, const key_t key, const weight_t weight) {
			if(element > max_n)
			{
				max_n = 2 * max_n;

//				vector<element_t> elementsNull(elements.size(), NULLINDEX);
//				elements.insert(elements.end(), elementsNull.begin(), elementsNull.end());
				elements.resize(max_n);
				vector<node_t> positionNull(position.size(), NULLINDEX);
				position.insert(position.end(), positionNull.begin(), positionNull.end());

			}
            if (position[element] == NULLINDEX) {
                element_t &back = elements[n];
                back.key = key;
                back.element = element;
                back.weight = weight;
                position[element] = n;
                sift_up(n++);
            } else {
                node_t el_pos = position[element];
                element_t &el = elements[el_pos];
                if (key > el.key) {
                    el.key = key;
                    sift_down(el_pos);
                } else {
                    el.key = key;
                    sift_up(el_pos);
                }
            }
        }


        // Clear the heap.
        inline void clear() {
            for (node_t i = 0; i < n; ++i) {
                position[elements[i].element] = NULLINDEX;
            }
            n = 0;
        }

        // Cheaper clear.
        inline void clear(node_t v) {
            position[v] = NULLINDEX;
        }

        inline void clear_n() {
            n = 0;
        }


        // Test whether an element is contained in the heap.
        inline bool contains(const node_t element) const {
            return position[element] != NULLINDEX;
        }

        inline key_t eleValue(const node_t element) const{//for the element in heap, return the key value
            return elements[position[element]].key;
        }


    protected:

        // Sift up an element.
        inline void sift_up(node_t i) {
            assert(i < n);
            node_t cur_i = i;
            while (cur_i > 0) {
                node_t parent_i = (cur_i-1) >> log_k;
                if (elements[parent_i].key > elements[cur_i].key)
                    swap(cur_i, parent_i);
                else
                    break;
                cur_i = parent_i;
            }
        }

        // Sift down an element.
        inline void sift_down(node_t i) {
            assert(i < n);

            while (true) {
                node_t min_ind = i;
                key_t min_key = elements[i].key;

                node_t child_ind_l = (i << log_k) + 1;
                node_t child_ind_u = std::min(child_ind_l + k, n);

                for (node_t j = child_ind_l; j < child_ind_u; ++j) {
                    if (elements[j].key < min_key) {
                        min_ind = j;
                        min_key = elements[j].key;
                    }
                }

                // Exchange?
                if (min_ind != i) {
                    swap(i, min_ind);
                    i = min_ind;
                } else {
                    break;
                }
            }
        }

        // Swap two elements in the heap.
        inline void swap(const node_t i, const node_t j) {
            element_t &el_i = elements[i];
            element_t &el_j = elements[j];

            // Exchange positions
            position[el_i.element] = j;
            position[el_j.element] = i;

            // Exchange elements
            element_t temp = el_i;
            el_i = el_j;
            el_j = temp;
        }



    private:

        // Number of elements in the heap.
        node_t n;

        // Number of maximal elements.
        node_t max_n;

        // Array of length heap_elements.
        vector<element_t> elements;

        // An array of positions for all elements.
        vector<node_t> position;

    };
}

#endif /* HEAD_H_ */
