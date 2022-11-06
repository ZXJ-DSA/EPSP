/*
 * head.h
 *
 *  Created on: 6 Jul 2022
 *      Author: zhangmengxuan
 */

#ifndef HEAD_H_
#define HEAD_H_

#include <stdio.h>
#include <math.h>
#include <vector>
#include <set>
#include <map>
#include<boost/algorithm/string/split.hpp>
#include<boost/algorithm/string/classification.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include<iostream>
#include<fstream>
#include<math.h>
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include <boost/thread/thread.hpp>
#include "Semaphore.h"

#define INF 999999999

using namespace std;
using namespace boost;

class Graph{
public:
	int edgenum;
	int nodenum;
	int pnum;//partition number
	int threadnum;//thread number
    string partiName;

//Whole graph
	//vector<set<int>> E;
	map<pair<int,int>,int> EdgeWei;
	vector<vector<pair<int,int>>> Neighbors;
	vector<int> NodeOrder;//nodeID order
	vector<int> vNodeOrder;//order nodeID

	vector<pair<pair<int,int>,int>> CutEdges;//the cut edges
	vector<vector<int>> BoundVer;//the boundary vertex in each partition
	vector<set<int>> BoundVerSet;
	set<int> TotalBoundSet;//all boundary vertex
	vector<int> VtoParID;//vertex to partition ID; (0,1,...,k-1; k means multiple partitions; -1 means no partition)
	vector<unordered_map<int,int>> EtoParID;//edge to partition ID; (0,1,...,k-1; k means cut edge)

	vector<unordered_map<int,int>> Label;
	vector<unordered_map<int,vector<int>>> PruningPointNew;//v {c,{u}}
	set<pair<int,int>> NoSupportedPair;

//Partitioned sub-graph
	vector<vector<vector<pair<int,int>>>> NeighborsParti;
	vector<vector<unordered_map<int,int>>> LabelParti;
	vector<vector<unordered_map<int,vector<int>>>> PruningPointParti;
	vector<set<pair<int,int>>> NoSupportedParti;

//Overlay graph
	//void OverlayGraphConstructPre();//construct before subgraph index construction
	void OverlayGraphConstructPost();//construct after subgraph index construction
	vector<map<int,vector<int>>> BedgeAllPID;//Boundary edges' all possible supported partition ID
	vector<map<int,pair<int,set<int>>>> BedgePID;//Boundary edges' supported partition ID
	//void OverlayGraphConstructPrune();//construct after subgraph index construction with pruning techniques
	void CorrectnessCheckOverlay();//Correctness check of overlay graph
	vector<vector<pair<int,int>>> NeighborsOverlay;
	vector<unordered_map<int,int>> LabelOverlay;
	vector<unordered_map<int,vector<int>>> PruningPointOverlay;//v {c,{u}}
	set<pair<int,int>> NoSupportedOverlay;

//PLL Index construction
	void IndexConst();
	void DijksPrune(int nodeID, vector<pair<int,int>>& vp);
//    void DijksPruneP(int nodeID, vector<pair<int,int>>& vp, vector<unordered_map<int,vector<int>>> & PruningPointNew);
	int ShortestDisQuery(int ID1,int ID2,vector<int>& SupNode, int& d);
/*//Index construction multiple threads
	void IndexConst1();
	void IndexConstMT(pair<int,int> p);
	void DijksPruneMT(int nodeID, vector<pair<int,int>>& vp,int partiID);
	int ShortestDisQueryMT(int ID1,int ID2,vector<int>& SupNode, int& d,int partiID);*/

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

//Query processing
	int Dijkstra(int ID1, int ID2, vector<vector<pair<int,int>>> &Neighbors);
	int HopQueryLocal(int ID1, int ID2, vector<unordered_map<int,int>> &Label);//within one partition or in overlay graph
	int HopQuery(int ID1, int ID2);//include all cases
	int HopQueryInOut(int in, int out);//one within partition, one is the boundary
	int HopQueryOutOut(int ID1, int ID2);//both are boundary
	int HopQueryInIn(int ID1, int ID2);//both are within partition
	void CorrectnessCheck();//Correctness check of query processing

//Index Update
	void DecreaseSingle(int a, int b, int oldW, int newW);//process one update
	void IncreaseSingle(int a, int b, int oldW, int newW);//process one update
	void IncreasePLL(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label,vector<unordered_map<int,vector<int>>> &PruningPointNew,set<pair<int,int>> &NoSupportedPair);
	void DecreasePLL(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label);

/*//distance function supporting increase
	int DisQueryVally(int ID1, int ID2, vector<vector<pair<int,int>>> &Neighbors, vector<unordered_map<int,int>> &Label);
	int DisQueryPeak(int ID1, int ID2, vector<unordered_map<int,int>> &Label);
	int DisQueryLower1(int ID1, int ID2, vector<vector<pair<int,int>>> &Neighbors, vector<unordered_map<int,int>> &Label);*/

//Graph Read
    void UpdateRead(string filename,vector<pair<pair<int,int>,int>>& TestData);
	void GraphRead(string filename);
	void OrderRead(string filename);
	void GraphPartitionRead(string filename);

//Efficiency check
	void EffiCheck(string filename, int runtimes);//read the random OD pairs
    void EffiCheck();

//Index size
	void Indexsize();

    //Post-boundary strategy
//The difference between no-boundary and post-boundary is the additional subgraph index L_i^', the query efficiency of query OD within one partition, update of L_i^'
    vector<vector<vector<pair<int,int>>>> NeighborsPartiPost;
    vector<vector<unordered_map<int,int>>> LabelPartiPost;
    vector<vector<unordered_map<int,vector<int>>>> PruningPointPartiPost;
    vector<set<pair<int,int>>> NoSupportedPartiPost;
    //query processing
    void CorrectnessCheckPost();
    int HopQueryPost(int ID1, int ID2);
    //index size
    void IndexsizePost();
    //efficiency check
    void EffiCheckPost(string filename,int runtimes);
    void EffiCheckPost1();
    //index construction
    void PartitionUpdate();
    //index update
    void insertion(int a, int b, int newW, vector<vector<pair<int,int>>> &Neighbors, vector<unordered_map<int,int>> &Label);//edge insertion
    //void DecreasePost(int a, int b, int oldW, int newW);//process one update
    //void IncreasePost(int a, int b, int oldW, int newW);//process one update

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

	// Constructor of the heap.
	heap(node_t n) : n(0), max_n(n), elements(n), position(n, NULLINDEX) {
	}

	heap() {

	}

	// Size of the heap.
	inline node_t size() const {
		return n;
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
