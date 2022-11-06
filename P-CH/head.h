/*
 * head.h
 *
 *  Created on: 19 Oct 2022
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
    string partiName;
	int edgenum;
	int nodenum;
	int pnum;//partition number
	int threadnum;//thread number

//Whole graph
	//vector<set<int>> E;
	map<pair<int,int>,int> EdgeWei;
	vector<vector<pair<int,int>>> Neighbor;
	vector<int> NodeOrder;//nodeID order
	vector<int> vNodeOrder;//order nodeID

	vector<pair<pair<int,int>,int>> CutEdges;//the cut edges
	vector<vector<int>> BoundVer;//the boundary vertex in each partition
	vector<set<int>> BoundVerSet;
	set<int> TotalBoundSet;//all boundary vertex
	vector<int> VtoParID;//vertex to partition ID; (0,1,...,k-1; k means multiple partitions; -1 means no partition)
	vector<unordered_map<int,int>> EtoParID;//edge to partition ID; (0,1,...,k-1; k means cut edge)

//Overlay graph
	//void OverlayGraphConstructPre();//construct before subgraph index construction
	void OverlayGraphConstructPost();//construct after subgraph index construction
	vector<unordered_map<int,int>> BedgePID;//Boundary edges' supported partition ID
	//void OverlayGraphConstructPrune();//construct after subgraph index construction with pruning techniques
	void CorrectnessCheckOverlay();//Correctness check of overlay graph

//CH index on a whole graph
	//Index construction
	void CHindex();
	Semaphore* sm = new Semaphore(threadnum);
	void CHindexMT();//multiple thread CH index construction
		//intermediate variable and function used in the H2H index construction
		vector<vector<pair<int,pair<int,int>>>> NeighborCon;
		vector<map<int, vector<int>>> SCconNodesMT;//multiple thread of SCconNodes
		vector<map<int,pair<int,int>>> E;

		void insertEorder(int u,int v,int w);
		void deleteEorder(int u,int v);
		void NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x);
		void insertEMTorder(int u,int v,int w);
		//on the whole partitioned graph
		vector<vector<vector<pair<int,int>>>> NeighborsParti;
		vector<vector<vector<pair<int,pair<int,int>>>>> NeighborCons;
		vector<vector<map<int, vector<int>>>> SCconNodesMTs;
		//on the overlay graph
		vector<vector<pair<int,int>>> NeighborsOverlay;
		vector<vector<pair<int,pair<int,int>>>> NeighborConOverlay;
		vector<map<int, vector<int>>> SCconNodesMTOverlay;

	//Query process
	int	QueryCH(int ID1, int ID2, vector<vector<pair<int,pair<int,int>>>> &NeighborCon);
	//Index update
	void CHdec(int a,int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbor,vector<vector<pair<int,pair<int,int>>>> &NeighborCon);//decrease
	void CHinc(int a,int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbor,vector<vector<pair<int,pair<int,int>>>> &NeighborCon,vector<map<int, vector<int>>> &SCconNodesMT);//increase
	void CHdecBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbor,vector<vector<pair<int,pair<int,int>>>> &NeighborCon);//decrease batch
	void CHincBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbor,vector<vector<pair<int,pair<int,int>>>> &NeighborCon,vector<map<int, vector<int>>> &SCconNodesMT);//increase batch
	vector<pair<pair<int,int>,pair<int,int>>> CHdecNew(int a,int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbor,vector<vector<pair<int,pair<int,int>>>> &NeighborCon);
	vector<pair<pair<int,int>,pair<int,int>>> CHincNew(int a,int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbor,vector<vector<pair<int,pair<int,int>>>> &NeighborCon,vector<map<int, vector<int>>> &SCconNodesMT);
//index update on partitioned graph
	void DecreaseSingle(int a, int b, int oldW, int newW);//process one update
	void IncreaseSingle(int a, int b, int oldW, int newW);//process one update

//Query processing on partitioned graph
	int Dijkstra(int ID1, int ID2);
	int DijkstraOverlay(int ID1, int ID2);
	int Query(int ID1, int ID2);
	int QueryCHInOut(int in, int out, int PID1);//one within partition, one is the boundary
	int QueryCHInIn(int ID1, int ID2, int PID);//in the same partition
	int QueryCHOutOut(int ID1, int ID2, int PID1, int PID2);//In different partitions
	void CorrectnessCheck();//Correctness check of query processing

//Graph Read
    void UpdateRead(string filename,vector<pair<pair<int,int>,int>>& TestData);
	void GraphRead(string filename);
	void OrderRead(string filename);
	void GraphPartitionRead(string filename);
	void EffiCheck(string filename, int runtimes);//read the random OD pairs
//Index Size
    void EffiCheck();
	void Indexsize();

    //Post Boundary
    void PartitionUpdate();
    void EffiCheckPost1();
    void EffiCheckPost(string filename,int runtimes);
    int QueryPost(int ID1, int ID2);
    void IndexsizePost();
    vector<vector<vector<pair<int,int>>>> NeighborsPartiPost;
    vector<vector<vector<pair<int,pair<int,int>>>>> NeighborConsPost;
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
