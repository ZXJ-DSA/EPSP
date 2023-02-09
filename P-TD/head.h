/*
 * head.h
 *
 *  Created on: 19 Sep 2022
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
#include "Timer.h"

#define INF 999999999
#define INCREASE 0
#define DECREASE 1
#define MIX 2

using namespace std;
using namespace boost;



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
	//vector<set<int>> FromNode;
	set<int> changedPos;
	vector<bool> FN;//another succint way of FromNode
	set<int> DisRe;
	vector<int> ch;
	int height, hdepth;//hdepty is the deepest node that a vertex still exists
	int pa;//parent
	int uniqueVertex;
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

class Graph{
public:
    string dataset;
    string algoName;
	int edgenum;
	int nodenum;
	int pnum;//partition number
	int threadnum=30;//thread number

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
	vector<unordered_map<int,vector<int>>> PruningPoint;//v {c,{u}}

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

//H2H index on a whole graph
	//Index construction
	void H2Hindex();
	void CHindex();
	Semaphore* sm = new Semaphore(threadnum);
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
    vector<vector<int>> RMQIndex;
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
    vector<vector<map<int, vector<int>>>> SCconNodesMTs;
    vector<vector<vector<int>>> VidtoTNids;
    //on the overlay graph
    vector<Node> TreeOverlay;
    vector<int> toRMQOverlay;
    vector<vector<int>> RMQIndexOverlay;
    vector<int> rankOverlay;
    int heightMaxOverlay;
    map<pair<int,int>,vector<int>> SCconNodesOverlay;
    vector<map<int, vector<int>>> SCconNodesOverlayMT;
    vector<vector<int>> VidtoTNidOverlay;

    vector<pair<pair<int,int>,int>> updateEdges;

	//Query process
	int QueryH2H(int ID1,int ID2);//shortest distance query with no partition
	int LCAQuery(int _p, int _q);
	//Index update
	void Decrease(int a,int b, int newW, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax);
	void EachNodeProBDis5(int child,vector<int>& line,set<int>& vertexIDChL, map<int,int>& checkedDis, vector<Node> &Tree, vector<int> &rank);
	void Increase(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<int>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid);
	void eachNodeProcessIncrease1(int children, vector<int>& line, int& changelabel, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid);
	void DecreaseBatch(vector<pair<pair<int,int>,int>>& wBatch, vector<vector<pair<int,int>>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax);//batch decrease for overlay graph
	void IncreaseBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<int>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid);
	void DecreaseSingle(int a, int b, int newW);//process one update
	void IncreaseSingle(int a, int b, int oldW, int newW);//process one update

//Query processing on partitioned graph
	int Dijkstra(int ID1, int ID2);
	int DijkstraOverlay(int ID1, int ID2);
	int QueryH2HPartition(int ID1, int ID2, int PID);//local process within one partition
	int LCAQueryPartition(int _p, int _q, int PID);
	int HopQueryOverlay(int ID1, int ID2);//on the overlay graph
	int LCAQueryOverlay(int _p, int _q);
	int HopQuery(int ID1, int ID2);//include all cases
	int HopQueryInOut(int in, int out);//one within partition, one is the boundary
	int HopQueryOutOut(int ID1, int ID2);//In different partitions
	int HopQueryInIn(int ID1, int ID2);//In same partition
	void CorrectnessCheck(int runtimes);//Correctness check of query processing

//Graph Read
	void GraphRead(string filename);
	void OrderRead(string filename);
	void GraphPartitionRead(string filename);
    void ReadUpdates(string filename);
	void EffiCheck(string filename, int runtimes);//read the random OD pairs
    void EffiCheck();

    //Index size
    void Indexsize();

    //Post boundary
    vector<vector<vector<pair<int,int>>>> NeighborsPartiPost;
    void PartitionUpdate();
    int HopQueryPost(int ID1, int ID2);
    void EffiCheckPost(string filename, int runtimes);
    void EffiCheckPost1();
    void IndexsizePost();

    //Pre boundary
    vector<int> dijkstra_candidate( int s, vector<int> &cands, vector<vector<pair<int,int>>> &graph );
    void PreBoundaryCompute(bool ifParallel);
    void preBoundaryPair(pair<int,int> p);

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
