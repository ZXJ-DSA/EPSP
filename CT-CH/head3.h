/*
 * head.h
 *
 *  Created on: 13 Oct 2022
 *      Author: zhangmengxuan
 */

#ifndef HEAD3_H_
#define HEAD3_H_

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
#include <stack>
#include <boost/thread/thread.hpp>
#include <chrono>
#include <string>
#include "Semaphore.h"


#define INF 99999999

using namespace std;
using namespace boost;

struct Nei{
	int nid;
	int w;
	int c;
};

struct Node{//tree node
	vector<pair<int,pair<int,int>>> vert;//neighID/weight/count(how many ways can lead to this super edge weight)
	vector<pair<int,Nei>> neighInf;//posID,neighID,weight,count(for shortcut information maintenance)
	vector<int> pos;
//    int posNum;//the number of pos, i.e., how many adjacent vertices are ancestors
	vector<int> dis, cnt;//the distance value and corresponding count number (i.e., record how many path has the shortest distance)
    vector<int> vAncestor;//the ancestors, which is corresponding to dis
//    vector<int> disInf;//the distances from this vertex to the interface vertices
    map<int,int> disInf;//the distances from this vertex to the interface vertices
	//vector<set<int>> FromNode;
	set<int> changedPos;
	vector<bool> FN;//another succinct way of FromNode, whether this distance label is directly obtained from shortcuts (vert)
    map<bool,bool> FNInf;//whether the interface distance is obtained from shortcuts (vert)
	set<int> DisRe;//record the vertex id that the distance label should be updated
    set<int> DisReInf;//record the vertex id that the interface label should be updated
	vector<int> ch;
	int height, hdepth;//hdepty is the deepest node that a vertex still exists
	int pa;//parent, the pa of root vertex is 0
	int uniqueVertex;//?vertex id of this tree node?
	vector<int> piv;//pivot vetex, used in path retrieval
    int treeroot;//the tree id of subtree root, i.e., rank[x]
	Node(){
		vert.clear();
		neighInf.clear();
		pos.clear();
		dis.clear();
        disInf.clear();
		cnt.clear();
        vAncestor.clear();
		ch.clear();
		pa = -1;
		uniqueVertex = -1;
		height = 0;
		hdepth = 0;
		changedPos.clear();
		FN.clear();
        FNInf.clear();
		DisRe.clear();
        DisReInf.clear();
		piv.clear();
        treeroot=-1;
	}
};

class Graph{
public:
    string dataset;
	int nodenum;
	int edgenum;
	vector<vector<pair<int,int>>> Neighbor;//original graph
	vector<int> DD,DD2;//intermediate variable in Contraction
	void ReadGraph(string filename);
	int threadnum;

	//core-periphery decomposition
	int Width;
//	vector<map<int,int>> Emap;//map structure for copied graph
	vector<int> NodeOrder;//nodeID order
	vector<int> vNodeOrder;//order nodeID
	vector<bool> existCore;
	int HighestOrder;
//	vector<NodeCore> TreeCore;
	vector<int> EulerSeqCore;
	vector<int> toRMQCore;
	vector<vector<int>> RMQIndexCore;
	void H2HconCore();
	void CHconsCore();
	void deleteECore(int u,int v);
	void insertECore(int u,int v,int w);
	void makeTreeCore();
	int matchCore(int x,vector<pair<int,pair<int,int>>> &vert);//vector<pair<int,int>> &vert

	//graph partition write & read
	int partiNum;
	void PartitionPreProcess();
	void WritePartition(string filename);
	void ReadPartition(string filename);
	void PartitionPostProcess();
	void OverlayGraph();
    void OverlayGraphProcess();//get the partition info and overlay graph
	//****************boundary vertex in each partition****************//
	vector<vector<int>> BoundVertex;//the interface vertex + root vertex
	vector<set<int>> BoundVertexSet;//only interface vertex, excluding root vertex
	//****************vertex tag for query answering & vertex order****************//
	vector<int> CoreTag;//-1 indicates core vertex and root vertex without children, i>=0 indicates non-core vertex (i.e., the root and inner-partition vertex) and which partition it belongs to
    vector<pair<bool,set<int>>> BoundTag;//first element: 1 (true) indicates boundary vertex (i.e., interface vertex), 0 (false) indicates non-boundary vertex; second element contains the partitions it belongs to
//    vector<bool> BoundTag;//1 (true) indicates boundary vertex (i.e., interface vertex), 0 (false) indicates non-boundary vertex
	//****************adjacent list for partitions and core****************//
	vector<vector<pair<int,int>>> AdjaCore;
	vector<map<int,int>> AdjaCoreMap;
    vector<map<int,int>> AdjaCoreMapOld;
	//****************Supported partition for boundary pair****************//
	vector<map<int,map<int,int>>> SuppPartiID;//<ID1,<ID2,<pid,weight>>>, the set of partitions that contains these two interface vertex
	vector<map<int,pair<int,set<int>>>> SuppPartiIDReal;//ID1,<ID2,<weight,set<pid>>>>, the set of partitions that support the super edge of these two interface vertex

    set<int> CoreVertex;

	//Dijkstra in partition & core
	int Dijkstra(int ID1, int ID2,vector<vector<pair<int,int>>> &Neighbor);
	int DijkstraPath(int ID1, int ID2);
    int DijkstraPathCore(int ID1, int ID2);
	int DijkstraParti(int ID1, int ID2, int pid);
	int DijkstraCore(int ID1, int ID2);
	int DijkstraCorePath(int ID1, int ID2);


	//PLL index construction
	vector<unordered_map<int,int>> Label;
    vector<unordered_map<int,set<int>>> PruningPoint;
//	vector<unordered_map<int,vector<int>>> PruningPointNew;//v {c,{u}}
//	set<pair<int,int>> NoSupportedPair;
	void IndexConst1();
	void DijksPrune1(int nodeID, vector<pair<int,int>>& vp);
	void IncreasePLL(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label,vector<unordered_map<int,vector<int>>> &PruningPointNew,set<pair<int,int>> &NoSupportedPair);
	void DecreasePLL(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label);
    void CoreGraphDebug(string graphfile);

	//WPSL index construction
	void IndexConstructMThread2New();
	bool DhopLableRefreshMulti2New(int step);
	void labelMultiThread2New(vector<unordered_map<int,int>>& newDhop, vector<int>& p,int step);
	void threadDistribute(vector<vector<int>>& processID);
	int ShortestDisQuery1(int ID1,int ID2,vector<int>& SupNode, int& d);
    void IncreasePSL(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label,vector<unordered_map<int,set<int>>> &PruningPointNew);
//	void IncreasePSL(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label,vector<unordered_map<int,vector<int>>> &PruningPointNew,set<pair<int,int>> &NoSupportedPair);
	void DecreasePSL(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label);
	int ShortestDisQuery(int ID1,int ID2,vector<unordered_map<int,int>> &Label);
    int DisQueryVally2(int ID1, int ID2, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label);
	int DisQueryVally(int ID1, int ID2, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label);
    pair<int,int> DisQueryPeak(int ID1, int ID2,vector<unordered_map<int,int>> &Label);
//	int DisQueryPeak(int ID1, int ID2,vector<unordered_map<int,int>> &Label);
	int DisQueryLower1(int ID1, int ID2, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label);
	Semaphore* sm = new Semaphore(threadnum);
	Semaphore* sm1 = new Semaphore(1);
	unordered_map<int, Semaphore*> mSm;
    vector<Semaphore*> vSm;
//	vector<vector<unordered_map<int,int>>> LabelStep;
	vector<unordered_map<int,int>> Dhop;
	vector<int> Dvectex;
	vector<bool> DvertexNew;
//	vector<unordered_map<int,unordered_map<int,int>>> PruningPointStepNew;//v {c,{u, step}}
	struct tri{int s;
			   int c;
			   int t;
	};

	//H2H index construction
    void H2HindexParallel(bool ifParallel);
    void H2HindexP(pair<int,int> p, bool ifParallel);
	void H2Hindex(bool ifParallel);

    //intermediate variable and function used in the H2H index construction
    vector<vector<pair<int,pair<int,int>>>> NeighborCon;
    vector<vector<pair<int,pair<int,int>>>> NeighborConCH;//contracted periphery graph for CH
//    vector<map<int, vector<int>>> SCconNodesMT;//supportive vertex, multiple thread of SCconNodes
    vector<map<int, vector<pair<int,int>>>> SCconNodesMT;//<ID1,<ID2,<x,weight>>>supportive vertex, multiple thread of SCconNodes
    vector<map<int,pair<int,int>>> E;
    vector<vector<int>> VidtoTNid;//record the child tree nodes whose vert neighbors contain this tree node (nodeID--->tree node rank)
    vector<int> rank;//rank[v]>0 indicates non-core vertex
    int heightMax;
    vector<Node> Tree;
//    vector<NodeCore> Tree;
    vector<int> EulerSeq;//prepare for the LCA calculation, EulerSeq is the Euler tour, i.e., E[1,...,2n-1]
    vector<int> toRMQ;//toRMQ[p] records the first occurrence of node p in the Euler tour, i.e., R[p]
    vector<vector<int>> RMQIndex;//?
    void insertEorder(int u,int v,int w);
    void deleteEorder(int u,int v);
    int match(int x,vector<pair<int,pair<int,int>>> &vert);
//    void makeRMQ();
    void makeRMQCore();
    void makeRMQDFS(int p, int height);
    void makeRMQDFSCore(int p, int height);
    void makeIndexDFSCore(int p, vector<int> &list, vector<int> & interface, map<int,unordered_map<int,int>> & disInfs);
//    void makeIndexDFS(int p, vector<int> &list);
    void NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x);
    void insertEMTorder(int u,int v,int w);


	int QueryH2H(int ID1,int ID2);//shortest distance query with no partition
	int LCAQuery(int _p, int _q);
	int QueryPeripheryTree(int ID1, int ID2, int PID);//query within partition
	int LCAQueryPartition(int _p, int _q, int PID);
    void DecreaseH2HNew(int a,int b, int newW, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifParallel);//either a or b is periphery vertex
    void InterfacePropagate(int child, vector<int>& interfaces, vector<Node> &Tree, bool ifIncrease);
    void InterfacePropagateParallel(pair<int,int> pRange, vector<int>& pids, bool ifIncrease);
	void EachNodeProBDis5New(int child,vector<int>& line, vector<int>& interfaces, set<int>& vertexIDChL, map<int,int>& checkedDis, vector<Node> &Tree, vector<int> &rank);
    void IncreaseH2HNew(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, bool ifParallel);
    void eachNodeProcessIncrease1New(int children, vector<int>& line, vector<int>& interfaces, int& changelabel, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid, int lowestH);


	//index size
	void indexsizeCTDijk();
	void indexsizeCTH2H();
	void EffiCheck(string filename,int runtimes);

	//Query processing
    int QueryDebug(int ID1, int ID2);
	int Query(int ID1, int ID2);
	int QueryPartiCore(int ID1, int ID2);
    int QueryPartiCoreDebug(int ID1, int ID2);
    int QuerySameParti(int ID1, int ID2);
	int QueryPartiParti(int ID1, int ID2);
	int QueryCore(int ID1, int ID2);
    int QueryCoreDebug(int ID1, int ID2);

    /// Extension query
    bool ifParallel = true;
    bool extUpdate = false;
    vector<unordered_set<int>> PartiVertex;//only insert in-partition vertex (i.e., periphery vertex, excluding the interface vertex) in this case
    vector<unordered_map<int,int>> IndexExt;//extended 2-hop label of vertex in periphery
    vector<bool> PartiUpdateExt;//indicates whether the periphery's extended index should be updated
    void ExtensionIndex(pair<int,int> pidRange, bool ifIncrease);
    void ExtensionIndex2(pair<int,int> pidRange, bool ifIncrease, vector<int>& partiForUpdate);
    void ExtensionIndexConstruct(bool ifParallel, bool ifIncrease);
    void ExtensionIndexUpdate(bool ifParallel, bool ifIncrease, vector<int>& partiForUpdate);
    int QueryPartiCoreExt(int ID1, int ID2);
    int QueryPartiCoreExtDebug(int ID1, int ID2);
    int QueryPartiPartiExt(int ID1, int ID2);
    int QueryPartiPartiExtDebug(int ID1, int ID2);

    /// New implementation
    void MDEContract();//contraction according to MDE
    void Create_tree();
    void TreeLabelCompute(pair<int,int> pidRange, vector<int> & pidRanks);
    void Compute_tree_label(bool ifParallel);
    void Construct_tree(bool ifParallel);
    void Construct_core();

	//Correctness Check
	void CorrectnessCheck(int runtimes);
	void CorrectnessCheckCore();
    void CorrectnessCheckCoreAll();

	//Index update
	void Decrease(int a, int b, int oldW, int newW);
	void Increase(int a, int b, int oldW, int newW);

	void TestDataRead(string filename, double Ratio, vector<pair<pair<int,int>,pair<int,int>>>& Data);
    void ReadUpdate(string filename,vector<pair<pair<int,int>,int>>& TestData);
    void ReadUpdate2(string filename,vector<pair<pair<int,int>,pair<int,int>>>& TestData);

	//Big Graph Processing
	void StainingMethod(int ID);
	void ODGene(int num, string filename);
	void UpdateGene(int num, string filename);
    void ReadLabels(string file);
    void WriteLabels(string file);
	void WriteCoreGraph(string graphfile);
	void ReadCoreGraph(string filename);
    vector<int> DFS_CC(vector<map<int,int>> & Edges, set<int> set_A, set<int> & set_B, int nodenum);

    int MinSpanTree(vector<vector<pair<int,int>>> & Neighbors);

    //New added functions
    void QueryPeriphery_CH(int ID1, set<int> todo, map<int,int> & results);//single-source CH query
    int QueryPeriphery_CH(int ID1, int ID2);
    int QuerySameParti_CH(int ID1, int ID2);
    int QueryPartiParti_CH(int ID1, int ID2);
    int QueryPartiCore_CH(int ID1, int ID2);
    int Query_CH(int ID1, int ID2);
    void EffiCheck_CH(string filename,int runtimes);
    void indexsizeCHH2H();
    void Create_Partitions();//get the partition
    void CorrectnessCheck_CH(int runtimes);

    void CHdecreaseBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch);
    void CHincreaseBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch);
    void Decrease_CH(int a, int b, int oldW, int newW);
    void Increase_CH(int a, int b, int oldW, int newW);
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

#endif /* HEAD3_H_ */
