#ifndef BATCHHLWW_H_
#define BATCHHLWW_H_

#include <sys/time.h>
#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <map>
#include <unordered_map>
#include <set>
#include <stack>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <boost/thread/thread.hpp>

#include "two_layer_queue.h"
#include "Timer.h"

//
// Implementation for weighted and undirected graphs.
//

#define INF 99999999
#define INCREASE 0
#define DECREASE 1
#define MIX 2

bool ifDebug = false;//

using namespace std;

class HighwayLabellingWW {
public:
  // Constructs labelling for a graph, given as a list of edges.
//  HighwayLabellingWW(std::string filename, int k);
  HighwayLabellingWW();
  ~HighwayLabellingWW();

  void ReadGraph(string filename, int k);
    void ConstructHighwayLabellingNewP(pair<int,int> vi, vector<int> & topk);
  void ConstructHighwayLabellingNew(int i, vector<int> & topk);
  void ConstructHighwayLabelling(int i, vector<int> & topk);
  void BuildIndex(vector<int> & topk, bool ifParallel);

  void UpdateLabelling(std::string filename, int m);
  void UpdateLabellingW(std::string filename, int m, int updateType, int updateBatch, int updateVolume,bool ifParallel);


  void BHL_Plus(std::vector<std::pair<std::string, std::pair<int, int> > > & updates);
  void BHL_PlusW(std::vector<std::pair<std::string, std::pair<int, int> > > & updates);
  void BatchSearch(vector<pair<string, pair<int, int> > > & updates, int rid , vector<int> & d_G, vector<int> & A, set<int> & V_AFF);
  void BatchSearch(vector<pair<pair<int,int>,int>> & updates, int rid , vector<int> & d_G, vector<int> & A, set<int> & V_AFF);
  void BHL(std::vector<std::pair<std::string, std::pair<int, int> > > & updates);
  void BHLW(std::vector<std::pair<std::string, std::pair<int, int> > > & updates);
  void BHLW(vector<pair<pair<int,int>,int>> & updates);
    void BHLWP(vector<pair<pair<int,int>,int>> & updates, bool ifParallel);
    void UpdateLandmark(pair<int,int> landmarkRange, vector<pair<pair<int,int>,int>> & updates);

  bool prunable(int i, int u, int *temp, int *A);
  bool prunable(int i, int u, vector<int> & temp, vector<int> & A);
  bool prunableW(int i, int u, vector<int> & temp, vector<int> & A);
  int query(int r, int v);
  int ldPair(int d, bool landmark);
  int ldDist(int ld);
  int ldAdd(int ld, bool landmark);

  void allocate();

    vector<int> OrderRead(string filename);
  void SelectLandmarks_HD(vector<int> & topk, string graphFile);
  int LabellingSize();

  int min(int a, int b);
  void RemoveLandmarks(vector<int> & topk);
  void HC_UB_naiveP(int s, int t, bool ifFull, int & dis);
  int HC_UB_naive(int s, int t);
  int HC_UB_naive_Full(int s, int t);
//  void QueryDistance(std::string pairs, std::string output, int runtimes);
  void EfficencyTest(string pairs, int runtimes);

  void storeLabelling(std::string filename);
  void loadLabelling_Full(std::string filename, vector<int> & topk);
  void loadLabelling_Pruned(std::string filename);

  void SaveLandmarks(string filename, vector<int> & topk);
  void LoadLandmarks(string filename, vector<int> & topk);
  int BFS(int ID1, int ID2, vector<vector<int> > & graph);
  int BFS(int ID1, int ID2, vector<unordered_map<int,int> > & graph);
  int Bi_BFS(int ID1, int ID2, vector<vector<int> > & graph);
  int Dijkstra(int ID1, int ID2, vector<unordered_map<int,int> > & Nodes);
  int BiDijkstra(int ID1, int ID2, vector<unordered_map<int,int> > & Nodes);
    int BiDijkstraBounded(int ID1, int ID2, vector<unordered_map<int,int>> & Nodes, int dis_upper);
  void BiDijkstraP(int ID1, int ID2, vector<unordered_map<int,int>> & Nodes, int & dis);
  void CorrectnessCheck(int runtimes,bool ifFull);
  void CorrectnessCheckW(int runtimes,bool ifFull);
  void DFS_CC(vector<unordered_map<int,int>> & Edges, set<int> & set_A);
  void DFS_CC(vector<vector<int>> & Edges, set<int> & set_A);
  void CheckCC();
//  int BatchHL_Query(int s,int t, vector<vector<int> > & graph,bool ifFull);
  int BatchHL_Query(int s,int t, vector<unordered_map<int,int> > & graph,bool ifFull);
  int BatchHL_QueryW(int s,int t, vector<unordered_map<int,int> > & graph,bool ifFull);
  void ReadUpdates(string filename);
  void verticesUpdate();


//private:
  int V;  // total number of vertices
  long long int E; // total number of edges
  int K; // total number of landmarks
  string dataset;
  int threadnum = 10;

  vector<vector<int>> distances, distances_1;
  vector<vector<int>> highway, highway_1;
  vector<vector<int>> vertices;//mapped vertex id of landmark
  vector<int> C;
  std::vector<unordered_map<int,int> > adj;//adjacent lists for weighted graph (ID2,weight)
  std::vector<unordered_map<int,int> > adjOld;//adjacent lists for weighted graph
  vector<unordered_map<int,int> > adjOldO;//adjacent lists for weighted graph
  std::map<int, int> landmarks;//map vertex id of landmark to simple id
  set<int> landmarkSet;//set of landmarks (vertex id)
  vector<int> landmarkTopK;//map from simple id to landmark vertex id
  vector<pair<pair<int,int>,int>> updateEdges;//vector of update edges


  double GetCurrentTimeSec() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
  }
  
  long GetCurrentTimeMilliSec() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * 1000LL + tv.tv_usec / 1000;
  }

  // Statistics
  double time_, time_querying_sec_, time_querying_millisec_;
};

HighwayLabellingWW::HighwayLabellingWW() { }

HighwayLabellingWW::~HighwayLabellingWW() { }

//created by Mengxuan, modified by Xinjie
namespace benchmark {

#define NULLINDEX 0xFFFFFFFF

    template<int log_k, typename id_t, typename k_t >//id,value
    class heap {
    public:

        // Expose types.
        typedef k_t key_t;
        typedef id_t node_t;

        // Some constants regarding the elements.
        //static const node_t NULLINDEX = 0xFFFFFFFF;
        static const node_t k = 1 << log_k;//equals k = 1*2^log_k, usually log_k = 2, so k = 4

        // A struct defining a heap element.
        struct element_t {
            key_t key;
            node_t element;

            element_t() : key(0), element(0) {}

            element_t(const key_t k, const node_t e) : key(k), element(e) {}
        };


        //public:

        // Constructor of the heap.
        heap(node_t n) : n(0), max_n(n), elements(n), position(n, NULLINDEX) {//n is the number of elements in current
            // state, max_n is the size of heap
        }

        heap(): n(0), max_n(0), elements(0), position(0, NULLINDEX) {}

        ~heap(){}

        // Risize the heap
        inline void resize(node_t a){
            n = 0; max_n = a;
            elements.resize(a);
            position.resize(a, NULLINDEX);
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
            --n;//n=n-1
            if (!empty()) {
                front = elements[n];//elements[n] is the top element
                position[front.element] = 0;//make its position valid, it is also the smallest one
                sift_down(0);
            }
        }

        inline key_t top_key() {//get the key, i.e. minimal cost
            assert(!empty());

            element_t &front = elements[0];

            return front.key;

        }

        inline node_t top_value() {//get the value, i.e. id number of minimal cost

            assert(!empty());

            element_t &front = elements[0];

            return front.element;
        }

        // Update an element of the heap.
        inline void update(const node_t element, const key_t key) {

            if (position[element] == NULLINDEX) {//if originally NULL
                element_t &back = elements[n];//add new element to position n
                back.key = key;
                back.element = element;
                position[element] = n;//set position id to n
                sift_up(n++);
            } else {//if already valid, update the value
                node_t el_pos = position[element];//position information
                element_t &el = elements[el_pos];//get the element
                if (key > el.key) {//update the element
//                if (key > el.key || (key <= el.key && element > el.element)) {//update the element || (elements[parent_i].key <= elements[cur_i].key && elements[parent_i].element > elements[cur_i].element)
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

        // Cheaper erase.
        inline void erase(node_t v) {
            position[v] = NULLINDEX;
        }

        inline void clear_n() {
            n = 0;
        }

        // Test whether an element is contained in the heap.
        inline bool contains(const node_t element) const {
            return position[element] != NULLINDEX;
        }

        //return current elements information
        void get_elements(std::vector<std::pair<int,int>> &e_vector){
            std::pair<int,int> temp_pair;

            for(int i=0;i<n;i++){
                temp_pair.first = elements[i].key;
                temp_pair.second = elements[i].element;
                e_vector.push_back(temp_pair);
            }
        }

    protected:

        // Sift up an element.
        inline void sift_up(node_t i) {
            assert(i < n);
            node_t cur_i = i;
            while (cur_i > 0) {
                node_t parent_i = (cur_i - 1) >> log_k;//equals (cur_i - 1)/(2^log_k)
                if (elements[parent_i].key > elements[cur_i].key)//compare with parent node, if smaller, then swap
//                if (elements[parent_i].key > elements[cur_i].key || (elements[parent_i].key <= elements[cur_i].key && elements[parent_i].element > elements[cur_i].element))//compare with parent node, if smaller, then swap
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

                node_t child_ind_l = (i << log_k) + 1;//equals i*2^log_k + 1
                node_t child_ind_u = std::min(child_ind_l + k, n);//equals min(child_ind_l+4,n)

                for (node_t j = child_ind_l; j < child_ind_u; ++j) {
                    if (elements[j].key < min_key) {
//                    if (elements[j].key < min_key || (elements[j].key >= min_key && elements[j].element < elements[i].element)) {
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
        std::vector<element_t> elements;

        // An array of positions for all elements.
        std::vector<node_t> position;
    };
}

#endif  // BATCHHLWW_H_
