//
// Created by Xinjie ZHOU on 25/8/2022.
//

#ifndef GTREE_H
#define GTREE_H

#include <iostream>
#include <fstream>
#include <sstream>

#include<stdio.h>
#include<metis.h>
#include<vector>
#include<stdlib.h>
#include<memory.h>
#include<unordered_map>
#include<map>
#include<unordered_set>
#include<set>
#include<deque>
#include<stack>
#include <list>
#include <queue>
#include<algorithm>
#include<sys/time.h>
#include <string>
#include <cassert>
#include <limits.h>
#include "Timer.h"
#include <boost/thread/thread.hpp>

using namespace std;
#define INF INT32_MAX
#define DECREASE 0
#define INCREASE 1
#define MIX 2

typedef int NodeId;
typedef unsigned int Distance;//long long

//string DataPath = "/Users/zhouxj/Documents/1-Research/1-Papers/0-My_Papers/2-EPSP/GTree-master/src/gtree/";
string DataPath = "/Users/zhouxj/Documents/1-Research/Datasets/";
//string DataPath = "/media/TraminerData/s4451682/GraphDataforPartition/Processed/";
int task_type = 1; //1: G*-Tree build and query; 2: G*-Tree query; 3: generate OD pairs. default: 1
//bool ifMac = true;
//bool ifMac = false;
int percent = 100;//percentage of shortcuts (comparing to edge number), default 100
int top_level = 1;
//bool ifDebug = false;//
string dirname;//intermediate file directory

// MACRO for timing
struct timeval tv;
long long ts, te;
#define TIME_TICK_START gettimeofday( &tv, NULL ); ts = tv.tv_sec * 100000 + tv.tv_usec / 10;
#define TIME_TICK_END gettimeofday( &tv, NULL ); te = tv.tv_sec * 100000 + tv.tv_usec / 10;
#define TIME_TICK_PRINT(T) printf("%s RESULT: %lld (0.01MS)\r\n", (#T), te - ts );
// ----------

// set all edge weight to 1 (unweighted graph)
#define ADJWEIGHT_SET_TO_ALL_ONE false //true
// we assume edge weight is integer, thus (input edge) * WEIGHT_INFLATE_FACTOR = (our edge weight)
#define WEIGHT_INFLATE_FACTOR 100000
// gtree fanout
#define PARTITION_PART 4
// gtree leaf node capacity = tau(in paper)
//#define LEAF_CAP 64 // maximum vertex number in leaf node
int LEAF_CAP = 64;
//#define THREAD_LIMIT 120//40 150

//const char *FILE_NODE = "cal.cnode";//const
//const char *FILE_EDGE = "cal.cedge";
//// gtree index disk storage
//const char *FILE_NODES_GTREE_PATH = "cal.paths.bin";
//const char *FILE_GTREE = "cal.gtree.bin";
//const char *FILE_ONTREE_MIND = "cal.minds.bin";
//// input
//const char *FILE_OBJECT = "cal.object";
//const char *FILE_SHORTCUT = "cal.shortcuts.bin";
//const char *FILE_QUERY = "cal.query";

typedef struct{
    double x,y;
    vector<int> adjnodes;   //adjacent vertex
    vector<int> adjweight;  //edge weight
    bool isborder;          //if border
    vector<int> gtreepath; // this is used to do sub-graph locating
    /// for G*-Tree
    int deep;//the number of levels from its leaf node to root node
    int inleaf;//the leaf node that it locates at
    int inleafpos;//the location position in leaf node
}Node;
// struct for tree node
typedef struct{
    vector<int> borders;    //border vertex
    vector<int> children;   //id of child nodes
    bool isleaf;            //if leaf node
    vector<int> leafnodes;  //vertex in the leaf node
    int father;             //id of father node
// ----- min dis -----
    vector<int> union_borders; // for non leaf node, it contains all borders of children; for leaf node, it contains all vertices of this node
    vector<int> mind; // min dis, row by row of union_borders
// ----- for pre query init, OCCURENCE LIST in paper -----
    vector<int> nonleafinvlist;
    vector<int> leafinvlist;
    vector<int> up_pos;         // the position IDs of the borders in parent node
    vector<int> current_pos;    // the position IDs of the borders of this tree node
/// for G*-Tree
// ----- for caching distances
    vector<vector<int>> cache;//store the distances from the borders of lca node to the borders of this node, dimension 1 is lca node
    bool is_cached;
// ----- for knn and range query-----
    unordered_set<int> oclist;
    bool is_visited;
// ----- for caching distances
    vector<int> cache_q;
}TreeNode;

// init status struct
typedef struct{
    int tnid; // tree node id
    set<int> nset; // node set
}Status;

class Gtree {
public:
    string dataset; //name of dataset
    unsigned int node_num = 0;      //the number of vertices
    unsigned int edge_num = 0;      //the number of edges
    int thread_num=10;
    int parti_num = 0;
    int noe = -1; // number of edges
    vector<Node> Nodes; // vector of vertex
    vector<TreeNode> GTree; //vector of tree node
    vector<TreeNode> GTreeO;
    vector<Node> NodesO;
    string FILE_NODE = ".cnode";//const
    string FILE_EDGE = ".cedge";
// gtree index disk storage
    string FILE_NODES_GTREE_PATH = ".paths.bin";
    string FILE_GTREE = ".gtree.bin";
    string FILE_ONTREE_MIND = ".minds.bin";
// input
    string FILE_OBJECT = ".object";
    string FILE_SHORTCUT = ".shortcuts";
    string FILE_QUERY = ".query";
    string FILE_UPDATE = ".update";
    vector<pair<pair<int,int>,int>> updateEdges;

//    // use for metis
////idx_t = int64_t / real_t = double
//    idx_t nvtxs; // |vertices|
//    idx_t ncon; // number of weight per vertex
//    idx_t *xadj; // array of adjacency of indices
//    idx_t *adjncy; // array of adjacency nodes
//    idx_t *vwgt; // array of weight of nodes
//    idx_t *adjwgt; // array of weight of edges in adjncy
//    idx_t nparts; // number of parts to partition
//    idx_t objval; // edge cut for partitioning solution
//    idx_t *part; // array of partition vector
//    idx_t options[METIS_NOPTIONS]; // option array

    void PathInit(bool ifMac);
    void options_setting();
    void load_graph();
    void ReadGraph_W();
    void data_transform_init( set<int> &nset );
    void init();
    void finalize();
    unordered_map<int,int> graph_partition( set<int> &nset );
    void build();
    void gtree_save();
    void load_gtreeQ();
    void build_up_and_down_pos();
    vector<int> dijkstra_candidate( int s, vector<int> &cands, vector<Node> &graph );

    void hierarchy_shortest_path_compute(int i, vector<int> & tn_v, vector<Node> & graph, vector<vector<unordered_map<int, unordered_map<int,int> >>> & vertex_pairsVV, int thread_i);
    void hierarchy_shortest_path_calculation(bool ifParallel);
    void hierarchy_shortest_path_compute_update(int i, vector<int> tn_v, vector<Node> & graph, vector<vector<unordered_map<int, unordered_map<int,int> >>> & vertex_pairsVV, int thread_i);
    void hierarchy_shortest_path_calculate(bool ifParallel);
    void hierarchy_shortest_path_save();
    void load_minds();
    int gtree_build(bool ifParallel); // build G-Tree
    void ODpairGenerate(int times);    //used to generate random OD pairs
    void UpdateGenerate(int times); //used to generate update pairs
    Distance Dijkstra(NodeId s, NodeId t, vector<Node> & Nodes);
    void ReadUpdates(string filename);
    bool CheckIfAncestor(int ID2, int nid);
    void init_update();
    void init_borders();
    void init_query(vector<TreeNode> & GTree);
    inline int find_LCA_pos(int src, int dst);
    bool dijkstra_candidate_update( int ID1, int ID1_pos, vector<int> &cands, vector<Node> &graph);//Dijkstra for leaf node update
    bool dijkstra_candidate_update( int ID1, int ID1_pos, int nid, vector<int> &cands, vector<Node> &graph, unordered_set<int> & borderSet);//Dijkstra for non-leaf node update
    void LeafNodeContract(int lnID, vector<Node> & graph, unordered_set<int> & leafsToCheck, unordered_set<int> & updatedLeafs);//unordered_set<int> & updatedLeafs,
    void NonLeafNodeContract(int lnID, int level_i, vector<Node> & graph, unordered_set<int> & nodesToCheck, unordered_set<int> & updatedNodes, vector< vector<int> > & treenodelevel, unordered_set<int> & borderSet);
    void GtreeUpdate(int updateType,int updateVolume, bool ifPrue);
    void LeafLevelUpdate(vector<Node> & graph, unordered_set<int> & leafsToCheck, unordered_set<int> & updatedNodes, unordered_set<int> & borderAffectedNodes, bool ifParallel);//update of the leaf level nodes
    void NonLeafLevelUpdate(vector<Node> & graph, int level_i, vector< vector<int> > & treenodelevel, unordered_set<int> & nodesToCheck, unordered_set<int> & updatedNodes, unordered_set<int> & borderAffectedNodes, bool ifParallel);//update of non-leaf level nodes
    void UpdateUpwardsPropagate(vector<Node> & graph, unordered_set<int> & updatedNodes, unordered_set<int> & borderAffectedNodes, bool ifParallel, int lca, unordered_set<int> & nodesToCheck);//propagate the update from current level to the highest affected level
    void GtreeUpdateParalel(int updateType,int updateVolume, int updateBatch, bool ifParallel);
    void UpdateDownwardsPropagate();//propagate the update from current level to the lowest affected level
    void LeafUpdate(int lnID, vector<Node> & graph, unordered_set<int> & leafsToCheck, unordered_set<int> & updatedLeafs,unordered_set<int> & borderAffectedNodes, unordered_map<int, bool> & flag_update);
    void NonLeafUpdate(int nid, vector<Node> & graph, unordered_map<int, bool> & flag_update);
    void UpdatePropagateToUpperLevels(vector<Node> & graph, unordered_set<int> & updatedLeafs, bool ifLeaf);
    void UpdatePropagateToUpperLevelsPrune(vector<Node> & graph, unordered_set<int> & updatedLeafs, unordered_set<int> & borderAffectedNodes);
    void UpdatePropagateToUpperLevelsPruneP(vector<Node> & graph, unordered_set<int> & updatedLeafs, unordered_set<int> & borderAffectedNodes);
    void GetChildNodes(int lca, unordered_set<int> & nodesToCheck, unordered_set<int> & leafsToCheck, unordered_set<int> & updatedLeafs);
    void GetChildLeafs(int lca, unordered_set<int> & leafsToCheck, unordered_set<int> & updatedLeafs);
    int ComputeDisByTree(int src, int dst, vector<TreeNode> & GTree);
    int dijkstra_p2p(int s, int t);
    int Distance_query(int src, int dst);
    void CorrectnessCheck(int times);
};

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
#endif //GTREE_H
