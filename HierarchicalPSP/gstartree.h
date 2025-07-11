//
// Created by Xinjie ZHOU on 25/8/2022.
//

#ifndef GSTARTREE_H
#define GSTARTREE_H

#include <queue>
#include "Semaphore.h"
#include "gtree.h"

struct leaf_node_pair {
    int x, y;
    double v;

    leaf_node_pair(int _x, int _y, double _v) {
        x = _x;
        y = _y;
        v = _v;
    }

    bool operator<(const leaf_node_pair &rhs) const {
        return v < rhs.v;   // top item is largest
    }
};

vector<int> _DD_;//true degree, temporal degree ,_DD2_
struct DegComp1{
    int x;
    DegComp1(int _x){
        x=_x;
    }
    bool operator< (const DegComp1 d) const{
        if(_DD_[x]!=_DD_[d.x])
            return _DD_[x]<_DD_[d.x];
//        if(_DD2_[x]!=_DD2_[x])
//            return _DD2_[x]<_DD2_[d.x];
        return x<d.x;
    }
};


struct Nei{
    int nid;
    int w;
    int c;
};

struct TDNode{//tree node
    vector<pair<int,pair<int,int>>> vert;//neighID/weight/count
    vector<int> pos;
    vector<int> dis, cnt;//the distance value and corresponding count number
    vector<int> vAncestor;//the ancestors, which is corresponding to dis
    vector<bool> FN;//another succint way of FromNode
    set<int> DisRe;
    vector<int> ch;
    int height, hdepth;//hdepty is the deepest node that a vertex still exists
    int pa;//parent
    int uniqueVertex;
    vector<int> piv;//pivot vetex, used in path retrieval
    TDNode(){
        vert.clear();
        pos.clear();
        dis.clear();
        vAncestor.clear();
        cnt.clear();
        ch.clear();
        pa = -1;
        uniqueVertex = -1;
        height = 0;
        hdepth = 0;
//        changedPos.clear();
        FN.clear();
        DisRe.clear();
        piv.clear();
    }
};

vector<int> NodeOrder_;
struct OrderCompMin{//prior to return the vertex with smaller order
    int x;
    OrderCompMin(int _x){
        x=_x;
    }
    bool operator< (const OrderCompMin& d) const{
        if(x==d.x){//avoid the redundant
            return false;
        }else{
            if(x!=d.x)
                return NodeOrder_[x]<NodeOrder_[d.x];
        }
    }
};

struct OrderCompMax{//prior to return the vertex with higher order
    int x;
    OrderCompMax(int _x){
        x=_x;
    }
    bool operator< (const OrderCompMax& d) const{
        if(x==d.x){//avoid the redundant
            return false;
        }else{
            if(x!=d.x)
                return NodeOrder_[x]>NodeOrder_[d.x];
        }
    }
};

struct OrderComp{//the edge with lower order first
    int x;
    int y;//order(x)<order(y)
    OrderComp(int _x, int _y){
        x=_x; y=_y;
    }
    bool operator< (const OrderComp& d) const{
        if(x==d.x && y==d.y){//avoid the redundant
            return false;
        }else{
            if(x!=d.x)
                return NodeOrder_[x]<NodeOrder_[d.x];
            if(y!=d.y)
                return NodeOrder_[y]<NodeOrder_[d.y];
        }
    }
};

class Gstartree:public Gtree{
public:
    unsigned long SHORTCUT_THRESHOLD = 0;
    unordered_map<int, unordered_map<int, vector<int> > > shortcuts;
    vector<int> query_objects;
//    int leafNodePairNum=0;
    vector<pair<int,int>> leafNodePairs;

    void IndexConstruction();
    void IndexMaintenance(int updateType, int updateBatch, int updateVolume, bool ifRead);
    void GstartreeIndexUpdate(vector<pair<pair<int,int>,pair<int,int> > > & updates, bool ifParallel);

//    void load_gtree();
    void load_shortcuts(string filename);
//    void build_up_and_down_pos();
    double compute_value_of_leaf_node_pair(int i, int j);//inline
    bool check_leaf_adjacent(TreeNode &li, TreeNode &lj);
    pair<int,int> compute_lca_level(int x, int y);
    void build_shortcuts_for_nodes(int x, int y, int top_level_x,int top_level_y);//inline
    void build_shortcuts_for_nodes(int x, int y);//inline
    void build_shortcuts();
    void save_shortcuts(string FS);

    int gstartree_build(); // build G*-Tree


    void init_ReadIndex();
    vector<int> load_objects();
//    int dijkstra_p2p(int s, int t);
//    vector<int> dijkstra_candidate(int s, unordered_set<int> &cands);
//    vector<int> dijkstra_candidate(int s, vector<int> &cands);
    inline bool check_shortcut(int ns, int nt);
//    inline int find_LCA_pos(int src, int dst);
    inline int dist_query(int src, int dst);
    void EfficiencyTest(int run_times); // query based on G*-Tree
    void CorrectnessCheck(int times);
    void QueryProcessingTest(int runtimes);
    int Query_TGTree(int ID1, int ID2);
    int Query_TGTreeDebug(int ID1, int ID2);
    int QueryPartiBoundary(int ID1,int ID2);
    int QueryPartiBoundaryDebug(int ID1,int ID2);
    int QueryPartiParti(int ID1,int ID2);
    int QueryH2H(int ID1,int ID2);
    int QueryH2HDebug(int ID1,int ID2);

    /// For N-TS-HP
    vector<int> NodeOrder;
    vector<int> vNodeOrder;
    vector<vector<int>> boundaryLevel;//hierarchical tree of boundary vertex
//    vector<bool> BoundaryTag;
    vector<vector<int>> VidtoTNid;//one vertex exist in those tree nodes (nodeID--->tree node rank)
    vector<int> EulerSeq;//prepare for the LCA calculation
    vector<int> toRMQ;
    vector<vector<int>> RMQIndex;
    vector<map<int, vector<int>>> SCconNodesMT;
//    map<int,vector<pair<int,int>>> OverlayGraph;
    vector<unordered_map<int,int>> OverlayGraph;
    int overlayNodeNum=0;
//    map<int,unordered_map<int,vector<pair<int,int>>>> PartiGraph;
    vector<unordered_map<int,pair<int,int>>> E;//ID1,ID2,weight,count
    vector<vector<pair<int,pair<int,int>>>> NeighborCon;
    int treeWidth=0;

    vector<Semaphore*> vSmNode;
    vector<int> rank;
    vector<TDNode> Tree;
    int heightMax;
    bool ifHierarchicalOrdering=false;//whether to use hierarchical vertex ordering

    void TDGTreeIndexBuild();
//    void leafNodeMatrix_save();
    void H2HIndexSave();
    void graphInfoSave(string filename);
    void load_graphOrdering(string filename);
    void init_TGTreeIndex();
    void IndexSizeH2H();
    void getOverlayLayers();//get boundary layers
    void StratifiedMDEOrdering(vector<vector<int>> & boundaryLevel);
    void getOriginalOverlayGraph();
    void getPartiGraph();
    void hierarchy_shortest_path_calculation_NoBoundary(bool ifParallel);
    void hierarchy_shortest_path_calculation_PreBoundary(bool ifParallel);
    vector<int> dijkstra_candidate_No( int s, vector<int> &cands, vector<Node> &graph);
    void leafNodeDisMatrixCompute(pair<int,int> p, vector<int> & leafNodes, bool ifGlobal);
//    void leafNodeAllPair_No(pair<int,int> p, vector<int> & leafNodes, map<int,unordered_map<int,vector<pair<int,int>>>> & graphs);
    void H2HconOrderMT(bool ifMDE);
    void MDEContract();
    void VertexContraction();
    void makeTree();
    int match(int x,vector<pair<int,pair<int,int>>> &vert);
    void makeRMQDFS(int p, int height);
    void makeRMQ();
    int LCAQuery(int _p, int _q);
    void makeIndex();//make index
    void makeIndexDFS(int p, vector<int> &list);

    void deleteE(int u,int v);
    void insertE(int u,int v,int w);
    void insertEMTorder(int u,int v,int w);
    void NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x);

    void TGTreeIndexUpdate(vector<pair<pair<int,int>,pair<int,int> > > & updates, bool ifParallel, int updateType);
    void PartitionUpdate_No(int tn,map<pair<int,int>,pair<int,int>> & affectedBPairs, bool ifParallel);
    void boundaryVUpdate_No(pair<int,int> p, int tn, vector<int> & cands, map<pair<int,int>,pair<int,int>> & affectedBPairs);
    void H2HdecBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch);//decrease
    void H2HincBatMT(vector<pair<pair<int,int>,pair<int,int>>>& wBatch);//increase
    void EachNodeProBDis5(int child,vector<int>& line,set<int>& vertexIDChL, map<int,int>& checkedDis);
    void eachNodeProcessIncrease1(int children, vector<int>& line, int& changelabel);


    /// LG-Tree
    vector<unordered_map<int,int>> Label;
    vector<vector<int>> vertexLevels;//vertex levels
    vector<int> levelTag;//level tags of all vertices
    void LGTreeIndexBuild();
    void VertexHierarchyBuild();
    void BorderLabelBuild();
    void IndexSizePLL();
    int Query_LGTree(int ID1, int ID2);
    int Query_LGTreeDebug(int ID1, int ID2);
    int QueryPLL(int ID1, int ID2);
    int QueryPLLDebug(int ID1, int ID2);
    int QueryPartiBoundaryPLL(int ID1,int ID2);
    int QueryPartiPartiPLL(int ID1,int ID2);
    int QueryPartiPartiPLLDebug(int ID1,int ID2);

    void init_LGTreeIndex();
    void overlayIndexSave(string filename);
    void load_overlayIndex(string filename);

    void LGTreeIndexUpdate(vector<pair<pair<int,int>,pair<int,int> > > & updates, bool ifParallel, int updateType);
    void PartitionUpdate_Pre(int tn,map<pair<int,int>,pair<int,int>> & affectedBPairs, bool ifParallel);
    void boundaryVUpdate_Pre(pair<int,int> p, int tn, vector<int> & cands, map<pair<int,int>,pair<int,int>> & affectedBPairs);
    void BorderLabelUpdateDec(vector<pair<pair<int,int>,pair<int,int>>>& wBatch);//decrease
    void BorderLabelUpdateInc(vector<pair<pair<int,int>,pair<int,int>>>& wBatch);//increase
    void CHdecBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& updateEdges);
    void CHincBatMT(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& updateEdges);
    int DisQueryVally(int ID1, int ID2);
};

#endif //GSTARTREE_H
