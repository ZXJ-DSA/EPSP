/*
 * head.h
 *
 *  Created on: 13 Oct 2022
 *      Author: zhangmengxuan, Xinjie ZHOU
 */

#ifndef HEAD_H_
#define HEAD_H_

#include <stdio.h>
#include <math.h>
#include <vector>
#include <map>
#include <set>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <iostream>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <stack>
#include <boost/thread/thread.hpp>
#include <chrono>
#include <string>
#include "Heap.h"
#include "labeling.hpp"
#include <omp.h>

#define INF 99999999
#define PreBoundary 1
#define NoBoundary 2
#define PostBoundary 3
//typedef unsigned int vertex;
typedef int vertex;

using namespace std;
using namespace boost;


extern vector<int> NodeOrder_;//nodeID order
extern vector<int> _DD_;

struct Nei{
	int nid;
	int w;
	int c;
};

struct OrderComp{
    int ID;
    OrderComp(){ID=0;}
    OrderComp(int _ID){
        ID=_ID;
    }
    bool operator< (const OrderComp d) const{
        if(NodeOrder_[ID]!=NodeOrder_[d.ID])
            return NodeOrder_[ID]>NodeOrder_[d.ID];
        return ID>d.ID;
    }
};

struct OrderComp2{
    int x;
    int y;//order(x)<order(y)
    OrderComp2(int _x, int _y){
        x=_x; y=_y;
    }
    bool operator< (const OrderComp2& d) const{
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

//struct OrderComp3{
//    int ID1, ID2;
//    OrderComp3(){
//        ID1=0, ID2=0;
//    }
//    OrderComp3(int _ID1, int _ID2){
//        ID1=_ID1;
//        ID2=_ID2;
//    }
//    bool operator< (const OrderComp3 d) const{//return the larger order vertex
//        return (NodeOrder_[ID1] > NodeOrder_[d.ID1]) || ((NodeOrder_[ID1] <= NodeOrder_[d.ID1]) && (NodeOrder_[ID2] > NodeOrder_[d.ID2]) );
//    }
//};

//return smallest Dis with largest-order vertex
struct MinComp{
    int s, t, Dis;//s is the source vertex while t is the target vertex
    MinComp(){
        s=0, t=0, Dis=0;
    }
    MinComp(int _ID1, int _ID2, int _Dis){
        s=_ID1; t=_ID2; Dis=_Dis;
    }
    bool operator< (const MinComp d) const{
        if(Dis != d.Dis){
            return Dis<d.Dis;
        }else{
            if(NodeOrder_[s]!=NodeOrder_[d.s]){
                return NodeOrder_[s]>NodeOrder_[d.s];
            }
            return NodeOrder_[t]>NodeOrder_[d.t];
        }
    }
};


struct DegComp1{//min-first
    int x;
    DegComp1(int _x){
        x=_x;
    }
    bool operator < (const DegComp1 d) const{
        if(_DD_[x]!=_DD_[d.x])
            return _DD_[x]<_DD_[d.x];
        return x<d.x;
    }
//    bool operator () (const DegComp1& a, const DegComp1& b) const{
//        return (_DD_[a.x] < _DD_[b.x]) || (_DD_[b.x] >= _DD_[a.x] && (a.x < b.x));
//    }

};

struct Node{//tree node
	vector<pair<int,pair<int,int>>> vert;//neighID/weight/count(how many ways can lead to this super edge weight)
//	vector<pair<int,Nei>> neighInf;//posID,neighID,weight,count(for shortcut information maintenance)
	vector<int> pos;
	vector<int> dis, cnt;//the distance value and corresponding count number (i.e., record how many path has the shortest distance)
    vector<int> vAncestor;//the ancestors, which is corresponding to dis
//    vector<int> disInfV;//the distances from this vertex to the interface vertices, for query efficiency test
    map<int,int> disInf;//the distances from this vertex to the interface vertices
	//vector<set<int>> FromNode;
//	set<int> changedPos;
	vector<bool> FN;//another succinct way of FromNode, whether this distance label is directly obtained from shortcuts (vert)
    map<int,bool> FNInf;//whether the interface distance is obtained from shortcuts (vert)
	set<int> DisRe;//record the vertex id that the distance label should be updated
    set<int> DisReInf;//record the vertex id that the interface label should be updated
	vector<int> ch;
	int height, hdepth;//hdepty is the deepest node that a vertex still exists
	int pa;//parent, the pa of root vertex is 0
	int uniqueVertex;//?vertex id of this tree node?
//	vector<int> piv;//pivot vetex, used in path retrieval
    int treeroot;//the tree id of subtree root, i.e., rank[x]
	Node(){
		vert.clear();
//		neighInf.clear();
		pos.clear();
		dis.clear();
//        disInfV.clear();
        disInf.clear();
		cnt.clear();
        vAncestor.clear();
		ch.clear();
		pa = -1;
		uniqueVertex = -1;
		height = 0;
		hdepth = 0;
//		changedPos.clear();
		FN.clear();
        FNInf.clear();
		DisRe.clear();
        DisReInf.clear();
//		piv.clear();
        treeroot=-1;
	}
};

class Graph{
public:
    string sourcePath;
    string dataset;
	int node_num;    //vertex number
	int edge_num;    //edge number
	vector<vector<pair<vertex,int>>> Neighbor;//original graph
    int partiNum;   //partition number
    int nodeNumCore;    //vertex number of core
	vector<int> DD; //intermediate variable in Contraction, DD2
	int threadnum=15;  //thread number
    int batchsize=15;
//    int indexType=0; //core index type
    int strategy = 2;//PSP strategy, 1: Pre-boundary; 2: No-boundary; 3: Post-boundary
//    bool ifOpt=false; // whether to use query-orient optimization
//    bool ifExtension=false; // whether to use extension optimization
    int algoCoreU=0;//algorithm for core update, (0: PDPLL; 1: SDPLL), default: 0
    int algoCoreC=0;//algorithm for core construction, (0: BPCL; 1: PCL; 2: PLL; 3: WPSL; 4: GLL; 5: Read)
    int algoTree=1;//algorithm for tree, (0: CH; 1: TD), default: 1

	//core-periphery decomposition
	int bandWidth;
	vector<int> NodeOrder;//nodeID order
	vector<int> vNodeOrder;//order nodeID
	vector<bool> existCore;
	int HighestOrder;
	vector<int> EulerSeqCore;
	vector<int> toRMQCore;
	vector<vector<int>> RMQIndexCore;

    bool ifVectorQ=false;//whether to use vector-based labels for querying
    vector<unordered_map<vertex,int>> Label;
    hl::Labeling LabelV;
    hl::PPR PPRV;
    vector<unordered_map<vertex,unordered_set<vertex>>> PruningPointSet;//{(v,c),u}
    vector<unordered_map<vertex,vertex>> PruningPointSet2;//<v,u,{c}>
    //	vector<vector<pair<int,int>>> PruningPointVector;//{v,(c,u)}
    //    vector<unordered_map<int,set<OrderComp>>> PruningPointSetOrder;//{(v,c),u}

    vector<unordered_set<int>> ChangedLabels;
    set<pair<int,int>> NoSupportedPair;


	//****************boundary vertex in each partition****************//
	vector<vector<int>> BoundVertex;//the interface vertex + root vertex
	vector<set<int>> BoundVertexSet;//only interface vertex, excluding root vertex
	//****************vertex tag for query answering & vertex order****************//
	vector<int> CoreTag;//-1 indicates core vertex and root vertex without children, i>=0 indicates non-core vertex (i.e., the root and inner-partition vertex) and which partition it belongs to
    vector<pair<bool,set<int>>> BoundTag;//first element: 1 (true) indicates boundary vertex (i.e., interface vertex), 0 (false) indicates non-boundary vertex; second element contains the partitions it belongs to
//    vector<bool> BoundTag;//1 (true) indicates boundary vertex (i.e., interface vertex), 0 (false) indicates non-boundary vertex
	//****************adjacent list for partitions and core****************//
	vector<vector<pair<vertex,int>>> AdjaCore;
	vector<map<int,int>> AdjaCoreMap;
    vector<map<int,int>> AdjaCoreMapOld;
//    unordered_map<int,unordered_map<int,int>> changedInterfaceEdges;
    set<int> CoreVertex;
	//****************Supported partition for boundary pair****************//
	vector<map<int,map<int,int>>> SuppPartiID;//<ID1,<ID2,<pid,weight>>>, record the partition and its supportive value for a given interface edge
	vector<map<int,pair<int,set<int>>>> SuppPartiIDReal;//ID1,<ID2,<weight,set<pid>>>>, //record the partitions that really support a given interface edge


    /// Index Construction
//    vector<omp_lock_t> oml;
    unordered_map<int, Semaphore*> mSm;
    vector<Semaphore*> vSm;
//	vector<vector<unordered_map<int,int>>> LabelStep;
    vector<unordered_map<vertex,int>> Dhop;
    vector<vector<pair<vertex,int>>> DhopV;
    vector<int> Dvectex;
    vector<bool> DvertexNew;
//    vector<pair<int,int>> RedundantLabels;//use to record the redundant labels caused by PSL
    vector<unordered_set<int>> RedundantLabels;//use to record the redundant labels caused by PSL

    //H2H index construction
    //intermediate variable and function used in the H2H index construction
    vector<vector<pair<int,pair<int,int>>>> NeighborCon;//ID1,ID2,(weight,count)
//    vector<map<int, vector<int>>> SCconNodesMT;//supportive vertex, multiple thread of SCconNodes
    vector<map<int, vector<pair<int,int>>>> SCconNodesMT;//<ID1,<ID2,<x,weight>>> supportive vertex, x is the contracted vertex
    vector<map<int,pair<int,int>>> E;//ID1,ID2,(weight,count)
    vector<vector<int>> VidtoTNid;//record the child tree nodes whose vert neighbors contain this tree node (nodeID--->tree node rank)
    vector<int> rank;//rank[v]>0 indicates non-core vertex
    int heightMax;
    vector<Node> Tree;
    vector<int> EulerSeq;//prepare for the LCA calculation, EulerSeq is the Euler tour, i.e., E[1,...,2n-1]
    vector<int> toRMQ;//toRMQ[p] records the first occurrence of node p in the Euler tour, i.e., R[p]
    vector<vector<int>> RMQIndex;//?

    vector<int> tRoots;//root vertex id of trees

    ~Graph(){
        clear();
    }
    void clear(){
        Neighbor.clear();
        AdjaCore.clear();
        AdjaCoreMap.clear();
        AdjaCoreMapOld.clear();
        Label.clear();
        Tree.clear();
        LabelV.clear();
        PPRV.clear();
        vSm.clear();

    }

    void IndexConstruction();
    void CTIndexConstruct(); //Core-tree index construction, tree index is TD
    void CTIndexConstructCH(); //Core-tree index construction, tree index is CH
    void MDEContract();
    void Create_partitions();//obtain the partitions for TD-CH
    void Create_tree();
    void Compute_tree_label(bool ifParallel);
    void Construct_tree(bool ifParallel);
    double Construct_core(int indexType);//0 indicates using PLL, 1 indicates using PSL

    int BoundaryDijkstra(int ID1, vector<int> & IDs, vector<int>& Dis, vector<vector<pair<vertex,int>>> &Neighbor);//For pre-boundary strategy

    void makeRMQCore();
    void makeRMQDFSCore(int p, int height);
    void TreeLabelCompute(pair<int,int> pidRange, vector<int> & pidRanks);
    void TreeLabelCompute2(pair<int,int> pidRange, vector<int> & pidRanks);
    void makeTreeIndexDFS(int p, vector<int> &ancestor, vector<int> & interface);
    void makeTreeIndexDFS2(int p, vector<int> &ancestor, vector<int> & interface, map<int,unordered_map<int,int>> & disInfs);
//    void makeIndexDFS(int p, vector<int> &list);
    void NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x);
    void insertEMTorder(int u,int v,int w);
    void deleteECore(int u,int v);
    void insertECore(int u,int v,int w);
    int matchCore(int x,vector<pair<int,pair<int,int>>> &vert);//vector<pair<int,int>> &vert
    int matchCoreCH(int x,vector<pair<int,pair<int,int>>> &vert);//vector<pair<int,int>> &vert
    void IndexsizeCTH2H();  //Core-tree index size computation
    void IndexsizeCTCH();  //Core-tree index size computation

	//PLL index construction
	void PLLIndexConstruct();
	void DijksPrune1(int nodeID, vector<pair<int,int>>& vp);

    double PLLConstructV(vector<vector<pair<vertex,int>>>& Neighbor);
    void DijksPrune1V(int nodeID,vector<vector<pair<vertex,int>>>& Neighbor);
    int PLLDisQuery1V(int ID1,unordered_map<vertex,int>& Lh,vector<int>& SupNode);


    //GLL
    void LCCIndexConstruct(vector<vector<pair<vertex,int>>> &Neighbor);
    void DijksPrune2(int nodeID, vector<vector<pair<vertex,int>>> &Neighbor);
    void DQ_Clean(vertex ID1);
    void DQ_Cleans(vector<vertex>& IDs);
    void GLLIndexConstruct(vector<vector<pair<vertex,int>>> &Neighbor);
    void CommitLabel(int ID, vector<pair<int,int>>& vp);
    void DijksPrune3(int nodeID, vector<pair<int,int>>& vp, vector<vector<pair<vertex,int>>> &Neighbor);

    double GLLConstructV(vector<vector<pair<vertex,int>>> &Neighbor);
    void DijksPrune4V(vertex nodeID, vector<vector<pair<vertex,int>>> &Neighbor);
    int ShortestDisQuery1V(int ID1,unordered_map<vertex,int>& Lh,vector<int>& SupNode, int& d);
    
    void BatchClean(vector<vertex> &bNodes, unordered_set<vertex>& setNodes, bool ifParallel);
    void LabelClean(vertex ID1, unordered_set<vertex>& setNodes);

    //PCL
    void PCLIndexConstruct(vector<vector<pair<vertex,int>>> &Neighbor);
    void PCLDijk(vertex nodeID, vector<vector<pair<vertex,int>>> &Neighbor);
    void PCLDijks(vector<vertex> & IDs, vector<vector<pair<vertex,int>>> &Neighbor);
    double PCLConstructV(vector<vector<pair<vertex,int>>> &Neighbor);
    void PCLDijkV(vertex nodeID, vector<vector<pair<vertex,int>>> &Neighbor);
    void PCLDijkV2(vector<vertex>& ProcessID, vector<vector<pair<vertex,int>>> &Neighbor);

    //BPCL
    void BPCLIndexConstruct(vector<vector<pair<vertex,int>>> &Neighbor);
    void BatchPCLDijk(vertex nodeID, unordered_set<vertex>& setNodes, vertex hID, vector<vector<pair<vertex,int>>> &Neighbor);
    void BatchPCLDijk2(vector<vertex>& p, unordered_set<vertex>& setNodes, vertex hID, vector<vector<pair<vertex,int>>> &Neighbor);
    double BPCLConstructV(vector<vector<pair<vertex,int>>> &Neighbor);
    void BatchPCLDijkV(vertex nodeID, unordered_set<vertex>& setNodes, vertex hID, vector<vector<pair<vertex,int>>> &Neighbor);
    void BatchPCLDijk2V(vector<vertex>& ProcessID,  unordered_set<vertex>& setNodes, vertex hID, vector<vector<pair<vertex,int>>> &Neighbor);

    void PCLDijkPPR(int nodeID);
    void ThreadDistribute(vector<vertex>& vertices, vector<vector<vertex>>& processID);
    void PPRConstruction(vertex nodeID, vector<vector<pair<vertex,int>>> &Neighbor);
    void PPRConstruction2(vector<vertex> & p, vector<vector<pair<vertex,int>>> &Neighbor);
    void PruningPointBuild(bool ifParallel, vector<vector<vertex>> & processID, vector<vector<pair<vertex,int>>> &Neighbor);
    void PPRConstructionV(vertex nodeID, vector<vector<pair<vertex,int>>> &Neighbor);
    void PPRConstructionV2(vector<vertex>& ProcessID, vector<vector<pair<vertex,int>>> &Neighbor);
    int ShortestDisQueryPeakV(int ID1,int ID2, unordered_map<vertex,int>& Lh, vector<int>& SupNode, int& d);

	//WPSL index construction
	void PSLIndexConstruct(vector<vector<pair<vertex,int>>> &Neighbor);
	bool DhopLableRefreshMulti2New(int step, vector<vector<pair<vertex,int>>> &Neighbor);
	void labelMultiThread2New(vector<unordered_map<vertex,int>>& newDhop, vector<int>& p,int step, vector<vector<pair<vertex,int>>> &Neighbor);
	void threadDistribute(vector<vector<int>>& processID);
	int ShortestDisQuery1(int ID1,int ID2,vector<int>& SupNode, int& d);
    int ShortestDisQuery2(int ID1,int ID2,vector<int>& SupNode, int& d);

    double PSLConstructV(vector<vector<pair<vertex,int>>> &Neighbor);
    bool DhopLableRefreshStepV(int step, vector<vector<pair<vertex,int>>> &Neighbor);
    void labelMultiThread2NewV(vector<vector<pair<vertex,int>>>& newDhop, int nodeID, vector<vector<pair<vertex,int>>> &Neighbor);
    void labelMultiThread2NewV2(vector<vector<pair<vertex,int>>>& newDhop, vector<int>& p, vector<vector<pair<vertex,int>>> &Neighbor);


    void CompareLabel(string filename1, string filename2);//compare two label, using filename1 as benchmark
    void CleanLabel(vector<unordered_map<vertex,int>> &Label);//clean all redundant labels
    void CleanLabel(vector<unordered_map<vertex,int>> &Label, string filename);//clean all redundant labels
    void CleanPruningPoints(vector<vector<pair<vertex,int>>> &Neighbor,vector<unordered_map<vertex,int>> &Label, vector<unordered_map<int,set<OrderComp>>> &PruningPointNew);//clean pruning points

    double timeCore=0;

	///Query processing
    void EffiCheck(string filename,int runtimes);
    void EffiCheckQueryTypeTest(int runtimes);
    void EffiCheckQueryType(string filename,int runtimes, int type);
	int Query(int ID1, int ID2);
    int QueryDebug(int ID1, int ID2);
	int QueryPartiCore(int ID1, int ID2);
    int QueryPartiCoreDebug(int ID1, int ID2);
    int QuerySameParti(int ID1, int ID2);
    int QuerySameParti2(int ID1, int ID2);
	int QueryPartiParti(int ID1, int ID2);
	int QueryCore(int ID1, int ID2);
    int QueryCoreDebug(int ID1, int ID2);
    int QueryH2H(int ID1,int ID2);//shortest distance query with no partition
    int LCAQuery(int _p, int _q);
    int QueryPeripheryTree(int ID1, int ID2, int PID);//
    int QueryPeripheryTree2(int ID1, int ID2, int PID);//
    int LCAQueryPartition(int _p, int _q, int PID);// query within partition

    int Query_CH(int ID1, int ID2);
    int QueryPartiCore_CH(int ID1, int ID2);
    int QueryPartiParti_CH(int ID1, int ID2);
    int QuerySameParti_CH(int ID1, int ID2);
    void QueryPeriphery_CH(int ID1, set<int> todo, map<int,int> & results);//single-source CH query
    int QueryPeriphery_CH(int ID1, int ID2);

    //Correctness Check
    void CorrectnessCheck(int runtimes);
    void CorrectnessCheckCore();
    void CorrectnessCheckCoreAll();
    void CoreGraphDebug(string graphfile);//core index debug
    void SameTreeQueryTest(string filename,int runtimes);
    void SameTreeQueryGen(vector<int>& tRoots, int times);
    void DFSTree(vector<int>& tNodes, int id);
    void SameTreeUpdateGen(vector<int>& tRoots, int times);

    //Dijkstra
    int Dijkstra(int ID1, int ID2,vector<vector<pair<vertex,int>>> &Neighbor);
    void RetrievePath(int ID1, int ID2, vector<int> & prece);
    int DijkstraCore(int ID1, int ID2);
    pair<int,int> DijkstraCoreDebug(int ID1, int ID2);

    /// Extension query
    bool ifParallel = true;
//    bool extUpdate = false;
    vector<unordered_set<int>> PartiVertex;//only insert in-partition vertex (i.e., periphery vertex, excluding the interface vertex) in this case
//    vector<unordered_map<int,int>> IndexExt;//extended 2-hop label of vertex in periphery
//    vector<bool> PartiUpdateExt;//indicates whether the periphery's extended index should be updated
//    void ExtensionIndex(pair<int,int> pidRange, bool ifIncrease);
//    void ExtensionIndex2(pair<int,int> pidRange, bool ifIncrease, vector<int>& partiForUpdate);
//    void ExtensionIndexConstruct(bool ifParallel, bool ifIncrease);
//    void ExtensionIndexUpdate(bool ifParallel, bool ifIncrease, vector<int>& partiForUpdate);
//    int QueryPartiCoreExt(int ID1, int ID2);
//    int QueryPartiCoreExtDebug(int ID1, int ID2);
//    int QueryPartiPartiExt(int ID1, int ID2);
//    int QueryPartiPartiExtDebug(int ID1, int ID2);

	/// Index update
    void IndexMaintenance(string updateFile, int updateType, int updateBatch);
    void IndexMaintenanceTypeTest(int updateType, int updateBatch);
    void IndexMaintenanceType(string updateFile, int updateType, int updateBatch, int type);
	void Decrease(int a, int b, int oldW, int newW);
	void Increase(int a, int b, int oldW, int newW);

    void QueryGenerationSameParti();
    void ODGeneSameParti(int num, string filename);

    void DecreasePSL(int a, int b, int oldW, int newW, vector<vector<pair<vertex,int>>> &Neighbors,vector<unordered_map<vertex,int>> &Label);
    void IncreasePSLNew(int a, int b, int oldW, int newW, vector<vector<pair<vertex,int>>> &Neighbor,vector<unordered_map<vertex,int>> &Label, vector<unordered_map<vertex,unordered_set<vertex>>> &PruningPointNew);//set version with NoSupportedPair
    void IncreasePSL(int a, int b, int oldW, int newW, vector<vector<pair<vertex,int>>> &Neighbors,vector<unordered_map<vertex,int>> &Label,vector<unordered_map<int,vector<int>>> &PruningPointNew,set<pair<int,int>> &NoSupportedPair);//original version, vector with NoSupportedPair

    // Core Update
    void CoarseUpdate(int LID, int HID, int oldW, queue<pair<int,pair<int,int>>>& WaitPro, queue<pair<int,pair<int,int>>>& WaitProP, vector<pair<int,int>>& AL1, vector<pair<int,int>>& AL2, vector<pair<int,int>>& AL2Check, vector<vector<pair<vertex,int>>> &Neighbor,vector<unordered_map<vertex,int>> &Label,bool ifDebug, int lid, int hid);//queue version
    void RefineUpdate(queue<pair<int,pair<int,int>>>& WaitPro, queue<pair<int,pair<int,int>>>& WaitProP, vector<pair<int,int>>& AL1, vector<pair<int,int>>& AL2, vector<pair<int,int>>& AL2Check, set<tuple<vertex,vertex,vertex>>& outdatedPruning, map<pair<vertex,vertex>,vertex>& newPruningPoints, vector<vector<pair<vertex,int>>> &Neighbor,vector<unordered_map<vertex,int>> &Label,vector<unordered_map<vertex,unordered_set<vertex>>> &PruningPointNew, bool ifDebug, int lid, int hid);//queue version
    bool PPRCheck(int curID, int hubID, vector<vector<pair<vertex,int>>> &Neighbor,vector<unordered_map<vertex,int>> &Label,vector<unordered_map<vertex,unordered_set<vertex>>> &PruningPointNew, queue<pair<int,pair<int,int>>>& WaitProPTem, vector<pair<int,int>> &AL2, map<pair<vertex,vertex>,vertex> &newPruningPoints, set<tuple<vertex,vertex,vertex>> &outdatedPruning, bool ifDebug, int lid, int hid);//queue version
    void PPRClean(vector<vector<pair<vertex,int>>> &Neighbors, map<pair<vertex,vertex>,vertex> &newPruningPoints, set<tuple<vertex,vertex,vertex>> &outdatedPruning, bool ifDebug, int lid, int hid);

    int PLLQuery(int ID1,int ID2,vector<unordered_map<vertex,int>> &Label);
    int DisQueryVallyNew(int ID1, int ID2, vector<vector<pair<vertex,int>>> &Neighbors,vector<unordered_map<vertex,int>> &Label);
    int DisQueryVally(int ID1, int ID2, vector<vector<pair<vertex,int>>> &Neighbors,vector<unordered_map<vertex,int>> &Label);
    int DisQueryVallyVert(int ID1, int ID2, vector<vector<pair<vertex,int>>> &Neighbors,vector<unordered_map<vertex,int>> &Label);
    pair<int,int> DisQueryVallyVert2(int ID1, int ID2, vector<vector<pair<vertex,int>>> &Neighbors,vector<unordered_map<vertex,int>> &Label);
    pair<int,int> DisQueryVally2(int ID1, int ID2, vector<vector<pair<vertex,int>>> &Neighbors,vector<unordered_map<vertex,int>> &Label);
    int DisQueryVallyDebug(int ID1, int ID2, vector<vector<pair<vertex,int>>> &Neighbors,vector<unordered_map<vertex,int>> &Label);
    pair<int,int> DisQueryPeak2(int ID1, int ID2,vector<unordered_map<vertex,int>> &Label);
    int DisQueryPeak(int ID1, int ID2,vector<unordered_map<vertex,int>> &Label);
    int DisQueryLower1(int ID1, int ID2, vector<vector<pair<vertex,int>>> &Neighbors,vector<unordered_map<vertex,int>> &Label);

    //Tree update
    void DecreaseH2HNew(int a,int b, int newW, vector<vector<pair<vertex,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifParallel);//either a or b is periphery vertex
    void IncreaseH2HNew(int a, int b, int oldW, int newW, vector<vector<pair<vertex,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, bool ifParallel);
    void InterfacePropagate(int child, vector<int>& interfaces, vector<Node> &Tree, bool ifIncrease);
    void InterfacePropagateParallel(pair<int,int> pRange, vector<int>& pids, bool ifIncrease);
    void AncestorEntryDecreaseUpdate(int child,vector<int>& line, vector<int>& interfaces, set<int>& vertexIDChL, map<int,int>& checkedDis, vector<Node> &Tree, vector<int> &rank);
    void AncestorEntryIncreaseUpdate(int children, vector<int>& line, vector<int>& interfaces, int& changelabel, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid, int lowestH);

    void CHdecreaseBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch);
    void CHincreaseBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch);
    void DecreaseCHNew(int a,int b, int newW, vector<vector<pair<vertex,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax);
    void IncreaseCHNew(int a, int b, int oldW, int newW, vector<vector<pair<vertex,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid);

    ///SDPLL
    void MDEOrdering();
    void DegreeOrdering();
    void BPCLConstruction();

    void PLLdec(int a,int b, int oldW, int newW, vector<vector<pair<vertex,int>>> &Neighbor);
    void PLLinc(int u,int v, int oldW, int newW, vector<vector<pair<vertex,int>>> &Neighbor);
    void DijkstraPrune1(int NodeID, int ID, int d, vector<vector<pair<vertex,int>>> &Neighbor);
    int DijkstraList(int ID1, vector<int>& distance, vector<vector<pair<vertex,int>>> &Neighbor);
    void PLLaffect(int u, int oldw, vector<int>& PAu, vector<int>& RAu, vector<int> disu, vector<int> disv, vector<int> disvNew, vector<vector<pair<vertex,int>>> &Neighbor);
    void PLLinvalid(int u, vector<int>& ILu, vector<int>& MLu, vector<int> RAu, vector<int> RAv, vector<int> PAv);
    void PLLremove(int u, vector<int> RAu, vector<int> ILu);
    void PLLPrune(int r, set<int> Affect, vector<vector<pair<vertex,int>>> &Neighbor);
    int ShortestDisQuery(int ID1,int ID2);
    int PrefixalDisQuery(int ID1, int ID2);

	/// Graph Preprocessing
    void ReadGraph(string filename);
    void ReadUpdate(string filename,vector<pair<pair<int,int>,int>>& TestData);
    void ReadUpdate2(string filename,vector<pair<pair<int,int>,pair<int,int>>>& TestData);
    void ReadUpdate3(string filename,vector<pair<pair<int,int>,tuple<int,int,int>>>& TestData);
	void StainingMethod(int ID);
	void ODGene(int num, string filename);
	void UpdateGene(int num, string filename);
    void CTQueryUpdateGen(int num);

    void WriteCoreIndex(string file);
    void ReadCoreIndex(string file);
    void WriteTreeIndex(string filename);
    void ReadTreeIndex(string file);
    void WriteTreeIndex2(string filename);
    void ReadTreeIndex2(string file);
    void WriteGraph(string graphfile);
	void WriteCoreGraph(string graphfile);
	void ReadCoreGraph(string filename);
    void WriteOrder(string filename);
    void ReadOrder(string filename);
    void CompareOrder(string filename1, string filename2);


    vector<int> DFS_CC(vector<map<int,int>> & Edges, set<int> set_A, set<int> & set_B, int nodenum);

};

#endif // HEAD_H_
