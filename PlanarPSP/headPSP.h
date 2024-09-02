/*
 * head.h
 *
 *  Created on: 24 August 2022
 *      Author: Xinjie ZHOU
 */

#ifndef HEADPSP_H_
#define HEADPSP_H_

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
#include <sys/stat.h>


#define INF 99999999
#define Dijk 0
#define CH 1
#define H2H 2
#define PLL 3
#define PreBoundary 1
#define NoBoundary 2
#define PostBoundary 3
#define EdgeCut 1
#define VertexCut 2
#define EdgeDecrease 1
#define EdgeIncrease 2
#define EdgeInsertion 3
#define VertexInsertion 4
#define EdgeDeletion 5
#define VertexDeletion 6

//typedef unsigned int vertex;
typedef int vertex;

using namespace std;
//using namespace boost;


extern vector<int> NodeOrder_;//nodeID order
extern vector<int> _DD_;
extern vector<int> _DD2_;

struct Nei{
	int nid;
	int w;
	int c;
};

struct OrderCompMax{// maximum-first, Higher-order first
    int ID;
    OrderCompMax(){ID=0;}
    OrderCompMax(int _ID){
        ID=_ID;
    }
    bool operator< (const OrderCompMax d) const{
        if(NodeOrder_[ID]!=NodeOrder_[d.ID])
            return NodeOrder_[ID]>NodeOrder_[d.ID];
        return ID>d.ID;
    }
};

struct OrderCompMin{//prior to reture the vertex with smaller order
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

struct OrderCompp{//prior to return the vertex with smaller order
    int x;
    OrderCompp(int _x){
        x=_x;
    }
    bool operator< (const OrderCompp& d) const{
        if(x==d.x){//avoid the redundant
            return false;
        }else{
            if(x!=d.x)
                return NodeOrder_[x]<NodeOrder_[d.x];
        }
    }
};


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


struct DegComp{//min-first
    int x;
    DegComp(int _x){
        x=_x;
    }
    bool operator < (const DegComp d) const{
        if(_DD_[x]!=_DD_[d.x])
            return _DD_[x]<_DD_[d.x];
        return x<d.x;
    }
//    bool operator () (const DegComp1& a, const DegComp1& b) const{
//        return (_DD_[a.x] < _DD_[b.x]) || (_DD_[b.x] >= _DD_[a.x] && (a.x < b.x));
//    }

};

struct DegComp2{//min-first
    int x;
    DegComp2(int _x){
        x=_x;
    }
    bool operator < (const DegComp2 d) const{
        if(_DD2_[x]!=_DD2_[d.x])
            return _DD2_[x]<_DD2_[d.x];
        return x<d.x;
    }
//    bool operator () (const DegComp1& a, const DegComp1& b) const{
//        return (_DD_[a.x] < _DD_[b.x]) || (_DD_[b.x] >= _DD_[a.x] && (a.x < b.x));
//    }

};

struct Node{//tree node
	vector<pair<int,pair<int,int>>> vert;//neighID/weight/count(how many ways can lead to this super edge weight)
	vector<int> pos;//the position index of neighbor vert. last element is the position index of this node
	vector<int> dis, cnt;//the distance value and corresponding count number (i.e., record how many path has the shortest distance)
    vector<int> disPost, cntPost;//the distance value of post-boundary strategy and corresponding count number
    vector<int> disInf;//the distance from this vertex to the boundary vertices and corresponding count number, cntInf
    vector<int> vAncestor;//the ancestors, which is corresponding to dis
    int rootPos;//index position of the root vertex of this partition in vAncestor
	vector<bool> FN;//another succinct way of FromNode, whether this distance label is directly obtained from shortcuts (vert)
    vector<bool> FNPost;//another succinct way of FromNode, whether this distance label is directly obtained from shortcuts (vert)
    vector<bool> FNInf;//whether the interface distance is obtained from shortcuts (vert)
	set<int> DisRe;//record the vertex id that the distance label should be updated
    set<int> DisRePost;//record the vertex id that the interface label should be updated
	vector<int> ch;//position index of children
	int height, hdepth;//hdepth is the deepest node that a vertex still exists
	int pa;//parent, the pa of root vertex is 0, position index
	int uniqueVertex;//vertex id of this tree node
//	vector<int> piv;//pivot vetex, used in path retrieval
//    int treeroot;//the tree id of subtree root, i.e., rank[x]
	Node(){
		vert.clear();
		pos.clear();
		dis.clear(); cnt.clear();
        disPost.clear(); cntPost.clear();
        disInf.clear();
        vAncestor.clear();
		ch.clear();
		pa = -1;
		uniqueVertex = -1;
		height = 0;
		hdepth = 0;
		FN.clear(); FNInf.clear();
		DisRe.clear();
        DisRePost.clear();
//		piv.clear();
//        treeroot=-1;
	}
};

class Graph{
public:
    string sourcePath;
    string dataset;
	int node_num=0;    //vertex number
	unsigned long long edge_num=0;    //edge number
	vector<vector<pair<vertex,int>>> Neighbor;//original graph
    vector<vector<pair<vertex,int>>> NeighborsParti;//<node_number,<in-partition adjacency lists>>
    vector<unordered_map<vertex,int>> NeighborsPartiPost;//adjacency lists of post-boundary partitions
    vector<vector<pair<vertex,int>>> NeighborsPartiPostV;//adjacency lists of post-boundary partitions
    vector<unordered_map<vertex,int>> NeighborsOverlay;//<node_number,<adjacency lists of overlay graph>>
    vector<vector<pair<vertex,int>>> NeighborsOverlayV;//<node_number,<adjacency lists of overlay graph>>
//    vector<unordered_map<int,int>> BoundaryShortcuts;
    vector<pair<int,bool>> PartiTag;//<node_number,<partition_id,if_boundary>>, for edge-cut
    vector<pair<int,set<int>>> PartiTagVC;//<node_number,<flag,set<partition_id>>>, flag=-1 means boundary vertex, flag>=0 means partition id, for vertex-cut
    vector<vector<vertex>> PartiVertex;//<partition_number,<in-partition vertices>>, in increasing vertex order, higher-rank vertex first
    vector<vector<vertex>> BoundVertex;//boundary vertices of each partition
//    vector<unordered_set<int>> BoundVertexSet;//boundary vertices of each partition
    vector<unordered_map<int,int>> BoundVertexMap;//map from boundary vertex id to its position id in BoundVertex
    vector<vertex> OverlayVertex;//overlay vertex in decreasing vertex order
//    vector<unordered_map<vertex,pair<int,int>>> repairShortcuts;//<ID1,ID2,<overlay weight, in-partition weight>>
    vector<map<int,int>> BoundShortcuts;//<ID1,ID2,weight>, source=-1 means stemming from core while source>=0 means stemming from specific partition
    int partiNum;   //actual partition number
    int pNum;   //nominated partition number
    bool ifParallel = true;
    bool ifIncrease = false;
    int cutType = 1;//1: edge-cut, 2: vertex-cut
    int nodeNumOverlay = 0;//vertex number of overlay graph
    vector<double> stageDurations;//durations of querying stages
    vector<pair<int,int>> GraphLocation;//Longitude and Latitude, already times 1000,000

    vector<bool> fullyConnected;//indicate whether this boundary vertex is fully-connected vertex
    bool ifFullOpt=false;//whether to use the overlay simplification optimization
    bool ifTreeOpt=false;//whether to use the tree decomposition optimization

    int HighestOrder;
    vector<bool> existOverlay;
    vector<int> partiRoots;//vertex id of partition root
    vector<map<int,map<int,int>>> SuppPartiID;//<ID1,<ID2,<pid,weight>>>, record the partition and its supportive value for a given interface edge
    vector<map<int,pair<int,set<int>>>> SuppPartiIDReal;//ID1,<ID2,<weight,set<pid>>>>, //record the partitions that really support a given interface edge
    set<int> affectedParti;//record the affected partitions
    vector<int> ProBeginIDV;
    vector<int> childNums;//children number of each vertex

	vector<int> DD; //intermediate variable in Contraction, DD2
	int threadnum=15;  //thread number
    string algoParti="NC";
    int algoChoice=2;//SP index. 1: CH; 2: H2H; 3: PLL
    int PSPStrategy=2;//PSP strategy. 1: Pre; 2: No; 3: Post
    int updateType=0;

	//vertex order
	vector<int> NodeOrder;//nodeID order, from 0
	vector<int> vNodeOrder;//order nodeID
    int treewidth=0;//treewidth
    int percent=0;

    /// Index Construction
//    vector<omp_lock_t> oml;
    unordered_map<int, Semaphore*> mSm;
    vector<Semaphore*> vSm;
    Semaphore* sm = new Semaphore(1);;// = new Semaphore(threadnum);

    //H2H index construction
    //intermediate variable and function used in the H2H index construction
    vector<map<int,pair<int,int>>> E;//ID1,ID2,(weight,count)
    vector<vector<pair<int,pair<int,int>>>> NeighborCon;//ID1,ID2,(weight,count)
    //for overlay graph
    vector<Node> Tree;
    vector<int> toRMQ;//toRMQ[p] records the first occurrence of node p in the Euler tour, i.e., R[p]
    vector<vector<int>> RMQIndex;//?
    vector<int> rank;//rank[v]>0 indicates non-core vertex
    int heightMax;
    vector<map<int, vector<pair<int,int>>>> SCconNodesMT;//<ID1,<ID2,<x,weight>>> where ID1<ID2, record the supportive vertices of a shortcut, only record edge once
    vector<vector<int>> VidtoTNid;//record the child tree nodes whose vert neighbors contain this tree node (nodeID--->tree node rank), for overlay

    //for no-boundary partitions
    vector<vertex> IDMap;//map the old id to new id in partitions
    vector<vector<Node>> Trees;//Trees for no-boundary
    vector<vector<int>> toRMQs;
    vector<vector<vector<int>>> RMQIndexs;
    vector<vector<int>> ranks;
    vector<int> heightMaxs;
    vector<map<int, vector<pair<int,int>>>> SCconNodesMTP;//for partitions. <ID1,<ID2,<x,weight>>> where ID1<ID2, record the supportive vertices of a shortcut, only record edge once
    vector<vector<int>> VidtoTNidP;//record the child tree nodes whose vert neighbors contain this tree node (nodeID--->tree node rank), for partition

    //for post-boundary partitions
//    vector<vertex> IDMapPost;//map the old id to new id in partitions
    vector<vector<Node>> TreesPost;//Trees for post-boundary
    vector<vector<int>> toRMQsPost;
    vector<vector<vector<int>>> RMQIndexsPost;
    vector<vector<int>> ranksPost;
    vector<int> heightMaxsPost;
    vector<map<int, vector<pair<int,int>>>> SCconNodesMTPost;//for partitions. <ID1,<ID2,<x,weight>>> where ID1<ID2, record the supportive vertices of a shortcut, only record edge once
    vector<vector<int>> VidtoTNidPost;//record the child tree nodes whose vert neighbors contain this tree node (nodeID--->tree node rank), for partition

    vector<bool> ifRepaired;

    vector<int> ProBeginVertexSetOverlay;
    set<int> vertexIDChLOverlay;
    vector<vector<int>> ProBeginVertexSetParti;
    vector<set<int>> vertexIDChLParti;


    vector<bool> vUpdated;// flag of whether the label of vertex has been updated
    vector<int> bHeights;
    int minWeight=INF;
    int maxWeight=0;
    int maxDegree=0;
    unsigned long long int qBNum=0;

    ~Graph(){
        clear();
    }
    void clear(){
        Neighbor.clear();
        Tree.clear();
        vSm.clear();
    }

    /// For non-partition SP index
    void HybridSPIndexConstruct();
    void MDEContraction(string orderfile);
    void NeighborComOrderGenerate(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p);
    void insertEMTOrderGenerate(int u,int v,int w);
    //for contraction
    void deleteEorder(int u,int v);
    void insertEorder(int u,int v,int w);
    int match(int x,vector<pair<int,pair<int,int>>> &vert);
    void makeTree();
    void makeIndex();
    void IndexsizeH2H();  //Index size computation

    /// For partitioned SP index
    void PCHIndexConstruct(int strategy);
    void PH2HIndexConstruct(int strategy);
    void PPLLIndexConstruct(int strategy);

    void PPLLIndexConstructVertexCut(int strategy);

    /// Pre-boundary Strategy
    void AllPairBoundaryDisCompute(bool ifParallel);
    void AllPairBoundaryDisUpdate(bool ifParallel, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>>& partiBatch, vector<pair<pair<int, int>, pair<int, int>>>& overlayBatch);
    void AllPairBoundaryDisUpdate(bool ifParallel, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>>& partiBatch, map<pair<int, int>, pair<int, int>>& overlayBatch);
    void AllPairBoundaryDisComputeParti(int pid);
    void AllPairBoundaryDisUpdateParti(int pid, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>>& partiBatch, vector<pair<pair<int, int>, pair<int, int>>>& overlayBatch);
    void AllPairBoundaryDisUpdatePartiMap(int pid, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>>& partiBatch, map<pair<int, int>, pair<int, int>>& overlayBatch);
    void AllPairBoundaryDisComputePartiV(vector<int>& p);
    void AllPairBoundaryDisUpdatePartiV(vector<int>& p, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>>& partiBatch, vector<pair<pair<int, int>, pair<int, int>>>& overlayBatch);
    void AllPairBoundaryDisUpdatePartiMapV(vector<int>& p, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>>& partiBatch, map<pair<int, int>, pair<int, int>>& overlayBatch);
    int BoundaryDijkstra(int ID1, vector<int> & IDs, vector<int>& Dis, vector<vector<pair<vertex,int>>> &Neighbor);


    ///PPLL
    bool ifVectorQ=false;//whether to use vector-based labels for querying
    int BPCLBSize=threadnum*10;//batch size of BPCL
    //overlay index
    vector<unordered_map<vertex,int>> Label;
    hl::Labeling LabelV;
    hl::PPR PPRV;
    vector<unordered_map<vertex,unordered_set<vertex>>> PruningPointSet;//{(v,c),u}
    vector<unordered_map<vertex,vertex>> PruningPointSet2;//<v,u,{c}>
//    vector<unordered_set<int>> ChangedLabels;
//    set<pair<int,int>> NoSupportedPair;



    //partition index
    vector<hl::Labeling> LabelVs;
    vector<hl::PPR> PPRVs;
    vector<unordered_map<vertex,int>> Labels;
//    vector<map<vertex,int>> Labels;
    vector<unordered_map<vertex,unordered_set<vertex>>> PruningPointSetP;//{(v,c),u}
    vector<unordered_map<vertex,vertex>> PruningPointSetP2;//<v,u,{c}>
    vector<unordered_map<vertex,int>> LabelsPost;
    vector<unordered_map<vertex,unordered_set<vertex>>> PruningPointSetPost;//{(v,c),u}
    vector<unordered_map<vertex,vertex>> PruningPointSetPost2;//<v,u,{c}>
    double coaseUpdateT=0;
    double refineUpdateT=0;
    double cleanUpdateT=0;
    unsigned long long int queryNum=0;
    unsigned long long int labelChangeNum=0;

//    vector<vector<unordered_map<vertex,int>>> LabelsBoundary;//PLL labels of boundary vertex for all partitions, the vertex ID is mapped by IDMap

    void ConstructPLL_OverlayIndex(vector<vector<pair<vertex,int>>> &Neighbor);
    void ConstructPLL_PartiIndex(bool ifParallel);
    void ConstructPLL_PartiIndexPost(bool ifParallel);
    void BPCLIndexConstructPartiVMap(vector<int>& p, vector<vector<pair<vertex, int>>> &Neighbor, vector<map<vertex,int>>& Label, vector<unordered_map<vertex,unordered_set<vertex>>>& PruningPointSet, vector<unordered_map<vertex,vertex>>& PruningPointSet2);
    void BPCLIndexConstructPartiV(vector<int>& p, vector<vector<pair<vertex, int>>> &Neighbor, vector<unordered_map<vertex,int>>& Label, vector<unordered_map<vertex,unordered_set<vertex>>>& PruningPointSet, vector<unordered_map<vertex,vertex>>& PruningPointSet2);
    void BPCLIndexConstructPartiMap(int pid, vector<vector<pair<vertex,int>>> &Neighbor, vector<map<vertex,int>>& Label, vector<unordered_map<vertex,unordered_set<vertex>>>& PruningPointSet, vector<unordered_map<vertex,vertex>>& PruningPointSet2);
    void BPCLIndexConstructParti(int pid, vector<vector<pair<vertex,int>>> &Neighbor, vector<unordered_map<vertex,int>>& Label, vector<unordered_map<vertex,unordered_set<vertex>>>& PruningPointSet, vector<unordered_map<vertex,vertex>>& PruningPointSet2);
    void BatchPCLDijk(vertex nodeID, unordered_set<vertex>& setNodes, vertex hID, vector<vector<pair<vertex,int>>> &Neighbor, vector<unordered_map<vertex,int>>& Label);
    void BatchPCLDijkNew(vertex nodeID, unordered_set<vertex>& setNodes, vertex hID, vector<vector<pair<vertex,int>>> &Neighbor, vector<unordered_map<vertex,int>>& Label, vector<tuple<int,int,int>> &labelTemp);
    void BatchPCLDijk2New(vector<vertex>& p, unordered_set<vertex>& setNodes, vertex hID, vector<vector<pair<vertex,int>>> &Neighbor, vector<unordered_map<vertex,int>>& Label, vector<tuple<int,int,int>> &labelTemp);
    void BatchPCLDijkNewMap(vertex nodeID, unordered_set<vertex>& setNodes, vertex hID, vector<vector<pair<vertex,int>>> &Neighbor, vector<map<vertex,int>>& Label, vector<tuple<int,int,int>> &labelTemp);
    void BatchPCLDijk2NewMap(vector<vertex>& p, unordered_set<vertex>& setNodes, vertex hID, vector<vector<pair<vertex,int>>> &Neighbor, vector<map<vertex,int>>& Label, vector<tuple<int,int,int>> &labelTemp);
    void BatchPCLDijk2(vector<vertex>& p, unordered_set<vertex>& setNodes, vertex hID, vector<vector<pair<vertex,int>>> &Neighbor, vector<unordered_map<vertex,int>>& Label);
    int ShortestDisQuery1(int ID1,int ID2,vector<int>& SupNode, int& d, vector<unordered_map<vertex,int>>& Label);
    int ShortestDisQuery1(int ID1,int ID2,vector<int>& SupNode, int& d, vector<map<vertex,int>>& Label);
    int ShortestDisQuery2(int ID1,int ID2,vector<int>& SupNode, int& d, vector<unordered_map<vertex,int>>& Label);
    int ShortestDisQuery2(int ID1,int ID2,vector<int>& SupNode, int& d, vector<map<vertex,int>>& Label);
    void PPRConstruction(vertex nodeID, vector<vector<pair<vertex,int>>> &Neighbor, vector<unordered_map<vertex,int>>& Label, vector<unordered_map<vertex,unordered_set<vertex>>>& PruningPointSet, vector<unordered_map<vertex,vertex>>& PruningPointSet2);
    void PPRConstruction(vertex nodeID, vector<vector<pair<vertex,int>>> &Neighbor, vector<map<vertex,int>>& Label, vector<unordered_map<vertex,unordered_set<vertex>>>& PruningPointSet, vector<unordered_map<vertex,vertex>>& PruningPointSet2);
    void PPRConstruction2(vector<vertex> & p, vector<vector<pair<vertex,int>>> &Neighbor, vector<unordered_map<vertex,int>>& Label, vector<unordered_map<vertex,unordered_set<vertex>>>& PruningPointSet, vector<unordered_map<vertex,vertex>>& PruningPointSet2);
    void PPRConstruction2Map(vector<vertex> & p, vector<vector<pair<vertex,int>>> &Neighbor, vector<map<vertex,int>>& Label, vector<unordered_map<vertex,unordered_set<vertex>>>& PruningPointSet, vector<unordered_map<vertex,vertex>>& PruningPointSet2);
    void PruningPointBuild(bool ifParallel, vector<vector<vertex>> & processID, vector<vector<pair<vertex,int>>> &Neighbor, vector<unordered_map<vertex,int>>& Label,  vector<unordered_map<vertex,unordered_set<vertex>>>& PruningPointSet, vector<unordered_map<vertex,vertex>>& PruningPointSet2);
    void PruningPointBuild(bool ifParallel, vector<vector<vertex>> & processID, vector<vector<pair<vertex,int>>> &Neighbor, vector<map<vertex,int>>& Label,  vector<unordered_map<vertex,unordered_set<vertex>>>& PruningPointSet, vector<unordered_map<vertex,vertex>>& PruningPointSet2);
    void PLLConstructV(vector<vector<pair<vertex,int>>>& Neighbor);
    void DijksPrune1V(int nodeID,vector<vector<pair<vertex,int>>>& Neighbor);
    int PLLDisQuery1V(int ID1,unordered_map<vertex,int>& Lh,vector<int>& SupNode);
    void ConstructPLL_OverlayGraph(bool ifParallel);
    void ConstructBoundaryShortcutPLL(int pid);
    void ConstructBoundaryShortcutPLLV(vector<int>& p);
    void ConstructPartitionPostPLL(bool ifParallel);
    void ConstructPostPartiPLL(int pid);
    void ConstructPostPartiPLLV(vector<int>& p);
    void IndexSizePPLL();

    void DecreasePartiBatchUpdatePLL(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay, bool ifOverlay);
    void PPLLBatchUpdateDec(vector<pair<pair<int, int>, pair<int, int>>> &wBatch);
    void PPLLBatchUpdateDecPre(vector<pair<pair<int, int>, pair<int, int>>> &wBatch);
    void DecreasePSL(int a, int b, int oldW, int newW, vector<vector<pair<vertex,int>>> &Neighbors,vector<unordered_map<vertex,int>> &Label);
    void DecreasePSL(int a, int b, int oldW, int newW, vector<vector<pair<vertex,int>>> &Neighbors,vector<map<vertex,int>> &Label);
    int PLLQuery(int ID1,int ID2,vector<unordered_map<vertex,int>> &Label);
    int PLLQuery(int ID1,int ID2,vector<map<vertex,int>> &Label);
    void Repair_PartiIndexPLL(bool ifParallel, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);
    void RepairPartitionIndexDecreasePLL(int pid, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);
    void RepairPartitionIndexPLLV(vector<int>& p, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);
    void RepairPartitionIndexIncreasePLL(int pid, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);
    void RepairPartitionIndexEdgeDeletePLL(int pid, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);

    void EdgeDeletePartiBatchUpdatePLL(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay, bool ifOverlay);
    void IncreasePartiBatchUpdatePLL(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay, bool ifOverlay);
    void PPLLBatchUpdateEdgeDelete(vector<pair<pair<int, int>, pair<int, int>>> &wBatch);
    void PPLLBatchUpdateEdgeDeletePre(vector<pair<pair<int, int>, pair<int, int>>> &wBatch);
    void PPLLBatchUpdateInc(vector<pair<pair<int, int>, pair<int, int>>> &wBatch);
    void PPLLBatchUpdateIncPre(vector<pair<pair<int, int>, pair<int, int>>> &wBatch);
    void EdgeDeletePSL(int a, int b, int oldW, vector<vector<pair<vertex,int>>> &Neighbor,vector<unordered_map<vertex,int>> &Label, vector<unordered_map<vertex,unordered_set<vertex>>> &PruningPointNew, vector<unordered_map<vertex,vertex>>& PruningPointSet2);
    void EdgeDeletePSL(int a, int b, int oldW, vector<vector<pair<vertex,int>>> &Neighbor,vector<map<vertex,int>> &Label, vector<unordered_map<vertex,unordered_set<vertex>>> &PruningPointNew, vector<unordered_map<vertex,vertex>>& PruningPointSet2);
    void IncreasePSL(int a, int b, int oldW, int newW, vector<vector<pair<vertex,int>>> &Neighbor,vector<map<vertex,int>> &Label, vector<unordered_map<vertex,unordered_set<vertex>>> &PruningPointNew, vector<unordered_map<vertex,vertex>>& PruningPointSet2);//set version with NoSupportedPair
    void IncreasePSL(int a, int b, int oldW, int newW, vector<vector<pair<vertex,int>>> &Neighbor,vector<unordered_map<vertex,int>> &Label, vector<unordered_map<vertex,unordered_set<vertex>>> &PruningPointNew, vector<unordered_map<vertex,vertex>>& PruningPointSet2);//set version with NoSupportedPair
    void CoarseUpdate(int LID,int HID, int oldW, queue<pair<int,pair<int,int>>>& WaitPro, queue<pair<int,pair<int,int>>>& WaitProP, vector<pair<int,int>>& AL1, vector<pair<int,int>>& AL2, vector<pair<int,int>>& AL2Check, vector<vector<pair<vertex,int>>> &Neighbor,vector<unordered_map<vertex,int>> &Label, bool ifDebug, int lid, int hid, vector<unordered_set<int>>& ChangedLabels);
    void RefineUpdate(queue<pair<int,pair<int,int>>>& WaitPro, queue<pair<int,pair<int,int>>>& WaitProP, vector<pair<int,int>>& AL1, vector<pair<int,int>>& AL2, vector<pair<int,int>>& AL2Check, set<tuple<vertex,vertex,vertex>>& outdatedPruning, map<pair<vertex,vertex>,vertex>& newPruningPoints, vector<unordered_map<vertex,vertex>>& PruningPointSet2, vector<vector<pair<vertex,int>>> &Neighbor,vector<unordered_map<vertex,int>> &Label,vector<unordered_map<vertex,unordered_set<vertex>>> &PruningPointNew, bool ifDebug, int lid, int hid, vector<unordered_set<int>>& ChangedLabels, set<pair<int,int>>& NoSupportedPair);//queue version
    bool PPRCheck(int curID, int hubID, vector<vector<pair<vertex,int>>> &Neighbor, vector<unordered_map<vertex,int>> &Label, vector<unordered_map<vertex,unordered_set<vertex>>> &PruningPointNew, queue<pair<int,pair<int,int>>>& WaitProPTem, vector<pair<int,int>> &AL2, map<pair<vertex,vertex>,vertex> &newPruningPoints, vector<unordered_map<vertex,vertex>>& PruningPointSet2, set<tuple<vertex,vertex,vertex>> &outdatedPruning, bool ifDebug, int lid, int hid, vector<unordered_set<int>>& ChangedLabels, set<pair<int,int>>& NoSupportedPair);//queue version
    void PPRClean(vector<vector<pair<vertex,int>>> &Neighbors, map<pair<vertex,vertex>,vertex> &newPruningPoints, set<tuple<vertex,vertex,vertex>> &outdatedPruning, bool ifDebug, int lid, int hid, vector<map<vertex,int>> &Label, vector<unordered_map<vertex,unordered_set<vertex>>>& PruningPointSet, vector<unordered_map<vertex,vertex>>& PruningPointSet2, vector<unordered_set<int>>& ChangedLabels);
    void CoarseUpdate(int LID,int HID, int oldW, queue<pair<int,pair<int,int>>>& WaitPro, queue<pair<int,pair<int,int>>>& WaitProP, vector<pair<int,int>>& AL1, vector<pair<int,int>>& AL2, vector<pair<int,int>>& AL2Check, vector<vector<pair<vertex,int>>> &Neighbor,vector<map<vertex,int>> &Label, bool ifDebug, int lid, int hid, vector<unordered_set<int>>& ChangedLabels);
    void RefineUpdate(queue<pair<int,pair<int,int>>>& WaitPro, queue<pair<int,pair<int,int>>>& WaitProP, vector<pair<int,int>>& AL1, vector<pair<int,int>>& AL2, vector<pair<int,int>>& AL2Check, set<tuple<vertex,vertex,vertex>>& outdatedPruning, map<pair<vertex,vertex>,vertex>& newPruningPoints, vector<unordered_map<vertex,vertex>>& PruningPointSet2, vector<vector<pair<vertex,int>>> &Neighbor,vector<map<vertex,int>> &Label,vector<unordered_map<vertex,unordered_set<vertex>>> &PruningPointNew, bool ifDebug, int lid, int hid, vector<unordered_set<int>>& ChangedLabels, set<pair<int,int>>& NoSupportedPair);//queue version
    bool PPRCheck(int curID, int hubID, vector<vector<pair<vertex,int>>> &Neighbor, vector<map<vertex,int>> &Label, vector<unordered_map<vertex,unordered_set<vertex>>> &PruningPointNew, queue<pair<int,pair<int,int>>>& WaitProPTem, vector<pair<int,int>> &AL2, map<pair<vertex,vertex>,vertex> &newPruningPoints, vector<unordered_map<vertex,vertex>>& PruningPointSet2, set<tuple<vertex,vertex,vertex>> &outdatedPruning, bool ifDebug, int lid, int hid, vector<unordered_set<int>>& ChangedLabels, set<pair<int,int>>& NoSupportedPair);//queue version
    void PPRClean(vector<vector<pair<vertex,int>>> &Neighbors, map<pair<vertex,vertex>,vertex> &newPruningPoints, set<tuple<vertex,vertex,vertex>> &outdatedPruning, bool ifDebug, int lid, int hid, vector<unordered_map<vertex,int>> &Label, vector<unordered_map<vertex,unordered_set<vertex>>>& PruningPointSet, vector<unordered_map<vertex,vertex>>& PruningPointSet2, vector<unordered_set<int>>& ChangedLabels);
    int DisQueryLower1(int ID1, int ID2, vector<vector<pair<vertex,int>>> &Neighbors,vector<unordered_map<vertex,int>> &Label);
    int DisQueryVally(int ID1, int ID2, vector<vector<pair<vertex,int>>> &Neighbors,vector<unordered_map<vertex,int>> &Label);
    pair<int,int> DisQueryVally2(int ID1, int ID2, vector<vector<pair<vertex,int>>> &Neighbors,vector<unordered_map<vertex,int>> &Label);
//    pair<int,int> DisQueryVallyVert2(int ID1, int ID2, vector<vector<pair<vertex,int>>> &Neighbors,vector<unordered_map<vertex,int>> &Label);
    int DisQueryPeak(int ID1, int ID2,vector<unordered_map<vertex,int>> &Label);
    pair<int,int> DisQueryPeak2(int ID1, int ID2,vector<unordered_map<vertex,int>> &Label);
    int DisQueryLower1(int ID1, int ID2, vector<vector<pair<vertex,int>>> &Neighbors,vector<map<vertex,int>> &Label);
    int DisQueryVally(int ID1, int ID2, vector<vector<pair<vertex,int>>> &Neighbors,vector<map<vertex,int>> &Label);
    pair<int,int> DisQueryVally2(int ID1, int ID2, vector<vector<pair<vertex,int>>> &Neighbors,vector<map<vertex,int>> &Label);
    int DisQueryPeak(int ID1, int ID2,vector<map<vertex,int>> &Label);
    pair<int,int> DisQueryPeak2(int ID1, int ID2,vector<map<vertex,int>> &Label);

    //For edge insertion
    void EdgeInsertionPLL(int ID1, int ID2, int newW, vector<vector<pair<vertex,int>>> &Neighbors,vector<unordered_map<vertex,int>> &Label);
    void EdgeInsertionPLL(int ID1, int ID2, int newW, vector<vector<pair<vertex,int>>> &Neighbors,vector<map<vertex,int>> &Label);
    void VertexInsertPartiBatchUpdatePLL(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<vertex,int>>> &Neighbors);
    void EdgeInsertPartiBatchUpdatePLL(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay, bool ifOverlay, vector<pair<pair<int,int>,int>>& updatedSC);
    void PPLLBatchUpdateEdgeInsert(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1);
    void PPLLBatchUpdateEdgeInsertPre(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1);
    void PPLLBatchUpdateVertexInsert(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1);
    void PPLLBatchUpdateVertexInsertPre(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1);

    /// For PH2H
    void PH2HIndexConstruct(); //PH2H index construction
    void ConstructBoundaryShortcutV(vector<int> & p, bool ifAllPair, bool ifCH);
    void ConstructBoundaryShortcut(int pid, bool ifCH);
    void ConstructBoundaryShortcutNoAllPair(int pid);
    void Construct_PartiIndex(bool ifParallel, bool ifLabelC);
    void Construct_OverlayGraph(bool ifParallel, bool ifCH);
    void Construct_OverlayGraphNoAllPair(bool ifParallel, bool ifCH);
    void Construct_OverlayIndex(bool ifLabelC);

    void PreConstructAllPairs(bool ifParallel);
    void PreConstructAllPairsPartiV(vector<int> & p);
    void PreConstructAllPairsParti(int pid);

    /// Post-boundary index
    void RefreshBoundaryEdgesAndLabelingPartiV(vector<int>& p);
    void RefreshBoundaryEdgesAndLabelingParti(int pid);
    void ConstructPartitionPost(bool ifParallel, bool ifCH);
    void ConstructPostParti(int pid, bool ifCH);
    void ConstructPostPartiV(vector<int>& p, bool ifCH);
    void ConstructPartitionPostIndex(bool ifParallel, bool ifLabelU);
    void ConstructPartitionPostIndexOpt(bool ifParallel);
    void Repair_PartiIndex(bool ifParallel, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);
    void Repair_PartiIndexPostMHLPost(bool ifParallel, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch, double & runT);
    void Repair_PartiIndexForOpt(bool ifParallel, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);
    void RepairPartitionIndexV(vector<int>& p, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch, bool ifOpt);
    void RepairPartitionIndexDecrease(int pid, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);
    void RepairPartitionIndexDecreaseForOpt(int pid, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);
    void RepairPartitionIndexIncrease(int pid, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);
    void RepairPartitionIndexIncreaseForOpt(int pid, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);


    void ThreadDistribute(vector<int>& vertices, vector<vector<int>>& processID);
    void ConstructPH2H_PartiV(vector<int> & P, vector<vector<Node>>& Trees, vector<vector<int>>& ranks, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs);
    void ConstructPH2H_PartiVCH(vector<int> & P, vector<vector<Node>>& Trees, vector<vector<int>>& ranks, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs);
    void ConstructPH2H_Parti(int pid, vector<vector<Node>>& Trees, vector<vector<int>>& ranks, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs);
    void ConstructPH2H_PartiLabel(int pid, vector<vector<Node>>& Trees, vector<vector<int>>& ranks, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs);
    void H2HCreateTree_Parti(int pid, vector<Node>& TreeP, vector<int>& rankP, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs);//Create tree for partition
    void H2HCreateIndex_Parti(int pid, vector<Node>& TreeP, vector<int>& rankP);//Create labels for partition
    void H2HCreateTree_Overlay();
    void H2HCreateIndex_Overlay();
    int ShortcutDisCheck(int ID1, int ID2);
    void makeTreeIndexDFSP(int p, vector<int>& list,  vector<Node>& TreeP, vector<int>& rankP);
    void makeRMQCoreP(int pid, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs, vector<vector<Node>>& Trees);
    void makeRMQDFSCoreP(int pid, int p, int height, vector<int>& EulerSeqP, vector<vector<int>>& toRMQs, vector<vector<Node>>& Trees);
    void makeRMQCore();
    void makeRMQDFSCore(int p, int height, vector<int>& EulerSeq);
    void makeIndexDFS(int p, vector<int> &list);
    void makeRMQ(vector<int>& toRMQ, vector<vector<int>>& RMQIndex, vector<Node>& Tree);
    void makeRMQDFS(int p, int height, vector<int>& EulerSeq, vector<int>& toRMQ, vector<Node>& Tree);
    void deleteECore(int u,int v);
    void insertECore(int u,int v,int w);
    void insertECoreMT(int u,int v,int w);
    int matchCore(int x,vector<pair<int,pair<int,int>>> &vert, vector<int>& rank);//vector<pair<int,int>> &vert
    int matchCoreParti(int x,vector<pair<int,pair<int,int>>> &vert, vector<int>& rank);//vector<pair<int,int>> &vert
    void IndexSizePH2H();  //Core-tree index size computation
    void NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x);
    void NeighborComorderH2H(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x);
    void insertEMTorder(int u,int v,int w);


    void VertexCutVertexOrdering();
    void EdgeCutVertexOrdering(int type);
    void SketchGraphBuild();
    void OverlayOrderingBuild();
    void OverlayOrderingBuildBoundaryFirst(int nodenum, vector<vector<int>> Neighbor);
    void PartitionOrderingBuildMDE(bool ifParallel);
    void PartitionOrderingBuildDegree(bool ifParallel);
    void OrderingAssemblyMDEBoundaryFirst(string filename);
    void OrderingAssemblyMDE(int pNum);
    void OrderingAssemblyBoundaryFirst(int pNum);
    void SketchOrder(vector<vector<pair<int,int>>> Neighbor, vector<int> &vNodeOrderSketch);
    void PartitionOrderingMDEV(vector<int>& p);
    void PartitionOrderingMDE(int pid);
    void PartitionOrderingDegreeV(vector<int>& p);
    void PartitionOrderingDegree(int pid);
    void deleteEOrderGenerate(int u,int v);
    void insertEOrderGenerate(int u,int v,int w);

//    vector<pair<pair<int,int>,int>> CutEdges;//the cut edges
    vector<vector<int>> NeighborSketch;
    vector<set<int>> NeighborSketchS;
//    vector<map<int,int>> vNodeOrderParti;
    vector<vector<int>> vNodeOrderParti;
    vector<int> vNodeOrderOverlay;


	///Query processing
    int QueryPCH(int ID1, int ID2, vector<vector<Node>>& Trees);
    int QueryPCHDebug(int ID1, int ID2, vector<vector<Node>>& Trees);
    int QueryPH2H(int ID1, int ID2);
    int QueryPPLL(int ID1, int ID2);
    int QueryPPLLDebug(int ID1, int ID2);

	//PH2H
    int QueryCore(int ID1, int ID2);
    int QueryH2HPartition(int ID1, int ID2, int PID);
    int QueryH2HPartitionDebug(int ID1, int ID2, int PID);
    int QueryH2HPartitionPost(int ID1, int ID2, int PID);
    int QueryPartiCore(int ID1, int ID2);
    int QueryPartiCoreDebug(int ID1, int ID2);
    int QuerySameParti(int ID1, int ID2);
    int QuerySamePartiPost(int ID1, int ID2);
    int QuerySamePartiPre(int ID1, int ID2);
    int QuerySamePartiPostOpt(int ID1, int ID2);
    int QueryPartiParti(int ID1, int ID2);



    //PCH
    int QueryCoreCH(int ID1, int ID2);
    int QueryCHPartition(int ID1, int ID2, int PID, vector<vector<Node>>& Trees);
    int QueryPartiCoreCH(int ID1, int ID2);
//    int QuerySamePartiCH(int ID1, int ID2);
    int QueryPartiPartiCH(int ID1, int ID2);


    //PPLL
    int QueryOverlayPLL(int ID1, int ID2);
    int QueryOverlayPLLDebug(int ID1, int ID2);
    int QuerySamePartiPLL(int ID1, int ID2, vector<unordered_map<vertex,int>>& Labels);
    int QuerySamePartiPLL(int ID1, int ID2, vector<map<vertex,int>>& Labels);
    int QueryPartiCorePLL(int ID1, int ID2);
    int QuerySamePartiPLLNo(int ID1, int ID2);
    int QueryPartiPartiPLL(int ID1, int ID2);
    int QueryPartiPartiPLLDebug(int ID1, int ID2);

    void EffiCheck(string filename,int runtimes);
    int QueryDebug(int ID1, int ID2);
    int QueryCoreDebug(int ID1, int ID2);
    int LCAQuery(int _p, int _q);
    int LCAQueryPartition(int _p, int _q, int PID);// query within partition
    int LCAQueryPartitionPost(int _p, int _q, int PID);// query within partition
    int LCAQueryOverlay(int _p, int _q);
    int LCAQuery(int _p, int _q, vector<int>& toRMQ, vector<vector<int>>& RMQIndex, vector<Node>& Tree);

    //Correctness Check
    void CorrectnessCheck(int runtimes);
    void CorrectnessCheckCore(int runtimes);
    void DFSTree(vector<int>& tNodes, int id);


    //Dijkstra
    int Dijkstra(int ID1, int ID2,vector<vector<pair<vertex,int>>> &Neighbor);
    int BiDijkstra(int ID1, int ID2,vector<vector<pair<vertex,int>>> &Neighbor);
    int Astar(int ID1, int ID2,vector<vector<pair<vertex,int>>> &Neighbor);
    int EuclideanDis(int s, int t);
    double EuclideanDis(pair<double,double> s, pair<double,double> t);
    void RetrievePath(int ID1, int ID2, vector<int> & prece,vector<vector<pair<vertex,int>>> &Neighbor);
    void RetrievePath(int ID1, int ID2, vector<int> & prece,vector<unordered_map<vertex,int>> &Neighbor);
    int DijkstraCore(int ID1, int ID2);
    int DijkstraPath(int ID1, int ID2,vector<vector<pair<vertex,int>>> &Neighbor);
    int DijkstraPath(int ID1, int ID2,vector<unordered_map<vertex,int>> &Neighbor);

    void EachNodeProBDis5H2H(int child,vector<int>& line,set<int>& vertexIDChL, map<int,int>& checkedDis);
    void eachNodeProcessIncrease1H2H(int children, vector<int>& line, int& changelabel);
    void PCHBatchUpdateDec(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    void PCHBatchUpdateEdgeInsert(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    void PCHBatchUpdateEdgeInsertPre(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    void PCHBatchUpdateVertexInsert(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    void PCHBatchUpdateVertexInsertPre(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    void PCHBatchUpdateDecPre(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    void PH2HBatchUpdateDec(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    void PH2HBatchUpdateDecPre(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    void PH2HBatchUpdateVertexInsert(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    void PH2HBatchUpdateVertexInsertPre(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    void PH2HBatchUpdateEdgeInsert(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    void PH2HBatchUpdateEdgeInsertPre(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    void PCHBatchUpdateInc(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    void PCHBatchUpdateEdgeDelete(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1);
    void PCHBatchUpdateEdgeDeletePre(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    void PCHBatchUpdateIncPre(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    void PH2HBatchUpdateInc(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    void PH2HBatchUpdateIncPre(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    void PH2HBatchUpdateEdgeDelete(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    void PH2HBatchUpdateEdgeDeletePre(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    void IndexConstruction();

	/// Index update
    void IndexMaintenance(int updateType, int batchNum, int batchSize, double changeRatio);
    void DecreaseOverlay(int a,int b, int newW, vector<unordered_map<vertex,int>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax);
    void DecreaseParti(int a,int b, int newW, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax);
    void EachNodeProBDis5(int child,vector<int>& line,set<int>& vertexIDChL, vector<Node> &Tree, vector<int> &rank);//map<int,int>& checkedDis,
    void EachNodeProBDis5PostMHLOverlay(int child,vector<int>& line,set<int>& vertexIDChL, vector<Node> &Tree, vector<int> &rank);//map<int,int>& checkedDis,
    void EachNodeProBDis5Parti(int child,vector<int>& line,set<int>& vertexIDChL, vector<Node> &Tree, vector<int> &rank);//map<int,int>& checkedDis, map<int,int>& checkedDis,
    void IncreaseOverlay(int a, int b, int oldW, int newW, vector<unordered_map<int,int>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid);
    void eachNodeProcessIncrease1(int children, vector<int>& line, int& changelabel, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid);
    void IncreaseParti(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid);
    void eachNodeProcessIncrease1Parti(int children, vector<int>& line, int& changelabel, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid);
    void DecreaseOverlayBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU);//batch decrease for overlay graph
    void EdgeInsertOverlayBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU);
    void EdgeInsertOverlayBatchCH(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax);
    void DecreaseH2HBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU);
    void DecreaseOverlayBatchLabel(vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int>& ProBeginVertexSet,set<int>& vertexIDChL);
    void DecreasePartiBatchUpdateCheck(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay);
    void DecreasePartiBatchUpdateCheckCH(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay, bool ifOpt, vector<pair<pair<int,int>,int>>& updatedSC);
    void EdgeInsertPartiBatchUpdateCheck(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay, bool ifLabel, vector<pair<pair<int,int>,int>>& updatedSC);
    void BottomUpNewShortcutInsertParti(int lid, int hid, int weight, vector<Node> &Tree, vector<int> &rank, int &newSCNum);
    void BottomUpNewShortcutInsertPartiCH(int lid, int hid, int weight, vector<Node> &Tree, vector<int> &rank, int &newSCNum);
    void BottomUpNewShortcutInsert(int lid, int hid, int weight, vector<Node> &Tree, vector<int> &rank, int &newSCNum);
    void BottomUpNewShortcutInsertCH(int lid, int hid, int weight, vector<Node> &Tree, vector<int> &rank, int &newSCNum);
    void EdgeInsertPartiBatchCH(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>>& Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<pair<pair<int,int>,int>>& updatedSC);
    void EdgeInsertPartiBatchH2H(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>>& Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<pair<pair<int,int>,int>>& updatedSC, bool ifLabelU);
    void EdgeInsertPartiBatchH2HPre(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>>& Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU);
    void EdgeInsertPartiBatchH2HPost(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<vertex,int>>& Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU);
    void VertexInsertPartiBatchUpdateCH(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch);
    void VertexInsertPartiBatchUpdateH2H(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch);
    void DecreasePartiBatch(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<pair<pair<int,int>,int>>& updatedSC, bool ifLabelU);
    void DecreasePartiBatchForOpt(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<pair<pair<int,int>,int>>& updatedSC, bool ifLabelU, bool ifConstruct);
    void DecreasePartiBatchPost(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>>& Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU);
    void DecreasePartiBatchPre(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>>& Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU);
    void DecreasePartiBatchLabel(vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int>& ProBeginVertexSet,set<int>& vertexIDChL);
    void IncreaseOverlayBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, bool ifLabelU);
    void IncreaseH2HBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, bool ifLabelU);
    void IncreaseOverlayBatchLabel(vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int>& ProBeginVertexSet, vector<vector<int>> &VidtoTNid);
    void IncreasePartiBatchUpdateCheck(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay, bool ifLabel,vector<pair<pair<int,int>,int>>& updatedSC);
    void IncreasePartiBatch(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, vector<pair<pair<int,int>,int>>& updatedSC, bool ifLabelU);
    void IncreasePartiBatchForOpt(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, vector<pair<pair<int,int>,int>>& updatedSC, bool ifLabelU);
    void IncreasePartiBatchPost(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, bool ifLabelU);
    void IncreasePartiBatchPre(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU);
    void EdgeDeletePartiBatchPre(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU);
    void IncreasePartiBatchLabel(vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int>& ProBeginVertexSet, vector<vector<int>> &VidtoTNid);

	/// Graph Preprocessing
    void ReadGraph(string filename);
    void ReadCoordinate(string filename);
    void ReadUpdate(string filename,vector<pair<pair<int,int>,int>>& TestData);
    void ReadUpdate2(string filename,vector<pair<pair<int,int>,pair<int,int>>>& TestData);
    void ReadUpdate3(string filename,vector<pair<pair<int,int>,tuple<int,int,int>>>& TestData);
	void StainingMethod(int ID);
	void ODGene(int num, string filename);
    void ODGeneParti(int num, string filename, double portion);
    void ODGeneSameParti(int num, string filename);
    void ODGeneCrossParti(int num, string filename);
	void UpdateGene(int num, string filename);
    void QueryGenerationParti(bool ifSame);
    void StructuralUpdateGeneration(int number);

    void WriteTreeIndexOverlay(string filename);
    void ReadTreeIndex(string file);
    void WriteTreeIndexParti(string filename);
    void WriteGraph(string graphfile);
    void WriteOrder(string filename);
    void ReadOrder(string filename);
    void CompareOrder(string filename1, string filename2);
    void GraphPartitionRead(string filename);
    void GraphPartitionReadVertexCut(string filename);

    void WriteCoreIndex(string file);
    void ReadCoreIndex(string file);
    void WriteCoreGraph(string graphfile);
    void ReadCoreGraph(string filename);

    vector<int> DFS_CC(vector<map<int,int>> & Edges, set<int> set_A, set<int> & set_B, int nodenum);
    vector<int> DFS_CC(vector<vector<pair<int,int>>> & Edges, set<int> set_A, set<int> & set_B, int nodenum);

};

#endif // HEADPSP_H_
