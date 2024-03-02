/*
 * head.h
 *
 *  Created on: 24 August 2022
 *      Author: Xinjie ZHOU
 */

#ifndef HEADSP_H_
#define HEADSP_H_

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
#include <omp.h>

#define INF 99999999
//typedef unsigned int vertex;
typedef int vertex;

using namespace std;
//using namespace boost;


extern vector<int> NodeOrder_;//nodeID order
extern vector<int> _DD_;
extern vector<int> NodeOrders;

struct Nei{
	int nid;
	int w;
	int c;
};

struct OrderComp{// maximum-first, Higher-order first
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

struct OrderCompp{//prior to reture the vertex with smaller order
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


struct OrderCompCH{
    int x;
    int y;//order(x)<order(y)
    OrderCompCH(int _x, int _y){
        x=_x; y=_y;
    }
    bool operator< (const OrderCompCH& d) const{
        if(x==d.x && y==d.y){//avoid the redundant
            return false;
        }else{
            if(x!=d.x)
                return NodeOrders[x]<NodeOrders[d.x];
            if(y!=d.y)
                return NodeOrders[y]<NodeOrders[d.y];
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
};

struct Node{//tree node
	vector<pair<int,pair<int,int>>> vert;//neighID/weight/count(how many ways can lead to this super edge weight)
//    vector<pair<int,pair<int,int>>> vertNo;//neighID/weight/count(how many ways can lead to this super edge weight)
	vector<int> pos;
	vector<int> dis, cnt;//the distance value and corresponding count number (i.e., record how many path has the shortest distance)
//    vector<int> disNo;//the distance value of post-boundary strategy
    vector<int> vAncestor;//the ancestors, which is corresponding to dis
	//vector<set<int>> FromNode;
//	set<int> changedPos;
	vector<bool> FN;//another succinct way of FromNode, whether this distance label is directly obtained from shortcuts (vert), i.e., one-hop. there still be a case that the dis is also supported by neighbor's label even it is true.
	set<int> DisRe;//record the vertex id that the distance label should be updated
	vector<int> ch;
	int height=0;//tree height of a tree node
    int hdepth=0;//hdepty is the deepest node that a vertex still exists
	int pa;//parent, the pa of root vertex is 0
	int uniqueVertex;//?vertex id of this tree node?
//	vector<int> piv;//pivot vetex, used in path retrieval
//    int treeroot;//the tree id of subtree root, i.e., rank[x]
	Node(){
		vert.clear();
//		neighInf.clear();
		pos.clear();
		dis.clear();
		cnt.clear();
        vAncestor.clear();
		ch.clear();
		pa = -1;
		uniqueVertex = -1;
		height = 0;
		hdepth = 0;
//		changedPos.clear();
		FN.clear();
		DisRe.clear();
//		piv.clear();
//        treeroot=-1;
	}
};

class Graph{
public:
    string graphfile;
    string dataset;
	int node_num=0;    //vertex number
	unsigned long long edge_num=0;    //edge number
	vector<vector<pair<vertex,int>>> Neighbor;//original graph
    int partiNum=16;   //partition number
    int algoIndex=2;    //SP index, 1: CH; 2:H2H

	vector<int> DD; //intermediate variable in Contraction, DD2
	int threadnum=15;  //thread number
    string algoParti="NC";

	//vertex order
	vector<int> NodeOrder;//nodeID order
	vector<int> vNodeOrder;//order nodeID

    /// Index Construction
//    vector<omp_lock_t> oml;
    unordered_map<int, Semaphore*> mSm;
    vector<Semaphore*> vSm;
    Semaphore* sm = new Semaphore(1);// = new Semaphore(threadnum);

    //H2H index construction
    //intermediate variable and function used in the H2H index construction
    vector<vector<pair<int,pair<int,int>>>> NeighborCon;//ID1,ID2,(weight,count)
    vector<map<int, vector<int>>> SCconNodesMT;// <ID1,ID2,<supportive vertices>> where ID1<ID2, record the supportive vertices of a shortcut, only record edge once
    vector<map<int,pair<int,int>>> E;//ID1,ID2,(weight,count)
    vector<vector<int>> VidtoTNid;// (nodeID,vector<tree node rank>), record the tree node id whose unique vertex involves this vertex as neighbor
    vector<int> rank;//rank[v]>0 indicates non-core vertex
    int heightMax;
    vector<Node> Tree;
    vector<int> EulerSeq;//prepare for the LCA calculation, EulerSeq is the Euler tour, i.e., E[1,...,2n-1]
    vector<int> toRMQ;//toRMQ[p] records the first occurrence of node p in the Euler tour, i.e., R[p]
    vector<vector<int>> RMQIndex;//For LCA query


    ~Graph(){
        clear();
    }
    void clear(){
        Neighbor.clear();
        Tree.clear();
        vSm.clear();
    }
    //// Index Construction
    void IndexConstruction(int algo);

    //// For CH index construction
    void CHIndexConstruct();
    void IndexsizeCHWP();

    //// For H2H index construction
    void H2HIndexConstruct(); //H2H index construction
    /// MDE contraction
    void MDEContraction(string orderfile);
    //for order generation
    void deleteEOrderGenerate(int u,int v);
    void NeighborComOrderGenerate(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p);
    void insertEMTOrderGenerate(int u,int v,int w);
    //for contraction
    void deleteEorder(int u,int v);
    void insertEorder(int u,int v,int w);
    void insertEMTorder(int u,int v,int w);
    void NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x);
    void NeighborComorderH2H(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x);
    /// Tree building
    void makeTree();
    int match(int x,vector<pair<int,pair<int,int>>> &vert);
    /// Index building
    void makeIndex();
    void makeRMQ();
    void makeRMQDFS(int p, int height);
    void makeIndexDFS(int p, vector<int>& list);

    void IndexsizeH2H();  //Index size computation


    /// For PLL index construction
    vector<unordered_map<vertex,int>> Label;
    vector<unordered_map<vertex,unordered_set<vertex>>> PruningPointSet;//{(v,c),u}
    vector<unordered_map<vertex,vertex>> PruningPointSet2;//<v,u,{c}>

    void PLLIndexConstruction(int type);
    void IndexSizePLL();
    //PLL
    void PLLIndexConstruct(vector<vector<pair<vertex,int>>> &Neighbor);
    void DijksPrune1(int nodeID, vector<pair<int,int>>& vp, vector<vector<pair<vertex,int>>> &Neighbor);
    int ShortestDisQuery1(int ID1,int ID2,vector<int>& SupNode, int& d);

    //PSL
    vector<unordered_map<vertex,int>> Dhop;
    vector<bool> DvertexNew;
    void PSLIndexConstruct(vector<vector<pair<vertex,int>>> &Neighbor);
    bool DhopLableRefreshMulti2New(int step, vector<vector<pair<vertex,int>>> &Neighbor);
    void labelMultiThread2New(vector<unordered_map<vertex,int>>& newDhop, vector<int>& p,int step, vector<vector<pair<vertex,int>>> &Neighbor);
    int ShortestDisQuery2(int ID1,int ID2,vector<int>& SupNode, int& d);
    void threadDistribute(vector<vector<int>>& processID);

    //BPCL
    int BPCLBatchSize=threadnum*10;
    void BPCLIndexConstruct(vector<vector<pair<vertex,int>>> &Neighbor);
    void BatchPCLDijk(vertex nodeID, unordered_set<vertex>& setNodes, vertex hID, vector<vector<pair<vertex,int>>> &Neighbor);
    void BatchPCLDijk2(vector<vertex>& p, unordered_set<vertex>& setNodes, vertex hID, vector<vector<pair<vertex,int>>> &Neighbor);
    void PruningPointBuild(bool ifParallel, vector<vector<vertex>> & processID, vector<vector<pair<vertex,int>>> &Neighbor);
    void PPRConstruction(vertex nodeID, vector<vector<pair<vertex,int>>> &Neighbor);
    void PPRConstruction2(vector<vertex> & p, vector<vector<pair<vertex,int>>> &Neighbor);
    void ThreadDistribute(vector<vertex>& vertices, vector<vector<vertex>>& processID);

    //// Query processing
    void CorrectnessCheck(int runtimes);
    void EffiCheck(int runtimes);
    unsigned long long EffiCheckThroughput(vector<pair<int,int>>& ODpair, int runtimes, int batchInterval, double updateT);

    // For CH query processing
    int	QueryCHWP(int ID1, int ID2);

    /// For H2H query processing
    int QueryH2H(int ID1,int ID2);
    int LCAQuery(int _p, int _q);

    // For PLL query
    int QueryPLL(int ID1, int ID2);

    //Dijkstra
    int Dijkstra(int ID1, int ID2,vector<vector<pair<vertex,int>>> &Neighbor);
    void RetrievePath(int ID1, int ID2, vector<int> & prece);

    //// For Index Maintenance test
    void IndexMaintenance(int updateType, int batchNum, bool ifBatch, int batchSize);
    // For CH
    void CHincBatMT(vector<pair<pair<int,int>,pair<int,int>>>& wBatch);//CHW increase
    void CHdecBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch);//CHW decrease

    // For H2H
    void H2HdecBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch);//decrease
    void H2HincBatMT(vector<pair<pair<int,int>,pair<int,int>>>& wBatch);//increase
    void EachNodeProBDis5(int child,vector<int>& line,set<int>& vertexIDChL, map<int,int>& checkedDis);
    void eachNodeProcessIncrease1(int children, vector<int>& line, int& changelabel);

    // For PLL
    void DecreasePSL(int a, int b, int oldW, int newW, vector<vector<pair<vertex,int>>> &Neighbors,vector<unordered_map<vertex,int>> &Label);
    int PLLQuery(int ID1,int ID2,vector<unordered_map<vertex,int>> &Label);
    void IncreasePSL(int a, int b, int oldW, int newW, vector<vector<pair<vertex,int>>> &Neighbor,vector<unordered_map<vertex,int>> &Label, vector<unordered_map<vertex,unordered_set<vertex>>> &PruningPointNew, vector<unordered_map<vertex,vertex>>& PruningPointSet2);//set version with NoSupportedPair
    void CoarseUpdate(int LID,int HID, int oldW, queue<pair<int,pair<int,int>>>& WaitPro, queue<pair<int,pair<int,int>>>& WaitProP, vector<pair<int,int>>& AL1, vector<pair<int,int>>& AL2, vector<pair<int,int>>& AL2Check, vector<vector<pair<vertex,int>>> &Neighbor,vector<unordered_map<vertex,int>> &Label, bool ifDebug, int lid, int hid, vector<unordered_set<int>>& ChangedLabels);
    void RefineUpdate(queue<pair<int,pair<int,int>>>& WaitPro, queue<pair<int,pair<int,int>>>& WaitProP, vector<pair<int,int>>& AL1, vector<pair<int,int>>& AL2, vector<pair<int,int>>& AL2Check, set<tuple<vertex,vertex,vertex>>& outdatedPruning, map<pair<vertex,vertex>,vertex>& newPruningPoints, vector<unordered_map<vertex,vertex>>& PruningPointSet2, vector<vector<pair<vertex,int>>> &Neighbor,vector<unordered_map<vertex,int>> &Label,vector<unordered_map<vertex,unordered_set<vertex>>> &PruningPointNew, bool ifDebug, int lid, int hid, vector<unordered_set<int>>& ChangedLabels, set<pair<int,int>>& NoSupportedPair);//queue version
    bool PPRCheck(int curID, int hubID, vector<vector<pair<vertex,int>>> &Neighbor, vector<unordered_map<vertex,int>> &Label, vector<unordered_map<vertex,unordered_set<vertex>>> &PruningPointNew, queue<pair<int,pair<int,int>>>& WaitProPTem, vector<pair<int,int>> &AL2, map<pair<vertex,vertex>,vertex> &newPruningPoints, vector<unordered_map<vertex,vertex>>& PruningPointSet2, set<tuple<vertex,vertex,vertex>> &outdatedPruning, bool ifDebug, int lid, int hid, vector<unordered_set<int>>& ChangedLabels, set<pair<int,int>>& NoSupportedPair);//queue version
    void PPRClean(vector<vector<pair<vertex,int>>> &Neighbors, map<pair<vertex,vertex>,vertex> &newPruningPoints, set<tuple<vertex,vertex,vertex>> &outdatedPruning, bool ifDebug, int lid, int hid, vector<unordered_map<vertex,int>> &Label, vector<unordered_map<vertex,unordered_set<vertex>>>& PruningPointSet, vector<unordered_map<vertex,vertex>>& PruningPointSet2, vector<unordered_set<int>>& ChangedLabels);
    int DisQueryLower1(int ID1, int ID2, vector<vector<pair<vertex,int>>> &Neighbors,vector<unordered_map<vertex,int>> &Label);
    int DisQueryVally(int ID1, int ID2, vector<vector<pair<vertex,int>>> &Neighbors,vector<unordered_map<vertex,int>> &Label);
    pair<int,int> DisQueryVally2(int ID1, int ID2, vector<vector<pair<vertex,int>>> &Neighbors,vector<unordered_map<vertex,int>> &Label);
    int DisQueryPeak(int ID1, int ID2,vector<unordered_map<vertex,int>> &Label);
    pair<int,int> DisQueryPeak2(int ID1, int ID2,vector<unordered_map<vertex,int>> &Label);

	/// Graph Preprocessing
    void ReadGraph(string filename);
    void ReadUpdate(string filename,vector<pair<pair<int,int>,int>>& TestData);
    void ReadUpdate2(string filename,vector<pair<pair<int,int>,pair<int,int>>>& TestData);
    void ReadUpdate3(string filename,vector<pair<pair<int,int>,tuple<int,int,int>>>& TestData);
	void StainingMethod(int ID);
	void ODGene(int num, string filename);
	void UpdateGene(int num, string filename);

    void WriteOrder(string filename);
    void ReadOrder(string filename);

    vector<int> DFS_CC(vector<map<int,int>> & Edges, set<int> set_A, set<int> & set_B, int nodenum);
    vector<int> DFS_CC(vector<vector<pair<int,int>>> & Edges, set<int> set_A, set<int> & set_B, int nodenum);

};

#endif // HEADSP_H_
