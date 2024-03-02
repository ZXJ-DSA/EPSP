/*
 * Baselines.cpp
 *
 *  Created on: 16 June 2023
 *     Authors: Xinjie ZHOU
 */
#include "head.h"


///**** LCC ****/
//function for index construction of LCC
void Graph::LCCIndexConstruct(vector<vector<pair<vertex,int>>> &Neighbor){
    bool ifParallel = true;

    Label.assign(node_num,unordered_map<vertex,int>());

    PruningPointSet2.clear();
    PruningPointSet2.assign(node_num,unordered_map<vertex,vertex>());
    PruningPointSet.clear();
    PruningPointSet.assign(node_num,unordered_map<vertex,unordered_set<vertex>>());
//    PruningPointSetOrder.clear();
//    PruningPointSetOrder.assign(node_num,unordered_map<int,set<OrderComp>>());


//    BatchPrunedDijkNew(11449, 0);
//    BatchPrunedDijkNew(12516, 0);
//    exit(0);

    int ID;
    int cnt=0;
    int stepShow = ceil(node_num/100000)*1000;
    stepShow = max(stepShow,1000);

    double runT1, runT2;
    runT1=0, runT2=0;
    cout<<"Step for show: "<<stepShow<<endl;

    int batchSize = 512;
    batchSize = threadnum;
    if(batchSize>node_num){
        batchSize=node_num;
    }
//    batchSize = 3;
    cout<<"Batch size: "<<batchSize<<endl;

    vector<vector<vertex>> batches;
    vector<vertex> bNodes;

    Timer tt;
    double time=0;

    vector<vertex> vertices;

    for(int i=node_num-1;i>=0;i--){
        ID=vNodeOrder[i];

        if(CoreTag[ID] != -1){//if ID is not core vertex
            if(!bNodes.empty()){
                batches.push_back(bNodes);
                bNodes.clear();
            }
            break;
        }
        vertices.emplace_back(ID);

        if(bNodes.size()<batchSize){
            bNodes.emplace_back(ID);
            if(bNodes.size()==batchSize || i==0){
                batches.push_back(bNodes);
                bNodes.clear();
//                bNodes.emplace_back(ID);
            }
        }
    }

    int batchNum = batches.size();
    int batchShow = 1;
    if (batchShow < batchNum / 10){
        batchShow = batchNum / 10;
    }

    cout<<"Total batch number: "<<batchNum<<endl;
    cout<<"Batch gap for show: "<<batchShow<<endl;
    int lid, hid;
    lid=200291, hid=200093;

    Timer tt2;
    tt2.start();
    vertex hID = vNodeOrder[node_num-1];
    for(int b1=0;b1<batches.size();++b1){
        bNodes=batches[b1];
//        hID = bNodes[0];
        unordered_set<vertex> setNodes; setNodes.clear();
        setNodes.insert(bNodes.begin(),bNodes.end());

        vector<pair<int,int>> vp;
        // process each batch
        tt.start();
        if(ifParallel){//use multiple thread
//        if(false){

            boost::thread_group thread;

            for(int i=0;i<bNodes.size();i++){
                ID = bNodes[i];
                thread.add_thread(new boost::thread(&Graph::DijksPrune2, this, ID, boost::ref(Neighbor)));
            }
            thread.join_all();
        }else{//use single thread
            for(int i=0;i<bNodes.size();i++){
                ID = bNodes[i];
                DijksPrune2(ID, Neighbor);
//                for(int j=0;j<vp.size();j++){
//                    Label[vp[j].first].insert(make_pair(ID, vp[j].second));
//                    //cout<<vp[j].first<<" "<<vp[j].second<<endl;
//                }
            }
        }

        tt.stop();
        time+=tt.GetRuntime();
        if(b1%batchShow==0){
//        if(b1%batchShow==0 && b1<451){
            cout<<"Batch "<<b1<<" finish! Order("<<NodeOrder[bNodes[0]]<<","<<NodeOrder[bNodes[bNodes.size()-1]]<<"). "<<time<<" s."<<endl;
            time = 0;

        }

    }
    tt2.stop();
    runT1 = tt2.GetRuntime();

    vector<vector<vertex>> processID;
    ThreadDistribute(vertices, processID);

    cout<<"Batch size for label cleaning: "<<processID[0].size()<<endl;
    tt2.start();
    boost::thread_group thread;
    for(int i=0;i<processID.size();i++){
        thread.add_thread(new boost::thread(&Graph::DQ_Cleans, this, boost::ref(processID[i])));
    }
    thread.join_all();
    tt2.stop();
    runT2 = tt2.GetRuntime();
    cout<<"Time used for label construction: "<<runT1+runT2<<" s. time for label cleaning: "<<runT2<<" s."<< endl;

    tt2.start();
    cout<<"Begin to construct pruning points..."<<endl;
    PruningPointBuild(ifParallel,processID,Neighbor);
    tt2.stop();
    cout<<"Time used for pruning point construction: "<<tt2.GetRuntime()<<" s."<<endl;
}
//function of pruning Dijkstra of LCC
void Graph::DijksPrune2(int nodeID, vector<vector<pair<vertex,int>>> &Neighbor){
    benchmark::heap<2, int, int> pqueue(node_num);
    pqueue.update(nodeID,0);

    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);

    distance[nodeID]=0;
    int topNodeID, topNodeDis;

    int NNodeID,NWeigh;
    bool ifPrune;

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        closed[topNodeID]=true;
        ifPrune = false;
        int TempDis; vector<int> SupNode;

        if(NodeOrder[topNodeID] > NodeOrder[nodeID]){
            continue;
        }

        vSm[nodeID]->wait();
        if(topNodeID != nodeID){
            vSm[topNodeID]->wait();
        }

        ShortestDisQuery1(nodeID, topNodeID,SupNode,TempDis);

        vSm[nodeID]->notify();
        if(topNodeID != nodeID){
            vSm[topNodeID]->notify();
        }

        if(TempDis<=topNodeDis) {//if the peak path is longer than this vally path
            continue;
        }

        vSm[topNodeID]->wait();
        Label[topNodeID].insert(make_pair(nodeID, topNodeDis));
        vSm[topNodeID]->notify();


        //Label[topNodeID].insert(nodeID, topNodeDis);
        for(auto it=AdjaCore[topNodeID].begin();it!=AdjaCore[topNodeID].end();it++){
            NNodeID=(*it).first;
            NWeigh=(*it).second+topNodeDis;
            if(!closed[NNodeID]){
                if(distance[NNodeID]>NWeigh){
                    distance[NNodeID]=NWeigh;
                    pqueue.update(NNodeID, NWeigh);
                }
            }
        }

    }
}
//GLL
void Graph::GLLIndexConstruct(vector<vector<pair<vertex,int>>> &Neighbor){
    bool ifParallel = true;

    Label.assign(node_num,unordered_map<vertex,int>());

    PruningPointSet2.clear();
    PruningPointSet2.assign(node_num,unordered_map<vertex,vertex>());
//    PruningPointSetOrder.clear();
//    PruningPointSetOrder.assign(node_num,unordered_map<int,set<OrderComp>>());
    PruningPointSet.clear();
    PruningPointSet.assign(node_num,unordered_map<vertex,unordered_set<vertex>>());

//    BatchPrunedDijkNew(11449, 0);
//    BatchPrunedDijkNew(12516, 0);
//    exit(0);

    int ID;
    int cnt=0;
    int stepShow = ceil(node_num/100000)*1000;
    stepShow = max(stepShow,1000);

    cout<<"Step for show: "<<stepShow<<endl;

    int batchSize = 512;
    batchSize = threadnum;
    if(batchSize>node_num){
        batchSize=node_num;
    }
//    batchSize = 3;
    cout<<"Batch size: "<<batchSize<<endl;

    vector<vector<vertex>> batches;
    vector<vertex> bNodes;

    vector<vertex> vertices;

    Timer tt;
    double time=0;

    for(int i=node_num-1;i>=0;i--){
        ID=vNodeOrder[i];

        if(CoreTag[ID] != -1){//if ID is not core vertex
            if(!bNodes.empty()){
                batches.push_back(bNodes);
                bNodes.clear();
            }
            break;
        }
        vertices.emplace_back(ID);
        if(bNodes.size()<batchSize){
            bNodes.emplace_back(ID);
            if(bNodes.size()==batchSize || i==0){
                batches.push_back(bNodes);
                bNodes.clear();
//                bNodes.emplace_back(ID);
            }
        }
    }

    int batchNum = batches.size();
    int batchShow = 1;
    if (batchShow < batchNum / 10){
        batchShow = batchNum / 10;
    }

    cout<<"Total batch number: "<<batchNum<<endl;
    cout<<"Batch gap for show: "<<batchShow<<endl;
    int lid, hid;
    lid=200291, hid=200093;

    Timer tt2;
    tt2.start();
    int hID = vNodeOrder[node_num-1];
    for(int b1=0;b1<batches.size();++b1){
        bNodes=batches[b1];
//        hID = bNodes[0];
        unordered_set<vertex> setNodes; setNodes.clear();
        setNodes.insert(bNodes.begin(),bNodes.end());

        vector<vector<pair<int,int>>> localTable;
        localTable.assign(bNodes.size(),vector<pair<int,int>>());

        // process each batch
        tt.start();
        if(ifParallel){//use multiple thread
//        if(false){
            boost::thread_group thread;
            for(int i=0;i<bNodes.size();i++){
                ID = bNodes[i];
//                thread.add_thread(new boost::thread(&Graph::DijksPrune2, this, ID, boost::ref(vp)));
                thread.add_thread(new boost::thread(&Graph::DijksPrune3, this, ID,  boost::ref(localTable[i]), boost::ref(Neighbor)));
            }
            thread.join_all();


            boost::thread_group thread2;
            for(int i=0;i<bNodes.size();i++){
                ID = bNodes[i];
                thread2.add_thread(new boost::thread(&Graph::CommitLabel, this, ID,  boost::ref(localTable[i])));
            }
            thread2.join_all();
        }else{//use single thread
            for(int i=0;i<bNodes.size();i++){
                ID = bNodes[i];
//                DijksPrune2(ID);
                DijksPrune3(ID, localTable[i], Neighbor);
                CommitLabel(ID,localTable[i]);
            }

        }

        BatchClean(bNodes, setNodes, ifParallel);
//        BatchClean(bNodes, setNodes, false);

        tt.stop();
        time+=tt.GetRuntime();
        if(b1%batchShow==0){
//        if(b1%batchShow==0 && b1<451){
            cout<<"Batch "<<b1<<" finish! Order("<<NodeOrder[bNodes[0]]<<","<<NodeOrder[bNodes[bNodes.size()-1]]<<"). "<<time<<" s."<<endl;
            time = 0;

        }

    }
    tt2.stop();
    cout<<"Time used for label construction: "<<tt2.GetRuntime()<<" s."<<endl;

    vector<vector<vertex>> processID;
    processID.assign(threadnum, vector<vertex>());
    ThreadDistribute(vertices, processID);
    tt2.start();
    cout<<"Begin to construct pruning points..."<<endl;
    PruningPointBuild(ifParallel, processID, Neighbor);
    tt2.stop();
    cout<<"Time used for pruning point construction: "<<tt2.GetRuntime()<<" s."<<endl;
}

void Graph::CommitLabel(int ID, vector<pair<int,int>>& vp){
    int topNodeID;
    for(auto it=vp.begin();it!=vp.end();++it){
        topNodeID = it->first;
        vSm[topNodeID]->wait();
        Label[topNodeID].insert(make_pair(ID,it->second));
        vSm[topNodeID]->notify();
    }
}
//function of pruning Dijkstra of LCC
void Graph::DijksPrune3(int nodeID, vector<pair<int,int>>& vp, vector<vector<pair<vertex,int>>> &Neighbor){
    benchmark::heap<2, int, int> pqueue(node_num);
    pqueue.update(nodeID,0);

    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);

    distance[nodeID]=0;
    int topNodeID, topNodeDis;
    int NNodeID,NWeigh;
    bool ifPrune;

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        closed[topNodeID]=true;
        ifPrune = false;
        int TempDis; vector<int> SupNode;

        if(NodeOrder[topNodeID] > NodeOrder[nodeID]){
            continue;
        }


        ShortestDisQuery1(nodeID, topNodeID,SupNode,TempDis);


        if(TempDis<=topNodeDis) {//if the peak path is longer than this vally path
            continue;
        }

//        vSm[topNodeID]->wait();
//        Label[topNodeID].insert(make_pair(nodeID, topNodeDis));
//        vSm[topNodeID]->notify();

        vp.emplace_back(make_pair(topNodeID, topNodeDis));


        //Label[topNodeID].insert(nodeID, topNodeDis);
        for(auto it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
            NNodeID=(*it).first;
            NWeigh=(*it).second+topNodeDis;
            if(!closed[NNodeID]){
                if(distance[NNodeID]>NWeigh){
                    distance[NNodeID]=NWeigh;
                    pqueue.update(NNodeID, NWeigh);
                }
            }
        }

    }
}
//function to clean redundant labels of ID1
void Graph::DQ_Clean(vertex ID1){
    pair<int,int> peakPair;
    int dispeak, peakhub;
    int ID2, Dis;
    vector<pair<int,int>> eraseLabels;

    vSm[ID1]->wait();
    for(auto it2=Label[ID1].begin();it2!=Label[ID1].end();++it2){
        ID2=it2->first, Dis=it2->second;
        if(ID1 == ID2) continue;

        peakPair=DisQueryPeak2(ID1,ID2,Label);
        dispeak=peakPair.first, peakhub=peakPair.second;
        if(dispeak<=Dis){
            eraseLabels.emplace_back(ID1,ID2);
        }

    }
    for(auto it2=eraseLabels.begin();it2!=eraseLabels.end();++it2){
        ID1=it2->first, ID2=it2->second;
        Label[ID1].erase(ID2);
    }
    vSm[ID1]->notify();
}
//function to clean redundant labels of IDs
void Graph::DQ_Cleans(vector<vertex>& IDs) {
    for(int i=0;i<IDs.size();++i){
        DQ_Clean(IDs[i]);
    }
}
//function of clean redundant labels of ID1
void Graph::LabelClean(vertex ID1, unordered_set<vertex>& setNodes){
    pair<int,int> peakPair;
    int dispeak, peakhub;
    int ID2, Dis;
    vector<int> eraseHubs;
    vector<pair<vertex,int>> remainedLabels;

    vSm[ID1]->wait();
    for(auto it2=Label[ID1].begin();it2!=Label[ID1].end();++it2){
        ID2=it2->first, Dis=it2->second;
        if(ID1 == ID2) continue;
        if(setNodes.find(ID2) != setNodes.end()){//if found
            vSm[ID2]->wait();
            peakPair=DisQueryPeak2(ID1,ID2,Label);
            vSm[ID2]->notify();
            dispeak=peakPair.first, peakhub=peakPair.second;
            if(dispeak<=Dis){
                eraseHubs.emplace_back(ID2);
            }else{
                remainedLabels.emplace_back(ID2,Dis);
            }
        }else{
            peakPair=DisQueryPeak2(ID1,ID2,Label);
            dispeak=peakPair.first, peakhub=peakPair.second;
            if(dispeak<=Dis){
                eraseHubs.emplace_back(ID2);
            }else{
                remainedLabels.emplace_back(ID2,Dis);
            }
        }
    }
    for(auto it2=eraseHubs.begin();it2!=eraseHubs.end();++it2){
        ID2=*it2;
        Label[ID1].erase(ID2);
    }
    LabelV.Labels[ID1] = remainedLabels;
    vSm[ID1]->notify();
}
//function of cleaning redundant labels generated in current batch
void Graph::BatchClean(vector<vertex> &bNodes, unordered_set<vertex>& setNodes, bool ifParallel){
    int ID;
//    Timer tt;
//    tt.start();
    if(ifParallel){
        boost::thread_group thread;
        for(int i=0;i<bNodes.size();i++){
            ID = bNodes[i];
            thread.add_thread(new boost::thread(&Graph::LabelClean, this, ID, boost::ref(setNodes)));
        }
        thread.join_all();
    }else{
        for(auto it=bNodes.begin();it!=bNodes.end();++it){
            ID = *it;
            LabelClean(ID, setNodes);
        }
    }
//    tt.stop();
//    cout<<"Time used for label cleaning: "<<tt.GetRuntime()<<" s."<<endl;
}

double Graph::GLLConstructV(vector<vector<pair<vertex,int>>> &Neighbor){
    cout<<"Vector-based implementation. Boost-based Labeling."<<endl;
    bool ifParallel = true;
    double t1=0,t2=0;

    Label.assign(node_num,unordered_map<vertex,int>());
    LabelV.resize(node_num);

    int ID;
    int cnt=0;
    int stepShow = ceil(node_num/100000)*1000;
    stepShow = max(stepShow,1000);

    cout<<"Step for show: "<<stepShow<<endl;

    int batchSize = 512;
    batchSize = threadnum;
    if(batchSize>node_num){
        batchSize=node_num;
    }
//    batchSize = 3;
    cout<<"Batch size: "<<batchSize<<endl;

    vector<vector<vertex>> batches;
    vector<vertex> bNodes;

    vector<vertex> vertices; vertices.clear();

    Timer tt;
    double time=0;

    for(int i=node_num-1;i>=0;i--){
        ID=vNodeOrder[i];

        if(CoreTag[ID] != -1){//if ID is not core vertex
            if(!bNodes.empty()){
                batches.push_back(bNodes);
                bNodes.clear();
            }
            break;
        }
        vertices.emplace_back(ID);
        if(bNodes.size()<batchSize){
            bNodes.emplace_back(ID);
            if(bNodes.size()==batchSize || i==0){
                batches.push_back(bNodes);
                bNodes.clear();
//                bNodes.emplace_back(ID);
            }
        }
    }

    int batchNum = batches.size();
    int batchShow = 1;
    if (batchShow < batchNum / 10){
        batchShow = batchNum / 10;
    }

    cout<<"Total batch number: "<<batchNum<<endl;
    cout<<"Batch gap for show: "<<batchShow<<endl;
    int lid, hid;
    lid=200291, hid=200093;

    Timer tt2;
    tt2.start();
    int hID = vNodeOrder[node_num-1];
    for(int b1=0;b1<batches.size();++b1){
        bNodes=batches[b1];
//        hID = bNodes[0];
        unordered_set<vertex> setNodes; setNodes.clear();
        setNodes.insert(bNodes.begin(),bNodes.end());

        vector<vector<pair<vertex,int>>> localTable;
        localTable.assign(bNodes.size(),vector<pair<vertex,int>>());

        // process each batch
        tt.start();
        if(ifParallel){//use multiple thread
//        if(false){

            boost::thread_group thread;
            for(int i=0;i<bNodes.size();i++){
                ID = bNodes[i];
//                thread.add_thread(new boost::thread(&Graph::DijksPrune2, this, ID, boost::ref(vp)));
                thread.add_thread(new boost::thread(&Graph::DijksPrune4V, this, ID, boost::ref(Neighbor)));
            }
            thread.join_all();

        }else{//use single thread
            for(int i=0;i<bNodes.size();i++){
                ID = bNodes[i];
//                DijksPrune2(ID);
                DijksPrune4V(ID,  Neighbor);
//                CommitLabel(ID,localTable[i]);
            }

        }

        BatchClean(bNodes, setNodes, ifParallel);
//        BatchClean(bNodes, setNodes, false);


        tt.stop();
        time+=tt.GetRuntime();
        if(b1%batchShow==0){
//        if(b1%batchShow==0 && b1<451){
            cout<<"Batch "<<b1<<" finish! Order("<<NodeOrder[bNodes[0]]<<","<<NodeOrder[bNodes[bNodes.size()-1]]<<"). "<<time<<" s."<<endl;
            time = 0;

        }

    }
    tt2.stop();
    t1=tt2.GetRuntime();
    LabelV.sort(threadnum);
    cout<<"Time used for label construction: "<<t1<<" s."<<endl;

    //original version
//    vector<vector<int>> processID;
//    processID.assign(threadnum, vector<int>());
//    ThreadDistribute(vertices, processID);
//    tt2.start();
//    cout<<"Begin to construct pruning points..."<<endl;
//    PruningPointBuild(ifParallel, processID, Neighbor);
//    tt2.stop();
//    t2=tt2.GetRuntime();
//    cout<<"Time used for pruning point construction: "<<t2<<" s."<<endl;


    if(algoCoreU==0){
        /// boost-based based implementation, post-processing
        cout<<"Begin to construct pruning points... boost-based implementation with post-processing."<<endl;
        PPRV.resize(node_num);
        vector<vector<vertex>> ProcessID;
        ThreadDistribute(vertices,ProcessID);
        tt2.start();
        boost::thread_group thread;
        for(int i=0;i<ProcessID.size();i++){
            thread.add_thread(new boost::thread(&Graph::PPRConstructionV2, this, boost::ref(ProcessID[i]),  boost::ref(Neighbor)));
        }
        thread.join_all();
        tt2.stop();
        t2=tt2.GetRuntime();
        cout<<"Time used for pruning point construction: "<<t2<<" s."<<endl;

        cout<<"Core index construction time (labeling + PPR): "<<t1+t2<<" s."<<endl;

        PPRV.postProcess(PruningPointSet,PruningPointSet2,vertices,threadnum);
    }


    return t1+t2;
}
//vector-based implementation
void Graph::DijksPrune4V(vertex nodeID, vector<vector<pair<vertex,int>>> &Neighbor){
    benchmark::heap<2, int, int> pqueue(node_num);
    pqueue.update(nodeID,0);

//    unordered_map<vertex,int> Lh;
//    for(auto it=LabelV.Labels[nodeID].begin();it!=LabelV.Labels[nodeID].end();++it){
//        Lh.insert({it->first,it->second});
//    }

    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);

    distance[nodeID]=0;
    int topNodeID, topNodeDis;
    int NNodeID,NWeigh;
    bool ifPrune;

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        closed[topNodeID]=true;
        ifPrune = false;
        int TempDis; vector<int> SupNode;

        if(NodeOrder[topNodeID] > NodeOrder[nodeID]){
            continue;
        }


        vSm[nodeID]->wait();
        if(topNodeID != nodeID){
            vSm[topNodeID]->wait();
        }

        ShortestDisQuery1V(topNodeID, Label[nodeID], SupNode,TempDis);

        vSm[nodeID]->notify();
        if(topNodeID != nodeID){
            vSm[topNodeID]->notify();
        }


        if(TempDis<=topNodeDis) {//if the peak path is shorter than this vally path
            continue;
        }

        vSm[topNodeID]->wait();
        Label[topNodeID].insert(make_pair(nodeID, topNodeDis));
        LabelV.add(topNodeID,nodeID,topNodeDis);
        vSm[topNodeID]->notify();

//        vp.emplace_back(make_pair(topNodeID, topNodeDis));


        //Label[topNodeID].insert(nodeID, topNodeDis);
        for(auto it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
            NNodeID=(*it).first;
            NWeigh=(*it).second+topNodeDis;
            if(!closed[NNodeID]){
                if(distance[NNodeID]>NWeigh){
                    distance[NNodeID]=NWeigh;
                    pqueue.update(NNodeID, NWeigh);
                }
            }
        }

    }
}
//query by current labels
int Graph::ShortestDisQuery1V(int ID1,unordered_map<vertex,int>& Lh,vector<int>& SupNode, int& d){
    d=INF;
    int hub, dis1, dis2;

    for(auto it=LabelV.Labels[ID1].begin();it!=LabelV.Labels[ID1].end();++it){
        hub=(*it).first;
        dis1=(*it).second;
        if(Lh.find(hub)!=Lh.end()){//if found
            dis2=Lh[hub];
            if(dis1+dis2<d){
                d=dis1+dis2;
                SupNode.clear();
                SupNode.push_back(hub);
            }else if(dis1+dis2==d){
                SupNode.push_back(hub);
            }
        }
    }

    return d;
}

/// Index Maintenance

/// SDPLL algorithm

void Graph::MDEOrdering(){
    set<DegComp1> Deg;//min first
    int degree;

    //initialize E
    map<int,pair<int,int>> m;
    E.assign(node_num,m);
    for(int i=0;i<Neighbor.size();i++){
        for(int j=0;j<Neighbor[i].size();j++)
            E[i].insert(make_pair(Neighbor[i][j].first,make_pair(Neighbor[i][j].second,1)));
    }

    _DD_.assign(node_num,0); //_DD2_.assign(node_num,0);
    DD.assign(node_num,0); //DD2.assign(node_num,0);

    for(int i=0;i<node_num;i++){
        degree=Neighbor[i].size();

        if(degree > 0){//get degree
            _DD_[i]=degree;
            DD[i]=degree;
            Deg.insert(DegComp1(i));
        }else{
            cout<<"Wrong!! Degree of "<<i<<" is "<<degree<<endl;
            exit(1);
        }
    }

    int ID1,ID2;
    vNodeOrder.clear();
    //vector<bool> exist;
    existCore.assign(node_num,true);//if in the core, all vertices is originally in core
    vector<bool> change;
    change.assign(node_num,false);//whether the neighbor (degree) has changed

    int x;
    int rank=0;
    while(!Deg.empty()) {
        int x=(*Deg.begin()).x;//minimum degree first
        while(change[x]){//update the degree if it is changed
            Deg.erase(DegComp1(x));
            _DD_[x]=DD[x];
//            _DD2_[x]=DD2[x];
            Deg.insert(DegComp1(x));
            change[x]=false;
            x=(*Deg.begin()).x;
        }

        vNodeOrder.push_back(x);//least important vertex first
        Deg.erase(Deg.begin());

        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();
        for(auto it=E[x].begin();it!=E[x].end();it++){
            if(existCore[(*it).first]){
                Neigh.push_back(*it);
            }
        }

        //delete the star
        for(int i=0;i<Neigh.size();i++){
            int y=Neigh[i].first;
            deleteECore(x,y);//delete x from y's adjacency list
            change[y]=true;
        }
        //add all-pair neighbors
        for(int i=0;i<Neigh.size();i++){
            ID1=Neigh[i].first;
            for(int j=i+1;j<Neigh.size();j++){
                ID2=Neigh[j].first;
                insertECore(ID1,ID2,Neigh[i].second.first+Neigh[j].second.first);


                change[ID1]=true;
                change[ID2]=true;
            }
        }
    }

    NodeOrder.assign(node_num,-1);
    for(int k=0;k<vNodeOrder.size();k++){
        NodeOrder[vNodeOrder[k]]=k;
    }

    NodeOrder_ = NodeOrder;

//    for(int i=node_num-1;i>node_num-5;--i){
//        cout<<vNodeOrder[i]<<" "<<Neighbor[vNodeOrder[i]].size()<<endl;
//    }
}

void Graph::DegreeOrdering(){
    set<DegComp1> Deg;//min first
    int degree;

    vNodeOrder.assign(node_num,-1);
    NodeOrder.assign(node_num,-1);
    _DD_.assign(node_num,0); //_DD2_.assign(node_num,0);

    for(int i=0;i<node_num;i++){
        degree=Neighbor[i].size();

        if(degree > 0){//get degree
            _DD_[i]=degree;
            Deg.insert(DegComp1(i));
        }else{
            cout<<"Wrong!! Degree of "<<i<<" is "<<degree<<endl;
            exit(1);
        }
    }


    int ID;
    int rank=0;
    while(!Deg.empty()) {
        ID = (*Deg.begin()).x;//minimum degree first
        NodeOrder[ID] = rank;
        vNodeOrder[rank]=ID;
        Deg.erase(Deg.begin());
        rank++;
    }

    for(int i=node_num-1;i>node_num-5;--i){
        cout<<vNodeOrder[i]<<" "<<Neighbor[vNodeOrder[i]].size()<<endl;
    }

    for(int i=0;i<5;++i){
        cout<<vNodeOrder[i]<<" "<<Neighbor[vNodeOrder[i]].size()<<endl;
    }
}

// Index Construction
void Graph::BPCLConstruction() {
    cout<<"SDPLL"<<endl;

//    MDEOrdering();
    DegreeOrdering();

    CoreTag.assign(node_num,-1);
//    vSm.reserve(node_num);
//    for(int i = 0; i < node_num; i++)
//    {
//        Semaphore* s = new Semaphore(1);
//        vSm.push_back(s);
//    }
//    BPCLIndexConstruct(Neighbor);//use batch optimization

    cout<<"PSL"<<endl;
    vSm.reserve(node_num);
    for(int i = 0; i < node_num; i++)
    {
        Semaphore* s = new Semaphore(1);
        vSm.push_back(s);
    }
    PSLIndexConstruct(Neighbor);
}




// Index Update
void Graph::PLLdec(int a,int b, int oldW, int newW, vector<vector<pair<vertex,int>>> &Neighbor){
    for(int i=0;i<Neighbor[a].size();i++){
        if(Neighbor[a][i].first==b){
            Neighbor[a][i].second=newW;
            break;
        }
    }
    for(int j=0;j<Neighbor[b].size();j++){
        if(Neighbor[b][j].first==a){
            Neighbor[b][j].second=newW;
            break;
        }
    }

    benchmark::heap<2, int, int> Q0(node_num);

    int hubID;
    for(auto it1=Label[a].begin();it1!=Label[a].end();it1++){
        hubID=(*it1).first;
        Q0.update(hubID, node_num-NodeOrder[hubID]);
        //Q0.update(hubID, NodeOrder[hubID]);
    }
    for(auto it2=Label[b].begin();it2!=Label[b].end();it2++){
        hubID=(*it2).first;
        Q0.update(hubID, node_num-NodeOrder[hubID]);
        //Q0.update(hubID, NodeOrder[hubID]);
    }

    int NodeID, NodeR;
    while(!Q0.empty()){
        Q0.extract_min(NodeID, NodeR);
        if(Label[a].find(NodeID)!=Label[a].end()){
            DijkstraPrune1(NodeID, b, Label[a][NodeID]+newW, Neighbor);
        }
        if(Label[b].find(NodeID)!=Label[b].end()){
            DijkstraPrune1(NodeID, a, Label[b][NodeID]+newW, Neighbor);
        }
    }
}

void Graph::PLLinc(int u,int v, int oldW, int newW, vector<vector<pair<vertex,int>>> &Neighbor){
    vector<int> Disu, Disv;
    vector<int> DisuNew, DisvNew;
    DijkstraList(u, Disu, Neighbor);
    DijkstraList(v, Disv, Neighbor);

    for(int i=0;i<Neighbor[u].size();i++){
        if(Neighbor[u][i].first==v){
            Neighbor[u][i].second=newW;
            break;
        }
    }
    for(int j=0;j<Neighbor[v].size();j++){
        if(Neighbor[v][j].first==u){
            Neighbor[v][j].second=newW;
            break;
        }
    }

    DijkstraList(u, DisuNew, Neighbor);
    DijkstraList(v, DisvNew, Neighbor);

    //identify the possible & real affected labels
    vector<int> PAu, PAv, RAu, RAv;
    PLLaffect(u, oldW, PAu, RAu, Disu, Disv, DisvNew, Neighbor);
    PLLaffect(v, oldW, PAv, RAv, Disv, Disu, DisuNew, Neighbor);

    //identify the invalid and missing labels
    vector<int> ILu,ILv,MLu,MLv;
    PLLinvalid(u, ILu, MLu, RAu, RAv, PAv);
    PLLinvalid(v, ILv, MLv, RAv, RAu, PAu);

    //remove the invalid label
    PLLremove(u, RAu, ILu);
    PLLremove(v, RAv, ILv);

    //connected or not after edge deletion
    int disuv=DijkstraCore(u,v);
    if(disuv!=INF){//still connected
        //collect IL&ML in decreasing vertex order
        //collect RA in a set
        set<int> Affect;
        for(int i=0;i<RAu.size();i++)
            Affect.insert(RAu[i]);
        for(int i=0;i<RAv.size();i++)
            Affect.insert(RAv[i]);

        map<int,int> Rmap;
        for(int i=0;i<ILu.size();i++)
            Rmap.insert(make_pair(NodeOrder[ILu[i]],ILu[i]));
        for(int i=0;i<ILv.size();i++)
            Rmap.insert(make_pair(NodeOrder[ILv[i]],ILv[i]));
        for(int i=0;i<MLu.size();i++)
            Rmap.insert(make_pair(NodeOrder[MLu[i]],MLu[i]));
        for(int i=0;i<MLv.size();i++)
            Rmap.insert(make_pair(NodeOrder[MLv[i]],MLv[i]));

        int rOrder, rID;
        while(Rmap.size()>0){
            rOrder=(*Rmap.rbegin()).first;
            rID=(*Rmap.rbegin()).second;
            PLLPrune(rID, Affect, Neighbor);
            Rmap.erase(rOrder);
        }
    }
}

void Graph::DijkstraPrune1(int NodeID, int ID, int d, vector<vector<pair<vertex,int>>> &Neighbor){
    benchmark::heap<2, int, int> pqueue(node_num);
    pqueue.update(ID,d);

    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);

    distance[ID]=d;
    int topNodeID, topNodeDis;
    int NNodeID,NWeigh;

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        closed[topNodeID]=true;

        int TempDis=PrefixalDisQuery(NodeID, topNodeID);
        if(TempDis<=topNodeDis)
            continue;

        Label[topNodeID][NodeID]=topNodeDis;
        for(auto it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){//Neighbor
            NNodeID=(*it).first;
            NWeigh=(*it).second+topNodeDis;
            if(!closed[NNodeID]){
                if(distance[NNodeID]>NWeigh){
                    distance[NNodeID]=NWeigh;
                    pqueue.update(NNodeID, NWeigh);
                }
            }
        }
    }
}

int Graph::DijkstraList(int ID1, vector<int>& distance, vector<vector<pair<vertex,int>>> &Neighbor){
    benchmark::heap<2, int, int> pqueue(node_num);
    pqueue.update(ID1,0);

    vector<bool> closed(node_num, false);
    distance.assign(node_num, INF);
//	vector<int> prece(node_num, 0);
    distance[ID1]=0;
    int topNodeID, topNodeDis;
    int NNodeID,NWeigh;

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        closed[topNodeID]=true;

        for(auto it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
            NNodeID=(*it).first;
            NWeigh=(*it).second+topNodeDis;
            if(!closed[NNodeID]){
                if(distance[NNodeID]>NWeigh){
                    distance[NNodeID]=NWeigh;
                    pqueue.update(NNodeID, NWeigh);
                    //	prece[NNodeID]=topNodeID;
                }
            }
        }
    }
    return 0;
}

void Graph::PLLaffect(int u, int oldw, vector<int>& PAu, vector<int>& RAu, vector<int> disu, vector<int> disv, vector<int> disvNew, vector<vector<pair<vertex,int>>> &Neighbor){
    vector<bool> Flag(node_num, false);
    Flag[u]=true;
    queue<int> Q;
    Q.push(u);
    int top,r;
    while(!Q.empty()){
        top=Q.front();
        Q.pop();
        for(int nid=0; nid<Neighbor[top].size(); nid++){
            r=Neighbor[top][nid].first;
            //cout<<"r "<<r<<endl;
            if(!Flag[r]){
                if(disv[r]==disu[r]+oldw){
                    Q.push(r);
                    PAu.push_back(r);
                    //cout<<"PA push"<<endl;
                    if(disvNew[r]!=disu[r]+oldw){
                        RAu.push_back(r);
                    }
                    Flag[r]=true;
                }
            }
        }
    }
}

void Graph::PLLinvalid(int u, vector<int>& ILu, vector<int>& MLu, vector<int> RAu, vector<int> RAv, vector<int> PAv){
    int t,r,w,ID;
    for(int i=0;i<RAu.size();i++){
        t=RAu[i];
        for(int j=0;j<RAv.size();j++){
            r=RAv[j];
            if(Label[t].find(r)!=Label[t].end()){
                ILu.push_back(r);

                //find ri
                int order=0;
                int ri=-1;
                for(auto it=Label[t].begin();it!=Label[t].end();it++){
                    ID=(*it).first;
                    if(NodeOrder[ID]<NodeOrder[r]){
                        if(NodeOrder[ID]>order){
                            order=NodeOrder[ID];
                            ri=ID;
                        }
                    }
                }
                if(ri!=-1){
                    set<int> PAvmap;
                    for(int k=0;k<PAv.size();k++)
                        PAvmap.insert(PAv[k]);

                    for(int Order=NodeOrder[ri]+1;Order<NodeOrder[r];Order++){
                        w=vNodeOrder[Order];
                        if(PAvmap.find(w)!=PAvmap.end()){
                            int disWT=ShortestDisQuery(w,t);
                            if(Label[w].find(r)!=Label[w].end() && Label[w][r]+Label[t][r]==disWT){
                                MLu.push_back(w);
                            }
                        }
                    }
                }

            }
        }
    }

}

void Graph::PLLremove(int u, vector<int> RAu, vector<int> ILu){
    int s,r;
    set<int> ILuMap;
    for(int k=0;k<ILu.size();k++)
        ILuMap.insert(ILu[k]);

    for(int i=0;i<RAu.size();i++){
        s=RAu[i];
        vector<int> R;

        for(auto it=Label[s].begin();it!=Label[s].end();it++){
            r=(*it).first;
            if(ILuMap.find(r)!=ILuMap.end()){
                R.push_back(r);
            }
        }
        for(int j=0;j<R.size();j++){
            Label[s].erase(R[j]);
        }
    }
}

void Graph::PLLPrune(int r, set<int> Affect, vector<vector<pair<vertex,int>>> &Neighbor){
    benchmark::heap<2, int, int> pqueue(node_num);
    pqueue.update(r,0);

    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);

    distance[r]=0;
    int topNodeID, topNodeDis;
    int NNodeID,NWeigh;

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        closed[topNodeID]=true;

        int TempDis=ShortestDisQuery(r, topNodeID);

        if(TempDis>topNodeDis){
            /*if(Affect.find(topNodeID)!=Affect.end()){
                if(Label[topNodeID][r]!=topNodeDis)
                    cout<<"inconsistent///////////////////////////"<<endl;
            }*/
            Label[topNodeID][r]=topNodeDis;
            for(auto it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
                NNodeID=(*it).first;
                NWeigh=(*it).second+topNodeDis;
                if(!closed[NNodeID]){
                    if(distance[NNodeID]>NWeigh){
                        distance[NNodeID]=NWeigh;
                        pqueue.update(NNodeID, NWeigh);
                    }
                }
            }
        }

    }
}

int Graph::ShortestDisQuery(int ID1,int ID2){
    if(ID1==ID2) return 0;
    if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
    int d=INF;
    int hub, dis1, dis2;
    for(auto it=Label[ID1].begin();it!=Label[ID1].end();it++){
        hub=(*it).first;
        dis1=(*it).second;
        if(Label[ID2].find(hub)!=Label[ID2].end()){
            dis2=Label[ID2][hub];
            if(dis1+dis2<d){
                d=dis1+dis2;
                //cout<<"details "<<hub<<" "<<dis1<<" "<<dis2<<" "<<d<<endl;
            }
        }
    }
    return d;
}

int Graph::PrefixalDisQuery(int ID1, int ID2){
    int d=INF;
    int hub, dis1, dis2;
    for(auto it=Label[ID1].begin();it!=Label[ID1].end();it++){
        hub=(*it).first;
        dis1=(*it).second;
        if(NodeOrder[hub]<=NodeOrder[ID1] && Label[ID2].find(hub)!=Label[ID2].end()){
            dis2=Label[ID2][hub];
            if(dis1+dis2<d){
                d=dis1+dis2;
                //cout<<"details "<<hub<<" "<<dis1<<" "<<dis2<<" "<<d<<endl;
            }
        }
    }
    return d;
}
