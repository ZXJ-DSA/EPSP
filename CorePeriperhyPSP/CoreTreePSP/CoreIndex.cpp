/*
 * CoreIndex.cpp
 *
 *  Created on: 16 June 2023
 *     Authors: Xinjie ZHOU
 */
#include "head.h"

vector<int> NodeOrder_;//nodeID order

//// Index Construction
//Function of constructing core index
double Graph::Construct_core(int indexType) {
    double rT=0;
    Timer tt;
    //construct the index of Core
    cout<<"Begin Core's Index Construction..."<<endl;
    if(indexType == 2){
        //****************PLL construction***************************
//        this->indexType=2;
        cout<<"PLL"<<endl;
//        PLLIndexConstruct();
        rT = PLLConstructV(AdjaCore);
        //***********************************************************
    }else if(indexType == 3){
        //****************PSL construction***************************
//        this->indexType=3;
        cout<<"PSL"<<endl;
        vSm.reserve(node_num);
        for(int i = 0; i < node_num; i++)
        {
            Semaphore* s = new Semaphore(1);
            vSm.push_back(s);
        }
//        PSLIndexConstruct(AdjaCore);
        rT = PSLConstructV(AdjaCore);
        //***********************************************************
    }else if(indexType == 4){
//        this->indexType = 4;
        cout<<"GLL"<<endl;
        vSm.reserve(node_num);
        for(int i = 0; i < node_num; i++)
        {
            Semaphore* s = new Semaphore(1);
            vSm.push_back(s);
        }
//        LCCIndexConstruct();//LCC
//        GLLIndexConstruct(AdjaCore);//GLL
        rT = GLLConstructV(AdjaCore);
    }else if(indexType == 1){
//        this->indexType=1;
        cout<<"PCL"<<endl;
        vSm.reserve(node_num);
        for(int i = 0; i < node_num; i++)
        {
            Semaphore* s = new Semaphore(1);
            vSm.push_back(s);
        }
//        PCLIndexConstruct(AdjaCore);//original version
        rT = PCLConstructV(AdjaCore);
    }else if(indexType == 0){
//        this->indexType=0;
        cout<<"BPCL"<<endl;
        vSm.reserve(node_num);
        for(int i = 0; i < node_num; i++)
        {
            Semaphore* s = new Semaphore(1);
            vSm.push_back(s);
        }
        tt.start();
        BPCLIndexConstruct(AdjaCore);//use batch optimization
        tt.stop();
        rT=tt.GetRuntime();
//        rT = BPCLConstructV(AdjaCore);
    }else if(indexType == 5){
        ReadCoreIndex(graphfile);//read PLL
    }

//    CorrectnessCheckCore();
    return rT;
}

///**** PLL ****/
//function for PLL construction
void Graph::PLLIndexConstruct(){

    Label.assign(node_num,unordered_map<vertex,int>());

//    PruningPointVector.clear();
//    PruningPointVector.assign(node_num,vector<pair<int,int>>());

    PruningPointSet.clear();
    PruningPointSet.assign(node_num,unordered_map<vertex,unordered_set<vertex>>());

    PruningPointSet2.clear();
    PruningPointSet2.assign(node_num,unordered_map<vertex,vertex>());
//    PruningPointSetOrder.clear();
//    PruningPointSetOrder.assign(node_num,unordered_map<int,set<OrderComp>>());

    int ID;
    int cnt=0;
    int stepShow = ceil(node_num/100000)*10000;
    stepShow = max(stepShow,1000);
    cout<<"Step for show: "<<stepShow<<endl;

    double time = 0;
    for(int i=node_num-1;i>=0;i--){
        ID=vNodeOrder[i];

//        if(CoreVertex.find(ID) == CoreVertex.end()){//if not found
//            cout<<"Finished."<<endl;
//            break;
//        }
        if(CoreTag[ID] != -1){//if not found
            break;
//            continue;
        }
        Timer tt1;
        tt1.start();
        vector<pair<int,int>> vp;
        DijksPrune1(ID,vp);
        tt1.stop();
        time+=tt1.GetRuntime();
        if(cnt%stepShow==0){
            cout<<"Node "<<cnt<<": "<<ID<<" ; vp.size: "<<vp.size()<<". Time: "<<time<<" s."<<endl;
            time = 0;
            //cout<<"ID "<<ID<<" vp.size "<<vp.size()<<endl;
        }
        cnt+=1;
//        for(int j=0;j<vp.size();j++){
//            Label[vp[j].first].insert(make_pair(ID, vp[j].second));
//            //cout<<vp[j].first<<" "<<vp[j].second<<endl;
//        }
    }
}
//original version
void Graph::DijksPrune1(int nodeID, vector<pair<int,int>>& vp){
    benchmark::heap<2, int, int> pqueue(node_num);
    pqueue.update(nodeID,0);

    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);

    distance[nodeID]=0;
    int topNodeID, topNodeDis;
    int NNodeID,NWeigh;

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        closed[topNodeID]=true;

        int TempDis; vector<int> SupNode;
        ShortestDisQuery1(nodeID, topNodeID,SupNode,TempDis);//the vertex with NodeOrder[topNodeID] > NodeOrder[nodeID] is implicitly continued.
        if(TempDis<=topNodeDis){//if dispeak <= disvally

            if(NodeOrder[topNodeID]<NodeOrder[nodeID]){
                for(int k=0;k<SupNode.size();k++){
                    int supn=SupNode[k];
                    unordered_map<int,vector<int>> map1;
//                    PruningPointVector[topNodeID].emplace_back(supn,nodeID);
//                    PruningPointVector[nodeID].emplace_back(supn,topNodeID);
                    
//                    PruningPointSetOrder[topNodeID][supn].insert(nodeID);
//                    PruningPointSetOrder[nodeID][supn].insert(topNodeID);
                    PruningPointSet[topNodeID][supn].insert(nodeID);
                    PruningPointSet[nodeID][supn].insert(topNodeID);
                    PruningPointSet2[topNodeID][nodeID]=supn;
                    PruningPointSet2[nodeID][topNodeID]=supn;

                    //cout<<"pruning point "<<topNodeID<<" "<<nodeID<<" "<<supn<<endl;
                }
            }
            continue;
        }

        Label[topNodeID].insert({nodeID, topNodeDis});
//        vp.push_back(make_pair(topNodeID,topNodeDis));
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
//vector-based implementation
double Graph::PLLConstructV(vector<vector<pair<vertex,int>>>& Neighbor){
    double t1;
    cout<<"Vector-based construction."<<endl;

    LabelV.resize(node_num);

    PPRV.resize(node_num);

    int ID;
    int cnt=0;
    int stepShow = ceil(node_num/10000)*2000;
    stepShow = max(stepShow,1000);
    cout<<"Step for show: "<<stepShow<<endl;
    vector<vertex> vertices;

    Timer tt1;
    tt1.start();
    Timer tt;
    tt.start();
    for(int i=node_num-1;i>=0;i--){//start from the least important vertex

        ID=vNodeOrder[i];
        if(CoreTag[ID]!=-1){
            break;
        }

        vertices.emplace_back(ID);
        DijksPrune1V(ID,Neighbor);
        if(cnt%stepShow==0){
            tt.stop();
            cout<<"Node "<<cnt<<": "<<ID<<" ; time: "<<tt.GetRuntime()<<" s."<<endl;
            tt.start();
            //cout<<"ID "<<ID<<" vp.size "<<vp.size()<<endl;
        }
        cnt+=1;
    }

    LabelV.sort(threadnum);
    tt1.stop();
    t1=tt1.GetRuntime();
    cout<<"Core index construction time (labeling + PPR): "<<t1<<" s."<<endl;
    LabelV.postProcess(Label,threadnum);
    PPRV.postProcess(PruningPointSet,PruningPointSet2,vertices,threadnum);
    return t1;
}
//vector-based implementation
void Graph::DijksPrune1V(int nodeID, vector<vector<pair<vertex,int>>>& Neighbor){
    benchmark::heap<2, int, int> pqueue(node_num);
    pqueue.update(nodeID,0);

    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);

    unordered_map<vertex,int> Lh;
    for(auto it=LabelV.Labels[nodeID].begin();it!=LabelV.Labels[nodeID].end();++it){
        Lh.insert({it->first,it->second});
    }

    distance[nodeID]=0;
    int topNodeID, topNodeDis;
    int NNodeID,NWeigh;

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        closed[topNodeID]=true;

        int TempDis; vector<int> SupNode;
//        PLLDisQuery1(nodeID, topNodeID,SupNode,TempDis);
        TempDis = PLLDisQuery1V(topNodeID,Lh,SupNode);

        if(TempDis<=topNodeDis){
            //for index update
            if(NodeOrder[topNodeID]<NodeOrder[nodeID]){
                for(int k=0;k<SupNode.size();k++){
                    int supn=SupNode[k];

                    PPRV.add(nodeID,supn,topNodeID);
                }
            }
            continue;
        }

        LabelV.add(topNodeID,nodeID,topNodeDis);
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
//vector-based implementation
int Graph::PLLDisQuery1V(int ID1,unordered_map<vertex,int>& Lh,vector<int>& SupNode){
    int d=INF;

    unsigned int hub;
    unsigned long long dis1, dis2;
    for(auto it=LabelV.Labels[ID1].begin();it!=LabelV.Labels[ID1].end();it++){
        hub=(*it).first;
        dis1=(*it).second;
        if(Lh.find(hub)!=Lh.end()){
            dis2=Lh[hub];
            if(dis1+dis2<d){
                d=dis1+dis2;
                SupNode.clear();
                SupNode.push_back(hub);
            }
            else if(dis1+dis2==d){
                SupNode.push_back(hub);
            }
        }
    }

    return d;
}

///**** PSL ****/
//function for PSL construction
void Graph::PSLIndexConstruct(vector<vector<pair<vertex,int>>> &Neighbor){
    cout<<"Hash table-based implementation."<<endl;
    Label.assign(node_num, unordered_map<vertex,int>());
    Dhop.assign(node_num, unordered_map<vertex,int>());

    PruningPointSet.clear();
    PruningPointSet.assign(node_num,unordered_map<vertex,unordered_set<vertex>>());

//    PruningPointSetOrder.clear();
//    PruningPointSetOrder.assign(node_num,unordered_map<int,set<OrderComp>>());
    PruningPointSet2.clear();
    PruningPointSet2.assign(node_num,unordered_map<vertex,vertex>());

    RedundantLabels.clear();
    RedundantLabels.assign(node_num,unordered_set<int>());

//	PruningPointStepNew.clear();
//	unordered_map<int,unordered_map<int,int>> map2;
//	PruningPointStepNew.assign(node_num, map2);
//	NoSupportedPair.clear();

    DvertexNew.assign(node_num, true);

    for(int i=0;i<node_num;i++){
        Label[i].insert(make_pair(i,0));
        Dhop[i].insert(make_pair(i,0));
        //Dvectex.push_back(i);
    }
    cout<<"step 0 finish!"<<endl;
    //LabelStep.push_back(Dhop);

    bool flag=true;
    int step=1;
    while(flag){

//		LabelStep.push_back(Dhop);
        flag=DhopLableRefreshMulti2New(step+1, Neighbor);
        cout<<"step "<<step<<" finish!"<<endl;
        step+=1;
    }
    cout<<"Index finish construction"<<endl;
    Dhop.clear();
    DvertexNew.clear();
}
//function for each step
bool Graph::DhopLableRefreshMulti2New(int step, vector<vector<pair<vertex,int>>> &Neighbor){
    bool flag=false;
    vector<unordered_map<vertex,int>> newDhop;
    newDhop.assign(node_num,unordered_map<vertex,int>());

    if(true){//use multiple thread
//    if(false){
        boost::thread_group thread;
        //int stepthr=Dvectex.size()/threadnum; //cout<<" Dvectex size "<<Dvectex.size()<<endl;
        vector<vector<int>> ProcessID;
        vector<int> vvv; ProcessID.assign(threadnum,vvv);
        threadDistribute(ProcessID);//distribute the pending-for-check vertices to different threads

        for(int i=0;i<ProcessID.size();i++){
            vector<int> p=ProcessID[i];
            thread.add_thread(new boost::thread(&Graph::labelMultiThread2New, this, boost::ref(newDhop), p,step, boost::ref(Neighbor)));
        }
        thread.join_all();
    }else{//use single thread
        vector<vector<int>> ProcessID;
        vector<int> vvv; ProcessID.assign(threadnum,vvv);
        threadDistribute(ProcessID);

        for(int i=0;i<ProcessID.size();i++){
            vector<int> p=ProcessID[i];
            labelMultiThread2New(newDhop,p,step, Neighbor);
        }
    }


    //cout<<"lllllllllllllllllllllllllllllll"<<endl;


    int zerolabel=0;
    Dhop.assign(newDhop.begin(), newDhop.end());
    Dvectex.clear(); set<int> Dset;
    DvertexNew.assign(node_num,false);
    for(int nodeID=0;nodeID<node_num;nodeID++){
        if(Dhop[nodeID].size()>0){//check labels of the last step, if nodeID gains new label in the last step
            flag=true;
            for(auto it=Dhop[nodeID].begin();it!=Dhop[nodeID].end();it++){
                if(Label[nodeID].find((*it).first)!=Label[nodeID].end()){//if found
                    Label[nodeID][(*it).first]=(*it).second;
                }else{//if not found
                    Label[nodeID].insert(*it);
                }
            }
            for(int i=0;i<Neighbor[nodeID].size();i++){//propagate label to Neighbors
                if(Dset.find(Neighbor[nodeID][i].first)==Dset.end()){//if not found
                    Dset.insert(Neighbor[nodeID][i].first);
                    //Dvectex.push_back(Neighbors[nodeID][i].first);
                    DvertexNew[Neighbor[nodeID][i].first]=true;//identify the neighbors for checking in next step
                }
            }
        }else if(Dhop[nodeID].size()==0)
            zerolabel+=1;
    }

    //cout<<"zero "<<zerolabel<<" Vertex to change "<<Dvectex.size()<<endl;
    return flag;
}
//function for propagation
void Graph::labelMultiThread2New(vector<unordered_map<vertex,int>>& newDhop, vector<int>& p,int step, vector<vector<pair<vertex,int>>> &Neighbor){
    for(int i=0;i<p.size();i++){
        int nodeID=p[i];
        unordered_map<vertex,int> Dhop0; Dhop0.clear();
        int neighID, neighW;
        for(int Index=0;Index<Neighbor[nodeID].size();Index++){
            neighID=Neighbor[nodeID][Index].first;
            //cout<<"neighID "<<neighID<<" Dhop.size "<<Dhop[neighID].size()<<endl;
            if(Dhop[neighID].size()>0){
                neighW=Neighbor[nodeID][Index].second;

                auto it=Dhop[neighID].begin();
                int hub, dis, d;
                for(;it!=Dhop[neighID].end();it++){
                    hub=(*it).first; dis=(*it).second; d=neighW+dis;
                    if(NodeOrder[hub]>NodeOrder[nodeID]){//if r(hub) > r(nodeID), i.e., r(w) > r(v)
                        int TempDis; vector<int> SupNode;

                        ShortestDisQuery1(nodeID, hub,SupNode,TempDis); //query by label, TempDis is the result

                        if(TempDis>d){//if existing peak path is longer than d
                            if(Dhop0.find(hub)!=Dhop0.end()){//if found
                                if(Dhop0[hub]>d)
                                    Dhop0[hub]=d;
                            }else{
                                Dhop0.insert(make_pair(hub,d));
                            }
//                            cout<<"L("<<nodeID<<","<<hub<<"): "<<TempDis<<" "<<d<<" "<<SupNode.size()<< endl;
                        }
                        else{//if the existing peak path is shorter than d. although d is not the shortest distance, we still need to record it. if Q(hub,nodeID,L<h) <= e(nodeID,neighID) + Lh-1(neighID)[hub], i.e., Q(w,v,L<h) <= e(v,u) + Lh-1(u)[w]
                            for(int k=0;k<SupNode.size();k++){
                                int supn=SupNode[k];

                                if(supn!=nodeID && supn!=hub){
                                    vSm[nodeID]->wait();
                                    PruningPointSet[nodeID][supn].insert(hub);
                                    PruningPointSet2[nodeID][hub]=supn;
//                                    PruningPointSetOrder[nodeID][supn].insert(hub);
                                    vSm[nodeID]->notify();

                                    vSm[hub]->wait();
                                    PruningPointSet[hub][supn].insert(nodeID);
                                    PruningPointSet2[hub][nodeID]=supn;
//                                    PruningPointSetOrder[hub][supn].insert(nodeID);
                                    vSm[hub]->notify();
                                }
                            }
                        }

                    }
                }
            }
        }
        newDhop[nodeID]=Dhop0;
    }
    //cout<<"one thread finish running!"<<endl;
}
//query by current labels
int Graph::ShortestDisQuery1(int ID1,int ID2,vector<int>& SupNode, int& d){
    d=INF;
    int hub, dis1, dis2;

    for(auto it=Label[ID1].begin();it!=Label[ID1].end();it++){
        hub=(*it).first;
        dis1=(*it).second;
        if(Label[ID2].find(hub)!=Label[ID2].end()){
            dis2=Label[ID2][hub];
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
//function of querying by current label, query peak
int Graph::ShortestDisQuery2(int ID1,int ID2,vector<int>& SupNode, int& d){
    d=INF;

    int hub, dis1, dis2;
    int finalHub=-1;

    for(auto it=Label[ID1].begin();it!=Label[ID1].end();it++){
        hub=(*it).first;
        dis1=(*it).second;
        if(hub!=ID2 && Label[ID2].find(hub)!=Label[ID2].end()){
            dis2=Label[ID2][hub];
            if(dis1+dis2<d){
                finalHub=hub;
                d=dis1+dis2;
                SupNode.clear();
                SupNode.push_back(hub);
            }else if(dis1+dis2==d){
                SupNode.push_back(hub);
            }
        }
    }

    return finalHub;
}
//thread allocation
void Graph::threadDistribute(vector<vector<int>>& processID){
    int ID;
    int cnt=0;
    int threadOrder;
    int a;
    for(int r=0;r<node_num;r++){
        ID=vNodeOrder[r];
        if(DvertexNew[ID]){//whether to
            a=cnt%(2*threadnum);
            if(a>=threadnum){
                threadOrder=2*threadnum-1-a;
            }else{
                threadOrder=a;
            }
            //cout<<"a "<<a<<" threadOrder "<<threadOrder<<endl;
            processID[threadOrder].push_back(ID);
            cnt+=1;
        }
    }
}

// vector-based label
double Graph::PSLConstructV(vector<vector<pair<vertex,int>>> &Neighbor){
    double t1=0,t2=0;
    cout<<"Vector-based implementation.  Boost-based labeling."<<endl;
    unordered_map<int,int> map0; map0.clear();
//    Label.assign(node_num, map0);
//    Dhop.assign(node_num, map0);
    LabelV.resize(node_num);
    PPRV.resize(node_num);
    DhopV.assign(node_num,vector<pair<vertex,int>>());

    DvertexNew.assign(node_num, true);
    vector<vertex> vertices;

    Timer tt1;
    tt1.start();
//    #pragma omp parallel for num_threads(threadnum)
    for(int i=node_num-1;i>=0;i--){
        int id=vNodeOrder[i];
        if(CoreTag[id]!=-1){
            break;
        }
        LabelV.add(id,id,0);
        DhopV[id].emplace_back(id,0);
        vertices.emplace_back(id);
    }

    //LabelStep.push_back(Dhop);
    Timer tt;
    bool flag=true;
    int step=0;
    cout<<"step "<<step<<" finish!"<<endl;
    step+=1;
    while(flag){
//		LabelStep.push_back(Dhop);
        tt.start();
        flag=DhopLableRefreshStepV(step+1, Neighbor);
        tt.stop();
        cout<<"step "<<step<<" finish! "<<tt.GetRuntime()<<" s."<< endl;
        step+=1;
    }
    LabelV.sort(threadnum);
    tt1.stop();
    t1=tt1.GetRuntime();
//    cout<<"Time for labeling construction: "<<t1<<" s."<<endl;
    cout<<"Core index construction time (labeling + PPR): "<<t1<<" s."<<endl;
    LabelV.postProcess(Label,threadnum);
    PPRV.postProcess(PruningPointSet,PruningPointSet2,vertices,threadnum);
    return t1;
}
// vector-based implementation
bool Graph::DhopLableRefreshStepV(int step, vector<vector<pair<vertex,int>>> &Neighbor){
    bool flag=false;

    vector<vector<pair<vertex,int>>> newDhopV;
    newDhopV.assign(node_num,vector<pair<vertex,int>>());

    /// boost-based implementation
    vector<int> ProcessID;
    for(int r=node_num-1;r>=0;r--){
        int ID=vNodeOrder[r];
        if(DvertexNew[ID] && CoreTag[ID]==-1){
            ProcessID.push_back(ID);
        }
    }
    vector<vector<int>> ProcessIDs;
    ProcessIDs.assign(threadnum,vector<int>() );
    for(int i=0;i<ProcessID.size();++i){
        int p=i%threadnum;
        ProcessIDs[p].push_back(ProcessID[i]);
    }
    boost::thread_group thread;
    for(int i=0;i<ProcessIDs.size();i++){
        thread.add_thread(new boost::thread(&Graph::labelMultiThread2NewV2, this, boost::ref(newDhopV), boost::ref(ProcessIDs[i]), boost::ref(Neighbor)));
    }
    thread.join_all();

//    /// openmp-based implementation
//    vector<int> ProcessID;
//    for(int r=node_num-1;r>=0;r--){
//        int ID=vNodeOrder[r];
//        if(DvertexNew[ID] && CoreTag[ID]==-1){
//            ProcessID.push_back(ID);
//        }
//    }
//    #pragma omp parallel for num_threads(threadnum) schedule (dynamic, threadnum)
//    for (int i = 0; i < ProcessID.size(); ++i) {
//        int ID=ProcessID[i];
//        labelMultiThread2NewV(newDhopV,ID,Neighbor);
//    }


    DhopV.assign(newDhopV.begin(), newDhopV.end());
    DvertexNew.assign(node_num,false);

    int hub;
    for(int nodeID=0;nodeID<node_num;nodeID++){
        if(DhopV[nodeID].size()>0){
            flag=true;
            unordered_map<int,int> imap; imap.clear();
            for(int i=0;i<LabelV.Labels[nodeID].size();++i){
                hub=LabelV.Labels[nodeID][i].first;
                imap.insert({hub,i});
            }
            for(auto it=DhopV[nodeID].begin();it!=DhopV[nodeID].end();++it){
                hub=it->first;
                if(imap.find(hub) != imap.end()){//if found
                    LabelV.Labels[nodeID][imap[hub]].second=it->second;
                }else{
                    LabelV.add(nodeID,it->first,it->second);
                }

            }
            for(int i=0;i<Neighbor[nodeID].size();i++){//Neighbors
                int neiV = Neighbor[nodeID][i].first;
                if(neiV<0 || neiV >= node_num){
                    cout<<"Wrong "<<neiV<<endl;
                }

                if(!DvertexNew[neiV]){
                    DvertexNew[neiV]=true;
                }
            }
        }
    }

    //cout<<"zero "<<zerolabel<<" Vertex to change "<<Dvectex.size()<<endl;
    return flag;
}
//vector-based implementation
void Graph::labelMultiThread2NewV(vector<vector<pair<vertex,int>>>& newDhop, int nodeID, vector<vector<pair<vertex,int>>> &Neighbor){
//	sm->wait();
    unordered_map<vertex,int> Lh;
    for(auto it=LabelV.Labels[nodeID].begin();it!=LabelV.Labels[nodeID].end();++it){
        Lh.insert({it->first,it->second});
    }

    unordered_map<int,int> Dhop0; Dhop0.clear();
    int neighID, neighW;

    for(int Index=0;Index<Neighbor[nodeID].size();Index++){//Neighbors
        neighID=Neighbor[nodeID][Index].first;
        //cout<<"neighID "<<neighID<<" Dhop.size "<<Dhop[neighID].size()<<endl;
        if(DhopV[neighID].size()>0){//if neighbor is active
            neighW=Neighbor[nodeID][Index].second;

            int hub, dis, d;
            for(auto it=DhopV[neighID].begin();it!=DhopV[neighID].end();++it){
                hub=(*it).first; dis=(*it).second; d=neighW+dis;
                if(NodeOrder[hub]>NodeOrder[nodeID]){
                    int TempDis; vector<int> SupNode;
//                        ShortestDisQuery1(nodeID, hub,SupNode,TempDis); //cout<<"nodeID "<<nodeID<<" TempDis "<<TempDis<<" d "<<d<<endl;
//                    TempDis = ShortestDisQuery2(nodeID, hub);
                    TempDis = PLLDisQuery1V(hub, Lh, SupNode);
                    if(TempDis>d){
                        if(Dhop0.find(hub)!=Dhop0.end()){
                            if(Dhop0[hub]>d) Dhop0[hub]=d;
                        }else{
                            Dhop0.insert(make_pair(hub,d));
                        }
                    }
                    else{
                        for(int k=0;k<SupNode.size();k++){
                            int supn=SupNode[k];

                            if(supn!=nodeID && supn!=hub){
                                PPRV.add(nodeID,supn,hub);
                            }
                        }
                    }

                }
            }
        }
    }


    for(auto it=Dhop0.begin();it!=Dhop0.end();++it){
        newDhop[nodeID].emplace_back(it->first,it->second);
    }

}

void Graph::labelMultiThread2NewV2(vector<vector<pair<vertex,int>>>& newDhop, vector<int>& p, vector<vector<pair<vertex,int>>> &Neighbor){
    for(int i=0;i<p.size();++i){
        int ID=p[i];
        labelMultiThread2NewV(newDhop,ID,Neighbor);
    }
}

///**** PCL ****/
//function for index construction of batch PLL
void Graph::PCLIndexConstruct(vector<vector<pair<vertex,int>>> &Neighbor){
    bool ifParallel = true;
    Label.assign(node_num,unordered_map<vertex,int>());

    PruningPointSet2.clear();
    PruningPointSet2.assign(node_num,unordered_map<vertex,vertex>());
//    PruningPointSetOrder.clear();
//    PruningPointSetOrder.assign(node_num,unordered_map<int,set<OrderComp>>());
    PruningPointSet.clear();
    PruningPointSet.assign(node_num,unordered_map<vertex,unordered_set<vertex>>());

//    PruningPointVector.clear();
//    PruningPointVector.assign(node_num,vector<pair<int,int>>());


    vertex ID;
    vector<vertex> vertices;

    Timer tt;
    double time=0;

    for(int i=node_num-1;i>=0;i--){
        ID=vNodeOrder[i];

        if(CoreTag[ID] != -1){//if ID is not core vertex
            break;
        }

        vertices.emplace_back(ID);
    }

    vector<vector<vertex>> processID;
    ThreadDistribute(vertices, processID);

    cout<<"Batch size: "<<processID[0].size()<<endl;

    tt.start();
    if(ifParallel){//use multiple thread
//        if(false){
        boost::thread_group thread;
        for(int i=0;i<processID.size();i++){
            thread.add_thread(new boost::thread(&Graph::PCLDijks, this, processID[i], boost::ref(Neighbor)));

        }
        thread.join_all();
    }else{//use single thread
        for(int i=0;i<processID.size();i++){
            PCLDijks(processID[i], Neighbor);
        }
    }

    tt.stop();
    cout<<"Time used for label construction: "<<tt.GetRuntime()<<" s."<<endl;

    tt.start();
    cout<<"Begin to construct pruning points..."<<endl;
    PruningPointBuild(ifParallel,processID, Neighbor);
//    PruningPointBuild(false);
    tt.stop();
    cout<<"Time used for pruning point construction: "<<tt.GetRuntime()<<" s."<<endl;
}
//function of original Dijkstra of PCL, 2023-06-16, pure search-based
void Graph::PCLDijk(vertex nodeID, vector<vector<pair<vertex,int>>> &Neighbor){
    benchmark::heap<2, int, int> pqueue(node_num);
    pqueue.update(nodeID,0);

    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);
    vector<int> peak(node_num,-1);//record the peak hub from nodeID to topNodeID
    vector<int> pre(node_num,-1);//record the predecessor vertex in the shortest path
    vector<bool> pruned(node_num,false);//record the PPR pruning information

    pre[nodeID]=nodeID;
    peak[nodeID]=nodeID;
    distance[nodeID]=0;
    int topNodeID, topNodeDis;
    int NNodeID,NWeigh;
    bool ifPrune;

    bool ifDebug = false;
//    ifDebug = true;

    int lid, hid;
    lid=181451, hid=190161;
    lid=12516, hid=11449;

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        closed[topNodeID]=true;
        if(NodeOrder[peak[topNodeID]] < NodeOrder[peak[pre[topNodeID]]]) {
            peak[topNodeID] = peak[pre[topNodeID]];
        }
        if(pruned[pre[topNodeID]]){
            pruned[topNodeID]= true;
        }
        ifPrune = false;
        int TempDis; vector<int> SupNode;

        //Query of dispeak
        if(NodeOrder[nodeID] <= NodeOrder[topNodeID] && NodeOrder[peak[topNodeID]] <= NodeOrder[topNodeID]){//if valley path
//            vSm[nodeID]->wait();
            Label[nodeID].insert(make_pair(topNodeID, topNodeDis));
//            vSm[nodeID]->notify();
        }

        //Label[topNodeID].insert(nodeID, topNodeDis);
        for(auto it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
            NNodeID=(*it).first;
            NWeigh=(*it).second+topNodeDis;

            if(!closed[NNodeID]){
                if(distance[NNodeID]>NWeigh){

                    pre[NNodeID]=topNodeID;
                    distance[NNodeID]=NWeigh;
                    pqueue.update(NNodeID, NWeigh);

                    if(peak[NNodeID] == -1){
                        peak[NNodeID] = NNodeID;
                    }

                }else if(distance[NNodeID]==NWeigh){//deal with equal distance scenario

                    if(NodeOrder[peak[topNodeID]] > NodeOrder[peak[pre[NNodeID]]]){//record the highest vertex among all shortest paths
//                        cout<<"change pre. L("<<nodeID<<","<<NNodeID<<"): "<<topNodeID<<"("<<NodeOrder[topNodeID]<<") "<<pre[NNodeID]<<"("<<NodeOrder[pre[NNodeID]]<<")"<<endl;
                        pre[NNodeID]=topNodeID;
                    }

                }
            }
//            if(NNodeID==lid){
//                cout<<topNodeID<< " "<<NNodeID<<" "<<pre[NNodeID]<<"("<<NodeOrder[pre[NNodeID]]<<")"<<endl;
//            }
        }
    }
}
//function of original Dijkstra of PCL, 2023-06-16, pure search-based
void Graph::PCLDijks(vector<vertex> & IDs, vector<vector<pair<vertex,int>>> &Neighbor){
    for(int i=0;i<IDs.size();++i){
        PCLDijk(IDs[i], Neighbor);
    }
}
//function for index construction of batch PLL, vector-based implementation
double Graph::PCLConstructV(vector<vector<pair<vertex,int>>> &Neighbor){
    cout<<"Vector-based implementation of PCL.  Boost-based labeling."<<endl;
    bool ifParallel = true;
    double t1,t2;

    LabelV.resize(node_num);


    vertex ID;
    vector<vertex> vertices;

    Timer tt;
    double time=0;

    for(int i=node_num-1;i>=0;i--){
        ID=vNodeOrder[i];

        if(CoreTag[ID] != -1){//if ID is not core vertex
            break;
        }
        vertices.emplace_back(ID);
    }

    vector<vector<vertex>> ProcessID;
    ThreadDistribute(vertices, ProcessID);

    cout<<"Batch size: "<<ProcessID[0].size()<<endl;

    /// boost-based implementation, slower than OpenMP
    tt.start();
    boost::thread_group thread1;
    for(int i=0;i<ProcessID.size();i++){
        thread1.add_thread(new boost::thread(&Graph::PCLDijkV2, this, boost::ref(ProcessID[i]), boost::ref(Neighbor)));
    }
    thread1.join_all();
    LabelV.sort(threadnum);
    tt.stop();
    t1=tt.GetRuntime();
    cout<<"Time used for boost-based label construction: "<<t1<<" s."<<endl;


//    /// openmp-based implementation
//    tt.start();
//#pragma omp parallel for num_threads(threadnum) schedule (dynamic, threadnum)
//    for(int i=0;i<vertices.size();++i){
//        int ID=vertices[i];
//        PCLDijkV(ID,Neighbor);
//    }
//    LabelV.sort(threadnum);
//    tt.stop();
//    t1=tt.GetRuntime();
//    cout<<"Time used for openmp-based label construction: "<<t1<<" s."<<endl;


/// boost-based based implementation, post-processing
    cout<<"Begin to construct pruning points... boost-based implementation with post-processing."<<endl;
    PPRV.resize(node_num);
    Timer tt2;
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

    LabelV.postProcess(Label,threadnum);
    PPRV.postProcess(PruningPointSet,PruningPointSet2,vertices,threadnum);

    return t1+t2;
}
//function of original Dijkstra of PCL, 2023-06-16, pure search-based, vector-based implementation
void Graph::PCLDijkV(vertex nodeID, vector<vector<pair<vertex,int>>> &Neighbor){
    benchmark::heap<2, int, int> pqueue(node_num);
    pqueue.update(nodeID,0);

    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);
    vector<int> peak(node_num,-1);//record the peak hub from nodeID to topNodeID
    vector<int> pre(node_num,-1);//record the predecessor vertex in the shortest path
    vector<bool> pruned(node_num,false);//record the PPR pruning information

    pre[nodeID]=nodeID;
    peak[nodeID]=nodeID;
    distance[nodeID]=0;
    int topNodeID, topNodeDis;
    int NNodeID,NWeigh;
    bool ifPrune;

    bool ifDebug = false;
//    ifDebug = true;

    int lid, hid;
    lid=181451, hid=190161;
    lid=12516, hid=11449;

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        closed[topNodeID]=true;
        if(NodeOrder[peak[topNodeID]] < NodeOrder[peak[pre[topNodeID]]]) {
            peak[topNodeID] = peak[pre[topNodeID]];
        }
        if(pruned[pre[topNodeID]]){
            pruned[topNodeID]= true;
        }
        ifPrune = false;
        int TempDis; vector<int> SupNode;

        //Query of dispeak
        if(NodeOrder[nodeID] <= NodeOrder[topNodeID] && NodeOrder[peak[topNodeID]] <= NodeOrder[topNodeID]){//if valley path
//            vSm[nodeID]->wait();
//            Label[nodeID].insert(make_pair(topNodeID, topNodeDis));
            LabelV.add_lockfree(nodeID,topNodeID,topNodeDis);
//            vSm[nodeID]->notify();
        }

        //Label[topNodeID].insert(nodeID, topNodeDis);
        for(auto it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
            NNodeID=(*it).first;
            NWeigh=(*it).second+topNodeDis;

            if(!closed[NNodeID]){
                if(distance[NNodeID]>NWeigh){

                    pre[NNodeID]=topNodeID;
                    distance[NNodeID]=NWeigh;
                    pqueue.update(NNodeID, NWeigh);

                    if(peak[NNodeID] == -1){
                        peak[NNodeID] = NNodeID;
                    }

                }else if(distance[NNodeID]==NWeigh){//deal with equal distance scenario

                    if(NodeOrder[peak[topNodeID]] > NodeOrder[peak[pre[NNodeID]]]){//record the highest vertex among all shortest paths
//                        cout<<"change pre. L("<<nodeID<<","<<NNodeID<<"): "<<topNodeID<<"("<<NodeOrder[topNodeID]<<") "<<pre[NNodeID]<<"("<<NodeOrder[pre[NNodeID]]<<")"<<endl;
                        pre[NNodeID]=topNodeID;
                    }

                }
            }
//            if(NNodeID==lid){
//                cout<<topNodeID<< " "<<NNodeID<<" "<<pre[NNodeID]<<"("<<NodeOrder[pre[NNodeID]]<<")"<<endl;
//            }
        }
    }
}
void Graph::PCLDijkV2(vector<vertex>& ProcessID, vector<vector<pair<vertex,int>>> &Neighbor){
    for(int i=0;i<ProcessID.size();++i){
        int id=ProcessID[i];
        PCLDijkV(id,Neighbor);
    }
}
///**** BPCL ****/
//function for index construction of batch PLL
void Graph::BPCLIndexConstruct(vector<vector<pair<vertex,int>>> &Neighbor){
    bool ifParallel = true;
    Label.assign(node_num,unordered_map<vertex,int>());

    PruningPointSet2.clear();
    PruningPointSet2.assign(node_num,unordered_map<vertex,vertex>());
//    PruningPointSetOrder.clear();
//    PruningPointSetOrder.assign(node_num,unordered_map<int,set<OrderComp>>());
    PruningPointSet.clear();
    PruningPointSet.assign(node_num,unordered_map<vertex,unordered_set<vertex>>());

//    PruningPointVector.clear();
//    PruningPointVector.assign(node_num,vector<pair<int,int>>());


    int ID;
    int cnt=0;
    int stepShow = ceil(node_num/100000)*10000;
    stepShow = max(stepShow,1000);

    cout<<"Step for show: "<<stepShow<<endl;

    int batchSize = 512;
    batchSize = batchsize;
    if(batchSize>node_num){
        batchSize=node_num;
    }
//    batchSize = 3;
    cout<<"Batch size: "<<batchSize<<endl;

    vector<vector<vertex>> batches;
    vector<vertex> bNodes;//vertices for current batch
    vector<vector<vertex>> ProcessID;


    vector<vertex> vertices;

    Timer tt;
    double time=0;
    int a=0;

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
        hID = bNodes[0];
        unordered_set<vertex> setNodes; setNodes.clear();
        setNodes.insert(bNodes.begin(),bNodes.end());

        vector<pair<int,int>> vp;
        // process each batch
        tt.start();
        /// openmp-based implementation


        /// boost-based implementation
        if(ifParallel){//use multiple thread
//        if(false){
            if(batchSize > threadnum){
                ProcessID.assign(threadnum,vector<vertex>());
                for(int j=0;j<bNodes.size();++j){
                    a = j%threadnum;
                    ProcessID[a].emplace_back(bNodes[j]);
                }
                boost::thread_group thread;
                for(int i=0;i<ProcessID.size();i++){
                    thread.add_thread(new boost::thread(&Graph::BatchPCLDijk2, this, boost::ref(ProcessID[i]), boost::ref(setNodes), hID, boost::ref(Neighbor)));
                }
                thread.join_all();
            }else if(batchSize == threadnum){
                boost::thread_group thread;
                for(int i=0;i<bNodes.size();i++){
                    ID = bNodes[i];
                    thread.add_thread(new boost::thread(&Graph::BatchPCLDijk, this, ID, boost::ref(setNodes), hID, boost::ref(Neighbor)));
                }
                thread.join_all();
            }

        }else{//use single thread
            for(int i=0;i<bNodes.size();i++){
                ID = bNodes[i];
                BatchPCLDijk(ID, setNodes, hID, Neighbor);
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
    cout<<"Time used for label construction: "<<tt2.GetRuntime()<<" s."<<endl;

    vector<vector<vertex>> processID;
    processID.assign(threadnum, vector<vertex>());

    ThreadDistribute(vertices, processID);

    tt2.start();
    cout<<"Begin to construct pruning points..."<<endl;
    PruningPointBuild(ifParallel, processID, Neighbor);
//    PruningPointBuild(false);
    tt2.stop();
    cout<<"Time used for pruning point construction: "<<tt2.GetRuntime()<<" s."<<endl;
}
//function of pruning Dijkstra of batch PLL, 2023-06-07, optimal correct version 1
void Graph::BatchPCLDijk(vertex nodeID, unordered_set<vertex>& setNodes, vertex hID, vector<vector<pair<vertex,int>>> &Neighbor){
    benchmark::heap<2, int, int> pqueue(node_num);
    pqueue.update(nodeID,0);

    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);
    vector<int> peak(node_num,-1);//record the peak hub from nodeID to topNodeID
    vector<int> pre(node_num,-1);//record the predecessor vertex in the shortest path

    pre[nodeID]=nodeID;
    peak[nodeID]=nodeID;
    distance[nodeID]=0;
    int topNodeID, topNodeDis;
    int NNodeID,NWeigh;
    bool ifPrune;

    bool ifDebug = false;
//    ifDebug = true;

    int lid, hid;
    lid=181451, hid=190161;
    lid=200291, hid=200093;

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        closed[topNodeID]=true;
        if(NodeOrder[peak[topNodeID]] < NodeOrder[peak[pre[topNodeID]]]) {
            peak[topNodeID] = peak[pre[topNodeID]];
        }
        ifPrune = false;
        int TempDis; vector<int> SupNode;

        if(NodeOrder[topNodeID] > NodeOrder[hID]) {
            continue;
        }
        else{//if NodeOrder[topNodeID] <= NodeOrder[hID]
            //Query of dispeak
            if(NodeOrder[nodeID] >= NodeOrder[topNodeID]){
                vSm[nodeID]->wait();
                if(topNodeID != nodeID){
                    vSm[topNodeID]->wait();
                }
                ShortestDisQuery1(topNodeID, nodeID, SupNode,TempDis);
                vSm[nodeID]->notify();
                if(topNodeID != nodeID){
                    vSm[topNodeID]->notify();
                }

                if(ifDebug){
                    if((nodeID==lid && topNodeID==hid) || (nodeID==hid && topNodeID==lid)){
                        pair<int,int> tempP = DijkstraCoreDebug(nodeID,topNodeID);
                        cout<<"Find 1. "<<nodeID<<"("<<NodeOrder[nodeID]<<") "<<topNodeID<<"("<<NodeOrder[topNodeID]<<"): "<<topNodeDis<<" "<<TempDis<<", peak in vally path: "<<peak[topNodeID]<<"("<<NodeOrder[peak[topNodeID]]<<"). true length: "<<DijkstraCore(nodeID,peak[topNodeID])+DijkstraCore(peak[topNodeID],topNodeID)<<" "<<tempP.first<<"("<<tempP.second<<","<<NodeOrder[tempP.second]<<")"<<endl;
                    }
                }

                if(TempDis>topNodeDis) {//if the peak path is longer than this vally path
                    if((NodeOrder[peak[topNodeID]] <= NodeOrder[topNodeID]) || (NodeOrder[peak[topNodeID]] <= NodeOrder[nodeID])) {

                        vSm[topNodeID]->wait();
                        if (Label[topNodeID].find(nodeID) == Label[topNodeID].end()) {//if not found
                            Label[topNodeID].insert(make_pair(nodeID, topNodeDis));
                        }
                        vSm[topNodeID]->notify();
                    }
//                    else {//if the peakhub is higher than topNodeID
//                        cout<<"Skip L1 "<<nodeID<<"("<<NodeOrder[nodeID]<<"),"<<topNodeID<<"("<<NodeOrder[topNodeID]<<"): "<<peakhubs[topNodeID]<<"("<<NodeOrder[peakhubs[topNodeID]]<<")"<<endl;
//                    }
                }
                else{//if TempDis<=topNodeDis
                    int highestHub = -1;
                    for (int k = 0; k < SupNode.size(); k++) {
                        int supn = SupNode[k];
                        if(highestHub == -1){
                            highestHub = supn;
                        }else if(NodeOrder[supn]>NodeOrder[highestHub]){
                            highestHub = supn;
                        }

                        ifPrune = true;
                        //cout<<"pruning point "<<topNodeID<<" "<<nodeID<<" "<<supn<<endl;
                    }

                    if(NodeOrder[highestHub] > NodeOrder[hID]){
                        continue;
                    }
//                    else{
//                        cout<<"Not continue 1. "<<nodeID<<" "<<topNodeID<<": "<<highestHub<<"("<<NodeOrder[highestHub]<<") "<<hID<<"("<<NodeOrder[hID]<<")"<<endl;
//                    }

                }
            }
            else{//if NodeOrder[nodeID] < NodeOrder[topNodeID]
                vSm[topNodeID]->wait();
                if(topNodeID != nodeID){
                    vSm[nodeID]->wait();
                }
                ShortestDisQuery1(nodeID, topNodeID,SupNode,TempDis);
                vSm[topNodeID]->notify();
                if(topNodeID != nodeID){
                    vSm[nodeID]->notify();
                }
                if(ifDebug){
                    if((nodeID==lid && topNodeID==hid) || (nodeID==hid && topNodeID==lid)){
                        pair<int,int> tempP = DijkstraCoreDebug(nodeID,topNodeID);
                        cout<<"Find 2. "<<nodeID<<"("<<NodeOrder[nodeID]<<") "<<topNodeID<<"("<<NodeOrder[topNodeID]<<"): "<<topNodeDis<<" "<<TempDis<<", peak in vally path: "<<peak[topNodeID]<<"("<<NodeOrder[peak[topNodeID]]<<"). true length: "<<DijkstraCore(nodeID,peak[topNodeID])+DijkstraCore(peak[topNodeID],topNodeID)<<" "<<tempP.first<<"("<<tempP.second<<","<<NodeOrder[tempP.second]<<")"<<endl;
                    }
                }

                if(TempDis>topNodeDis) {//if the peak path is longer than this vally path
                    if((NodeOrder[peak[topNodeID]] <= NodeOrder[topNodeID]) || (NodeOrder[peak[topNodeID]] <= NodeOrder[nodeID])) {//if the topNodeID is the highest vertex in the path, reversely insert it.

                        bool update=false;
                        vSm[nodeID]->wait();
                        if (Label[nodeID].find(topNodeID) == Label[nodeID].end()) {//if not found
                            Label[nodeID].insert(make_pair(topNodeID, topNodeDis));
                            update=true;
                        }else if(Label[nodeID][topNodeID] > topNodeDis){
                            cout<<"Larger label!! "<<Label[nodeID][topNodeID]<<" "<<topNodeDis<<endl;
                            Label[nodeID][topNodeID] = topNodeDis;
                            update=true;
                        }
                        vSm[nodeID]->notify();

//                        if(ifDebug){
//                            if(update){
//                                pair<int,int> DijkPair= DijkstraCoreDebug(nodeID,topNodeID);
//                                if(topNodeDis != DijkPair.first){
//                                    cout<<"reverse insertion. L "<<nodeID<<"("<<NodeOrder[nodeID]<<"),"<<topNodeID<<"("<<NodeOrder[topNodeID]<<"): "<<Label[nodeID][topNodeID]<<" "<<topNodeDis<<" "<<TempDis<<" "<<DijkPair.first<<"("<<NodeOrder[DijkPair.second]<<")"<<endl;
//                                }
//                            }
//                        }
                    }
//                    else{//if the peakhub is higher than topNodeID
//                        cout<<"Skip L2 "<<nodeID<<"("<<NodeOrder[nodeID]<<"),"<<topNodeID<<"("<<NodeOrder[topNodeID]<<"): "<<peakhubs[topNodeID]<<"("<<NodeOrder[peakhubs[topNodeID]]<<")"<<endl;
//                    }
                }
                else{//if TempDis<=topNodeDis
                    int highestHub = -1;
                    for (int k = 0; k < SupNode.size(); k++) {
                        int supn = SupNode[k];
                        if(highestHub == -1){
                            highestHub = supn;
                        }else if(NodeOrder[supn]>NodeOrder[highestHub]){//>
                            highestHub = supn;
                        }

                        ifPrune = true;
                        //cout<<"pruning point "<<topNodeID<<" "<<nodeID<<" "<<supn<<endl;
                    }

                    if(NodeOrder[highestHub] > NodeOrder[hID]){//if the dispeak is sourced from higher vertex
                        continue;
                    }
//                    else{//else we have to propagate the information to neighbors to ensure the vally paths under hID can be recovered.
//                        cout<<"Not continue 2. "<<nodeID<<" "<<topNodeID<<": "<<highestHub<<"("<<NodeOrder[highestHub]<<") "<<hID<<"("<<NodeOrder[hID]<<")"<<endl;
//                    }

                }
            }
        }

        //Label[topNodeID].insert(nodeID, topNodeDis);
        for(auto it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
            NNodeID=(*it).first;
            NWeigh=(*it).second+topNodeDis;
            if(!closed[NNodeID]){
                if(distance[NNodeID]>NWeigh){
                    if(ifDebug){
                        if(nodeID == hid && NNodeID == lid){
                            cout<<"here! "<<nodeID<<" "<<topNodeID<<"("<<peak[topNodeID]<<","<<NodeOrder[peak[topNodeID]]<<"), "<<NNodeID<<"("<<peak[NNodeID]<<","<<NodeOrder[peak[NNodeID]]<<"): "<<distance[NNodeID]<<" "<<NWeigh<<" "<<endl;
                        }
                    }

                    pre[NNodeID]=topNodeID;
                    distance[NNodeID]=NWeigh;
                    pqueue.update(NNodeID, NWeigh);

                    if(peak[NNodeID] == -1){
                        peak[NNodeID] = NNodeID;
                    }

                }else if(distance[NNodeID]==NWeigh){//deal with equal distance scenario
                    if(NodeOrder[peak[topNodeID]] > NodeOrder[peak[pre[NNodeID]]]){//record the highest vertex among all shortest paths
//                        cout<<"change pre. L("<<nodeID<<","<<NNodeID<<"): "<<topNodeID<<"("<<NodeOrder[topNodeID]<<") "<<pre[NNodeID]<<"("<<NodeOrder[pre[NNodeID]]<<")"<<endl;
                        pre[NNodeID]=topNodeID;
                    }

                }
            }
        }
    }
}
//function of pruning Dijkstra of batch PLL, 2023-06-07, optimal correct version 1
void Graph::BatchPCLDijk2(vector<vertex>& p, unordered_set<vertex>& setNodes, vertex hID, vector<vector<pair<vertex,int>>> &Neighbor){
    for(auto it=p.begin();it!=p.end();++it){
        int nodeID = *it;
        BatchPCLDijk(nodeID, setNodes, hID, Neighbor);
    }
}

// vector-based implementation
double Graph::BPCLConstructV(vector<vector<pair<vertex,int>>> &Neighbor){
    double t1,t2;

    cout<<"Vector-based implementation. Boost-based labeling."<<endl;
    bool ifParallel = true;

    Label.assign(node_num,unordered_map<vertex,int>());
    LabelV.resize(node_num);

    int ID;
    int cnt=0;
//    int stepShow = ceil(node_num/100000)*10000;
//    stepShow = max(stepShow,1000);
//    cout<<"Step for show: "<<stepShow<<endl;

    int batchSize = 512;
    batchSize = batchsize;
    if(batchSize>node_num){
        batchSize=node_num;
    }
//    batchSize = 3;
    cout<<"Batch size: "<<batchSize<<endl;

    vector<vector<vertex>> batches;
    vector<vertex> bNodes;//vertices for current batch
    vector<vector<vertex>> ProcessID;


    vector<vertex> vertices;

    Timer tt;
    double time=0;
    int a=0;

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
        hID = bNodes[0];
        unordered_set<vertex> setNodes; setNodes.clear();
        setNodes.insert(bNodes.begin(),bNodes.end());

        vector<pair<int,int>> vp;
        // process each batch
        tt.start();

        if(ifParallel){//use multiple thread
//        if(false){
            /// boost-based implementation
            if(batchSize > threadnum){
                ProcessID.assign(threadnum,vector<vertex>());
                for(int j=0;j<bNodes.size();++j){
                    a = j%threadnum;
                    ProcessID[a].emplace_back(bNodes[j]);
                }
                boost::thread_group thread;
                for(int i=0;i<ProcessID.size();i++){
                    thread.add_thread(new boost::thread(&Graph::BatchPCLDijk2V, this, boost::ref(ProcessID[i]), boost::ref(setNodes), hID, boost::ref(Neighbor)));
                }
                thread.join_all();


            }else if(batchSize == threadnum){
                boost::thread_group thread;
                for(int i=0;i<bNodes.size();i++){
                    ID = bNodes[i];
                    thread.add_thread(new boost::thread(&Graph::BatchPCLDijkV, this, ID, boost::ref(setNodes), hID, boost::ref(Neighbor)));
                }
                thread.join_all();
            }

//            /// openmp-based implementation
//                #pragma omp parallel for num_threads(threadnum) schedule (dynamic, threadnum)
//                for(int j=0;j<bNodes.size();++j){
//                    int id=bNodes[j];
//                    BatchPCLDijkV(id, setNodes, hID, Neighbor);
//                }

        }else{//use single thread
            for(int i=0;i<bNodes.size();i++){
                ID = bNodes[i];
                BatchPCLDijkV(ID, setNodes, hID, Neighbor);
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
    t1=tt2.GetRuntime();
    LabelV.sort(threadnum);
    cout<<"Time used for label construction: "<<t1<<" s."<<endl;

    if(algoCoreU==0){
        /// boost-based based implementation, post-processing
        cout<<"Begin to construct pruning points... boost-based implementation with post-processing."<<endl;
        PPRV.resize(node_num);
        tt2.start();
        ThreadDistribute(vertices,ProcessID);
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
//    PPRV.postProcess(PruningPointSet,PruningPointSet2,vertices);
    }


//    WriteLabels(graphfile);
    return t1+t2;
}
// batch pruned Dijkstra
void Graph::BatchPCLDijkV(vertex nodeID, unordered_set<vertex>& setNodes, vertex hID, vector<vector<pair<vertex,int>>> &Neighbor){
    benchmark::heap<2, int, int> pqueue(node_num);
    pqueue.update(nodeID,0);

//    unordered_map<vertex,uint> Lh;
//    for(auto it=LabelV.Labels[nodeID].begin();it!=LabelV.Labels[nodeID].end();++it){
//        Lh.insert({it->first,it->second});
//    }

    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);
    vector<int> peak(node_num,-1);//record the peak hub from nodeID to topNodeID
    vector<int> pre(node_num,-1);//record the predecessor vertex in the shortest path

    pre[nodeID]=nodeID;
    peak[nodeID]=nodeID;
    distance[nodeID]=0;
    int topNodeID, topNodeDis;
    int NNodeID,NWeigh;
    bool ifPrune;

    bool ifDebug = false;
//    ifDebug = true;

    int lid, hid;
    lid=181451, hid=190161;
    lid=200291, hid=200093;

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        closed[topNodeID]=true;
        if(NodeOrder[peak[topNodeID]] < NodeOrder[peak[pre[topNodeID]]]) {
            peak[topNodeID] = peak[pre[topNodeID]];
        }
        ifPrune = false;
        int TempDis; vector<int> SupNode;

        if(NodeOrder[topNodeID] > NodeOrder[hID]) {
            continue;
        }
        else{//if NodeOrder[topNodeID] <= NodeOrder[hID]
            //Query of dispeak
            if(NodeOrder[nodeID] >= NodeOrder[topNodeID]){
                vSm[nodeID]->wait();
                if(topNodeID != nodeID){
                    vSm[topNodeID]->wait();
                }
                ShortestDisQuery1V(topNodeID, Label[nodeID], SupNode,TempDis);
                vSm[nodeID]->notify();
                if(topNodeID != nodeID){
                    vSm[topNodeID]->notify();
                }


                if(TempDis>topNodeDis) {//if the peak path is longer than this vally path
                    if((NodeOrder[peak[topNodeID]] <= NodeOrder[topNodeID]) || (NodeOrder[peak[topNodeID]] <= NodeOrder[nodeID])) {

                        vSm[topNodeID]->wait();
                        if (Label[topNodeID].find(nodeID) == Label[topNodeID].end()) {//if not found
                            Label[topNodeID].insert(make_pair(nodeID, topNodeDis));
                            LabelV.add(topNodeID,nodeID,topNodeDis);
                        }

                        vSm[topNodeID]->notify();
                    }

                }
                else{//if TempDis<=topNodeDis
                    int highestHub = -1;
                    for (int k = 0; k < SupNode.size(); k++) {
                        int supn = SupNode[k];
                        if(highestHub == -1){
                            highestHub = supn;
                        }else if(NodeOrder[supn]>NodeOrder[highestHub]){
                            highestHub = supn;
                        }

                        ifPrune = true;
                        //cout<<"pruning point "<<topNodeID<<" "<<nodeID<<" "<<supn<<endl;
                    }

                    if(NodeOrder[highestHub] > NodeOrder[hID]){
                        continue;
                    }


                }
            }
            else{//if NodeOrder[nodeID] < NodeOrder[topNodeID]
                vSm[topNodeID]->wait();
                if(topNodeID != nodeID){
                    vSm[nodeID]->wait();
                }
                ShortestDisQuery1V(topNodeID, Label[nodeID],SupNode,TempDis);
                vSm[topNodeID]->notify();
                if(topNodeID != nodeID){
                    vSm[nodeID]->notify();
                }

                if(TempDis>topNodeDis) {//if the peak path is longer than this vally path
                    if((NodeOrder[peak[topNodeID]] <= NodeOrder[topNodeID]) || (NodeOrder[peak[topNodeID]] <= NodeOrder[nodeID])) {//if the topNodeID is the highest vertex in the path, reversely insert it.

                        bool update=false;
                        vSm[nodeID]->wait();
                        if (Label[nodeID].find(topNodeID) == Label[nodeID].end()) {//if not found
                            Label[nodeID].insert(make_pair(topNodeID, topNodeDis));
                            LabelV.add(nodeID,topNodeID,topNodeDis);
                            update=true;
                        }else if(Label[nodeID][topNodeID] > topNodeDis){
                            cout<<"Larger label!! "<<Label[nodeID][topNodeID]<<" "<<topNodeDis<<endl;
                            Label[nodeID][topNodeID] = topNodeDis;
                            update=true;
                            exit(1);
                        }
                        vSm[nodeID]->notify();

                    }
                }
                else{//if TempDis<=topNodeDis
                    int highestHub = -1;
                    for (int k = 0; k < SupNode.size(); k++) {
                        int supn = SupNode[k];
                        if(highestHub == -1){
                            highestHub = supn;
                        }else if(NodeOrder[supn]>NodeOrder[highestHub]){//>
                            highestHub = supn;
                        }

                        ifPrune = true;
                        //cout<<"pruning point "<<topNodeID<<" "<<nodeID<<" "<<supn<<endl;
                    }

                    if(NodeOrder[highestHub] > NodeOrder[hID]){//if the dispeak is sourced from higher vertex
                        continue;
                    }

                }
            }
        }

        //Label[topNodeID].insert(nodeID, topNodeDis);
        for(auto it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
            NNodeID=(*it).first;
            NWeigh=(*it).second+topNodeDis;
            if(!closed[NNodeID]){
                if(distance[NNodeID]>NWeigh){
                    if(ifDebug){
                        if(nodeID == hid && NNodeID == lid){
                            cout<<"here! "<<nodeID<<" "<<topNodeID<<"("<<peak[topNodeID]<<","<<NodeOrder[peak[topNodeID]]<<"), "<<NNodeID<<"("<<peak[NNodeID]<<","<<NodeOrder[peak[NNodeID]]<<"): "<<distance[NNodeID]<<" "<<NWeigh<<" "<<endl;
                        }
                    }

                    pre[NNodeID]=topNodeID;
                    distance[NNodeID]=NWeigh;
                    pqueue.update(NNodeID, NWeigh);

                    if(peak[NNodeID] == -1){
                        peak[NNodeID] = NNodeID;
                    }

                }else if(distance[NNodeID]==NWeigh){//deal with equal distance scenario
                    if(NodeOrder[peak[topNodeID]] > NodeOrder[peak[pre[NNodeID]]]){//record the highest vertex among all shortest paths
//                        cout<<"change pre. L("<<nodeID<<","<<NNodeID<<"): "<<topNodeID<<"("<<NodeOrder[topNodeID]<<") "<<pre[NNodeID]<<"("<<NodeOrder[pre[NNodeID]]<<")"<<endl;
                        pre[NNodeID]=topNodeID;
                    }

                }
            }
        }
    }
}
void Graph::BatchPCLDijk2V(vector<vertex>& ProcessID,  unordered_set<vertex>& setNodes, vertex hID, vector<vector<pair<vertex,int>>> &Neighbor){
    for(int i=0;i<ProcessID.size();++i){
        BatchPCLDijkV(ProcessID[i],setNodes,hID,Neighbor);
    }
}

void Graph::PCLDijkPPR(int nodeID){
    benchmark::heap<2, int, int> pqueue(node_num);
    pqueue.update(nodeID,0);

    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);
    vector<int> peak(node_num,-1);//record the peak hub from nodeID to topNodeID
    vector<int> pre(node_num,-1);//record the predecessor vertex in the shortest path
    vector<bool> pruned(node_num,false);//record the PPR pruning information
    vector<vector<int>> pres(node_num,vector<int>());

    pres[nodeID].emplace_back(nodeID);
    pre[nodeID]=nodeID;
    peak[nodeID]=nodeID;
    distance[nodeID]=0;
    int topNodeID, topNodeDis;
    int NNodeID,NWeigh;
    bool ifPrune;

    bool ifDebug = false;
//    ifDebug = true;

    int lid, hid;
    lid=181451, hid=190161;
    lid=12516, hid=11449;

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        closed[topNodeID]=true;
        if(NodeOrder[peak[topNodeID]] < NodeOrder[peak[pre[topNodeID]]]) {
            peak[topNodeID] = peak[pre[topNodeID]];
        }
        if(pruned[pre[topNodeID]]){
            pruned[topNodeID]= true;
        }
        ifPrune = false;
        int TempDis; vector<int> SupNode;

        //Query of dispeak
        if(NodeOrder[nodeID] >= NodeOrder[topNodeID]){
            if(NodeOrder[peak[topNodeID]] <= NodeOrder[nodeID]) {//if valley path
                vSm[topNodeID]->wait();
                if (Label[topNodeID].find(nodeID) == Label[topNodeID].end()) {//if not found
                    Label[topNodeID].insert(make_pair(nodeID, topNodeDis));
                }
                vSm[topNodeID]->notify();
            }else{//if peak path
                ifPrune = false;
                if(!pruned[topNodeID]) {//if the predecessor of topNodeID has not been pruned
//                        cout<<"here!!"<<endl;
                    for (auto it = pres[topNodeID].begin(); it != pres[topNodeID].end(); ++it) {
                        int supn = peak[*it];
                        if ((NodeOrder[topNodeID] < NodeOrder[supn]) &&
                            (NodeOrder[nodeID] < NodeOrder[supn])) {
                            vSm[nodeID]->wait();
//                            PruningPointSetOrder[nodeID][supn].insert(topNodeID);
                            PruningPointSet[nodeID][supn].insert(topNodeID);
                            PruningPointSet2[nodeID][topNodeID] = supn;
                            vSm[nodeID]->notify();

                            vSm[topNodeID]->wait();
//                            PruningPointSetOrder[topNodeID][supn].insert(nodeID);
                            PruningPointSet[topNodeID][supn].insert(nodeID);
                            PruningPointSet2[topNodeID][nodeID] = supn;
                            vSm[topNodeID]->notify();
                            ifPrune = true;
                        }

                    }

                }


                if(ifPrune){
                    pruned[topNodeID] = true;
                }
            }

        }
        else{//if NodeOrder[nodeID] < NodeOrder[topNodeID]
            if(NodeOrder[peak[topNodeID]] <= NodeOrder[topNodeID]) {//if the topNodeID is the highest vertex in the path, reversely insert it.

                bool update=false;
                vSm[nodeID]->wait();
                if (Label[nodeID].find(topNodeID) == Label[nodeID].end()) {//if not found
                    Label[nodeID].insert(make_pair(topNodeID, topNodeDis));
                    update=true;
                }
                vSm[nodeID]->notify();
            }
            else{
                ifPrune = false;
                if(!pruned[topNodeID]) {//if the predecessor of topNodeID has not been pruned

                    for (auto it = pres[topNodeID].begin(); it != pres[topNodeID].end(); ++it) {
                        int supn = peak[*it];
                        if ((NodeOrder[topNodeID] < NodeOrder[supn]) &&
                            (NodeOrder[nodeID] < NodeOrder[supn])) {
                            vSm[nodeID]->wait();
//                            PruningPointSetOrder[nodeID][supn].insert(topNodeID);
                            PruningPointSet[nodeID][supn].insert(topNodeID);
                            PruningPointSet2[nodeID][topNodeID] = supn;
                            vSm[nodeID]->notify();

                            vSm[topNodeID]->wait();
//                            PruningPointSetOrder[topNodeID][supn].insert(nodeID);
                            PruningPointSet[topNodeID][supn].insert(nodeID);
                            PruningPointSet2[topNodeID][nodeID] = supn;
                            vSm[topNodeID]->notify();
                            ifPrune = true;
                        }

                    }

                }

                if(ifPrune){
                    pruned[topNodeID] = true;
                }
            }

        }


        //Label[topNodeID].insert(nodeID, topNodeDis);
        for(auto it=AdjaCore[topNodeID].begin();it!=AdjaCore[topNodeID].end();it++){
            NNodeID=(*it).first;
            NWeigh=(*it).second+topNodeDis;

            if(!closed[NNodeID]){
                if(distance[NNodeID]>NWeigh){

                    pre[NNodeID]=topNodeID;
                    distance[NNodeID]=NWeigh;
                    pqueue.update(NNodeID, NWeigh);

                    pres[NNodeID].clear();
                    pres[NNodeID].emplace_back(topNodeID);

                    if(peak[NNodeID] == -1){
                        peak[NNodeID] = NNodeID;
                    }

                }else if(distance[NNodeID]==NWeigh){//deal with equal distance scenario

                    if(NodeOrder[peak[topNodeID]] > NodeOrder[peak[pre[NNodeID]]]){//record the highest vertex among all shortest paths
//                        cout<<"change pre. L("<<nodeID<<","<<NNodeID<<"): "<<topNodeID<<"("<<NodeOrder[topNodeID]<<") "<<pre[NNodeID]<<"("<<NodeOrder[pre[NNodeID]]<<")"<<endl;
                        pre[NNodeID]=topNodeID;
                    }

                    pres[NNodeID].emplace_back(topNodeID);

                }
            }
//            if(NNodeID==lid){
//                cout<<topNodeID<< " "<<NNodeID<<" "<<pre[NNodeID]<<"("<<NodeOrder[pre[NNodeID]]<<")"<<endl;
//            }
        }
    }
}
//function of vertex allocation
void Graph::ThreadDistribute(vector<vertex>& vertices, vector<vector<vertex>>& processID){
    processID.assign(threadnum, vector<vertex>());
    int pid=0;
    for(int i=0;i<vertices.size();++i){
        pid=i%threadnum;
        processID[pid].emplace_back(vertices[i]);
    }
}
//function of building the pruning point records for a given vertex
void Graph::PPRConstruction(vertex nodeID, vector<vector<pair<vertex,int>>> &Neighbor){
    benchmark::heap<2, int, int> pqueue(node_num);
    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);

    pqueue.update(nodeID,0);
    distance[nodeID]=0;
    int topNodeID, topNodeDis;
    int NNodeID,NWeigh;

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        closed[topNodeID]=true;

        if(NodeOrder[topNodeID] > NodeOrder[nodeID]){
            continue;
        }

        int TempDis; vector<int> SupNode;
        //if NodeOrder[topNodeID] <= NodeOrder[nodeID]
        ShortestDisQuery2(topNodeID, nodeID, SupNode,TempDis);//DisQueryPeak
//        ShortestDisQuery1(topNodeID, nodeID, SupNode,TempDis);

        if(TempDis<=topNodeDis){//if dispeak <= disvally
            for(int k=0;k<SupNode.size();k++){
                int supn=SupNode[k];
                unordered_map<int,vector<int>> map1;

                vSm[nodeID]->wait();
//                PruningPointVector[nodeID].emplace_back(supn,topNodeID);
//                PruningPointSetOrder[nodeID][supn].insert(topNodeID);
                PruningPointSet[nodeID][supn].insert(topNodeID);
                PruningPointSet2[nodeID][topNodeID] = supn;
                vSm[nodeID]->notify();

                vSm[topNodeID]->wait();
//                PruningPointVector[topNodeID].emplace_back(supn,nodeID);
//                PruningPointSetOrder[topNodeID][supn].insert(nodeID);
                PruningPointSet[topNodeID][supn].insert(nodeID);
                PruningPointSet2[topNodeID][nodeID] = supn;
                vSm[topNodeID]->notify();

                //cout<<"pruning point "<<topNodeID<<" "<<nodeID<<" "<<supn<<endl;
            }
            continue;

        }

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

void Graph::PPRConstruction2(vector<vertex> & p, vector<vector<pair<vertex,int>>> &Neighbor){
//    cout<<"1"<<endl;
    for(int i=0;i<p.size();++i){
        PPRConstruction(p[i], Neighbor);
    }
}
//function of building the pruning points
void Graph::PruningPointBuild(bool ifParallel, vector<vector<vertex>> & processID, vector<vector<pair<vertex,int>>> &Neighbor){

    if(ifParallel){
//    if(false){
        cout<<"Batch size: "<<processID[0].size()<<endl;
        boost::thread_group thread;
        for(auto j=0;j<processID.size();++j){
            thread.add_thread(new boost::thread(&Graph::PPRConstruction2, this, boost::ref(processID[j]), boost::ref(Neighbor) ));
        }
        thread.join_all();
    }else{
        for(auto j=0;j<processID.size();++j){
            PPRConstruction2(processID[j], Neighbor);
        }
    }
}
//vector-based implementation
void Graph::PPRConstructionV(vertex nodeID, vector<vector<pair<vertex,int>>> &Neighbor){
    benchmark::heap<2, int, int> pqueue(node_num);
    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);

    pqueue.update(nodeID,0);
    distance[nodeID]=0;
    int topNodeID, topNodeDis;
    int NNodeID,NWeigh;

    unordered_map<vertex,int> Lh;
    for(int i=0;i<LabelV.Labels[nodeID].size();++i){
        Lh.insert({LabelV.Labels[nodeID][i].first,LabelV.Labels[nodeID][i].second});
    }

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        closed[topNodeID]=true;

        if(NodeOrder[topNodeID] > NodeOrder[nodeID]){
            continue;
        }

        int TempDis; vector<int> SupNode;
        //if NodeOrder[topNodeID] <= NodeOrder[nodeID]
        ShortestDisQueryPeakV(topNodeID, nodeID, Lh, SupNode,TempDis);//DisQueryPeak
//        ShortestDisQueryPeakV(topNodeID, nodeID, Label[nodeID], SupNode,TempDis);//DisQueryPeak
//        ShortestDisQueryPeak(topNodeID, nodeID, SupNode,TempDis);//DisQueryPeak

        if(TempDis<=topNodeDis){//if dispeak <= disvally
            for(int k=0;k<SupNode.size();k++){
                int supn=SupNode[k];
                unordered_map<int,vector<int>> map1;

                PPRV.add_lockfree(nodeID,supn,topNodeID);

                //cout<<"pruning point "<<topNodeID<<" "<<nodeID<<" "<<supn<<endl;
            }
            continue;

        }

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
void Graph::PPRConstructionV2(vector<vertex>& ProcessID, vector<vector<pair<vertex,int>>> &Neighbor){
    for(int i=0;i<ProcessID.size();++i){
        int id=ProcessID[i];
        PPRConstructionV(id,Neighbor);
    }
}
int Graph::ShortestDisQueryPeakV(int ID1,int ID2, unordered_map<vertex,int>& Lh, vector<int>& SupNode, int& d){
    d=INF;
    int hub, dis1, dis2;
    int finalHub=-1;

    for(auto it=LabelV.Labels[ID1].begin();it!=LabelV.Labels[ID1].end();++it){
        hub=(*it).first;
        dis1=(*it).second;
        if(hub!=ID2 && Lh.find(hub)!=Lh.end()){
            dis2=Lh[hub];
            if(dis1+dis2<d){
                finalHub=hub;
                d=dis1+dis2;
                SupNode.clear();
                SupNode.push_back(hub);
            }else if(dis1+dis2==d){
                SupNode.push_back(hub);
            }
        }
    }

    return finalHub;
}

//// Index Maintenance
//function for maintaining the edge weight decrease updates of 2-hop labeling
void Graph::DecreasePSL(int a, int b, int oldW, int newW, vector<vector<pair<vertex,int>>> &Neighbors,vector<unordered_map<vertex,int>> &Label){
    //update the edges on original graph
    for(int i=0;i<Neighbors[a].size();i++){
        if(Neighbors[a][i].first==b){
            Neighbors[a][i].second=newW;
            break;
        }
    }
    for(int i=0;i<Neighbors[b].size();i++){
        if(Neighbors[b][i].first==a){
            Neighbors[b][i].second=newW;
            break;
        }
    }

    //check the dis(a,b)
    int Dab=PLLQuery(a,b,Label);//query by PSL label

    int LID,HID;
    if(NodeOrder[a]>NodeOrder[b]){
        LID=b; HID=a;
    }else{
        LID=a; HID=b;
    }

    if(Dab>newW){//the index change is triggered
        vector<vector<pair<int,int>>> Change;
        vector<pair<int,int>> vec;
        Change.assign(node_num,vec);
        set<int> WaitPro;

        Label[LID][HID]=newW;
//        if(LabelV.Labels[LID].size()!=Label[LID].size()){
//            LabelV.Labels[LID].emplace_back(HID,newW);
//        }
        Change[LID].push_back(make_pair(HID,newW));
        WaitPro.insert(LID);

        //check the label of a,b
        int hubid, hubdis;

        for(auto it1=Label[LID].begin();it1!=Label[LID].end();it1++){
//        for(auto it1=LabelV.Labels[LID].begin();it1!=LabelV.Labels[LID].end();++it1){
            hubid=(*it1).first;
            hubdis=(*it1).second;
//            if(Label[LID].find(hubid)==Label[LID].end()){
//                cout<<"Wrong! Not found!"<<endl; exit(1);
//            }
            hubdis=Label[LID][hubid];
            if(NodeOrder[hubid]>NodeOrder[HID] && newW+hubdis<PLLQuery(HID,hubid,Label)){
                Label[HID][hubid]=newW+hubdis;
//                if(LabelV.Labels[HID].size()!=Label[HID].size()){
//                    LabelV.Labels[HID].emplace_back(hubid,newW+hubdis);
//                }
                Change[HID].push_back(make_pair(hubid, newW+hubdis));
                WaitPro.insert(HID);
            }
        }

        for(auto it2=Label[HID].begin();it2!=Label[HID].end();it2++){
//        for(auto it2=LabelV.Labels[HID].begin();it2!=LabelV.Labels[HID].end();++it2){
            hubid=(*it2).first;
            hubdis=(*it2).second;
//            if(Label[HID].find(hubid)==Label[HID].end()){
//                cout<<"Wrong! Not found!"<<endl; exit(1);
//            }
            hubdis=Label[HID][hubid];
            if(newW+hubdis<PLLQuery(LID, hubid,Label)){
                Label[LID][hubid]=newW+hubdis;
//                if(LabelV.Labels[LID].size()!=Label[LID].size()){
//                    LabelV.Labels[LID].emplace_back(hubid,newW+hubdis);
//                }
                Change[LID].push_back(make_pair(hubid, newW+hubdis));
                WaitPro.insert(LID);
            }
        }

        //check the label of their neighbors step by step
        while(WaitPro.size()>0){
            set<int> WaitProTem;
            vector<vector<pair<int,int>>> ChangeTem;
            vector<pair<int,int>> vec;
            ChangeTem.assign(node_num,vec);

            for(set<int>::iterator it=WaitPro.begin();it!=WaitPro.end();it++){
                int curID=*it;
                vector<pair<int,int>> curChange=Change[curID];
                int neiID, neiDis, hID, hDis;
                for(int i=0;i<Neighbors[curID].size();i++){
                    neiID=Neighbors[curID][i].first;
                    neiDis=Neighbors[curID][i].second;

                    for(int j=0;j<curChange.size();j++){//for each decreased label, we check all its neighbors
                        hID=curChange[j].first; hDis=curChange[j].second;
                        if(NodeOrder[hID]>NodeOrder[neiID] && PLLQuery(neiID, hID,Label)>neiDis+hDis){
                            Label[neiID][hID]=neiDis+hDis;
//                            if(LabelV.Labels[neiID].size()!=Label[neiID].size()){
//                                LabelV.Labels[neiID].emplace_back(hID,neiDis+hDis);
//                            }
                            WaitProTem.insert(neiID);
                            ChangeTem[neiID].push_back(make_pair(hID, neiDis+hDis));
                        }
                    }
                }
            }

            WaitPro=WaitProTem;
            Change=ChangeTem;
        }
    }
}
//new version: set version with NoSupportedPair, clean label, 2023-05-14, iterative update, queue version, correct
void Graph::IncreasePSLNew(int a, int b, int oldW, int newW, vector<vector<pair<vertex,int>>> &Neighbor,vector<unordered_map<vertex,int>> &Label, vector<unordered_map<vertex,unordered_set<vertex>>> &PruningPointNew){//
    for(int i=0;i<Neighbor[a].size();i++){
        if(Neighbor[a][i].first==b){
//            cout<<"Core "<<a<<" "<<b<<" "<<Neighbor[a][i].second<<" "<<newW<<endl;
            Neighbor[a][i].second=newW;
            break;
        }
    }
    for(int i=0;i<Neighbor[b].size();i++){
        if(Neighbor[b][i].first==a){
//            cout<<"Core "<<b<<" "<<a<<" "<<Neighbor[b][i].second<<" "<<newW<<endl;
            Neighbor[b][i].second=newW;
            break;
        }
    }

    int LID, HID;
    if(NodeOrder[a]>NodeOrder[b]){
        LID=b; HID=a;
    }else{
        LID=a; HID=b;
    }

    ///for debug
    bool ifDebug = false;
//    ifDebug = true;

    int lid, hid;
    lid=735, hid=294;//core 29-1

    int lid2,hid2;
    lid2=147082, hid2=196720;

//    cout<<"LID: "<<LID<<" ("<<NodeOrder[LID]<<") ; HID: "<<HID<<" ("<<NodeOrder[HID]<<")"<<endl;
    int dislower,disvally,dispeak,peakhub,vallyID;
    pair<int,int> peakPair;//distance, hubID
    pair<int,int> vallyPair;//distance, vallyID
    set<tuple<int,int,int>> outdatedPruning;//<nodeID,supportNode,prunedID>
//    map<pair<int,int>,int> outdatedPruning;//<<nodeID,prunedID>,supportNode>
    map<pair<int,int>,int> newPruningPoints;//<<nodeID,prunedID>,supportNode>
    outdatedPruning.clear();
    newPruningPoints.clear();
    //activate or not
    dislower = DisQueryLower1(LID,HID,Neighbor,Label);
    dispeak = DisQueryPeak(LID,HID,Label);
//    cout<<"dis "<<dis<<" "<<oldW<<endl;

    if(dispeak<=oldW){
        if(ifDebug){
            cout<<"Not triggered! "<<LID<<" "<<HID<<": "<<dislower<<" "<<dispeak<<" "<<oldW<<" "<<newW<<endl;
        }
        return;
    }
    if((Label[LID].find(HID)!=Label[LID].end()) && (dislower>oldW)){//index update triggered
        if(ifDebug){
            cout<<"Triggered! "<<Label[LID][HID]<<endl;
        }

        if(ifDebug){
            if(Label[lid].find(hid) != Label[lid].end()){//if found
                cout<<"Before update. Label("<<lid<<","<<hid<<"): "<<Label[lid][hid]<<" "<<DisQueryVally(lid,hid,Neighbor,Label)<<" "<<DisQueryPeak(lid,hid,Label)<<" "<<DijkstraCore(lid,hid)<<endl;
                if(PruningPointSet2[lid].find(hid) != PruningPointSet2[lid].end()){
                    cout<<"Pruning point "<<lid<<" "<<hid<<" "<<PruningPointSet2[lid][hid]<<endl;
                }
            }else{//if not found
                cout<<"Label("<<lid<<","<<hid<<"): -1 "<<DisQueryVally(lid,hid,Neighbor,Label)<<" "<<DisQueryPeak(lid,hid,Label)<<" "<<DijkstraCore(lid,hid)<<endl;
                if(PruningPointSet2[lid].find(hid) != PruningPointSet2[lid].end()){
                    cout<<"Pruning point "<<lid<<" "<<hid<<" "<<PruningPointSet2[lid][hid]<<endl;
                }
            }
        }

        set<tuple<vertex,vertex,vertex>> outdatedPruning;//<nodeID,supportNode,prunedID>
        map<pair<vertex,vertex>,vertex> newPruningPoints;//<<nodeID,prunedID>,supportNode>
        outdatedPruning.clear();
        newPruningPoints.clear();
        NoSupportedPair.clear();

        queue<pair<int,pair<int,int>>> WaitPro;
        queue<pair<int,pair<int,int>>> WaitProP;
        vector<pair<int,int>> AL1; AL1.clear();
        vector<pair<int,int>> AL2; AL2.clear();
        vector<pair<int,int>> AL2Check; AL2Check.clear();
        ChangedLabels.assign(node_num,unordered_set<int>());

        int curID, hubID, hubDis;
        int dis, cnt;



        /// weight change source 1, Label(a,b)
//        vallyPair = DisQueryVallyVert2(LID,HID,Neighbor,Label);
        vallyPair = DisQueryVally2(LID,HID,Neighbor,Label);
        disvally=vallyPair.first, vallyID=vallyPair.second;
        peakPair = DisQueryPeak2(LID,HID,Label);
        dispeak=peakPair.first, peakhub=peakPair.second;

//        if(disvally >= dispeak){
//            cout<<"! dispeak <= disvally for L("<<LID<<","<<HID<<"): "<<disvally<<" "<<dispeak<<endl;
//        }

        if(ifDebug){
            int Dijk= DijkstraCore(LID,HID);
            if(Dijk != disvally){
                cout<<"ab "<<LID<<" "<<HID<<": "<<Label[LID][HID]<<"("<<oldW<<") "<<disvally<<"("<<vallyID<<") "<<dispeak<<"("<<peakhub<<") "<<Dijk<<endl;
            }
        }

        if(Label[LID][HID] < disvally){
            AL1.emplace_back(LID,HID);
//            WaitPro.push(make_pair(LID, make_pair(HID,Label[LID][HID])));
//            Label[LID][HID]=disvally;//correct to the new value
//            ChangedLabels[LID].insert(HID);
        }

//        cout<<"Begin iterative propagation..."<<endl;
        int iterations=1;

        while(!AL1.empty()){
            if(iterations > 1){
                cout<<"Iteration "<<iterations<<endl;
            }
            iterations++;

            //update AL1: AL1->AL1
            CoarseUpdate(LID, HID, oldW, WaitPro, WaitProP, AL1, AL2, AL2Check, Neighbor, Label, ifDebug, lid, hid);
//            cout<<"After AL1: "<<AL1.size()<<" "<<AL2Check.size()<<" "<<WaitPro.size()<<" "<<WaitProP.size()<< endl;

            if(ifDebug){
                vallyPair= DisQueryVally2(lid,hid,Neighbor,Label);
                disvally=vallyPair.first, vallyID=vallyPair.second;
                peakPair= DisQueryPeak2(lid,hid,Label);
                dispeak=peakPair.first, peakhub=peakPair.second;
                if(Label[lid].find(hid) != Label[lid].end()){
                    cout<<"Between AL1 and AL2. "<<lid<<" "<<hid<<": "<<Label[lid][hid]<<" "<<disvally<<"("<<vallyID<<") "<< dispeak<<"("<<peakhub<<") "<< DijkstraCore(lid,hid)<< endl;
                }else{
                    cout<<"Between AL1 and AL2. "<<lid<<" "<<hid<<": -1 "<<disvally<<"("<<vallyID<<") "<< dispeak<<"("<<peakhub<<") "<< DijkstraCore(lid,hid)<< endl;
                }
            }

            //update AL2: AL1->AL2, AL2->AL2
            RefineUpdate( WaitPro, WaitProP, AL1, AL2, AL2Check, outdatedPruning, newPruningPoints, Neighbor, Label, PruningPointNew, ifDebug, lid, hid);
//            cout<<"After AL2: "<<AL1.size()<<" "<<AL2Check.size()<<" "<<WaitPro.size()<<" "<<WaitProP.size()<<endl;
        }
//        cout<<"Done."<<endl;

//        cout<<"Remove and add pruning point."<<endl;
        PPRClean(Neighbor, newPruningPoints, outdatedPruning, ifDebug, lid, hid);
/*        int id1, hub, id2;
        ///remove old pruning point
        for(auto it2=outdatedPruning.begin();it2!=outdatedPruning.end();++it2){
            id1 = get<0>(*it2); id2 = get<2>(*it2); hub = get<1>(*it2);

            disvally=DisQueryVally(id1, id2,Neighbor,Label);
            peakPair=DisQueryPeak2(id1, id2,Label);
            dispeak=peakPair.first; peakhub=peakPair.second;

            if((disvally < dispeak)){
                if(ifDebug){
                    if(lid == id1 && hid == id2){
                        cout<<"Remove pruning point: ("<<id1<<","<<id2<<"): hub "<<hub<<", "<<dispeak<<" "<<disvally<<endl;
                    }
                }
                if(PruningPointSet2[id1].find(id2) != PruningPointSet2[id1].end()){
                    PruningPointSet2[id1].erase(id2);
                }
                if(PruningPointSet2[id2].find(id1) != PruningPointSet2[id2].end()){
                    PruningPointSet2[id2].erase(id1);
                }
                if(PruningPointNew[id1].find(hub) != PruningPointNew[id1].end()){//if found
                    if(PruningPointNew[id1][hub].find(id2) != PruningPointNew[id1][hub].end() ){//if found
                        PruningPointNew[id1][hub].erase(id2);
                        if(PruningPointNew[id1][hub].empty()){
                            PruningPointNew[id1].erase(hub);
                        }
                    }
                }
                if(PruningPointNew[id2].find(hub) != PruningPointNew[id2].end()){//if found
                    if(PruningPointNew[id2][hub].find(id1) != PruningPointNew[id2][hub].end() ){//if found
                        PruningPointNew[id2][hub].erase(id1);
                        if(PruningPointNew[id2][hub].empty()){
                            PruningPointNew[id2].erase(hub);
                        }
                    }
                }
            }else {//if disvally >= dispeak
//                if(hub != peakhub){
//                    cout<<"hub is inconsistent for outdatedPruning point ("<<id1<<" "<<id2<<") ! "<<hub<<" "<<peakhub<<endl;
//                }
                PruningPointNew[id1][peakhub].insert(id2);
                PruningPointNew[id2][peakhub].insert(id1);
                PruningPointSet2[id1][id2]=peakhub;
                PruningPointSet2[id2][id1]=peakhub;
                if(Label[id1].find(id2) != Label[id1].end()){//if found
//                    cout<<"Outdated Pruning point erase wrong. "<<id1<<" "<<id2<<" "<<Label[id1][id2]<<" "<<disvally<<" "<<dispeak<<"("<<peakhub<<","<<hub<<") "<<DijkstraCore(id1,id2)<<endl;
                    Label[id1].erase(id2);
                }
            }
        }
        outdatedPruning.clear();
        ///add new pruning point
        for(auto it=newPruningPoints.begin();it!=newPruningPoints.end();++it){
            id1 = it->first.first; id2 = it->first.second; hub = it->second;
            vallyPair=DisQueryVally2(id1, id2,Neighbor,Label);
            disvally=vallyPair.first; vallyID=vallyPair.second;
            peakPair=DisQueryPeak2(id1, id2,Label);
            dispeak=peakPair.first; peakhub=peakPair.second;
            if(disvally < dispeak){//it is normal for this case as we add all peak pairs to newPruningPoint which may not be correct peak pairs
                if((Label[id1].find(id2) == Label[id1].end()) || (Label[id1].find(id2) != Label[id1].end() && Label[id1][id2] != disvally)) {//if not found or found but incorrect
                    if(ifDebug){
                        int disDijk = Dijkstra(id1,id2,Neighbor);
                        if(disvally == disDijk){
                            if(Label[id1].find(id2) != Label[id1].end()){//if found
                                cout<<"add label by new pruning point. L("<<id1<<","<<id2<<") "<<Label[id1][id2]<<" "<<disvally<<endl;
                            }else{//if not found
                                cout<<"add label by new pruning point. L("<<id1<<","<<id2<<") -1 "<<disvally<<endl;
                            }
                            //                    Label[id1][id2] = disvally;
                        }else if(dispeak == disDijk){
                            cout<<"!! dispeak is correct. "<<Label[id1][id2]<<" "<<disvally<<" "<<dispeak;
                            if(Label[id1].find(id2) != Label[id1].end()) {//if found
                                cout<<" !! erase Label("<<id1<<","<<id2<<")";
                                Label[id1].erase(id2);
                            }
                            cout<<endl;
                        }else{

                            if(Label[id1].find(id2) != Label[id1].end()){
                                cout<<"!!!!!!!! Totally wrong !!!!!!!! L("<<id1<<","<<id2<<"): " <<Label[id1][id2]<<" "<<disvally<<" "<<dispeak<<" "<<disDijk<<endl;
                            }else{
                                cout<<"!!!!!!!! Totally wrong !!!!!!!! L("<<id1<<","<<id2<<"): -1 "<<disvally<<" "<<dispeak<<" "<<disDijk<<endl;
                            }
                        }
                    }
                    if(Label[id1].find(id2) == Label[id1].end()) {//if not found or found but incorrect
                        cout<<"missing label of outdatedPruning. L("<<id1<<" "<<id2<<"): -1 "<<disvally<<" "<<dispeak<<" "<<endl;
                    }else if(Label[id1].find(id2) != Label[id1].end() && Label[id1][id2] != disvally){
                        cout<<"incorrect label of outdatedPruning. L("<<id1<<" "<<id2<<"): "<<Label[id1][id2]<<" "<<disvally<<" "<<dispeak<<" "<<endl;
                    }
                }

            }else{//if disvally >= dispeak
//                if(hub != peakhub){
//                    cout<<"hub is inconsistent for new pruning point ("<<id1<<" "<<id2<<") ! "<<hub<<" "<<peakhub<<endl;
//                }
                PruningPointNew[id1][peakhub].insert(id2);
                PruningPointNew[id2][peakhub].insert(id1);
                PruningPointSet2[id1][id2]=peakhub;
                PruningPointSet2[id2][id1]=peakhub;
                if(Label[id1].find(id2) != Label[id1].end()) {//if found
//                    cout<<"New Pruning point erase wrong. "<<id1<<" "<<id2<<" "<<Label[id1][id2]<<" "<<disvally<<" "<<dispeak<<"("<<peakhub<<","<<hub<<") "<<DijkstraCore(id1,id2)<<endl;
                    Label[id1].erase(id2);
                }

            }

        }
        /// remove redundant label
        //cout<<"remove redundant labels..."<<endl;
        //        CleanLabel(Label);
        vector<pair<int,int>> eraseLabels; eraseLabels.clear();
        for(int i=node_num-1;i>=0;--i){
            id1 = vNodeOrder[i];
            if(CoreTag[id1] != -1){
//            cout<<"break "<<i<<" "<<id1<<endl;
                break;
            }
            if(!ChangedLabels[id1].empty()){
                for(auto it=ChangedLabels[id1].begin();it!=ChangedLabels[id1].end();++it){
                    id2 = *it;
                    vallyPair=DisQueryVally2(id1, id2,Neighbor,Label);
                    disvally=vallyPair.first; vallyID=vallyPair.second;
                    peakPair=DisQueryPeak2(id1, id2,Label);
                    dispeak=peakPair.first; peakhub=peakPair.second;
                    if(dispeak <= disvally){
                        PruningPointNew[id1][peakhub].insert(id2);
                        PruningPointNew[id2][peakhub].insert(id1);
                        PruningPointSet2[id1][id2]=peakhub;
                        PruningPointSet2[id2][id1]=peakhub;
                        eraseLabels.emplace_back(id1,id2);
//                    if(Label[id1].find(id2) != Label[id1].end()){//if found
////                        cout<<"erase redundant label ("<<id1<<" "<<id2<<"): "<<Label[id1][id2]<<" "<<dispeak<<endl;
//                        Label[id1].erase(id2);
//                    }
                    }
                    else{//if disvally < dispeak
                        if(PruningPointSet2[id1].find(id2) != PruningPointSet2[id1].end()){
                            PruningPointSet2[id1].erase(id2);
                        }
                        if(PruningPointSet2[id2].find(id1) != PruningPointSet2[id2].end()){
                            PruningPointSet2[id2].erase(id1);
                        }
                        if(PruningPointNew[id1].find(hub) != PruningPointNew[id1].end()){//if found
                            if(PruningPointNew[id1][hub].find(id2) != PruningPointNew[id1][hub].end() ){//if found
                                PruningPointNew[id1][hub].erase(id2);
                                if(PruningPointNew[id1][hub].empty()){
                                    PruningPointNew[id1].erase(hub);
                                }
                            }
                        }
                        if(PruningPointNew[id2].find(hub) != PruningPointNew[id2].end()){//if found
                            if(PruningPointNew[id2][hub].find(id1) != PruningPointNew[id2][hub].end() ){//if found
                                PruningPointNew[id2][hub].erase(id1);
                                if(PruningPointNew[id2][hub].empty()){
                                    PruningPointNew[id2].erase(hub);
                                }
                            }
                        }
                        if(Label[id1].find(id2) == Label[id1].end()){//if not found
                            cout<<"missing label! L("<<id1<<","<<id2<<"): -1 "<<disvally<<endl;
//                        Label[id1][id2]=disvally;
                        }else if(Label[id1][id2] != disvally){
                            cout<<"incorrect vally label! L("<<id1<<","<<id2<<"): "<<Label[id1][id2]<<" "<<disvally<<endl;
//                        Label[id1][id2]=disvally;
                        }
                    }
                }
            }
        }
        if(!eraseLabels.empty()){
            cout<<"eraseLabels is not empty! ";
            bool flag_erase=false;
            for(auto it2=eraseLabels.begin();it2!=eraseLabels.end();++it2){
                id1 = it2->first, id2 = it2->second;
                if(Label[id1].find(id2) != Label[id1].end()){//if found
                    Label[id1].erase(id2);
                    flag_erase = true;
                }
            }
//            cout<<"label erased."<<endl;
        }*/

    }
    else{
        if(ifDebug){
            cout<<"Not triggered! "<<LID<<" "<<HID<<": "<<dislower<<" "<<dispeak<<" "<<oldW<<" "<<newW<<endl;
        }
    }

}

//function of propagating the AL1 labels, output the updated labels for AL2 check, queue version, 2023-03-26, new
void Graph::CoarseUpdate(int LID,int HID, int oldW, queue<pair<int,pair<int,int>>>& WaitPro, queue<pair<int,pair<int,int>>>& WaitProP, vector<pair<int,int>>& AL1, vector<pair<int,int>>& AL2, vector<pair<int,int>>& AL2Check, vector<vector<pair<vertex,int>>> &Neighbor,vector<unordered_map<vertex,int>> &Label, bool ifDebug, int lid, int hid){
    int disvally,dispeak,peakhub,vallyID;
    pair<int,int> peakPair;//distance, hubID
    pair<int,int> vallyPair;//distance, vallyID

    int curID,hubID,hubDis;
    int neiID, neiDis;

    int lid2,hid2;
    lid2=147082, hid2=196720;

    /// weight change source 1, Label(a,b)
//    vallyPair = DisQueryVallyVert2(LID,HID,Neighbor,Label);
    vallyPair = DisQueryVally2(LID,HID,Neighbor,Label);
    disvally=vallyPair.first, vallyID=vallyPair.second;
    peakPair = DisQueryPeak2(LID,HID,Label);
    dispeak=peakPair.first, peakhub=peakPair.second;

//    if(disvally >= dispeak){
//        cout<<"! dispeak <= disvally for L("<<LID<<","<<HID<<"): "<<disvally<<" "<<dispeak<<endl;
//    }

    if(ifDebug){
        int Dijk= DijkstraCore(LID,HID);
        if(Dijk != disvally){
            cout<<"ab "<<LID<<" "<<HID<<": "<<Label[LID][HID]<<"("<<oldW<<") "<<disvally<<"("<<vallyID<<") "<<dispeak<<"("<<peakhub<<") "<<Dijk<<endl;
        }
    }

    if(Label[LID][HID] < disvally){
//        AL1.emplace_back(LID,HID);
        WaitPro.push(make_pair(LID, make_pair(HID,Label[LID][HID])));
        Label[LID][HID]=disvally;//correct to the new value
        ChangedLabels[LID].insert(HID);
        AL2Check.emplace_back(LID,HID);
    }

    /// weight change source 2, the labels that are directly affected by the e(a,b), using asynchronous propagation
    //Detect the affected labels of HID
    for(auto it=Label[LID].begin();it!=Label[LID].end();it++){
        hubID=(*it).first; hubDis=(*it).second;
        if(Label[HID].find(hubID)!=Label[HID].end() && oldW+hubDis==Label[HID][hubID]){
            AL1.emplace_back(HID,hubID);
            vallyPair= DisQueryVally2(HID,hubID,Neighbor,Label);
            disvally=vallyPair.first, vallyID= vallyPair.second;

            if(Label[HID][hubID] < disvally){
                if(ifDebug){
                    if((lid == HID && hid ==hubID) || (lid2 == HID && hid2 == hubID)){
                        cout<<"Find. "<<HID<<" "<<hubID<<" "<<Label[HID][hubID]<<" "<<disvally<<endl;
                    }
                }
                WaitPro.push(make_pair(HID, make_pair(hubID,Label[HID][hubID])));

                Label[HID][hubID]=disvally;//correct to the new value
                ChangedLabels[HID].insert(hubID);
                //AL1->AL2, check PPR
                AL2Check.emplace_back(HID,hubID);
            }
        }
    }
    //Detect the affected labels of LID
    for(auto it=Label[HID].begin();it!=Label[HID].end();it++){
        hubID=(*it).first; hubDis=(*it).second;
        if(hubID!=HID && Label[LID].find(hubID)!=Label[LID].end() && oldW+hubDis==Label[LID][hubID]){
            AL1.emplace_back(LID,hubID);
            vallyPair= DisQueryVally2(LID,hubID,Neighbor,Label);
            disvally=vallyPair.first, vallyID= vallyPair.second;

            if(Label[LID][hubID] < disvally){
                if(ifDebug){
                    if((lid == LID && hid ==hubID) || (lid2 == LID && hid2 == hubID)){
                        cout<<"Find. "<<LID<<" "<<hubID<<" "<<Label[LID][hubID]<<" "<<disvally<<endl;
                    }
                }
                WaitPro.push(make_pair(LID, make_pair(hubID,Label[LID][hubID])));

                Label[LID][hubID]=disvally;//correct to the new value
                ChangedLabels[LID].insert(hubID);
                //AL1->AL2, check PPR
                AL2Check.emplace_back(LID,hubID);
            }
        }
    }

    //Update the source labels of AL1
//    for(auto it=AL1.begin();it!=AL1.end();++it){
//        curID = it->first, hubID = it->second;
//        vallyPair= DisQueryVally2(curID,hubID,Neighbor,Label);
//        disvally=vallyPair.first, vallyID= vallyPair.second;
//
//        if((Label[curID].find(hubID) != Label[curID].end()) &&  (Label[curID][hubID] < disvally)){
//            if(ifDebug){
//                if((lid == curID && hid ==hubID) || (lid2 == curID && hid2 == hubID)){
//                    cout<<"Find. "<<curID<<" "<<hubID<<" "<<Label[curID][hubID]<<" "<<disvally<<endl;
//                }
//            }
//            WaitPro.push(make_pair(curID, make_pair(hubID,Label[curID][hubID])));
//
//            Label[curID][hubID]=disvally;//correct to the new value
//            ChangedLabels[curID].insert(hubID);
//            //AL1->AL2, check PPR
//            AL2Check.emplace_back(curID,hubID);
//        }
//    }
    AL1.clear();

    while(!WaitPro.empty()){
        queue<pair<int,pair<int,int>>> WaitProTem;
        while(!WaitPro.empty()){
            auto front = WaitPro.front();
            curID = front.first, hubID = front.second.first, hubDis = front.second.second;
            WaitPro.pop();
//            cout<<"AL("<<curID<<","<<hubID<<"): "<<hubDis<<endl;

            //lower-order vertex
            for (int k = 0; k < Neighbor[curID].size(); k++) {//for each neighbor of curID which has the hub hID
                neiID = Neighbor[curID][k].first;
                neiDis = Neighbor[curID][k].second;
//                if (ifDebug) {
//                        if (neiID == lid && hubID == hid) {
//                            cout << "AL1->AL1 detection. curID " << curID << " L(" << neiID << "," << hubID << ")"
//                                 << ": ";
//                            if (Label[neiID].find(hubID) != Label[neiID].end()) {//if found
//                                cout << Label[neiID][hubID] << " ";
//                            }
//                            cout << DisQueryVallyDebug(neiID, hubID, Neighbor) << " "<<DijkstraCore(neiID,hubID)<< endl;
//                        }
//                }
                //detect and update the affected label of neighbor
                if ((Label[neiID].find(hubID) != Label[neiID].end()) && (neiDis + hubDis == Label[neiID][hubID])) {
//                        AL1.insert(OrderComp3(neiID, hubID));
                    vallyPair = DisQueryVally2(neiID, hubID, Neighbor, Label);
                    disvally = vallyPair.first, vallyID = vallyPair.second;
                    peakPair = DisQueryPeak2(neiID, hubID, Label);
                    dispeak = peakPair.first, peakhub = peakPair.second;

                    if (ifDebug) {
                        if ((neiID == lid && hubID == hid) || (neiID == lid2 && hubID == hid2)) {
                            if (Label[neiID].find(hubID) != Label[neiID].end()) {//if found
                                cout << "!!!!!! AL1->AL1 lower. L(" << neiID << "," << hubID << "): " << Label[neiID][hubID] << " "<<disvally << "(" << vallyID << ") " << dispeak << "(" << peakhub << ") "<< DijkstraCore(neiID, hubID) << endl;
                            }else{
                                cout << "!!!!!! AL1->AL1 lower. L(" << neiID << "," << hubID << "): -1 "<<disvally << "(" << vallyID << ") " << dispeak << "(" << peakhub << ") "<< DijkstraCore(neiID, hubID) << endl;
                            }
                        }
                    }
                    //update the label
                    if (Label[neiID][hubID] < disvally) {//original setting
                        WaitProTem.push(make_pair(neiID, make_pair(hubID,Label[neiID][hubID])));

                        Label[neiID][hubID] = disvally;
                        ChangedLabels[neiID].insert(hubID);
                        //AL1->AL2, check PPR
//                        AL2Check.emplace_back(neiID,hubID);
                    }
                    //AL1->AL2, check PPR
                    AL2Check.emplace_back(neiID,hubID);
                }
            }

        }

        WaitPro = WaitProTem;
    }

}
//function of propagating the AL2 labels, output the updated labels for AL1 check, queue version, new
void Graph::RefineUpdate(queue<pair<int,pair<int,int>>>& WaitPro, queue<pair<int,pair<int,int>>>& WaitProP, vector<pair<int,int>>& AL1, vector<pair<int,int>>& AL2, vector<pair<int,int>>& AL2Check, set<tuple<vertex,vertex,vertex>>& outdatedPruning, map<pair<vertex,vertex>,vertex>& newPruningPoints, vector<vector<pair<vertex,int>>> &Neighbor,vector<unordered_map<vertex,int>> &Label, vector<unordered_map<vertex,unordered_set<vertex>>> &PruningPointNew, bool ifDebug, int lid, int hid){
    int disvally,dispeak,peakhub,vallyID;
    pair<int,int> peakPair;//distance, hubID
    pair<int,int> vallyPair;//distance, vallyID

    int curID,hubID,hubDis;

    bool ifUpdate = false;

    //check and get the source AL2, AL1->AL2
    for(int i=0;i<AL2Check.size();++i){//for each increased label
        curID = AL2Check[i].first, hubID = AL2Check[i].second;

        bool temp_bool = PPRCheck(curID,hubID,Neighbor,Label,PruningPointNew,WaitProP,AL2,newPruningPoints,outdatedPruning,ifDebug,lid,hid);///
        if(temp_bool){
            ifUpdate = true;
        }

//        WaitProP.push(make_pair(curID, make_pair(hubID,-1)));///

//        if(ifUpdate){//if update
//            vallyPair=DisQueryVally2(curID,hubID,Neighbor,Label);
//            disvally=vallyPair.first, vallyID=vallyPair.second;
//            peakPair=DisQueryPeak2(curID,hubID,Label);
//            dispeak=peakPair.first, peakhub=peakPair.second;
//
//            if(disvally < dispeak){
//                int Dijk = DijkstraCore(curID,hubID);
//                if(Dijk != Label[curID][hubID]){
//                    cout<<"Incorrect L("<<curID<<","<<hubID<<"): "<<Label[curID][hubID]<<" "<<disvally<<"("<<vallyID<<") "<< dispeak<<"("<<peakhub<<") "<< Dijk<<endl;
//                }
//            }
//            cout<<"new insertion caused by L("<<curID<<","<<hubID<<") "<<Label[curID][hubID]<<" "<<Dijk<<endl;
//        }

    }

    AL2Check.clear();
    AL1.clear();
    if(ifDebug){
        if(ifUpdate){
            cout<<"!! Have new labels inserted!"<<endl;
        }
        if(Label[lid].find(hid) != Label[lid].end()){
            cout<<lid<<" "<<hid<<": "<<Label[lid][hid]<<endl;
        }
        cout<<"AL2->AL2"<<endl;
    }

    int round=0;
    int u,neiid,neiw,cDis;

    while(!WaitProP.empty()){
        queue<pair<int,pair<int,int>>> WaitProPTem;
        while(!WaitProP.empty()){
//        cout<<"round "<<round<<endl;
//        round++;
            auto front = WaitProP.front();
            curID = front.first, u = front.second.first; //hubID=front.second.second;
            WaitProP.pop();

            //detect the affected neighbors of AL2
            //Lower-order vertex
            for(int l = 0; l < Neighbor[curID].size(); l++) {
                neiid = Neighbor[curID][l].first;
                neiw = Neighbor[curID][l].second;
                if (NodeOrder[neiid] < NodeOrder[u]) {
                    vallyPair = DisQueryVally2(neiid, u, Neighbor, Label);
                    disvally = vallyPair.first; vallyID = vallyPair.second;
                    peakPair = DisQueryPeak2(neiid, u, Label);
                    dispeak = peakPair.first; peakhub = peakPair.second;//dispeak may be incorrect for the time being
                    bool flag = false;
                    if (PruningPointSet2[neiid].find(u) != PruningPointSet2[neiid].end()) {
                        hubID = PruningPointSet2[neiid][u];//original peak hub
                        flag = true;
                    } else {
                        hubID = -1;
                    }

                    if (ifDebug) {
                        if (u == hid && neiid == lid) {
                            if (Label[neiid].find(u) != Label[neiid].end()) {//if found
                                cout << "!!! here lower. " << curID << " L(" << neiid << "," << u << "): " << Label[neiid][u] << " "<<disvally << "(" << vallyID << ") " << dispeak << "(" << peakhub << ") " << DijkstraCore(neiid, u) << endl;
                            }else{
                                cout << "!!! here lower. " << curID << " L(" << neiid << "," << u << "): -1 "<<disvally << "(" << vallyID << ") " << dispeak << "(" << peakhub << ") " << DijkstraCore(neiid, u) << endl;
                            }

                        }
                    }

                    if (Label[neiid].find(u) == Label[neiid].end()) {// if not found
                        if (disvally < dispeak) {//if disvally is smaller, add new label
//                        WaitProPTem.insert(neiid);
//                        ChangePTem[neiid].insert(u);
                            WaitProPTem.push(make_pair(neiid, make_pair(u,-1)));

                            AL2.emplace_back(neiid, u);
                            Label[neiid][u] = disvally;
                            ChangedLabels[neiid].insert(u);
                            NoSupportedPair.insert(make_pair(neiid, u));///
//                       AL2Checked.insert(make_pair(neiid,u));
//                       cout<<"AL2->AL2 new insertion. L("<<neiid<<","<<u<<"): "<<Label[neiid][u]<<" "<<disvally<<"("<<vallyID<<") "<<dispeak<<"("<<peakhub<<") "<<DijkstraCore(neiid,u)<<endl;

                            if (hubID != -1) {
                                outdatedPruning.insert(make_tuple(neiid, hubID, u));//
                            } else {//hub == -1, if there is no peak hub in PruningPointSet2, but the label from u to neiid may truly pruned
                                outdatedPruning.insert(make_tuple(neiid, peakhub, u));//
                            }

                            if (PruningPointSet2[neiid].find(u) != PruningPointSet2[neiid].end()) {///
                                PruningPointSet2[neiid].erase(u);
                            }
                            if (PruningPointSet2[u].find(neiid) != PruningPointSet2[u].end()) {
                                PruningPointSet2[u].erase(neiid);
                            }
                        }
                        else {//if dispeak<=disvally
                            if (peakhub != -1 && peakhub != hubID) {//may add redundant PruningPoints
                                newPruningPoints[make_pair(neiid,u)] = peakhub;// this may not be true as the dispeak may be wrong
                                if (hubID != -1) {
                                    outdatedPruning.insert(make_tuple(neiid, hubID, u));//
                                }
                                PruningPointSet2[neiid][u] = peakhub;
                                PruningPointSet2[u][neiid] = peakhub;
                                PruningPointNew[neiid][peakhub].insert(u);
                                PruningPointNew[u][peakhub].insert(neiid);
                            }
                        }

                    }
                    else if ((Label[neiid].find(u) != Label[neiid].end()) && (Label[neiid][u] != disvally)){//if found, check whether the new label L(curID,u) will affect L(neiid,u)
                        //old solution
//                                    if(hubID != -1){
//                                        AL2.insert(OrderComp3(neiid,u));
//                                        WaitProPTem.insert(neiid);///AL2
//                                        ChangePTem[neiid].insert(u);
//                                    }else{
//                                        WaitProTem.insert(neiid);///AL1
//                                        ChangeTem[neiid].insert(OrderComp2(u, Label[neiid][u]));
//                                    }
//                                    Label[neiid][u]=disvally;
//                                    ChangedLabels[neiid].insert(u);

//                        if(ifDebug){
//                            cout<<"AL2->AL2 lower existing label. "<<curID<<"("<<NodeOrder[curID]<<") L("<<neiid<<","<<u<<"): "<<Label[neiid][u]<<" "<<disvally<<"("<<vallyID<<") "<<dispeak<<"("<<peakhub<<") "<<DijkstraCore(neiid,u)<<endl;
//                        }

                        if(Label[neiid][u] > disvally){//for decrease update
                            AL2.emplace_back(neiid, u);
                            WaitProPTem.push(make_pair(neiid, make_pair(u,-1)));///AL2
                            Label[neiid][u]=disvally;
                            ChangedLabels[neiid].insert(u);
                        }else if(Label[neiid][u] < disvally){//for increase update
                            AL1.emplace_back(neiid,u);
                            cout<<"!!! Caused by original redundant labels !! L<d_1 for Refine update (lower)! L("<<neiid<<","<<u<<"): "<<Label[neiid][u]<<" "<<disvally<<" "<<DijkstraCore(neiid,u)<< endl;
                        }

                    }

                }
            }

            //Higher-order vertex, not necessary
            /*for(int l=0; l<Neighbor[u].size();++l){
                neiid = Neighbor[u][l].first;
                neiw = Neighbor[u][l].second;
                if (NodeOrder[curID] < NodeOrder[neiid]) {
                    vallyPair = DisQueryVally2(curID, neiid, Neighbor, Label);
                    disvally = vallyPair.first; vallyID = vallyPair.second;
                    peakPair = DisQueryPeak2(curID, neiid, Label);
                    dispeak = peakPair.first; peakhub = peakPair.second;//dispeak may be incorrect for the time being
                    bool flag = false;
                    if (PruningPointSet2[curID].find(neiid) != PruningPointSet2[curID].end()) {
                        hubID = PruningPointSet2[curID][neiid];//original peak hub
                        flag = true;
                    } else {
                        hubID = -1;
                    }

                    if (ifDebug) {
                        if (neiid == hid && curID == lid) {
                            if (Label[curID].find(neiid) != Label[curID].end()) {//if found
                                cout << "!!! here higher. " << curID << " L(" << curID << "," << neiid << "): " << Label[curID][neiid] << " "<< disvally << "(" << vallyID << ") " << dispeak << "(" << peakhub << ") " << DijkstraCore(curID, neiid) << endl;
                            }else{
                                cout << "!!! here higher. " << curID << " L(" << curID << "," << neiid << "): -1 "<< disvally << "(" << vallyID << ") " << dispeak << "(" << peakhub << ") " << DijkstraCore(curID, neiid) << endl;
                            }
                        }
                    }

                    if (Label[curID].find(neiid) == Label[curID].end()) {// if not found
                        if (disvally < dispeak) {//if disvally is smaller, add new label
                            WaitProPTem.push(make_pair(curID, make_pair(neiid,-1)));

                            AL2.emplace_back(curID,neiid);
                            Label[curID][neiid] = disvally;
                            ChangedLabels[curID].insert(neiid);
                            NoSupportedPair.insert(make_pair(curID,neiid));

                            if (hubID != -1) {
                                outdatedPruning.insert(make_tuple(curID, hubID, neiid));//
                            } else {//hub == -1, if there is no peak hub in PruningPointSet2, but the label from u to neiid may truly pruned
                                outdatedPruning.insert(make_tuple(curID, peakhub, neiid));//
                            }

                            if (PruningPointSet2[neiid].find(curID) != PruningPointSet2[neiid].end()) {///
                                PruningPointSet2[neiid].erase(curID);
                            }
                            if (PruningPointSet2[curID].find(neiid) != PruningPointSet2[curID].end()) {
                                PruningPointSet2[curID].erase(neiid);
                            }
                        }
                        else {//if dispeak<=disvally
                            if (peakhub != -1 && peakhub != hubID) {//may add redundant PruningPoints
                                newPruningPoints[make_pair(curID,neiid)] = peakhub;// this may not be true as the dispeak may be wrong
                                if (hubID != -1) {
                                    outdatedPruning.insert(make_tuple(curID, hubID, neiid));//
                                }
                                PruningPointSet2[neiid][curID] = peakhub;
                                PruningPointSet2[curID][neiid] = peakhub;
                                PruningPointNew[neiid][peakhub].insert(curID);
                                PruningPointNew[curID][peakhub].insert(neiid);
                            }
                        }

                    }
                    else {//if found, check whether the new label L(curID,u) will affect L(neiid,u)

                        if (Label[curID][neiid] != disvally) {
//                            if(ifDebug){
//                                cout<<"AL2->AL2 higher existing label. "<<curID<<"("<<NodeOrder[curID]<<") L("<<neiid<<","<<u<<"): "<<Label[curID][neiid]<<" "<<disvally<<"("<<vallyID<<") "<<dispeak<<"("<<peakhub<<") "<<DijkstraCore(curID,neiid)<<endl;
//                            }

                            if(Label[curID][neiid] > disvally){//for decrease update
                                AL2.emplace_back(curID,neiid);
                                WaitProPTem.push(make_pair(curID, make_pair(neiid,-1)));///AL2
                                Label[curID][neiid]=disvally;
                                ChangedLabels[curID].insert(neiid);
                            }else if(Label[curID][neiid] < disvally){//for increase update
                                AL1.emplace_back(curID,neiid);
                                cout<<"!!! Caused by original redundant labels !! L<d_1 for Refine update (higher)! L("<<curID<<","<<neiid<<"): "<<Label[curID][neiid]<<" "<<disvally<<" "<<DijkstraCore(curID,neiid)<< endl;

                            }
                        }
                    }

                }
            }*/

        }
        WaitProP=WaitProPTem;
    }

}
//function of checking the pruned label by PPR, pair version, queue version, correct
bool Graph::PPRCheck(int curID, int hubID, vector<vector<pair<vertex,int>>> &Neighbor,vector<unordered_map<vertex,int>> &Label, vector<unordered_map<vertex,unordered_set<vertex>>> &PruningPointNew, queue<pair<int,pair<int,int>>>& WaitProPTem, vector<pair<int,int>> &AL2, map<pair<vertex,vertex>,vertex> &newPruningPoints, set<tuple<vertex,vertex,vertex>> &outdatedPruning, bool ifDebug, int lid, int hid){
    int disvally,dispeak,peakhub,vallyID;
    pair<int,int> peakPair;//distance, hubID
    pair<int,int> vallyPair;//distance, vallyID
    bool flag = false;

//    if(ifDebug){
//        if(curID==211099 && hubID==211425){
//            cout<<"PPR source find. L("<<curID<<","<<hubID<<"): "<<Label[curID][hubID]<<endl;
//        }
//    }

    if(PruningPointNew[curID].find(hubID)!=PruningPointNew[curID].end()){
        for(auto snum=PruningPointNew[curID][hubID].begin();snum!=PruningPointNew[curID][hubID].end();++snum){//for each pruned vertex
            int s=*snum;//snum->ID;
            if((NodeOrder[s]<NodeOrder[curID]) && (NoSupportedPair.find(make_pair(s,curID))==NoSupportedPair.end())){//not found
//            if(NodeOrder[s]<NodeOrder[curID]){//not found
//                AL2.insert(OrderComp3(s,curID));///
                vallyPair=DisQueryVally2(s,curID,Neighbor,Label);
                disvally=vallyPair.first; vallyID=vallyPair.second;//disvally may not be correct if there exist neighbor which is also the pruning point of (curID,hID)
                peakPair=DisQueryPeak2(s,curID,Label);//the original dispeak is smaller than disvally
                dispeak=peakPair.first; peakhub=peakPair.second;

                if(ifDebug){
                    if(s == lid && curID == hid){
                        cout<<"!!!!!!!!!!!!!!!!!!!!!!!! target pruned label L("<<s<<","<<curID<<") by "<<hubID<<": "<<disvally<<"("<<vallyID<<") "<<dispeak<<"("<<peakhub<<") "<<DijkstraCore(s,curID)<<endl;
                    }
                }

                if(disvally<dispeak){//the pruned label should be inserted back
//                    WaitProPTem.insert(s);
//                    ChangePTem[s].insert(curID);
                    WaitProPTem.push(make_pair(s, make_pair(curID,hubID)));

                    if(ifDebug){
//                        int Dijk= DijkstraCore(s,curID);
//                        if(Dijk != disvally){
//                            if(Label[s].find(curID) != Label[s].end()){//if found
//                                cout<<"AL1->AL2. L("<<s<<","<<curID<<") "<<Label[s][curID]<<" "<<disvally<<"("<<vallyID<<") "<<dispeak<<"("<<peakhub<<") "<<Dijk<<endl;
//                            }else{//if not found
//                                cout<<"AL1->AL2. L("<<s<<","<<curID<<") -1 "<<disvally<<"("<<vallyID<<") "<<dispeak<<"("<<peakhub<<") "<<Dijk<<endl;
//                            }
//                        }
                    }

                    AL2.emplace_back(s,curID);
                    Label[s][curID]=disvally;
                    ChangedLabels[s].insert(curID);
                    flag = true;

                    outdatedPruning.insert(make_tuple(s,hubID,curID));
                    NoSupportedPair.insert(make_pair(s,curID));

                    if(PruningPointSet2[curID].find(s) != PruningPointSet2[curID].end()){///
                        PruningPointSet2[curID].erase(s);
                    }
                    if(PruningPointSet2[s].find(curID) != PruningPointSet2[s].end()){
                        PruningPointSet2[s].erase(curID);
                    }

                }
                else {//if dispeak<=disvally

                    if(peakhub != -1 && peakhub != hubID){
                        newPruningPoints[make_pair(s,curID)] = peakhub;
                        outdatedPruning.insert(make_tuple(s,hubID,curID));
                        PruningPointSet2[curID][s]=peakhub;
                        PruningPointSet2[s][curID]=peakhub;
                        PruningPointNew[curID][peakhub].insert(s);
                        PruningPointNew[s][peakhub].insert(curID);
                    }
                }
            }else if(NodeOrder[curID]<NodeOrder[s] && NoSupportedPair.find(make_pair(curID,s))==NoSupportedPair.end()){
//            } else if(NodeOrder[s]>NodeOrder[curID]){
//                AL2.insert(OrderComp3(curID,s));///
                vallyPair=DisQueryVally2(curID,s,Neighbor,Label);
                disvally=vallyPair.first; vallyID=vallyPair.second;
                peakPair=DisQueryPeak2(curID,s,Label);//the old P2 distance dominates P1 distance, but the increase of d(curID,hID)
                dispeak=peakPair.first; peakhub=peakPair.second;

                if(ifDebug){
                    if(s == hid && curID == lid){
                        cout<<"!!!!!!!!!!!!!!!!!!!!!!!! target pruned label L("<<curID<<","<<s<<") by "<<hubID<<": "<<disvally<<"("<<vallyID<<") "<<dispeak<<"("<<peakhub<<") "<<DijkstraCore(s,curID)<<endl;
                    }
                }

                if(disvally<dispeak){//disvally may be incorrect
//                    WaitProPTem.insert(curID);///
//                    ChangePTem[curID].insert(s);
                    WaitProPTem.push(make_pair(curID, make_pair(s,hubID)));

                    if(ifDebug){
//                        int Dijk= DijkstraCore(s,curID);
//                        if(Dijk != disvally){
//                            if(Label[curID].find(s) != Label[curID].end()){//if found
//                                cout<<"AL1->AL2. L("<<curID<<","<<s<<") "<<Label[curID][s]<<" "<<disvally<<"("<<vallyID<<") "<<dispeak<<"("<<peakhub<<") "<<Dijk<<endl;
//                            }else{//if not found
//                                cout<<"AL1->AL2. L("<<curID<<","<<s<<") -1 "<<disvally<<"("<<vallyID<<") "<<dispeak<<"("<<peakhub<<") "<<Dijk<<endl;
//                            }
//                        }
                    }

                    AL2.emplace_back(curID,s);
                    Label[curID][s]=disvally;//the label here may be incorrect for the time being
                    ChangedLabels[curID].insert(s);
                    flag = true;

                    outdatedPruning.insert(make_tuple(curID,hubID,s));
                    NoSupportedPair.insert(make_pair(curID,s));

                    if(PruningPointSet2[curID].find(s) != PruningPointSet2[curID].end()){
                        PruningPointSet2[curID].erase(s);
                    }
                    if(PruningPointSet2[s].find(curID) != PruningPointSet2[s].end()){
                        PruningPointSet2[s].erase(curID);
                    }
                }
                else {//if dispeak<=disvally

                    if(peakhub != -1 && peakhub != hubID) {
                        newPruningPoints[make_pair(curID,s)] = peakhub;
                        outdatedPruning.insert(make_tuple(curID, hubID, s));///
                        PruningPointSet2[curID][s]=peakhub;
                        PruningPointSet2[s][curID]=peakhub;
                        PruningPointNew[curID][peakhub].insert(s);
                        PruningPointNew[s][peakhub].insert(curID);
                    }
                }
            }
        }
    }
    return flag;
}
//function for cleaning redundant labels and PPR
void Graph::PPRClean(vector<vector<pair<vertex,int>>> &Neighbor, map<pair<vertex,vertex>,vertex> &newPruningPoints, set<tuple<vertex,vertex,vertex>> &outdatedPruning, bool ifDebug, int lid, int hid){
    int dislower,disvally,dispeak;
    vertex peakhub,vallyID;
    pair<int,vertex> peakPair;//distance, hubID
    pair<int,vertex> vallyPair;//distance, vallyID
    int id1, hub, id2;
    ///remove old pruning point
    for(auto it2=outdatedPruning.begin();it2!=outdatedPruning.end();++it2){
        id1 = get<0>(*it2); id2 = get<2>(*it2); hub = get<1>(*it2);

        disvally=DisQueryVally(id1, id2,Neighbor,Label);
        peakPair=DisQueryPeak2(id1, id2,Label);
        dispeak=peakPair.first; peakhub=peakPair.second;

        if((disvally < dispeak)){
            if(ifDebug){
                if(lid == id1 && hid == id2){
                    cout<<"Remove pruning point: ("<<id1<<","<<id2<<"): hub "<<hub<<", "<<dispeak<<" "<<disvally<<endl;
                }
            }
            if(PruningPointSet2[id1].find(id2) != PruningPointSet2[id1].end()){
                PruningPointSet2[id1].erase(id2);
            }
            if(PruningPointSet2[id2].find(id1) != PruningPointSet2[id2].end()){
                PruningPointSet2[id2].erase(id1);
            }
            if(PruningPointSet[id1].find(hub) != PruningPointSet[id1].end()){//if found
                if(PruningPointSet[id1][hub].find(id2) != PruningPointSet[id1][hub].end() ){//if found
                    PruningPointSet[id1][hub].erase(id2);
                    if(PruningPointSet[id1][hub].empty()){
                        PruningPointSet[id1].erase(hub);
                    }
                }
            }
            if(PruningPointSet[id2].find(hub) != PruningPointSet[id2].end()){//if found
                if(PruningPointSet[id2][hub].find(id1) != PruningPointSet[id2][hub].end() ){//if found
                    PruningPointSet[id2][hub].erase(id1);
                    if(PruningPointSet[id2][hub].empty()){
                        PruningPointSet[id2].erase(hub);
                    }
                }
            }
        }else {//if disvally >= dispeak
//                if(hub != peakhub){
//                    cout<<"hub is inconsistent for outdatedPruning point ("<<id1<<" "<<id2<<") ! "<<hub<<" "<<peakhub<<endl;
//                }
            PruningPointSet[id1][peakhub].insert(id2);
            PruningPointSet[id2][peakhub].insert(id1);
            PruningPointSet2[id1][id2]=peakhub;
            PruningPointSet2[id2][id1]=peakhub;
            if(Label[id1].find(id2) != Label[id1].end()){//if found
//                    cout<<"Outdated Pruning point erase wrong. "<<id1<<" "<<id2<<" "<<Label[id1][id2]<<" "<<disvally<<" "<<dispeak<<"("<<peakhub<<","<<hub<<") "<<DijkstraCore(id1,id2)<<endl;
                Label[id1].erase(id2);
            }
        }
    }
    outdatedPruning.clear();
    ///add new pruning point
    for(auto it=newPruningPoints.begin();it!=newPruningPoints.end();++it){
        id1 = it->first.first; id2 = it->first.second; hub = it->second;
        vallyPair=DisQueryVally2(id1, id2,Neighbor,Label);
        disvally=vallyPair.first; vallyID=vallyPair.second;
        peakPair=DisQueryPeak2(id1, id2,Label);
        dispeak=peakPair.first; peakhub=peakPair.second;
        if(disvally < dispeak){//it is normal for this case as we add all peak pairs to newPruningPoint which may not be correct peak pairs
            if((Label[id1].find(id2) == Label[id1].end()) || (Label[id1].find(id2) != Label[id1].end() && Label[id1][id2] != disvally)) {//if not found or found but incorrect
                if(ifDebug){
                    int disDijk = Dijkstra(id1,id2,Neighbor);
                    if(disvally == disDijk){
                        if(Label[id1].find(id2) != Label[id1].end()){//if found
                            cout<<"add label by new pruning point. L("<<id1<<","<<id2<<") "<<Label[id1][id2]<<" "<<disvally<<endl;
                        }else{//if not found
                            cout<<"add label by new pruning point. L("<<id1<<","<<id2<<") -1 "<<disvally<<endl;
                        }
                        //                    Label[id1][id2] = disvally;
                    }else if(dispeak == disDijk){
                        cout<<"!! dispeak is correct. "<<Label[id1][id2]<<" "<<disvally<<" "<<dispeak;
                        if(Label[id1].find(id2) != Label[id1].end()) {//if found
                            cout<<" !! erase Label("<<id1<<","<<id2<<")";
                            Label[id1].erase(id2);
                        }
                        cout<<endl;
                    }else{

                        if(Label[id1].find(id2) != Label[id1].end()){
                            cout<<"!!!!!!!! Totally wrong !!!!!!!! L("<<id1<<","<<id2<<"): " <<Label[id1][id2]<<" "<<disvally<<" "<<dispeak<<" "<<disDijk<<endl;
                        }else{
                            cout<<"!!!!!!!! Totally wrong !!!!!!!! L("<<id1<<","<<id2<<"): -1 "<<disvally<<" "<<dispeak<<" "<<disDijk<<endl;
                        }
                    }
                }
                if(Label[id1].find(id2) == Label[id1].end()) {//if not found or found but incorrect
                    cout<<"missing label of outdatedPruning. L("<<id1<<" "<<id2<<"): -1 "<<disvally<<" "<<dispeak<<" "<<endl;
                }else if(Label[id1].find(id2) != Label[id1].end() && Label[id1][id2] != disvally){
                    cout<<"incorrect label of outdatedPruning. L("<<id1<<" "<<id2<<"): "<<Label[id1][id2]<<" "<<disvally<<" "<<dispeak<<" "<<endl;
                }
            }

        }else{//if disvally >= dispeak
//                if(hub != peakhub){
//                    cout<<"hub is inconsistent for new pruning point ("<<id1<<" "<<id2<<") ! "<<hub<<" "<<peakhub<<endl;
//                }
            PruningPointSet[id1][peakhub].insert(id2);
            PruningPointSet[id2][peakhub].insert(id1);
            PruningPointSet2[id1][id2]=peakhub;
            PruningPointSet2[id2][id1]=peakhub;
            if(Label[id1].find(id2) != Label[id1].end()) {//if found
//                    cout<<"New Pruning point erase wrong. "<<id1<<" "<<id2<<" "<<Label[id1][id2]<<" "<<disvally<<" "<<dispeak<<"("<<peakhub<<","<<hub<<") "<<DijkstraCore(id1,id2)<<endl;
                Label[id1].erase(id2);
            }

        }

    }
    /// remove redundant label
    //cout<<"remove redundant labels..."<<endl;
    //        CleanLabel(Label);
    vector<pair<int,int>> eraseLabels; eraseLabels.clear();
    for(int i=node_num-1;i>=0;--i){
        id1 = vNodeOrder[i];
        if(!ChangedLabels[id1].empty()){
            for(auto it=ChangedLabels[id1].begin();it!=ChangedLabels[id1].end();++it){
                id2 = *it;
                vallyPair=DisQueryVally2(id1, id2,Neighbor,Label);
                disvally=vallyPair.first; vallyID=vallyPair.second;
                peakPair=DisQueryPeak2(id1, id2,Label);
                dispeak=peakPair.first; peakhub=peakPair.second;
                if(dispeak <= disvally){
                    PruningPointSet[id1][peakhub].insert(id2);
                    PruningPointSet[id2][peakhub].insert(id1);
                    PruningPointSet2[id1][id2]=peakhub;
                    PruningPointSet2[id2][id1]=peakhub;
                    eraseLabels.emplace_back(id1,id2);
//                    if(Label[id1].find(id2) != Label[id1].end()){//if found
////                        cout<<"erase redundant label ("<<id1<<" "<<id2<<"): "<<Label[id1][id2]<<" "<<dispeak<<endl;
//                        Label[id1].erase(id2);
//                    }
                }
                else{//if disvally < dispeak
                    if(PruningPointSet2[id1].find(id2) != PruningPointSet2[id1].end()){
                        PruningPointSet2[id1].erase(id2);
                    }
                    if(PruningPointSet2[id2].find(id1) != PruningPointSet2[id2].end()){
                        PruningPointSet2[id2].erase(id1);
                    }
                    if(PruningPointSet[id1].find(hub) != PruningPointSet[id1].end()){//if found
                        if(PruningPointSet[id1][hub].find(id2) != PruningPointSet[id1][hub].end() ){//if found
                            PruningPointSet[id1][hub].erase(id2);
                            if(PruningPointSet[id1][hub].empty()){
                                PruningPointSet[id1].erase(hub);
                            }
                        }
                    }
                    if(PruningPointSet[id2].find(hub) != PruningPointSet[id2].end()){//if found
                        if(PruningPointSet[id2][hub].find(id1) != PruningPointSet[id2][hub].end() ){//if found
                            PruningPointSet[id2][hub].erase(id1);
                            if(PruningPointSet[id2][hub].empty()){
                                PruningPointSet[id2].erase(hub);
                            }
                        }
                    }
                    if(Label[id1].find(id2) == Label[id1].end()){//if not found
                        cout<<"missing label! L("<<id1<<","<<id2<<"): -1 "<<disvally<<endl;
//                        Label[id1][id2]=disvally;
                    }else if(Label[id1][id2] != disvally){
                        cout<<"incorrect vally label! L("<<id1<<","<<id2<<"): "<<Label[id1][id2]<<" "<<disvally<<endl;
//                        Label[id1][id2]=disvally;
                    }
                }
            }
        }
    }
    if(!eraseLabels.empty()){
//        cout<<"eraseLabels is not empty! ";
        bool flag_erase=false;
        for(auto it2=eraseLabels.begin();it2!=eraseLabels.end();++it2){
            id1 = it2->first, id2 = it2->second;
            if(Label[id1].find(id2) != Label[id1].end()){//if found
                Label[id1].erase(id2);
                flag_erase = true;
            }
        }
//            cout<<"label erased."<<endl;
    }
}

//original version of WPSL, incorrect
void Graph::IncreasePSL(int a, int b, int oldW, int newW, vector<vector<pair<vertex,int>>> &Neighbors,vector<unordered_map<vertex,int>> &Label,vector<unordered_map<int,vector<int>>> &PruningPointNew,set<pair<int,int>> &NoSupportedPair){
    for(int i=0;i<Neighbors[a].size();i++){
        if(Neighbors[a][i].first==b){
//            cout<<i<<" "<<a<<" "<<b<<" "<<Neighbors[a][i].second<<" "<<newW<<endl;
            if(oldW != Neighbors[a][i].second){
                cout<<"Old edge weight is incorrect! "<<a<<" "<<b<<": "<<oldW<<" "<<Neighbors[a][i].second<<endl;
                oldW = Neighbors[a][i].second;
            }
            Neighbors[a][i].second=newW;
            break;
        }
    }
    for(int i=0;i<Neighbors[b].size();i++){
        if(Neighbors[b][i].first==a){
//            cout<<i<<" "<<a<<" "<<b<<" "<<Neighbors[b][i].second<<" "<<newW<<endl;
            if(oldW != Neighbors[b][i].second){
                cout<<"Old edge weight is incorrect! "<<b<<" "<<a<<": "<<oldW<<" "<<Neighbors[b][i].second<<endl;
                oldW = Neighbors[b][i].second;
            }
            Neighbors[b][i].second=newW;
            break;
        }
    }

    bool ifDebug = false;//false;
//    ifDebug = true;
    int lid,hid,lid2,hid2;
    lid=71972, lid2=72339, hid=27296;
    lid=223905, lid2=223904, hid=230142;
    lid=54337, lid2=54394, hid=27784, hid2=230142;

    int LID,HID;
    if(NodeOrder[a]>NodeOrder[b]){
        LID=b; HID=a;
    }else{
        LID=a; HID=b;
    }
//    cout<<LID<<"("<<NodeOrder[LID]<<") "<<HID<<"("<<NodeOrder[HID]<<")"<<endl;
    pair<int,int> peakPair, valleyPair;
    int dis,disvally,valleyID,dispeak,peakhub;

    //activate or not
    dis=DisQueryLower1(LID,HID,Neighbors,Label);//d1', via a neighbor of LID (which has lower order than HID) to reach HID
    dispeak=DisQueryPeak(LID,HID,Label);//d2, via the common hub vertex of LID and HID
    if(dispeak<=oldW)//if d2 is lower than oldW, the increase of oldW will not affect the shortest distance between LID and HID
        return;
    if(ifDebug){
        if(dis == INF)
            cout<<DisQueryLower1(LID,HID,Neighbors,Label)<<endl;
        if(dis==oldW){
            cout<<dis<<" "<<oldW<<endl;
        }
    }

    if(dis>oldW){//if d1' is larger than oldW, it indicates the shortest distance between LID and HID is equal to oldW, index update triggered
        vector<vector<pair<int,int>>> Change;//the label that is changed
        vector<pair<int,int>> vec;
        Change.assign(node_num,vec);
        set<int> WaitPro;//the vertex waited for label propagation, i.e., AL1
        vector<vector<int>> ChangeP;
        ChangeP.assign(node_num,vector<int>());
        set<int> WaitProP;//the vertex waited for label propagation, i.e., AL2

        WaitPro.insert(LID);
        Change[LID].push_back(make_pair(HID, oldW));
        disvally=DisQueryVally(LID,HID,Neighbors,Label);//shortest distance query through the neighbor's label
        Label[LID][HID]=disvally;//may not be the final correct value at this moment

        //cout<<"start "<<LID<<" "<<HID<<" "<<disvally<<endl;

        //affected by the w(a,b)
        int hubID, hDis;
        int dis, cnt;
        for(auto it=Label[HID].begin();it!=Label[HID].end();++it){//check the vertex has higher order than HID and update LID's label (except HID)
            hubID=(*it).first; hDis=(*it).second;
            if(hubID!=HID && Label[LID].find(hubID)!=Label[LID].end()){// && oldW+hDis==Label[LID][hubID]
                if(oldW+hDis==Label[LID][hubID]){//if the shortest path pass through e(a,b)
                    disvally=DisQueryVally(LID,hubID,Neighbors,Label);
                    if(Label[LID][hubID]<disvally){//if the new distance of P1 is higher, update it
                        WaitPro.insert(LID);
                        Change[LID].push_back(make_pair(hubID, oldW+hDis));//record the changed label of LID with the old distance
                        Label[LID][hubID]=disvally;//should be correct
                        //cout<<"weight "<<LID<<" "<<hubID<<" "<<disvally<<endl;
                    }
                }

            }
        }

        for(auto it=Label[LID].begin();it!=Label[LID].end();it++){//update HID's label
            hubID=(*it).first; hDis=(*it).second;
            if(Label[HID].find(hubID)!=Label[HID].end()){//&& oldW+hDis==Label[HID][hubID]
                if(oldW+hDis==Label[HID][hubID]){
                    disvally=DisQueryVally(HID,hubID,Neighbors,Label);
                    if(Label[HID][hubID]<disvally){
                        WaitPro.insert(HID);
                        Change[HID].push_back(make_pair(hubID, oldW+hDis));//record the changed label
                        Label[HID][hubID]=disvally;//should be correct
                        //cout<<"weight "<<HID<<" "<<hubID<<" "<<disvally<<endl;
                    }
                }
            }
        }

        while(WaitProP.size()>0 || WaitPro.size()>0){
            set<int> WaitProTem;
            vector<vector<pair<int,int>>> ChangeTem;
            vector<pair<int,int>> vec;
            ChangeTem.assign(node_num,vec);
            set<int> WaitProPTem;
            vector<int> vecint;
            vector<vector<int>> ChangePTem;
            ChangePTem.assign(node_num,vecint);

            //Change->Change & ChangeP: AL1->AL1 and AL1->AL2
            for(auto it=WaitPro.begin();it!=WaitPro.end();it++){
                int curID=*it;
                vector<pair<int,int>> curChange=Change[curID];
                int neiID, neiDis, hID, hDis;


                for(int j=0;j<curChange.size();j++){
                    hID=curChange[j].first; hDis=curChange[j].second;
                    //Change->Change: AL1->AL1
                    for(int k=0;k<Neighbors[curID].size();k++){
                        neiID=Neighbors[curID][k].first; neiDis=Neighbors[curID][k].second;
                        if(Label[neiID].find(hID)!=Label[neiID].end() ){//&& neiDis+hDis==Label[neiID][hID]
                            if(neiDis+hDis==Label[neiID][hID]){

                                valleyPair=DisQueryVally2(neiID,hID,Neighbors,Label);
                                disvally=valleyPair.first, valleyID=valleyPair.second;
                                if(Label[neiID][hID]<disvally){
                                    if(ifDebug){
                                        if((neiID == lid || neiID == lid2) && hID == hid){
                                            cout<<"AL1->AL1: "<<neiID<<" "<<hID<<" "<<Label[neiID][hID]<<" "<<disvally<<"("<<valleyID<<") "<<DijkstraCore(neiID,hID)<<endl;
                                        }
                                    }

                                    WaitProTem.insert(neiID);
                                    ChangeTem[neiID].push_back(make_pair(hID, neiDis+hDis));
                                    Label[neiID][hID]=disvally;//may not be the final correct value at this moment


                                    //cout<<"1--1 "<<neiID<<" "<<hID<<" "<<disvally<<endl;
                                }
                            }

                        }
                    }

                    if(ifDebug){
                        if((curID == lid || curID == lid2 ) && hID == hid){
                            if(PruningPointNew[curID].find(hID)!=PruningPointNew[curID].end()){
                                cout<<"Affected Label find! "<<curID<<" "<<hID<<" "<<PruningPointNew[curID][hID].size()<<endl;
                            }else {
                                cout << "Affected Label find! " << curID << " " << hID << endl;
                            }
                        }
                    }

                    //Change->ChangeP: AL1->AL2
                    if(PruningPointNew[curID].find(hID)!=PruningPointNew[curID].end()){//check whether the pruned vertex should be inserted to curID again

                        for(int snum=0;snum<PruningPointNew[curID][hID].size();snum++){
                            int s=PruningPointNew[curID][hID][snum];

                            if(NodeOrder[s]<NodeOrder[curID] && NoSupportedPair.find(make_pair(s,curID))==NoSupportedPair.end()){//it is not in NoSupportedPair
                                disvally=DisQueryVally(s,curID,Neighbors,Label);

                                peakPair=DisQueryPeak2(s,curID,Label);
                                dispeak=peakPair.first; //peakhub=peakPair.second;
                                if(dispeak>disvally){
                                    WaitProPTem.insert(s);
                                    ChangePTem[s].push_back(curID);
                                    Label[s][curID]=disvally;
                                    NoSupportedPair.insert(make_pair(s,curID));
                                    if(ifDebug){
                                        if((s == lid || s==lid2) && curID == hid){
                                            cout<<"AL1->AL2 1: "<<s<<" "<<curID<<" "<<disvally<<endl;
                                        }
                                    }


                                    //cout<<"1--2 "<<s<<" "<<curID<<" "<<disvally<<endl;
                                }
                            }else if(NodeOrder[s]>NodeOrder[curID] && NoSupportedPair.find(make_pair(curID,s))==NoSupportedPair.end()){
                                disvally=DisQueryVally(curID,s,Neighbors,Label);
                                peakPair=DisQueryPeak2(curID,s,Label);//the old P2 distance dominates P1 distance, but the increase of d(curID,hID) may violate the previous balance
                                dispeak=peakPair.first; //peakhub=peakPair.second;
                                if(dispeak>disvally){
                                    WaitProPTem.insert(curID);
                                    ChangePTem[curID].push_back(s);
                                    Label[curID][s]=disvally;//should be the final correct value
                                    NoSupportedPair.insert(make_pair(curID,s));//insert to indicate that there is no pruning between curID and s
                                    if(ifDebug){
                                        if((curID == lid || curID==lid2) && s == hid){
                                            cout<<"AL1->AL2 2: "<<curID<<" "<<s<<" "<<disvally<<endl;
                                        }
                                    }

                                    //cout<<"1--2 "<<curID<<" "<<s<<" "<<disvally<<endl;
                                }
                            }
                        }
                    }
                }
            }

            //ChangeP->CHangeP: AL2->AL2
            int v,u,neiid,neiw;
            for(auto itp=WaitProP.begin();itp!=WaitProP.end();itp++){
                v=*itp;
                for(int k=0;k<ChangeP[v].size();k++){
                    u=ChangeP[v][k];
                    for(int l=0;l<Neighbors[v].size();l++){
                        neiid=Neighbors[v][l].first; neiw=Neighbors[v][l].second;
                        if(NodeOrder[neiid]<NodeOrder[u]){
                            valleyPair=DisQueryVally2(neiid, u,Neighbors,Label);
                            peakPair=DisQueryPeak2(neiid, u,Label);
                            disvally=valleyPair.first, valleyID=valleyPair.second;
                            dispeak=peakPair.first, peakhub=peakPair.second;

                            if(ifDebug){
                                if((neiid == lid || neiid==lid2) && u == hid){
                                    if(Label[neiid].find(u)!=Label[neiid].end()){
                                        cout<<"AL2->AL2: "<<neiid<<" "<<u<<" "<<Label[neiid][u]<<" "<<disvally<<"("<<valleyID<<") "<<dispeak<<"("<<peakhub<<") "<<DijkstraCore(neiid,u)<<endl;
                                    }else{
                                        cout<<"AL2->AL2: "<<neiid<<" "<<u<<" "<<disvally<<"("<<valleyID<<") "<<dispeak<<"("<<peakhub<<") "<<DijkstraCore(neiid,u)<<endl;
                                    }

                                }
                            }

                            if(Label[neiid].find(u)==Label[neiid].end()){
                                if(disvally<dispeak) {
                                    WaitProPTem.insert(neiid);
                                    ChangePTem[neiid].push_back(u);
                                    Label[neiid][u]=disvally;//should be the final correct value

                                    NoSupportedPair.insert(make_pair(neiid,u));///
                                }

                            }else if(Label[neiid].find(u)!=Label[neiid].end() && Label[neiid][u]!=disvally){
                                WaitProPTem.insert(neiid);
                                ChangePTem[neiid].push_back(u);
                                Label[neiid][u]=disvally;//should be the final correct value
                            }

                        }
                    }
                }
            }

            WaitPro=WaitProTem;
            Change=ChangeTem;
            WaitProP=WaitProPTem;
            ChangeP=ChangePTem;
        }

    }
}
//function for querying the distance by label. old version
int Graph::PLLQuery(int ID1,int ID2,vector<unordered_map<vertex,int>> &Label){
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
//function for querying the distance by label. vector-based version
/*int Graph::PLLQuery(int ID1,int ID2,vector<unordered_map<vertex,int>> &Label){
    int d=INF;

    int hub, dis1, dis2;
    if(LabelV.Labels[ID1].size() != Label[ID1].size()){
        cout<<"Size inconsistent between Labels and LabelV! "<<LabelV.Labels[ID1].size()<<" "<<Label[ID1].size() <<endl; exit(1);
    }
    for(auto it=LabelV.Labels[ID1].begin();it!=LabelV.Labels[ID1].end();++it){
//    for(auto it=Label[ID1].begin();it!=Label[ID1].end();it++){
        hub=(*it).first;
        if(Label[ID1].find(hub)==Label[ID1].end()){//if not found
            cout<<"Hub inconsistent between Labels and LabelV! "<<endl; exit(1);
        }
        dis1=Label[ID1][hub];
//        dis1=(*it).second;
        if(Label[ID2].find(hub)!=Label[ID2].end()){
            dis2=Label[ID2][hub];
            if(dis1+dis2<d){
                d=dis1+dis2;
                //cout<<"details "<<hub<<" "<<dis1<<" "<<dis2<<" "<<d<<endl;
            }
        }
    }
    return d;
}*/
//function of computing distance from ID1(LID) to ID2(HID), via a neighbor of ID1 (which has lower order than HID) on updated graph
int Graph::DisQueryLower1(int ID1, int ID2, vector<vector<pair<vertex,int>>> &Neighbors,vector<unordered_map<vertex,int>> &Label){
    int neiID,neiDis;
    int d=INF;
    for(int i=0;i<Neighbors[ID1].size();i++){
        neiID=Neighbors[ID1][i].first;
        neiDis=Neighbors[ID1][i].second;
        if(NodeOrder[neiID]<NodeOrder[ID2] && Label[neiID].find(ID2)!=Label[neiID].end()){//if neiID has lower order and ID2 is a hub of neiID
            if(neiDis+Label[neiID][ID2]<d){
                d=neiDis+Label[neiID][ID2];
            }
        }
    }
    return d;
}
//old version. function of computing the shortest distance through the label of ID1's neighbors, i.e., d1.
int Graph::DisQueryVally(int ID1, int ID2, vector<vector<pair<vertex,int>>> &Neighbors,vector<unordered_map<vertex,int>> &Label){
    int neiID,neiDis;
    int d=INF;
    for(int i=0;i<Neighbors[ID1].size();i++){
        neiID=Neighbors[ID1][i].first;
        neiDis=Neighbors[ID1][i].second;
        if(NodeOrder[neiID]<=NodeOrder[ID2] && Label[neiID].find(ID2)!=Label[neiID].end()){//if found
            if(neiDis+Label[neiID][ID2]<d){
                d=neiDis+Label[neiID][ID2];
            }
        }
    }
    return d;
}
//function of computing vally distance, return distance and neighbor
pair<int,int> Graph::DisQueryVally2(int ID1, int ID2, vector<vector<pair<vertex,int>>> &Neighbors,vector<unordered_map<vertex,int>> &Label){
    int neiID,neiDis;
    int d=INF;
    int finalNeigh = -1;
    for(int i=0;i<Neighbors[ID1].size();i++){
        neiID=Neighbors[ID1][i].first;
//        if(neiID == 204062){
//            cout<<ID1<<" "<<neiID<<" "<<ID2<<endl;
//            if(PruningPointSet2[neiID].find(ID2)!=PruningPointSet2[neiID].end()){
//                cout<<"Pruned "<<PruningPointSet2[neiID][ID2]<<endl;
//            }
//        }
        neiDis=Neighbors[ID1][i].second;
        if(NodeOrder[neiID]<=NodeOrder[ID2] && Label[neiID].find(ID2)!=Label[neiID].end()){
            if(neiDis+Label[neiID][ID2]<d){
                d=neiDis+Label[neiID][ID2];
                finalNeigh = neiID;
//                cout<<neiID<<" "<<d<<endl;
            }
        }
    }
    return make_pair(d,finalNeigh);
}
//function of updating the label with the changed edge, 2023-03-10
int Graph::DisQueryVallyVert(int ID1, int ID2, vector<vector<pair<vertex,int>>> &Neighbors,vector<unordered_map<vertex,int>> &Label){
    int neiID,neiDis;
    int d=INF;

    for(int i=0;i<Neighbors[ID1].size();i++){
        neiID=Neighbors[ID1][i].first;
        neiDis=Neighbors[ID1][i].second;
        if(NodeOrder[neiID]<=NodeOrder[ID2] && Label[neiID].find(ID2)!=Label[neiID].end()){//if found
            if(neiDis+Label[ID1][ID2] != Label[neiID][ID2]){//avoid the case that the shortest path between neiID and ID2 pass L(ID1,ID2)
                if(neiDis+Label[neiID][ID2]<d){
                    d=neiDis+Label[neiID][ID2];
                }
            }else{//deal with special case
                int dnei=INF;
                int neiID2,neiDis2;
                for(int j=0;j<Neighbors[neiID].size();++j){
                    neiID2=Neighbors[neiID][j].first;
                    if(neiID2 != ID1 && NodeOrder[neiID2] <= NodeOrder[ID2] && Label[neiID2].find(ID2)!=Label[neiID2].end()){
                        neiDis2=Neighbors[neiID][j].second;
                        if(neiDis2+Label[neiID2][ID2]<dnei){
                            dnei=neiDis2+Label[neiID2][ID2];
                        }
                    }
                }
                if(neiDis+dnei<d){
                    d=neiDis+dnei;
                }
            }
        }
    }
    return d;
}
//function of updating the label with the changed edge, 2023-03-10
pair<int,int> Graph::DisQueryVallyVert2(int ID1, int ID2, vector<vector<pair<vertex,int>>> &Neighbors,vector<unordered_map<vertex,int>> &Label){
    int neiID,neiDis;
    int d=INF;
    int finalNeigh = -1;
    for(int i=0;i<Neighbors[ID1].size();i++){
        neiID=Neighbors[ID1][i].first;
        neiDis=Neighbors[ID1][i].second;
        if(NodeOrder[neiID]<=NodeOrder[ID2] && Label[neiID].find(ID2)!=Label[neiID].end()){//if found
            if(neiDis+Label[ID1][ID2] != Label[neiID][ID2]){//avoid the case that the shortest path between neiID and ID2 pass L(ID1,ID2)
                if(neiDis+Label[neiID][ID2]<d){
                    d=neiDis+Label[neiID][ID2];
                    finalNeigh = neiID;
                }
            }else{
                int dnei=INF;
                int neiID2,neiDis2;
                for(int j=0;j<Neighbors[neiID].size();++j){
                    neiID2=Neighbors[neiID][j].first;
                    if(neiID2 != ID1 && NodeOrder[neiID2] <= NodeOrder[ID2] && Label[neiID2].find(ID2)!=Label[neiID2].end()){
                        neiDis2=Neighbors[neiID][j].second;
                        if(neiDis2+Label[neiID2][ID2]<dnei){
                            dnei=neiDis2+Label[neiID2][ID2];
                        }
                    }
                }
                if(neiDis+dnei<d){
                    d=neiDis+dnei;
                    finalNeigh=neiID;
                }
            }
        }
    }
    return make_pair(d,finalNeigh);
}
//debug version. function of computing the shortest distance through the label of ID1's neighbors, i.e., d1.
int Graph::DisQueryVallyDebug(int ID1, int ID2, vector<vector<pair<vertex,int>>> &Neighbors,vector<unordered_map<vertex,int>> &Label){
    int neiID,neiDis;
    int d=INF;
    for(int i=0;i<Neighbors[ID1].size();i++){
        neiID=Neighbors[ID1][i].first;
        neiDis=Neighbors[ID1][i].second;
        if(NodeOrder[neiID]<=NodeOrder[ID2] && Label[neiID].find(ID2)!=Label[neiID].end()){//if found
            if(neiDis+Label[neiID][ID2]<d){
                d=neiDis+Label[neiID][ID2];
//                cout<<ID1<<" "<<ID2<<": "<<neiID<<" "<<d<<" "<<neiDis<<" "<<Label[neiID][ID2]<<endl;
            }
        }
    }
    return d;
}

//old version: function of computing d2
int Graph::DisQueryPeak(int ID1, int ID2,vector<unordered_map<vertex,int>> &Label){
    int d=INF;

    int hub, dis1;
    for(auto it=Label[ID1].begin();it!=Label[ID1].end();it++){
        hub=(*it).first;
        dis1=(*it).second;
        if(hub!=ID2 && Label[ID2].find(hub)!=Label[ID2].end()){//hub is the hub of ID1 and ID2 and hub != ID2
            if(dis1+Label[ID2][hub]<d){
                d=dis1+Label[ID2][hub];
            }
        }
    }
    return d;
}
//function of computing d2, new version
pair<int,int> Graph::DisQueryPeak2(int ID1, int ID2,vector<unordered_map<vertex,int>> &Label){
    int d=INF;

    int hub, dis1, finalHub=-1;

    for(auto it=Label[ID1].begin();it!=Label[ID1].end();it++){
        hub=(*it).first;
        dis1=(*it).second;
        if(hub!=ID2 && Label[ID2].find(hub)!=Label[ID2].end()){//hub is the hub of ID1 and ID2 and hub != ID2
            if(dis1+Label[ID2][hub]<d){
                d=dis1+Label[ID2][hub];
                finalHub = hub;
            }

        }
    }
    return make_pair(d,finalHub);
}


//function for comparing the PLL label with PSL label
void Graph::CompareLabel(string filename1, string filename2){
    cout<<"Comparing labels..."<<endl;
    ifstream IF1(filename1);
    if(!IF1){
        cout<<"Cannot open file "<<filename1<<endl;
        exit(1);
    }
    ifstream IF2(filename2);
    if(!IF2){
        cout<<"Cannot open file "<<filename2<<endl;
        exit(1);
    }
    vector<unordered_map<vertex,int>> Label1, Label2;//Label1 is the ground truth
    Label1.assign(node_num,unordered_map<vertex,int>());
    Label2.assign(node_num,unordered_map<vertex,int>());
    string line;
    int ID1,ID2,weight;

    //read label 1
    getline(IF1,line);
    vector<string> vs;
    boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
    int node_num=stoi(vs[0]);
    assert(node_num == node_num);
    getline(IF1,line);
    while(!line.empty()){
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
        ID1=stoi(vs[0]), ID2=stoi(vs[1]), weight=stoi(vs[2]);
        Label1[ID1][ID2] = weight;
        if(IF1.eof())
            break;
        getline(IF1,line);
    }
    IF1.close();

    //read label 2
    getline(IF2,line);
    vs.clear();
    boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
    node_num=stoi(vs[0]);
    assert(node_num == node_num);
    getline(IF2,line);
    while(!line.empty()){
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
        ID1=stoi(vs[0]), ID2=stoi(vs[1]), weight=stoi(vs[2]);
        Label2[ID1][ID2] = weight;
        if(IF2.eof())
            break;
        getline(IF2,line);
    }
    IF2.close();

    pair<int,int> peakpair, vallypair;
    int dispeak, peakhub, disvally, vallyid;

    for(int i=node_num-1;i>=0;--i){
        ID1 = vNodeOrder[i];
        if(CoreTag[ID1] != -1){
            break;
        }
        if(Label1[ID1].size() != Label2[ID1].size()){
            cout<<"Inconsistent node "<<ID1<<"("<<NodeOrder[ID1]<<") : "<<Label2[ID1].size()<<" "<<Label1[ID1].size()<<endl;

            if(Label1[ID1].size() < Label2[ID1].size()){
                for(auto it=Label2[ID1].begin();it!=Label2[ID1].end();++it){
                    ID2 = it->first;
                    if(Label1[ID1].find(ID2) == Label1[ID1].end()){//if not found
                        peakpair = DisQueryPeak2(ID1,ID2,Label2);
                        dispeak = peakpair.first, peakhub = peakpair.second;
                        vallypair = DisQueryVally2(ID1,ID2,AdjaCore,Label2);
                        disvally = vallypair.first, vallyid = vallypair.second;
                        int disDijk = DijkstraCore(ID1,ID2);
                        cout<<ID2<<"("<<NodeOrder[ID2]<<") "<<Label2[ID1][ID2]<<" "<<disvally<<"("<<vallyid<<","<<NodeOrder[vallyid]<<") "<<dispeak<<"("<<peakhub<<","<<NodeOrder[peakhub]<<") "<<disDijk<<endl;
                    }
                }
            }else{// if(Label2[ID1].size() <= Label1[ID1].size())
                cout<<"!!!!! Ground truth is larger!"<<endl;
                for(auto it=Label1[ID1].begin();it!=Label1[ID1].end();++it){
                    ID2 = it->first;
                    if(Label2[ID1].find(ID2) == Label2[ID1].end()){//if not found
                        cout<<ID2<<"("<<NodeOrder[ID2]<<","<<Label1[ID1][ID2]<<") ";
                    }
                }
                cout<<endl;
            }

        }
    }
    cout<<"Done."<<endl;
}
//function of cleaning the redundant PLL labels
void Graph::CleanLabel(vector<unordered_map<vertex,int>> &Label){
    int HID,dis;
    int ID;
    pair<int,int> peakpair, vallypair;
    int dispeak, disvally, peakhub, vallyid;
    cout<<"Cleaning label..."<<endl;
    Timer tt;
    tt.start();

    for(int i=node_num-2;i>=0;i--){
        ID = vNodeOrder[i];
        if(CoreTag[ID] != -1){
            break;
        }
        vector<pair<int,int>> eraseIDs;
        for(auto it1=Label[ID].begin();it1!=Label[ID].end();++it1){
            HID = it1->first; dis = it1->second;
            peakpair = DisQueryPeak2(ID,HID,Label);
            dispeak = peakpair.first, peakhub = peakpair.second;
            if(dispeak <= dis){
                eraseIDs.emplace_back(HID,dispeak);
            }
        }
        if(!eraseIDs.empty()){
            for(auto it2=eraseIDs.begin();it2!=eraseIDs.end();++it2){
                HID = it2->first, dispeak = it2->second;
                cout<<"Order "<<i<<": erase label L("<<ID<<","<<HID<<") "<<Label[ID][HID]<<" "<<dispeak<<endl;
                Label[ID].erase(HID);

            }
        }
    }
    tt.stop();
    cout<<"Time used for cleaning label: "<<tt.GetRuntime()<<" s."<<endl;
}
//function of cleaning the redundant PLL label and write the results to filename
void Graph::CleanLabel(vector<unordered_map<vertex,int>> &Label, string filename){
    int HID,dis;
    int ID;
    pair<int,int> peakpair, vallypair;
    int dispeak, disvally, peakhub, vallyid;
    cout<<"Cleaning label..."<<endl;
    Timer tt;
    tt.start();
    ofstream OF(filename+".clean");

    for(int i=node_num-2;i>=0;i--){
        ID = vNodeOrder[i];
        if(CoreTag[ID] != -1){
            break;
        }
        vector<pair<int,int>> eraseIDs;
        for(auto it1=Label[ID].begin();it1!=Label[ID].end();++it1){
            HID = it1->first; dis = it1->second;
            peakpair = DisQueryPeak2(ID,HID,Label);
            dispeak = peakpair.first, peakhub = peakpair.second;
            if(dispeak <= dis){
                eraseIDs.emplace_back(HID,dispeak);
            }
        }
        if(!eraseIDs.empty()){
            for(auto it2=eraseIDs.begin();it2!=eraseIDs.end();++it2){
                HID = it2->first, dispeak = it2->second;
//                cout<<"Order "<<i<<": erase label L("<<ID<<","<<HID<<") "<<Label[ID][HID]<<" "<<dispeak<<endl;
                OF<<"Order "<<i<<": erase label L("<<ID<<","<<HID<<") "<<Label[ID][HID]<<" "<<dispeak<<endl;
                Label[ID].erase(HID);

            }
        }
    }
    OF.close();
    tt.stop();
    cout<<"Time used for cleaning label: "<<tt.GetRuntime()<<" s."<<endl;
}