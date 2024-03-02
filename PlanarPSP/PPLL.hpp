/*
 * PPLL.hpp
 *
 *  Created on: 16 June 2023
 *      Author: Xinjie ZHOU
 */
#include "headPSP.h"
#include "labeling.hpp"


//function for index construction of batch PLL
void Graph::ConstructPLL_OverlayIndex(vector<vector<pair<vertex,int>>> &Neighbor){
    bool ifParallel = true;
    Label.assign(node_num,unordered_map<vertex,int>());

    PruningPointSet2.clear();
    PruningPointSet2.assign(node_num,unordered_map<vertex,vertex>());
    PruningPointSet.clear();
    PruningPointSet.assign(node_num,unordered_map<vertex,unordered_set<vertex>>());

    int ID;
    int cnt=0;
    int stepShow = ceil(nodeNumOverlay/100000)*10000;
    stepShow = max(stepShow,1000);

//    cout<<"Step for show: "<<stepShow<<endl;

    int batchSize = 512;
    batchSize = BPCLBSize;
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

    for(int i=0;i<OverlayVertex.size();++i){
        ID=OverlayVertex[i];
        vertices.emplace_back(ID);
        if(bNodes.size()<batchSize){
            bNodes.emplace_back(ID);
            if(bNodes.size()==batchSize){
                batches.push_back(bNodes);
                bNodes.clear();
            }
        }
    }
    if(!bNodes.empty()){
        batches.push_back(bNodes);
        bNodes.clear();
    }

//    for(int i=node_num-1;i>=0;i--){
//        ID=vNodeOrder[i];
//        if(!PartiTag[ID].second){//if ID is not core vertex
//            if(!bNodes.empty()){
//                batches.push_back(bNodes);
//                bNodes.clear();
//            }
//            break;
//        }
//        vertices.emplace_back(ID);
//        if(bNodes.size()<batchSize){
//            bNodes.emplace_back(ID);
//            if(bNodes.size()==batchSize || i==0){
//                batches.push_back(bNodes);
//                bNodes.clear();
//            }
//        }
//    }

    int batchNum = batches.size();
    int batchShow = 1;
    if (batchShow < batchNum / 10){
        batchShow = batchNum / 10;
    }

    cout<<"Total batch number: "<<batchNum<<endl;
//    cout<<"Batch gap for show: "<<batchShow<<endl;
    int lid, hid;

    Timer tt2;
    tt2.start();
    vertex hID = vNodeOrder[node_num-1];
    for(int b1=0;b1<batches.size();++b1){
        bNodes=batches[b1];
        hID = bNodes[0];
//        cout<<"Batch "<<b1<<": "<<hID<<"("<<NodeOrder[hID]<<") "<<PartiTag[hID].second<<endl;
        unordered_set<vertex> setNodes; setNodes.clear();
        setNodes.insert(bNodes.begin(),bNodes.end());

        vector<pair<int,int>> vp;
        // process each batch
        tt.start();

        /// boost-based implementation
//        if(ifParallel){//use multiple thread
        if(false){
            if(batchSize > threadnum){
                ProcessID.assign(threadnum,vector<vertex>());
                for(int j=0;j<bNodes.size();++j){
                    a = j%threadnum;
                    ProcessID[a].emplace_back(bNodes[j]);
                }
                boost::thread_group thread;
                for(int i=0;i<ProcessID.size();i++){
                    thread.add_thread(new boost::thread(&Graph::BatchPCLDijk2, this, boost::ref(ProcessID[i]), boost::ref(setNodes), hID, boost::ref(Neighbor), boost::ref(Label)));
                }
                thread.join_all();
            }else if(batchSize == threadnum){
                boost::thread_group thread;
                for(int i=0;i<bNodes.size();i++){
                    ID = bNodes[i];
                    thread.add_thread(new boost::thread(&Graph::BatchPCLDijk, this, ID, boost::ref(setNodes), hID, boost::ref(Neighbor), boost::ref(Label)));
                }
                thread.join_all();
            }

        }else{//use single thread
            for(int i=0;i<bNodes.size();i++){
                ID = bNodes[i];
                BatchPCLDijk(ID, setNodes, hID, Neighbor, Label);
//                for(int j=0;j<vp.size();j++){
//                    Label[vp[j].first].insert(make_pair(ID, vp[j].second));
//                    //cout<<vp[j].first<<" "<<vp[j].second<<endl;
//                }
            }
        }

        tt.stop();
        time+=tt.GetRuntime();
//        if(b1%batchShow==0){
//            cout<<"Batch "<<b1<<" finish! Order("<<NodeOrder[bNodes[0]]<<","<<NodeOrder[bNodes[bNodes.size()-1]]<<"). "<<time<<" s."<<endl;
//            time = 0;
//        }

    }
    tt2.stop();
    cout<<"Time used for label construction: "<<tt2.GetRuntime()<<" s."<<endl;

    vector<vector<vertex>> processID;
    processID.assign(threadnum, vector<vertex>());

    ThreadDistribute(vertices, processID);

    tt2.start();
//    cout<<"Begin to construct pruning points..."<<endl;
    PruningPointBuild(ifParallel, processID, Neighbor, Label, PruningPointSet, PruningPointSet2);
//    PruningPointBuild(false);
    tt2.stop();
    cout<<"Time used for pruning point construction: "<<tt2.GetRuntime()<<" s."<<endl;
}

//function for partition index construction for PPLL
void Graph::ConstructPLL_PartiIndex(bool ifParallel) {
    Labels.assign(node_num,unordered_map<vertex,int>());

    PruningPointSetP2.clear();
    PruningPointSetP2.assign(node_num,unordered_map<vertex,vertex>());
    PruningPointSetP.clear();
    PruningPointSetP.assign(node_num,unordered_map<vertex,unordered_set<vertex>>());
    cout<<"Batch size: "<<BPCLBSize<<endl;
    if(ifParallel){
        cout<<"Multiple thread computation for partition index construction!"<<endl;
        //multi-thread
        if(threadnum<partiNum){
            vector<vector<int>> processID;
            processID.assign(threadnum, vector<int>());
            vector<int> vertices;
            for(int pid=0;pid<partiNum;++pid){
                vertices.emplace_back(pid);
            }
            ThreadDistribute(vertices, processID);
            cout<<"Batch number: "<<processID[0].size()<<endl;
            boost::thread_group thread;

            for(auto j=0;j<processID.size();++j){
                thread.add_thread(new boost::thread(&Graph::BPCLIndexConstructPartiV, this, boost::ref(processID[j]), boost::ref(NeighborsParti), boost::ref(Labels), boost::ref(PruningPointSetP), boost::ref(PruningPointSetP2)));
            }


            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                thread.add_thread(new boost::thread(&Graph::BPCLIndexConstructParti, this, j, boost::ref(NeighborsParti), boost::ref(Labels), boost::ref(PruningPointSetP), boost::ref(PruningPointSetP2)));
            }
            thread.join_all();
        }

    }
    else{
        cout<<"Single thread computation!"<<endl;
        //single-thread
        for(int pid=0;pid<partiNum;++pid){
//            cout<<"Partition "<<pid<<endl;
            BPCLIndexConstructParti(pid, NeighborsParti, Labels, PruningPointSetP, PruningPointSetP2);
        }
    }
}

//function for partition index construction for PPLL
void Graph::ConstructPLL_PartiIndexPost(bool ifParallel) {
    LabelsPost.assign(node_num,unordered_map<vertex,int>());

    PruningPointSetPost2.clear();
    PruningPointSetPost2.assign(node_num,unordered_map<vertex,vertex>());
    PruningPointSetPost.clear();
    PruningPointSetPost.assign(node_num,unordered_map<vertex,unordered_set<vertex>>());
//    cout<<"Batch size: "<<BPCLBSize<<endl;
    if(ifParallel){
        cout<<"Multiple thread computation for partition index construction!"<<endl;
        //multi-thread
        if(threadnum<partiNum){
            vector<vector<int>> processID;
            processID.assign(threadnum, vector<int>());
            vector<int> vertices;
            for(int pid=0;pid<partiNum;++pid){
                vertices.emplace_back(pid);
            }
            ThreadDistribute(vertices, processID);
            cout<<"Batch number: "<<processID[0].size()<<endl;
            boost::thread_group thread;

            for(auto j=0;j<processID.size();++j){
                thread.add_thread(new boost::thread(&Graph::BPCLIndexConstructPartiV, this, boost::ref(processID[j]), boost::ref(NeighborsPartiPostV), boost::ref(LabelsPost), boost::ref(PruningPointSetPost), boost::ref(PruningPointSetPost2)));
            }


            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                thread.add_thread(new boost::thread(&Graph::BPCLIndexConstructParti, this, j, boost::ref(NeighborsPartiPostV), boost::ref(LabelsPost), boost::ref(PruningPointSetPost), boost::ref(PruningPointSetPost2)));
            }
            thread.join_all();
        }

    }
    else{
        cout<<"Single thread computation!"<<endl;
        //single-thread
        for(int pid=0;pid<partiNum;++pid){
//            cout<<"Partition "<<pid<<endl;
            BPCLIndexConstructParti(pid, NeighborsPartiPostV, LabelsPost, PruningPointSetPost, PruningPointSetPost2);
        }
    }
}

void Graph::ConstructPLL_OverlayGraph(bool ifParallel){
    if(ifParallel){
        //multiple threads
        if(threadnum<partiNum){
            vector<vector<int>> processID;
            processID.assign(threadnum, vector<int>());
            vector<int> vertices;
            for(int pid=0;pid<partiNum;++pid){
                vertices.emplace_back(pid);
            }
            ThreadDistribute(vertices, processID);
            cout<<"Batch number: "<<processID[0].size()<<endl;
            boost::thread_group thread;
            for(auto j=0;j<processID.size();++j){
                thread.add_thread(new boost::thread(&Graph::ConstructBoundaryShortcutPLLV, this, boost::ref(processID[j])));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                thread.add_thread(new boost::thread(&Graph::ConstructBoundaryShortcutPLL, this, j));
            }
            thread.join_all();
        }

    }
    else{
        //single thread
        for(int k=0;k<partiNum;k++){
            ConstructBoundaryShortcutPLL(k);
        }
    }
}

void Graph::ConstructBoundaryShortcutPLLV(vector<int>& p){
    int pid;
    for(int i=0;i<p.size();++i){
        pid=p[i];
        ConstructBoundaryShortcutPLL(pid);
    }
}

void Graph::ConstructBoundaryShortcutPLL(int pid){
    //boundary edges
    int ID1,ID2,weight;

    for(int i=0;i<BoundVertex[pid].size();i++){
        ID1=BoundVertex[pid][i];
        for(int j=i+1;j<BoundVertex[pid].size();j++){
            ID2=BoundVertex[pid][j];
//            cout<<ID1<<" "<<ID2<<endl;
            weight=QuerySamePartiPLL(ID1,ID2,Labels);
            if(weight>=INF){
                continue;
            }
            NeighborsOverlay[ID1][ID2]=weight;
            NeighborsOverlay[ID2][ID1]=weight;
            bool ifExist=false;
            for(int p=0;p<NeighborsOverlayV[ID1].size();++p){
                if(NeighborsOverlayV[ID1][p].first==ID2){
                    if(NeighborsOverlayV[ID1][p].second>weight){
                        NeighborsOverlayV[ID1][p].second==weight;
                    }
                    ifExist=true;
                    break;
                }
            }
            if(!ifExist){//if not found
                NeighborsOverlayV[ID1].emplace_back(ID2,weight);
                NeighborsOverlayV[ID2].emplace_back(ID1,weight);
            }else{//if found
                ifExist= false;
                for(int p=0;p<NeighborsOverlayV[ID2].size();++p){
                    if(NeighborsOverlayV[ID2][p].first==ID1){
                        if(NeighborsOverlayV[ID2][p].second>weight){
                            NeighborsOverlayV[ID2][p].second==weight;
                        }
                        ifExist=true;
                        break;
                    }
                }
                if(!ifExist){
                    cout<<"Wrong! Inconsistent."<<endl; exit(1);
                }
            }
        }
    }


}

void Graph::ConstructPartitionPostPLL(bool ifParallel){
    NeighborsPartiPost.assign(node_num,unordered_map<int,int>());

    if(ifParallel){
        // multi-thread
//        cout<<"Multi-thread computation!"<<endl;
        //multi-thread
        if(threadnum<partiNum){
            vector<vector<int>> processID;
            processID.assign(threadnum, vector<int>());
            vector<int> vertices;
            for(int pid=0;pid<partiNum;++pid){
                vertices.emplace_back(pid);
            }
            ThreadDistribute(vertices, processID);
//            cout<<"Batch number: "<<processID[0].size()<<endl;
            boost::thread_group thread;
            for(auto j=0;j<processID.size();++j){
                thread.add_thread(new boost::thread(&Graph::ConstructPostPartiPLLV, this, boost::ref(processID[j])));

            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                thread.add_thread(new boost::thread(&Graph::ConstructPostPartiPLL, this, j));

            }
            thread.join_all();
        }
    }
    else{
        // single thread
        for(int k=0;k<partiNum;k++){
            cout<<"Repairing partition "<<k<<endl;
            ConstructPostPartiPLL(k);

        }
    }
}

void Graph::ConstructPostPartiPLLV(vector<int>& p){
    for(int i=0;i<p.size();++i){
        ConstructPostPartiPLL(p[i]);
    }
}

void Graph::ConstructPostPartiPLL(int pid){
    int ID;
    for(int i=0;i<PartiVertex[pid].size();++i){
        ID=PartiVertex[pid][i];
        NeighborsPartiPost[ID].insert(NeighborsParti[ID].begin(),NeighborsParti[ID].end());
    }
    int ID1,ID2,weight=-1;

    for(int i=0;i<BoundVertex[pid].size();++i){
        ID1=BoundVertex[pid][i];
        for(int j=i+1;j<BoundVertex[pid].size();++j){
            ID2=BoundVertex[pid][j];
//            assert(NeighborsOverlay[ID1].find(ID2)!=NeighborsOverlay[ID1].end());//may not be true for no all-pair strategy
//            weight=NeighborsOverlay[ID1][ID2];
//            if(NeighborsPartiPost[ID1].find(ID2)==NeighborsPartiPost[ID1].end()){//if not found
//                NeighborsPartiPost[ID1].insert({ID2,weight});
//                NeighborsPartiPost[ID2].insert({ID1,weight});
//            }
            weight= QueryOverlayPLL(ID1,ID2);
            NeighborsPartiPost[ID1][ID2]=weight;
            NeighborsPartiPost[ID2][ID1]=weight;
            bool ifExist=false;
            for(int p=0;p<NeighborsPartiPostV[ID1].size();++p){
                if(NeighborsPartiPostV[ID1][p].first==ID2){
                    NeighborsPartiPostV[ID1][p].second==weight;
                    ifExist=true;
                    break;
                }
            }
            if(!ifExist){//if not found
                NeighborsPartiPostV[ID1].emplace_back(ID2,weight);
                NeighborsPartiPostV[ID2].emplace_back(ID1,weight);
            }else{
                ifExist= false;
                for(int p=0;p<NeighborsPartiPostV[ID2].size();++p){
                    if(NeighborsPartiPostV[ID2][p].first==ID1){
                        NeighborsPartiPostV[ID2][p].second==weight;
                        ifExist=true;
                        break;
                    }
                }
                if(!ifExist){
                    cout<<"Wrong! Inconsistent."<<endl; exit(1);
                }
            }
        }
    }


}

void Graph::BPCLIndexConstructPartiV(vector<int>& p, vector<vector<pair<vertex, int>>> &Neighbor, vector<unordered_map<vertex,int>>& Label, vector<unordered_map<vertex,unordered_set<vertex>>>& PruningPointSet, vector<unordered_map<vertex,vertex>>& PruningPointSet2) {
    int pid;
    for(int i=0;i<p.size();++i){
        pid=p[i];
        BPCLIndexConstructParti(pid, Neighbor, Label, PruningPointSet, PruningPointSet2);
    }
}

void Graph::BPCLIndexConstructParti(int pid, vector<vector<pair<vertex, int>>> &Neighbor, vector<unordered_map<vertex,int>>& Label, vector<unordered_map<vertex,unordered_set<vertex>>>& PruningPointSet, vector<unordered_map<vertex,vertex>>& PruningPointSet2) {
    bool ifParallel = true;

    int ID;
    int cnt=0;
    int stepShow = ceil(node_num/100000)*10000;
    stepShow = max(stepShow,1000);

//    cout<<"Step for show: "<<stepShow<<endl;

    int batchSize = 512;
    batchSize = BPCLBSize;
    if(batchSize>node_num){
        batchSize=node_num;
    }
//    batchSize = 3;
//    cout<<"Batch size: "<<batchSize<<endl;

    vector<vector<vertex>> batches;
    vector<vertex> bNodes;//vertices for current batch
    vector<vector<vertex>> ProcessID;
    vector<vertex> vertices;

    Timer tt;
    double time=0;
    int a=0;

    for(int i=0;i<PartiVertex[pid].size();++i){
        ID=PartiVertex[pid][i];
//        cout<<i<<": "<<ID<<"("<<NodeOrder[ID]<<") "<<PartiTag[ID].second<<endl;
        vertices.emplace_back(ID);
        if(bNodes.size()<batchSize){
            bNodes.emplace_back(ID);
            if(bNodes.size()==batchSize){
                batches.push_back(bNodes);
                bNodes.clear();
            }
        }
    }
    if(!bNodes.empty()){
        batches.push_back(bNodes);
        bNodes.clear();
    }

    int batchNum = batches.size();
    int batchShow = 1;
    if (batchShow < batchNum / 3){
        batchShow = batchNum / 3;
    }

//    cout<<"Total batch number: "<<batchNum<<endl;
//    cout<<"Batch gap for show: "<<batchShow<<endl;
    int lid, hid;
    int partiThread=1;
    partiThread=max(partiThread, threadnum/partiNum);

    Timer tt2;
    tt2.start();
    vertex hID = vNodeOrder[node_num-1];
    for(int b1=0;b1<batches.size();++b1){
        bNodes=batches[b1];
        hID = bNodes[0];
//        cout<<"Batch "<<b1<<": "<<hID<<"("<<NodeOrder[hID]<<"). bNodes size: "<<bNodes.size()<<endl;
        unordered_set<vertex> setNodes; setNodes.clear();
        setNodes.insert(bNodes.begin(),bNodes.end());

        vector<pair<int,int>> vp;
        // process each batch
        tt.start();

        /// boost-based implementation
        if(ifParallel && partiThread>1){//use multiple thread
//        if(false){
            if(batchSize > partiThread){
                ProcessID.assign(partiThread,vector<vertex>());
                for(int j=0;j<bNodes.size();++j){
                    a = j%partiThread;
                    ProcessID[a].emplace_back(bNodes[j]);
                }
                boost::thread_group thread;
                for(int i=0;i<ProcessID.size();i++){
                    thread.add_thread(new boost::thread(&Graph::BatchPCLDijk2, this, boost::ref(ProcessID[i]), boost::ref(setNodes), hID, boost::ref(Neighbor), boost::ref(Label)));
                }
                thread.join_all();
            }else if(batchSize == partiThread){
                boost::thread_group thread;
                for(int i=0;i<bNodes.size();i++){
                    ID = bNodes[i];
                    thread.add_thread(new boost::thread(&Graph::BatchPCLDijk, this, ID, boost::ref(setNodes), hID, boost::ref(Neighbor), boost::ref(Label)));
                }
                thread.join_all();
            }

        }
        else{//use single thread
            for(int i=0;i<bNodes.size();i++){
                ID = bNodes[i];
                BatchPCLDijk(ID, setNodes, hID, Neighbor, Label);
//                for(int j=0;j<vp.size();j++){
//                    Label[vp[j].first].insert(make_pair(ID, vp[j].second));
//                    //cout<<vp[j].first<<" "<<vp[j].second<<endl;
//                }
            }
        }

        tt.stop();
        time+=tt.GetRuntime();
//        if(b1%batchShow==0){
//            cout<<"Batch "<<b1<<" finish! Order("<<NodeOrder[bNodes[0]]<<","<<NodeOrder[bNodes[bNodes.size()-1]]<<"). "<<time<<" s."<<endl;
//            time = 0;
//        }

    }
    tt2.stop();
//    cout<<"Time used for label construction: "<<tt2.GetRuntime()<<" s."<<endl;

    vector<vector<vertex>> processID;
    processID.assign(threadnum, vector<vertex>());

    ThreadDistribute(vertices, processID);

    tt2.start();
//    cout<<"Begin to construct pruning points..."<<endl;
    PruningPointBuild(ifParallel, processID, Neighbor, Label, PruningPointSet, PruningPointSet2);
//    PruningPointBuild(false);
    tt2.stop();
//    cout<<"Time used for pruning point construction: "<<tt2.GetRuntime()<<" s."<<endl;
}

//function of building the pruning points
void Graph::PruningPointBuild(bool ifParallel, vector<vector<vertex>> & processID, vector<vector<pair<vertex,int>>> &Neighbor, vector<unordered_map<vertex,int>>& Label, vector<unordered_map<vertex,unordered_set<vertex>>>& PruningPointSet, vector<unordered_map<vertex,vertex>>& PruningPointSet2){

    if(ifParallel){
//    if(false){
//        cout<<"Batch size: "<<processID[0].size()<<endl;
        boost::thread_group thread;
        for(auto j=0;j<processID.size();++j){
            thread.add_thread(new boost::thread(&Graph::PPRConstruction2, this, boost::ref(processID[j]), boost::ref(Neighbor), boost::ref(Label), boost::ref(PruningPointSet), boost::ref(PruningPointSet2) ));
        }
        thread.join_all();
    }else{
        for(auto j=0;j<processID.size();++j){
            PPRConstruction2(processID[j], Neighbor, Label, PruningPointSet, PruningPointSet2);
        }
    }
}

void Graph::PPRConstruction2(vector<vertex> & p, vector<vector<pair<vertex,int>>> &Neighbor, vector<unordered_map<vertex,int>>& Label, vector<unordered_map<vertex,unordered_set<vertex>>>& PruningPointSet, vector<unordered_map<vertex,vertex>>& PruningPointSet2){
//    cout<<"1"<<endl;
    for(int i=0;i<p.size();++i){
        PPRConstruction(p[i], Neighbor, Label, PruningPointSet, PruningPointSet2);
    }
}

//function of building the pruning point records for a given vertex
void Graph::PPRConstruction(vertex nodeID, vector<vector<pair<vertex,int>>> &Neighbor, vector<unordered_map<vertex,int>>& Label,  vector<unordered_map<vertex,unordered_set<vertex>>>& PruningPointSet, vector<unordered_map<vertex,vertex>>& PruningPointSet2){
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
        ShortestDisQuery2(topNodeID, nodeID, SupNode,TempDis, Label);//DisQueryPeak
//        ShortestDisQuery1(topNodeID, nodeID, SupNode,TempDis);

        if(TempDis<=topNodeDis){//if dispeak <= disvally
            for(int k=0;k<SupNode.size();k++){
                int supn=SupNode[k];
                unordered_map<int,vector<int>> map1;

                vSm[nodeID]->wait();
                PruningPointSet[nodeID][supn].insert(topNodeID);
                PruningPointSet2[nodeID][topNodeID] = supn;
                vSm[nodeID]->notify();

                vSm[topNodeID]->wait();
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

//function of pruning Dijkstra of batch PLL, 2023-06-07, optimal correct version 1
void Graph::BatchPCLDijk(vertex nodeID, unordered_set<vertex>& setNodes, vertex hID, vector<vector<pair<vertex,int>>> &Neighbor, vector<unordered_map<vertex,int>>& Label){
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
                ShortestDisQuery1(topNodeID, nodeID, SupNode,TempDis,Label);
                vSm[nodeID]->notify();
                if(topNodeID != nodeID){
                    vSm[topNodeID]->notify();
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
                ShortestDisQuery1(nodeID, topNodeID,SupNode,TempDis,Label);
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
void Graph::BatchPCLDijk2(vector<vertex>& p, unordered_set<vertex>& setNodes, vertex hID, vector<vector<pair<vertex,int>>> &Neighbor, vector<unordered_map<vertex,int>>& Label){
    for(auto it=p.begin();it!=p.end();++it){
        int nodeID = *it;
        BatchPCLDijk(nodeID, setNodes, hID, Neighbor, Label);
    }
}

//query by current labels
int Graph::ShortestDisQuery1(int ID1,int ID2,vector<int>& SupNode, int& d, vector<unordered_map<vertex,int>>& Label){
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
int Graph::ShortestDisQuery2(int ID1,int ID2,vector<int>& SupNode, int& d, vector<unordered_map<vertex,int>>& Label){
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

//vector-based implementation
void Graph::PLLConstructV(vector<vector<pair<vertex,int>>>& Neighbor){
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
        if(PartiTag[ID].second){
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
//    PPRV.postProcess(PruningPointSet,PruningPointSet2,vertices,threadnum);
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

//function for Query processing, new
int Graph::QueryPPLL(int ID1, int ID2){
    int dis=INF;

    if(PartiTag[ID1].first==PartiTag[ID2].first){//if same partition
        if(PSPStrategy==PreBoundary){
//            cout<<"Same-parti"<<endl;
            dis= QuerySamePartiPLL(ID1,ID2,Labels);
        }
        else if(PSPStrategy==PostBoundary){
            dis= QuerySamePartiPLL(ID1,ID2,LabelsPost);
        }
        else if(PSPStrategy==NoBoundary){//no-boundary
            if(PartiTag[ID1].second && PartiTag[ID2].second){//Case 1: both in overlay graph
//                cout<<"Same-parti: Core-Core"<<endl;
                dis=QueryOverlayPLL(ID1, ID2);
            }
            else if(PartiTag[ID1].second && !PartiTag[ID2].second){//Case 2: ID2 in partition, ID1 in core
//                cout<<"Same-parti: Core-Parti"<<endl;
                dis=QueryPartiCorePLL(ID2, ID1);
            }else if(!PartiTag[ID1].second && PartiTag[ID2].second){//Case 2: ID1 in partition, ID2 in core
//                cout<<"Same-parti: Parti-Core"<<endl;
                dis = QueryPartiCorePLL(ID1, ID2);
            }else{//Case 3: Same partition
//                cout<<"Same-parti: Non-boundary"<<endl;
                dis= QuerySamePartiPLLNo(ID1,ID2);
            }
        }
    }
    else{//if different partition
        if(PartiTag[ID1].second && PartiTag[ID2].second){//Case 1: both in overlay graph
//            cout<<"Different Parti: Core-Core"<<endl;
            dis=QueryOverlayPLL(ID1, ID2);
        }
        else if(PartiTag[ID1].second && !PartiTag[ID2].second){//Case 2: ID2 in partition, ID1 in core
//            cout<<"Different Parti: Core-Parti"<<endl;
            dis=QueryPartiCorePLL(ID2, ID1);
        }else if(!PartiTag[ID1].second && PartiTag[ID2].second){//Case 2: ID1 in partition, ID2 in core
//            cout<<"Different Parti: Parti-Core"<<endl;
            dis = QueryPartiCorePLL(ID1, ID2);
        }
        else{//Case 4: Different Partitions
//            cout<<"Different Parti: Parti-Parti"<<endl;
            dis = QueryPartiPartiPLL(ID1, ID2);
        }
    }

    return dis;
}

int Graph::QueryOverlayPLL(int ID1, int ID2){
    if(!PartiTag[ID1].second || !PartiTag[ID2].second){
        cout<<"Not overlay vertex. "<<ID1<<"("<<PartiTag[ID1].second<<") "<<ID2<<"("<<PartiTag[ID2].second<<")"<<endl; exit(1);
    }

    int hub=-1, dis=INF, dis1, dis2;
    bool ifFind=false;
    int lid,hid;
    if(NodeOrder[ID1]>NodeOrder[ID2]){
        hid=ID1, lid=ID2;
    }else{
        hid=ID2, lid=ID1;
    }
    if(Label[lid].find(hid)!=Label[lid].end()){//if found
        dis=Label[lid][hid];
        ifFind=true;
    }else{
        for(auto it=Label[ID1].begin();it!=Label[ID1].end();++it){
            hub=it->first; dis1=it->second;
            if(Label[ID2].find(hub)!=Label[ID2].end()){//if found
                dis2=Label[ID2][hub];
                if(dis>dis1+dis2){
                    dis=dis1+dis2;
                }
                ifFind=true;
            }
        }
    }

    if(!ifFind){
        cout<<"No common hub! "<<ID1<<" "<<ID2<<endl; exit(1);
    }
    return dis;
}

int Graph::QuerySamePartiPLL(int ID1, int ID2, vector<unordered_map<vertex,int>>& Labels){
    if(PartiTag[ID1].first != PartiTag[ID2].first){
        cout<<"Not same partition. "<<ID1<<"("<<PartiTag[ID1].first<<") "<<ID2<<"("<<PartiTag[ID2].first<<")"<<endl; exit(1);
    }
    int pid=PartiTag[ID1].first;
    int hub=-1, dis=INF, dis1, dis2;
    bool ifFind=false;
    int lid,hid;
    if(NodeOrder[ID1]>NodeOrder[ID2]){
        hid=ID1, lid=ID2;
    }else{
        hid=ID2, lid=ID1;
    }
    if(Labels[lid].find(hid)!=Labels[lid].end()){//if found
        dis=Labels[lid][hid];
        ifFind=true;
    }else{
        for(auto it=Labels[hid].begin();it!=Labels[hid].end();++it){
            hub=it->first; dis1=it->second;
            if(Labels[lid].find(hub)!=Labels[lid].end()){//if found
                dis2=Labels[lid][hub];
                if(dis>dis1+dis2){
                    dis=dis1+dis2;
                }
                ifFind=true;
            }
        }
    }

//    if(!ifFind){
//        cout<<"No common hub for "<<ID1<<" "<<ID2<<"; pid: "<<pid<<endl; exit(1);
//    }
    return dis;
}

//Case 4: Same tree, for original version
int Graph::QuerySamePartiPLLNo(int ID1, int ID2){//both are within the partition
    int d=INF;

    int pid1=PartiTag[ID1].first;
    int pid2=PartiTag[ID2].first;
    if(pid1==pid2){//if in the same partition
//        cout<<"Same-Parti"<<endl;
        int temp_dis = QuerySamePartiPLL(ID1,ID2,Labels);/// d2 may be wrong sometimes
        if(temp_dis<d)//QueryH2HPartition(ID1,ID2,pid1)
            d=temp_dis;//QueryH2HPartition(ID1,ID2,pid1);
        vector<int> B=BoundVertex[pid1];
        map<int,int> m1,m2;
        m1.clear();
        m2.clear();
        vector<int> B1,B2;
        B1.clear();
        B2.clear();
        int bID,d1,d2;
        for(int i=0;i<B.size();i++){
            bID=B[i];
            d1=QuerySamePartiPLL(ID1,bID,Labels);
            d2=QuerySamePartiPLL(ID2,bID,Labels);

            if(d1<d){
                B1.push_back(bID);
                m1.insert(make_pair(bID,d1));
            }
            if(d2<d){
                B2.push_back(bID);
                m2.insert(make_pair(bID,d2));
            }
        }

        int bID1, bID2, tempdis;
        if(!B1.empty() && !B2.empty()){
            for(int k=0;k<B1.size();k++){
                bID1=B1[k];
                if(m1[bID1]>d)
                    continue;
                for(int z=0;z<B2.size();z++){
                    bID2=B2[z];
                    if(m2[bID2]>d)
                        continue;
                    tempdis=m1[bID1]+QueryOverlayPLL(bID1,bID2)+m2[bID2];
                    if(tempdis<d)
                        d=tempdis;
                }
            }
        }

    }else{//if in different partitions
        cout<<"Wrong for same partition query!"<<endl;
        exit(1);
    }

    return d;
}

//Case 2: one core, one tree
int Graph::QueryPartiCorePLL(int ID1, int ID2){//ID1 partition, ID2 core
    int d=INF;

    int pid=PartiTag[ID1].first;
    int bid;
    int dis1,dis2;
    if(PSPStrategy==NoBoundary || PSPStrategy==PreBoundary){
        for(auto it=BoundVertex[pid].begin();it!=BoundVertex[pid].end();++it){
            bid=*it;
            dis1= QuerySamePartiPLL(ID1,bid,Labels);
            dis2= QueryOverlayPLL(bid,ID2);
            if(d>dis1+dis2)
                d=dis1+dis2;
        }
    }
    else if(PSPStrategy==PostBoundary){
        for(auto it=BoundVertex[pid].begin();it!=BoundVertex[pid].end();++it){
            bid=*it;
            dis1= QuerySamePartiPLL(ID1,bid,LabelsPost);
            dis2= QueryOverlayPLL(bid,ID2);
            if(d>dis1+dis2)
                d=dis1+dis2;
        }
    }

    return d;
}

//Case 3: Different trees
int Graph::QueryPartiPartiPLL(int ID1, int ID2){//both are within the partition
    int d=INF;

    int pid1=PartiTag[ID1].first;
    int pid2=PartiTag[ID2].first;
    if(pid1==pid2){//if in the same partition
        cout<<"Wrong for partition-partition query!"<<endl;
        exit(1);

    }else{//if in different partitions
//        cout<<"Parti-Parti: "<<pid1<<" "<<pid2<<endl;
        vector<int> B1=BoundVertex[pid1];
        vector<int> B2=BoundVertex[pid2];

        map<int,int> m1,m2;
        m1.clear();
        m2.clear();
        int bID1, bID2, tempdis;
        int b1,b2,d1,d2;

        if(PSPStrategy==NoBoundary || PSPStrategy==PreBoundary){
            for(int i=0;i<B1.size();i++){
                bID1=B1[i];
                m1.insert(make_pair(bID1, QuerySamePartiPLL(ID1,bID1,Labels)));
            }
            for(int j=0;j<B2.size();j++){
                bID2=B2[j];
                m2.insert(make_pair(bID2,QuerySamePartiPLL(ID2,bID2,Labels)));
            }
        }else if(PSPStrategy==PostBoundary){
            for(int i=0;i<B1.size();i++){
                bID1=B1[i];
                m1.insert(make_pair(bID1, QuerySamePartiPLL(ID1,bID1,LabelsPost)));
            }
            for(int j=0;j<B2.size();j++){
                bID2=B2[j];
                m2.insert(make_pair(bID2,QuerySamePartiPLL(ID2,bID2,LabelsPost)));
            }
        }


        for(int k=0;k<B1.size();k++){
            bID1=B1[k];

            if(m1[bID1]>d)
                continue;

            for(int z=0;z<B2.size();z++){
                bID2=B2[z];

                if(m2[bID2]>d)
                    continue;

                tempdis=m1[bID1]+QueryOverlayPLL(bID1,bID2)+m2[bID2];
                if(tempdis<d){
                    d=tempdis;
                    d1=m1[bID1]; d2=m2[bID2];
                    b1=bID1; b2=bID2;
                }

            }
        }

//        cout<<"b1, b2, d1, d2: "<<b1<<" "<<b2<<" "<<d1<<" "<<d2<<endl;
    }

    return d;
}

/// index maintenance
//function for throughput test of decrease update
void Graph::PPLLBatchUpdateDec(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1) {
    Timer tt;
    Timer tt2;
    tt.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }

    vUpdated.assign(node_num, false);
    if(!partiBatch.empty()){
        if(partiBatch.size()>threadnum){
            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
        }
//        cout<<"Update Partition number: "<<partiBatch.size()<<endl;
        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::DecreasePartiBatchUpdatePLL, this, pid, boost::ref(it->second), boost::ref(overlayBatch), true));
        }
        thread.join_all();
    }

    if(!overlayBatch.empty()){
//        cout<<"OverlayBatch size: "<<overlayBatch.size()<<endl;
        int a,b,oldW,newW;
        for(int i=0;i<overlayBatch.size();++i){
            a=overlayBatch[i].first.first; b=overlayBatch[i].first.second;
            oldW=overlayBatch[i].second.first; newW=overlayBatch[i].second.second;
            if(NeighborsOverlay[a].find(b)!=NeighborsOverlay[a].end()){//if found
                NeighborsOverlay[a][b]=newW;
            }else{
                cout<<"Wrong. Not found "<<a<<" "<<b<<endl; exit(1);
            }
            if(NeighborsOverlay[b].find(a)!=NeighborsOverlay[b].end()){//if found
                NeighborsOverlay[b][a]=newW;
            }else{
                cout<<"Wrong. Not found "<<b<<" "<<a<<endl; exit(1);
            }
            DecreasePSL(a,b,oldW,newW,NeighborsOverlayV,Label);
        }
    }

    // repair the partition index
    if(PSPStrategy>=PostBoundary){
        tt2.start();
        Repair_PartiIndexPLL(true, false, partiBatch);//post
        tt2.stop();

    }

    tt.stop();
}

//function for throughput test of decrease update
void Graph::PPLLBatchUpdateDecPre(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1) {
    Timer tt;
    Timer tt2;
    tt.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }

    AllPairBoundaryDisUpdate(true, false, partiBatch, overlayBatch);

    vUpdated.assign(node_num, false);
    if(!partiBatch.empty()){
        if(partiBatch.size()>threadnum){
            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
        }
//        cout<<"Update Partition number: "<<partiBatch.size()<<endl;
        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::DecreasePartiBatchUpdatePLL, this, pid, boost::ref(it->second), boost::ref(overlayBatch), false));
        }
        thread.join_all();
    }

    if(!overlayBatch.empty()){
        cout<<"OverlayBatch size: "<<overlayBatch.size()<<endl;
        int a,b,oldW,newW;
        for(int i=0;i<overlayBatch.size();++i){
            a=overlayBatch[i].first.first; b=overlayBatch[i].first.second;
            oldW=overlayBatch[i].second.first; newW=overlayBatch[i].second.second;
            if(NeighborsOverlay[a].find(b)!=NeighborsOverlay[a].end()){//if found
                NeighborsOverlay[a][b]=newW;
            }else{
                cout<<"Wrong. Not found "<<a<<" "<<b<<endl; exit(1);
            }
            if(NeighborsOverlay[b].find(a)!=NeighborsOverlay[b].end()){//if found
                NeighborsOverlay[b][a]=newW;
            }else{
                cout<<"Wrong. Not found "<<b<<" "<<a<<endl; exit(1);
            }
            DecreasePSL(a,b,oldW,newW,NeighborsOverlayV,Label);
        }
    }

    tt.stop();
}

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

//partition update PCH
void Graph::DecreasePartiBatchUpdatePLL(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay, bool ifOverlay){
    //partition batch decrease update
//    vector<pair<pair<int,int>,int>> updatedSC;
    int ID1,ID2,oldW,newW;
    for(int i=0;i<wBatch.size();++i){
        ID1=wBatch[i].first.first; ID2=wBatch[i].first.second;
        oldW=wBatch[i].second.first; newW=wBatch[i].second.second;
        DecreasePSL(ID1,ID2,oldW,newW,NeighborsParti,Labels);
    }

    if(ifOverlay){
        //check the boundary edge within partition
        int bid1,bid2,olddis,newdis;
        for(int i=0;i<BoundVertex[pid].size();++i){
            bid1=BoundVertex[pid][i];
            for(int j=i+1;j<BoundVertex[pid].size();++j){
                bid2=BoundVertex[pid][j];
                newdis= QuerySamePartiPLL(bid1,bid2,Labels);
                if(newdis>=INF){
                    continue;
                }
                if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
                    olddis=NeighborsOverlay[bid1][bid2];
                }else{//if not found
                    cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
                }

                if(newdis<olddis){
//            cout<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
//                NeighborsOverlay[bid1][bid2]=newdis;
//                NeighborsOverlay[bid2][bid1]=newdis;
                    sm->wait();
                    weightOverlay.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));//weightOverlay collect the changed edges on overlay graph
                    sm->notify();
                }else if(newdis>olddis){
                    cout<<"Wrong boundary shortcut value! "<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl; exit(1);
                }
            }
        }
    }

}

void Graph::IncreasePartiBatchUpdatePLL(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay, bool ifOverlay){
    //partition batch decrease update
//    vector<pair<pair<int,int>,int>> updatedSC;
    int ID1,ID2,oldW,newW;
    for(int i=0;i<wBatch.size();++i){
        ID1=wBatch[i].first.first; ID2=wBatch[i].first.second;
        oldW=wBatch[i].second.first; newW=wBatch[i].second.second;
        IncreasePSL(ID1,ID2,oldW,newW,NeighborsParti,Labels,PruningPointSetP,PruningPointSetP2);
    }

    if(ifOverlay){
        //check the boundary edge within partition
        int bid1,bid2,olddis,newdis;
        for(int i=0;i<BoundVertex[pid].size();++i){
            bid1=BoundVertex[pid][i];
            for(int j=i+1;j<BoundVertex[pid].size();++j){
                bid2=BoundVertex[pid][j];
                newdis= QuerySamePartiPLL(bid1,bid2,Labels);
                if(newdis>=INF){
                    continue;
                }
                if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
                    olddis=NeighborsOverlay[bid1][bid2];
                }else{//if not found
                    cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
                }

                if(newdis>olddis){
//            cout<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
//                NeighborsOverlay[bid1][bid2]=newdis;
//                NeighborsOverlay[bid2][bid1]=newdis;
                    sm->wait();
                    weightOverlay.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));//weightOverlay collect the changed edges on overlay graph
                    sm->notify();
                }else if(newdis<olddis){
                    cout<<"Wrong boundary shortcut value! "<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl; exit(1);
                }
            }
        }
    }

}

//Function of repair the partition index
void Graph::Repair_PartiIndexPLL(bool ifParallel, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch){
//    repairShortcuts.assign(node_num, unordered_map<vertex,pair<int,int>>());
    ifRepaired.assign(partiNum, false);
    if(ifParallel){
        // multi-thread
//        cout<<"Multi-thread computation!"<<endl;
        //multi-thread
        if(threadnum<partiNum){
            vector<vector<int>> processID;
            processID.assign(threadnum, vector<int>());
            vector<int> vertices;
            for(int pid=0;pid<partiNum;++pid){
                vertices.emplace_back(pid);
            }
            ThreadDistribute(vertices, processID);
//            cout<<"Batch number: "<<processID[0].size()<<endl;
            boost::thread_group thread;
            for(auto j=0;j<processID.size();++j){
                thread.add_thread(new boost::thread(&Graph::RepairPartitionIndexPLLV, this, boost::ref(processID[j]), ifIncrease, boost::ref(partiBatch)));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            if(ifIncrease){//increase update
                boost::thread_group thread;
                for(auto j=0;j<partiNum;++j){
                    thread.add_thread(new boost::thread(&Graph::RepairPartitionIndexIncreasePLL, this, j, boost::ref(partiBatch)));
                }
                thread.join_all();
            }
            else{//decrease update
                boost::thread_group thread;
                for(auto j=0;j<partiNum;++j){
                    thread.add_thread(new boost::thread(&Graph::RepairPartitionIndexDecreasePLL, this, j, boost::ref(partiBatch)));
                }
                thread.join_all();
            }

        }
    }
    else{
        // single thread
        if(ifIncrease){
            for(int k=0;k<partiNum;k++){
//                cout<<"Repairing partition "<<k<<endl;
                RepairPartitionIndexIncreasePLL(k,partiBatch);
            }
        }
        else{
            for(int k=0;k<partiNum;k++){
//                cout<<"Repairing partition "<<k<<endl;
                RepairPartitionIndexDecreasePLL(k,partiBatch);
            }
        }

    }

//    int pNum=0;
//    for(auto it=ifRepaired.begin();it!=ifRepaired.end();++it){
//        if(*it){
//            ++pNum;
//        }
//    }
//    cout<<"Repaired partition number: "<<pNum<<endl;

}

void Graph::RepairPartitionIndexDecreasePLL(int pid, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch) {
    int ID1,ID2,woverlay;
    int wlocal=INF;
    //boundary edges
    vector<pair<pair<int,int>,pair<int,int>>> weightsParti;//collect the changed edges on overlay graph
    weightsParti.clear();

    if(partiBatch.find(pid)!=partiBatch.end()){//if found
        weightsParti=partiBatch[pid];
    }

    for(int i=0;i<BoundVertex[pid].size();i++){
        ID1=BoundVertex[pid][i];
        for(int j=i+1;j<BoundVertex[pid].size();j++){
            ID2=BoundVertex[pid][j];

            if(NeighborsPartiPost[ID1].find(ID2)!=NeighborsPartiPost[ID1].end()){//if found
                wlocal=NeighborsPartiPost[ID1][ID2];
            }else{
                cout<<"Post. Not found edge e("<<ID1<<","<<ID2<<") in overlay graph!"<<endl; exit(1);
            }
            woverlay=QueryOverlayPLL(ID1,ID2);
            bool found=false;//whether the boundary edge exist or not
            int wei;
            if(woverlay<wlocal){
                if(ID1<ID2){
                    weightsParti.emplace_back(make_pair(ID1,ID2),make_pair(wlocal,woverlay));
                }else{
                    weightsParti.emplace_back(make_pair(ID2,ID1),make_pair(wlocal,woverlay));
                }

//                DecreaseParti(ID1,ID2, woverlay, NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid]);
            }else if(woverlay>wlocal){
                cout<<"Something wrong: shortest path in the overlay graph rather than in the subgraph. "<<ID1<<" "<<ID2<<" "<< wlocal<<" "<<woverlay<< endl; exit(1);
            }

        }
    }

    if(!weightsParti.empty()){
//        cout<<"Size of weightsParti of partition "<<pid<<" : "<<weightsParti.size()<<endl;
        int a,b,oldW, newW;
        for(int i=0;i<weightsParti.size();++i){
            a=weightsParti[i].first.first, b=weightsParti[i].first.second, oldW=weightsParti[i].second.first, newW=weightsParti[i].second.second;
            if(NeighborsPartiPost[a].find(b)!=NeighborsPartiPost[a].end()){//if found
                NeighborsPartiPost[a][b]=newW;
            }else{
                cout<<"Wrong. Not found "<<a<<" "<<b<<endl; exit(1);
            }
            if(NeighborsPartiPost[b].find(a)!=NeighborsPartiPost[b].end()){//if found
                NeighborsPartiPost[b][a]=newW;
            }else{
                cout<<"Wrong. Not found "<<b<<" "<<a<<endl; exit(1);
            }
            DecreasePSL(a,b,oldW,newW, NeighborsPartiPostV, LabelsPost);
        }

        ifRepaired[pid]=true;
    }

}

void Graph::RepairPartitionIndexPLLV(vector<int>& p, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch) {

    if(ifIncrease){
        for(int i=0;i<p.size();++i){
            RepairPartitionIndexIncreasePLL(p[i], partiBatch);
        }
    }
    else{
        for(int i=0;i<p.size();++i){
            RepairPartitionIndexDecreasePLL(p[i], partiBatch);
        }
    }

}

void Graph::RepairPartitionIndexIncreasePLL(int pid, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch) {
    int ID1,ID2,woverlay;
    int wlocal=INF;
    //boundary edges
    vector<pair<pair<int,int>,pair<int,int>>> weightsParti;//collect the changed edges on overlay graph
    weightsParti.clear();
    map<pair<int,int>,pair<int,int>> updatesSet;

    if(partiBatch.find(pid)!=partiBatch.end()){//if found
//        cout<<"Find for "<<pid<<" "<<partiBatch[pid].size()<<endl;
//        weightsParti=partiBatch[pid];
        for(auto it=partiBatch[pid].begin();it!=partiBatch[pid].end();++it){
            ID1=it->first.first, ID2=it->first.second;
            if(ID1<ID2){
                if(updatesSet.find(make_pair(ID1,ID2))==updatesSet.end()){//if not found
                    updatesSet.insert(*it);
                }else{
                    cout<<"Original edge. Already exists! "<<ID1<<" "<<ID2<<" "<<updatesSet[make_pair(ID1,ID2)].first<<" "<<updatesSet[make_pair(ID1,ID2)].second<<endl; exit(1);
                }
//                    weightsParti.emplace_back(make_pair(ID1,ID2), make_pair(wlocal, woverlay));//oldW,newW
            }else{
                if(updatesSet.find(make_pair(ID2,ID1))==updatesSet.end()){//if not found
                    updatesSet.insert(*it);
                }else{
                    cout<<"Original edge. Already exists! "<<ID2<<" "<<ID1<<" "<<updatesSet[make_pair(ID2,ID1)].first<<" "<<updatesSet[make_pair(ID2,ID1)].second<<endl; exit(1);
                }
//                    weightsParti.emplace_back(make_pair(ID2,ID1), make_pair(wlocal, woverlay));//oldW,newW
            }
        }
    }

    for(int i=0;i<BoundVertex[pid].size();i++){
        ID1=BoundVertex[pid][i];
        for(int j=i+1;j<BoundVertex[pid].size();j++){
            ID2=BoundVertex[pid][j];
            if(NeighborsPartiPost[ID1].find(ID2)!=NeighborsPartiPost[ID1].end()){//if found
                wlocal=NeighborsPartiPost[ID1][ID2];
            }else{
                cout<<"Post. Not found edge e("<<ID1<<","<<ID2<<") in overlay graph!"<<endl; exit(1);
            }
            woverlay=QueryOverlayPLL(ID1,ID2);
            bool found=false;//whether the boundary edge exist or not
            int wei;
//            if(pid==49){
//                cout<<pid<<": "<<ID1<<" "<<ID2<<" "<<wlocal<<" "<<woverlay<<endl;
//            }

            if(woverlay>wlocal){
                // update partition index
                if(ID1<ID2){
                    if(updatesSet.find(make_pair(ID1,ID2))==updatesSet.end()){//if not found
                        updatesSet.insert({make_pair(ID1,ID2),make_pair(wlocal, woverlay)});
                    }else{
                        cout<<"All-pair check. Already exists! "<<ID1<<" "<<ID2<<" "<<updatesSet[make_pair(ID1,ID2)].first<<" "<<updatesSet[make_pair(ID1,ID2)].second<<endl;
                        updatesSet[make_pair(ID1,ID2)].first=wlocal,updatesSet[make_pair(ID1,ID2)].second=woverlay;
//                        exit(1);
                    }
//                    weightsParti.emplace_back(make_pair(ID1,ID2), make_pair(wlocal, woverlay));//oldW,newW
                }else{
                    if(updatesSet.find(make_pair(ID2,ID1))==updatesSet.end()){//if not found
                        updatesSet.insert({make_pair(ID2,ID1),make_pair(wlocal, woverlay)});
                    }else{
                        cout<<"All-pair check. Already exists! "<<ID2<<" "<<ID1<<" "<<updatesSet[make_pair(ID2,ID1)].first<<" "<<updatesSet[make_pair(ID2,ID1)].second<<endl;
                        updatesSet[make_pair(ID2,ID1)].first=wlocal,updatesSet[make_pair(ID1,ID2)].second=woverlay;
//                        exit(1);
                    }
//                    weightsParti.emplace_back(make_pair(ID2,ID1), make_pair(wlocal, woverlay));//oldW,newW
                }

//                DecreaseParti(ID1,ID2, woverlay, NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid]);
            }else if(woverlay<wlocal){
                cout<<"Something wrong: shortest path in the overlay graph rather than in the subgraph. "<<ID1<<" "<<ID2<<" "<< wlocal<<" "<<woverlay<< endl; exit(1);
            }

        }
    }

    for(auto it=updatesSet.begin();it!=updatesSet.end();++it){
        weightsParti.emplace_back(*it);
    }

    if(!weightsParti.empty()){
//        cout<<"Size of weightsParti of partition "<<pid<<" : "<<weightsParti.size()<<endl;
        int a,b,oldW, newW;
        for(int i=0;i<weightsParti.size();++i) {
            a = weightsParti[i].first.first, b = weightsParti[i].first.second, oldW = weightsParti[i].second.first, newW = weightsParti[i].second.second;
            if (NeighborsPartiPost[a].find(b) != NeighborsPartiPost[a].end()) {//if found
                NeighborsPartiPost[a][b] = newW;
            } else {
                cout << "Wrong. Not found " << a << " " << b << endl;
                exit(1);
            }
            if (NeighborsPartiPost[b].find(a) != NeighborsPartiPost[b].end()) {//if found
                NeighborsPartiPost[b][a] = newW;
            } else {
                cout << "Wrong. Not found " << b << " " << a << endl;
                exit(1);
            }
            IncreasePSL(a,b,oldW,newW, NeighborsPartiPostV, LabelsPost, PruningPointSetPost, PruningPointSetPost2);
        }
        ifRepaired[pid]=true;
    }

}

//function for throughput test of increase update
void Graph::PPLLBatchUpdateInc(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1) {
    Timer tt;
    Timer tt2;
    tt.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }

    vUpdated.assign(node_num, false);
    if(!partiBatch.empty()){
        if(partiBatch.size()>threadnum){
            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
        }
//        cout<<"Update Partition number: "<<partiBatch.size()<<endl;
        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::IncreasePartiBatchUpdatePLL, this, pid, boost::ref(it->second), boost::ref(overlayBatch), true));
        }
        thread.join_all();
    }

    if(!overlayBatch.empty()){
//        cout<<"OverlayBatch size: "<<overlayBatch.size()<<endl;
        int a,b,oldW,newW;
        for(int i=0;i<overlayBatch.size();++i){
            a=overlayBatch[i].first.first; b=overlayBatch[i].first.second; oldW=overlayBatch[i].second.first; newW=overlayBatch[i].second.second;
            if(NeighborsOverlay[a].find(b)!=NeighborsOverlay[a].end()){//if found
                NeighborsOverlay[a][b]=newW;
            }else{
                cout<<"Wrong. Not found "<<a<<" "<<b<<endl; exit(1);
            }
            if(NeighborsOverlay[b].find(a)!=NeighborsOverlay[b].end()){//if found
                NeighborsOverlay[b][a]=newW;
            }else{
                cout<<"Wrong. Not found "<<b<<" "<<a<<endl; exit(1);
            }
            IncreasePSL(a,b,oldW,newW,NeighborsOverlayV,Label,PruningPointSet,PruningPointSet2);
        }
    }


    // repair the partition index
    if(PSPStrategy>=PostBoundary){
        Repair_PartiIndexPLL(true, true, partiBatch);//post
    }

    tt.stop();
//    CorrectnessCheck(100);
}

//function for throughput test of increase update
void Graph::PPLLBatchUpdateIncPre(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1) {
    Timer tt;
    Timer tt2;
    tt.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }

    vUpdated.assign(node_num, false);
    AllPairBoundaryDisUpdate(true, true, partiBatch, overlayBatch);

    vUpdated.assign(node_num, false);
    if(!partiBatch.empty()){
        if(partiBatch.size()>threadnum){
            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
        }
//        cout<<"Update Partition number: "<<partiBatch.size()<<endl;
        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::IncreasePartiBatchUpdatePLL, this, pid, boost::ref(it->second), boost::ref(overlayBatch), false));
        }
        thread.join_all();
    }

    if(!overlayBatch.empty()){
//        cout<<"OverlayBatch size: "<<overlayBatch.size()<<endl;
        int a,b,oldW,newW;
        for(int i=0;i<overlayBatch.size();++i){
            a=overlayBatch[i].first.first; b=overlayBatch[i].first.second;
            oldW=overlayBatch[i].second.first; newW=overlayBatch[i].second.second;
            if(NeighborsOverlay[a].find(b)!=NeighborsOverlay[a].end()){//if found
                NeighborsOverlay[a][b]=newW;
            }else{
                cout<<"Wrong. Not found "<<a<<" "<<b<<endl; exit(1);
            }
            if(NeighborsOverlay[b].find(a)!=NeighborsOverlay[b].end()){//if found
                NeighborsOverlay[b][a]=newW;
            }else{
                cout<<"Wrong. Not found "<<b<<" "<<a<<endl; exit(1);
            }
            IncreasePSL(a,b,oldW,newW,NeighborsOverlayV,Label,PruningPointSet,PruningPointSet2);
        }
    }
    tt.stop();
//    CorrectnessCheck(100);
}

//new version: set version with NoSupportedPair, clean label, 2023-05-14, iterative update, queue version, correct
void Graph::IncreasePSL(int a, int b, int oldW, int newW, vector<vector<pair<vertex,int>>> &Neighbor,vector<unordered_map<vertex,int>> &Label, vector<unordered_map<vertex,unordered_set<vertex>>> &PruningPointNew, vector<unordered_map<vertex,vertex>>& PruningPointSet2){//
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
        vector<unordered_set<int>> ChangedLabels;
        set<pair<int,int>> NoSupportedPair;
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
            CoarseUpdate(LID, HID, oldW, WaitPro, WaitProP, AL1, AL2, AL2Check, Neighbor, Label, ifDebug, lid, hid, ChangedLabels);
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
            RefineUpdate( WaitPro, WaitProP, AL1, AL2, AL2Check, outdatedPruning, newPruningPoints, PruningPointSet2, Neighbor, Label, PruningPointNew, ifDebug, lid, hid, ChangedLabels, NoSupportedPair);
//            cout<<"After AL2: "<<AL1.size()<<" "<<AL2Check.size()<<" "<<WaitPro.size()<<" "<<WaitProP.size()<<endl;
        }
//        cout<<"Done."<<endl;

//        cout<<"Remove and add pruning point."<<endl;
        PPRClean(Neighbor, newPruningPoints, outdatedPruning, ifDebug, lid, hid, Label, PruningPointNew, PruningPointSet2, ChangedLabels);

    }
    else{
        if(ifDebug){
            cout<<"Not triggered! "<<LID<<" "<<HID<<": "<<dislower<<" "<<dispeak<<" "<<oldW<<" "<<newW<<endl;
        }
    }

}
//function of propagating the AL1 labels, output the updated labels for AL2 check, queue version, 2023-03-26, new
void Graph::CoarseUpdate(int LID,int HID, int oldW, queue<pair<int,pair<int,int>>>& WaitPro, queue<pair<int,pair<int,int>>>& WaitProP, vector<pair<int,int>>& AL1, vector<pair<int,int>>& AL2, vector<pair<int,int>>& AL2Check, vector<vector<pair<vertex,int>>> &Neighbor,vector<unordered_map<vertex,int>> &Label, bool ifDebug, int lid, int hid, vector<unordered_set<int>>& ChangedLabels){
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
void Graph::RefineUpdate(queue<pair<int,pair<int,int>>>& WaitPro, queue<pair<int,pair<int,int>>>& WaitProP, vector<pair<int,int>>& AL1, vector<pair<int,int>>& AL2, vector<pair<int,int>>& AL2Check, set<tuple<vertex,vertex,vertex>>& outdatedPruning, map<pair<vertex,vertex>,vertex>& newPruningPoints, vector<unordered_map<vertex,vertex>>& PruningPointSet2, vector<vector<pair<vertex,int>>> &Neighbor,vector<unordered_map<vertex,int>> &Label, vector<unordered_map<vertex,unordered_set<vertex>>> &PruningPointNew, bool ifDebug, int lid, int hid, vector<unordered_set<int>>& ChangedLabels, set<pair<int,int>>& NoSupportedPair){
    int disvally,dispeak,peakhub,vallyID;
    pair<int,int> peakPair;//distance, hubID
    pair<int,int> vallyPair;//distance, vallyID

    int curID,hubID,hubDis;

    bool ifUpdate = false;

    //check and get the source AL2, AL1->AL2
    for(int i=0;i<AL2Check.size();++i){//for each increased label
        curID = AL2Check[i].first, hubID = AL2Check[i].second;

        bool temp_bool = PPRCheck(curID,hubID,Neighbor,Label,PruningPointNew,WaitProP,AL2,newPruningPoints, PruningPointSet2, outdatedPruning,ifDebug,lid,hid, ChangedLabels, NoSupportedPair);///
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

        }
        WaitProP=WaitProPTem;
    }

}
//function of checking the pruned label by PPR, pair version, queue version, correct
bool Graph::PPRCheck(int curID, int hubID, vector<vector<pair<vertex,int>>> &Neighbor,vector<unordered_map<vertex,int>> &Label, vector<unordered_map<vertex,unordered_set<vertex>>> &PruningPointNew, queue<pair<int,pair<int,int>>>& WaitProPTem, vector<pair<int,int>> &AL2, map<pair<vertex,vertex>,vertex> &newPruningPoints, vector<unordered_map<vertex,vertex>>& PruningPointSet2, set<tuple<vertex,vertex,vertex>> &outdatedPruning, bool ifDebug, int lid, int hid, vector<unordered_set<int>>& ChangedLabels, set<pair<int,int>>& NoSupportedPair){
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
void Graph::PPRClean(vector<vector<pair<vertex,int>>> &Neighbor, map<pair<vertex,vertex>,vertex> &newPruningPoints, set<tuple<vertex,vertex,vertex>> &outdatedPruning, bool ifDebug, int lid, int hid, vector<unordered_map<vertex,int>> &Label, vector<unordered_map<vertex,unordered_set<vertex>>>& PruningPointSet, vector<unordered_map<vertex,vertex>>& PruningPointSet2, vector<unordered_set<int>>& ChangedLabels){
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