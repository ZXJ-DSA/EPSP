/*
 * TreeIndex.cpp
 *
 *  Created on: 16 June 2023
 *      Author: Xinjie ZHOU
 */
#include "headPSP.h"

vector<int> NodeOrder_;//nodeID order
vector<int> _DD_;//true degree, temporal degree ,_DD2_
vector<int> _DD2_;//count

//// Index Construction
//Function of constructing tree index for partitions
void Graph::AllPairBoundaryDisCompute(bool ifParallel) {

    if(ifParallel){
        cout<<"Multiple thread computation for all-pair boundary distance computation!"<<endl;
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
                thread.add_thread(new boost::thread(&Graph::AllPairBoundaryDisComputePartiV, this, boost::ref(processID[j])));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                thread.add_thread(new boost::thread(&Graph::AllPairBoundaryDisComputeParti, this, j));
            }
            thread.join_all();
        }

    }
    else{
        cout<<"Single thread computation!"<<endl;
        //single-thread

        for(int pid=0;pid<partiNum;++pid){
//            cout<<"Partition "<<pid<<endl;
            AllPairBoundaryDisComputeParti(pid);
        }


    }

}

//Function of constructing tree index for partitions
void Graph::AllPairBoundaryDisUpdate(bool ifParallel, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>>& partiBatch, vector<pair<pair<int, int>, pair<int, int>>>& overlayBatch) {

    if(ifParallel){
//        cout<<"Multiple thread computation for all-pair boundary distance computation!"<<endl;
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
                thread.add_thread(new boost::thread(&Graph::AllPairBoundaryDisUpdatePartiV, this, boost::ref(processID[j]), ifIncrease, boost::ref(partiBatch), boost::ref(overlayBatch)));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                thread.add_thread(new boost::thread(&Graph::AllPairBoundaryDisUpdateParti, this, j, ifIncrease, boost::ref(partiBatch), boost::ref(overlayBatch)));
            }
            thread.join_all();
        }

    }
    else{
        cout<<"Single thread computation!"<<endl;
        //single-thread
        for(int pid=0;pid<partiNum;++pid){
//            cout<<"Partition "<<pid<<endl;
            AllPairBoundaryDisUpdateParti(pid,ifIncrease,partiBatch, overlayBatch);
        }
    }

}

void Graph::AllPairBoundaryDisUpdate(bool ifParallel, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>>& partiBatch, map<pair<int, int>, pair<int, int>>& overlayBatch) {

    if(ifParallel){
//        cout<<"Multiple thread computation for all-pair boundary distance computation!"<<endl;
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
                thread.add_thread(new boost::thread(&Graph::AllPairBoundaryDisUpdatePartiMapV, this, boost::ref(processID[j]), ifIncrease, boost::ref(partiBatch), boost::ref(overlayBatch)));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                thread.add_thread(new boost::thread(&Graph::AllPairBoundaryDisUpdatePartiMap, this, j, ifIncrease, boost::ref(partiBatch), boost::ref(overlayBatch)));
            }
            thread.join_all();
        }

    }
    else{
        cout<<"Single thread computation!"<<endl;
        //single-thread
        for(int pid=0;pid<partiNum;++pid){
//            cout<<"Partition "<<pid<<endl;
            AllPairBoundaryDisUpdatePartiMap(pid,ifIncrease,partiBatch, overlayBatch);
        }
    }

}

void Graph::AllPairBoundaryDisComputePartiV(vector<int>& p) {
    for(auto it=p.begin();it!=p.end();++it){
        AllPairBoundaryDisComputeParti(*it);
    }
}

void Graph::AllPairBoundaryDisUpdatePartiMapV(vector<int>& p, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>>& partiBatch, map<pair<int, int>, pair<int, int>>& overlayBatch) {
    for(auto it=p.begin();it!=p.end();++it){
        AllPairBoundaryDisUpdatePartiMap(*it, ifIncrease, partiBatch, overlayBatch);
    }
}

void Graph::AllPairBoundaryDisUpdatePartiV(vector<int>& p, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>>& partiBatch, vector<pair<pair<int, int>, pair<int, int>>>& overlayBatch) {
    for(auto it=p.begin();it!=p.end();++it){
        AllPairBoundaryDisUpdateParti(*it, ifIncrease, partiBatch, overlayBatch);
    }
}

void Graph::AllPairBoundaryDisComputeParti(int pid) {
    int bid1,bid2;
    for(int i=0;i<BoundVertex[pid].size();++i){
        bid1=BoundVertex[pid][i];
        vector<int> IDs;
        for(int j=i+1;j<BoundVertex[pid].size();++j){
            bid2=BoundVertex[pid][j];
            IDs.push_back(bid2);
        }
        vector<int> Dis;
        BoundaryDijkstra(bid1, IDs, Dis, Neighbor);
        // Insert to overlay graph
        for(int j=0;j<IDs.size();++j){
            bid2=IDs[j];
            if(NeighborsOverlay[bid1].find(bid2)!=NeighborsOverlay[bid1].end()){//if found
                if(NeighborsOverlay[bid1][bid2]>Dis[j]){
                    NeighborsOverlay[bid1][bid2]=Dis[j];
                }
                else if(NeighborsOverlay[bid1][bid2]<Dis[j]){
                    cout<<"Wrong overlay! "<<bid1<<" "<<bid2<<" "<<NeighborsOverlay[bid1][bid2]<<" "<<Dis[j]<<endl; exit(1);
                }

                for(int p=0;p<NeighborsOverlayV[bid1].size();++p){
                    if(bid2 == NeighborsOverlayV[bid1][p].first){
                        if(NeighborsOverlayV[bid1][p].second > Dis[j]){
                            NeighborsOverlayV[bid1][p].second = Dis[j];
                        }else if(NeighborsOverlayV[bid1][p].second < Dis[j]){
                            cout<<"Wrong! "<<bid1<<" "<<bid2<<" "<<NeighborsOverlayV[bid1][p].second<<" "<<Dis[j]<<endl; exit(1);
                        }
                    }
                }
            }else{
                NeighborsOverlay[bid1][bid2]=Dis[j];
                NeighborsOverlayV[bid1].push_back(make_pair(bid2,Dis[j]));
            }
            if(NeighborsOverlay[bid2].find(bid1)!=NeighborsOverlay[bid2].end()){
                if(NeighborsOverlay[bid2][bid1]>Dis[j]){
                    NeighborsOverlay[bid2][bid1]=Dis[j];
                }else if(NeighborsOverlay[bid2][bid1]<Dis[j]){
                    cout<<"Wrong overlay! "<<bid1<<" "<<bid2<<" "<<NeighborsOverlay[bid2][bid1]<<" "<<Dis[j]<<endl; exit(1);
                }

                for(int p=0;p<NeighborsOverlayV[bid2].size();++p){
                    if(bid1 == NeighborsOverlayV[bid2][p].first){
                        if(NeighborsOverlayV[bid2][p].second > Dis[j]){
                            NeighborsOverlayV[bid2][p].second = Dis[j];
                        }else if(NeighborsOverlayV[bid2][p].second < Dis[j]){
                            cout<<"Wrong! "<<bid1<<" "<<bid2<<" "<<NeighborsOverlayV[bid2][p].second<<" "<<Dis[j]<<endl; exit(1);
                        }
                    }
                }
            }else{
                NeighborsOverlay[bid2][bid1]=Dis[j];
                NeighborsOverlayV[bid2].push_back(make_pair(bid1,Dis[j]));
            }
        }
        // Insert to partition graph
        bool ifExist=false;
        for(int j=0;j<IDs.size();++j){
            bid2=IDs[j];
            ifExist=false;
            for(int p=0;p<NeighborsParti[bid1].size();++p){
                if(bid2 == NeighborsParti[bid1][p].first){
                    if(NeighborsParti[bid1][p].second > Dis[j]){
                        NeighborsParti[bid1][p].second = Dis[j];
                    }else if(NeighborsParti[bid1][p].second < Dis[j]){
                        cout<<"Wrong parti! "<<bid1<<" "<<bid2<<" "<<NeighborsParti[bid1][p].second<<" "<<Dis[j]<<endl; exit(1);
                    }
                    ifExist=true;
                }
            }
            if(!ifExist){
                NeighborsParti[bid1].emplace_back(bid2,Dis[j]);
            }
            ifExist=false;
            for(int p=0;p<NeighborsParti[bid2].size();++p){
                if(bid1 == NeighborsParti[bid2][p].first){
                    if(NeighborsParti[bid2][p].second > Dis[j]){
                        NeighborsParti[bid2][p].second = Dis[j];
                    }else if(NeighborsParti[bid2][p].second < Dis[j]){
                        cout<<"Wrong parti! "<<bid1<<" "<<bid2<<" "<<NeighborsParti[bid2][p].second<<" "<<Dis[j]<<endl; exit(1);
                    }
                    ifExist=true;
                }
            }
            if(!ifExist){
                NeighborsParti[bid2].emplace_back(bid1,Dis[j]);
            }
        }
    }
}

void Graph::AllPairBoundaryDisUpdateParti(int pid, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>>& partiBatch, vector<pair<pair<int, int>, pair<int, int>>>& overlayBatch) {
    int bid1,bid2;

    for(int i=0;i<BoundVertex[pid].size();++i){
        bid1=BoundVertex[pid][i];
        vector<int> IDs;
        for(int j=i+1;j<BoundVertex[pid].size();++j){
            bid2=BoundVertex[pid][j];
            IDs.push_back(bid2);
        }
        vector<int> Dis;
        BoundaryDijkstra(bid1, IDs, Dis, Neighbor);
        if(ifIncrease){
            // Insert to overlay graph
            for(int j=0;j<IDs.size();++j){
                bid2=IDs[j];
                if(NeighborsOverlay[bid1].find(bid2)!=NeighborsOverlay[bid1].end()){
                    if(NeighborsOverlay[bid1][bid2]<Dis[j]){
                        sm->wait();
                        overlayBatch.push_back({make_pair(bid1,bid2), make_pair(NeighborsOverlay[bid1][bid2],Dis[j])});
                        sm->notify();
//                        NeighborsOverlay[bid1][bid2]=Dis[j];
                    }
                    else if(NeighborsOverlay[bid1][bid2]>Dis[j]){
                        cout<<"Wrong overlay! "<<bid1<<" "<<bid2<<" "<<NeighborsOverlay[bid1][bid2]<<" "<<Dis[j]<<endl; exit(1);
                    }
                }else{
                    cout<<"Wrong overlay! Not found. "<<bid1<<" "<<bid2<<" "<<endl; exit(1);
                }
//                if(NeighborsOverlay[bid2].find(bid1)!=NeighborsOverlay[bid2].end()){
//                    if(NeighborsOverlay[bid2][bid1]<Dis[j]){
////                        NeighborsOverlay[bid2][bid1]=Dis[j];
//                    }else if(NeighborsOverlay[bid2][bid1]>Dis[j]){
//                        cout<<"Wrong overlay! "<<bid1<<" "<<bid2<<" "<<NeighborsOverlay[bid2][bid1]<<" "<<Dis[j]<<endl; exit(1);
//                    }
//
//                }else{
//                    cout<<"Wrong overlay! Not found. "<<bid1<<" "<<bid2<<" "<<endl; exit(1);
//                }
            }
            // Insert to partition graph
            for(int j=0;j<IDs.size();++j){
                bid2=IDs[j];
                for(int p=0;p<NeighborsParti[bid1].size();++p){
                    if(bid2 == NeighborsParti[bid1][p].first){
                        if(NeighborsParti[bid1][p].second < Dis[j]){
                            if(partiBatch.find(pid)==partiBatch.end()){//if not found
                                sm->wait();
                                partiBatch.insert({pid, vector<pair<pair<int, int>, pair<int, int>>>()});
                                partiBatch[pid].emplace_back(make_pair(bid1,bid2), make_pair(NeighborsParti[bid1][p].second,Dis[j]));
                                sm->notify();
                            }else{//if found
                                partiBatch[pid].emplace_back(make_pair(bid1,bid2), make_pair(NeighborsParti[bid1][p].second,Dis[j]));
                            }
//                            NeighborsParti[bid1][p].second = Dis[j];
                        }else if(NeighborsParti[bid1][p].second > Dis[j]){
                            cout<<"Wrong parti! "<<bid1<<" "<<bid2<<" "<<NeighborsParti[bid1][p].second<<" "<<Dis[j]<<endl; exit(1);
                        }
                    }
                }
//                for(int p=0;p<NeighborsParti[bid2].size();++p){
//                    if(bid1 == NeighborsParti[bid2][p].first){
//                        if(NeighborsParti[bid2][p].second < Dis[j]){
//                            NeighborsParti[bid2][p].second = Dis[j];
//                        }else if(NeighborsParti[bid2][p].second > Dis[j]){
//                            cout<<"Wrong parti! "<<bid1<<" "<<bid2<<" "<<NeighborsParti[bid2][p].second<<" "<<Dis[j]<<endl; exit(1);
//                        }
//                    }
//                }
            }
        }else{//if decrease update
            // Insert to overlay graph
            for(int j=0;j<IDs.size();++j){
                bid2=IDs[j];
                if(NeighborsOverlay[bid1].find(bid2)!=NeighborsOverlay[bid1].end()){
                    if(NeighborsOverlay[bid1][bid2]>Dis[j]){
                        sm->wait();///
                        overlayBatch.push_back({make_pair(bid1,bid2), make_pair(NeighborsOverlay[bid1][bid2],Dis[j])});
                        sm->notify();
//                        NeighborsOverlay[bid1][bid2]=Dis[j];
                    }
                    else if(NeighborsOverlay[bid1][bid2]<Dis[j]){
                        cout<<"Wrong overlay! "<<bid1<<" "<<bid2<<" "<<NeighborsOverlay[bid1][bid2]<<" "<<Dis[j]<<endl; exit(1);
                    }

                }else{
                    cout<<"Wrong overlay! Not found. "<<bid1<<" "<<bid2<<" "<<endl; exit(1);
                }
//                if(NeighborsOverlay[bid2].find(bid1)!=NeighborsOverlay[bid2].end()){
//                    if(NeighborsOverlay[bid2][bid1]>Dis[j]){
//                        NeighborsOverlay[bid2][bid1]=Dis[j];
//                    }else if(NeighborsOverlay[bid2][bid1]<Dis[j]){
//                        cout<<"Wrong overlay! "<<bid1<<" "<<bid2<<" "<<NeighborsOverlay[bid2][bid1]<<" "<<Dis[j]<<endl; exit(1);
//                    }
//                }else{
//                    cout<<"Wrong overlay! Not found. "<<bid1<<" "<<bid2<<" "<<endl; exit(1);
//                }
            }
            // Insert to partition graph
            for(int j=0;j<IDs.size();++j){
                bid2=IDs[j];
                for(int p=0;p<NeighborsParti[bid1].size();++p){
                    if(bid2 == NeighborsParti[bid1][p].first){
                        if(NeighborsParti[bid1][p].second > Dis[j]){
                            if(partiBatch.find(pid)==partiBatch.end()){//if not found
                                sm->wait();
                                partiBatch.insert({pid, vector<pair<pair<int, int>, pair<int, int>>>()});
                                partiBatch[pid].emplace_back(make_pair(bid1,bid2), make_pair(NeighborsParti[bid1][p].second,Dis[j]));
                                sm->notify();
                            }else{//if found
                                partiBatch[pid].emplace_back(make_pair(bid1,bid2), make_pair(NeighborsParti[bid1][p].second,Dis[j]));
                            }
//                            NeighborsParti[bid1][p].second = Dis[j];
                        }else if(NeighborsParti[bid1][p].second < Dis[j]){
                            cout<<"Wrong parti! "<<bid1<<" "<<bid2<<" "<<NeighborsParti[bid1][p].second<<" "<<Dis[j]<<endl; exit(1);
                        }
                    }
                }
//                for(int p=0;p<NeighborsParti[bid2].size();++p){
//                    if(bid1 == NeighborsParti[bid2][p].first){
//                        if(NeighborsParti[bid2][p].second > Dis[j]){
//                            NeighborsParti[bid2][p].second = Dis[j];
//                        }else if(NeighborsParti[bid2][p].second < Dis[j]){
//                            cout<<"Wrong parti! "<<bid1<<" "<<bid2<<" "<<NeighborsParti[bid2][p].second<<" "<<Dis[j]<<endl; exit(1);
//                        }
//                    }
//                }
            }
        }
    }
}

void Graph::AllPairBoundaryDisUpdatePartiMap(int pid, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>>& partiBatch, map<pair<int, int>, pair<int, int>>& overlayBatch) {
    int bid1,bid2;

    for(int i=0;i<BoundVertex[pid].size();++i){
        bid1=BoundVertex[pid][i];
        vector<int> IDs;
        for(int j=i+1;j<BoundVertex[pid].size();++j){
            bid2=BoundVertex[pid][j];
            IDs.push_back(bid2);
        }
        vector<int> Dis;
        BoundaryDijkstra(bid1, IDs, Dis, Neighbor);
        if(ifIncrease){
            // Insert to overlay graph
            for(int j=0;j<IDs.size();++j){
                bid2=IDs[j];
                if(NeighborsOverlay[bid1].find(bid2)!=NeighborsOverlay[bid1].end()){
                    if(NeighborsOverlay[bid1][bid2]<Dis[j]){
                        sm->wait();
//                        overlayBatch.insert({make_pair(bid1,bid2), make_pair(NeighborsOverlay[bid1][bid2],Dis[j])});
                        if(bid1<=bid2){
                            if(overlayBatch.find(make_pair(bid1,bid2))==overlayBatch.end()){//if not found
                                overlayBatch.insert({make_pair(bid1,bid2), make_pair(NeighborsOverlay[bid1][bid2],Dis[j])});
                            }else {
                                cout<<"exist update. "<<bid1<<" "<<bid2<<" "<<overlayBatch[make_pair(bid1,bid2)].second<<" "<<Dis[j]<<endl;
                                if(overlayBatch[make_pair(bid1,bid2)].first!=NeighborsOverlay[bid1][bid2] || overlayBatch[make_pair(bid1,bid2)].second>Dis[j]){
                                    overlayBatch[make_pair(bid1,bid2)]=make_pair(NeighborsOverlay[bid1][bid2],Dis[j]);
                                }

                            }
                        }
                        else{
                            if(overlayBatch.find(make_pair(bid2,bid1))==overlayBatch.end()){//if not found
                                overlayBatch.insert({make_pair(bid2,bid1), make_pair(NeighborsOverlay[bid1][bid2],Dis[j])});
                            }else {
                                cout<<"exist update. "<<bid2<<" "<<bid1<<" "<<overlayBatch[make_pair(bid2,bid1)].second<<" "<<Dis[j]<<endl;
                                if(overlayBatch[make_pair(bid2,bid1)].first!=NeighborsOverlay[bid1][bid2] || overlayBatch[make_pair(bid2,bid1)].second>Dis[j]){
                                    overlayBatch[make_pair(bid2,bid1)]=make_pair(NeighborsOverlay[bid1][bid2],Dis[j]);
                                }

                            }
                        }
                        sm->notify();
//                        NeighborsOverlay[bid1][bid2]=Dis[j];
                    }
                    else if(NeighborsOverlay[bid1][bid2]>Dis[j]){
                        cout<<"Wrong overlay! "<<bid1<<" "<<bid2<<" "<<NeighborsOverlay[bid1][bid2]<<" "<<Dis[j]<<endl; exit(1);
                    }
                }else{
                    cout<<"Wrong overlay! Not found. "<<bid1<<" "<<bid2<<" "<<endl; exit(1);
                }
//                if(NeighborsOverlay[bid2].find(bid1)!=NeighborsOverlay[bid2].end()){
//                    if(NeighborsOverlay[bid2][bid1]<Dis[j]){
////                        NeighborsOverlay[bid2][bid1]=Dis[j];
//                    }else if(NeighborsOverlay[bid2][bid1]>Dis[j]){
//                        cout<<"Wrong overlay! "<<bid1<<" "<<bid2<<" "<<NeighborsOverlay[bid2][bid1]<<" "<<Dis[j]<<endl; exit(1);
//                    }
//
//                }else{
//                    cout<<"Wrong overlay! Not found. "<<bid1<<" "<<bid2<<" "<<endl; exit(1);
//                }
            }
            // Insert to partition graph
            for(int j=0;j<IDs.size();++j){
                bid2=IDs[j];
                for(int p=0;p<NeighborsParti[bid1].size();++p){
                    if(bid2 == NeighborsParti[bid1][p].first){
                        if(NeighborsParti[bid1][p].second < Dis[j]){
                            if(partiBatch.find(pid)==partiBatch.end()){//if not found
                                sm->wait();
                                partiBatch.insert({pid, vector<pair<pair<int, int>, pair<int, int>>>()});
                                partiBatch[pid].emplace_back(make_pair(bid1,bid2), make_pair(NeighborsParti[bid1][p].second,Dis[j]));
                                sm->notify();
                            }else{//if found
                                partiBatch[pid].emplace_back(make_pair(bid1,bid2), make_pair(NeighborsParti[bid1][p].second,Dis[j]));
                            }
//                            NeighborsParti[bid1][p].second = Dis[j];
                        }else if(NeighborsParti[bid1][p].second > Dis[j]){
                            cout<<"Wrong parti! "<<bid1<<" "<<bid2<<" "<<NeighborsParti[bid1][p].second<<" "<<Dis[j]<<endl; exit(1);
                        }
                    }
                }
//                for(int p=0;p<NeighborsParti[bid2].size();++p){
//                    if(bid1 == NeighborsParti[bid2][p].first){
//                        if(NeighborsParti[bid2][p].second < Dis[j]){
//                            NeighborsParti[bid2][p].second = Dis[j];
//                        }else if(NeighborsParti[bid2][p].second > Dis[j]){
//                            cout<<"Wrong parti! "<<bid1<<" "<<bid2<<" "<<NeighborsParti[bid2][p].second<<" "<<Dis[j]<<endl; exit(1);
//                        }
//                    }
//                }
            }
        }else{//if decrease update
            // Insert to overlay graph
            for(int j=0;j<IDs.size();++j){
                bid2=IDs[j];
                if(NeighborsOverlay[bid1].find(bid2)!=NeighborsOverlay[bid1].end()){
                    if(NeighborsOverlay[bid1][bid2]>Dis[j]){
                        sm->wait();///
                        if(bid1<=bid2){
                            if(overlayBatch.find(make_pair(bid1,bid2))==overlayBatch.end()){//if not found
                                overlayBatch.insert({make_pair(bid1,bid2), make_pair(NeighborsOverlay[bid1][bid2],Dis[j])});
                            }else {
                                cout<<"exist update. "<<bid1<<" "<<bid2<<" "<<overlayBatch[make_pair(bid1,bid2)].second<<" "<<Dis[j]<<endl;
                                if(overlayBatch[make_pair(bid1,bid2)].first!=NeighborsOverlay[bid1][bid2] || overlayBatch[make_pair(bid1,bid2)].second>Dis[j]){
                                    overlayBatch[make_pair(bid1,bid2)]=make_pair(NeighborsOverlay[bid1][bid2],Dis[j]);
                                }

                            }
                        }
                        else{
                            if(overlayBatch.find(make_pair(bid2,bid1))==overlayBatch.end()){//if not found
                                overlayBatch.insert({make_pair(bid2,bid1), make_pair(NeighborsOverlay[bid1][bid2],Dis[j])});
                            }else {
                                cout<<"exist update. "<<bid2<<" "<<bid1<<" "<<overlayBatch[make_pair(bid2,bid1)].second<<" "<<Dis[j]<<endl;
                                if(overlayBatch[make_pair(bid2,bid1)].first!=NeighborsOverlay[bid1][bid2] || overlayBatch[make_pair(bid2,bid1)].second>Dis[j]){
                                    overlayBatch[make_pair(bid2,bid1)]=make_pair(NeighborsOverlay[bid1][bid2],Dis[j]);
                                }

                            }
                        }
                        sm->notify();
//                        NeighborsOverlay[bid1][bid2]=Dis[j];
                    }
                    else if(NeighborsOverlay[bid1][bid2]<Dis[j]){
                        cout<<"Wrong overlay! "<<bid1<<" "<<bid2<<" "<<NeighborsOverlay[bid1][bid2]<<" "<<Dis[j]<<endl; exit(1);
                    }

                }else{
                    cout<<"Wrong overlay! Not found. "<<bid1<<" "<<bid2<<" "<<endl; exit(1);
                }
//                if(NeighborsOverlay[bid2].find(bid1)!=NeighborsOverlay[bid2].end()){
//                    if(NeighborsOverlay[bid2][bid1]>Dis[j]){
//                        NeighborsOverlay[bid2][bid1]=Dis[j];
//                    }else if(NeighborsOverlay[bid2][bid1]<Dis[j]){
//                        cout<<"Wrong overlay! "<<bid1<<" "<<bid2<<" "<<NeighborsOverlay[bid2][bid1]<<" "<<Dis[j]<<endl; exit(1);
//                    }
//                }else{
//                    cout<<"Wrong overlay! Not found. "<<bid1<<" "<<bid2<<" "<<endl; exit(1);
//                }
            }
            // Insert to partition graph
            for(int j=0;j<IDs.size();++j){
                bid2=IDs[j];
                for(int p=0;p<NeighborsParti[bid1].size();++p){
                    if(bid2 == NeighborsParti[bid1][p].first){
                        if(NeighborsParti[bid1][p].second > Dis[j]){
                            if(partiBatch.find(pid)==partiBatch.end()){//if not found
                                sm->wait();
                                partiBatch.insert({pid, vector<pair<pair<int, int>, pair<int, int>>>()});
                                partiBatch[pid].emplace_back(make_pair(bid1,bid2), make_pair(NeighborsParti[bid1][p].second,Dis[j]));
                                sm->notify();
                            }else{//if found
                                partiBatch[pid].emplace_back(make_pair(bid1,bid2), make_pair(NeighborsParti[bid1][p].second,Dis[j]));
                            }
//                            NeighborsParti[bid1][p].second = Dis[j];
                        }else if(NeighborsParti[bid1][p].second < Dis[j]){
                            cout<<"Wrong parti! "<<bid1<<" "<<bid2<<" "<<NeighborsParti[bid1][p].second<<" "<<Dis[j]<<endl; exit(1);
                        }
                    }
                }
//                for(int p=0;p<NeighborsParti[bid2].size();++p){
//                    if(bid1 == NeighborsParti[bid2][p].first){
//                        if(NeighborsParti[bid2][p].second > Dis[j]){
//                            NeighborsParti[bid2][p].second = Dis[j];
//                        }else if(NeighborsParti[bid2][p].second < Dis[j]){
//                            cout<<"Wrong parti! "<<bid1<<" "<<bid2<<" "<<NeighborsParti[bid2][p].second<<" "<<Dis[j]<<endl; exit(1);
//                        }
//                    }
//                }
            }
        }
    }
}

//Dijkstra's algorithm
int Graph::BoundaryDijkstra(int ID1, vector<int> & IDs, vector<int>& Dis, vector<vector<pair<vertex,int>>> &Neighbor){
    benchmark::heap<2, int, int> pqueue(node_num);
    pqueue.update(ID1,0);

    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);
    vector<int> prece(node_num, 0);
    distance[ID1]=0;
    int topNodeID, topNodeDis;
    int NNodeID,NWeigh;
    int ID2;
    unordered_map<int,int> results;
    unordered_set<int> targets;
    targets.insert(IDs.begin(),IDs.end());

    int d=INF;//initialize d to infinite for the unreachable case

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        if(targets.find(topNodeID)!=targets.end()){//if found
            results.insert({topNodeID,distance[topNodeID]});
            targets.erase(topNodeID);
            if(targets.empty()){
                break;
            }
        }
        closed[topNodeID]=true;

        for(auto it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
            NNodeID=(*it).first;
            NWeigh=(*it).second+topNodeDis;
            if(!closed[NNodeID]){
                if(distance[NNodeID]>NWeigh){
                    distance[NNodeID]=NWeigh;
                    pqueue.update(NNodeID, NWeigh);
                    prece[NNodeID]=topNodeID;
                }
            }
        }
    }
    if(!targets.empty()){
        cout<<"Not find all targets! "<<targets.size()<<endl; exit(1);
    }
    for(int i=0;i<IDs.size();++i){
        if(results.find(IDs[i])!=results.end()){
            Dis.push_back(results[IDs[i]]);
        }else{
            cout<<"Not find vertex "<<IDs[i]<<endl; exit(1);
        }
    }
    //retrieve path
//    RetrievePath(ID1, ID2, prece);

    return d;
}

//Function of constructing tree index for partitions
void Graph::Construct_PartiIndex(bool ifParallel, bool ifLabelC){
    //for H2H update
    SCconNodesMTP.assign(node_num, map<int, vector<pair<int,int>>>());
    VidtoTNidP.assign(node_num,vector<int>());
    NeighborCon.assign(node_num, vector<pair<int,pair<int,int>>>());
//    _DD_.assign(node_num,0); //_DD2_.assign(node_num,0);
//    DD.assign(node_num,0); //DD2.assign(node_num,0);
    Trees.assign(partiNum,vector<Node>());
    toRMQs.assign(partiNum,vector<int>());
    RMQIndexs.assign(partiNum,vector<vector<int>>());
    ranks.assign(partiNum,vector<int>());
    heightMaxs.assign(partiNum,0);
    ProBeginVertexSetParti.assign(partiNum,vector<int>());
    vertexIDChLParti.assign(partiNum,set<int>());
//    BoundEdges.assign(node_num,map<int,pair<int,int>>());

    //initialize E
    map<int,pair<int,int>> m;
    E.assign(node_num,m);
    for(int i=0;i<NeighborsParti.size();i++){
        for(auto it=NeighborsParti[i].begin();it!=NeighborsParti[i].end();++it){
            E[i].insert(make_pair(it->first,make_pair(it->second,1)));
        }
    }
    for(int pid=0;pid<partiNum;++pid){
        ranks[pid].assign(PartiVertex[pid].size(),-1);
    }

    if(ifParallel){
//        cout<<"Multiple thread computation for partition index construction!"<<endl;
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
            if(ifLabelC){
                for(auto j=0;j<processID.size();++j){
                    thread.add_thread(new boost::thread(&Graph::ConstructPH2H_PartiV, this, boost::ref(processID[j]), boost::ref(Trees), boost::ref(ranks), boost::ref(SCconNodesMTP), boost::ref(VidtoTNidP), boost::ref(heightMaxs), boost::ref(toRMQs), boost::ref(RMQIndexs)));
                }
            }else{
                for(auto j=0;j<processID.size();++j){
                    thread.add_thread(new boost::thread(&Graph::ConstructPH2H_PartiVCH, this, boost::ref(processID[j]), boost::ref(Trees), boost::ref(ranks), boost::ref(SCconNodesMTP), boost::ref(VidtoTNidP), boost::ref(heightMaxs), boost::ref(toRMQs), boost::ref(RMQIndexs)));
                }
            }

            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                if(ifLabelC){
                    thread.add_thread(new boost::thread(&Graph::ConstructPH2H_Parti, this, j, boost::ref(Trees), boost::ref(ranks), boost::ref(SCconNodesMTP), boost::ref(VidtoTNidP), boost::ref(heightMaxs), boost::ref(toRMQs), boost::ref(RMQIndexs)));
                }else{
                    thread.add_thread(new boost::thread(&Graph::H2HCreateTree_Parti, this, j, boost::ref(Trees[j]), boost::ref(ranks[j]), boost::ref(SCconNodesMTP), boost::ref(VidtoTNidP), boost::ref(heightMaxs)));
                }

            }
            thread.join_all();
        }

    }
    else{
        cout<<"Single thread computation!"<<endl;
        //single-thread
        for(int pid=0;pid<partiNum;++pid){
            cout<<"Partition "<<pid<<endl;
            if(ifLabelC){
                ConstructPH2H_Parti(pid, Trees, ranks, SCconNodesMTP, VidtoTNidP, heightMaxs, toRMQs, RMQIndexs);
            }else{
                H2HCreateTree_Parti(pid, Trees[pid], ranks[pid], SCconNodesMTP, VidtoTNidP, heightMaxs);
            }

        }
    }

    vector<int> treeSize;
    int aveHeight=0;
    for(int i=0;i<partiNum;++i){
        treeSize.emplace_back(Trees[i].size());
        aveHeight+=heightMaxs[i];
    }
    cout<<"Partition graph! Maximum tree node number: "<< *max_element(treeSize.begin(),treeSize.end()) <<" ; Maximum tree height: "<< *max_element(heightMaxs.begin(),heightMaxs.end())<<" ; Average tree height: "<< aveHeight/partiNum<< endl;
}


//function of vertex allocation
void Graph::ThreadDistribute(vector<int>& vertices, vector<vector<int>>& processID){
//        processID.assign(threadnum, vector<NodeId>());
    int pid=0;
    for(int i=0;i<vertices.size();++i){
        pid=i%processID.size();
        processID[pid].emplace_back(vertices[i]);
    }
}
void Graph::ConstructBoundaryShortcut(int pid, bool ifCH){
    //boundary edges
    int ID1,ID2,weight;
    if(ifFullOpt){//use overlay simplification optimization
        if(ifCH){
            for(int i=0;i<BoundVertex[pid].size();i++){
                ID1=BoundVertex[pid][i];
                if(fullyConnected[ID1]) continue;
                for(int j=i+1;j<BoundVertex[pid].size();j++){
                    ID2=BoundVertex[pid][j];
                    if(fullyConnected[ID2]) continue;
                    weight=QueryCHPartition(ID1,ID2,pid,Trees);
                    NeighborsOverlay[ID1][ID2]=weight;
                    NeighborsOverlay[ID2][ID1]=weight;
                }
            }
        }
        else{
            for(int i=0;i<BoundVertex[pid].size();i++){
                ID1=BoundVertex[pid][i];
                if(fullyConnected[ID1]) continue;
                for(int j=i+1;j<BoundVertex[pid].size();j++){
                    ID2=BoundVertex[pid][j];
                    if(fullyConnected[ID2]) continue;
                    weight=QueryH2HPartition(ID1,ID2,pid);
                    NeighborsOverlay[ID1][ID2]=weight;
                    NeighborsOverlay[ID2][ID1]=weight;
                }
            }
        }
    }
    else{//not use overlay simplification optimization
        if(ifCH){
            for(int i=0;i<BoundVertex[pid].size();i++){
                ID1=BoundVertex[pid][i];
                for(int j=i+1;j<BoundVertex[pid].size();j++){
                    ID2=BoundVertex[pid][j];
                    weight=QueryCHPartition(ID1,ID2,pid,Trees);
                    NeighborsOverlay[ID1][ID2]=weight;
                    NeighborsOverlay[ID2][ID1]=weight;
                }
            }
        }
        else{
            for(int i=0;i<BoundVertex[pid].size();i++){
                ID1=BoundVertex[pid][i];
                for(int j=i+1;j<BoundVertex[pid].size();j++){
                    ID2=BoundVertex[pid][j];
                    weight=QueryH2HPartition(ID1,ID2,pid);
                    NeighborsOverlay[ID1][ID2]=weight;
                    NeighborsOverlay[ID2][ID1]=weight;
                }
            }
        }
    }


}

void Graph::ConstructBoundaryShortcutNoAllPair(int pid){
    //boundary edges
    int ID1,ID2,weight;
    if(ifFullOpt){//with overlay optimization
        for(int i=0;i<BoundVertex[pid].size();i++){
            ID1=BoundVertex[pid][i];
            if(fullyConnected[ID1]) continue;
            if(!NeighborCon[ID1].empty()){
                for(auto it=NeighborCon[ID1].begin();it!=NeighborCon[ID1].end();++it){
                    ID2=it->first; weight=it->second.first;
                    if(fullyConnected[ID2]) continue;
                    if(!PartiTag[ID2].second){//if ID2 is not boundary vertex
                        cout<<"Wrong for this shortcut! "<<ID1<<"("<<PartiTag[ID1].first<<","<<PartiTag[ID1].second<<") "<<ID2<<"("<<PartiTag[ID2].first<<","<<PartiTag[ID2].second<<")"<<endl; exit(1);//?
                    }
                    if(NeighborsOverlay[ID1].find(ID2)==NeighborsOverlay[ID1].end()){//if not found
                        NeighborsOverlay[ID1][ID2]=weight;
                        NeighborsOverlay[ID2][ID1]=weight;
                    }else if(NeighborsOverlay[ID1][ID2]>weight){
                        NeighborsOverlay[ID1][ID2]=weight;
                        NeighborsOverlay[ID2][ID1]=weight;
                    }
                }
            }
        }
    }else{
        for(int i=0;i<BoundVertex[pid].size();i++){
            ID1=BoundVertex[pid][i];
            if(!NeighborCon[ID1].empty()){
                for(auto it=NeighborCon[ID1].begin();it!=NeighborCon[ID1].end();++it){
                    ID2=it->first; weight=it->second.first;
                    if(!PartiTag[ID2].second){//if ID2 is not boundary vertex
                        cout<<"Wrong for this shortcut! "<<ID1<<"("<<PartiTag[ID1].first<<","<<PartiTag[ID1].second<<") "<<ID2<<"("<<PartiTag[ID2].first<<","<<PartiTag[ID2].second<<")"<<endl; exit(1);//?
                    }
                    if(NeighborsOverlay[ID1].find(ID2)==NeighborsOverlay[ID1].end()){//if not found
                        NeighborsOverlay[ID1][ID2]=weight;
                        NeighborsOverlay[ID2][ID1]=weight;
                    }else if(NeighborsOverlay[ID1][ID2]>weight){
                        NeighborsOverlay[ID1][ID2]=weight;
                        NeighborsOverlay[ID2][ID1]=weight;
                    }
                }
            }
        }
    }

}
void Graph::ConstructBoundaryShortcutV(vector<int> & p, bool ifAllPair, bool ifCH){
    if(ifAllPair){
        for(int i=0;i<p.size();++i){
            ConstructBoundaryShortcut(p[i], ifCH);
        }
    }else{
        for(int i=0;i<p.size();++i){
            ConstructBoundaryShortcutNoAllPair(p[i]);
        }
    }

}

void Graph::Construct_OverlayGraph(bool ifParallel, bool ifCH){
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
                thread.add_thread(new boost::thread(&Graph::ConstructBoundaryShortcutV, this, boost::ref(processID[j]) ,true, ifCH ));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                thread.add_thread(new boost::thread(&Graph::ConstructBoundaryShortcut, this, j, ifCH));
            }
            thread.join_all();
        }

    }
    else{
        //single thread
        for(int k=0;k<partiNum;k++){
            ConstructBoundaryShortcut(k, ifCH);
        }
    }


}

void Graph::Construct_OverlayGraphNoAllPair(bool ifParallel, bool ifCH){
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
                thread.add_thread(new boost::thread(&Graph::ConstructBoundaryShortcutV, this, boost::ref(processID[j]), false, false ));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                thread.add_thread(new boost::thread(&Graph::ConstructBoundaryShortcutNoAllPair, this, j));
            }
            thread.join_all();
        }

    }
    else{
        //single thread
        for(int k=0;k<partiNum;k++){
            ConstructBoundaryShortcutNoAllPair(k);
        }
    }


}

//Function of constructing tree index for overlay grpah
void Graph::Construct_OverlayIndex(bool ifLabelC){
    //Create tree for partition
    H2HCreateTree_Overlay();
    if(ifLabelC){
        //Create labels for partition
        H2HCreateIndex_Overlay();
    }
}

void Graph::PreConstructAllPairs(bool ifParallel){
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
                thread.add_thread(new boost::thread(&Graph::PreConstructAllPairsPartiV, this, boost::ref(processID[j]) ));
            }

            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                thread.add_thread(new boost::thread(&Graph::PreConstructAllPairsParti, this, j));
            }
            thread.join_all();
        }

    }
    else{
        cout<<"Single thread computation!"<<endl;
        //single-thread
        for(int pid=0;pid<partiNum;++pid){
//            cout<<"Partition "<<pid<<endl;
            PreConstructAllPairsParti(pid);
        }
    }
}

void Graph::PreConstructAllPairsPartiV(vector<int> & p){
    for(int i=0;i<p.size();++i){
        PreConstructAllPairsParti(p[i]);
    }
}

void Graph::PreConstructAllPairsParti(int pid){
    int ID1,ID2;
    bool flagEdge;
    int newENum=0;
    for(int i=0;i<BoundVertex[pid].size();++i){
        ID1=BoundVertex[pid][i];
        for(int j=i+1;j<BoundVertex[pid].size();++j){
            ID2=BoundVertex[pid][j];
            flagEdge=false;
            for(auto it=NeighborsParti[ID1].begin();it!=NeighborsParti[ID1].end();++it){
                if(ID2==it->first){
                    flagEdge=true;
                    break;
                }
            }
            if(!flagEdge){//if not exist edge
                NeighborsParti[ID1].emplace_back(ID2,INF/3);//insert edge with weight equals INF
                NeighborsParti[ID2].emplace_back(ID1,INF/3);//insert edge with weight equals INF
                newENum+=2;
            }
        }
    }
//    cout<<"New edge number of partition "<<pid<<" : "<<newENum<<endl;
}


void Graph::ConstructPartitionPostIndexOpt(bool ifParallel){

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

                thread.add_thread(new boost::thread(&Graph::RefreshBoundaryEdgesAndLabelingPartiV, this, boost::ref(processID[j])));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                thread.add_thread(new boost::thread(&Graph::RefreshBoundaryEdgesAndLabelingParti, this, j));
            }
            thread.join_all();
        }
    }
    else{
        // single thread
        for(int k=0;k<partiNum;k++){
            cout<<"Repairing partition "<<k<<endl;
            RefreshBoundaryEdgesAndLabelingParti(k);
        }
    }

}

void Graph::RefreshBoundaryEdgesAndLabelingPartiV(vector<int>& p){
    for(int i=0;i<p.size();++i){
        RefreshBoundaryEdgesAndLabelingParti(p[i]);
    }
}

void Graph::RefreshBoundaryEdgesAndLabelingParti(int pid) {
    int ID;
    vector<pair<pair<int,int>,pair<int,int>>> wBatch;

    for(int i=0;i<PartiVertex[pid].size();++i){
        ID=PartiVertex[pid][i];
        NeighborsPartiPost[ID].insert(NeighborsParti[ID].begin(),NeighborsParti[ID].end());
    }

    int ID1,ID2;
    int wlocal, woverlay;
    int k;
    for(int i=0;i<BoundVertex[pid].size();++i){
        ID1=BoundVertex[pid][i];
        for(int j=i+1;j<BoundVertex[pid].size();++j){
            ID2=BoundVertex[pid][j];
            if(NeighborsPartiPost[ID1].find(ID2)==NeighborsPartiPost[ID1].end()){//if not found
                cout<<pid<<": boundary edge between "<<ID1<<" and "<<ID2<<" does not exist! "<<endl; exit(1);
            }
            wlocal=NeighborsPartiPost[ID1][ID2];
//            wlocal=-1;
//            for(k=0;k<NeighborsParti[ID1].size();++k){
//                if(NeighborsParti[ID1][k].first==ID2){
//                    wlocal=NeighborsParti[ID1][k].second;
//                    break;
//                }
//            }
            woverlay= QueryCore(ID1,ID2);
//            int d2= Dijkstra(ID1,ID2,Neighbor);
//            if(d2!=weight){
//                cout<<"Incorrect! "<<ID1<<" "<<ID2<<" "<<weight<<" "<<d2<<endl; exit(1);
//            }
            if(woverlay<wlocal){
//                BoundEdges[ID1].insert({ID2, make_pair(woverlay,-1)});
//                BoundEdges[ID2].insert({ID1, make_pair(woverlay,-1)});
//                NeighborsParti[ID1][k].second=woverlay;
                wBatch.emplace_back(make_pair(ID1,ID2), make_pair(wlocal,woverlay));
            }else if(woverlay>wlocal){
                cout<<"Something wrong happens. woverlay is larger! "<<wlocal<<" "<<woverlay<<endl; exit(1);
            }
//            else{
//                BoundEdges[ID1].insert({ID2, make_pair(woverlay,pid)});
//                BoundEdges[ID2].insert({ID1, make_pair(woverlay,pid)});
//            }

        }
    }
//    cout<<"wBatch size: "<<wBatch.size()<<endl;

//    E.assign(node_num,map<int,pair<int,int>>());
//    for(int i=0;i<PartiVertex[pid].size();i++){
//        int id=PartiVertex[pid][i];
//        for(auto it=NeighborsParti[id].begin();it!=NeighborsParti[id].end();++it){
//            E[id].insert(make_pair(it->first,make_pair(it->second,1)));
//        }
//    }
//    H2HCreateTree_Parti(pid, Trees[pid], ranks[pid], SCconNodesMTP, VidtoTNidP, heightMaxs);
    vector<pair<pair<int,int>,int>> updatedSC;
    //update shortcuts
    DecreasePartiBatchForOpt(pid, wBatch, NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid], updatedSC, false, true);
    //Create LCA index
    makeRMQCoreP(pid, toRMQs, RMQIndexs, Trees);
    //Create labels for partition
    H2HCreateIndex_Parti(pid, Trees[pid], ranks[pid]);
}

void Graph::ConstructPartitionPost(bool ifParallel, bool ifCH){
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
                thread.add_thread(new boost::thread(&Graph::ConstructPostPartiV, this, boost::ref(processID[j]), ifCH));

            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                thread.add_thread(new boost::thread(&Graph::ConstructPostParti, this, j, ifCH));

            }
            thread.join_all();
        }
    }
    else{
        // single thread
        for(int k=0;k<partiNum;k++){
            cout<<"Repairing partition "<<k<<endl;
            ConstructPostParti(k,ifCH);

        }
    }
}

void Graph::ConstructPostParti(int pid, bool ifCH){
    int ID;
    for(int i=0;i<PartiVertex[pid].size();++i){
        ID=PartiVertex[pid][i];
        NeighborsPartiPost[ID].insert(NeighborsParti[ID].begin(),NeighborsParti[ID].end());
    }
    int ID1,ID2,weight=-1;
    if(ifCH){
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
                weight= QueryCoreCH(ID1,ID2);
                NeighborsPartiPost[ID1][ID2]=weight;
                NeighborsPartiPost[ID2][ID1]=weight;
            }
        }
    }
    else{
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
                weight= QueryCore(ID1,ID2);
                NeighborsPartiPost[ID1][ID2]=weight;
                NeighborsPartiPost[ID2][ID1]=weight;
            }
        }
    }

}

void Graph::ConstructPostPartiV(vector<int>& p, bool ifCH){
    for(int i=0;i<p.size();++i){
        ConstructPostParti(p[i], ifCH);
    }
}
//old version
void Graph::ConstructPartitionPostIndex(bool ifParallel, bool ifLabelC){
    SCconNodesMTPost.assign(node_num, map<int, vector<pair<int,int>>>());
    VidtoTNidPost.assign(node_num,vector<int>());
    NeighborCon.assign(node_num, vector<pair<int,pair<int,int>>>());
//    _DD_.assign(node_num,0); //_DD2_.assign(node_num,0);
//    DD.assign(node_num,0); //DD2.assign(node_num,0);
    TreesPost.assign(partiNum,vector<Node>());
    toRMQsPost.assign(partiNum,vector<int>());
    RMQIndexsPost.assign(partiNum,vector<vector<int>>());
    ranksPost.assign(partiNum,vector<int>());
    heightMaxsPost.assign(partiNum,0);

    //initialize E
    map<int,pair<int,int>> m;
    E.assign(node_num,m);
    for(int i=0;i<NeighborsPartiPost.size();i++){
        for(auto it=NeighborsPartiPost[i].begin();it!=NeighborsPartiPost[i].end();++it){
            E[i].insert(make_pair(it->first,make_pair(it->second,1)));
        }
    }
    for(int pid=0;pid<partiNum;++pid){
        ranksPost[pid].assign(PartiVertex[pid].size(),-1);
    }

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
                if(ifLabelC){
                    thread.add_thread(new boost::thread(&Graph::ConstructPH2H_PartiV, this, boost::ref(processID[j]), boost::ref(TreesPost), boost::ref(ranksPost), boost::ref(SCconNodesMTPost), boost::ref(VidtoTNidPost), boost::ref(heightMaxsPost), boost::ref(toRMQsPost), boost::ref(RMQIndexsPost) ));
                }else{
                    thread.add_thread(new boost::thread(&Graph::ConstructPH2H_PartiVCH, this, boost::ref(processID[j]), boost::ref(TreesPost), boost::ref(ranksPost), boost::ref(SCconNodesMTPost), boost::ref(VidtoTNidPost), boost::ref(heightMaxsPost), boost::ref(toRMQsPost), boost::ref(RMQIndexsPost) ));
                }

            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                if(ifLabelC){
                    thread.add_thread(new boost::thread(&Graph::ConstructPH2H_Parti, this, j, boost::ref(TreesPost), boost::ref(ranksPost), boost::ref(SCconNodesMTPost), boost::ref(VidtoTNidPost), boost::ref(heightMaxsPost), boost::ref(toRMQsPost), boost::ref(RMQIndexsPost)));
                }else{
                    thread.add_thread(new boost::thread(&Graph::H2HCreateTree_Parti, this, j, boost::ref(TreesPost[j]), boost::ref(ranksPost[j]), boost::ref(SCconNodesMTPost), boost::ref(VidtoTNidPost), boost::ref(heightMaxsPost)));
                }

            }
            thread.join_all();
        }
    }
    else{
        cout<<"Single thread computation!"<<endl;
        //single-thread
        for(int pid=0;pid<partiNum;++pid){
//            cout<<"Partition "<<pid<<endl;
            if(ifLabelC){
                ConstructPH2H_Parti(pid,TreesPost, ranksPost, SCconNodesMTPost, VidtoTNidPost, heightMaxsPost, toRMQsPost, RMQIndexsPost);
            }else{
                H2HCreateTree_Parti(pid,TreesPost[pid], ranksPost[pid], SCconNodesMTPost, VidtoTNidPost, heightMaxsPost);
            }
        }
    }
    vector<int> treeSize;
    int aveHeight=0;
    for(int i=0;i<partiNum;++i){
        treeSize.emplace_back(TreesPost[i].size());
        aveHeight+=heightMaxsPost[i];
    }
    cout<<"Post-boundary partition graph! Maximum tree node number: "<< *max_element(treeSize.begin(),treeSize.end()) <<" ; Maximum tree height: "<< *max_element(heightMaxsPost.begin(),heightMaxsPost.end())<<" ; Average tree height: "<< aveHeight/partiNum<< endl;
}

//Function of repair the partition index
void Graph::Repair_PartiIndex(bool ifParallel, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch){
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
                thread.add_thread(new boost::thread(&Graph::RepairPartitionIndexV, this, boost::ref(processID[j]), ifIncrease, boost::ref(partiBatch), false));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            if(ifIncrease){//increase update
                boost::thread_group thread;
                for(auto j=0;j<partiNum;++j){
                    thread.add_thread(new boost::thread(&Graph::RepairPartitionIndexIncrease, this, j, boost::ref(partiBatch)));
                }
                thread.join_all();
            }
            else{//decrease update
                boost::thread_group thread;
                for(auto j=0;j<partiNum;++j){
                    thread.add_thread(new boost::thread(&Graph::RepairPartitionIndexDecrease, this, j, boost::ref(partiBatch)));
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
                RepairPartitionIndexIncrease(k,partiBatch);
            }
        }
        else{
            for(int k=0;k<partiNum;k++){
//                cout<<"Repairing partition "<<k<<endl;
                RepairPartitionIndexDecrease(k,partiBatch);
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

void Graph::Repair_PartiIndexForOpt(bool ifParallel, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch){
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
                thread.add_thread(new boost::thread(&Graph::RepairPartitionIndexV, this, boost::ref(processID[j]), ifIncrease, boost::ref(partiBatch), true));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            if(ifIncrease){//increase update
                boost::thread_group thread;
                for(auto j=0;j<partiNum;++j){
                    thread.add_thread(new boost::thread(&Graph::RepairPartitionIndexIncreaseForOpt, this, j, boost::ref(partiBatch)));
                }
                thread.join_all();
            }
            else{//decrease update
                boost::thread_group thread;
                for(auto j=0;j<partiNum;++j){
                    thread.add_thread(new boost::thread(&Graph::RepairPartitionIndexDecreaseForOpt, this, j, boost::ref(partiBatch)));
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
                RepairPartitionIndexIncreaseForOpt(k,partiBatch);
            }
        }
        else{
            for(int k=0;k<partiNum;k++){
//                cout<<"Repairing partition "<<k<<endl;
                RepairPartitionIndexDecreaseForOpt(k,partiBatch);
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

//Function of repair the partition index
void Graph::RepairPartitionIndexV(vector<int>& p, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch, bool ifOpt) {
    if(ifOpt){
        if(ifIncrease){
            for(int i=0;i<p.size();++i){
                RepairPartitionIndexIncreaseForOpt(p[i], partiBatch);
            }
        }
        else{
            for(int i=0;i<p.size();++i){
                RepairPartitionIndexDecreaseForOpt(p[i], partiBatch);
            }
        }
    }else{
        if(ifIncrease){
            for(int i=0;i<p.size();++i){
                RepairPartitionIndexIncrease(p[i], partiBatch);
            }
        }
        else{
            for(int i=0;i<p.size();++i){
                RepairPartitionIndexDecrease(p[i], partiBatch);
            }
        }
    }


}
void Graph::RepairPartitionIndexDecrease(int pid, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch) {
    int ID1,ID2,woverlay;
    int wlocal=INF;
    //boundary edges
    vector<pair<pair<int,int>,pair<int,int>>> weightsParti;//collect the changed edges on overlay graph
    weightsParti.clear();
    vector<pair<pair<int,int>,pair<int,int>>> weightsPartiInsert;//collect the changed edges on overlay graph

    if(partiBatch.find(pid)!=partiBatch.end()){//if found
        for(int i=0;i<partiBatch[pid].size();++i){
            ID1=partiBatch[pid][i].first.first, ID2=partiBatch[pid][i].first.second;
            if(NeighborsPartiPost[ID1].find(ID2)==NeighborsPartiPost[ID1].end()){//if not found
                weightsPartiInsert.emplace_back(partiBatch[pid][i]);
            }else{
                weightsParti.emplace_back(partiBatch[pid][i]);
            }
        }
//        weightsParti=partiBatch[pid];
    }
    if(algoChoice==CH){
        for(int i=0;i<BoundVertex[pid].size();i++){
            ID1=BoundVertex[pid][i];
            for(int j=i+1;j<BoundVertex[pid].size();j++){
                ID2=BoundVertex[pid][j];

                if(NeighborsPartiPost[ID1].find(ID2)!=NeighborsPartiPost[ID1].end()){//if found
                    wlocal=NeighborsPartiPost[ID1][ID2];
                }else{
                    cout<<"Not found edge e("<<ID1<<","<<ID2<<") in overlay graph!"<<endl; exit(1);
                }
                woverlay=QueryCoreCH(ID1,ID2);
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
    }else if(algoChoice==H2H){
        for(int i=0;i<BoundVertex[pid].size();i++){
            ID1=BoundVertex[pid][i];
            for(int j=i+1;j<BoundVertex[pid].size();j++){
                ID2=BoundVertex[pid][j];

                if(NeighborsPartiPost[ID1].find(ID2)!=NeighborsPartiPost[ID1].end()){//if found
                    wlocal=NeighborsPartiPost[ID1][ID2];
                }else{
                    cout<<"Not found edge e("<<ID1<<","<<ID2<<") in overlay graph!"<<endl; exit(1);
                }
                woverlay=QueryCore(ID1,ID2);
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
    }

    if(!weightsPartiInsert.empty()){
//        cout<<"weightsPartiInsert size: "<<weightsPartiInsert.size()<<endl;
        bool ifLabel=false;
        if(algoChoice==H2H){
            ifLabel=true;
        }
        EdgeInsertPartiBatchH2HPost(pid,weightsPartiInsert, NeighborsPartiPost, TreesPost[pid], ranksPost[pid], heightMaxsPost[pid], ifLabel);
    }
    if(!weightsParti.empty()){
//        cout<<"Size of weightsParti of partition "<<pid<<" : "<<weightsParti.size()<<endl;
        bool ifLabel=false;
        if(algoChoice==H2H){
            ifLabel=true;
        }
        DecreasePartiBatchPost(pid,weightsParti, NeighborsPartiPost, TreesPost[pid], ranksPost[pid], heightMaxsPost[pid], ifLabel);
        ifRepaired[pid]=true;
    }

}

void Graph::RepairPartitionIndexDecreaseForOpt(int pid, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch) {
    int ID1,ID2,woverlay;
    int wlocal=INF;
    //boundary edges
    vector<pair<pair<int,int>,pair<int,int>>> weightsParti;//collect the changed edges on overlay graph
    weightsParti.clear();

//    if(pid==12){
//        cout<<pid<<endl;
//    }

//    if(partiBatch.find(pid)!=partiBatch.end()){//if found
//        weightsParti=partiBatch[pid];
//    }

    for(int i=0;i<BoundVertex[pid].size();i++){
        ID1=BoundVertex[pid][i];
        for(int j=i+1;j<BoundVertex[pid].size();j++){
            ID2=BoundVertex[pid][j];
//            wlocal=-1;
//            for(auto it=NeighborsParti[ID1].begin();it!=NeighborsParti[ID1].end();++it){
//                if(it->first==ID2){
//                    wlocal=it->second;
//                    break;
//                }
//            }
            if(NeighborsPartiPost[ID1].find(ID2)==NeighborsPartiPost[ID1].end()){
                cout<<"Not found! "<<ID1<<" "<<ID2<<endl; exit(1);
            }
            wlocal=NeighborsPartiPost[ID1][ID2];
//            wlocal= QueryH2HPartition(ID1,ID2,pid);
            if(wlocal==-1){//if found
                cout<<"Not found edge e("<<ID1<<","<<ID2<<") in overlay graph!"<<endl; exit(1);
            }
            woverlay=QueryCore(ID1,ID2);
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
//        cout<<"Size of weightsParti (with new) of partition "<<pid<<" : "<<weightsParti.size()<<endl;
        bool ifPost=false;
        if(PSPStrategy>=PostBoundary){
            ifPost=true;
        }
        vector<pair<pair<int,int>,int>> updatedSC;
        DecreasePartiBatchForOpt(pid,weightsParti, NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid], updatedSC, ifPost, false);
        ifRepaired[pid]=true;
    }else if(partiBatch.find(pid)!=partiBatch.end()){
//        cout<<"Size of weightsParti of partition "<<pid<<" : "<<weightsParti.size()<<endl;
        DecreasePartiBatchLabel(TreesPost[pid], ranks[pid], heightMaxs[pid], ProBeginVertexSetParti[pid], vertexIDChLParti[pid]);
        ifRepaired[pid]=true;
    }

}

void Graph::RepairPartitionIndexIncrease(int pid, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch) {
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
    if(algoChoice==CH){
        for(int i=0;i<BoundVertex[pid].size();i++){
            ID1=BoundVertex[pid][i];
            for(int j=i+1;j<BoundVertex[pid].size();j++){
                ID2=BoundVertex[pid][j];
                if(NeighborsPartiPost[ID1].find(ID2)!=NeighborsPartiPost[ID1].end()){//if found
                    wlocal=NeighborsPartiPost[ID1][ID2];
                }else{
                    cout<<"Not found edge e("<<ID1<<","<<ID2<<") in overlay graph!"<<endl; exit(1);
                }
                woverlay=QueryCoreCH(ID1,ID2);
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
    }else if(algoChoice==H2H){
        for(int i=0;i<BoundVertex[pid].size();i++){
            ID1=BoundVertex[pid][i];
            for(int j=i+1;j<BoundVertex[pid].size();j++){
                ID2=BoundVertex[pid][j];
                if(NeighborsPartiPost[ID1].find(ID2)!=NeighborsPartiPost[ID1].end()){//if found
                    wlocal=NeighborsPartiPost[ID1][ID2];
                }else{
                    cout<<"Not found edge e("<<ID1<<","<<ID2<<") in overlay graph!"<<endl; exit(1);
                }
                woverlay=QueryCore(ID1,ID2);
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
    }
    for(auto it=updatesSet.begin();it!=updatesSet.end();++it){
        weightsParti.emplace_back(*it);
    }

    if(!weightsParti.empty()){
//        cout<<"Size of weightsParti of partition "<<pid<<" : "<<weightsParti.size()<<endl;
        bool ifPost=false;
        if(algoChoice==H2H){
            ifPost=true;
        }
        IncreasePartiBatchPost(pid, weightsParti, NeighborsPartiPost, TreesPost[pid], ranksPost[pid], heightMaxsPost[pid],SCconNodesMTPost,VidtoTNidPost, ifPost);
        ifRepaired[pid]=true;
    }

}

void Graph::RepairPartitionIndexIncreaseForOpt(int pid, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch) {
    int ID1,ID2,woverlay;
    int wlocal=INF;
    //boundary edges
    vector<pair<pair<int,int>,pair<int,int>>> weightsParti;//collect the changed edges on overlay graph
    weightsParti.clear();
    map<pair<int,int>,pair<int,int>> updatesSet;

//    if(partiBatch.find(pid)!=partiBatch.end()){//if found
////        cout<<"Find for "<<pid<<" "<<partiBatch[pid].size()<<endl;
////        weightsParti=partiBatch[pid];
//        for(auto it=partiBatch[pid].begin();it!=partiBatch[pid].end();++it){
//            ID1=it->first.first, ID2=it->first.second;
//            if(ID1<ID2){
//                if(updatesSet.find(make_pair(ID1,ID2))==updatesSet.end()){//if not found
//                    updatesSet.insert(*it);
//                }else{
//                    cout<<"Already exists! "<<ID1<<" "<<ID2<<endl; exit(1);
//                }
////                    weightsParti.emplace_back(make_pair(ID1,ID2), make_pair(wlocal, woverlay));//oldW,newW
//            }else{
//                if(updatesSet.find(make_pair(ID2,ID1))==updatesSet.end()){//if not found
//                    updatesSet.insert(*it);
//                }else{
//                    cout<<"Already exists! "<<ID1<<" "<<ID2<<endl; exit(1);
//                }
////                    weightsParti.emplace_back(make_pair(ID2,ID1), make_pair(wlocal, woverlay));//oldW,newW
//            }
//        }
//    }
    for(int i=0;i<BoundVertex[pid].size();i++){
        ID1=BoundVertex[pid][i];
        for(int j=i+1;j<BoundVertex[pid].size();j++){
            ID2=BoundVertex[pid][j];
            wlocal=-1;
            for(auto it=NeighborsParti[ID1].begin();it!=NeighborsParti[ID1].end();++it){
                if(it->first==ID2){
                    wlocal=it->second;
                    break;
                }
            }
            if(wlocal==-1){//if found
                cout<<"Not found edge e("<<ID1<<","<<ID2<<") in overlay graph!"<<endl; exit(1);
            }
            woverlay=QueryCore(ID1,ID2);
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
                        cout<<"Already exists! "<<ID1<<" "<<ID2<<" "<<updatesSet[make_pair(ID1,ID2)].first<<" "<<updatesSet[make_pair(ID1,ID2)].second<<endl; exit(1);
                    }
//                    weightsParti.emplace_back(make_pair(ID1,ID2), make_pair(wlocal, woverlay));//oldW,newW
                }else{
                    if(updatesSet.find(make_pair(ID2,ID1))==updatesSet.end()){//if not found
                        updatesSet.insert({make_pair(ID2,ID1),make_pair(wlocal, woverlay)});
                    }else{
                        cout<<"Already exists! "<<ID2<<" "<<ID1<<" "<<updatesSet[make_pair(ID2,ID1)].first<<" "<<updatesSet[make_pair(ID2,ID1)].second<<endl; exit(1);
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
//        cout<<"Size of weightsParti (with new) of partition "<<pid<<" : "<<weightsParti.size()<<endl;
        bool ifPost=false;
        if(PSPStrategy>=PostBoundary){
            ifPost=true;
        }
        vector<pair<pair<int,int>,int>> updatedSC;
        IncreasePartiBatchForOpt(pid, weightsParti, NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid],SCconNodesMTP,VidtoTNidP,updatedSC, ifPost);
        ifRepaired[pid]=true;
    }else if(partiBatch.find(pid)!=partiBatch.end()){
//        cout<<"Size of weightsParti of partition "<<pid<<" : "<<weightsParti.size()<<endl;
//        if(pid==49){
//            cout<<pid<<endl;
//            for(auto it=ProBeginVertexSetParti[pid].begin();it!=ProBeginVertexSetParti[pid].end();++it){
//                cout<<*it<<"("<<NodeOrder[*it]<<")"<<endl;
//            }
//        }
        IncreasePartiBatchLabel(Trees[pid], ranks[pid], heightMaxs[pid], ProBeginVertexSetParti[pid], VidtoTNidP);
        ifRepaired[pid]=true;
    }

}


void Graph::ConstructPH2H_PartiV(vector<int> & P, vector<vector<Node>>& Trees, vector<vector<int>>& ranks, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs){
    int PID;
    for(int i=0;i<P.size();++i){
        PID=P[i];
        ConstructPH2H_Parti(PID, Trees, ranks, SCconNodesMTP, VidtoTNidP, heightMaxs, toRMQs, RMQIndexs);
    }
}

void Graph::ConstructPH2H_PartiVCH(vector<int> & P, vector<vector<Node>>& Trees, vector<vector<int>>& ranks, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs){
    int pid;
    for(int i=0;i<P.size();++i){
        pid=P[i];
        H2HCreateTree_Parti(pid, Trees[pid], ranks[pid], SCconNodesMTP, VidtoTNidP, heightMaxs);
//        ConstructPH2H_Parti(PID, Trees, ranks, SCconNodesMTP, VidtoTNidP, heightMaxs, toRMQs, RMQIndexs, ifLabelC);
    }
}

void Graph::ConstructPH2H_Parti(int pid, vector<vector<Node>>& Trees, vector<vector<int>>& ranks, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs){
    //Create tree for partition
    H2HCreateTree_Parti(pid, Trees[pid], ranks[pid], SCconNodesMTP, VidtoTNidP, heightMaxs);
    //Create LCA index
    makeRMQCoreP(pid, toRMQs, RMQIndexs, Trees);
    //Create labels for partition
    H2HCreateIndex_Parti(pid, Trees[pid], ranks[pid]);
}
//void Graph::ConstructPH2H_PartiNoLabel(int pid, vector<vector<Node>>& Trees, vector<vector<int>>& ranks, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs){
//    //Create tree for partition
//    H2HCreateTree_Parti(pid, Trees[pid], ranks[pid], SCconNodesMTP, VidtoTNidP, heightMaxs);
//}
void Graph::ConstructPH2H_PartiLabel(int pid, vector<vector<Node>>& Trees, vector<vector<int>>& ranks, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs){
    //Create LCA index
    makeRMQCoreP(pid, toRMQs, RMQIndexs, Trees);
    //Create labels for partition
    H2HCreateIndex_Parti(pid, Trees[pid], ranks[pid]);
}

//Function of Creating tree for partition
void Graph::H2HCreateTree_Parti(int pid, vector<Node>& TreeP, vector<int>& rankP, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs) {
//    cout<<"Partition "<<pid<<": boundary number is "<<BoundVertex[pid].size()<<endl;
    /// Contraction
    int degree;
    int ID, ID1, ID2;
    unordered_map<vertex,bool> existCore; existCore.clear();
    for(auto it=PartiVertex[pid].begin();it!=PartiVertex[pid].end();++it){
        existCore.insert({*it,true});
    }

//    unordered_map<vertex,int> rankP; rankP.clear();
    for(int id=PartiVertex[pid].size()-1;id>=0;--id){
        ID = PartiVertex[pid][id];
//        rankP.insert({ID,-1});
        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();
        for(auto it=E[ID].begin();it!=E[ID].end();it++){
            if(existCore.find(it->first)==existCore.end()){// not found
                cout<<"Wrong neighbor! "<<ID<<"("<<PartiTag[ID].first<<") "<< it->first<<"("<<PartiTag[it->first].first<<")"<<endl; exit(1);
            }
            else{// if found
                if(existCore[it->first]){
                    Neigh.emplace_back(*it);
                }else{
                    cout<<"Not in core!"<<it->first<<endl; exit(1);
                }
            }
        }
        NeighborCon[ID].assign(Neigh.begin(),Neigh.end());
//        cout<<ID<<" "<<NeighborCon[ID].size()<<" "<<PartiTag[ID].second<<endl;

        existCore[ID]=false;
        //delete the star
        for(int i=0;i<Neigh.size();i++){
            int y=Neigh[i].first;
            deleteECore(ID,y);//delete ID from y's adjacency list
        }
        //add all-pair neighbors
        for(int i=0;i<Neigh.size();i++){
            ID1=Neigh[i].first;
            for(int j=i+1;j<Neigh.size();j++){
                ID2=Neigh[j].first;
                insertECore(ID1,ID2,Neigh[i].second.first+Neigh[j].second.first);
                /// For TD update
                if(ID1<ID2){
                    SCconNodesMTP[ID1][ID2].emplace_back(ID,Neigh[i].second.first+Neigh[j].second.first);//supportive vertex, no direction, may contain the supportive vertices for shortcuts between interface vertices
                }
                else{
                    SCconNodesMTP[ID2][ID1].emplace_back(ID,Neigh[i].second.first+Neigh[j].second.first);
                }

            }
        }


    }

    /// Create Tree
    ID=PartiVertex[pid][0];
    Node root;//virtual root node
    if(NeighborCon[ID].empty()){
//        cout<<"There exist non-virtual root!"<<endl;
        root.uniqueVertex=ID;
    }else{
        cout<<"Wrong!"<<endl; exit(1);
    }
    root.height=1;
    TreeP.push_back(root);
//    cout<<"0 "<<ID<<" "<<IDMap[ID]<<endl;
    rankP[IDMap[ID]] = 0;
//    rankP[ID] = 0;
    int CCNum=1;
    for(int id=1;id<PartiVertex[pid].size();++id){
        ID = PartiVertex[pid][id];
//        cout<<id<<" "<<ID<<" "<<IDMap[ID]<<endl;
        int nn;
        if(existCore[ID]){
            cout<<"Wrong: should be out of core"<<endl; exit(1);
        }
        Node nod;
        nod.vert=NeighborCon[ID];//
        nod.uniqueVertex=ID;
        int pa=matchCoreParti(ID,NeighborCon[ID], rankP);

//        cout<<"pa "<<pa<<" "<<TreeP[pa].height<<endl;
        if(pa!=-1){
            TreeP[pa].ch.push_back(TreeP.size());
            nod.pa=pa;
            nod.height=TreeP[pa].height+1;
            /// for update
            nod.hdepth=TreeP[pa].height+1;
            for(int i=0;i<NeighborCon[ID].size();i++){//for the neighbors which have higher order
                nn=NeighborCon[ID][i].first;
                if(PartiTag[nn].first != pid){
                    cout<<"Wrong nn! "<<PartiTag[nn].first <<" "<< pid<<endl; exit(1);
                }
                VidtoTNidP[nn].emplace_back(TreeP.size());
//            if(rankP.find(nn)==rankP.end()){
//                cout<<"Not found "<<nn<<" in rank!"<<endl; exit(1);
//            }
                if(TreeP[rankP[IDMap[nn]]].hdepth<TreeP[pa].height+1){
                    TreeP[rankP[IDMap[nn]]].hdepth=TreeP[pa].height+1;
                }

            }
        }else{
            CCNum+=1;
//            cout<<"PID "<<pid<<": "<<ID<<" is a solitary vertex."<<endl;
            nod.pa=-1;
            nod.height=1;
            /// for update
            nod.hdepth=1;
            for(int i=0;i<NeighborCon[ID].size();i++){//for the neighbors which have higher order
                nn=NeighborCon[ID][i].first;
                if(PartiTag[nn].first != pid){
                    cout<<"Wrong nn! "<<PartiTag[nn].first <<" "<< pid<<endl; exit(1);
                }
                VidtoTNidP[nn].emplace_back(TreeP.size());
//            if(rankP.find(nn)==rankP.end()){
//                cout<<"Not found "<<nn<<" in rank!"<<endl; exit(1);
//            }
                if(TreeP[rankP[IDMap[nn]]].hdepth<TreeP[pa].height+1){
                    TreeP[rankP[IDMap[nn]]].hdepth=TreeP[pa].height+1;
                }

            }
        }

        if(nod.height>heightMaxs[pid]){
            heightMaxs[pid]=nod.height;
        }

//        if(rankP.find(ID)==rankP.end()){
//            cout<<"Not found "<<ID<<" in rank!"<<endl; exit(1);
//        }
        rankP[IDMap[ID]]=TreeP.size();//the position of tree, higher-order vertex has lower rank
        TreeP.push_back(nod);

    }

    if(CCNum>1){
        cout<<"Partition "<<pid<<" is not a connected component. "<<CCNum<<endl;
        if(algoChoice==H2H){
            cout<<"P-TD requires the partition to be fully connected!"<<endl; exit(1);
        }
    }
//    cout<<pid<<"'s tree node number: "<<TreeP.size()<<" ; tree height: "<< heightMaxs[pid]<<endl;
}
//Function of tree-label index construction for partition
void Graph::H2HCreateIndex_Parti(int pid, vector<Node>& TreeP, vector<int>& rankP){
//    cout<<"Computing Tree Label for partition "<<pid<<endl;
    //initialize
    vector<int> list; //list.clear();
    list.push_back(TreeP[0].uniqueVertex);
    TreeP[0].pos.clear();
    TreeP[0].pos.push_back(0);

    for(int i=0;i<TreeP[0].ch.size();i++){
        makeTreeIndexDFSP(TreeP[0].ch[i],list,TreeP, rankP);
    }

}

//Function of Creating tree for overlay graph
void Graph::H2HCreateTree_Overlay() {
    //for H2H update
    SCconNodesMT.assign(node_num, map<int, vector<pair<int,int>>>());
    NeighborCon.assign(node_num, vector<pair<int,pair<int,int>>>());
//    _DD_.assign(node_num,0); //_DD2_.assign(node_num,0);
    DD.assign(node_num,0); //DD2.assign(node_num,0);
    VidtoTNid.assign(node_num,vector<int>());
    //initialize E
    map<int,pair<int,int>> m;
    E.assign(node_num,m);
    for(int i=0;i<NeighborsOverlay.size();i++){
        if(!NeighborsOverlay[i].empty()){
            for(auto it=NeighborsOverlay[i].begin();it!=NeighborsOverlay[i].end();++it){
                E[i].insert(make_pair(it->first,make_pair(it->second,1)));
            }
        }
    }
    /// Contraction
    int degree;
    int ID, ID1, ID2;
    vector<bool> existCore(node_num,true);
    rank.assign(node_num,-1);

    bool flagAdd = true;
    for(int id=OverlayVertex.size()-1;id>=0;--id){
        ID = OverlayVertex[id];
//        cout<<ID<<" "<<NodeOrder[ID]<<endl;
        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();
        for(auto it=E[ID].begin();it!=E[ID].end();it++){
            if(existCore[it->first]){
                Neigh.emplace_back(*it);
            }else{
                cout<<"Not in core!"<<it->first<<endl; exit(1);
            }
        }
        NeighborCon[ID].assign(Neigh.begin(),Neigh.end());

        existCore[ID]=false;
        //delete the star
        for(int i=0;i<Neigh.size();i++){
            int y=Neigh[i].first;
            deleteECore(ID,y);//delete ID from y's adjacency list
        }

        if(Neigh.size()<=100){
            //single thread
            for(int i=0;i<Neigh.size();i++){
                ID1=Neigh[i].first;
                for(int j=i+1;j<Neigh.size();j++){
                    ID2=Neigh[j].first;
                    insertECore(ID1,ID2,Neigh[i].second.first+Neigh[j].second.first);
                    /// For TD update
                    if(ID1<ID2){
                        SCconNodesMT[ID1][ID2].emplace_back(ID,Neigh[i].second.first+Neigh[j].second.first);//supportive vertex, no direction, may contain the supportive vertices for shortcuts between interface vertices
                    }else{
                        SCconNodesMT[ID2][ID1].emplace_back(ID,Neigh[i].second.first+Neigh[j].second.first);
                    }
                }
            }
        }else{
//            cout<<"Multiple thread for contraction. "<<ID<<" "<<Neigh.size()<<endl;
            //multiple thread
            if(Neigh.size()>threadnum){
                int step=Neigh.size()/threadnum;
                boost::thread_group thread;
                for(int i=0;i<threadnum;i++){
                    pair<int,int> p;
                    p.first=i*step;
                    if(i==threadnum-1)
                        p.second=Neigh.size();
                    else
                        p.second=(i+1)*step;
                    thread.add_thread(new boost::thread(&Graph::NeighborComorder, this, boost::ref(Neigh), p, ID));
                }
                thread.join_all();
            }else{
                boost::thread_group thread;
                for(int i=0;i<Neigh.size();i++){
                    pair<int,int> p;
                    p.first=i; p.second=(i+1);
                    thread.add_thread(new boost::thread(&Graph::NeighborComorder, this, boost::ref(Neigh), p, ID));
                }
                thread.join_all();
            }
        }

        //add all-pair neighbors
//        for(int i=0;i<Neigh.size();i++){
//            ID1=Neigh[i].first;
//            for(int j=i+1;j<Neigh.size();j++){
//                ID2=Neigh[j].first;
//                insertECore(ID1,ID2,Neigh[i].second.first+Neigh[j].second.first);
//                /// For TD update
//                SCconNodesMT[ID1][ID2].emplace_back(ID,Neigh[i].second.first+Neigh[j].second.first);//supportive vertex, no direction, may contain the supportive vertices for shortcuts between interface vertices
//                SCconNodesMT[ID2][ID1].emplace_back(ID,Neigh[i].second.first+Neigh[j].second.first);
//            }
//        }


    }
//    cout<<"Flag 1"<<endl;
    /// Create Tree
    ID=OverlayVertex[0];
    Node root;//virtual root node
    if(NeighborCon[ID].empty()){
//        cout<<"There exist non-virtual root!"<<endl;
        root.uniqueVertex=ID;
    }else{
        cout<<"Wrong!"<<endl; exit(1);
    }
    root.height=1;
    Tree.push_back(root);
//    rankP[IDMap[ID]] = 0;
    rank[ID] = 0;
    heightMax=0;

    for(int id=1;id<OverlayVertex.size();++id){
        ID = OverlayVertex[id];
//        cout<<ID<<" "<<NodeOrder[ID]<<endl;
        int nn;
        if(existCore[ID]){
            cout<<"Wrong: should be out of core"<<endl; exit(1);
        }
        Node nod;
        nod.vert=NeighborCon[ID];//
        nod.uniqueVertex=ID;
        int pa=matchCore(ID,NeighborCon[ID], rank);

        //cout<<"pa "<<pa<<endl;

        Tree[pa].ch.push_back(Tree.size());
        nod.pa=pa;
        nod.height=Tree[pa].height+1;
        /// for update
        nod.hdepth=Tree[pa].height+1;
        for(int i=0;i<NeighborCon[ID].size();i++){//for the neighbors which have higher order
            nn=NeighborCon[ID][i].first;
            VidtoTNid[nn].emplace_back(Tree.size());
            if(Tree[rank[nn]].hdepth<Tree[pa].height+1){
                Tree[rank[nn]].hdepth=Tree[pa].height+1;
            }

        }
        if(nod.height>heightMax){
            heightMax=nod.height;
        }

        rank[ID]=Tree.size();//the position of tree, higher-order vertex has lower rank
        Tree.push_back(nod);

    }

    /// LCA index
    makeRMQCore();//build LCA index

    cout<<"Overlay graph! Tree node number: "<<Tree.size()<<" ; Tree height: "<< heightMax<<endl;
}

void Graph::NeighborComorderH2H(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x){
//    sm->wait();
    int ID1, w1;
    int ID2, w2;
    for(int k=p.first;k<p.second;k++){
        ID1=Neighvec[k].first;
        w1=Neighvec[k].second.first;
        for(int h=0;h<Neighvec.size();h++){
            ID2=Neighvec[h].first;
            w2=Neighvec[h].second.first;
            if(ID1==ID2){
                continue;
            }
            insertEMTorder(ID1, ID2, w1+w2);
            if(ID1<ID2){
                if(SCconNodesMT[ID1].find(ID2)==SCconNodesMT[ID1].end()){//if not found
                    SCconNodesMT[ID1].insert({ID2,vector<pair<int,int>>()});
                }
                SCconNodesMT[ID1][ID2].emplace_back(x,w1+w2);
            }

        }
    }
//    sm->notify();
}

void Graph::NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x){
//    sm->wait();
    int ID1, w1;
    int ID2, w2;
    for(int k=p.first;k<p.second;k++){
        ID1=Neighvec[k].first;
        w1=Neighvec[k].second.first;
        for(int h=0;h<Neighvec.size();h++){
            ID2=Neighvec[h].first;
            w2=Neighvec[h].second.first;
            if(ID1==ID2){
                continue;
            }
            insertECoreMT(ID1,ID2,w1+w2);
            /// For TD update
            if(ID1<ID2){
                SCconNodesMT[ID1][ID2].emplace_back(x,w1+w2);//supportive vertex, no direction, may contain the supportive vertices for shortcuts between interface vertices
            }
        }
    }
//    sm->notify();
}

//Function of tree-label index construction for partition
void Graph::H2HCreateIndex_Overlay(){
    //initialize
    vector<int> list; //list.clear();
    list.push_back(Tree[0].uniqueVertex);
    Tree[0].pos.clear();
    Tree[0].pos.push_back(0);

    for(int i=0;i<Tree[0].ch.size();i++){
        makeIndexDFS(Tree[0].ch[i],list);
    }

}



//function of computing the H2H label of peripheries: original version
void Graph::makeTreeIndexDFSP(int p, vector<int>& list,  vector<Node>& TreeP, vector<int>& rankP){
//initialize
//    cout<<"Map "<<p<<" "<<IDMap[p]<<endl;
//    p=IDMap[p];
    int NeiNum=TreeP[p].vert.size();
    TreeP[p].pos.assign(NeiNum+1,0);
    TreeP[p].dis.assign(list.size(),INF);
    TreeP[p].cnt.assign(list.size(),0);
    TreeP[p].FN.assign(list.size(),true);
    TreeP[p].vAncestor=list;

    //pos
    //map<int,Nei> Nmap; Nmap.clear();//shortcut infor ordered by the pos ID
    for(int i=0;i<NeiNum;i++){
        for(int j=0;j<list.size();j++){
            if(TreeP[p].vert[i].first==list[j]){
                TreeP[p].pos[i]=j;//record the position of neighbors
                TreeP[p].dis[j]=TreeP[p].vert[i].second.first;
                TreeP[p].cnt[j]=1;
                break;
            }
        }
    }
    TreeP[p].pos[NeiNum]=list.size();


    //dis
    for(int i=0;i<NeiNum;i++){
        int x=TreeP[p].vert[i].first;
        int disvb=TreeP[p].vert[i].second.first;
        int k=TreeP[p].pos[i];//the kth ancestor is x

        for(int j=0;j<list.size();j++){//check the distance to the j-th ancestor could be updated by neighbors, including the valley path and peak path
            int y=list[j];//the jth ancestor is y

            int z;//the distance from x to y
            if(k!=j){
                if(k<j)//x is the ancestor of y, peak path
                    z=TreeP[rankP[IDMap[y]]].dis[k];
                else if(k>j)//y is the ancestor of x, valley path
                    z=TreeP[rankP[IDMap[x]]].dis[j];

                if(TreeP[p].dis[j]>z+disvb){
                    TreeP[p].dis[j]=z+disvb;
                    TreeP[p].FN[j]=false;
                    TreeP[p].cnt[j]=1;
                }else if(TreeP[p].dis[j]==z+disvb){
                    TreeP[p].cnt[j]+=1;
                }
            }
        }
    }

    //nested loop
    list.push_back(TreeP[p].uniqueVertex);
    for(int i=0;i<TreeP[p].ch.size();i++){
        makeTreeIndexDFSP(TreeP[p].ch[i],list, TreeP, rankP);
    }
    list.pop_back();
}

void Graph::makeIndexDFS(int p, vector<int>& list){
    //initialize
    int NeiNum=Tree[p].vert.size();
    Tree[p].pos.assign(NeiNum+1,0);
    Tree[p].dis.assign(list.size(),INF);
    Tree[p].cnt.assign(list.size(),0);
    Tree[p].FN.assign(list.size(),true);
    Tree[p].vAncestor=list;

    //pos
    //map<int,Nei> Nmap; Nmap.clear();//shortcut infor ordered by the pos ID
    for(int i=0;i<NeiNum;i++){
        for(int j=0;j<list.size();j++){
            if(Tree[p].vert[i].first==list[j]){
                Tree[p].pos[i]=j;
                Tree[p].dis[j]=Tree[p].vert[i].second.first;
                Tree[p].cnt[j]=1;
                break;
            }
        }
    }
    Tree[p].pos[NeiNum]=list.size();

    //dis
    for(int i=0;i<NeiNum;i++){
        int x=Tree[p].vert[i].first;
        int disvb=Tree[p].vert[i].second.first;
        int k=Tree[p].pos[i];//the kth ancestor is x

        for(int j=0;j<list.size();j++){//check the distance to the j-th ancestor could be updated by neighbors, including the valley path and peak path
            int y=list[j];//the jth ancestor is y

            int z;//the distance from x to y
            if(k!=j){
                if(k<j)//x is the ancestor of y, peak path
                    z=Tree[rank[y]].dis[k];
                else if(k>j)//y is the ancestor of x, valley path
                    z=Tree[rank[x]].dis[j];

                if(Tree[p].dis[j]>z+disvb){
                    Tree[p].dis[j]=z+disvb;
                    Tree[p].FN[j]=false;
                    Tree[p].cnt[j]=1;
                }else if(Tree[p].dis[j]==z+disvb){
                    Tree[p].cnt[j]+=1;
                }
            }
        }
    }

    //nested loop
    list.push_back(Tree[p].uniqueVertex);
    for(int i=0;i<Tree[p].ch.size();i++){
        makeIndexDFS(Tree[p].ch[i],list);
    }
    list.pop_back();
}



//function of erasing edge (u,v), i.e., erase u from v's adjacency list.
void Graph::deleteECore(int u,int v){
//	if(Emap[u].find(v)!=Emap[u].end()){
//		Emap[u].erase(Emap[u].find(v));
//		DD[u]--;
//	}

    if(E[v].find(u)!=E[v].end()){
        E[v].erase(E[v].find(u));
        DD[v]--;
    }
}
//function of inserting edge (u,v)
void Graph::insertECore(int u,int v,int w){
    if(E[u].find(v)==E[u].end()){//if not found
        E[u].insert(make_pair(v,make_pair(w,1)));
        DD[u]++;
//		DD2[u]++;
    }
    else{//if found
        if(E[u][v].first>w)
            E[u][v]= make_pair(w,1);
        else if(E[u][v].first==w)
            E[u][v].second+=1;
    }

    if(E[v].find(u)==E[v].end()){
        E[v].insert(make_pair(u,make_pair(w,1)));
        DD[v]++;
//		DD2[v]++;
    }
    else{
        if(E[v][u].first>w)
            E[v][u]=make_pair(w,1);
        else if(E[v][u].first==w)
            E[v][u].second+=1;
    }
}
void Graph::insertECoreMT(int u,int v,int w){
    if(E[u].find(v)==E[u].end()){//if not found
        E[u].insert(make_pair(v,make_pair(w,1)));
        DD[u]++;
//		DD2[u]++;
    }
    else{//if found
        if(E[u][v].first>w)
            E[u][v]= make_pair(w,1);
        else if(E[u][v].first==w)
            E[u][v].second+=1;
    }

}

//compute the father tree node
int Graph::matchCore(int x,vector<pair<int,pair<int,int>>> &vert, vector<int>& rank){
    int nearest=vert[0].first;
    for(int i=1;i<vert.size();i++){
        if(NodeOrder[vert[i].first]<NodeOrder[nearest])//get the least node order
            nearest=vert[i].first;
    }
    //cout<<nearest<<" "<<rankCore[nearest]<<endl;
    return rank[nearest];
}


//compute the father tree node
int Graph::matchCoreParti(int x,vector<pair<int,pair<int,int>>> &vert, vector<int>& rank){
    if(vert.empty()){//deal with solitary vertex
        return -1;
    }
    int nearest=vert[0].first;
    for(int i=1;i<vert.size();i++){
        if(NodeOrder[vert[i].first]<NodeOrder[nearest])//get the least node order
            nearest=vert[i].first;
    }
    //cout<<nearest<<" "<<rankCore[nearest]<<endl;
//    if(rank.find(nearest)==rank.end()){//if not found
//        cout<<"Not found "<<nearest<<" in rank!"<<endl; exit(1);
//    }
    return rank[IDMap[nearest]];
}

//construct RMQ index
void Graph::makeRMQCoreP(int pid, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs, vector<vector<Node>>& Trees){
    //EulerSeq.clear();
    toRMQs[pid].assign(PartiVertex[pid].size(),0);
    vector<int> EulerSeqP;
    //RMQIndex.clear();
    makeRMQDFSCoreP(pid, 0, 1, EulerSeqP, toRMQs, Trees);
    RMQIndexs[pid].push_back(EulerSeqP);

    int m = EulerSeqP.size();
//    cout<<"m: "<<m<<endl;
    for (int i = 2, k = 1; i < m; i = i * 2, k++){
        vector<int> tmp;
        //tmp.clear();
        tmp.assign(m,0);
        for (int j = 0; j < m - i; j++){
            int x = RMQIndexs[pid][k - 1][j], y = RMQIndexs[pid][k - 1][j + i / 2];
//            cout<<"x and y: "<<x<<" "<<y<<endl;
            if (Trees[pid][x].height < Trees[pid][y].height)
                tmp[j] = x;
            else tmp[j] = y;
        }
        RMQIndexs[pid].push_back(tmp);
    }
}

void Graph::makeRMQDFSCoreP(int pid, int p, int height, vector<int>& EulerSeqP, vector<vector<int>>& toRMQs, vector<vector<Node>>& Trees){
    toRMQs[pid][p] = EulerSeqP.size();//toRMQ[p] records the first occurrence of node p in the Euler tour, i.e., R[p]
    EulerSeqP.push_back(p);//EulerSeq is the Euler tour, i.e., E[1,...,2n-1]
    for (int i = 0; i < Trees[pid][p].ch.size(); i++){
        makeRMQDFSCoreP(pid,Trees[pid][p].ch[i], height + 1, EulerSeqP, toRMQs, Trees);
        EulerSeqP.push_back(p);
    }
}

void Graph::makeRMQ(vector<int>& toRMQ, vector<vector<int>>& RMQIndex, vector<Node>& Tree){
    vector<int> EulerSeq;
    EulerSeq.clear();
    toRMQ.assign(node_num,0);
    //RMQIndex.clear();
    makeRMQDFS(0, 1, EulerSeq, toRMQ, Tree);
    RMQIndex.push_back(EulerSeq);

    int m = EulerSeq.size();
    for (int i = 2, k = 1; i < m; i = i * 2, k++){
        vector<int> tmp;
        //tmp.clear();
        tmp.assign(m,0);
        for (int j = 0; j < m - i; j++){
            int x = RMQIndex[k - 1][j], y = RMQIndex[k - 1][j + i / 2];
            if (Tree[x].height < Tree[y].height)
                tmp[j] = x;
            else tmp[j] = y;
        }
        RMQIndex.push_back(tmp);
    }
}

void Graph::makeRMQDFS(int p, int height, vector<int>& EulerSeq, vector<int>& toRMQ, vector<Node>& Tree){
    toRMQ[p] = EulerSeq.size();//toRMQ[p] records the first occurrence of node p in the Euler tour, i.e., R[p]
    EulerSeq.push_back(p);//EulerSeq is the Euler tour, i.e., E[1,...,2n-1]
    for (int i = 0; i < Tree[p].ch.size(); i++){
        makeRMQDFS(Tree[p].ch[i], height + 1, EulerSeq, toRMQ, Tree);
        EulerSeq.push_back(p);
    }
}

//construct RMQ index
void Graph::makeRMQCore(){
    vector<int> EulerSeq;
    EulerSeq.clear();
    toRMQ.assign(node_num,0);
    //RMQIndex.clear();
    makeRMQDFSCore(0, 1, EulerSeq);
    RMQIndex.push_back(EulerSeq);

    int m = EulerSeq.size();
    for (int i = 2, k = 1; i < m; i = i * 2, k++){
        vector<int> tmp;
        //tmp.clear();
        tmp.assign(m,0);
        for (int j = 0; j < m - i; j++){
            int x = RMQIndex[k - 1][j], y = RMQIndex[k - 1][j + i / 2];
            if (Tree[x].height < Tree[y].height)
                tmp[j] = x;
            else tmp[j] = y;
        }
        RMQIndex.push_back(tmp);
    }
}

void Graph::makeRMQDFSCore(int p, int height, vector<int>& EulerSeq){
    toRMQ[p] = EulerSeq.size();//toRMQ[p] records the first occurrence of node p in the Euler tour, i.e., R[p]
    EulerSeq.push_back(p);//EulerSeq is the Euler tour, i.e., E[1,...,2n-1]
    for (int i = 0; i < Tree[p].ch.size(); i++){
        makeRMQDFSCore(Tree[p].ch[i], height + 1, EulerSeq);
        EulerSeq.push_back(p);
    }
}


/// Query Processing

int Graph::LCAQuery(int _p, int _q){
    int p = toRMQ[_p], q = toRMQ[_q];
    if (p > q){
        int x = p;
        p = q;
        q = x;
    }
    int len = q - p + 1;
    int i = 1, k = 0;
    while (i * 2 < len){
        i *= 2;
        k++;
    }
    q = q - i + 1;
    if (Tree[RMQIndex[k][p]].height < Tree[RMQIndex[k][q]].height)
        return RMQIndex[k][p];
    else return RMQIndex[k][q];
}


//function for Query processing, new
int Graph::QueryPCH(int ID1, int ID2, vector<vector<Node>>& Trees){
    int dis=INF;

    if(PartiTag[ID1].first==PartiTag[ID2].first){//if same partition
        if(PSPStrategy==PreBoundary || PSPStrategy==PostBoundary){
//            cout<<"Same-parti"<<endl;
            dis= QueryCHPartition(ID1,ID2,PartiTag[ID1].first,Trees);
        }else if(PSPStrategy==NoBoundary){
            if(PartiTag[ID1].second && PartiTag[ID2].second){//Case 1: both in overlay graph
//                cout<<"Same-parti: Core-Core"<<endl;
                dis=QueryCoreCH(ID1, ID2);
            }
            else if(PartiTag[ID1].second && !PartiTag[ID2].second){//Case 2: ID2 in partition, ID1 in core
//                cout<<"Same-parti: Core-Parti"<<endl;
                dis=QueryPartiCoreCH(ID2, ID1);
            }else if(!PartiTag[ID1].second && PartiTag[ID2].second){//Case 2: ID1 in partition, ID2 in core
//                cout<<"Same-parti: Parti-Core"<<endl;
                dis = QueryPartiCoreCH(ID1, ID2);
            }else{//Case 3: Same partition
//                cout<<"Same-parti: same-parti"<<endl;
                dis= QueryPartiPartiCH(ID1,ID2);
            }
        }
    }
    else{//if different partition
        if(PartiTag[ID1].second && PartiTag[ID2].second){//Case 1: both in overlay graph
//            cout<<"Different Parti: Core-Core"<<endl;
            dis=QueryCoreCH(ID1, ID2);
        }
        else if(PartiTag[ID1].second && !PartiTag[ID2].second){//Case 2: ID2 in partition, ID1 in core
//            cout<<"Different Parti: Core-Parti"<<endl;
            dis=QueryPartiCoreCH(ID2, ID1);
        }else if(!PartiTag[ID1].second && PartiTag[ID2].second){//Case 2: ID1 in partition, ID2 in core
//            cout<<"Different Parti: Parti-Core"<<endl;
            dis=QueryPartiCoreCH(ID1, ID2);
        }
        else{//Case 4: Different Partitions
//            cout<<"Different Parti: Parti-Parti"<<endl;
            dis = QueryPartiPartiCH(ID1, ID2);
        }

    }

    return dis;
}

//function for Query processing, new
int Graph::QueryPCHDebug(int ID1, int ID2, vector<vector<Node>>& Trees){
    int dis=INF;

    if(PartiTag[ID1].first==PartiTag[ID2].first){//if same partition
        if(PSPStrategy==PreBoundary || PSPStrategy==PostBoundary){
            cout<<"Same-parti"<<endl;
            dis= QueryCHPartition(ID1,ID2,PartiTag[ID1].first,Trees);
        }else if(PSPStrategy==NoBoundary){
            if(PartiTag[ID1].second && PartiTag[ID2].second){//Case 1: both in overlay graph
                cout<<"Same-parti: Core-Core"<<endl;
                dis=QueryCoreCH(ID1, ID2);
            }
            else if(PartiTag[ID1].second && !PartiTag[ID2].second){//Case 2: ID2 in partition, ID1 in core
                cout<<"Same-parti: Core-Parti"<<endl;
                dis=QueryPartiCoreCH(ID2, ID1);
            }else if(!PartiTag[ID1].second && PartiTag[ID2].second){//Case 2: ID1 in partition, ID2 in core
                cout<<"Same-parti: Parti-Core"<<endl;
                dis = QueryPartiCoreCH(ID1, ID2);
            }else{//Case 3: Same partition
                cout<<"Same-parti: same-parti"<<endl;
                dis= QueryPartiPartiCH(ID1,ID2);
            }
        }
    }
    else{//if different partition
        if(PartiTag[ID1].second && PartiTag[ID2].second){//Case 1: both in overlay graph
            cout<<"Different Parti: Core-Core"<<endl;
            dis=QueryCoreCH(ID1, ID2);
        }
        else if(PartiTag[ID1].second && !PartiTag[ID2].second){//Case 2: ID2 in partition, ID1 in core
            cout<<"Different Parti: Core-Parti"<<endl;
            dis=QueryPartiCoreCH(ID2, ID1);
        }else if(!PartiTag[ID1].second && PartiTag[ID2].second){//Case 2: ID1 in partition, ID2 in core
            cout<<"Different Parti: Parti-Core"<<endl;
            dis=QueryPartiCoreCH(ID1, ID2);
        }
        else{//Case 4: Different Partitions
            cout<<"Different Parti: Parti-Parti"<<endl;
            dis = QueryPartiPartiCH(ID1, ID2);
        }

    }

    return dis;
}

//function for Query processing, new
int Graph::QueryPH2H(int ID1, int ID2){
    int dis=INF;


    if(PartiTag[ID1].first==PartiTag[ID2].first){//if same partition
        if(PSPStrategy==PreBoundary){
//            cout<<"Same-parti"<<endl;
            dis= QuerySamePartiPre(ID1,ID2);
        }else if(PSPStrategy==PostBoundary){
//            cout<<"Same-parti"<<endl;
            dis= QuerySamePartiPost(ID1,ID2);
        }else if(PSPStrategy==NoBoundary){//no-boundary
            if(PartiTag[ID1].second && PartiTag[ID2].second){//Case 1: both in overlay graph
//                cout<<"Same-parti: Core-Core"<<endl;
                dis=QueryCore(ID1, ID2);
            }
            else if(PartiTag[ID1].second && !PartiTag[ID2].second){//Case 2: ID2 in partition, ID1 in core
//                cout<<"Same-parti: Core-Parti"<<endl;
                dis=QueryPartiCore(ID2, ID1);
            }else if(!PartiTag[ID1].second && PartiTag[ID2].second){//Case 2: ID1 in partition, ID2 in core
//                cout<<"Same-parti: Parti-Core"<<endl;
                dis = QueryPartiCore(ID1, ID2);
            }else{//Case 3: Same partition
//                cout<<"Same-parti: Non-boundary"<<endl;
                dis= QuerySameParti(ID1,ID2);
            }
        }
    }
    else{//if different partition
        if(PartiTag[ID1].second && PartiTag[ID2].second){//Case 1: both in overlay graph
//            cout<<"Different Parti: Core-Core"<<endl;
            dis=QueryCore(ID1, ID2);
        }
        else if(PartiTag[ID1].second && !PartiTag[ID2].second){//Case 2: ID2 in partition, ID1 in core
//            cout<<"Different Parti: Core-Parti"<<endl;
            dis=QueryPartiCore(ID2, ID1);
        }else if(!PartiTag[ID1].second && PartiTag[ID2].second){//Case 2: ID1 in partition, ID2 in core
//            cout<<"Different Parti: Parti-Core"<<endl;
            dis = QueryPartiCore(ID1, ID2);
        }
        else{//Case 4: Different Partitions
//            cout<<"Different Parti: Parti-Parti"<<endl;
            dis = QueryPartiParti(ID1, ID2);
        }
    }

    return dis;
}



//Case 1: query on overlay graph
int Graph::QueryCore(int ID1, int ID2){
    if(!PartiTag[ID1].second || !PartiTag[ID2].second){
        cout<<"Not overlay vertex! "<<ID1<<"("<<PartiTag[ID1].second<<") "<<ID2<<"("<<PartiTag[ID2].second<<")"<<endl; exit(1);
    }
    if(ID1==ID2) return 0;
    if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
    int r1=rank[ID1], r2=rank[ID2];
    int LCA=LCAQueryOverlay(r1,r2);

    if(LCA==r1)
        return Tree[r2].dis[Tree[r1].pos.back()];
    else if(LCA==r2)
        return Tree[r1].dis[Tree[r2].pos.back()];
    else{
        int tmp=INF;
        for(int i=0;i<Tree[LCA].pos.size();i++){
//            if(Tree[LCA].pos[i]>=Tree[r1].dis.size() || Tree[LCA].pos[i]>=Tree[r2].dis.size()){
//                cout<<ID1<<"("<<r1<<","<<PartiTag[ID1].second<<") "<<ID2<<"("<<r2<<","<<PartiTag[ID2].second<<") "<<LCA<<" "<<Tree.size()<<": "<<Tree[LCA].pos[i]<<" "<<Tree[r1].dis.size()<<" "<<Tree[r2].dis.size()<<endl;
//            }
            if(tmp>Tree[r1].dis[Tree[LCA].pos[i]]+Tree[r2].dis[Tree[LCA].pos[i]])
                tmp=Tree[r1].dis[Tree[LCA].pos[i]]+Tree[r2].dis[Tree[LCA].pos[i]];
        }
        return tmp;
    }
}
int Graph::QueryCoreDebug(int ID1, int ID2){
    if(!PartiTag[ID1].second || !PartiTag[ID2].second){
        cout<<"Not overlay vertex! "<<ID1<<"("<<PartiTag[ID1].second<<") "<<ID2<<"("<<PartiTag[ID2].second<<")"<<endl; exit(1);
    }
    if(ID1==ID2) return 0;
    if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
    int r1=rank[ID1], r2=rank[ID2];
    int LCA=LCAQueryOverlay(r1,r2);
    int d1,d2,ancestor1,ancestor2;
    int dis=INF;

    cout<<"Overlay Query."<<endl;

    if(LCA==r1){
        cout<<"r1. Hub of QueryCore: "<<Tree[LCA].uniqueVertex<<endl;
        dis=Tree[r2].dis[Tree[r1].pos.back()];
        cout<<dis<<" "<<DijkstraCore(ID1,ID2)<<" "<<Dijkstra(ID1,ID2,Neighbor)<<endl;
        DijkstraPath(ID1,ID2,NeighborsOverlay);
        DijkstraPath(ID1,ID2,Neighbor);
    }
    else if(LCA==r2){
        cout<<"r2. Hub of QueryCore: "<<Tree[LCA].uniqueVertex<<endl;
        dis=Tree[r1].dis[Tree[r2].pos.back()];
        cout<<dis<<" "<<DijkstraCore(ID1,ID2)<<" "<<Dijkstra(ID1,ID2,Neighbor)<<endl;
        DijkstraPath(ID1,ID2,NeighborsOverlay);
        DijkstraPath(ID1,ID2,Neighbor);
    }
    else{
        cout<<"LCA="<<LCA<<" "<<Tree[LCA].uniqueVertex<<endl;
        for(int i=0;i<Tree[LCA].pos.size();i++){
            if(dis>Tree[r1].dis[Tree[LCA].pos[i]]+Tree[r2].dis[Tree[LCA].pos[i]]){
                d1=Tree[r1].dis[Tree[LCA].pos[i]], d2=Tree[r2].dis[Tree[LCA].pos[i]];
                ancestor1=Tree[r1].vAncestor[Tree[LCA].pos[i]], ancestor2=Tree[r2].vAncestor[Tree[LCA].pos[i]];
                dis=d1+d2;
                cout<<ancestor1<<"("<<NodeOrder[ancestor1]<<") "<<ancestor2<<"("<<NodeOrder[ancestor2]<<") "<<dis<<endl;

            }
        }
        cout<<"Hub of QueryCore: "<<ancestor1<<"("<<NodeOrder[ancestor1]<<") "<<ancestor2<<"("<<NodeOrder[ancestor2]<<")"<< endl;
        DijkstraPath(ID1,ID2,NeighborsOverlay);
        DijkstraPath(ID1,ID2,Neighbor);
    }
    return dis;
}

//Case 2: one core, one tree
int Graph::QueryPartiCore(int ID1, int ID2){//ID1 partition, ID2 core
    int d=INF;

    int pid=PartiTag[ID1].first;
    int bid;
    int dis1,dis2;
    if(PSPStrategy==NoBoundary || PSPStrategy==PreBoundary){
        for(auto it=BoundVertex[pid].begin();it!=BoundVertex[pid].end();++it){
            bid=*it;
            dis1= QueryH2HPartition(ID1,bid,pid);
            dis2= QueryCore(bid,ID2);
//            cout<<ID1<<" "<<bid<<" "<<ID2<<": "<<dis1<<" "<<dis2<<" "<<dis1+dis2<<" "<<d<<endl;
            if(d>dis1+dis2)
                d=dis1+dis2;
        }
    }
    else if(PSPStrategy==PostBoundary){
        for(auto it=BoundVertex[pid].begin();it!=BoundVertex[pid].end();++it){
            bid=*it;
            dis1= QueryH2HPartitionPost(ID1,bid,pid);
            dis2= QueryCore(bid,ID2);
            if(d>dis1+dis2)
                d=dis1+dis2;
        }
    }

    return d;
}

int Graph::QueryPartiCoreDebug(int ID1, int ID2){//ID1 partition, ID2 core
    int d=INF;

    int pid=PartiTag[ID1].first;
    int bid;
    int dis1,dis2;
    int finalbid,finaldis1,finaldis2;
    if(PSPStrategy==NoBoundary || PSPStrategy==PreBoundary){
        for(auto it=BoundVertex[pid].begin();it!=BoundVertex[pid].end();++it){
            bid=*it;
//            dis1= QueryH2HPartitionDebug(ID1,bid,pid);
            dis1= QueryH2HPartition(ID1,bid,pid);
            dis2= QueryCore(bid,ID2);
            cout<<bid<<": "<<dis1<<" "<<dis2<<" "<<d<<endl;
            if(d>dis1+dis2){
                d=dis1+dis2;
//                cout<<bid<<": "<<dis1<<" "<<dis2<<" "<<d<<endl;
                finalbid=bid, finaldis1=dis1, finaldis2=dis2;
            }

        }
    }
    else if(PSPStrategy==PostBoundary){
        for(auto it=BoundVertex[pid].begin();it!=BoundVertex[pid].end();++it){
            bid=*it;
            dis1= QueryH2HPartitionPost(ID1,bid,pid);
            dis2= QueryCore(bid,ID2);
            if(d>dis1+dis2){
                d=dis1+dis2;
                cout<<bid<<": "<<dis1<<" "<<dis2<<" "<<d<<endl;
                finalbid=bid, finaldis1=dis1, finaldis2=dis2;
            }
        }
    }

    int dDijk_s=Dijkstra(ID1,finalbid,Neighbor), dDijk_t=Dijkstra(finalbid,ID2,Neighbor);
    cout<<ID1<<" "<<finalbid<<"("<<NodeOrder[finalbid]<<","<<PartiTag[finalbid].first<<","<<PartiTag[finalbid].second<<") "<<ID2<<" : "<<finaldis1<<" "<<finaldis2<<" "<<d<<" ; "<<dDijk_s<<" "<<dDijk_t<<" "<<Dijkstra(ID1,ID2,Neighbor)<<endl;
    if(finaldis2==0){
        QueryH2HPartitionDebug(ID1,ID2,pid);
    }
    DijkstraPath(ID1,ID2,Neighbor);
    return d;
}

//Case 3: Different trees
int Graph::QueryPartiParti(int ID1, int ID2){//both are within the partition
    int d=INF;

    int pid1=PartiTag[ID1].first;
    int pid2=PartiTag[ID2].first;
    if(pid1==pid2){//if in the same partition
        cout<<"Wrong for partition-partition query!"<<endl;
        exit(1);

    }else{//if in different partitions
//        cout<<"Parti-Parti: "<<pid1<<" "<<pid2<<endl;
//        vector<int> B1=BoundVertex[pid1];
//        vector<int> B2=BoundVertex[pid2];

        map<int,int> m1,m2;
        m1.clear();
        m2.clear();
        int bID1, bID2, tempdis;
        int b1,b2,d1,d2;

//        qBNum=qBNum+BoundVertex[pid1].size()+BoundVertex[pid2].size();
        if(PSPStrategy==NoBoundary || PSPStrategy==PreBoundary){
            for(int i=0;i<BoundVertex[pid1].size();i++){
                bID1=BoundVertex[pid1][i];
                m1.insert(make_pair(bID1, QueryH2HPartition(ID1,bID1,pid1)));
            }
            for(int j=0;j<BoundVertex[pid2].size();j++){
                bID2=BoundVertex[pid2][j];
                m2.insert(make_pair(bID2,QueryH2HPartition(ID2,bID2,pid2)));
            }
        }else if(PSPStrategy==PostBoundary){
            for(int i=0;i<BoundVertex[pid1].size();i++){
                bID1=BoundVertex[pid1][i];
                m1.insert(make_pair(bID1, QueryH2HPartitionPost(ID1,bID1,pid1)));
            }
            for(int j=0;j<BoundVertex[pid2].size();j++){
                bID2=BoundVertex[pid2][j];
                m2.insert(make_pair(bID2,QueryH2HPartitionPost(ID2,bID2,pid2)));
            }
        }

        qBNum=qBNum+BoundVertex[pid1].size()*BoundVertex[pid2].size();
        for(int k=0;k<BoundVertex[pid1].size();k++){
            bID1=BoundVertex[pid1][k];

            if(m1[bID1]>d)
                continue;

            for(int z=0;z<BoundVertex[pid2].size();z++){
                bID2=BoundVertex[pid2][z];

                if(m2[bID2]>d)
                    continue;

//                qBNum++;
                tempdis=m1[bID1]+QueryCore(bID1,bID2)+m2[bID2];
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


//Case 4: Same tree, for original version
int Graph::QuerySameParti(int ID1, int ID2){//both are within the partition
    int d=INF;

    int pid1=PartiTag[ID1].first;
    int pid2=PartiTag[ID2].first;
    if(pid1==pid2){//if in the same partition
//        cout<<"Same-Parti"<<endl;
        int temp_dis = QueryH2HPartition(ID1,ID2,pid1);/// d2 may be wrong sometimes
        if(temp_dis<d)//QueryH2HPartition(ID1,ID2,pid1)
            d=temp_dis;//QueryH2HPartition(ID1,ID2,pid1);
//        vector<int> B=BoundVertex[pid1];
        map<int,int> m1,m2;
        m1.clear();
        m2.clear();
        vector<int> B1,B2;
        B1.clear();
        B2.clear();
        int bID,d1,d2;
        for(int i=0;i<BoundVertex[pid1].size();i++){
            bID=BoundVertex[pid1][i];
            d1=QueryH2HPartition(ID1,bID,pid1);
            d2=QueryH2HPartition(ID2,bID,pid1);

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
                    tempdis=m1[bID1]+QueryCore(bID1,bID2)+m2[bID2];
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

//Case 4: Same tree, for query-orient version
int Graph::QuerySamePartiPost(int ID1, int ID2){//both are within the partition
    int d=INF;

    int pid1=PartiTag[ID1].first;
    int pid2=PartiTag[ID2].first;
    if(pid1==pid2){//if in the same partition
//        cout<<"Same-Parti"<<endl;
        d = QueryH2HPartitionPost(ID1,ID2,pid1);

    }else{//if in different partitions
        cout<<"Wrong for same partition query!"<<endl;
        exit(1);
    }

    return d;
}

//Case 4: Same tree, for query-orient version
int Graph::QuerySamePartiPre(int ID1, int ID2){//both are within the partition
    int d=INF;

    int pid1=PartiTag[ID1].first;
    int pid2=PartiTag[ID2].first;
    if(pid1==pid2){//if in the same partition
//        cout<<"Same-Parti"<<endl;
        d = QueryH2HPartition(ID1,ID2,pid1);

    }else{//if in different partitions
        cout<<"Wrong for same partition query!"<<endl;
        exit(1);
    }

    return d;
}

int Graph::QuerySamePartiPostOpt(int ID1, int ID2){//both are within the partition
    int d=INF;

    int pid1=PartiTag[ID1].first;
    int pid2=PartiTag[ID2].first;
    if(pid1==pid2){//if in the same partition
//        cout<<"Same-Parti"<<endl;
        d = QueryH2HPartition(ID1,ID2,pid1);

    }else{//if in different partitions
        cout<<"Wrong for same partition query!"<<endl;
        exit(1);
    }

    return d;
}

int Graph::QueryCoreCH(int ID1, int ID2){
    if(ID1==ID2) return 0;
    int d=INF;
    benchmark::heap<2,int,int> fHeapForward(node_num);
    benchmark::heap<2, int, int> fHeapBackward(node_num);

    //closed or not
    vector<bool> vVisitedF(node_num, false);
    vector<bool> vVisitedB(node_num, false);
    //the existing shortest distance
    vector<int>	vDistanceForward(node_num, INF);
    vector<int>	vDistanceBackward(node_num, INF);
    //stop search or not
    bool bF = false;
    bool bB = false;
    vDistanceForward[ID1] = 0;
    vDistanceBackward[ID2] = 0;
    fHeapForward.update(ID1,0);
    fHeapBackward.update(ID2,0);
    int rForward, rBackward;

    int topNodeIDForward, topNodeIDBackward,topDisForward,topDisBackward, neighborNodeID, neighborLength;

    while(!fHeapForward.empty() || !fHeapBackward.empty() )
    {
        if(bF && bB)
            break;
        if(bF && fHeapBackward.empty())
            break;
        if(bB && fHeapForward.empty())
            break;
        //Forward Search
        if(!fHeapForward.empty() && !bF)
        {
            fHeapForward.extract_min(topNodeIDForward, topDisForward);
            //cout<<topNodeIDForward<<" "<<topDisForward<<" "<<(NodeOrder[topNodeIDForward]>NodeOrder[137099])<<endl;

            if(vDistanceForward[topNodeIDForward] > d)
                bF = true;

            vVisitedF[topNodeIDForward] = true;

            if(vVisitedB[topNodeIDForward]){
                int distTmp=topDisForward+vDistanceBackward[topNodeIDForward];
                if(distTmp<d){
                    d=distTmp;
                    //cout<<"forwardtopvalue "<<topDisForward<<" "<<vDistanceBackward[topNodeIDForward]<<" "<<d<<"meet "<<topNodeIDForward<<endl;
                }
            }

//            for(auto out=NeighborCon[topNodeIDForward].begin();out!=NeighborCon[topNodeIDForward].end();out++){
            rForward=rank[topNodeIDForward];
            for(auto out=Tree[rForward].vert.begin();out!=Tree[rForward].vert.end();out++){
                neighborNodeID = (*out).first;
                neighborLength = (*out).second.first;

                int df = vDistanceForward[topNodeIDForward] + neighborLength;
                if(!vVisitedF[neighborNodeID]){
                    if(vDistanceForward[neighborNodeID] > df){
                        //if(neighborNodeID==37783) cout<<"neighdis "<<vDistanceForward[neighborNodeID]<<" adddis "<<df<<endl;
                        vDistanceForward[neighborNodeID] = df;
                        fHeapForward.update(neighborNodeID, df);
                    }
                }
            }
        }

        //Backward Search
        if(!fHeapBackward.empty() && !bB)
        {
            fHeapBackward.extract_min(topNodeIDBackward, topDisBackward);

            if(vDistanceBackward[topNodeIDBackward] > d)
                bB = true;

            vVisitedB[topNodeIDBackward] = true;

            if(vVisitedF[topNodeIDBackward]){
                int distTmp=topDisBackward+vDistanceForward[topNodeIDBackward];
                if(distTmp<d){
                    d=distTmp;
                    //cout<<"backtopvalue "<<topDisBackward<<" "<<vDistanceForward[topNodeIDBackward]<<" "<<d<<"meet "<<topNodeIDBackward<<endl;
                }
            }

//            for(auto in=NeighborCon[topNodeIDBackward].begin();in!=NeighborCon[topNodeIDBackward].end();in++){
            rBackward=rank[topNodeIDBackward];
            for(auto in=Tree[rBackward].vert.begin();in!=Tree[rBackward].vert.end();in++){
                neighborNodeID = (*in).first;
                neighborLength = (*in).second.first;

                int db = vDistanceBackward[topNodeIDBackward] + neighborLength;
                if(!vVisitedB[neighborNodeID]){
                    if(vDistanceBackward[neighborNodeID]>db){
                        vDistanceBackward[neighborNodeID] = db;
                        fHeapBackward.update(neighborNodeID, db);
                    }
                }
            }
        }
    }
    return d;
}

//Query within one partition, PCH
int Graph::QueryCHPartition(int ID1, int ID2, int PID, vector<vector<Node>>& Trees){
    if(ID1==ID2) return 0;
    int d=INF;
    benchmark::heap<2,int,int> fHeapForward(node_num);
    benchmark::heap<2, int, int> fHeapBackward(node_num);

    //closed or not
    vector<bool> vVisitedF(node_num, false);
    vector<bool> vVisitedB(node_num, false);
    //the existing shortest distance
    vector<int>	vDistanceForward(node_num, INF);
    vector<int>	vDistanceBackward(node_num, INF);
    //stop search or not
    bool bF = false;
    bool bB = false;
    vDistanceForward[ID1] = 0;
    vDistanceBackward[ID2] = 0;
    fHeapForward.update(ID1,0);
    fHeapBackward.update(ID2,0);
    int rForward, rBackward;

    int topNodeIDForward, topNodeIDBackward,topDisForward,topDisBackward, neighborNodeID, neighborLength;

    while(!fHeapForward.empty() || !fHeapBackward.empty() )
    {
        if(bF && bB)
            break;
        if(bF && fHeapBackward.empty())
            break;
        if(bB && fHeapForward.empty())
            break;
        //Forward Search
        if(!fHeapForward.empty() && !bF)
        {
            fHeapForward.extract_min(topNodeIDForward, topDisForward);
            //cout<<topNodeIDForward<<" "<<topDisForward<<" "<<(NodeOrder[topNodeIDForward]>NodeOrder[137099])<<endl;

            if(vDistanceForward[topNodeIDForward] > d)
                bF = true;

            vVisitedF[topNodeIDForward] = true;

            if(vVisitedB[topNodeIDForward]){
                int distTmp=topDisForward+vDistanceBackward[topNodeIDForward];
                if(distTmp<d){
                    d=distTmp;
                    //cout<<"forwardtopvalue "<<topDisForward<<" "<<vDistanceBackward[topNodeIDForward]<<" "<<d<<"meet "<<topNodeIDForward<<endl;
                }
            }

//            for(auto out=NeighborCon[topNodeIDForward].begin();out!=NeighborCon[topNodeIDForward].end();out++){
            rForward=ranks[PID][IDMap[topNodeIDForward]];
            for(auto out=Trees[PID][rForward].vert.begin();out!=Trees[PID][rForward].vert.end();out++){
                neighborNodeID = (*out).first;
                neighborLength = (*out).second.first;

                int df = vDistanceForward[topNodeIDForward] + neighborLength;
                if(!vVisitedF[neighborNodeID]){
                    if(vDistanceForward[neighborNodeID] > df){
                        //if(neighborNodeID==37783) cout<<"neighdis "<<vDistanceForward[neighborNodeID]<<" adddis "<<df<<endl;
                        vDistanceForward[neighborNodeID] = df;
                        fHeapForward.update(neighborNodeID, df);
                    }
                }
            }
        }

        //Backward Search
        if(!fHeapBackward.empty() && !bB)
        {
            fHeapBackward.extract_min(topNodeIDBackward, topDisBackward);

            if(vDistanceBackward[topNodeIDBackward] > d)
                bB = true;

            vVisitedB[topNodeIDBackward] = true;

            if(vVisitedF[topNodeIDBackward]){
                int distTmp=topDisBackward+vDistanceForward[topNodeIDBackward];
                if(distTmp<d){
                    d=distTmp;
                    //cout<<"backtopvalue "<<topDisBackward<<" "<<vDistanceForward[topNodeIDBackward]<<" "<<d<<"meet "<<topNodeIDBackward<<endl;
                }
            }

//            for(auto in=NeighborCon[topNodeIDBackward].begin();in!=NeighborCon[topNodeIDBackward].end();in++){
            rBackward=ranks[PID][IDMap[topNodeIDBackward]];
            for(auto in=Trees[PID][rBackward].vert.begin();in!=Trees[PID][rBackward].vert.end();in++){
                neighborNodeID = (*in).first;
                neighborLength = (*in).second.first;

                int db = vDistanceBackward[topNodeIDBackward] + neighborLength;
                if(!vVisitedB[neighborNodeID]){
                    if(vDistanceBackward[neighborNodeID]>db){
                        vDistanceBackward[neighborNodeID] = db;
                        fHeapBackward.update(neighborNodeID, db);
                    }
                }
            }
        }
    }
    return d;
}


//Query from partition to core, no-boundary of PCH
int Graph::QueryPartiCoreCH(int ID1, int ID2){//ID1: partition vertex
    if(ID1==ID2) return 0;
    if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
    int d=INF;
    benchmark::heap<2,int,int> fHeapForward(node_num);
    benchmark::heap<2, int, int> fHeapBackward(node_num);

    //closed or not
    vector<bool> vVisitedF(node_num, false);
    vector<bool> vVisitedB(node_num, false);
    //the existing shortest distance
    vector<int>	vDistanceForward(node_num, INF);
    vector<int>	vDistanceBackward(node_num, INF);
    //stop search or not
    bool bF = false;
    bool bB = false;
    vDistanceForward[ID1] = 0;
    vDistanceBackward[ID2] = 0;
    fHeapForward.update(ID1,0);
    fHeapBackward.update(ID2,0);

    int rForward, rBackward;
    int PID;
    int topNodeIDForward, topNodeIDBackward,topDisForward,topDisBackward, neighborNodeID, neighborLength;

    while(!fHeapForward.empty() || !fHeapBackward.empty() )
    {
        if(bF && bB)
            break;
        if(bF && fHeapBackward.empty())
            break;
        if(bB && fHeapForward.empty())
            break;
        //Forward Search
        if(!fHeapForward.empty() && !bF)
        {
            fHeapForward.extract_min(topNodeIDForward, topDisForward);
            //cout<<topNodeIDForward<<" "<<topDisForward<<" "<<(NodeOrder[topNodeIDForward]>NodeOrder[137099])<<endl;

            if(vDistanceForward[topNodeIDForward] > d)
                bF = true;

            vVisitedF[topNodeIDForward] = true;

            if(vVisitedB[topNodeIDForward]){
                int distTmp=topDisForward+vDistanceBackward[topNodeIDForward];
                if(distTmp<d){
                    d=distTmp;
                    //cout<<"forwardtopvalue "<<topDisForward<<" "<<vDistanceBackward[topNodeIDForward]<<" "<<d<<"meet "<<topNodeIDForward<<endl;
                }
            }

            //if(VtoParID[topNodeIDForward]==PID1){
//            for(auto out=NeighborCons[PID1][topNodeIDForward].begin();out!=NeighborCons[PID1][topNodeIDForward].end();out++){
            PID=PartiTag[topNodeIDForward].first;
            rForward=ranks[PID][IDMap[topNodeIDForward]];
            for(auto out=Trees[PID][rForward].vert.begin();out!=Trees[PID][rForward].vert.end();out++){
                neighborNodeID = (*out).first;
                neighborLength = (*out).second.first;

                int df = vDistanceForward[topNodeIDForward] + neighborLength;
                if(!vVisitedF[neighborNodeID]){
                    if(vDistanceForward[neighborNodeID] > df){
                        //if(neighborNodeID==37783) cout<<"neighdis "<<vDistanceForward[neighborNodeID]<<" adddis "<<df<<endl;
                        vDistanceForward[neighborNodeID] = df;
                        fHeapForward.update(neighborNodeID, df);
                    }
                }
            }
            //}else{
//            for(auto out=NeighborConOverlay[topNodeIDForward].begin();out!=NeighborConOverlay[topNodeIDForward].end();out++){
            if(PartiTag[topNodeIDForward].second){
                rForward=rank[topNodeIDForward];
                if(rForward!=-1){
                    for(auto out=Tree[rForward].vert.begin();out!=Tree[rForward].vert.end();out++){
                        neighborNodeID = (*out).first;
                        neighborLength = (*out).second.first;

                        int df = vDistanceForward[topNodeIDForward] + neighborLength;
                        if(!vVisitedF[neighborNodeID]){
                            if(vDistanceForward[neighborNodeID] > df){
                                //if(neighborNodeID==37783) cout<<"neighdis "<<vDistanceForward[neighborNodeID]<<" adddis "<<df<<endl;
                                vDistanceForward[neighborNodeID] = df;
                                fHeapForward.update(neighborNodeID, df);
                            }
                        }
                    }
                }
            }


            //}

        }

        //Backward Search
        if(!fHeapBackward.empty() && !bB)
        {
            fHeapBackward.extract_min(topNodeIDBackward, topDisBackward);

            if(vDistanceBackward[topNodeIDBackward] > d)
                bB = true;

            vVisitedB[topNodeIDBackward] = true;

            if(vVisitedF[topNodeIDBackward]){
                int distTmp=topDisBackward+vDistanceForward[topNodeIDBackward];
                if(distTmp<d){
                    d=distTmp;
                    //cout<<"backtopvalue "<<topDisBackward<<" "<<vDistanceForward[topNodeIDBackward]<<" "<<d<<"meet "<<topNodeIDBackward<<endl;
                }
            }

//            for(auto in=NeighborConOverlay[topNodeIDBackward].begin();in!=NeighborConOverlay[topNodeIDBackward].end();in++){
            rBackward=rank[topNodeIDBackward];
            for(auto in=Tree[rBackward].vert.begin();in!=Tree[rBackward].vert.end();in++){
                neighborNodeID = (*in).first;
                neighborLength = (*in).second.first;

                int db = vDistanceBackward[topNodeIDBackward] + neighborLength;
                if(!vVisitedB[neighborNodeID]){
                    if(vDistanceBackward[neighborNodeID]>db){
                        vDistanceBackward[neighborNodeID] = db;
                        fHeapBackward.update(neighborNodeID, db);
                    }
                }
            }
        }
    }
    return d;
}

//Query from partition to core, no-boundary of PCH
int Graph::QueryPartiPartiCH(int ID1, int ID2) {//ID1 and ID2 partition vertices
    if(ID1==ID2) return 0;
    if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
    int d=INF;
    benchmark::heap<2,int,int> fHeapForward(node_num);
    benchmark::heap<2, int, int> fHeapBackward(node_num);

    //closed or not
    vector<bool> vVisitedF(node_num, false);
    vector<bool> vVisitedB(node_num, false);
    //the existing shortest distance
    vector<int>	vDistanceForward(node_num, INF);
    vector<int>	vDistanceBackward(node_num, INF);
    //stop search or not
    bool bF = false;
    bool bB = false;
    vDistanceForward[ID1] = 0;
    vDistanceBackward[ID2] = 0;
    fHeapForward.update(ID1,0);
    fHeapBackward.update(ID2,0);

    int rForward, rBackward;
    int PID;
    int topNodeIDForward, topNodeIDBackward,topDisForward,topDisBackward, neighborNodeID, neighborLength;

    while(!fHeapForward.empty() || !fHeapBackward.empty() )
    {
        if(bF && bB)
            break;
        if(bF && fHeapBackward.empty())
            break;
        if(bB && fHeapForward.empty())
            break;
        //Forward Search
        if(!fHeapForward.empty() && !bF)
        {
            fHeapForward.extract_min(topNodeIDForward, topDisForward);
            //cout<<topNodeIDForward<<" "<<topDisForward<<" "<<(NodeOrder[topNodeIDForward]>NodeOrder[137099])<<endl;

            if(vDistanceForward[topNodeIDForward] > d)
                bF = true;

            vVisitedF[topNodeIDForward] = true;

            if(vVisitedB[topNodeIDForward]){
                int distTmp=topDisForward+vDistanceBackward[topNodeIDForward];
                if(distTmp<d){
                    d=distTmp;
                    //cout<<"forwardtopvalue "<<topDisForward<<" "<<vDistanceBackward[topNodeIDForward]<<" "<<d<<"meet "<<topNodeIDForward<<endl;
                }
            }

            //if(VtoParID[topNodeIDForward]==PID1){
//            for(auto out=NeighborCons[PID1][topNodeIDForward].begin();out!=NeighborCons[PID1][topNodeIDForward].end();out++){
            PID=PartiTag[topNodeIDForward].first;
            rForward=ranks[PID][IDMap[topNodeIDForward]];
            for(auto out=Trees[PID][rForward].vert.begin();out!=Trees[PID][rForward].vert.end();out++){
                neighborNodeID = (*out).first;
                neighborLength = (*out).second.first;

                int df = vDistanceForward[topNodeIDForward] + neighborLength;
                if(!vVisitedF[neighborNodeID]){
                    if(vDistanceForward[neighborNodeID] > df){
                        //if(neighborNodeID==37783) cout<<"neighdis "<<vDistanceForward[neighborNodeID]<<" adddis "<<df<<endl;
                        vDistanceForward[neighborNodeID] = df;
                        fHeapForward.update(neighborNodeID, df);
                    }
                }
            }
            //}else{
//            for(auto out=NeighborConOverlay[topNodeIDForward].begin();out!=NeighborConOverlay[topNodeIDForward].end();out++){
            if(PartiTag[topNodeIDForward].second){
                rForward=rank[topNodeIDForward];
                if(rForward!=-1){
                    for(auto out=Tree[rForward].vert.begin();out!=Tree[rForward].vert.end();out++){
                        neighborNodeID = (*out).first;
                        neighborLength = (*out).second.first;

                        int df = vDistanceForward[topNodeIDForward] + neighborLength;
                        if(!vVisitedF[neighborNodeID]){
                            if(vDistanceForward[neighborNodeID] > df){
                                //if(neighborNodeID==37783) cout<<"neighdis "<<vDistanceForward[neighborNodeID]<<" adddis "<<df<<endl;
                                vDistanceForward[neighborNodeID] = df;
                                fHeapForward.update(neighborNodeID, df);
                            }
                        }
                    }
                }
            }


            //}

        }

        //Backward Search
        if(!fHeapBackward.empty() && !bB)
        {
            fHeapBackward.extract_min(topNodeIDBackward, topDisBackward);

            if(vDistanceBackward[topNodeIDBackward] > d)
                bB = true;

            vVisitedB[topNodeIDBackward] = true;

            if(vVisitedF[topNodeIDBackward]){
                int distTmp=topDisBackward+vDistanceForward[topNodeIDBackward];
                if(distTmp<d){
                    d=distTmp;
                    //cout<<"backtopvalue "<<topDisBackward<<" "<<vDistanceForward[topNodeIDBackward]<<" "<<d<<"meet "<<topNodeIDBackward<<endl;
                }
            }

//            for(auto in=NeighborConOverlay[topNodeIDBackward].begin();in!=NeighborConOverlay[topNodeIDBackward].end();in++){
            PID=PartiTag[topNodeIDBackward].first;
            rBackward=ranks[PID][IDMap[topNodeIDBackward]];
            for(auto in=Trees[PID][rBackward].vert.begin();in!=Trees[PID][rBackward].vert.end();in++){
                neighborNodeID = (*in).first;
                neighborLength = (*in).second.first;

                int db = vDistanceBackward[topNodeIDBackward] + neighborLength;
                if(!vVisitedB[neighborNodeID]){
                    if(vDistanceBackward[neighborNodeID]>db){
                        vDistanceBackward[neighborNodeID] = db;
                        fHeapBackward.update(neighborNodeID, db);
                    }
                }
            }
            if(PartiTag[topNodeIDBackward].second){
                rBackward=rank[topNodeIDBackward];
                if(rBackward!=-1) {
                    for (auto in = Tree[rBackward].vert.begin(); in != Tree[rBackward].vert.end(); in++) {
                        neighborNodeID = (*in).first;
                        neighborLength = (*in).second.first;

                        int db = vDistanceBackward[topNodeIDBackward] + neighborLength;
                        if (!vVisitedB[neighborNodeID]) {
                            if (vDistanceBackward[neighborNodeID] > db) {
                                vDistanceBackward[neighborNodeID] = db;
                                fHeapBackward.update(neighborNodeID, db);
                            }
                        }
                    }
                }
            }

        }
    }
    return d;
}

//Query within one partition, no-boundary
int Graph::QueryH2HPartition(int ID1, int ID2, int PID){
    if(ID1==ID2) return 0;
    if(PartiTag[ID1].first!=PID || PartiTag[ID2].first!=PID){
        cout<<"Wrong! ID1 and ID2 are not in the same partition! "<<PartiTag[ID1].first<<" "<<PartiTag[ID2].first<<" "<<PID<<endl; exit(1);
    }
    int dis=INF;
    int dis1, dis2, dis1Final, dis2Final, hub;
    int r1=ranks[PID][IDMap[ID1]], r2=ranks[PID][IDMap[ID2]];
//    cout<<"ID1: "<<ID1<<" "<<IDMap[ID1]<<" "<<r1<<" ; ID2: "<<ID2<<" "<<IDMap[ID2]<<" "<<r2<<endl;
    int LCA=LCAQueryPartition(r1,r2,PID);
//    cout<<"LCA: "<<LCA<<endl;
    if(LCA==r1){
        dis = Trees[PID][r2].dis[Trees[PID][r1].pos.back()];
//        cout<<"r1 "<<dis<<" "<<Dijkstra(ID1,ID2,Neighbor)<<endl;
    }
    else if(LCA==r2){
        dis = Trees[PID][r1].dis[Trees[PID][r2].pos.back()];
//        cout<<"r2 "<<dis<<" "<<Dijkstra(ID1,ID2,Neighbor)<<endl;
    }
    else{

        for(int i=0;i<Trees[PID][LCA].pos.size();i++){
            dis1=Trees[PID][r1].dis[Trees[PID][LCA].pos[i]], dis2=Trees[PID][r2].dis[Trees[PID][LCA].pos[i]];
            if(dis>dis1+dis2){
                dis=dis1+dis2;
                dis1Final=dis1, dis2Final=dis2;
                hub=Trees[PID][r1].vAncestor[Trees[PID][LCA].pos[i]];
            }

        }
//        int dDijk_s=Dijkstra(ID1,hub,Neighbor), dDijk_t=Dijkstra(hub,ID2,Neighbor);
//        cout<<ID1<<" "<<hub<<"("<<NodeOrder[hub]<<","<<PartiTag[hub].first<<","<<PartiTag[hub].second<<") "<<ID2<<" : "<<dis1Final<<" "<<dis2Final<<" "<<dis<<" ; "<<dDijk_s<<" "<<dDijk_t<<" "<<Dijkstra(ID1,ID2,Neighbor)<<endl;
    }
    return dis;
}

int Graph::QueryH2HPartitionDebug(int ID1, int ID2, int PID){
    if(ID1==ID2) return 0;
    if(PartiTag[ID1].first!=PID || PartiTag[ID2].first!=PID){
        cout<<"Wrong! ID1 and ID2 are not in the same partition! "<<PartiTag[ID1].first<<" "<<PartiTag[ID2].first<<" "<<PID<<endl; exit(1);
    }
    int dis=INF;
    int dis1, dis2, dis1Final, dis2Final, hub;
    int r1=ranks[PID][IDMap[ID1]], r2=ranks[PID][IDMap[ID2]];
//    cout<<"ID1: "<<ID1<<" "<<IDMap[ID1]<<" "<<r1<<" ; ID2: "<<ID2<<" "<<IDMap[ID2]<<" "<<r2<<endl;
    int LCA=LCAQueryPartition(r1,r2,PID);
//    cout<<"LCA: "<<LCA<<endl;
    bool flag=false;
//    if(ID1==10 && ID2==768){
//        cout<<"Find. "<<ID1<<" "<<ID2<<endl;
        flag=true;
//    }
    if(LCA==r1){
        dis = Trees[PID][r2].dis[Trees[PID][r1].pos.back()];
        if(flag)
            cout<<"r1 "<<dis<<" "<<Dijkstra(ID1,ID2,Neighbor)<<endl;
    }
    else if(LCA==r2){
        dis = Trees[PID][r1].dis[Trees[PID][r2].pos.back()];
        if(flag)
            cout<<"r2 "<<dis<<" "<<Dijkstra(ID1,ID2,Neighbor)<<endl;
    }
    else{

        for(int i=0;i<Trees[PID][LCA].pos.size();i++){
            dis1=Trees[PID][r1].dis[Trees[PID][LCA].pos[i]], dis2=Trees[PID][r2].dis[Trees[PID][LCA].pos[i]];
            if(dis>dis1+dis2){
                dis=dis1+dis2;
                dis1Final=dis1, dis2Final=dis2;
                hub=Trees[PID][r1].vAncestor[Trees[PID][LCA].pos[i]];
            }

        }
        int dDijk_s=Dijkstra(ID1,hub,Neighbor), dDijk_t=Dijkstra(hub,ID2,Neighbor);
        cout<<ID1<<" "<<hub<<"("<<NodeOrder[hub]<<","<<PartiTag[hub].first<<","<<PartiTag[hub].second<<") "<<ID2<<" : "<<dis1Final<<" "<<dis2Final<<" "<<dis<<" ; "<<dDijk_s<<" "<<dDijk_t<<" "<<Dijkstra(ID1,ID2,Neighbor)<<endl;
//        DijkstraPath(ID1,ID2,Neighbor);
    }
    return dis;
}

//Query within one partition, no-boundary
int Graph::QueryH2HPartitionPost(int ID1, int ID2, int PID){
    if(ID1==ID2) return 0;
    if(PartiTag[ID1].first!=PID || PartiTag[ID2].first!=PID){
        cout<<"Wrong! ID1 and ID2 are not in the same partition! "<<PartiTag[ID1].first<<" "<<PartiTag[ID2].first<<" "<<PID<<endl; exit(1);
    }
    int r1=ranksPost[PID][IDMap[ID1]], r2=ranksPost[PID][IDMap[ID2]];
//    cout<<"ID1: "<<ID1<<" "<<IDMap[ID1]<<" "<<r1<<" ; ID2: "<<ID2<<" "<<IDMap[ID2]<<" "<<r2<<endl;
    int LCA=LCAQueryPartitionPost(r1,r2,PID);
//    cout<<"LCA: "<<LCA<<endl;
    if(LCA==r1)
        return TreesPost[PID][r2].dis[TreesPost[PID][r1].pos.back()];
    else if(LCA==r2)
        return TreesPost[PID][r1].dis[TreesPost[PID][r2].pos.back()];
    else{
        int tmp=INF;
        for(int i=0;i<TreesPost[PID][LCA].pos.size();i++){
            if(tmp>TreesPost[PID][r1].dis[TreesPost[PID][LCA].pos[i]]+TreesPost[PID][r2].dis[TreesPost[PID][LCA].pos[i]])
                tmp=TreesPost[PID][r1].dis[TreesPost[PID][LCA].pos[i]]+TreesPost[PID][r2].dis[TreesPost[PID][LCA].pos[i]];
        }
        return tmp;
    }
}

int Graph::LCAQueryPartition(int _p, int _q, int PID){
    int p = toRMQs[PID][_p], q = toRMQs[PID][_q];
    if (p > q){
        int x = p;
        p = q;
        q = x;
    }
    int len = q - p + 1;
    int i = 1, k = 0;
    while (i * 2 < len){
        i *= 2;
        k++;
    }
    q = q - i + 1;
    if (Trees[PID][RMQIndexs[PID][k][p]].height < Trees[PID][RMQIndexs[PID][k][q]].height)
        return RMQIndexs[PID][k][p];
    else return RMQIndexs[PID][k][q];
}

int Graph::LCAQueryPartitionPost(int _p, int _q, int PID){
    int p = toRMQsPost[PID][_p], q = toRMQsPost[PID][_q];
    if (p > q){
        int x = p;
        p = q;
        q = x;
    }
    int len = q - p + 1;
    int i = 1, k = 0;
    while (i * 2 < len){
        i *= 2;
        k++;
    }
    q = q - i + 1;
    if (TreesPost[PID][RMQIndexsPost[PID][k][p]].height < TreesPost[PID][RMQIndexsPost[PID][k][q]].height)
        return RMQIndexsPost[PID][k][p];
    else return RMQIndexsPost[PID][k][q];
}

int Graph::LCAQueryOverlay(int _p, int _q){
    int p = toRMQ[_p], q = toRMQ[_q];
    if (p > q){
        int x = p;
        p = q;
        q = x;
    }
    int len = q - p + 1;
    int i = 1, k = 0;
    while (i * 2 < len){
        i *= 2;
        k++;
    }
    q = q - i + 1;
    if (Tree[RMQIndex[k][p]].height < Tree[RMQIndex[k][q]].height)
        return RMQIndex[k][p];
    else return RMQIndex[k][q];
}

int Graph::LCAQuery(int _p, int _q, vector<int>& toRMQ, vector<vector<int>>& RMQIndex, vector<Node>& Tree){
    int p = toRMQ[_p], q = toRMQ[_q];
    if (p > q){
        int x = p;
        p = q;
        q = x;
    }
    int len = q - p + 1;
    int i = 1, k = 0;
    while (i * 2 < len){
        i *= 2;
        k++;
    }
    q = q - i + 1;
    if (Tree[RMQIndex[k][p]].height < Tree[RMQIndex[k][q]].height)
        return RMQIndex[k][p];
    else return RMQIndex[k][q];
}


/// Index maintenance
//H2H index update
void Graph::DecreaseOverlay(int a,int b, int newW, vector<unordered_map<vertex,int>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax){
//    map<int,int> checkedDis;//map<tree node ID, distance index>

    if(Neighbors[a].find(b)!=Neighbors[a].end()){
        Neighbors[a][b]=newW;
    }else{
        cout<<"Wrong for Neighbors!"<<endl; exit(1);
    }
    if(Neighbors[b].find(a)!=Neighbors[b].end()){
        Neighbors[b][a]=newW;
    }else{
        cout<<"Wrong for Neighbors!"<<endl; exit(1);
    }

    for(int i=0;i<Tree.size();i++){
        Tree[i].DisRe.clear();
    }

    int lid,hid;
    if(NodeOrder[a]<NodeOrder[b]){
        lid=a;hid=b;
    }else{
        lid=b;hid=a;
    }

    int IniH=Tree[rank[lid]].height;//the height where weight change begins
    int ProH=Tree[rank[lid]].height; int ProID=lid;
    vector<set<int>> SCre;//record the shortcut change in each height
    set<int> ss;//ss.clear();
    SCre.assign(ProH+1,ss);

    //map<int,set<int>> DisRe;//rankid; record the distance change caused by the shortcut in each height
    //DisRe.clear();

    int MinH;

    bool tri=false;
    for(int i=0;i<Tree[rank[lid]].vert.size();i++){
        if(Tree[rank[lid]].vert[i].first==hid){
            if(Tree[rank[lid]].vert[i].second.first>newW){
                Tree[rank[lid]].vert[i].second.first=newW;
                Tree[rank[lid]].vert[i].second.second=1;
                tri=true;
                SCre[ProH].insert(hid);
                MinH=IniH;
            }else if(Tree[rank[lid]].vert[i].second.first==newW){
                Tree[rank[lid]].vert[i].second.second+=1;
            }
            break;
        }
    }

    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed

    vector<int> ProIDRecord; ProIDRecord.assign(ProH+1,0);

    //int ProBeginH;
    int ProBeginID;
    if(tri){
        //cout<<"Bottom-up ;;;;;;;;;;;;;;;;;; "<<endl;
        while(ProH>=MinH){

            ProIDRecord[ProH]=ProID;
            vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
            bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
            for(auto it=SCre[ProH].begin();it!=SCre[ProH].end();it++){
                int Cid=*it; int Cw;//=OCdis[make_pair(ProID,Cid)];
                int cidH=Tree[rank[Cid]].height-1;

                map<int,int> Hnei; //Hnei.clear();
                vector<pair<int,int>> Lnei; //Lnei.clear();
                for(int j=0;j<Vert.size();j++){
                    if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                        Hnei[Vert[j].first]=Vert[j].second.first;
                    }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                        Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                    }else{
                        Cw=Vert[j].second.first;
                    }
                }

                if(Tree[rank[ProID]].dis[cidH]>=Cw){
                    Tree[rank[ProID]].dis[cidH]=Cw;
                    Tree[rank[ProID]].FN[cidH]=true;
                    ProIDdisCha=true;
                    Tree[rank[ProID]].DisRe.insert(Cid);
                    //DisRe[rank[ProID]].insert(Cid); //cout<<"dischange Cid "<<Cid<<endl;
                }

                int hid,hidHeight,lid,lidHeight,wsum;
                for(int j=0;j<Tree[rank[Cid]].vert.size();j++){
                    hid=Tree[rank[Cid]].vert[j].first;hidHeight=Tree[rank[hid]].height-1;
                    if(Hnei.find(hid)!=Hnei.end()){
                        wsum=Cw+Hnei[hid];
                        if(wsum<Tree[rank[Cid]].vert[j].second.first){
                            Tree[rank[Cid]].vert[j].second.first=wsum;
                            Tree[rank[Cid]].vert[j].second.second=1;
                            SCre[Tree[rank[Cid]].height].insert(hid);
                            if(Tree[rank[Cid]].height<MinH) MinH=Tree[rank[Cid]].height;

                        }else if(wsum==Tree[rank[Cid]].vert[j].second.first){
                            Tree[rank[Cid]].vert[j].second.second+=1;
                        }

                    }
                }
                for(int j=0;j<Lnei.size();j++){
                    lid=Lnei[j].first;lidHeight=Tree[rank[lid]].height-1;
                    for(int k=0;k<Tree[rank[lid]].vert.size();k++){
                        if(Tree[rank[lid]].vert[k].first==Cid){
                            wsum=Cw+Lnei[j].second;
                            if(Tree[rank[lid]].vert[k].second.first>wsum){
                                Tree[rank[lid]].vert[k].second.first=wsum;
                                Tree[rank[lid]].vert[k].second.second=1;
                                SCre[Tree[rank[lid]].height].insert(Cid);
                                if(Tree[rank[lid]].height<MinH) MinH=Tree[rank[lid]].height;

                            }else if(Tree[rank[lid]].vert[k].second.first==wsum){
                                Tree[rank[lid]].vert[k].second.second+=1;
                            }

                            break;
                        }
                    }
                }
            }

            if(ProIDdisCha){//if the distance labeling is detected changed
                vertexIDChL.insert(ProID);
                //ProBeginH=ProH;
                ProBeginID=ProID;
            }

            ProH-=1;
            ProID=Tree[Tree[rank[ProID]].pa].uniqueVertex;
        }

        vector<int> linee; //linee.clear();
        linee.reserve(heightMax);
        int pachidd=Tree[Tree[rank[ProBeginID]].pa].uniqueVertex;
        while(Tree[rank[pachidd]].height>1){
            linee.insert(linee.begin(),pachidd);
            pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
        }
        linee.insert(linee.begin(),pachidd);

        //top-down process
        EachNodeProBDis5(rank[ProBeginID], linee, vertexIDChL,  Tree, rank);
    }
    //return checkedDis.size();
}

void Graph::EachNodeProBDis5(int child,vector<int>& line,set<int>& vertexIDChL, vector<Node> &Tree, vector<int> &rank){
    bool ProIDdisCha=false;
    vector<int> cntNew(line.size(),0);
    vector<bool> flagUpdate(line.size(),false);
    if(Tree[child].DisRe.size()!=0){
        for(int k=0;k<Tree[child].vert.size();k++){
            int b=Tree[child].vert[k].first, bH=Tree[rank[b]].height-1,vbW=Tree[child].vert[k].second.first;
            if(Tree[child].FN[bH]){
                if(Tree[child].DisRe.find(b)!=Tree[child].DisRe.end()){//all ancestor check
                    for(int i=0;i<bH;i++){
//                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
                            Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
                            Tree[child].FN[i]=false;
                            Tree[child].cnt[i]=1;//new
                            ProIDdisCha=true;
                            flagUpdate[i]=true;
                            cntNew[i]=1;
                        }
                        else if(Tree[child].dis[i]==vbW+Tree[rank[b]].dis[i]){
                            cntNew[i]++;
                            if(flagUpdate[i]) {
                                Tree[child].cnt[i]+=1;//new
                            }
                            else if(cntNew[i]>Tree[child].cnt[i]){
                                Tree[child].cnt[i]=cntNew[i];//use cntNew to redress the cnt value since the edge decrease may lead to more ways for dis (i.e., increase the cnt)
                            }
                        }
                    }
                    for(int i=bH+1;i<line.size();i++){
//                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
                            Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
                            Tree[child].FN[i]=false;
                            Tree[child].cnt[i]=1;//new
                            ProIDdisCha=true;
                            flagUpdate[i]=true;
                            cntNew[i]=1;
                        }
                        else if(Tree[child].dis[i]==vbW+Tree[rank[line[i]]].dis[bH]){
                            cntNew[i]++;
                            if(flagUpdate[i]) {
                                Tree[child].cnt[i]+=1;//new
                            }
                            else if(cntNew[i]>Tree[child].cnt[i]){
                                Tree[child].cnt[i]=cntNew[i];
                            }
                        }
                    }

                }else{//partial ancestor check

                    if(vertexIDChL.find(b)!=vertexIDChL.end()){
                        for(int i=0;i<bH;i++){
//                            checkedDis.insert(make_pair(child,i));
                            if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
                                Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
                                Tree[child].FN[i]=false;
                                Tree[child].cnt[i]=1;//new
                                ProIDdisCha=true;
                                flagUpdate[i]=true;
                                cntNew[i]=1;
                            }
                            else if(Tree[child].dis[i]==vbW+Tree[rank[b]].dis[i]){
                                cntNew[i]++;
                                if(flagUpdate[i]) {
                                    Tree[child].cnt[i]+=1;//new
                                }
                                else if(cntNew[i]>Tree[child].cnt[i]){
                                    Tree[child].cnt[i]=cntNew[i];
                                }
                            }
                        }
                    }
                    for(int i=bH+1;i<line.size();i++){
//                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
                            Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
                            Tree[child].FN[i]=false;
                            Tree[child].cnt[i]=1;//new
                            ProIDdisCha=true;
                            flagUpdate[i]=true;
                            cntNew[i]=1;
                        }
                        else if(Tree[child].dis[i]==vbW+Tree[rank[line[i]]].dis[bH]){
                            cntNew[i]++;
                            if(flagUpdate[i]) {
                                Tree[child].cnt[i]+=1;//new
                            }
                            else if(cntNew[i]>Tree[child].cnt[i]){
                                Tree[child].cnt[i]=cntNew[i];
                            }
                        }
                    }

                }
            }
        }
    }else{
        for(int k=0;k<Tree[child].vert.size();k++){
            int b=Tree[child].vert[k].first, bH=Tree[rank[b]].height-1,vbW=Tree[child].vert[k].second.first;
            if(Tree[child].FN[bH]){
                if(vertexIDChL.find(b)!=vertexIDChL.end()){
                    for(int i=0;i<bH;i++){
//                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
                            Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
                            Tree[child].FN[i]=false;
                            Tree[child].cnt[i]=1;//new
                            ProIDdisCha=true;
                            flagUpdate[i]=true;
                            cntNew[i]=1;
                        }
                        else if(Tree[child].dis[i]==vbW+Tree[rank[b]].dis[i]){
                            cntNew[i]++;
                            if(flagUpdate[i]) {
                                Tree[child].cnt[i]+=1;//new
                            }
                            else if(cntNew[i]>Tree[child].cnt[i]){
                                Tree[child].cnt[i]=cntNew[i];
                            }
                        }
                    }
                }
                for(int i=bH+1;i<line.size();i++){
//                    checkedDis.insert(make_pair(child,i));
                    if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
                        Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
                        Tree[child].FN[i]=false;
                        Tree[child].cnt[i]=1;//new
                        ProIDdisCha=true;
                        flagUpdate[i]=true;
                        cntNew[i]=1;
                    }
                    else if(Tree[child].dis[i]==vbW+Tree[rank[line[i]]].dis[bH]){
                        cntNew[i]++;
                        if(flagUpdate[i]) {
                            Tree[child].cnt[i]+=1;//new
                        }
                        else if(cntNew[i]>Tree[child].cnt[i]){
                            Tree[child].cnt[i]=cntNew[i];
                        }
                    }
                }
            }
        }
    }

    if(ProIDdisCha){
//        sm->wait();
        vertexIDChL.insert(Tree[child].uniqueVertex);
//        vUpdated[Tree[child].uniqueVertex] = true;
//        sm->notify();
    }

    line.push_back(Tree[child].uniqueVertex);
    for(int i=0;i<Tree[child].ch.size();i++){
        EachNodeProBDis5(Tree[child].ch[i], line, vertexIDChL,Tree, rank);
    }
    line.pop_back();

}


void Graph::DecreaseParti(int a,int b, int newW, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax){
    map<int,int> checkedDis;//map<tree node ID, distance index>

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

    for(int i=0;i<Tree.size();i++){
        Tree[i].DisRe.clear();
    }

    int lid,hid;
    if(NodeOrder[a]<NodeOrder[b]){
        lid=a;hid=b;
    }else{
        lid=b;hid=a;
    }
    int lidM=IDMap[lid];
    int IniH=Tree[rank[lidM]].height;//the height where weight change begins
    int ProH=Tree[rank[lidM]].height; int ProID=lid;
    vector<set<int>> SCre;//record the shortcut change in each height
    set<int> ss;//ss.clear();
    SCre.assign(ProH+1,ss);

    //map<int,set<int>> DisRe;//rankid; record the distance change caused by the shortcut in each height
    //DisRe.clear();

    int MinH;

    bool tri=false;
    for(int i=0;i<Tree[rank[lidM]].vert.size();i++){
        if(Tree[rank[lidM]].vert[i].first==hid){
            if(Tree[rank[lidM]].vert[i].second.first>newW){
                Tree[rank[lidM]].vert[i].second.first=newW;
                Tree[rank[lidM]].vert[i].second.second=1;
                tri=true;
                SCre[ProH].insert(hid);
                MinH=IniH;
            }else if(Tree[rank[lidM]].vert[i].second.first==newW){
                Tree[rank[lidM]].vert[i].second.second+=1;
            }
            break;
        }
    }

    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed

    vector<int> ProIDRecord; ProIDRecord.assign(ProH+1,0);

    //int ProBeginH;
    int ProBeginID;
    if(tri){
//        cout<<"Bottom-up ;;;;;;;;;;;;;;;;;; "<<endl;
        while(ProH>=MinH){

            ProIDRecord[ProH]=ProID;
            int ProIDM=IDMap[ProID];
            vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProIDM]].vert;
            bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
            for(auto it=SCre[ProH].begin();it!=SCre[ProH].end();it++){
                int Cid=*it; int Cw;//=OCdis[make_pair(ProID,Cid)];
                int CidM=IDMap[Cid];
                int cidH=Tree[rank[CidM]].height-1;

                map<int,int> Hnei; //Hnei.clear();
                vector<pair<int,int>> Lnei; //Lnei.clear();
                for(int j=0;j<Vert.size();j++){
                    if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                        Hnei[Vert[j].first]=Vert[j].second.first;
                    }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                        Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                    }else{
                        Cw=Vert[j].second.first;
                    }
                }

                if(Tree[rank[ProIDM]].dis[cidH]>=Cw){
                    Tree[rank[ProIDM]].dis[cidH]=Cw;
                    Tree[rank[ProIDM]].FN[cidH]=true;
                    ProIDdisCha=true;
                    Tree[rank[ProIDM]].DisRe.insert(Cid);
                    //DisRe[rank[ProID]].insert(Cid); //cout<<"dischange Cid "<<Cid<<endl;
                }

                int hid,hidHeight,lid,lidHeight,wsum;
                for(int j=0;j<Tree[rank[CidM]].vert.size();j++){
                    hid=Tree[rank[CidM]].vert[j].first;hidHeight=Tree[rank[IDMap[hid]]].height-1;
                    if(Hnei.find(hid)!=Hnei.end()){
                        wsum=Cw+Hnei[hid];
                        if(wsum<Tree[rank[CidM]].vert[j].second.first){
                            Tree[rank[CidM]].vert[j].second.first=wsum;
                            Tree[rank[CidM]].vert[j].second.second=1;
                            SCre[Tree[rank[CidM]].height].insert(hid);
                            if(Tree[rank[CidM]].height<MinH) MinH=Tree[rank[CidM]].height;

                        }else if(wsum==Tree[rank[CidM]].vert[j].second.first){
                            Tree[rank[CidM]].vert[j].second.second+=1;
                        }

                    }
                }
                for(int j=0;j<Lnei.size();j++){
                    lid=Lnei[j].first;
                    int lidM=IDMap[lid];
                    lidHeight=Tree[rank[lidM]].height-1;
                    for(int k=0;k<Tree[rank[lidM]].vert.size();k++){
                        if(Tree[rank[lidM]].vert[k].first==Cid){
                            wsum=Cw+Lnei[j].second;
                            if(Tree[rank[lidM]].vert[k].second.first>wsum){
                                Tree[rank[lidM]].vert[k].second.first=wsum;
                                Tree[rank[lidM]].vert[k].second.second=1;
                                SCre[Tree[rank[lidM]].height].insert(Cid);
                                if(Tree[rank[lidM]].height<MinH) MinH=Tree[rank[lidM]].height;

                            }else if(Tree[rank[lidM]].vert[k].second.first==wsum){
                                Tree[rank[lidM]].vert[k].second.second+=1;
                            }

                            break;
                        }
                    }
                }
            }

            if(ProIDdisCha){//if the distance labeling is detected changed
                vertexIDChL.insert(ProID);
                //ProBeginH=ProH;
                ProBeginID=ProID;
            }

            ProH-=1;
            ProID=Tree[Tree[rank[ProIDM]].pa].uniqueVertex;
        }

        vector<int> linee; //linee.clear();
        linee.reserve(heightMax);
        int pachidd=Tree[Tree[rank[IDMap[ProBeginID]]].pa].uniqueVertex;
        while(Tree[rank[IDMap[pachidd]]].height>1){
            linee.insert(linee.begin(),pachidd);
            pachidd=Tree[Tree[rank[IDMap[pachidd]]].pa].uniqueVertex;
        }
        linee.insert(linee.begin(),pachidd);

        //top-down process
        EachNodeProBDis5Parti(rank[IDMap[ProBeginID]], linee, vertexIDChL, Tree, rank);
    }
    else{
//        cout<<"Not trigger update! "<<lid<<" "<<hid<<endl;
    }
    //return checkedDis.size();
}

void Graph::EachNodeProBDis5Parti(int child,vector<int>& line,set<int>& vertexIDChL, vector<Node> &Tree, vector<int> &rank){
    bool ProIDdisCha=false;
    int childID=Tree[child].uniqueVertex;
    int bHeight=BoundVertex[PartiTag[childID].first].size();
    vector<int> cntNew(line.size(),0);
    vector<bool> flagUpdate(line.size(),false);
    if(Tree[child].DisRe.size()!=0){
        for(int k=0;k<Tree[child].vert.size();k++){
            int b=Tree[child].vert[k].first, bH=Tree[rank[IDMap[b]]].height-1,vbW=Tree[child].vert[k].second.first;
            int bM=IDMap[b];
            if(Tree[child].FN[bH]){
                if(Tree[child].DisRe.find(b)!=Tree[child].DisRe.end()){//all ancestor check
                    for(int i=0;i<bH;i++){
//                        if(childID==1204 && Tree[child].vAncestor[i]==2218){
//                            cout<<"Find! "<<childID<<" "<<Tree[child].vAncestor[i]<<" "<<Tree[child].dis[i]<<endl;
//                        }
//                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[bM]].dis[i]){
                            Tree[child].dis[i]=vbW+Tree[rank[bM]].dis[i];
                            Tree[child].FN[i]=false;
                            Tree[child].cnt[i]=1;//new
                            ProIDdisCha=true;
                            flagUpdate[i]=true;
                            cntNew[i]=1;
                        }
                        else if(Tree[child].dis[i]==vbW+Tree[rank[bM]].dis[i]){
                            cntNew[i]++;
                            if(flagUpdate[i]) {
                                Tree[child].cnt[i]+=1;//new
                            }
                            else if(cntNew[i]>Tree[child].cnt[i]){
                                Tree[child].cnt[i]=cntNew[i];//use cntNew to redress the cnt value since the edge decrease may lead to more ways for dis (i.e., increase the cnt)
                            }
                        }
                    }
                    for(int i=bH+1;i<line.size();i++){
//                        if(childID==1204 && Tree[child].vAncestor[i]==2218){
//                            cout<<"Find! "<<childID<<" "<<Tree[child].vAncestor[i]<<" "<<Tree[child].dis[i]<<endl;
//                        }
//                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[IDMap[line[i]]]].dis[bH]){
                            Tree[child].dis[i]=vbW+Tree[rank[IDMap[line[i]]]].dis[bH];
                            Tree[child].FN[i]=false;
                            Tree[child].cnt[i]=1;//new
                            ProIDdisCha=true;
                            flagUpdate[i]=true;
                            cntNew[i]=1;
                        }
                        else if(Tree[child].dis[i]==vbW+Tree[rank[IDMap[line[i]]]].dis[bH]){
                            cntNew[i]++;
                            if(flagUpdate[i]) {
                                Tree[child].cnt[i]+=1;//new
                            }
                            else if(cntNew[i]>Tree[child].cnt[i]){
                                Tree[child].cnt[i]=cntNew[i];
                            }
                        }
                    }

                }else{//partial ancestor check

                    if(vertexIDChL.find(b)!=vertexIDChL.end()){
                        for(int i=0;i<bH;i++){
//                            if(childID==1204 && Tree[child].vAncestor[i]==2218){
//                                cout<<"Find2! "<<childID<<" "<<Tree[child].vAncestor[i]<<endl;
//                            }
//                            checkedDis.insert(make_pair(child,i));
                            if(Tree[child].dis[i]>vbW+Tree[rank[bM]].dis[i]){
                                Tree[child].dis[i]=vbW+Tree[rank[bM]].dis[i];
                                Tree[child].FN[i]=false;
                                Tree[child].cnt[i]=1;//new
                                ProIDdisCha=true;
                                flagUpdate[i]=true;
                                cntNew[i]=1;
                            }
                            else if(Tree[child].dis[i]==vbW+Tree[rank[bM]].dis[i]){
                                cntNew[i]++;
                                if(flagUpdate[i]) {
                                    Tree[child].cnt[i]+=1;//new
                                }
                                else if(cntNew[i]>Tree[child].cnt[i]){
                                    Tree[child].cnt[i]=cntNew[i];
                                }
                            }
                        }
                    }
                    for(int i=bH+1;i<line.size();i++){
//                        if(childID==1204 && Tree[child].vAncestor[i]==2218){
//                            cout<<"Find2! "<<childID<<" "<<Tree[child].vAncestor[i]<<endl;
//                        }
//                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[IDMap[line[i]]]].dis[bH]){
                            Tree[child].dis[i]=vbW+Tree[rank[IDMap[line[i]]]].dis[bH];
                            Tree[child].FN[i]=false;
                            Tree[child].cnt[i]=1;//new
                            ProIDdisCha=true;
                            flagUpdate[i]=true;
                            cntNew[i]=1;
                        }
                        else if(Tree[child].dis[i]==vbW+Tree[rank[IDMap[line[i]]]].dis[bH]){
                            cntNew[i]++;
                            if(flagUpdate[i]) {
                                Tree[child].cnt[i]+=1;//new
                            }
                            else if(cntNew[i]>Tree[child].cnt[i]){
                                Tree[child].cnt[i]=cntNew[i];
                            }
                        }
                    }

                }
            }
        }
    }else{
        for(int k=0;k<Tree[child].vert.size();k++){
            int b=Tree[child].vert[k].first, bH=Tree[rank[IDMap[b]]].height-1,vbW=Tree[child].vert[k].second.first;
            int bM=IDMap[b];
            if(Tree[child].FN[bH]){
                if(vertexIDChL.find(b)!=vertexIDChL.end()){
                    for(int i=0;i<bH;i++){
//                        if(childID==1204 && Tree[child].vAncestor[i]==2218){
//                            cout<<"Find3.1! "<<childID<<" "<<Tree[child].vAncestor[i]<<" "<<Tree[child].dis[i]<<" "<<vbW+Tree[rank[bM]].dis[i]<<endl;
//                        }
//                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[bM]].dis[i]){
                            Tree[child].dis[i]=vbW+Tree[rank[bM]].dis[i];
                            Tree[child].FN[i]=false;
                            Tree[child].cnt[i]=1;//new
                            ProIDdisCha=true;
                            flagUpdate[i]=true;
                            cntNew[i]=1;
                        }
                        else if(Tree[child].dis[i]==vbW+Tree[rank[bM]].dis[i]){
                            cntNew[i]++;
                            if(flagUpdate[i]) {
                                Tree[child].cnt[i]+=1;//new
                            }
                            else if(cntNew[i]>Tree[child].cnt[i]){
                                Tree[child].cnt[i]=cntNew[i];
                            }
                        }
                    }
                }
                for(int i=bH+1;i<line.size();i++){
//                    if(childID==1204 && Tree[child].vAncestor[i]==2218){
//                        cout<<"Find3.2! "<<childID<<" "<<Tree[child].vAncestor[i]<<" "<<Tree[child].dis[i]<<" "<<vbW+Tree[rank[IDMap[line[i]]]].dis[bH]<<endl;
//                    }
//                    checkedDis.insert(make_pair(child,i));
                    if(Tree[child].dis[i]>vbW+Tree[rank[IDMap[line[i]]]].dis[bH]){
                        Tree[child].dis[i]=vbW+Tree[rank[IDMap[line[i]]]].dis[bH];
                        Tree[child].FN[i]=false;
                        Tree[child].cnt[i]=1;//new
                        ProIDdisCha=true;
                        flagUpdate[i]=true;
                        cntNew[i]=1;
                    }
                    else if(Tree[child].dis[i]==vbW+Tree[rank[IDMap[line[i]]]].dis[bH]){
                        cntNew[i]++;
                        if(flagUpdate[i]) {
                            Tree[child].cnt[i]+=1;//new
                        }
                        else if(cntNew[i]>Tree[child].cnt[i]){
                            Tree[child].cnt[i]=cntNew[i];
                        }
                    }
                }
            }
        }
    }

    if(ProIDdisCha){
        vertexIDChL.insert(Tree[child].uniqueVertex);
        vUpdated[Tree[child].uniqueVertex]=true;
    }

    line.push_back(Tree[child].uniqueVertex);
    for(int i=0;i<Tree[child].ch.size();i++){
        EachNodeProBDis5Parti(Tree[child].ch[i], line, vertexIDChL, Tree, rank);
    }
    line.pop_back();

}

//batch update for overlay graph
void Graph::DecreaseOverlayBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU){
    map<int,int> checkedDis;

    for(int i=0;i<Tree.size();i++){
        Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
    }

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<OrderCompp>> SCre; //SCre.clear();
    SCre.assign(node_num,set<OrderCompp>());//{vertexID, set<int>}
    set<OrderCompMin> OC; //OC.clear();//vertexID in decreasing node order

    //    vector<int> ProBeginVertexSet; //ProBeginVertexSet.clear();
    vector<int> ProBeginVertexSetNew;
//    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed
    ProBeginVertexSetOverlay.clear(); vertexIDChLOverlay.clear();

    int a,b,newW,lid,hid;
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first; b=wBatch[k].first.second ;newW=wBatch[k].second.second;
        if(NodeOrder[a]<NodeOrder[b]){
            lid=a;hid=b;
        }else{
            lid=b;hid=a;
        }

        if(Neighbor[a].find(b)!=Neighbor[a].end()){
            Neighbor[a][b]=newW;
        }
//        else{
//            cout<<"Wrong for Neighbors! "<<a<<" "<<b<<endl;
//            exit(1);
//        }
        if(Neighbor[b].find(a)!=Neighbor[b].end()){
            Neighbor[b][a]=newW;
        }
//        else{
//            cout<<"Wrong for Neighbors! "<<b<<" "<<a<<endl;
//            exit(1);
//        }

        for(int i=0;i<Tree[rank[lid]].vert.size();i++){
            if(Tree[rank[lid]].vert[i].first==hid){
                if(Tree[rank[lid]].vert[i].second.first>newW){
                    Tree[rank[lid]].vert[i].second.first=newW;
                    Tree[rank[lid]].vert[i].second.second=1;
                    SCre[lid].insert(hid);
                    OC.insert(OrderCompMin(lid));
                }else if(Tree[rank[lid]].vert[i].second.first==newW){
                    Tree[rank[lid]].vert[i].second.second+=1;
                }
                break;
            }
        }

    }


    int ProID;
    //processing the stars
    while(!OC.empty()){
        ProID=(*OC.begin()).x;
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
        bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=it->x; int Cw;
            int cidH=Tree[rank[Cid]].height-1;

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }else{
                    Cw=Vert[j].second.first;
                }
            }

            if(algoChoice==H2H){
                if(Tree[rank[ProID]].dis[cidH]>Cw){
                    Tree[rank[ProID]].dis[cidH]=Cw;
                    Tree[rank[ProID]].FN[cidH]=true;
                    ProIDdisCha=true;
                    vUpdated[ProID]=true;
                    Tree[rank[ProID]].DisRe.insert(Cid);
                    Tree[rank[ProID]].cnt[cidH]=1;//new
                }else if(Tree[rank[ProID]].dis[cidH]==Cw){
                    Tree[rank[ProID]].FN[cidH]=true;
                    Tree[rank[ProID]].cnt[cidH]+=1;//new
                }

            }

            int hid,hidHeight,lid,lidHeight,wsum;
            for(int j=0;j<Tree[rank[Cid]].vert.size();j++){
                hid=Tree[rank[Cid]].vert[j].first;hidHeight=Tree[rank[hid]].height-1;
                if(Hnei.find(hid)!=Hnei.end()){
                    wsum=Cw+Hnei[hid];
                    if(wsum<Tree[rank[Cid]].vert[j].second.first){
                        Tree[rank[Cid]].vert[j].second.first=wsum;
                        Tree[rank[Cid]].vert[j].second.second=1;
                        SCre[Cid].insert(hid);
                        OC.insert(OrderCompMin(Cid));
                    }else if(wsum==Tree[rank[Cid]].vert[j].second.first){
                        Tree[rank[Cid]].vert[j].second.second+=1;
                    }

                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid=Lnei[j].first;lidHeight=Tree[rank[lid]].height-1;
                for(int k=0;k<Tree[rank[lid]].vert.size();k++){
                    if(Tree[rank[lid]].vert[k].first==Cid){
                        wsum=Cw+Lnei[j].second;
                        if(Tree[rank[lid]].vert[k].second.first>wsum){
                            Tree[rank[lid]].vert[k].second.first=wsum;
                            Tree[rank[lid]].vert[k].second.second=1;
                            SCre[lid].insert(Cid);
                            OC.insert(OrderCompMin(lid));
                        }else if(Tree[rank[lid]].vert[k].second.first==wsum){
                            if(SCre[ProID].find(lid)==SCre[ProID].end()) {//if not found, to avoid repeat count
                                Tree[rank[lid]].vert[k].second.second += 1;
                            }
                        }

                        break;
                    }
                }
            }
        }

        if(ProIDdisCha){//if the distance labeling is dectected changed
            vertexIDChLOverlay.insert(ProID);
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetOverlay.size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[ProID],r;
            for(int i=0;i<ProBeginVertexSetOverlay.size();i++){
                r=rank[ProBeginVertexSetOverlay[i]];
                if(LCAQueryOverlay(rnew,r)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetOverlay[i]);
                }
            }
            ProBeginVertexSetOverlay=ProBeginVertexSetNew;
        }
    }

    //cout<<"Finish bottom-up refresh"<<endl;

    if(ifLabelU){
        int ProBeginVertexID;
        for(int i=0;i<ProBeginVertexSetOverlay.size();i++){
            ProBeginVertexID=ProBeginVertexSetOverlay[i];
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
            while(Tree[rank[pachidd]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);
            EachNodeProBDis5(rank[ProBeginVertexID], linee, vertexIDChLOverlay,Tree,rank);
        }
    }

    //return checkedDis.size();
}

void Graph::EdgeInsertOverlayBatchCH(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax){
    map<int,int> checkedDis;

    for(int i=0;i<Tree.size();i++){
        Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
    }

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<OrderCompp>> SCre; //SCre.clear();
    SCre.assign(node_num,set<OrderCompp>());//{vertexID, set<int>}
    set<OrderCompMin> OC; //OC.clear();//vertexID in decreasing node order

    //    vector<int> ProBeginVertexSet; //ProBeginVertexSet.clear();
    vector<int> ProBeginVertexSetNew;
//    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed
    ProBeginVertexSetOverlay.clear(); vertexIDChLOverlay.clear();

    int a,b,newW,lid,hid;
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first; b=wBatch[k].first.second ;newW=wBatch[k].second.second;
        if(NodeOrder[a]<NodeOrder[b]){
            lid=a;hid=b;
        }else{
            lid=b;hid=a;
        }

        if(Neighbor[a].find(b)!=Neighbor[a].end()){
            cout<<"Already exist this edge. "<<a<<" "<<b<<" "<<Neighbor[a][b]<<endl;
            exit(1);
        }
        if(Neighbor[b].find(a)!=Neighbor[b].end()){
            cout<<"Already exist this edge. "<<a<<" "<<b<<" "<<Neighbor[b][a]<<endl;
            exit(1);
        }

        Neighbor[a].insert({b,newW});
        Neighbor[b].insert({a,newW});
        bool ifFind=false;
        for(int i=0;i<Tree[rank[lid]].vert.size();i++){
            if(Tree[rank[lid]].vert[i].first==hid){
                ifFind=true;
                if(Tree[rank[lid]].vert[i].second.first>newW){
                    Tree[rank[lid]].vert[i].second.first=newW;
                    Tree[rank[lid]].vert[i].second.second=1;
                    SCre[lid].insert(hid);
                    OC.insert(OrderCompMin(lid));
                }else if(Tree[rank[lid]].vert[i].second.first==newW){
                    Tree[rank[lid]].vert[i].second.second+=1;
                }
                break;
            }
        }

        if(!ifFind){
            int newSCNum=0;
            cout<<"Not found edge "<<lid<<" "<<hid<<" "<<newW<<endl;
            Tree[rank[lid]].vert.emplace_back(hid, make_pair(newW,1));
            newSCNum++;
            SCre[lid].insert(hid);
            OC.insert(OrderCompMin(lid));

            BottomUpNewShortcutInsertCH(lid,hid,newW,Tree,rank,newSCNum);
        }
    }


    int ProID;
    //processing the stars
    while(!OC.empty()){
        ProID=(*OC.begin()).x;
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
        bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=it->x; int Cw;
            int cidH=Tree[rank[Cid]].height-1;

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }else{
                    Cw=Vert[j].second.first;
                }
            }


            int hid,hidHeight,lid,lidHeight,wsum;
            for(int j=0;j<Tree[rank[Cid]].vert.size();j++){
                hid=Tree[rank[Cid]].vert[j].first;hidHeight=Tree[rank[hid]].height-1;
                if(Hnei.find(hid)!=Hnei.end()){
                    wsum=Cw+Hnei[hid];
                    if(wsum<Tree[rank[Cid]].vert[j].second.first){
                        Tree[rank[Cid]].vert[j].second.first=wsum;
                        Tree[rank[Cid]].vert[j].second.second=1;
                        SCre[Cid].insert(hid);
                        OC.insert(OrderCompMin(Cid));
                    }else if(wsum==Tree[rank[Cid]].vert[j].second.first){
                        Tree[rank[Cid]].vert[j].second.second+=1;
                    }

                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid=Lnei[j].first;lidHeight=Tree[rank[lid]].height-1;
                for(int k=0;k<Tree[rank[lid]].vert.size();k++){
                    if(Tree[rank[lid]].vert[k].first==Cid){
                        wsum=Cw+Lnei[j].second;
                        if(Tree[rank[lid]].vert[k].second.first>wsum){
                            Tree[rank[lid]].vert[k].second.first=wsum;
                            Tree[rank[lid]].vert[k].second.second=1;
                            SCre[lid].insert(Cid);
                            OC.insert(OrderCompMin(lid));
                        }else if(Tree[rank[lid]].vert[k].second.first==wsum){
                            if(SCre[ProID].find(lid)==SCre[ProID].end()) {//if not found, to avoid repeat count
                                Tree[rank[lid]].vert[k].second.second += 1;
                            }
                        }

                        break;
                    }
                }
            }
        }

        if(ProIDdisCha){//if the distance labeling is dectected changed
            vertexIDChLOverlay.insert(ProID);
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetOverlay.size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[ProID],r;
            for(int i=0;i<ProBeginVertexSetOverlay.size();i++){
                r=rank[ProBeginVertexSetOverlay[i]];
                if(LCAQueryOverlay(rnew,r)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetOverlay[i]);
                }
            }
            ProBeginVertexSetOverlay=ProBeginVertexSetNew;
        }
    }

    //cout<<"Finish bottom-up refresh"<<endl;



    //return checkedDis.size();
}

void Graph::EdgeInsertOverlayBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU){
    map<int,int> checkedDis;

    for(int i=0;i<Tree.size();i++){
        Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
    }

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<OrderCompp>> SCre; //SCre.clear();
    SCre.assign(node_num,set<OrderCompp>());//{vertexID, set<int>}
    set<OrderCompMin> OC; //OC.clear();//vertexID in decreasing node order

    //    vector<int> ProBeginVertexSet; //ProBeginVertexSet.clear();
    vector<int> ProBeginVertexSetNew;
//    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed
    ProBeginVertexSetOverlay.clear(); vertexIDChLOverlay.clear();

    int a,b,newW,lid,hid;
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first; b=wBatch[k].first.second ;newW=wBatch[k].second.second;
        if(NodeOrder[a]<NodeOrder[b]){
            lid=a;hid=b;
        }else{
            lid=b;hid=a;
        }

        if(Neighbor[a].find(b)!=Neighbor[a].end()){
            cout<<"Already exist this edge. "<<a<<" "<<b<<" "<<Neighbor[a][b]<<endl;
            exit(1);
        }
        if(Neighbor[b].find(a)!=Neighbor[b].end()){
            cout<<"Already exist this edge. "<<a<<" "<<b<<" "<<Neighbor[b][a]<<endl;
            exit(1);
        }

        Neighbor[a].insert({b,newW});
        Neighbor[b].insert({a,newW});
        bool ifFind=false;
        for(int i=0;i<Tree[rank[lid]].vert.size();i++){
            if(Tree[rank[lid]].vert[i].first==hid){
                ifFind=true;
                if(Tree[rank[lid]].vert[i].second.first>newW){
                    Tree[rank[lid]].vert[i].second.first=newW;
                    Tree[rank[lid]].vert[i].second.second=1;
                    SCre[lid].insert(hid);
                    OC.insert(OrderCompMin(lid));
                }else if(Tree[rank[lid]].vert[i].second.first==newW){
                    Tree[rank[lid]].vert[i].second.second+=1;
                }
                break;
            }
        }

        if(!ifFind){
            int newSCNum=0;
            cout<<"Not found edge "<<lid<<" "<<hid<<" "<<newW<<endl;
            ifFind=false;
            Tree[rank[lid]].vert.emplace_back(hid, make_pair(newW,1));
            newSCNum++;
            for(int i=0;i<Tree[rank[lid]].vAncestor.size();++i){
                if(Tree[rank[lid]].vAncestor[i]==hid){
                    Tree[rank[lid]].pos.push_back(i);
                    ifFind=true;
                }
            }
            SCre[lid].insert(hid);
            OC.insert(OrderCompMin(lid));
            if(!ifFind){//if false
                cout<<"Cross-branch edge! Should reorganize the tree structure! "<<lid<<" "<<hid<<" "<<newW<<endl;
                exit(1);
            }
            BottomUpNewShortcutInsert(lid,hid,newW,Tree,rank,newSCNum);
        }
    }


    int ProID;
    //processing the stars
    while(!OC.empty()){
        ProID=(*OC.begin()).x;
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
        bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=it->x; int Cw;
            int cidH=Tree[rank[Cid]].height-1;

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }else{
                    Cw=Vert[j].second.first;
                }
            }

            if(algoChoice==H2H){
                if(Tree[rank[ProID]].dis[cidH]>Cw){
                    Tree[rank[ProID]].dis[cidH]=Cw;
                    Tree[rank[ProID]].FN[cidH]=true;
                    ProIDdisCha=true;
                    vUpdated[ProID]=true;
                    Tree[rank[ProID]].DisRe.insert(Cid);
                    Tree[rank[ProID]].cnt[cidH]=1;//new
                }else if(Tree[rank[ProID]].dis[cidH]==Cw){
                    Tree[rank[ProID]].FN[cidH]=true;
                    Tree[rank[ProID]].cnt[cidH]+=1;//new
                }

            }

            int hid,hidHeight,lid,lidHeight,wsum;
            for(int j=0;j<Tree[rank[Cid]].vert.size();j++){
                hid=Tree[rank[Cid]].vert[j].first;hidHeight=Tree[rank[hid]].height-1;
                if(Hnei.find(hid)!=Hnei.end()){
                    wsum=Cw+Hnei[hid];
                    if(wsum<Tree[rank[Cid]].vert[j].second.first){
                        Tree[rank[Cid]].vert[j].second.first=wsum;
                        Tree[rank[Cid]].vert[j].second.second=1;
                        SCre[Cid].insert(hid);
                        OC.insert(OrderCompMin(Cid));
                    }else if(wsum==Tree[rank[Cid]].vert[j].second.first){
                        Tree[rank[Cid]].vert[j].second.second+=1;
                    }

                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid=Lnei[j].first;lidHeight=Tree[rank[lid]].height-1;
                for(int k=0;k<Tree[rank[lid]].vert.size();k++){
                    if(Tree[rank[lid]].vert[k].first==Cid){
                        wsum=Cw+Lnei[j].second;
                        if(Tree[rank[lid]].vert[k].second.first>wsum){
                            Tree[rank[lid]].vert[k].second.first=wsum;
                            Tree[rank[lid]].vert[k].second.second=1;
                            SCre[lid].insert(Cid);
                            OC.insert(OrderCompMin(lid));
                        }else if(Tree[rank[lid]].vert[k].second.first==wsum){
                            if(SCre[ProID].find(lid)==SCre[ProID].end()) {//if not found, to avoid repeat count
                                Tree[rank[lid]].vert[k].second.second += 1;
                            }
                        }

                        break;
                    }
                }
            }
        }

        if(ProIDdisCha){//if the distance labeling is dectected changed
            vertexIDChLOverlay.insert(ProID);
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetOverlay.size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[ProID],r;
            for(int i=0;i<ProBeginVertexSetOverlay.size();i++){
                r=rank[ProBeginVertexSetOverlay[i]];
                if(LCAQueryOverlay(rnew,r)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetOverlay[i]);
                }
            }
            ProBeginVertexSetOverlay=ProBeginVertexSetNew;
        }
    }

    //cout<<"Finish bottom-up refresh"<<endl;

    if(ifLabelU){
        int ProBeginVertexID;
        for(int i=0;i<ProBeginVertexSetOverlay.size();i++){
            ProBeginVertexID=ProBeginVertexSetOverlay[i];
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
            while(Tree[rank[pachidd]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);
            EachNodeProBDis5(rank[ProBeginVertexID], linee, vertexIDChLOverlay,Tree,rank);
        }
    }

    //return checkedDis.size();
}


void Graph::DecreaseOverlayBatchLabel(vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int> &ProBeginVertexSet, set<int> &vertexIDChL) {
    int ProBeginVertexID;
    for(int i=0;i<ProBeginVertexSet.size();i++){
        ProBeginVertexID=ProBeginVertexSet[i];
        vector<int> linee; //linee.clear();
        linee.reserve(heightMax);
        int pachidd=Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
        while(Tree[rank[pachidd]].height>1){
            linee.insert(linee.begin(),pachidd);
            pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
        }
        linee.insert(linee.begin(),pachidd);
        EachNodeProBDis5(rank[ProBeginVertexID], linee, vertexIDChL,Tree,rank);
    }
}

void Graph::BottomUpNewShortcutInsertParti(int lid, int hid, int weight, vector<Node> &Tree, vector<int> &rank, int &newSCNum){
    int lidM = IDMap[lid], hidM = IDMap[hid];
    int id1,id2,w1,newW;
    int fatherWeight=INF;
    bool ifFind=false;
    int fatherID=Tree[Tree[rank[lidM]].pa].uniqueVertex;
    for(int i=0;i<Tree[rank[lidM]].vert.size();++i){
        id1=Tree[rank[lidM]].vert[i].first; w1=Tree[rank[lidM]].vert[i].second.first;

        if(NodeOrder[id1]<NodeOrder[hid]){
            int id1M=IDMap[id1];
            ifFind=false;
            newW = weight+w1;
            newW = INF;
            for(int j=0;j<Tree[rank[id1M]].vert.size();++j){
                id2=Tree[rank[id1M]].vert[j].first;
                if(id2==hid){
//                    if(Tree[rank[id1M]].vert[j].second.first>newW){
//                        Tree[rank[id1M]].vert[j].second.first=newW;
//                        Tree[rank[id1M]].vert[j].second.second=1;
//                        cout<<"Update SC 1. "<<id1<<"("<<NodeOrder[id1]<<") "<<hid<<" "<<newW<<endl;
//                        if(id1==fatherID){
//                            fatherWeight=newW;
//                        }
//                    }else if(Tree[rank[id1M]].vert[j].second.first==newW){
//                        Tree[rank[id1M]].vert[j].second.second+=1;
//                        if(id1==fatherID){
//                            fatherWeight=newW;
//                        }
//                    }
                    ifFind=true;
                    break;
                }
            }
            if(!ifFind){//if not found
                Tree[rank[id1M]].vert.emplace_back(hid, make_pair(newW,1));
                ifFind=false;
                newSCNum++;
//                cout<<"New SC 1. "<<id1<<"("<<NodeOrder[id1]<<") "<<hid<<" "<<newW<<endl;
                for(int k=0;k<Tree[rank[id1M]].vAncestor.size();++k) {
                    if (Tree[rank[id1M]].vAncestor[k] == hid) {
//                        Tree[rank[id1M]].pos.push_back(k);
                        int temp=Tree[rank[id1M]].pos.back();
                        Tree[rank[id1M]].pos[Tree[rank[id1M]].pos.size()-1]=k;
                        Tree[rank[id1M]].pos.push_back(temp);
                        ifFind = true;
                        break;
                    }
                }
                if(!ifFind){
                    cout<<"Wrong. Not find hid. "<<id1<<" "<<hid<<endl; exit(1);
                }
            }
        }
        else if(NodeOrder[id1]>NodeOrder[hid]){
            ifFind=false;
            newW = weight+w1;
            newW = INF;
            for(int j=0;j<Tree[rank[hidM]].vert.size();++j){
                id2=Tree[rank[hidM]].vert[j].first;
                if(id2==id1){
//                    if(Tree[rank[hidM]].vert[j].second.first>newW){
//                        cout<<"Update SC 2. "<<hid<<"("<<NodeOrder[hid]<<") "<<id1<<"("<<NodeOrder[id1]<<") "<<newW<<endl;
//                        Tree[rank[hidM]].vert[j].second.first=newW;
//                        Tree[rank[hidM]].vert[j].second.second=1;
//                    }else if(Tree[rank[hidM]].vert[j].second.first==newW){
//                        Tree[rank[hidM]].vert[j].second.second+=1;
//                    }
                    ifFind=true;
                    break;
                }
            }
            if(!ifFind){//if not found
                Tree[rank[hidM]].vert.emplace_back(id1, make_pair(newW,1));
                ifFind=false;
                newSCNum++;
//                cout<<"New SC 2. "<<hid<<"("<<NodeOrder[hid]<<") "<<id1<<"("<<NodeOrder[id1]<<") "<<newW<<endl;
                for(int k=0;k<Tree[rank[hidM]].vAncestor.size();++k) {
                    if (Tree[rank[hidM]].vAncestor[k] == id1) {
//                        Tree[rank[hidM]].pos.push_back(k);
                        int temp=Tree[rank[hidM]].pos.back();
                        Tree[rank[hidM]].pos[Tree[rank[hidM]].pos.size()-1]=k;
                        Tree[rank[hidM]].pos.push_back(temp);
                        ifFind = true;
                        break;
                    }
                }
                if(!ifFind){
                    cout<<"Wrong. Not find id1. "<<id1<<" "<<hid<<endl; exit(1);
                }
            }
        }


    }

    if(Tree[rank[lidM]].height>Tree[rank[hidM]].height){
//        cout<<lid<<"("<<Tree[rank[lidM]].height<<") "<<hid<<"("<<Tree[rank[hidM]].height<<") "<<fatherWeight<<endl;
        BottomUpNewShortcutInsertParti(fatherID,hid,fatherWeight,Tree,rank,newSCNum);
    }
}

void Graph::BottomUpNewShortcutInsertPartiCH(int lid, int hid, int weight, vector<Node> &Tree, vector<int> &rank, int &newSCNum){
    int lidM = IDMap[lid], hidM = IDMap[hid];
    int id1,id2,w1,newW;
    int fatherWeight=INF;
    bool ifFind=false;
    int fatherID=Tree[Tree[rank[lidM]].pa].uniqueVertex;
    for(int i=0;i<Tree[rank[lidM]].vert.size();++i){
        id1=Tree[rank[lidM]].vert[i].first; w1=Tree[rank[lidM]].vert[i].second.first;

        if(NodeOrder[id1]<NodeOrder[hid]){
            int id1M=IDMap[id1];
            ifFind=false;
            newW = weight+w1;
            newW = INF;
            for(int j=0;j<Tree[rank[id1M]].vert.size();++j){
                id2=Tree[rank[id1M]].vert[j].first;
                if(id2==hid){
                    ifFind=true;
                    break;
                }
            }
            if(!ifFind){//if not found
                Tree[rank[id1M]].vert.emplace_back(hid, make_pair(newW,1));
                ifFind=false;
                newSCNum++;
            }
        }
        else if(NodeOrder[id1]>NodeOrder[hid]){
            ifFind=false;
            newW = weight+w1;
            newW = INF;
            for(int j=0;j<Tree[rank[hidM]].vert.size();++j){
                id2=Tree[rank[hidM]].vert[j].first;
                if(id2==id1){
                    ifFind=true;
                    break;
                }
            }
            if(!ifFind){//if not found
                Tree[rank[hidM]].vert.emplace_back(id1, make_pair(newW,1));
                ifFind=false;
                newSCNum++;

            }
        }


    }

    if(Tree[rank[lidM]].height>Tree[rank[hidM]].height){
//        cout<<lid<<"("<<Tree[rank[lidM]].height<<") "<<hid<<"("<<Tree[rank[hidM]].height<<") "<<fatherWeight<<endl;
        BottomUpNewShortcutInsertPartiCH(fatherID,hid,fatherWeight,Tree,rank,newSCNum);
    }
}

void Graph::BottomUpNewShortcutInsert(int lid, int hid, int weight, vector<Node> &Tree, vector<int> &rank, int &newSCNum){
    int id1,id2,w1,newW;
    int fatherWeight=INF;
    bool ifFind=false;
    int fatherID=Tree[Tree[rank[lid]].pa].uniqueVertex;
    for(int i=0;i<Tree[rank[lid]].vert.size();++i){
        id1=Tree[rank[lid]].vert[i].first; w1=Tree[rank[lid]].vert[i].second.first;

        if(NodeOrder[id1]<NodeOrder[hid]){
            ifFind=false;
            for(int j=0;j<Tree[rank[id1]].vert.size();++j){
                id2=Tree[rank[id1]].vert[j].first;
                if(id2==hid){
                    ifFind=true;
                    break;
                }
            }
            if(!ifFind){//if not found
                newW = weight+w1;
                newW = INF;
                Tree[rank[id1]].vert.emplace_back(hid, make_pair(newW,1));
                newSCNum++;
                ifFind=false;
                for(int k=0;k<Tree[rank[id1]].vAncestor.size();++k) {
                    if (Tree[rank[id1]].vAncestor[k] == hid) {
//                        Tree[rank[id1]].pos.push_back(k);
                        int temp=Tree[rank[id1]].pos.back();
                        Tree[rank[id1]].pos[Tree[rank[id1]].pos.size()-1]=k;
                        Tree[rank[id1]].pos.push_back(temp);
                        ifFind = true;
                        break;
                    }
                }
                if(!ifFind){
                    cout<<"Wrong. Not find hid. "<<id1<<" "<<hid<<endl; exit(1);
                }
//                if(id1==fatherID){
//                    fatherWeight=newW;
//                }
            }
        }
        else if(NodeOrder[id1]>NodeOrder[hid]){
            ifFind=false;
            for(int j=0;j<Tree[rank[hid]].vert.size();++j){
                id2=Tree[rank[hid]].vert[j].first;
                if(id2==id1){
                    ifFind=true; break;
                }
            }
            if(!ifFind){//if not found
                newW = weight+w1;
                newW = INF;
                Tree[rank[hid]].vert.emplace_back(id1, make_pair(newW,1));
                newSCNum++;
                ifFind=false;
                for(int k=0;k<Tree[rank[hid]].vAncestor.size();++k) {
                    if (Tree[rank[hid]].vAncestor[k] == id1) {
//                        Tree[rank[hid]].pos.push_back(k);
                        int temp=Tree[rank[hid]].pos.back();
                        Tree[rank[hid]].pos[Tree[rank[hid]].pos.size()-1]=k;
                        Tree[rank[hid]].pos.push_back(temp);
                        ifFind = true;
                        break;
                    }
                }
                if(!ifFind){
                    cout<<"Wrong. Not find id1. "<<id1<<" "<<hid<<endl; exit(1);
                }
            }
        }


    }

    if(Tree[rank[lid]].height>Tree[rank[hid]].height){
//        cout<<lid<<"("<<Tree[rank[lid]].height<<") "<<hid<<"("<<Tree[rank[hid]].height<<") "<<fatherWeight<<endl;
        BottomUpNewShortcutInsert(fatherID,hid,fatherWeight,Tree,rank,newSCNum);
    }
}

void Graph::BottomUpNewShortcutInsertCH(int lid, int hid, int weight, vector<Node> &Tree, vector<int> &rank, int &newSCNum){
    int id1,id2,w1,newW;
    int fatherWeight=INF;
    bool ifFind=false;
    int fatherID=Tree[Tree[rank[lid]].pa].uniqueVertex;
    for(int i=0;i<Tree[rank[lid]].vert.size();++i){
        id1=Tree[rank[lid]].vert[i].first; w1=Tree[rank[lid]].vert[i].second.first;

        if(NodeOrder[id1]<NodeOrder[hid]){
            ifFind=false;
            for(int j=0;j<Tree[rank[id1]].vert.size();++j){
                id2=Tree[rank[id1]].vert[j].first;
                if(id2==hid){
                    ifFind=true;
                    break;
                }
            }
            if(!ifFind){//if not found
                newW = weight+w1;
                newW = INF;
                Tree[rank[id1]].vert.emplace_back(hid, make_pair(newW,1));
                newSCNum++;
            }
        }
        else if(NodeOrder[id1]>NodeOrder[hid]){
            ifFind=false;
            for(int j=0;j<Tree[rank[hid]].vert.size();++j){
                id2=Tree[rank[hid]].vert[j].first;
                if(id2==id1){
                    ifFind=true;
                    break;
                }
            }
            if(!ifFind){//if not found
                newW = weight+w1;
                newW = INF;
                Tree[rank[hid]].vert.emplace_back(id1, make_pair(newW,1));
                newSCNum++;
            }
        }

    }

    if(Tree[rank[lid]].height>Tree[rank[hid]].height){
//        cout<<lid<<"("<<Tree[rank[lid]].height<<") "<<hid<<"("<<Tree[rank[hid]].height<<") "<<fatherWeight<<endl;
        BottomUpNewShortcutInsertCH(fatherID,hid,fatherWeight,Tree,rank,newSCNum);
    }
}

void Graph::EdgeInsertPartiBatchCH(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>>& Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<pair<pair<int,int>,int>>& updatedSC){
    map<int,int> checkedDis;

    for(int i=0;i<Tree.size();i++){
        Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
    }

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<OrderCompp>> SCre; //SCre.clear();
    SCre.assign(node_num,set<OrderCompp>());//{vertexID, set<int>}
    set<OrderCompMin> OC; //OC.clear();//vertexID in decreasing node order

    int a,b,newW,lid,hid,lidM,hidM;
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first; b=wBatch[k].first.second ;newW=wBatch[k].second.second;
        if(NodeOrder[a]<NodeOrder[b]){
            lid=a;hid=b;
        }else{
            lid=b;hid=a;
        }
        lidM = IDMap[lid]; hidM = IDMap[hid];

//        if(pid==23){
//            cout<<pid<<": "<<a<<" "<<b<<" "<<newW<<endl;
//        }
        bool ifFind=false;
        for(int i=0;i<Neighbors[a].size();i++){
            if(Neighbors[a][i].first==b){
//                cout<<"Already exist this edge. "<<a<<" "<<b<<" "<<Neighbors[a][i].second<<endl;
                if(Neighbors[a][i].second>newW){
                    Neighbors[a][i].second=newW;
                }else{
//                    cout<<"Wrong update. "<<a<<" "<<b<<" "<<Neighbors[a][i].second<<" "<<newW<<endl; exit(1);
                }

                ifFind=true; break;
//                exit(1);
            }
        }

        for(int i=0;i<Neighbors[b].size();i++){
            if(Neighbors[b][i].first==a){
//                cout<<"Already exist this edge. "<<a<<" "<<b<<" "<<Neighbors[b][i].second<<endl;
                if(Neighbors[b][i].second>newW){
                    Neighbors[b][i].second=newW;
                }else{
//                    cout<<"Wrong update. "<<a<<" "<<b<<" "<<Neighbors[b][i].second<<" "<<newW<<endl; exit(1);
                }

                ifFind=true; break;
//                exit(1);
            }
        }

        if(!ifFind) {
            Neighbors[a].emplace_back(b, newW);
            Neighbors[b].emplace_back(a, newW);
        }
        ifFind=false;
        for(int i=0;i<Tree[rank[lidM]].vert.size();i++){
//            cout<<"vert "<<lid<<" "<<Tree[rank[lidM]].vert[i].first<<endl;
            if(Tree[rank[lidM]].vert[i].first==hid){
                ifFind=true;
                if(Tree[rank[lidM]].vert[i].second.first>newW){
                    Tree[rank[lidM]].vert[i].second.first=newW;
                    Tree[rank[lidM]].vert[i].second.second=1;
                    SCre[lid].insert(hid);
                    OC.insert(OrderCompMin(lid));
                }else if(Tree[rank[lidM]].vert[i].second.first==newW){
                    Tree[rank[lidM]].vert[i].second.second+=1;
                }
                break;
            }
        }
        int newSCNum=0;
        if(!ifFind){
            //            cout<<"Not found edge "<<lid<<" "<<hid<<" "<<newW<<endl;
            ifFind=false;
            Tree[rank[lidM]].vert.emplace_back(hid, make_pair(newW,1));
            newSCNum++;
//            cout<<"New SC 1. "<<lid<<" "<<hid<<" "<<newW<<endl;

            SCre[lid].insert(hid);
            OC.insert(OrderCompMin(lid));
            BottomUpNewShortcutInsertPartiCH(lid,hid,newW,Tree,rank,newSCNum);
//            cout<<"newSCNum: "<<newSCNum<<endl;
        }

    }

//    vector<int> ProBeginVertexSet; //ProBeginVertexSet.clear();
//    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed
    ProBeginVertexSetParti[pid].clear(); vertexIDChLParti[pid].clear();
    vector<int> ProBeginVertexSetNew;
    int ProBeginVertexID;
    int ProID, ProIDM;
    //processing the stars
    while(!OC.empty()){
        ProID=(*OC.begin()).x;
        ProIDM=IDMap[ProID];
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProIDM]].vert;
        bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=it->x; int Cw; int CidM=IDMap[Cid];
            int cidH=Tree[rank[CidM]].height-1;

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }else{
                    Cw=Vert[j].second.first;
                }
            }

            if(PartiTag[ProID].second){//if boundary vertex
                updatedSC.emplace_back(make_pair(ProID,Cid),Cw);
            }



            int hid2,hidHeight2,lid2,lidHeight2,wsum,lid2M;
            for(int j=0;j<Tree[rank[CidM]].vert.size();j++){
                hid2=Tree[rank[CidM]].vert[j].first;hidHeight2=Tree[rank[IDMap[hid2]]].height-1;
                if(Hnei.find(hid2)!=Hnei.end()){
                    wsum=Cw+Hnei[hid2];
                    if(wsum<Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.first=wsum;
                        Tree[rank[CidM]].vert[j].second.second=1;
                        SCre[Cid].insert(hid2);
                        OC.insert(OrderCompMin(Cid));
                    }else if(wsum==Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.second+=1;
                    }

                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid2=Lnei[j].first; lid2M=IDMap[lid2];
                lidHeight2=Tree[rank[lid2M]].height-1;
                for(int k=0;k<Tree[rank[lid2M]].vert.size();k++){
                    if(Tree[rank[lid2M]].vert[k].first==Cid){
                        wsum=Cw+Lnei[j].second;
                        if(Tree[rank[lid2M]].vert[k].second.first>wsum){
                            Tree[rank[lid2M]].vert[k].second.first=wsum;
                            Tree[rank[lid2M]].vert[k].second.second=1;
                            SCre[lid2].insert(Cid);
                            OC.insert(OrderCompMin(lid2));
                        }else if(Tree[rank[lid2M]].vert[k].second.first==wsum){
                            if(SCre[ProID].find(OrderCompp(lid2))==SCre[ProID].end()) {//if not found, avoid repeated count
                                Tree[rank[lid2M]].vert[k].second.second += 1;
                            }
                        }

                        break;
                    }
                }
            }
        }

//        if(ProIDdisCha){//if the distance labeling is dectected changed
//            vertexIDChLParti[pid].insert(ProID);
//            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetParti[pid].size()+1);
//            ProBeginVertexSetNew.push_back(ProID);
//            int rnew=rank[IDMap[ProID]],r;
//            for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
//                r=rank[IDMap[ProBeginVertexSetParti[pid][i]]];
//                if(LCAQueryPartition(rnew,r,pid)!=rnew){
//                    ProBeginVertexSetNew.push_back(ProBeginVertexSetParti[pid][i]);
//                }
//            }
//            ProBeginVertexSetParti[pid]=ProBeginVertexSetNew;
//        }
    }

    //cout<<"Finish bottom-up refresh"<<endl;

}


//batch update for partition graph of PH2H, edge insertion
void Graph::EdgeInsertPartiBatchH2H(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>>& Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<pair<pair<int,int>,int>>& updatedSC, bool ifLabelU){
    map<int,int> checkedDis;

    for(int i=0;i<Tree.size();i++){
        Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
    }

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<OrderCompp>> SCre; //SCre.clear();
    SCre.assign(node_num,set<OrderCompp>());//{vertexID, set<int>}
    set<OrderCompMin> OC; //OC.clear();//vertexID in decreasing node order

    int a,b,newW,lid,hid,lidM,hidM;
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first; b=wBatch[k].first.second ;newW=wBatch[k].second.second;
        if(NodeOrder[a]<NodeOrder[b]){
            lid=a;hid=b;
        }else{
            lid=b;hid=a;
        }
        lidM = IDMap[lid]; hidM = IDMap[hid];

//        if(pid==23){
//            cout<<pid<<": "<<a<<" "<<b<<" "<<newW<<endl;
//        }

        for(int i=0;i<Neighbors[a].size();i++){
            if(Neighbors[a][i].first==b){
                cout<<"Already exist this edge. "<<a<<" "<<b<<" "<<Neighbors[a][i].second<<endl;
                exit(1);
            }
        }
        for(int i=0;i<Neighbors[b].size();i++){
            if(Neighbors[b][i].first==a){
                cout<<"Already exist this edge. "<<a<<" "<<b<<" "<<Neighbors[b][i].second<<endl;
                exit(1);
            }
        }

        Neighbors[a].emplace_back(b,newW);
        Neighbors[b].emplace_back(a,newW);
        bool ifFind=false;
        for(int i=0;i<Tree[rank[lidM]].vert.size();i++){
//            cout<<"vert "<<lid<<" "<<Tree[rank[lidM]].vert[i].first<<endl;
            if(Tree[rank[lidM]].vert[i].first==hid){
                ifFind=true;
                if(Tree[rank[lidM]].vert[i].second.first>newW){
                    Tree[rank[lidM]].vert[i].second.first=newW;
                    Tree[rank[lidM]].vert[i].second.second=1;
                    SCre[lid].insert(hid);
                    OC.insert(OrderCompMin(lid));
                }else if(Tree[rank[lidM]].vert[i].second.first==newW){
                    Tree[rank[lidM]].vert[i].second.second+=1;
                }
                break;
            }
        }
        int newSCNum=0;
        if(!ifFind){
//            cout<<"Not found edge "<<lid<<" "<<hid<<" "<<newW<<endl;
            ifFind=false;
            Tree[rank[lidM]].vert.emplace_back(hid, make_pair(newW,1));
            newSCNum++;
//            cout<<"New SC 1. "<<lid<<" "<<hid<<" "<<newW<<endl;
            for(int i=0;i<Tree[rank[lidM]].vAncestor.size();++i){
                if(Tree[rank[lidM]].vAncestor[i]==hid){
                    int temp=Tree[rank[lidM]].pos.back();
                    Tree[rank[lidM]].pos[Tree[rank[lidM]].pos.size()-1]=i;
                    Tree[rank[lidM]].pos.push_back(temp);
                    ifFind=true;
                }
            }
            SCre[lid].insert(hid);
            OC.insert(OrderCompMin(lid));
            if(!ifFind){//if false
                cout<<"Cross-branch edge! Should reorganize the tree structure! "<<lid<<" "<<hid<<" "<<newW<<endl;
                exit(1);
            }
            BottomUpNewShortcutInsertParti(lid,hid,newW,Tree,rank,newSCNum);
//            cout<<"newSCNum: "<<newSCNum<<endl;
        }
    }

//    vector<int> ProBeginVertexSet; //ProBeginVertexSet.clear();
//    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed
    ProBeginVertexSetParti[pid].clear(); vertexIDChLParti[pid].clear();
    vector<int> ProBeginVertexSetNew;
    int ProBeginVertexID;
    int ProID, ProIDM;
    //processing the stars
    while(!OC.empty()){
        ProID=(*OC.begin()).x;
        ProIDM=IDMap[ProID];
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProIDM]].vert;
        bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=it->x; int Cw; int CidM=IDMap[Cid];
            int cidH=Tree[rank[CidM]].height-1;

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }else{
                    Cw=Vert[j].second.first;
                }
            }

            if(PartiTag[ProID].second){//if boundary vertex
                updatedSC.emplace_back(make_pair(ProID,Cid),Cw);
            }

            if(algoChoice==H2H){
                if(Tree[rank[ProIDM]].dis[cidH]>Cw){
                    Tree[rank[ProIDM]].dis[cidH]=Cw;
                    Tree[rank[ProIDM]].FN[cidH]=true;
                    ProIDdisCha=true;
                    vUpdated[ProID]=true;
                    Tree[rank[ProIDM]].DisRe.insert(Cid);
                    Tree[rank[ProIDM]].cnt[cidH]=1;//new
                }else if(Tree[rank[ProIDM]].dis[cidH]==Cw){
                    Tree[rank[ProIDM]].FN[cidH]=true;
                    Tree[rank[ProIDM]].cnt[cidH]+=1;//new
                }
            }


            int hid2,hidHeight2,lid2,lidHeight2,wsum,lid2M;
            for(int j=0;j<Tree[rank[CidM]].vert.size();j++){
                hid2=Tree[rank[CidM]].vert[j].first;hidHeight2=Tree[rank[IDMap[hid2]]].height-1;
                if(Hnei.find(hid2)!=Hnei.end()){
                    wsum=Cw+Hnei[hid2];
                    if(wsum<Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.first=wsum;
                        Tree[rank[CidM]].vert[j].second.second=1;
                        SCre[Cid].insert(hid2);
                        OC.insert(OrderCompMin(Cid));
                    }else if(wsum==Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.second+=1;
                    }

                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid2=Lnei[j].first; lid2M=IDMap[lid2];
                lidHeight2=Tree[rank[lid2M]].height-1;
                for(int k=0;k<Tree[rank[lid2M]].vert.size();k++){
                    if(Tree[rank[lid2M]].vert[k].first==Cid){
                        wsum=Cw+Lnei[j].second;
                        if(Tree[rank[lid2M]].vert[k].second.first>wsum){
                            Tree[rank[lid2M]].vert[k].second.first=wsum;
                            Tree[rank[lid2M]].vert[k].second.second=1;
                            SCre[lid2].insert(Cid);
                            OC.insert(OrderCompMin(lid2));
                        }else if(Tree[rank[lid2M]].vert[k].second.first==wsum){
                            if(SCre[ProID].find(OrderCompp(lid2))==SCre[ProID].end()) {//if not found, avoid repeated count
                                Tree[rank[lid2M]].vert[k].second.second += 1;
                            }
                        }

                        break;
                    }
                }
            }
        }

        if(ProIDdisCha){//if the distance labeling is dectected changed
            vertexIDChLParti[pid].insert(ProID);
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetParti[pid].size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[IDMap[ProID]],r;
            for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
                r=rank[IDMap[ProBeginVertexSetParti[pid][i]]];
                if(LCAQueryPartition(rnew,r,pid)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetParti[pid][i]);
                }
            }
            ProBeginVertexSetParti[pid]=ProBeginVertexSetNew;
        }
    }

    //cout<<"Finish bottom-up refresh"<<endl;
    if(ifLabelU){
        for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
            ProBeginVertexID=ProBeginVertexSetParti[pid][i];
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[IDMap[ProBeginVertexID]]].pa].uniqueVertex;
            while(Tree[rank[IDMap[pachidd]]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[IDMap[pachidd]]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);
            EachNodeProBDis5Parti(rank[IDMap[ProBeginVertexID]], linee, vertexIDChLParti[pid],Tree,rank);
        }
    }

    //return checkedDis.size();
}

//batch update for partition graph of PH2H, edge insertion
void Graph::EdgeInsertPartiBatchH2HPre(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>>& Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU){
    map<int,int> checkedDis;

    for(int i=0;i<Tree.size();i++){
        Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
    }

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<OrderCompp>> SCre; //SCre.clear();
    SCre.assign(node_num,set<OrderCompp>());//{vertexID, set<int>}
    set<OrderCompMin> OC; //OC.clear();//vertexID in decreasing node order

    int a,b,newW,lid,hid,lidM,hidM;
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first; b=wBatch[k].first.second ;newW=wBatch[k].second.second;
        if(NodeOrder[a]<NodeOrder[b]){
            lid=a;hid=b;
        }else{
            lid=b;hid=a;
        }
        lidM = IDMap[lid]; hidM = IDMap[hid];
        bool ifFind=false;
        for(int i=0;i<Neighbors[a].size();i++){
            if(Neighbors[a][i].first==b){
//                cout<<"Already exist this edge. "<<a<<" "<<b<<" "<<Neighbors[a][i].second<<endl;
                if(Neighbors[a][i].second>newW){
                    Neighbors[a][i].second=newW;
                }else{
//                    cout<<"Wrong update. "<<a<<" "<<b<<" "<<Neighbors[a][i].second<<" "<<newW<<endl; exit(1);
                }

                ifFind=true; break;
//                exit(1);
            }
        }

        for(int i=0;i<Neighbors[b].size();i++){
            if(Neighbors[b][i].first==a){
//                cout<<"Already exist this edge. "<<a<<" "<<b<<" "<<Neighbors[b][i].second<<endl;
                if(Neighbors[b][i].second>newW){
                    Neighbors[b][i].second=newW;
                }else{
//                    cout<<"Wrong update. "<<a<<" "<<b<<" "<<Neighbors[b][i].second<<" "<<newW<<endl; exit(1);
                }

                ifFind=true; break;
//                exit(1);
            }
        }

        if(!ifFind) {
            Neighbors[a].emplace_back(b, newW);
            Neighbors[b].emplace_back(a, newW);
        }
        ifFind=false;
        for(int i=0;i<Tree[rank[lidM]].vert.size();i++){
//            cout<<"vert "<<lid<<" "<<Tree[rank[lidM]].vert[i].first<<endl;
            if(Tree[rank[lidM]].vert[i].first==hid){
                ifFind=true;
                if(Tree[rank[lidM]].vert[i].second.first>newW){
                    Tree[rank[lidM]].vert[i].second.first=newW;
                    Tree[rank[lidM]].vert[i].second.second=1;
                    SCre[lid].insert(hid);
                    OC.insert(OrderCompMin(lid));
                }else if(Tree[rank[lidM]].vert[i].second.first==newW){
                    Tree[rank[lidM]].vert[i].second.second+=1;
                }
                break;
            }
        }
        int newSCNum=0;
        if(!ifFind){
//            cout<<"Not found edge "<<lid<<" "<<hid<<" "<<newW<<endl;
            ifFind=false;
            Tree[rank[lidM]].vert.emplace_back(hid, make_pair(newW,1));
            newSCNum++;
//            cout<<"New SC 1. "<<lid<<" "<<hid<<" "<<newW<<endl;
            for(int i=0;i<Tree[rank[lidM]].vAncestor.size();++i){
                if(Tree[rank[lidM]].vAncestor[i]==hid){
                    int temp=Tree[rank[lidM]].pos.back();
                    Tree[rank[lidM]].pos[Tree[rank[lidM]].pos.size()-1]=i;
                    Tree[rank[lidM]].pos.push_back(temp);
                    ifFind=true;
                }
            }
            SCre[lid].insert(hid);
            OC.insert(OrderCompMin(lid));
            if(!ifFind){//if false
                cout<<"Cross-branch edge! Should reorganize the tree structure! "<<lid<<" "<<hid<<" "<<newW<<endl;
                exit(1);
            }
            BottomUpNewShortcutInsertParti(lid,hid,newW,Tree,rank,newSCNum);
//            cout<<"newSCNum: "<<newSCNum<<endl;
        }
    }

//    vector<int> ProBeginVertexSet; //ProBeginVertexSet.clear();
//    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed
    ProBeginVertexSetParti[pid].clear(); vertexIDChLParti[pid].clear();
    vector<int> ProBeginVertexSetNew;
    int ProBeginVertexID;
    int ProID, ProIDM;
    //processing the stars
    while(!OC.empty()){
        ProID=(*OC.begin()).x;
        ProIDM=IDMap[ProID];
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProIDM]].vert;
        bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=it->x; int Cw; int CidM=IDMap[Cid];
            int cidH=Tree[rank[CidM]].height-1;

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }else{
                    Cw=Vert[j].second.first;
                }
            }


            if(algoChoice==H2H){
                if(Tree[rank[ProIDM]].dis[cidH]>Cw){
                    Tree[rank[ProIDM]].dis[cidH]=Cw;
                    Tree[rank[ProIDM]].FN[cidH]=true;
                    ProIDdisCha=true;
                    vUpdated[ProID]=true;
                    Tree[rank[ProIDM]].DisRe.insert(Cid);
                    Tree[rank[ProIDM]].cnt[cidH]=1;//new
                }else if(Tree[rank[ProIDM]].dis[cidH]==Cw){
                    Tree[rank[ProIDM]].FN[cidH]=true;
                    Tree[rank[ProIDM]].cnt[cidH]+=1;//new
                }
            }


            int hid2,hidHeight2,lid2,lidHeight2,wsum,lid2M;
            for(int j=0;j<Tree[rank[CidM]].vert.size();j++){
                hid2=Tree[rank[CidM]].vert[j].first;hidHeight2=Tree[rank[IDMap[hid2]]].height-1;
                if(Hnei.find(hid2)!=Hnei.end()){
                    wsum=Cw+Hnei[hid2];
                    if(wsum<Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.first=wsum;
                        Tree[rank[CidM]].vert[j].second.second=1;
                        SCre[Cid].insert(hid2);
                        OC.insert(OrderCompMin(Cid));
                    }else if(wsum==Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.second+=1;
                    }

                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid2=Lnei[j].first; lid2M=IDMap[lid2];
                lidHeight2=Tree[rank[lid2M]].height-1;
                for(int k=0;k<Tree[rank[lid2M]].vert.size();k++){
                    if(Tree[rank[lid2M]].vert[k].first==Cid){
                        wsum=Cw+Lnei[j].second;
                        if(Tree[rank[lid2M]].vert[k].second.first>wsum){
                            Tree[rank[lid2M]].vert[k].second.first=wsum;
                            Tree[rank[lid2M]].vert[k].second.second=1;
                            SCre[lid2].insert(Cid);
                            OC.insert(OrderCompMin(lid2));
                        }else if(Tree[rank[lid2M]].vert[k].second.first==wsum){
                            if(SCre[ProID].find(OrderCompp(lid2))==SCre[ProID].end()) {//if not found, avoid repeated count
                                Tree[rank[lid2M]].vert[k].second.second += 1;
                            }
                        }

                        break;
                    }
                }
            }
        }

        if(ProIDdisCha){//if the distance labeling is dectected changed
            vertexIDChLParti[pid].insert(ProID);
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetParti[pid].size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[IDMap[ProID]],r;
            for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
                r=rank[IDMap[ProBeginVertexSetParti[pid][i]]];
                if(LCAQueryPartition(rnew,r,pid)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetParti[pid][i]);
                }
            }
            ProBeginVertexSetParti[pid]=ProBeginVertexSetNew;
        }
    }

    //cout<<"Finish bottom-up refresh"<<endl;
    if(ifLabelU){
        for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
            ProBeginVertexID=ProBeginVertexSetParti[pid][i];
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[IDMap[ProBeginVertexID]]].pa].uniqueVertex;
            while(Tree[rank[IDMap[pachidd]]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[IDMap[pachidd]]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);
            EachNodeProBDis5Parti(rank[IDMap[ProBeginVertexID]], linee, vertexIDChLParti[pid],Tree,rank);
        }
    }

    //return checkedDis.size();
}

//batch update for partition graph of PH2H, edge insertion
void Graph::EdgeInsertPartiBatchH2HPost(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<vertex,int>>& Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU){
    map<int,int> checkedDis;
    for(int i=0;i<Tree.size();i++){
        Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
    }

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<OrderCompp>> SCre; //SCre.clear();
    SCre.assign(node_num,set<OrderCompp>());//{vertexID, set<int>}
    set<OrderCompMin> OC; //OC.clear();//vertexID in decreasing node order

    int a,b,newW,lid,hid,lidM,hidM;
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first; b=wBatch[k].first.second ;newW=wBatch[k].second.second;
        if(NodeOrder[a]<NodeOrder[b]){
            lid=a;hid=b;
        }else{
            lid=b;hid=a;
        }
        lidM = IDMap[lid]; hidM = IDMap[hid];

//        if(pid==23){
//            cout<<pid<<": "<<a<<" "<<b<<" "<<newW<<endl;
//        }


        if(Neighbors[a].find(b)!=Neighbors[a].end()){
            cout<<"Already exist this edge. "<<a<<" "<<b<<" "<<Neighbors[a][b]<<endl;
            exit(1);
        }


        if(Neighbors[b].find(a)!=Neighbors[b].end()){
            cout<<"Already exist this edge. "<<a<<" "<<b<<" "<<Neighbors[b][a]<<endl;
            exit(1);
        }


        Neighbors[a].insert({b,newW});
        Neighbors[b].insert({a,newW});
        bool ifFind=false;
        for(int i=0;i<Tree[rank[lidM]].vert.size();i++){
//            cout<<"vert "<<lid<<" "<<Tree[rank[lidM]].vert[i].first<<endl;
            if(Tree[rank[lidM]].vert[i].first==hid){
                ifFind=true;
                if(Tree[rank[lidM]].vert[i].second.first>newW){
                    Tree[rank[lidM]].vert[i].second.first=newW;
                    Tree[rank[lidM]].vert[i].second.second=1;
                    SCre[lid].insert(hid);
                    OC.insert(OrderCompMin(lid));
                }else if(Tree[rank[lidM]].vert[i].second.first==newW){
                    Tree[rank[lidM]].vert[i].second.second+=1;
                }
                break;
            }
        }
        int newSCNum=0;
        if(!ifFind){
//            cout<<"Not found edge "<<lid<<" "<<hid<<" "<<newW<<endl;
            ifFind=false;
            Tree[rank[lidM]].vert.emplace_back(hid, make_pair(newW,1));
            newSCNum++;
            SCre[lid].insert(hid);
            OC.insert(OrderCompMin(lid));
//            cout<<"New SC 1. "<<lid<<" "<<hid<<" "<<newW<<endl;
            if(algoChoice==H2H){
                for(int i=0;i<Tree[rank[lidM]].vAncestor.size();++i){
                    if(Tree[rank[lidM]].vAncestor[i]==hid){
                        int temp=Tree[rank[lidM]].pos.back();
                        Tree[rank[lidM]].pos[Tree[rank[lidM]].pos.size()-1]=i;
                        Tree[rank[lidM]].pos.push_back(temp);
                        ifFind=true;
                    }
                }
                if(!ifFind){//if false
                    cout<<"Cross-branch edge! Should reorganize the tree structure! "<<lid<<" "<<hid<<" "<<newW<<endl;
                    exit(1);
                }
                BottomUpNewShortcutInsertParti(lid,hid,newW,Tree,rank,newSCNum);
            }
            else{
                BottomUpNewShortcutInsertPartiCH(lid,hid,newW,Tree,rank,newSCNum);
            }

//            cout<<"newSCNum: "<<newSCNum<<endl;
        }
    }

//    vector<int> ProBeginVertexSet; //ProBeginVertexSet.clear();
//    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed
    ProBeginVertexSetParti[pid].clear(); vertexIDChLParti[pid].clear();
    vector<int> ProBeginVertexSetNew;
    int ProBeginVertexID;
    int ProID, ProIDM;
    //processing the stars
    while(!OC.empty()){
        ProID=(*OC.begin()).x;
        ProIDM=IDMap[ProID];
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProIDM]].vert;
        bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=it->x; int Cw; int CidM=IDMap[Cid];
            int cidH=Tree[rank[CidM]].height-1;

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }else{
                    Cw=Vert[j].second.first;
                }
            }

//            if(ProID==1204){
//                cout<<ProID<<" "<<Cid<<" "<<Cw<<endl;
//            }

            if(algoChoice==H2H){
                if(Tree[rank[ProIDM]].dis[cidH]>Cw){
                    Tree[rank[ProIDM]].dis[cidH]=Cw;
                    Tree[rank[ProIDM]].FN[cidH]=true;
                    ProIDdisCha=true;
                    vUpdated[ProID]=true;
                    Tree[rank[ProIDM]].DisRe.insert(Cid);
                    Tree[rank[ProIDM]].cnt[cidH]=1;//new
                }else if(Tree[rank[ProIDM]].dis[cidH]==Cw){
                    Tree[rank[ProIDM]].FN[cidH]=true;
                    Tree[rank[ProIDM]].cnt[cidH]+=1;//new
                }
            }


            int hid2,hidHeight2,lid2,lidHeight2,wsum,lid2M;
            for(int j=0;j<Tree[rank[CidM]].vert.size();j++){
                hid2=Tree[rank[CidM]].vert[j].first;hidHeight2=Tree[rank[IDMap[hid2]]].height-1;
                if(Hnei.find(hid2)!=Hnei.end()){
                    wsum=Cw+Hnei[hid2];
                    if(wsum<Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.first=wsum;
                        Tree[rank[CidM]].vert[j].second.second=1;
                        SCre[Cid].insert(hid2);
                        OC.insert(OrderCompMin(Cid));
                    }else if(wsum==Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.second+=1;
                    }

                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid2=Lnei[j].first; lid2M=IDMap[lid2];
                lidHeight2=Tree[rank[lid2M]].height-1;
                for(int k=0;k<Tree[rank[lid2M]].vert.size();k++){
                    if(Tree[rank[lid2M]].vert[k].first==Cid){
                        wsum=Cw+Lnei[j].second;
                        if(Tree[rank[lid2M]].vert[k].second.first>wsum){
                            Tree[rank[lid2M]].vert[k].second.first=wsum;
                            Tree[rank[lid2M]].vert[k].second.second=1;
                            SCre[lid2].insert(Cid);
                            OC.insert(OrderCompMin(lid2));
                        }else if(Tree[rank[lid2M]].vert[k].second.first==wsum){
                            if(SCre[ProID].find(OrderCompp(lid2))==SCre[ProID].end()) {//if not found, avoid repeated count
                                Tree[rank[lid2M]].vert[k].second.second += 1;
                            }
                        }

                        break;
                    }
                }
            }
        }

        if(ProIDdisCha){//if the distance labeling is dectected changed
            vertexIDChLParti[pid].insert(ProID);
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetParti[pid].size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[IDMap[ProID]],r;
            for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
                r=rank[IDMap[ProBeginVertexSetParti[pid][i]]];
                if(LCAQueryPartition(rnew,r,pid)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetParti[pid][i]);
                }
            }
            ProBeginVertexSetParti[pid]=ProBeginVertexSetNew;
        }
    }

    //cout<<"Finish bottom-up refresh"<<endl;
    if(ifLabelU){
        for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
            ProBeginVertexID=ProBeginVertexSetParti[pid][i];
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[IDMap[ProBeginVertexID]]].pa].uniqueVertex;
            while(Tree[rank[IDMap[pachidd]]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[IDMap[pachidd]]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);
            EachNodeProBDis5Parti(rank[IDMap[ProBeginVertexID]], linee, vertexIDChLParti[pid],Tree,rank);
        }
    }

    //return checkedDis.size();
}

//batch update for partition graph of PH2H
void Graph::DecreasePartiBatch(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>>& Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<pair<pair<int,int>,int>>& updatedSC, bool ifLabelU){
    map<int,int> checkedDis;

    for(int i=0;i<Tree.size();i++){
        Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
    }

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<OrderCompp>> SCre; //SCre.clear();
    SCre.assign(node_num,set<OrderCompp>());//{vertexID, set<int>}
    set<OrderCompMin> OC; //OC.clear();//vertexID in decreasing node order

    int a,b,newW,lid,hid,lidM,hidM;
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first; b=wBatch[k].first.second ;newW=wBatch[k].second.second;
        if(NodeOrder[a]<NodeOrder[b]){
            lid=a;hid=b;
        }else{
            lid=b;hid=a;
        }
        lidM = IDMap[lid]; hidM = IDMap[hid];

//        if(pid==23){
//            cout<<pid<<": "<<a<<" "<<b<<" "<<newW<<endl;
//        }

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

        for(int i=0;i<Tree[rank[lidM]].vert.size();i++){
            if(Tree[rank[lidM]].vert[i].first==hid){
                if(Tree[rank[lidM]].vert[i].second.first>newW){
                    Tree[rank[lidM]].vert[i].second.first=newW;
                    Tree[rank[lidM]].vert[i].second.second=1;
                    SCre[lid].insert(hid);
                    OC.insert(OrderCompMin(lid));
                }else if(Tree[rank[lidM]].vert[i].second.first==newW){
                    Tree[rank[lidM]].vert[i].second.second+=1;
                }
                break;
            }
        }

    }

//    vector<int> ProBeginVertexSet; //ProBeginVertexSet.clear();
//    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed
    ProBeginVertexSetParti[pid].clear(); vertexIDChLParti[pid].clear();
    vector<int> ProBeginVertexSetNew;
    int ProBeginVertexID;
    int ProID, ProIDM;
    //processing the stars
    while(!OC.empty()){
        ProID=(*OC.begin()).x;
        ProIDM=IDMap[ProID];
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProIDM]].vert;
        bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=it->x; int Cw; int CidM=IDMap[Cid];
            int cidH=Tree[rank[CidM]].height-1;

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }else{
                    Cw=Vert[j].second.first;
                }
            }

            if(PartiTag[ProID].second){//if boundary vertex
                updatedSC.emplace_back(make_pair(ProID,Cid),Cw);
            }

            if(algoChoice==H2H){
                if(Tree[rank[ProIDM]].dis[cidH]>Cw){
                    Tree[rank[ProIDM]].dis[cidH]=Cw;
                    Tree[rank[ProIDM]].FN[cidH]=true;
                    ProIDdisCha=true;
                    vUpdated[ProID]=true;
                    Tree[rank[ProIDM]].DisRe.insert(Cid);
                    Tree[rank[ProIDM]].cnt[cidH]=1;//new
                }else if(Tree[rank[ProIDM]].dis[cidH]==Cw){
                    Tree[rank[ProIDM]].FN[cidH]=true;
                    Tree[rank[ProIDM]].cnt[cidH]+=1;//new
                }
            }


            int hid2,hidHeight2,lid2,lidHeight2,wsum,lid2M;
            for(int j=0;j<Tree[rank[CidM]].vert.size();j++){
                hid2=Tree[rank[CidM]].vert[j].first;hidHeight2=Tree[rank[IDMap[hid2]]].height-1;
                if(Hnei.find(hid2)!=Hnei.end()){
                    wsum=Cw+Hnei[hid2];
                    if(wsum<Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.first=wsum;
                        Tree[rank[CidM]].vert[j].second.second=1;
                        SCre[Cid].insert(hid2);
                        OC.insert(OrderCompMin(Cid));
                    }else if(wsum==Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.second+=1;
                    }

                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid2=Lnei[j].first; lid2M=IDMap[lid2];
                lidHeight2=Tree[rank[lid2M]].height-1;
                for(int k=0;k<Tree[rank[lid2M]].vert.size();k++){
                    if(Tree[rank[lid2M]].vert[k].first==Cid){
                        wsum=Cw+Lnei[j].second;
                        if(Tree[rank[lid2M]].vert[k].second.first>wsum){
                            Tree[rank[lid2M]].vert[k].second.first=wsum;
                            Tree[rank[lid2M]].vert[k].second.second=1;
                            SCre[lid2].insert(Cid);
                            OC.insert(OrderCompMin(lid2));
                        }else if(Tree[rank[lid2M]].vert[k].second.first==wsum){
                            if(SCre[ProID].find(OrderCompp(lid2))==SCre[ProID].end()) {//if not found, avoid repeated count
                                Tree[rank[lid2M]].vert[k].second.second += 1;
                            }
                        }

                        break;
                    }
                }
            }
        }

        if(ProIDdisCha){//if the distance labeling is dectected changed
            vertexIDChLParti[pid].insert(ProID);
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetParti[pid].size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[IDMap[ProID]],r;
            for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
                r=rank[IDMap[ProBeginVertexSetParti[pid][i]]];
                if(LCAQueryPartition(rnew,r,pid)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetParti[pid][i]);
                }
            }
            ProBeginVertexSetParti[pid]=ProBeginVertexSetNew;
        }
    }

    //cout<<"Finish bottom-up refresh"<<endl;
    if(ifLabelU){
        for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
            ProBeginVertexID=ProBeginVertexSetParti[pid][i];
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[IDMap[ProBeginVertexID]]].pa].uniqueVertex;
            while(Tree[rank[IDMap[pachidd]]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[IDMap[pachidd]]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);
            EachNodeProBDis5Parti(rank[IDMap[ProBeginVertexID]], linee, vertexIDChLParti[pid],Tree,rank);
        }
    }

    //return checkedDis.size();
}

void Graph::DecreasePartiBatchForOpt(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>>& Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<pair<pair<int,int>,int>>& updatedSC, bool ifLabelU, bool ifConstruct){
    map<int,int> checkedDis;

//    for(int i=0;i<Tree.size();i++){
//        Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
//    }

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<OrderCompp>> SCre; //SCre.clear();
    SCre.assign(node_num,set<OrderCompp>());//{vertexID, set<int>}
    set<OrderCompMin> OC; //OC.clear();//vertexID in decreasing node order

    int a,b,oldW,newW,lid,hid,lidM,hidM;
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first; b=wBatch[k].first.second ;newW=wBatch[k].second.second;oldW=wBatch[k].second.first;
        if(NodeOrder[a]<NodeOrder[b]){
            lid=a;hid=b;
        }else{
            lid=b;hid=a;
        }
        lidM = IDMap[lid]; hidM = IDMap[hid];

        for(int i=0;i<Neighbors[a].size();i++){
            if(Neighbors[a][i].first==b){
                if(Neighbors[a][i].second!=oldW){
                    cout<<"Old edge weight inconsistent! "<<a<<" "<<b<<" "<<oldW<<" "<<Neighbors[a][i].second<<endl; exit(1);
                }
                Neighbors[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbors[b].size();i++){
            if(Neighbors[b][i].first==a){
                if(Neighbors[b][i].second!=oldW){
                    cout<<"Old edge weight inconsistent! "<<b<<" "<<a<<" "<<oldW<<" "<<Neighbors[b][i].second<<endl; exit(1);
                }
                Neighbors[b][i].second=newW;
                break;
            }
        }

        for(int i=0;i<Tree[rank[lidM]].vert.size();i++){
            if(Tree[rank[lidM]].vert[i].first==hid){
                if(Tree[rank[lidM]].vert[i].second.first>newW){
                    Tree[rank[lidM]].vert[i].second.first=newW;
                    Tree[rank[lidM]].vert[i].second.second=1;
                    SCre[lid].insert(hid);
                    OC.insert(OrderCompMin(lid));
                }else if(Tree[rank[lidM]].vert[i].second.first==newW){
//                    if(pid==12){
//                        cout<<"Equal. "<<pid<<": "<<lid<<" "<<hid<<" "<<oldW<<" "<<newW<<" ("<<Tree[rank[lidM]].vert[i].second.second<<")"<<endl;
//                    }
                    if(!ifConstruct){//should not plus one if in index construction
//                        cout<<"Equal. "<<pid<<": "<<lid<<" "<<hid<<" "<<oldW<<" "<<newW<<" ("<<Tree[rank[lidM]].vert[i].second.second<<")"<<endl;
                        Tree[rank[lidM]].vert[i].second.second+=1;
                    }
                }else{
                    cout<<"Unexpected result! "<<Tree[rank[lidM]].vert[i].second.first<<" "<<newW<<endl; exit(1);
                }
                break;
            }
        }

    }

//    if(pid==12){
//        cout<<pid<<endl;
//    }

//    vector<int> ProBeginVertexSet; //ProBeginVertexSet.clear();
//    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed
    if(!ifLabelU){//for PH2H_Post update
        ProBeginVertexSetParti[pid].clear(); vertexIDChLParti[pid].clear();
        for(int i=0;i<Tree.size();i++){
            Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
        }
    }

    vector<int> ProBeginVertexSetNew;
    int ProBeginVertexID;
    int ProID, ProIDM;
    //processing the stars
    while(!OC.empty()){
        ProID=(*OC.begin()).x;
        ProIDM=IDMap[ProID];
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProIDM]].vert;
        bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=it->x; int Cw; int CidM=IDMap[Cid];
            int cidH=Tree[rank[CidM]].height-1;

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }else{
                    Cw=Vert[j].second.first;
                }
            }

            if(PartiTag[ProID].second){//if boundary vertex
                updatedSC.emplace_back(make_pair(ProID,Cid),Cw);
            }

//            if(algoUpdate!=PCH_No && !ifConstruct){
//                if(Tree[rank[ProIDM]].dis[cidH]>Cw){
//                    Tree[rank[ProIDM]].dis[cidH]=Cw;
//                    Tree[rank[ProIDM]].FN[cidH]=true;
//                    ProIDdisCha=true;
//                    vUpdated[ProID]=true;
//                    Tree[rank[ProIDM]].DisRe.insert(Cid);
//                }else if(Tree[rank[ProIDM]].dis[cidH]==Cw){
//                    Tree[rank[ProIDM]].FN[cidH]=true;
//                }
//            }


            int hid2,hidHeight2,lid2,lidHeight2,wsum,lid2M;
            for(int j=0;j<Tree[rank[CidM]].vert.size();j++){
                hid2=Tree[rank[CidM]].vert[j].first;hidHeight2=Tree[rank[IDMap[hid2]]].height-1;
                if(Hnei.find(hid2)!=Hnei.end()){
                    wsum=Cw+Hnei[hid2];
                    if(wsum<Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.first=wsum;
                        Tree[rank[CidM]].vert[j].second.second=1;
                        SCre[Cid].insert(hid2);
                        OC.insert(OrderCompMin(Cid));
                    }else if(wsum==Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.second+=1;
                    }

                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid2=Lnei[j].first; lid2M=IDMap[lid2];
                lidHeight2=Tree[rank[lid2M]].height-1;
                for(int k=0;k<Tree[rank[lid2M]].vert.size();k++){
                    if(Tree[rank[lid2M]].vert[k].first==Cid){
                        wsum=Cw+Lnei[j].second;
                        if(Tree[rank[lid2M]].vert[k].second.first>wsum){
                            Tree[rank[lid2M]].vert[k].second.first=wsum;
                            Tree[rank[lid2M]].vert[k].second.second=1;
                            SCre[lid2].insert(Cid);
                            OC.insert(OrderCompMin(lid2));
                        }else if(Tree[rank[lid2M]].vert[k].second.first==wsum){
                            if(SCre[ProID].find(OrderCompp(lid2))==SCre[ProID].end()) {//if not found, avoid repeated count
                                Tree[rank[lid2M]].vert[k].second.second += 1;
                            }
                        }

                        break;
                    }
                }
            }
        }

        if(ProIDdisCha){//if the distance labeling is dectected changed
            vertexIDChLParti[pid].insert(ProID);
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetParti[pid].size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[IDMap[ProID]],r;
            for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
                r=rank[IDMap[ProBeginVertexSetParti[pid][i]]];
                if(LCAQueryPartition(rnew,r,pid)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetParti[pid][i]);
                }
            }
            ProBeginVertexSetParti[pid]=ProBeginVertexSetNew;
        }
    }

    //cout<<"Finish bottom-up refresh"<<endl;

    if(ifLabelU){
//        cout<<pid<<": ProBeginVertexSet size: "<<ProBeginVertexSetParti[pid].size()<<endl;
        for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
            ProBeginVertexID=ProBeginVertexSetParti[pid][i];
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[IDMap[ProBeginVertexID]]].pa].uniqueVertex;
            while(Tree[rank[IDMap[pachidd]]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[IDMap[pachidd]]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);
            EachNodeProBDis5Parti(rank[IDMap[ProBeginVertexID]], linee, vertexIDChLParti[pid],Tree,rank);
        }
    }

    //return checkedDis.size();
}


void Graph::DecreasePartiBatchLabel(vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int>& ProBeginVertexSet, set<int>& vertexIDChL){
    int ProBeginVertexID;
    for(int i=0;i<ProBeginVertexSet.size();i++){
        ProBeginVertexID=ProBeginVertexSet[i];
        vector<int> linee; //linee.clear();
        linee.reserve(heightMax);
        int pachidd=Tree[Tree[rank[IDMap[ProBeginVertexID]]].pa].uniqueVertex;
        while(Tree[rank[IDMap[pachidd]]].height>1){
            linee.insert(linee.begin(),pachidd);
            pachidd=Tree[Tree[rank[IDMap[pachidd]]].pa].uniqueVertex;
        }
        linee.insert(linee.begin(),pachidd);
        EachNodeProBDis5Parti(rank[IDMap[ProBeginVertexID]], linee, vertexIDChL,Tree,rank);
    }
}

//batch update for partition graph
void Graph::DecreasePartiBatchPost(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>>& Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU){
    map<int,int> checkedDis;

    for(int i=0;i<Tree.size();i++){
        Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
    }

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<OrderCompp>> SCre; //SCre.clear();
    SCre.assign(node_num,set<OrderCompp>());//{vertexID, set<int>}
    set<OrderCompMin> OC; //OC.clear();//vertexID in decreasing node order

    int a,b,newW,lid,hid,lidM,hidM;
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first; b=wBatch[k].first.second ;newW=wBatch[k].second.second;
        if(NodeOrder[a]<NodeOrder[b]){
            lid=a;hid=b;
        }else{
            lid=b;hid=a;
        }
        lidM = IDMap[lid]; hidM = IDMap[hid];

        if(Neighbors[a].find(b)!=Neighbors[a].end()){
            Neighbors[a][b]=newW;
        }else{
            cout<<"Not found edge! "<<a<<" "<<b<<" "<<newW<<endl; exit(1);
        }

        if(Neighbors[b].find(a)!=Neighbors[b].end()){
            Neighbors[b][a]=newW;
        }else{
            cout<<"Not found edge! "<<b<<" "<<a<<" "<<newW<<endl; exit(1);
        }

        for(int i=0;i<Tree[rank[lidM]].vert.size();i++){
            if(Tree[rank[lidM]].vert[i].first==hid){
                if(Tree[rank[lidM]].vert[i].second.first>newW){
                    Tree[rank[lidM]].vert[i].second.first=newW;
                    Tree[rank[lidM]].vert[i].second.second=1;
                    SCre[lid].insert(hid);
                    OC.insert(OrderCompMin(lid));
                }else if(Tree[rank[lidM]].vert[i].second.first==newW){
                    Tree[rank[lidM]].vert[i].second.second+=1;
                }
                break;
            }
        }

    }

//    vector<int> ProBeginVertexSet; //ProBeginVertexSet.clear();
//    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distance labeling has changed
    ProBeginVertexSetParti[pid].clear(); vertexIDChLParti[pid].clear();
    vector<int> ProBeginVertexSetNew;
    int ProBeginVertexID;
    int ProID, ProIDM;
    //processing the stars
    while(!OC.empty()){
        ProID=(*OC.begin()).x;
        ProIDM=IDMap[ProID];
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProIDM]].vert;
        bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=it->x; int Cw; int CidM=IDMap[Cid];
            int cidH=Tree[rank[CidM]].height-1;

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }else{
                    Cw=Vert[j].second.first;
                }
            }

            if(algoChoice==H2H){
                if(Tree[rank[ProIDM]].dis[cidH]>Cw){
                    Tree[rank[ProIDM]].dis[cidH]=Cw;
                    Tree[rank[ProIDM]].FN[cidH]=true;
                    ProIDdisCha=true;
                    Tree[rank[ProIDM]].DisRe.insert(Cid);
                    Tree[rank[ProIDM]].cnt[cidH]=1;//new
                }else if(Tree[rank[ProIDM]].dis[cidH]==Cw){
                    Tree[rank[ProIDM]].FN[cidH]=true;
                    Tree[rank[ProIDM]].cnt[cidH]+=1;//new
                }
            }


            int hid2,hidHeight2,lid2,lidHeight2,wsum,lid2M;
            for(int j=0;j<Tree[rank[CidM]].vert.size();j++){
                hid2=Tree[rank[CidM]].vert[j].first;hidHeight2=Tree[rank[IDMap[hid2]]].height-1;
                if(Hnei.find(hid2)!=Hnei.end()){
                    wsum=Cw+Hnei[hid2];
                    if(wsum<Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.first=wsum;
                        Tree[rank[CidM]].vert[j].second.second=1;
                        SCre[Cid].insert(hid2);
                        OC.insert(OrderCompMin(Cid));
                    }else if(wsum==Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.second+=1;
                    }

                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid2=Lnei[j].first; lid2M=IDMap[lid2];
                lidHeight2=Tree[rank[lid2M]].height-1;
                for(int k=0;k<Tree[rank[lid2M]].vert.size();k++){
                    if(Tree[rank[lid2M]].vert[k].first==Cid){
                        wsum=Cw+Lnei[j].second;
                        if(Tree[rank[lid2M]].vert[k].second.first>wsum){
                            Tree[rank[lid2M]].vert[k].second.first=wsum;
                            Tree[rank[lid2M]].vert[k].second.second=1;
                            SCre[lid2].insert(Cid);
                            OC.insert(OrderCompMin(lid2));
                        }else if(Tree[rank[lid2M]].vert[k].second.first==wsum){
                            if(SCre[ProID].find(OrderCompp(lid2))==SCre[ProID].end()) {//if not found, avoid repeated count
                                Tree[rank[lid2M]].vert[k].second.second += 1;
                            }
                        }

                        break;
                    }
                }
            }
        }

        if(ProIDdisCha){//if the distance labeling is dectected changed
            vertexIDChLParti[pid].insert(ProID);
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetParti[pid].size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[IDMap[ProID]],r;
            for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
                r=rank[IDMap[ProBeginVertexSetParti[pid][i]]];
                if(LCAQueryPartition(rnew,r,pid)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetParti[pid][i]);
                }
            }
            ProBeginVertexSetParti[pid]=ProBeginVertexSetNew;
        }
    }

    //cout<<"Finish bottom-up refresh"<<endl;
    if(ifLabelU){
        for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
            ProBeginVertexID=ProBeginVertexSetParti[pid][i];
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[IDMap[ProBeginVertexID]]].pa].uniqueVertex;
            while(Tree[rank[IDMap[pachidd]]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[IDMap[pachidd]]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);
            EachNodeProBDis5Parti(rank[IDMap[ProBeginVertexID]], linee, vertexIDChLParti[pid],Tree,rank);
        }
    }

    //return checkedDis.size();
}

//batch update for partition graph
void Graph::DecreasePartiBatchPre(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>>& Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU){
    map<int,int> checkedDis;

    for(int i=0;i<Tree.size();i++){
        Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
    }

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<OrderCompp>> SCre; //SCre.clear();
    SCre.assign(node_num,set<OrderCompp>());//{vertexID, set<int>}
    set<OrderCompMin> OC; //OC.clear();//vertexID in decreasing node order

    int a,b,newW,lid,hid,lidM,hidM;
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first; b=wBatch[k].first.second ;newW=wBatch[k].second.second;
        if(NodeOrder[a]<NodeOrder[b]){
            lid=a;hid=b;
        }else{
            lid=b;hid=a;
        }
        lidM = IDMap[lid]; hidM = IDMap[hid];

//        if(pid==23){
//            cout<<pid<<": "<<a<<" "<<b<<" "<<newW<<endl;
//        }

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

        for(int i=0;i<Tree[rank[lidM]].vert.size();i++){
            if(Tree[rank[lidM]].vert[i].first==hid){
                if(Tree[rank[lidM]].vert[i].second.first>newW){
                    Tree[rank[lidM]].vert[i].second.first=newW;
                    Tree[rank[lidM]].vert[i].second.second=1;
                    SCre[lid].insert(hid);
                    OC.insert(OrderCompMin(lid));
                }else if(Tree[rank[lidM]].vert[i].second.first==newW){
                    Tree[rank[lidM]].vert[i].second.second+=1;
                }
                break;
            }
        }

    }

//    vector<int> ProBeginVertexSet; //ProBeginVertexSet.clear();
//    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distance labeling has changed
    ProBeginVertexSetParti[pid].clear(); vertexIDChLParti[pid].clear();
    vector<int> ProBeginVertexSetNew;
    int ProBeginVertexID;
    int ProID, ProIDM;
    //processing the stars
    while(!OC.empty()){
        ProID=(*OC.begin()).x;
        ProIDM=IDMap[ProID];
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProIDM]].vert;
        bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=it->x; int Cw; int CidM=IDMap[Cid];
            int cidH=Tree[rank[CidM]].height-1;

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }else{
                    Cw=Vert[j].second.first;
                }
            }

            if(algoChoice==H2H){
                if(Tree[rank[ProIDM]].dis[cidH]>Cw){
                    Tree[rank[ProIDM]].dis[cidH]=Cw;
                    Tree[rank[ProIDM]].FN[cidH]=true;
                    ProIDdisCha=true;
                    Tree[rank[ProIDM]].DisRe.insert(Cid);
                    Tree[rank[ProIDM]].cnt[cidH]=1;//new
                }else if(Tree[rank[ProIDM]].dis[cidH]==Cw){
                    Tree[rank[ProIDM]].FN[cidH]=true;
                    Tree[rank[ProIDM]].cnt[cidH]+=1;//new
                }
            }


            int hid2,hidHeight2,lid2,lidHeight2,wsum,lid2M;
            for(int j=0;j<Tree[rank[CidM]].vert.size();j++){
                hid2=Tree[rank[CidM]].vert[j].first;hidHeight2=Tree[rank[IDMap[hid2]]].height-1;
                if(Hnei.find(hid2)!=Hnei.end()){
                    wsum=Cw+Hnei[hid2];
                    if(wsum<Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.first=wsum;
                        Tree[rank[CidM]].vert[j].second.second=1;
                        SCre[Cid].insert(hid2);
                        OC.insert(OrderCompMin(Cid));
                    }else if(wsum==Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.second+=1;
                    }

                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid2=Lnei[j].first; lid2M=IDMap[lid2];
                lidHeight2=Tree[rank[lid2M]].height-1;
                for(int k=0;k<Tree[rank[lid2M]].vert.size();k++){
                    if(Tree[rank[lid2M]].vert[k].first==Cid){
                        wsum=Cw+Lnei[j].second;
                        if(Tree[rank[lid2M]].vert[k].second.first>wsum){
                            Tree[rank[lid2M]].vert[k].second.first=wsum;
                            Tree[rank[lid2M]].vert[k].second.second=1;
                            SCre[lid2].insert(Cid);
                            OC.insert(OrderCompMin(lid2));
                        }else if(Tree[rank[lid2M]].vert[k].second.first==wsum){
                            if(SCre[ProID].find(OrderCompp(lid2))==SCre[ProID].end()) {//if not found, avoid repeated count
                                Tree[rank[lid2M]].vert[k].second.second += 1;
                            }
                        }

                        break;
                    }
                }
            }
        }

        if(ProIDdisCha){//if the distance labeling is dectected changed
            vertexIDChLParti[pid].insert(ProID);
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetParti[pid].size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[IDMap[ProID]],r;
            for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
                r=rank[IDMap[ProBeginVertexSetParti[pid][i]]];
                if(LCAQueryPartition(rnew,r,pid)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetParti[pid][i]);
                }
            }
            ProBeginVertexSetParti[pid]=ProBeginVertexSetNew;
        }
    }

    //cout<<"Finish bottom-up refresh"<<endl;
    if(ifLabelU){
        for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
            ProBeginVertexID=ProBeginVertexSetParti[pid][i];
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[IDMap[ProBeginVertexID]]].pa].uniqueVertex;
            while(Tree[rank[IDMap[pachidd]]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[IDMap[pachidd]]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);
            EachNodeProBDis5Parti(rank[IDMap[ProBeginVertexID]], linee, vertexIDChLParti[pid],Tree,rank);
        }
    }

    //return checkedDis.size();
}

void Graph::IncreaseOverlay(int a,int b, int oldW, int newW, vector<unordered_map<int,int>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid){
    int ChangeNum=0;
    map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
    //OCdis.clear();

    if(Neighbors[a].find(b)!=Neighbors[a].end()){
        Neighbors[a][b]=newW;
    }else{
        cout<<"Wrong for Neighbors!"<<endl; exit(1);
    }
    if(Neighbors[b].find(a)!=Neighbors[b].end()){
        Neighbors[b][a]=newW;
    }else{
        cout<<"Wrong for Neighbors!"<<endl; exit(1);
    }

    int lid,hid;
    if(NodeOrder[a]<NodeOrder[b]){
        lid=a;hid=b;
    }else{
        lid=b;hid=a;
    }

    int IniH=Tree[rank[lid]].height;//the height where weight change begins
    int ProH=Tree[rank[lid]].height; int ProID=lid;
    vector<set<int>> SCre;//record the shortcut change in each height
    set<int> vec; //vec.clear();
    SCre.assign(IniH+1,vec);
    int MinH;

    vector<int> line; //line.clear();
    line.reserve(heightMax);
    int pachid=ProID;
    while(Tree[rank[pachid]].height>1){
        line.insert(line.begin(),pachid);
        pachid=Tree[Tree[rank[pachid]].pa].uniqueVertex;
    }
    line.insert(line.begin(),pachid);

    bool tri=false;
    for(int i=0;i<Tree[rank[lid]].vert.size();i++){
        if(Tree[rank[lid]].vert[i].first==hid){
            if(Tree[rank[lid]].vert[i].second.first==oldW){
                Tree[rank[lid]].vert[i].second.second-=1;
                if(Tree[rank[lid]].vert[i].second.second<1){
                    OCdis[make_pair(lid,hid)]=oldW;
                    SCre[ProH].insert(hid);
                    MinH=IniH;
                    tri=true;//cout<<"Trigger the Shortcut Change or not? "<<tri<<endl;
                }
            }
            break;
        }
    }

    bool influence; int ProBeginID;
    if(tri){
        //shortcut update
        while(ProH>=MinH){
            influence=false;
            vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
            for(auto it=SCre[ProH].begin();it!=SCre[ProH].end();it++){
                int Cid=*it; int Cw=OCdis[make_pair(ProID,Cid)];
                int cidH=Tree[rank[Cid]].height-1;

                map<int,int> Hnei; //Hnei.clear();
                vector<pair<int,int>> Lnei; //Lnei.clear();
                for(int j=0;j<Vert.size();j++){
                    if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                        Hnei[Vert[j].first]=Vert[j].second.first;
                    }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                        Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                    }
                }
                //check the affected shortcuts
                int hid,lid;
                for(int j=0;j<Tree[rank[Cid]].vert.size();j++){
                    hid=Tree[rank[Cid]].vert[j].first;
                    if(Hnei.find(hid)!=Hnei.end()){
                        if(Cw+Hnei[hid]==Tree[rank[Cid]].vert[j].second.first){
                            Tree[rank[Cid]].vert[j].second.second-=1;
                            if(Tree[rank[Cid]].vert[j].second.second<1){
                                SCre[Tree[rank[Cid]].height].insert(hid);
                                if(Tree[rank[Cid]].height<MinH) MinH=Tree[rank[Cid]].height;
                                OCdis[make_pair(Cid,hid)]=Cw+Hnei[hid];
                            }
                        }
                    }
                }
                for(int j=0;j<Lnei.size();j++){
                    lid=Lnei[j].first;
                    for(int k=0;k<Tree[rank[lid]].vert.size();k++){
                        if(Tree[rank[lid]].vert[k].first==Cid){
                            if(Tree[rank[lid]].vert[k].second.first==Cw+Lnei[j].second){
                                Tree[rank[lid]].vert[k].second.second-=1;
                                if(Tree[rank[lid]].vert[k].second.second<1){
                                    SCre[Tree[rank[lid]].height].insert(Cid);
                                    if(Tree[rank[lid]].height<MinH) MinH=Tree[rank[lid]].height;
                                    OCdis[make_pair(lid,Cid)]=Cw+Lnei[j].second;
                                }
                            }
                            break;
                        }
                    }
                }

                //before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
                if(Tree[rank[ProID]].FN[cidH]){
                    influence=true;
                    //higher than Cid
                    for(int i=0;i<cidH;i++){
                        if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[Cid]].dis[i]){
                            Tree[rank[ProID]].cnt[i]-=1;
                        }
                    }

                    //equal to Cid
                    Tree[rank[ProID]].FN[cidH]=false;
                    Tree[rank[ProID]].cnt[cidH]-=1;

                    //lower than Cid
                    for(int i=cidH+1;i<Tree[rank[ProID]].dis.size();i++){
                        if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[line[i]]].dis[cidH]){
                            Tree[rank[ProID]].cnt[i]-=1;
                        }
                    }
                }

                //get the new value of shortcut
                //	cout<<Cw<<" increase to ";
                Cw=INF; int countwt=0;

                for(auto it2=Neighbors[ProID].begin();it2!=Neighbors[ProID].end();++it2){
                    if(it2->first==Cid){
                        Cw=it2->second;//the weight value in the original graph
                        countwt=1;
                        break;
                    }
                }

                int ssw,wtt,wid;
                vector<pair<int,int>> Wnodes; //Wnodes.clear();
                /*if(ProID<Cid)
                    Wnodes=SCconNodes[make_pair(ProID,Cid)]; //cout<<"wid num "<<Wnodes.size()<<endl;
                else
                    Wnodes=SCconNodes[make_pair(Cid,ProID)];*/

                if(ProID<Cid)
                    Wnodes=SCconNodesMT[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
                else
                    Wnodes=SCconNodesMT[Cid][ProID];
                for(int i=0;i<Wnodes.size();i++){
                    wid=Wnodes[i].first;
                    for(int j=0;j<Tree[rank[wid]].vert.size();j++){
                        if(Tree[rank[wid]].vert[j].first==ProID){
                            ssw=Tree[rank[wid]].vert[j].second.first;
                        }
                        if(Tree[rank[wid]].vert[j].first==Cid){
                            wtt=Tree[rank[wid]].vert[j].second.first;
                        }
                    }

                    if(ssw+wtt<Cw){
                        Cw=ssw+wtt;
                        countwt=1;
                    }else if(ssw+wtt==Cw){
                        countwt+=1;
                    }
                }

                //cout<<Cw<<endl;
                //refresh the shortcut to the new value
                for(int i=0;i<Tree[rank[ProID]].vert.size();i++){
                    if(Tree[rank[ProID]].vert[i].first==Cid){
                        Tree[rank[ProID]].vert[i].second.first=Cw;
                        Tree[rank[ProID]].vert[i].second.second=countwt;
                        break;
                    }
                }
            }

            if(influence){
                ProBeginID=ProID;
            }

            ProH-=1;
            ProID=Tree[Tree[rank[ProID]].pa].uniqueVertex;
        }
    }

    vector<int> line1; //line1.clear();
    line1.reserve(heightMax);
    pachid=Tree[Tree[rank[ProBeginID]].pa].uniqueVertex;
    while(Tree[rank[pachid]].height>1){
        line1.insert(line1.begin(),pachid);
        pachid=Tree[Tree[rank[pachid]].pa].uniqueVertex;
    }
    line1.insert(line1.begin(),pachid);

    eachNodeProcessIncrease1(rank[ProBeginID],line1,ChangeNum,Tree,rank,VidtoTNid);

    //return ChangeNum;
}

void Graph::eachNodeProcessIncrease1(int children, vector<int>& line, int& changelabel, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid){
    int childID=Tree[children].uniqueVertex;
    int childH=Tree[children].height-1;
    for(int i=0;i<Tree[children].dis.size();i++){
//        if(childID==6410 && Tree[children].vAncestor[i]==6411){
//            cout<<"Find 3. "<<childID<<" "<<Tree[children].vAncestor[i]<<" "<<Tree[children].dis[i]<<endl;
//        }
        if(Tree[children].cnt[i]<=0){
            vUpdated[childID] = true;
            changelabel+=1;
            //firstly, check which dis can be infected
            int disBF=Tree[children].dis[i];
            int PID;
            //chidlID
            for(int k=0;k<VidtoTNid[childID].size();k++){
                PID=VidtoTNid[childID][k];
                if(Tree[PID].FN[childH] && Tree[PID].dis[i]==disBF+Tree[PID].dis[childH]){
                    Tree[PID].cnt[i]-=1;
                }
            }

            //line[i]
            for(int k=0;k<VidtoTNid[line[i]].size();k++){
                PID=VidtoTNid[line[i]][k];
                if(Tree[PID].height>Tree[children].height && Tree[PID].vAncestor[childH] == childID){///
                    if(Tree[PID].FN[i] && Tree[PID].dis[childH]==disBF+Tree[PID].dis[i]){
                        Tree[PID].cnt[childH]-=1;
                    }
                }
            }

            //secondly, calculate the actual distance
            int dis=INF; int count=0;
            int Dvb; int b,bH; int DDvb=INF;
            for(int j=0;j<Tree[children].vert.size();j++){
                Dvb=Tree[children].vert[j].second.first;
                b=Tree[children].vert[j].first;
                bH=Tree[rank[b]].height-1;
                if(bH<i){
                    if(Dvb+Tree[rank[line[i]]].dis[bH]<dis){
                        dis=Dvb+Tree[rank[line[i]]].dis[bH];
                        count=1;
                    }else if(Dvb+Tree[rank[line[i]]].dis[bH]==dis){
                        count+=1;
                    }
                }else if(bH==i){
                    DDvb=Dvb;
                    if(Dvb<dis){
                        dis=Dvb;
                        count=1;
                    }else if(Dvb==dis){
                        count+=1;
                    }
                }else{
                    if(Dvb+Tree[rank[b]].dis[i]<dis){
                        dis=Dvb+Tree[rank[b]].dis[i];
                        count=1;
                    }else if(Dvb+Tree[rank[b]].dis[i]==dis){
                        count+=1;
                    }
                }
            }
            if(DDvb==dis) Tree[children].FN[i]=true;
            Tree[children].dis[i]=dis;
            Tree[children].cnt[i]=count;
        }
    }

    line.push_back(childID);
    for(int i=0;i<Tree[children].ch.size();i++){
        eachNodeProcessIncrease1(Tree[children].ch[i],line,changelabel,Tree,rank,VidtoTNid);
    }
    line.pop_back();
}

void Graph::IncreaseParti(int a,int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid){
    int ChangeNum=0;
    map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
    //OCdis.clear();

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

    int lid,hid;
    if(NodeOrder[a]<NodeOrder[b]){
        lid=a;hid=b;
    }else{
        lid=b;hid=a;
    }
    int lidM=IDMap[lid];
    int IniH=Tree[rank[lidM]].height;//the height where weight change begins
    int ProH=Tree[rank[lidM]].height; int ProID=lid;
    vector<set<int>> SCre;//record the shortcut change in each height
    set<int> vec; //vec.clear();
    SCre.assign(IniH+1,vec);
    int MinH;

    vector<int> line; //line.clear();
    line.reserve(heightMax);
    int pachid=ProID;
    while(Tree[rank[IDMap[pachid]]].height>1){
        line.insert(line.begin(),pachid);
        pachid=Tree[Tree[rank[IDMap[pachid]]].pa].uniqueVertex;
    }
    line.insert(line.begin(),pachid);

    bool tri=false;
    for(int i=0;i<Tree[rank[lidM]].vert.size();i++){
        if(Tree[rank[lidM]].vert[i].first==hid){
            if(Tree[rank[lidM]].vert[i].second.first==oldW){
                Tree[rank[lidM]].vert[i].second.second-=1;
                if(Tree[rank[lidM]].vert[i].second.second<1){
                    OCdis[make_pair(lid,hid)]=oldW;
                    SCre[ProH].insert(hid);
                    MinH=IniH;
                    tri=true;//cout<<"Trigger the Shortcut Change or not? "<<tri<<endl;
                }
            }
            break;
        }
    }

    bool influence; int ProBeginID;
    if(tri){
        //shortcut update
        while(ProH>=MinH){
            influence=false;
            int ProIDM=IDMap[ProID];
            vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProIDM]].vert;
            for(auto it=SCre[ProH].begin();it!=SCre[ProH].end();it++){
                int Cid=*it; int Cw=OCdis[make_pair(ProID,Cid)];
                int cidH=Tree[rank[IDMap[Cid]]].height-1;
                int CidM=IDMap[Cid];

                map<int,int> Hnei; //Hnei.clear();
                vector<pair<int,int>> Lnei; //Lnei.clear();
                for(int j=0;j<Vert.size();j++){
                    if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                        Hnei[Vert[j].first]=Vert[j].second.first;
                    }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                        Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                    }
                }
                //check the affected shortcuts
                int hid,lid;
                for(int j=0;j<Tree[rank[CidM]].vert.size();j++){// for higher-order vertices
                    hid=Tree[rank[CidM]].vert[j].first;
                    if(Hnei.find(hid)!=Hnei.end()){
                        if(Cw+Hnei[hid]==Tree[rank[CidM]].vert[j].second.first){
                            Tree[rank[CidM]].vert[j].second.second-=1;
                            if(Tree[rank[CidM]].vert[j].second.second<1){
                                SCre[Tree[rank[CidM]].height].insert(hid);
                                if(Tree[rank[CidM]].height<MinH) MinH=Tree[rank[CidM]].height;
                                OCdis[make_pair(Cid,hid)]=Cw+Hnei[hid];
                            }
                        }
                    }
                }
                for(int j=0;j<Lnei.size();j++){//for lower-order vertices
                    lid=Lnei[j].first; lidM=IDMap[lid];
                    for(int k=0;k<Tree[rank[lidM]].vert.size();k++){
                        if(Tree[rank[lidM]].vert[k].first==Cid){
                            if(Tree[rank[lidM]].vert[k].second.first==Cw+Lnei[j].second){
                                Tree[rank[lidM]].vert[k].second.second-=1;
                                if(Tree[rank[lidM]].vert[k].second.second<1){
                                    SCre[Tree[rank[lidM]].height].insert(Cid);
                                    if(Tree[rank[lidM]].height<MinH) MinH=Tree[rank[lidM]].height;
                                    OCdis[make_pair(lid,Cid)]=Cw+Lnei[j].second;
                                }
                            }
                            break;
                        }
                    }
                }

                //before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
                if(Tree[rank[ProIDM]].FN[cidH]){//if the label is obtained from shortcut
                    influence=true;
                    //higher than Cid
                    for(int i=0;i<cidH;i++){
                        if(Tree[rank[ProIDM]].dis[i]==Cw+Tree[rank[IDMap[Cid]]].dis[i]){
                            Tree[rank[ProIDM]].cnt[i]-=1;
                        }
                    }

                    //equal to Cid
                    Tree[rank[ProIDM]].FN[cidH]=false;
                    Tree[rank[ProIDM]].cnt[cidH]-=1;

                    //lower than Cid
                    for(int i=cidH+1;i<Tree[rank[ProIDM]].dis.size();i++){
                        if(Tree[rank[ProIDM]].dis[i]==Cw+Tree[rank[IDMap[line[i]]]].dis[cidH]){
                            Tree[rank[ProIDM]].cnt[i]-=1;
                        }
                    }
                }

                //get the new value of shortcut
                //	cout<<Cw<<" increase to ";
                Cw=INF; int countwt=0;

                for(int i=0;i<Neighbors[ProID].size();i++){
                    if(Neighbors[ProID][i].first==Cid){
                        Cw=Neighbors[ProID][i].second;//the weight value in the original graph
                        countwt=1;
                        break;
                    }
                }

                int ssw,wtt,wid;
                vector<pair<int,int>> Wnodes; //Wnodes.clear();
                /*if(ProID<Cid)
                    Wnodes=SCconNodes[make_pair(ProID,Cid)]; //cout<<"wid num "<<Wnodes.size()<<endl;
                else
                    Wnodes=SCconNodes[make_pair(Cid,ProID)];*/

                if(ProID<Cid)
                    Wnodes=SCconNodesMT[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
                else
                    Wnodes=SCconNodesMT[Cid][ProID];
                for(int i=0;i<Wnodes.size();i++){//for each supportive vertex of this shortcut
                    wid=Wnodes[i].first;
                    int widM=IDMap[wid];
                    for(int j=0;j<Tree[rank[widM]].vert.size();j++){
                        if(Tree[rank[widM]].vert[j].first==ProID){
                            ssw=Tree[rank[widM]].vert[j].second.first;
                        }
                        if(Tree[rank[widM]].vert[j].first==Cid){
                            wtt=Tree[rank[widM]].vert[j].second.first;
                        }
                    }

                    if(ssw+wtt<Cw){
                        Cw=ssw+wtt;
                        countwt=1;
                    }else if(ssw+wtt==Cw){
                        countwt+=1;
                    }
                }

                //cout<<Cw<<endl;
                //refresh the shortcut to the new value
                for(int i=0;i<Tree[rank[ProIDM]].vert.size();i++){
                    if(Tree[rank[ProIDM]].vert[i].first==Cid){
                        Tree[rank[ProIDM]].vert[i].second.first=Cw;
                        Tree[rank[ProIDM]].vert[i].second.second=countwt;
                        break;
                    }
                }
            }

            if(influence){
                ProBeginID=ProID;
            }

            ProH-=1;
            ProID=Tree[Tree[rank[ProIDM]].pa].uniqueVertex;
        }
    }

    vector<int> line1; //line1.clear();
    line1.reserve(heightMax);
    pachid=Tree[Tree[rank[IDMap[ProBeginID]]].pa].uniqueVertex;
    while(Tree[rank[IDMap[pachid]]].height>1){
        line1.insert(line1.begin(),pachid);
        pachid=Tree[Tree[rank[IDMap[pachid]]].pa].uniqueVertex;
    }
    line1.insert(line1.begin(),pachid);

    eachNodeProcessIncrease1Parti(rank[IDMap[ProBeginID]],line1,ChangeNum,Tree,rank,VidtoTNid);

    //return ChangeNum;
}

void Graph::eachNodeProcessIncrease1Parti(int children, vector<int>& line, int& changelabel, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid){
    int childID=Tree[children].uniqueVertex;
    int childH=Tree[children].height-1;
//    int childM=IDMap[children];

//    cout<<childID<<" "<<childH<<endl;
    for(int i=0;i<Tree[children].dis.size();i++){
//        if(childID==7988){
//            cout<<childID<<": "<<i<<", "<<Tree[children].vAncestor[i]<< " "<<Tree[children].cnt[i]<< endl;
//        }
//        if(childID==7988 && Tree[children].vAncestor[i]==7984){
//            cout<<"Find 4. "<<childID<<" "<<Tree[children].vAncestor[i]<<" "<<Tree[children].dis[i]<<endl;
//        }

        if(Tree[children].cnt[i]<=0){//the distance from child to line[i] should be updated
            if(i<BoundVertex[PartiTag[childID].first].size()){
                vUpdated[childID] = true;
            }

            changelabel+=1;
            //firstly, check which dis can be infected
            int disBF=Tree[children].dis[i];
            int tid;
            //chidlID
            for(int k=0;k<VidtoTNid[childID].size();k++){//check the tree nodes containing child
                tid=VidtoTNid[childID][k];
//                cout<<k<<" "<<tid<<" "<<IDMap[tid]<<" "<<childH<<" "<<Tree[tid].dis.size()<<" "<<Tree[tid].FN.size()<<endl;
                if(Tree[tid].FN[childH] && Tree[tid].dis[i]==disBF+Tree[tid].dis[childH]){//if the distance from tid to line[i] sources from child, valley path
                    Tree[tid].cnt[i]-=1;
                }
            }
//            cout<<"Flag 2"<<endl;
            //line[i]
            for(int k=0;k<VidtoTNid[line[i]].size();k++){
                tid=VidtoTNid[line[i]][k];
                if(Tree[tid].height>Tree[children].height && Tree[tid].vAncestor[childH] == childID){//children is the ancestor of tid
                    if(Tree[tid].FN[i] && Tree[tid].dis[childH]==disBF+Tree[tid].dis[i]){//if the distance from tid to child sources from line[i], peak path, out of scope
                        Tree[tid].cnt[childH]-=1;
                    }
                }
            }

            //secondly, calculate the actual distance
            int dis=INF; int count=0;
            int Dvb; int b,bH; int DDvb=INF;
            for(int j=0;j<Tree[children].vert.size();j++){
                Dvb=Tree[children].vert[j].second.first;
                b=Tree[children].vert[j].first;
                bH=Tree[rank[IDMap[b]]].height-1;
                if(bH<i){//if b is the ancestor of line[i]
                    if(Dvb+Tree[rank[IDMap[line[i]]]].dis[bH]<dis){
                        dis=Dvb+Tree[rank[IDMap[line[i]]]].dis[bH];
                        count=1;
                    }else if(Dvb+Tree[rank[IDMap[line[i]]]].dis[bH]==dis){
                        count+=1;
                    }
                }else if(bH==i){
                    DDvb=Dvb;//shortcut
                    if(Dvb<dis){
                        dis=Dvb;
                        count=1;
                    }else if(Dvb==dis){
                        count+=1;
                    }
                }else{//if line[i] is the ancestor of b
                    if(Dvb+Tree[rank[IDMap[b]]].dis[i]<dis){
                        dis=Dvb+Tree[rank[IDMap[b]]].dis[i];
                        count=1;
                    }else if(Dvb+Tree[rank[IDMap[b]]].dis[i]==dis){
                        count+=1;
                    }
                }
            }
            if(DDvb==dis) Tree[children].FN[i]=true;
            Tree[children].dis[i]=dis;
            Tree[children].cnt[i]=count;
        }
    }

    line.push_back(childID);
    for(int i=0;i<Tree[children].ch.size();i++){
        eachNodeProcessIncrease1Parti(Tree[children].ch[i],line,changelabel,Tree,rank,VidtoTNid);
    }
    line.pop_back();
}

void Graph::IncreaseOverlayBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, bool ifLabelU){
    int checknum=0;
    map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
    OCdis.clear();

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<OrderCompp>> SCre; //SCre.clear();
    SCre.assign(node_num,set<OrderCompp>());//{vertexID, set<int>}
    set<OrderCompMin> OC; OC.clear();//vertexID in decreasing node order
    bool flag=false;
    for(int k=0;k<wBatch.size();k++){
        int a=wBatch[k].first.first;
        int b=wBatch[k].first.second;
        int oldW=wBatch[k].second.first;
        int newW=wBatch[k].second.second;

        if(oldW!=newW){
            if(Neighbor[a].find(b)!=Neighbor[a].end()){
                if(Neighbor[a][b]!=oldW){//only works for no-boundary
                    cout<<"Inconsistent! "<<Neighbor[a][b]<<" "<<oldW<<endl;
                    exit(1);
                }
                Neighbor[a][b]=newW;
            }
//            else{
//                cout<<"Wrong for Neighbors!"<<endl; exit(1);
//            }
            if(Neighbor[b].find(a)!=Neighbor[b].end()){
                if(Neighbor[b][a]!=oldW){
                    cout<<"Inconsistent! "<<Neighbor[b][a]<<" "<<oldW<<endl;
                    exit(1);
                }
                Neighbor[b][a]=newW;
            }
//            else{
//                cout<<"Wrong for Neighbors!"<<endl; exit(1);
//            }

            int lid,hid;
            if(NodeOrder[a]<NodeOrder[b]){
                lid=a;hid=b;
            }else{
                lid=b;hid=a;
            }
//            if(lid==6410 && hid==6411){
//                cout<<"Find. "<<lid<<" "<<hid<<" "<<oldW<<" "<<newW<<endl;
//            }
            for(int i=0;i<Tree[rank[lid]].vert.size();i++){
                if(Tree[rank[lid]].vert[i].first==hid){
                    if(Tree[rank[lid]].vert[i].second.first==oldW){
                        Tree[rank[lid]].vert[i].second.second-=1;
                        if(Tree[rank[lid]].vert[i].second.second<1){
                            OCdis[make_pair(lid,hid)]=oldW;
                            SCre[lid].insert(hid);
                            OC.insert(OrderCompMin(lid));
                        }
                    }
//                    else{
//                        cout<<"oldW is incorrect. "<<lid<<" "<<hid<<" "<<oldW<<" "<<Tree[rank[lid]].vert[i].second.first<<endl; exit(1);
//                    }
                    break;
                }
            }
        }
    }

//    vector<int> ProBeginVertexSet; ProBeginVertexSet.clear();
    ProBeginVertexSetOverlay.clear();
    vector<int> ProBeginVertexSetNew;
    bool influence;
    int ProID; vector<int> line;
    while(!OC.empty()){
        ProID=(*OC.begin()).x;
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
        influence=false;

        //each ProID corresponds to a line
        line.clear(); line.reserve(heightMax);
        int pachid=ProID;
        while(Tree[rank[pachid]].height>1){
            line.insert(line.begin(),pachid);
            pachid=Tree[Tree[rank[pachid]].pa].uniqueVertex;
        }
        line.insert(line.begin(),pachid);

        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=it->x; int Cw=OCdis[make_pair(ProID,Cid)];
            int cidH=Tree[rank[Cid]].height-1;

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }
            }

//            if(ProID==6410 && Cid==6411){
//                cout<<"Find 2. "<<ProID<<" "<<Cid<<" "<<Cw<<endl;
//                flag=true;
//            }

            //check the affected shortcuts
            int hid,lid;
            for(int j=0;j<Tree[rank[Cid]].vert.size();j++){
                hid=Tree[rank[Cid]].vert[j].first;
                if(Hnei.find(hid)!=Hnei.end()){
                    if(Cw+Hnei[hid]==Tree[rank[Cid]].vert[j].second.first){
                        Tree[rank[Cid]].vert[j].second.second-=1;
                        if(Tree[rank[Cid]].vert[j].second.second<1){
                            SCre[Cid].insert(hid);
                            OC.insert(OrderCompMin(Cid));
                            OCdis[make_pair(Cid,hid)]=Cw+Hnei[hid];
                        }
                    }
                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid=Lnei[j].first;
                for(int k=0;k<Tree[rank[lid]].vert.size();k++){
                    if(Tree[rank[lid]].vert[k].first==Cid){
                        if(Tree[rank[lid]].vert[k].second.first==Cw+Lnei[j].second){
                            if(SCre[ProID].find(lid)==SCre[ProID].end()) {//if not found
                                Tree[rank[lid]].vert[k].second.second -= 1;
                                if (Tree[rank[lid]].vert[k].second.second < 1) {
                                    SCre[lid].insert(Cid);
                                    OC.insert(OrderCompMin(lid));
                                    OCdis[make_pair(lid, Cid)] = Cw + Lnei[j].second;
                                }
                            }
                        }
                        break;
                    }
                }
            }

            //get the new value of shortcut
            //	cout<<Cw<<" increase to ";
            int newCw=INF; int countwt=0;

            for(auto it2=Neighbor[ProID].begin();it2!=Neighbor[ProID].end();++it2){
                if(it2->first==Cid){
                    newCw=it2->second;//the weight value in the original graph
                    countwt=1;
                    break;
                }
            }

            int ssw,wtt,wid;
            vector<pair<int,int>> Wnodes;
            Wnodes.clear();

            if(ProID<Cid)
                Wnodes=SCconNodesMT[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
            else
                Wnodes=SCconNodesMT[Cid][ProID];
            if(!Wnodes.empty()){
                for(int i=0;i<Wnodes.size();i++){
                    wid=Wnodes[i].first;
                    for(int j=0;j<Tree[rank[wid]].vert.size();j++){
                        if(Tree[rank[wid]].vert[j].first==ProID){
                            ssw=Tree[rank[wid]].vert[j].second.first;
                        }
                        if(Tree[rank[wid]].vert[j].first==Cid){
                            wtt=Tree[rank[wid]].vert[j].second.first;
                        }
                    }

                    if(ssw+wtt<newCw){
                        newCw=ssw+wtt;
                        countwt=1;
                    }else if(ssw+wtt==newCw){
                        countwt+=1;
                    }
                }
            }

            //cout<<Cw<<endl;
            //refresh the shortcut to the new value
            for(int i=0;i<Tree[rank[ProID]].vert.size();i++){
                if(Tree[rank[ProID]].vert[i].first==Cid){
                    Tree[rank[ProID]].vert[i].second.first=newCw;
                    Tree[rank[ProID]].vert[i].second.second=countwt;
                    break;
                }
            }


            if(algoChoice==H2H){
                if(newCw>Cw){
                    //before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
                    if(Tree[rank[ProID]].FN[cidH]){
                        influence=true;
                        //higher than Cid
                        for(int i=0;i<cidH;i++){
                            if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[Cid]].dis[i]){
                                Tree[rank[ProID]].cnt[i]-=1;
                            }
                        }

                        //equal to Cid
                        Tree[rank[ProID]].FN[cidH]=false;
                        Tree[rank[ProID]].cnt[cidH]-=1;
                        if(flag){
                            cout<<"Flag. "<<ProID<<" "<<Cid<<" "<<Cw<<" "<<newCw<<" "<<Tree[rank[ProID]].cnt[cidH]<<endl;
                            flag=false;
                        }

                        //lower than Cid
                        for(int i=cidH+1;i<Tree[rank[ProID]].dis.size();i++){
                            if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[line[i]]].dis[cidH]){
                                Tree[rank[ProID]].cnt[i]-=1;
                            }
                        }
                    }
                }

            }

        }

        if(influence){
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetOverlay.size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[ProID],r;
            for(int i=0;i<ProBeginVertexSetOverlay.size();i++){
                r=rank[ProBeginVertexSetOverlay[i]];
                if(LCAQueryOverlay(rnew,r)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetOverlay[i]);
                }
            }
            ProBeginVertexSetOverlay=ProBeginVertexSetNew;
        }

    }

    if(ifLabelU){
        int ProBeginVertexID;
        for(int i=0;i<ProBeginVertexSetOverlay.size();i++){
            ProBeginVertexID=ProBeginVertexSetOverlay[i];
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
            while(Tree[rank[pachidd]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);

            eachNodeProcessIncrease1(rank[ProBeginVertexID], linee,checknum,Tree,rank,VidtoTNid);
        }
    }

    //return checknum;
}


void Graph::IncreaseOverlayBatchLabel(vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int> &ProBeginVertexSet, vector<vector<int>> &VidtoTNid) {
    int ProBeginVertexID;
    int checknum=0;
    for(int i=0;i<ProBeginVertexSet.size();i++){
        ProBeginVertexID=ProBeginVertexSet[i];
        vector<int> linee; //linee.clear();
        linee.reserve(heightMax);
        int pachidd=Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
        while(Tree[rank[pachidd]].height>1){
            linee.insert(linee.begin(),pachidd);
            pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
        }
        linee.insert(linee.begin(),pachidd);

        eachNodeProcessIncrease1(rank[ProBeginVertexID], linee,checknum,Tree,rank,VidtoTNid);
    }
}

void Graph::IncreasePartiBatch(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, vector<pair<pair<int,int>,int>>& updatedSC, bool ifLabelU){
    int checknum=0;
    map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
    OCdis.clear();

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<OrderCompp>> SCre; //SCre.clear();
    SCre.assign(node_num,set<OrderCompp>());//{vertexID, set<int>}
    set<OrderCompMin> OC; OC.clear();//vertexID in decreasing node order

    for(int k=0;k<wBatch.size();k++){
        int a=wBatch[k].first.first;
        int b=wBatch[k].first.second;
        int oldW=wBatch[k].second.first;
        int newW=wBatch[k].second.second;

        if(oldW<newW){
            for(int i=0;i<Neighbors[a].size();i++){
                if(Neighbors[a][i].first==b){
                    if(Neighbors[a][i].second<newW){
                        Neighbors[a][i].second=newW;
                    }else{
                        cout<<"Wrong update. "<<Neighbors[a][i].second<<" "<<newW<<endl; exit(1);
                    }
                    break;
                }
            }
            for(int i=0;i<Neighbors[b].size();i++){
                if(Neighbors[b][i].first==a){
                    if(Neighbors[b][i].second<newW){
                        Neighbors[b][i].second=newW;
                    }else{
                        cout<<"Wrong update. "<<Neighbors[b][i].second<<" "<<newW<<endl; exit(1);
                    }
                    break;
                }
            }

            int lid,hid;
            if(NodeOrder[a]<NodeOrder[b]){
                lid=a;hid=b;
            }else{
                lid=b;hid=a;
            }
            int lidM=IDMap[lid];

            for(int i=0;i<Tree[rank[lidM]].vert.size();i++){
                if(Tree[rank[lidM]].vert[i].first==hid){
                    if(Tree[rank[lidM]].vert[i].second.first==oldW){
                        Tree[rank[lidM]].vert[i].second.second-=1;
                        if(Tree[rank[lidM]].vert[i].second.second<1){
                            OCdis[make_pair(lid,hid)]=oldW;
                            SCre[lid].insert(hid);
                            OC.insert(OrderCompMin(lid));
                        }
                    }
                    break;
                }
            }
        }
    }

//    vector<int> ProBeginVertexSet; ProBeginVertexSet.clear();
    ProBeginVertexSetParti[pid].clear();
    vector<int> ProBeginVertexSetNew;
    bool influence;
    int ProID, ProIDM; vector<int> line;
    while(!OC.empty()){
        ProID=(*OC.begin()).x; ProIDM=IDMap[ProID];
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProIDM]].vert;
        influence=false;

        //each ProID corresponds to a line
        line.clear(); line.reserve(heightMax);
        int pachid=ProID;
        while(Tree[rank[IDMap[pachid]]].height>1){
            line.insert(line.begin(),pachid);
            pachid=Tree[Tree[rank[IDMap[pachid]]].pa].uniqueVertex;
        }
        line.insert(line.begin(),pachid);

        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=it->x; int Cw=OCdis[make_pair(ProID,Cid)];
            int cidH=Tree[rank[IDMap[Cid]]].height-1;
            int CidM=IDMap[Cid];

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }
            }


            //check the affected shortcuts
            int hid2,lid2,lid2M;
            for(int j=0;j<Tree[rank[CidM]].vert.size();j++){
                hid2=Tree[rank[CidM]].vert[j].first;
                if(Hnei.find(hid2)!=Hnei.end()){
                    if(Cw+Hnei[hid2]==Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.second-=1;
                        if(Tree[rank[CidM]].vert[j].second.second<1){
                            SCre[Cid].insert(hid2);
                            OC.insert(OrderCompMin(Cid));
                            OCdis[make_pair(Cid,hid2)]=Cw+Hnei[hid2];
                        }
                    }
                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid2=Lnei[j].first; lid2M=IDMap[lid2];
                for(int k=0;k<Tree[rank[lid2M]].vert.size();k++){
                    if(Tree[rank[lid2M]].vert[k].first==Cid){
                        if(Tree[rank[lid2M]].vert[k].second.first==Cw+Lnei[j].second){
                            if(SCre[ProID].find(lid2)==SCre[ProID].end()) {//if not found
                                Tree[rank[lid2M]].vert[k].second.second -= 1;
                                if (Tree[rank[lid2M]].vert[k].second.second < 1) {
                                    SCre[lid2].insert(Cid);
                                    OC.insert(OrderCompMin(lid2));
                                    OCdis[make_pair(lid2, Cid)] = Cw + Lnei[j].second;
                                }
                            }
                        }
                        break;
                    }
                }
            }

            //get the new value of shortcut
            //	cout<<Cw<<" increase to ";
            int newCw=INF; int countwt=0;

            for(auto it2=Neighbors[ProID].begin();it2!=Neighbors[ProID].end();++it2){
                if(it2->first==Cid){
                    newCw=it2->second;//the weight value in the original graph
                    countwt=1;
                    break;
                }
            }

            int ssw,wtt,wid,widM;
            vector<pair<int,int>> Wnodes;
            Wnodes.clear();

            if(ProID<Cid)
                Wnodes=SCconNodesMT[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
            else
                Wnodes=SCconNodesMT[Cid][ProID];
            if(Wnodes.size()>0){
                for(int i=0;i<Wnodes.size();i++){
                    wid=Wnodes[i].first; widM=IDMap[wid];
                    for(int j=0;j<Tree[rank[widM]].vert.size();j++){
                        if(Tree[rank[widM]].vert[j].first==ProID){
                            ssw=Tree[rank[widM]].vert[j].second.first;
                        }
                        if(Tree[rank[widM]].vert[j].first==Cid){
                            wtt=Tree[rank[widM]].vert[j].second.first;
                        }
                    }

                    if(ssw+wtt<newCw){
                        newCw=ssw+wtt;
                        countwt=1;
                    }else if(ssw+wtt==newCw){
                        countwt+=1;
                    }
                }
            }

//            if(PartiTag[ProID].second){//if boundary vertex
//                updatedSC.emplace_back(make_pair(ProID,Cid),newCw);
//            }

            //cout<<Cw<<endl;
            //refresh the shortcut to the new value
            for(int i=0;i<Tree[rank[ProIDM]].vert.size();i++){
                if(Tree[rank[ProIDM]].vert[i].first==Cid){
                    Tree[rank[ProIDM]].vert[i].second.first=newCw;
                    Tree[rank[ProIDM]].vert[i].second.second=countwt;
                    break;
                }
            }

            if(algoChoice==H2H){
                if(newCw>Cw){
                    //before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
                    if(Tree[rank[ProIDM]].FN[cidH]){
                        influence=true;
                        //higher than Cid
                        for(int i=0;i<cidH;i++){
                            if(Tree[rank[ProIDM]].dis[i]==Cw+Tree[rank[CidM]].dis[i]){
                                Tree[rank[ProIDM]].cnt[i]-=1;
                            }
                        }

                        //equal to Cid
                        Tree[rank[ProIDM]].FN[cidH]=false;
                        Tree[rank[ProIDM]].cnt[cidH]-=1;

                        //lower than Cid
                        for(int i=cidH+1;i<Tree[rank[ProIDM]].dis.size();i++){
                            if(Tree[rank[ProIDM]].dis[i]==Cw+Tree[rank[IDMap[line[i]]]].dis[cidH]){
                                Tree[rank[ProIDM]].cnt[i]-=1;
                            }
                        }
                    }
                }

            }



        }

        if(influence){
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetParti[pid].size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[IDMap[ProID]],r;
            for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
                r=rank[IDMap[ProBeginVertexSetParti[pid][i]]];
                if(LCAQueryPartition(rnew,r,pid)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetParti[pid][i]);
                }
            }
            ProBeginVertexSetParti[pid]=ProBeginVertexSetNew;
        }

    }

    if(ifLabelU){
        int ProBeginVertexID;
        for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
            ProBeginVertexID=ProBeginVertexSetParti[pid][i];
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[IDMap[ProBeginVertexID]]].pa].uniqueVertex;
            while(Tree[rank[IDMap[pachidd]]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[IDMap[pachidd]]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);

            eachNodeProcessIncrease1Parti(rank[IDMap[ProBeginVertexID]], linee,checknum,Tree,rank,VidtoTNid);
        }
    }

    //return checknum;
}

void Graph::IncreasePartiBatchForOpt(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, vector<pair<pair<int,int>,int>>& updatedSC, bool ifLabelU){
    int checknum=0;
    map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
    OCdis.clear();

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<OrderCompp>> SCre; //SCre.clear();
    SCre.assign(node_num,set<OrderCompp>());//{vertexID, set<int>}
    set<OrderCompMin> OC; OC.clear();//vertexID in decreasing node order

    for(int k=0;k<wBatch.size();k++){
        int a=wBatch[k].first.first;
        int b=wBatch[k].first.second;
        int oldW=wBatch[k].second.first;
        int newW=wBatch[k].second.second;

        if(oldW!=newW){//oldW may be incorrect!
            for(int i=0;i<Neighbors[a].size();i++){
                if(Neighbors[a][i].first==b){
                    if(oldW!=Neighbors[a][i].second){
                        cout<<"Edge weight inconsistent! "<<a<<" "<<b<<" "<<Neighbors[a][i].second<<" "<<oldW<<endl;
                        oldW=Neighbors[a][i].second;
                    }
                    Neighbors[a][i].second=newW;
                    break;
                }
            }
            for(int i=0;i<Neighbors[b].size();i++){
                if(Neighbors[b][i].first==a){
                    if(oldW!=Neighbors[b][i].second){
                        cout<<"Edge weight inconsistent! "<<b<<" "<<a<<" "<<Neighbors[b][i].second<<" "<<oldW<<endl;
                        oldW=Neighbors[b][i].second;
                    }
                    Neighbors[b][i].second=newW;
                    break;
                }
            }

            int lid,hid;
            if(NodeOrder[a]<NodeOrder[b]){
                lid=a;hid=b;
            }else{
                lid=b;hid=a;
            }
            int lidM=IDMap[lid];

            for(int i=0;i<Tree[rank[lidM]].vert.size();i++){
                if(Tree[rank[lidM]].vert[i].first==hid){
                    if(Tree[rank[lidM]].vert[i].second.first==oldW){
                        Tree[rank[lidM]].vert[i].second.second-=1;
                        if(Tree[rank[lidM]].vert[i].second.second<1){
                            OCdis[make_pair(lid,hid)]=oldW;
                            SCre[lid].insert(hid);
                            OC.insert(OrderCompMin(lid));
                        }
                    }
                    break;
                }
            }
        }
    }

//    vector<int> ProBeginVertexSet; ProBeginVertexSet.clear();
    if(!ifLabelU){
        ProBeginVertexSetParti[pid].clear();
    }

    vector<int> ProBeginVertexSetNew;
    bool influence;
    int ProID, ProIDM; vector<int> line;
    while(!OC.empty()){
        ProID=(*OC.begin()).x; ProIDM=IDMap[ProID];
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProIDM]].vert;
        influence=false;

        //each ProID corresponds to a line
        line.clear(); line.reserve(heightMax);
        int pachid=ProID;
        while(Tree[rank[IDMap[pachid]]].height>1){
            line.insert(line.begin(),pachid);
            pachid=Tree[Tree[rank[IDMap[pachid]]].pa].uniqueVertex;
        }
        line.insert(line.begin(),pachid);

        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=it->x; int Cw=OCdis[make_pair(ProID,Cid)];
            int cidH=Tree[rank[IDMap[Cid]]].height-1;
            int CidM=IDMap[Cid];

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }
            }
            //check the affected shortcuts
            int hid2,lid2,lid2M;
            for(int j=0;j<Tree[rank[CidM]].vert.size();j++){
                hid2=Tree[rank[CidM]].vert[j].first;
                if(Hnei.find(hid2)!=Hnei.end()){
                    if(Cw+Hnei[hid2]==Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.second-=1;
                        if(Tree[rank[CidM]].vert[j].second.second<1){
                            SCre[Cid].insert(hid2);
                            OC.insert(OrderCompMin(Cid));
                            OCdis[make_pair(Cid,hid2)]=Cw+Hnei[hid2];
                        }
                    }
                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid2=Lnei[j].first; lid2M=IDMap[lid2];
                for(int k=0;k<Tree[rank[lid2M]].vert.size();k++){
                    if(Tree[rank[lid2M]].vert[k].first==Cid){
                        if(Tree[rank[lid2M]].vert[k].second.first==Cw+Lnei[j].second){
                            if(SCre[ProID].find(lid2)==SCre[ProID].end()) {//if not found
                                Tree[rank[lid2M]].vert[k].second.second -= 1;
                                if (Tree[rank[lid2M]].vert[k].second.second < 1) {
                                    SCre[lid2].insert(Cid);
                                    OC.insert(OrderCompMin(lid2));
                                    OCdis[make_pair(lid2, Cid)] = Cw + Lnei[j].second;
                                }
                            }
                        }
                        break;
                    }
                }
            }


            //get the new value of shortcut
            //	cout<<Cw<<" increase to ";
            int newCw=INF; int countwt=0;

            for(auto it2=Neighbors[ProID].begin();it2!=Neighbors[ProID].end();++it2){
                if(it2->first==Cid){
                    newCw=it2->second;//the weight value in the original graph
                    countwt=1;
                    break;
                }
            }

            int ssw,wtt,wid,widM;
            vector<pair<int,int>> Wnodes;
            Wnodes.clear();

            if(ProID<Cid)
                Wnodes=SCconNodesMT[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
            else
                Wnodes=SCconNodesMT[Cid][ProID];
            if(Wnodes.size()>0){
                for(int i=0;i<Wnodes.size();i++){
                    wid=Wnodes[i].first; widM=IDMap[wid];
                    for(int j=0;j<Tree[rank[widM]].vert.size();j++){
                        if(Tree[rank[widM]].vert[j].first==ProID){
                            ssw=Tree[rank[widM]].vert[j].second.first;
                        }
                        if(Tree[rank[widM]].vert[j].first==Cid){
                            wtt=Tree[rank[widM]].vert[j].second.first;
                        }
                    }

                    if(ssw+wtt<newCw){
                        newCw=ssw+wtt;
                        countwt=1;
                    }else if(ssw+wtt==newCw){
                        countwt+=1;
                    }
                }
            }

            if(PartiTag[ProID].second){//if boundary vertex
//                cout<<"Boundary shortcut update: "<<ProID<<" "<<Cid<<" "<<Cw<<" "<<newCw<<endl;
                updatedSC.emplace_back(make_pair(ProID,Cid),newCw);
            }

//            if(ProID==788 && Cid==835){
//                cout<<"Find. "<<ProID<<" "<<Cid<<" "<<newCw<<" "<<countwt<<endl;
//            }
            //cout<<Cw<<endl;
            //refresh the shortcut to the new value
            for(int i=0;i<Tree[rank[ProIDM]].vert.size();i++){
                if(Tree[rank[ProIDM]].vert[i].first==Cid){
                    Tree[rank[ProIDM]].vert[i].second.first=newCw;
                    Tree[rank[ProIDM]].vert[i].second.second=countwt;
                    break;
                }
            }

            if(algoChoice==H2H){
                if(newCw>Cw){
                    if(Tree[rank[ProIDM]].FN[cidH]){
                        influence=true;
                        //higher than Cid
                        for(int i=0;i<cidH;i++){
                            if(Tree[rank[ProIDM]].dis[i]==Cw+Tree[rank[CidM]].dis[i]){
                                Tree[rank[ProIDM]].cnt[i]-=1;
                            }
                        }

                        //equal to Cid
                        Tree[rank[ProIDM]].FN[cidH]=false;
                        Tree[rank[ProIDM]].cnt[cidH]-=1;

                        //lower than Cid
                        for(int i=cidH+1;i<Tree[rank[ProIDM]].dis.size();i++){
                            if(Tree[rank[ProIDM]].dis[i]==Cw+Tree[rank[IDMap[line[i]]]].dis[cidH]){
                                Tree[rank[ProIDM]].cnt[i]-=1;
                            }
                        }
                    }
                }
            }


        }

        if(influence){
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetParti[pid].size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[IDMap[ProID]],r;
            for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
                r=rank[IDMap[ProBeginVertexSetParti[pid][i]]];
                if(LCAQueryPartition(rnew,r,pid)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetParti[pid][i]);
                }
            }
            ProBeginVertexSetParti[pid]=ProBeginVertexSetNew;
        }

    }

    if(ifLabelU){
        int ProBeginVertexID;
        for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
            ProBeginVertexID=ProBeginVertexSetParti[pid][i];
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[IDMap[ProBeginVertexID]]].pa].uniqueVertex;
            while(Tree[rank[IDMap[pachidd]]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[IDMap[pachidd]]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);

            eachNodeProcessIncrease1Parti(rank[IDMap[ProBeginVertexID]], linee,checknum,Tree,rank,VidtoTNid);
        }
    }

    //return checknum;
}


void Graph::IncreasePartiBatchPost(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, bool ifLabelU){
    int checknum=0;
    map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
    OCdis.clear();

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<OrderCompp>> SCre; //SCre.clear();
    SCre.assign(node_num,set<OrderCompp>());//{vertexID, set<int>}
    set<OrderCompMin> OC; OC.clear();//vertexID in decreasing node order

    for(int k=0;k<wBatch.size();k++){
        int a=wBatch[k].first.first;
        int b=wBatch[k].first.second;
        int oldW=wBatch[k].second.first;
        int newW=wBatch[k].second.second;

        if(oldW!=newW){
            if(Neighbors[a].find(b)!=Neighbors[a].end()){
                Neighbors[a][b]=newW;
            }else{
                cout<<"Not found edge! "<<endl; exit(1);
            }

            if(Neighbors[b].find(a)!=Neighbors[b].end()){
                Neighbors[b][a]=newW;
            }else{
                cout<<"Not found edge! "<<endl; exit(1);
            }

            int lid,hid;
            if(NodeOrder[a]<NodeOrder[b]){
                lid=a;hid=b;
            }else{
                lid=b;hid=a;
            }
            int lidM=IDMap[lid];

            for(int i=0;i<Tree[rank[lidM]].vert.size();i++){
                if(Tree[rank[lidM]].vert[i].first==hid){
                    if(Tree[rank[lidM]].vert[i].second.first==oldW){
                        Tree[rank[lidM]].vert[i].second.second-=1;
                        if(Tree[rank[lidM]].vert[i].second.second<1){
                            OCdis[make_pair(lid,hid)]=oldW;
                            SCre[lid].insert(hid);
                            OC.insert(OrderCompMin(lid));
                        }
                    }
                    break;
                }
            }
        }
    }

//    vector<int> ProBeginVertexSet; ProBeginVertexSet.clear();
    ProBeginVertexSetParti[pid].clear();
    vector<int> ProBeginVertexSetNew;
    bool influence;
    int ProID, ProIDM; vector<int> line;
    while(!OC.empty()){
        ProID=(*OC.begin()).x; ProIDM=IDMap[ProID];
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProIDM]].vert;
        influence=false;

        //each ProID corresponds to a line
        line.clear(); line.reserve(heightMax);
        int pachid=ProID;
        while(Tree[rank[IDMap[pachid]]].height>1){
            line.insert(line.begin(),pachid);
            pachid=Tree[Tree[rank[IDMap[pachid]]].pa].uniqueVertex;
        }
        line.insert(line.begin(),pachid);

        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=it->x; int Cw=OCdis[make_pair(ProID,Cid)];
            int cidH=Tree[rank[IDMap[Cid]]].height-1;
            int CidM=IDMap[Cid];

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }
            }
            //check the affected shortcuts
            int hid2,lid2,lid2M;
            for(int j=0;j<Tree[rank[CidM]].vert.size();j++){
                hid2=Tree[rank[CidM]].vert[j].first;
                if(Hnei.find(hid2)!=Hnei.end()){
                    if(Cw+Hnei[hid2]==Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.second-=1;
                        if(Tree[rank[CidM]].vert[j].second.second<1){
                            SCre[Cid].insert(hid2);
                            OC.insert(OrderCompMin(Cid));
                            OCdis[make_pair(Cid,hid2)]=Cw+Hnei[hid2];
                        }
                    }
                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid2=Lnei[j].first; lid2M=IDMap[lid2];
                for(int k=0;k<Tree[rank[lid2M]].vert.size();k++){
                    if(Tree[rank[lid2M]].vert[k].first==Cid){
                        if(Tree[rank[lid2M]].vert[k].second.first==Cw+Lnei[j].second){
                            if(SCre[ProID].find(lid2)==SCre[ProID].end()) {//if not found
                                Tree[rank[lid2M]].vert[k].second.second -= 1;
                                if (Tree[rank[lid2M]].vert[k].second.second < 1) {
                                    SCre[lid2].insert(Cid);
                                    OC.insert(OrderCompMin(lid2));
                                    OCdis[make_pair(lid2, Cid)] = Cw + Lnei[j].second;
                                }
                            }
                        }
                        break;
                    }
                }
            }

            //get the new value of shortcut
            //	cout<<Cw<<" increase to ";
            int newCw=INF; int countwt=0;

            for(auto it2=Neighbors[ProID].begin();it2!=Neighbors[ProID].end();++it2){
                if(it2->first==Cid){
                    newCw=it2->second;//the weight value in the original graph
                    countwt=1;
                    break;
                }
            }

            int ssw,wtt,wid,widM;
            vector<pair<int,int>> Wnodes;
            Wnodes.clear();
            /*if(SCconNodes.find(make_pair(ProID,Cid))!=SCconNodes.end())
                Wnodes=SCconNodes[make_pair(ProID,Cid)];
            else if(SCconNodes.find(make_pair(Cid,ProID))!=SCconNodes.end())
                Wnodes=SCconNodes[make_pair(Cid,ProID)];*/
            if(ProID<Cid)
                Wnodes=SCconNodesMT[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
            else
                Wnodes=SCconNodesMT[Cid][ProID];
            if(Wnodes.size()>0){
                for(int i=0;i<Wnodes.size();i++){
                    wid=Wnodes[i].first; widM=IDMap[wid];
                    for(int j=0;j<Tree[rank[widM]].vert.size();j++){
                        if(Tree[rank[widM]].vert[j].first==ProID){
                            ssw=Tree[rank[widM]].vert[j].second.first;
                        }
                        if(Tree[rank[widM]].vert[j].first==Cid){
                            wtt=Tree[rank[widM]].vert[j].second.first;
                        }
                    }

                    if(ssw+wtt<newCw){
                        newCw=ssw+wtt;
                        countwt=1;
                    }else if(ssw+wtt==newCw){
                        countwt+=1;
                    }
                }
            }

            //cout<<Cw<<endl;
            //refresh the shortcut to the new value
            for(int i=0;i<Tree[rank[ProIDM]].vert.size();i++){
                if(Tree[rank[ProIDM]].vert[i].first==Cid){
                    Tree[rank[ProIDM]].vert[i].second.first=newCw;
                    Tree[rank[ProIDM]].vert[i].second.second=countwt;
                    break;
                }
            }

            if(algoChoice==H2H){
                if(newCw>Cw){
                    //before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
                    if(Tree[rank[ProIDM]].FN[cidH]){
                        influence=true;
                        //higher than Cid
                        for(int i=0;i<cidH;i++){
                            if(Tree[rank[ProIDM]].dis[i]==Cw+Tree[rank[CidM]].dis[i]){
                                Tree[rank[ProIDM]].cnt[i]-=1;
                            }
                        }

                        //equal to Cid
                        Tree[rank[ProIDM]].FN[cidH]=false;
                        Tree[rank[ProIDM]].cnt[cidH]-=1;

                        //lower than Cid
                        for(int i=cidH+1;i<Tree[rank[ProIDM]].dis.size();i++){
                            if(Tree[rank[ProIDM]].dis[i]==Cw+Tree[rank[IDMap[line[i]]]].dis[cidH]){
                                Tree[rank[ProIDM]].cnt[i]-=1;
                            }
                        }
                    }
                }

            }



        }

        if(influence){
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetParti[pid].size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[IDMap[ProID]],r;
            for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
                r=rank[IDMap[ProBeginVertexSetParti[pid][i]]];
                if(LCAQueryPartitionPost(rnew,r,pid)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetParti[pid][i]);
                }
            }
            ProBeginVertexSetParti[pid]=ProBeginVertexSetNew;
        }

    }

    if(ifLabelU){
        int ProBeginVertexID;
        for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
            ProBeginVertexID=ProBeginVertexSetParti[pid][i];
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[IDMap[ProBeginVertexID]]].pa].uniqueVertex;
            while(Tree[rank[IDMap[pachidd]]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[IDMap[pachidd]]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);

            eachNodeProcessIncrease1Parti(rank[IDMap[ProBeginVertexID]], linee,checknum,Tree,rank,VidtoTNid);
        }
    }

    //return checknum;
}

void Graph::IncreasePartiBatchPre(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU){
    int checknum=0;
    map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
    OCdis.clear();

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<OrderCompp>> SCre; //SCre.clear();
    SCre.assign(node_num,set<OrderCompp>());//{vertexID, set<int>}
    set<OrderCompMin> OC; OC.clear();//vertexID in decreasing node order
    bool flag=false;

    for(int k=0;k<wBatch.size();k++){
        int a=wBatch[k].first.first;
        int b=wBatch[k].first.second;
        int oldW=wBatch[k].second.first;
        int newW=wBatch[k].second.second;

        if(oldW!=newW){
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

            int lid,hid;
            if(NodeOrder[a]<NodeOrder[b]){
                lid=a;hid=b;
            }else{
                lid=b;hid=a;
            }
            int lidM=IDMap[lid];

            for(int i=0;i<Tree[rank[lidM]].vert.size();i++){
                if(Tree[rank[lidM]].vert[i].first==hid){
//                    if(lid==7988 && hid==7984){
//                        cout<<"Find 1. "<<lid<<" "<<hid<<" "<<Tree[rank[lidM]].vert[i].second.first<<" "<<Tree[rank[lidM]].vert[i].second.second <<" "<<oldW<<endl;
//                    }
                    if(Tree[rank[lidM]].vert[i].second.first==oldW){
                        Tree[rank[lidM]].vert[i].second.second-=1;
                        if(Tree[rank[lidM]].vert[i].second.second<1){
                            OCdis[make_pair(lid,hid)]=oldW;
                            SCre[lid].insert(hid);
                            OC.insert(OrderCompMin(lid));
                        }
                    }
                    break;
                }
            }
        }
    }

//    vector<int> ProBeginVertexSet; ProBeginVertexSet.clear();
    ProBeginVertexSetParti[pid].clear();
    vector<int> ProBeginVertexSetNew;
    bool influence;
    int ProID, ProIDM; vector<int> line;
    while(!OC.empty()){
        ProID=(*OC.begin()).x; ProIDM=IDMap[ProID];
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProIDM]].vert;
        influence=false;

        //each ProID corresponds to a line
        line.clear(); line.reserve(heightMax);
        int pachid=ProID;
        while(Tree[rank[IDMap[pachid]]].height>1){
            line.insert(line.begin(),pachid);
            pachid=Tree[Tree[rank[IDMap[pachid]]].pa].uniqueVertex;
        }
        line.insert(line.begin(),pachid);

        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=it->x; int Cw=OCdis[make_pair(ProID,Cid)];
            int cidH=Tree[rank[IDMap[Cid]]].height-1;
            int CidM=IDMap[Cid];

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }
            }



            //check the affected shortcuts
            int hid2,lid2,lid2M;
            for(int j=0;j<Tree[rank[CidM]].vert.size();j++){
                hid2=Tree[rank[CidM]].vert[j].first;
                if(Hnei.find(hid2)!=Hnei.end()){
                    if(Cw+Hnei[hid2]==Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.second-=1;
                        if(Tree[rank[CidM]].vert[j].second.second<1){
                            SCre[Cid].insert(hid2);
                            OC.insert(OrderCompMin(Cid));
                            OCdis[make_pair(Cid,hid2)]=Cw+Hnei[hid2];
                        }
                    }
                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid2=Lnei[j].first; lid2M=IDMap[lid2];
                for(int k=0;k<Tree[rank[lid2M]].vert.size();k++){
                    if(Tree[rank[lid2M]].vert[k].first==Cid){
                        if(Tree[rank[lid2M]].vert[k].second.first==Cw+Lnei[j].second){
                            if(SCre[ProID].find(lid2)==SCre[ProID].end()) {//if not found
                                Tree[rank[lid2M]].vert[k].second.second -= 1;
                                if (Tree[rank[lid2M]].vert[k].second.second < 1) {
                                    SCre[lid2].insert(Cid);
                                    OC.insert(OrderCompMin(lid2));
                                    OCdis[make_pair(lid2, Cid)] = Cw + Lnei[j].second;
                                }
                            }
                        }
                        break;
                    }
                }
            }

            //get the new value of shortcut
            //	cout<<Cw<<" increase to ";
            int newCw=INF; int countwt=0;
            for(auto it2=Neighbors[ProID].begin();it2!=Neighbors[ProID].end();++it2){
                if(it2->first==Cid){
                    newCw=it2->second;//the weight value in the original graph
                    countwt=1;
                    break;
                }
            }

            int ssw,wtt,wid,widM;
            vector<pair<int,int>> Wnodes;
            Wnodes.clear();
            /*if(SCconNodes.find(make_pair(ProID,Cid))!=SCconNodes.end())
                Wnodes=SCconNodes[make_pair(ProID,Cid)];
            else if(SCconNodes.find(make_pair(Cid,ProID))!=SCconNodes.end())
                Wnodes=SCconNodes[make_pair(Cid,ProID)];*/
            if(ProID<Cid)
                Wnodes=SCconNodesMTP[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
            else
                Wnodes=SCconNodesMTP[Cid][ProID];
            if(Wnodes.size()>0){
                for(int i=0;i<Wnodes.size();i++){
                    wid=Wnodes[i].first; widM=IDMap[wid];
                    for(int j=0;j<Tree[rank[widM]].vert.size();j++){
                        if(Tree[rank[widM]].vert[j].first==ProID){
                            ssw=Tree[rank[widM]].vert[j].second.first;
                        }
                        if(Tree[rank[widM]].vert[j].first==Cid){
                            wtt=Tree[rank[widM]].vert[j].second.first;
                        }
                    }

                    if(ssw+wtt<newCw){
                        newCw=ssw+wtt;
                        countwt=1;
                    }else if(ssw+wtt==newCw){
                        countwt+=1;
                    }
                }
            }

//            if(ProID==7988 && Cid==7984){
//                cout<<"Find 2. "<<ProID<<" "<<Cid<<" "<<Cw<< " "<<newCw<<endl;
//            }
            //cout<<Cw<<endl;
            //refresh the shortcut to the new value
            for(int i=0;i<Tree[rank[ProIDM]].vert.size();i++){
                if(Tree[rank[ProIDM]].vert[i].first==Cid){
                    Tree[rank[ProIDM]].vert[i].second.first=newCw;
                    Tree[rank[ProIDM]].vert[i].second.second=countwt;
                    break;
                }
            }

            if(algoChoice==H2H){
                if(newCw>Cw){
                    //before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
                    if(Tree[rank[ProIDM]].FN[cidH]){
                        influence=true;
                        //higher than Cid
                        for(int i=0;i<cidH;i++){
                            if(Tree[rank[ProIDM]].dis[i]==Cw+Tree[rank[CidM]].dis[i]){
                                Tree[rank[ProIDM]].cnt[i]-=1;
                            }
                        }

                        //equal to Cid
                        Tree[rank[ProIDM]].FN[cidH]=false;
                        Tree[rank[ProIDM]].cnt[cidH]-=1;
//                        if(ProID==7988 && Cid==7984){
//                            cout<<"Find 3. "<<ProID<<" "<<Cid<<" "<<Tree[rank[ProIDM]].cnt[cidH]<<" "<<Tree[rank[ProIDM]].height<<endl;
//                            flag=true;
//                        }
                        //lower than Cid
                        for(int i=cidH+1;i<Tree[rank[ProIDM]].dis.size();i++){
                            if(Tree[rank[ProIDM]].dis[i]==Cw+Tree[rank[IDMap[line[i]]]].dis[cidH]){
                                Tree[rank[ProIDM]].cnt[i]-=1;
                            }
                        }
                    }
                }

            }



        }

        if(influence){
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetParti[pid].size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[IDMap[ProID]],r;
            for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
                r=rank[IDMap[ProBeginVertexSetParti[pid][i]]];
                if(LCAQueryPartition(rnew,r,pid)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetParti[pid][i]);
                }
            }
            ProBeginVertexSetParti[pid]=ProBeginVertexSetNew;
//            if(flag){
//                cout<<"Flag "<<ProID;
//                for(int j=0;j<ProBeginVertexSetNew.size();++j){
//                    cout<<" "<<ProBeginVertexSetNew[j]<<"("<<Tree[rank[IDMap[ProBeginVertexSetNew[j]]]].height<<")";
//                }
//                cout<<endl;
////                flag=false;
//            }
        }

    }

    if(ifLabelU){
        int ProBeginVertexID;
        for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
            ProBeginVertexID=ProBeginVertexSetParti[pid][i];
//            cout<<"ProBeginVertex ID "<<ProBeginVertexID<<" "<<Tree[rank[IDMap[ProBeginVertexID]]].uniqueVertex<<" "<<Tree[rank[IDMap[ProBeginVertexID]]].height<<endl;
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[IDMap[ProBeginVertexID]]].pa].uniqueVertex;
            while(Tree[rank[IDMap[pachidd]]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[IDMap[pachidd]]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);

            eachNodeProcessIncrease1Parti(rank[IDMap[ProBeginVertexID]], linee,checknum,Tree,rank,VidtoTNidP);
        }
    }

    //return checknum;
}

void Graph::EdgeDeletePartiBatchPre(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU){
    int checknum=0;
    map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
    OCdis.clear();

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<OrderCompp>> SCre; //SCre.clear();
    SCre.assign(node_num,set<OrderCompp>());//{vertexID, set<int>}
    set<OrderCompMin> OC; OC.clear();//vertexID in decreasing node order

    for(int k=0;k<wBatch.size();k++){
        int a=wBatch[k].first.first;
        int b=wBatch[k].first.second;
        int oldW=wBatch[k].second.first;
        int newW=wBatch[k].second.second;

        if(oldW!=newW){
            for(int i=0;i<Neighbors[a].size();i++){
                if(Neighbors[a][i].first==b){
                    Neighbors[a][i].second=newW;
//                    Neighbors[a].erase(Neighbors[a].begin()+i);
                    break;
                }
            }
            for(int i=0;i<Neighbors[b].size();i++){
                if(Neighbors[b][i].first==a){
                    Neighbors[b][i].second=newW;
                    break;
                }
            }

            int lid,hid;
            if(NodeOrder[a]<NodeOrder[b]){
                lid=a;hid=b;
            }else{
                lid=b;hid=a;
            }
            int lidM=IDMap[lid];

            for(int i=0;i<Tree[rank[lidM]].vert.size();i++){
                if(Tree[rank[lidM]].vert[i].first==hid){
                    if(Tree[rank[lidM]].vert[i].second.first==oldW){
                        Tree[rank[lidM]].vert[i].second.second-=1;
                        if(Tree[rank[lidM]].vert[i].second.second<1){
                            OCdis[make_pair(lid,hid)]=oldW;
                            SCre[lid].insert(hid);
                            OC.insert(OrderCompMin(lid));
                        }
                    }
                    break;
                }
            }
        }
    }

//    vector<int> ProBeginVertexSet; ProBeginVertexSet.clear();
    ProBeginVertexSetParti[pid].clear();
    vector<int> ProBeginVertexSetNew;
    bool influence;
    int ProID, ProIDM; vector<int> line;
    while(!OC.empty()){
        ProID=(*OC.begin()).x; ProIDM=IDMap[ProID];
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProIDM]].vert;
        influence=false;

        //each ProID corresponds to a line
        line.clear(); line.reserve(heightMax);
        int pachid=ProID;
        while(Tree[rank[IDMap[pachid]]].height>1){
            line.insert(line.begin(),pachid);
            pachid=Tree[Tree[rank[IDMap[pachid]]].pa].uniqueVertex;
        }
        line.insert(line.begin(),pachid);

        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=it->x; int Cw=OCdis[make_pair(ProID,Cid)];
            int cidH=Tree[rank[IDMap[Cid]]].height-1;
            int CidM=IDMap[Cid];

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }
            }
            //check the affected shortcuts
            int hid2,lid2,lid2M;
            for(int j=0;j<Tree[rank[CidM]].vert.size();j++){
                hid2=Tree[rank[CidM]].vert[j].first;
                if(Hnei.find(hid2)!=Hnei.end()){
                    if(Cw+Hnei[hid2]==Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.second-=1;
                        if(Tree[rank[CidM]].vert[j].second.second<1){
                            SCre[Cid].insert(hid2);
                            OC.insert(OrderCompMin(Cid));
                            OCdis[make_pair(Cid,hid2)]=Cw+Hnei[hid2];
                        }
                    }
                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid2=Lnei[j].first; lid2M=IDMap[lid2];
                for(int k=0;k<Tree[rank[lid2M]].vert.size();k++){
                    if(Tree[rank[lid2M]].vert[k].first==Cid){
                        if(Tree[rank[lid2M]].vert[k].second.first==Cw+Lnei[j].second){
                            if(SCre[ProID].find(lid2)==SCre[ProID].end()) {//if not found
                                Tree[rank[lid2M]].vert[k].second.second -= 1;
                                if (Tree[rank[lid2M]].vert[k].second.second < 1) {
                                    SCre[lid2].insert(Cid);
                                    OC.insert(OrderCompMin(lid2));
                                    OCdis[make_pair(lid2, Cid)] = Cw + Lnei[j].second;
                                }
                            }
                        }
                        break;
                    }
                }
            }

            //get the new value of shortcut
            //	cout<<Cw<<" increase to ";
            int newCw=INF; int countwt=0;
            for(auto it2=Neighbors[ProID].begin();it2!=Neighbors[ProID].end();++it2){
                if(it2->first==Cid){
                    newCw=it2->second;//the weight value in the original graph
                    countwt=1;
                    break;
                }
            }

            int ssw,wtt,wid,widM;
            vector<pair<int,int>> Wnodes;
            Wnodes.clear();
            /*if(SCconNodes.find(make_pair(ProID,Cid))!=SCconNodes.end())
                Wnodes=SCconNodes[make_pair(ProID,Cid)];
            else if(SCconNodes.find(make_pair(Cid,ProID))!=SCconNodes.end())
                Wnodes=SCconNodes[make_pair(Cid,ProID)];*/
            if(ProID<Cid)
                Wnodes=SCconNodesMTP[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
            else
                Wnodes=SCconNodesMTP[Cid][ProID];
            if(Wnodes.size()>0){
                for(int i=0;i<Wnodes.size();i++){
                    wid=Wnodes[i].first; widM=IDMap[wid];
                    for(int j=0;j<Tree[rank[widM]].vert.size();j++){
                        if(Tree[rank[widM]].vert[j].first==ProID){
                            ssw=Tree[rank[widM]].vert[j].second.first;
                        }
                        if(Tree[rank[widM]].vert[j].first==Cid){
                            wtt=Tree[rank[widM]].vert[j].second.first;
                        }
                    }

                    if(ssw+wtt<newCw){
                        newCw=ssw+wtt;
                        countwt=1;
                    }else if(ssw+wtt==newCw){
                        countwt+=1;
                    }
                }
            }

            //cout<<Cw<<endl;
            //refresh the shortcut to the new value
            for(int i=0;i<Tree[rank[ProIDM]].vert.size();i++){
                if(Tree[rank[ProIDM]].vert[i].first==Cid){
                    Tree[rank[ProIDM]].vert[i].second.first=newCw;
                    Tree[rank[ProIDM]].vert[i].second.second=countwt;
                    break;
                }
            }

            if(algoChoice==H2H){
                if(newCw>Cw){
                    //before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
                    if(Tree[rank[ProIDM]].FN[cidH]){
                        influence=true;
                        //higher than Cid
                        for(int i=0;i<cidH;i++){
                            if(Tree[rank[ProIDM]].dis[i]==Cw+Tree[rank[CidM]].dis[i]){
                                Tree[rank[ProIDM]].cnt[i]-=1;
                            }
                        }

                        //equal to Cid
                        Tree[rank[ProIDM]].FN[cidH]=false;
                        Tree[rank[ProIDM]].cnt[cidH]-=1;

                        //lower than Cid
                        for(int i=cidH+1;i<Tree[rank[ProIDM]].dis.size();i++){
                            if(Tree[rank[ProIDM]].dis[i]==Cw+Tree[rank[IDMap[line[i]]]].dis[cidH]){
                                Tree[rank[ProIDM]].cnt[i]-=1;
                            }
                        }
                    }
                }

            }



        }

        if(influence){
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetParti[pid].size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[IDMap[ProID]],r;
            for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
                r=rank[IDMap[ProBeginVertexSetParti[pid][i]]];
                if(LCAQueryPartition(rnew,r,pid)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetParti[pid][i]);
                }
            }
            ProBeginVertexSetParti[pid]=ProBeginVertexSetNew;
        }

    }

    if(ifLabelU){
        int ProBeginVertexID;
        for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
            ProBeginVertexID=ProBeginVertexSetParti[pid][i];
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[IDMap[ProBeginVertexID]]].pa].uniqueVertex;
            while(Tree[rank[IDMap[pachidd]]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[IDMap[pachidd]]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);

            eachNodeProcessIncrease1Parti(rank[IDMap[ProBeginVertexID]], linee,checknum,Tree,rank,VidtoTNidP);
        }
    }

    //return checknum;
}

void Graph::IncreasePartiBatchLabel(vector<Node> &Tree, vector<int> &rank, int heightMax,  vector<int> &ProBeginVertexSet, vector<vector<int>> &VidtoTNid) {
    int ProBeginVertexID;
    int checknum=0;
    for(int i=0;i<ProBeginVertexSet.size();i++){
        ProBeginVertexID=ProBeginVertexSet[i];
        vector<int> linee; //linee.clear();
        linee.reserve(heightMax);
        int pachidd=Tree[Tree[rank[IDMap[ProBeginVertexID]]].pa].uniqueVertex;
        while(Tree[rank[IDMap[pachidd]]].height>1){
            linee.insert(linee.begin(),pachidd);
            pachidd=Tree[Tree[rank[IDMap[pachidd]]].pa].uniqueVertex;
        }
        linee.insert(linee.begin(),pachidd);

        eachNodeProcessIncrease1Parti(rank[IDMap[ProBeginVertexID]], linee,checknum,Tree,rank,VidtoTNidP);
    }
}

int Graph::ShortcutDisCheck(int ID1, int ID2){
    int d=INF;
    int ID;
    if(PartiTag[ID1].first!=PartiTag[ID2].first){
        cout<<ID1<<" and "<<ID2<<" are not in the same partition!"<<endl; exit(1);
    }
    int PID=PartiTag[ID1].first;
    int lid,hid;
    if(NodeOrder[ID1]>NodeOrder[ID2]){
        hid=ID1, lid=ID2;
    }else{
        hid=ID2, lid=ID1;
    }
    bool flag=false;
/*    for(int i=0;i<Trees[PID][ranks[PID][IDMap[lid]]].vert.size();i++){
        if(Trees[PID][ranks[PID][IDMap[lid]]].vert[i].first==hid){
            d=Trees[PID][ranks[PID][IDMap[lid]]].vert[i].second.first;
            flag=true;
            break;
        }
    }*/

    int Cw=INF;
    for(auto it2=NeighborsParti[lid].begin();it2!=NeighborsParti[lid].end();++it2){
        if(it2->first==hid){
            Cw=it2->second;//the weight value in the original graph
//            cout<<lid<<" "<<hid<<", Cw1: "<<Cw<<endl;
            flag=true;
            break;
        }
    }

    if(SCconNodesMTP[lid].find(hid) != SCconNodesMTP[lid].end()){//if found
        int widM;
        for(auto it=SCconNodesMTP[lid][hid].begin();it!=SCconNodesMTP[lid][hid].end();++it){
            ID=it->first;
            if(PartiTag[ID].second){//if boundary vertex
//                cout<<"Continue "<<ID<<endl;
                continue;
            }
            widM=IDMap[ID];
            int ssw=-1,wtt=-1;
            for(int j=0;j<Trees[PID][ranks[PID][widM]].vert.size();j++){
                if(Trees[PID][ranks[PID][widM]].vert[j].first==lid){
                    ssw=Trees[PID][ranks[PID][widM]].vert[j].second.first;
                    break;
                }

            }
            for(int j=0;j<Trees[PID][ranks[PID][widM]].vert.size();j++) {
                if (Trees[PID][ranks[PID][widM]].vert[j].first == hid) {
                    wtt = Trees[PID][ranks[PID][widM]].vert[j].second.first;
                    break;
                }
            }
            if(ssw==-1 || wtt==-1){
                cout<<"Wrong! "<<ssw<<" "<<wtt<<endl; exit(1);
            }

            if(ssw+wtt<Cw){
                Cw=ssw+wtt;
//                cout<<lid<<" "<<hid<<", Cw2: "<<Cw<<endl;
            }
        }
        //cout<<Cw<<endl;
        flag=true;
    }
//    else{
//        cout<<"Shortcut sc("<<ID1<<","<<ID2<<") does not exist!"<<endl;
//    }
    d=Cw;
    if(!flag){
        cout<<"Not found! "<<endl; exit(1);
    }
    return d;
}

void Graph::VertexInsertPartiBatchUpdateCH(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch){
    int a,b,oldW,newW;
    for(int i=0;i<wBatch.size();++i){
        a=wBatch[i].first.first; b=wBatch[i].first.second;
        oldW=wBatch[i].second.first; newW=wBatch[i].second.second;
//        cout<<a<<" "<<b<<" "<<oldW<<" "<<newW<<endl;
        // label inherit
        Node node;
        node.uniqueVertex=a;
        int bRank=ranks[pid][IDMap[b]];
//        cout<<Trees[pid][bRank].uniqueVertex<<"("<<Trees[pid][bRank].height<<"):";
//        for(int j=0;j<Trees[pid][bRank].vAncestor.size();++j){
//            cout<<" "<<Trees[pid][bRank].vAncestor[j]<<"("<<Trees[pid][bRank].dis[j]<<")";
//        }
//        cout<<endl;
        node.pa=bRank;
        node.height=Trees[pid][bRank].height+1;
        node.hdepth=node.hdepth;
        if(node.height>heightMaxs[pid]){
            heightMaxs[pid]=node.height;
        }
        node.vert.emplace_back(b,make_pair(newW,1));

        ranks[pid].push_back(Trees[pid].size());
        Trees[pid][bRank].ch.push_back(Trees[pid].size());
        Trees[pid].emplace_back(node);
//        int aRank=ranks[pid][IDMap[a]];
//        cout<<Trees[pid][aRank].uniqueVertex<<"("<<Trees[pid][aRank].height<<"):";
//        for(int j=0;j<Trees[pid][aRank].vAncestor.size();++j){
//            cout<<" "<<Trees[pid][aRank].vAncestor[j]<<"("<<Trees[pid][aRank].dis[j]<<")";
//        }
//        cout<<endl;
//        for(int j=0;j<Trees[pid][aRank].pos.size();++j){
//            cout<<Trees[pid][aRank].pos[j]<<" ";
//        }
//        cout<<endl;
        if(PSPStrategy>=PostBoundary){
            Node node;
            node.uniqueVertex=a;
            int bRank=ranksPost[pid][IDMap[b]];
            node.pa=bRank;
            node.height=TreesPost[pid][bRank].height+1;
            node.hdepth=node.hdepth;
            if(node.height>heightMaxsPost[pid]){
                heightMaxsPost[pid]=node.height;
            }
            node.vert.emplace_back(b,make_pair(newW,1));
            for(int j=0;j<TreesPost[pid][bRank].dis.size();++j){
                node.dis.push_back(TreesPost[pid][bRank].dis[j]+newW);
            }
            node.dis.push_back(newW);
            node.cnt=TreesPost[pid][bRank].cnt;
            node.cnt.push_back(1);
            node.vAncestor=TreesPost[pid][bRank].vAncestor;
            node.pos.push_back(node.vAncestor.size());
            node.vAncestor.push_back(b);
            node.pos.push_back(node.vAncestor.size());
            node.FN.assign(TreesPost[pid][bRank].FN.size(),false);
            node.FN.push_back(true);

            ranksPost[pid].push_back(TreesPost[pid].size());
            TreesPost[pid][bRank].ch.push_back(TreesPost[pid].size());
            TreesPost[pid].emplace_back(node);
        }
    }
//    cout<<"LCA index rebuild."<<endl;
    RMQIndexs[pid].clear();
    toRMQs[pid].clear();
    makeRMQCoreP(pid, toRMQs, RMQIndexs, Trees);
    if(PSPStrategy>=PostBoundary){
        RMQIndexsPost[pid].clear();
        toRMQsPost[pid].clear();
        makeRMQCoreP(pid, toRMQsPost, RMQIndexsPost, TreesPost);
    }
}

void Graph::VertexInsertPartiBatchUpdateH2H(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch){
    int a,b,oldW,newW;
    for(int i=0;i<wBatch.size();++i){
        a=wBatch[i].first.first; b=wBatch[i].first.second;
        oldW=wBatch[i].second.first; newW=wBatch[i].second.second;
//        cout<<a<<" "<<b<<" "<<oldW<<" "<<newW<<endl;
        // label inherit
        Node node;
        node.uniqueVertex=a;
        int bRank=ranks[pid][IDMap[b]];
//        cout<<Trees[pid][bRank].uniqueVertex<<"("<<Trees[pid][bRank].height<<"):";
//        for(int j=0;j<Trees[pid][bRank].vAncestor.size();++j){
//            cout<<" "<<Trees[pid][bRank].vAncestor[j]<<"("<<Trees[pid][bRank].dis[j]<<")";
//        }
//        cout<<endl;
        node.pa=bRank;
        node.height=Trees[pid][bRank].height+1;
        node.hdepth=node.hdepth;
        if(node.height>heightMaxs[pid]){
            heightMaxs[pid]=node.height;
        }
        node.vert.emplace_back(b,make_pair(newW,1));
        for(int j=0;j<Trees[pid][bRank].dis.size();++j){
            node.dis.push_back(Trees[pid][bRank].dis[j]+newW);
        }
        node.dis.push_back(newW);
        node.cnt=Trees[pid][bRank].cnt;
        node.cnt.push_back(1);
        node.vAncestor=Trees[pid][bRank].vAncestor;
        node.pos.push_back(node.vAncestor.size());
        node.vAncestor.push_back(b);
        node.pos.push_back(node.vAncestor.size());
        node.FN.assign(Trees[pid][bRank].FN.size(),false);
        node.FN.push_back(true);

        ranks[pid].push_back(Trees[pid].size());
        Trees[pid][bRank].ch.push_back(Trees[pid].size());
        Trees[pid].emplace_back(node);
//        int aRank=ranks[pid][IDMap[a]];
//        cout<<Trees[pid][aRank].uniqueVertex<<"("<<Trees[pid][aRank].height<<"):";
//        for(int j=0;j<Trees[pid][aRank].vAncestor.size();++j){
//            cout<<" "<<Trees[pid][aRank].vAncestor[j]<<"("<<Trees[pid][aRank].dis[j]<<")";
//        }
//        cout<<endl;
//        for(int j=0;j<Trees[pid][aRank].pos.size();++j){
//            cout<<Trees[pid][aRank].pos[j]<<" ";
//        }
//        cout<<endl;
        if(PSPStrategy>=PostBoundary){
            Node node;
            node.uniqueVertex=a;
            int bRank=ranksPost[pid][IDMap[b]];
            node.pa=bRank;
            node.height=TreesPost[pid][bRank].height+1;
            node.hdepth=node.hdepth;
            if(node.height>heightMaxsPost[pid]){
                heightMaxsPost[pid]=node.height;
            }
            node.vert.emplace_back(b,make_pair(newW,1));
            for(int j=0;j<TreesPost[pid][bRank].dis.size();++j){
                node.dis.push_back(TreesPost[pid][bRank].dis[j]+newW);
            }
            node.dis.push_back(newW);
            node.cnt=TreesPost[pid][bRank].cnt;
            node.cnt.push_back(1);
            node.vAncestor=TreesPost[pid][bRank].vAncestor;
            node.pos.push_back(node.vAncestor.size());
            node.vAncestor.push_back(b);
            node.pos.push_back(node.vAncestor.size());
            node.FN.assign(TreesPost[pid][bRank].FN.size(),false);
            node.FN.push_back(true);

            ranksPost[pid].push_back(TreesPost[pid].size());
            TreesPost[pid][bRank].ch.push_back(TreesPost[pid].size());
            TreesPost[pid].emplace_back(node);
        }
    }
//    cout<<"LCA index rebuild."<<endl;
    RMQIndexs[pid].clear();
    toRMQs[pid].clear();
    makeRMQCoreP(pid, toRMQs, RMQIndexs, Trees);
    if(PSPStrategy>=PostBoundary){
        RMQIndexsPost[pid].clear();
        toRMQsPost[pid].clear();
        makeRMQCoreP(pid, toRMQsPost, RMQIndexsPost, TreesPost);
    }
}
