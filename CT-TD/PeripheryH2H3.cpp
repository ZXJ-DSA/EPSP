/*
 * PeripheryH2H.cpp
 *
 *  Created on: 14 Oct 2022
 *      Author: zhangmengxuan
 */
#include "head3.h"

extern bool ifExtension;

//Function of constructing tree decomposition
void Graph::Construct_tree(bool ifParallel){
    H2HContract();//MDE-based contraction
    Create_tree();//Create tree
    Compute_tree_label(ifParallel);//Construct periphery index (H2H label + interface label)
}

void Graph::Construct_core(int indexType) {

    //construct the index of Core
    cout<<"Begin Core's Index Construction..."<<endl;
    if(indexType == 1){
        //****************PSL construction***************************
        this->indexType=1;
        cout<<"PSL"<<endl;
        vSm.reserve(nodenum);
        for(int i = 0; i < nodenum; i++)
        {
            Semaphore* s = new Semaphore(1);
            vSm.push_back(s);
        }
        IndexConstructMThread2New();

        //***********************************************************
    }else if(indexType == 0){
        //****************PLL construction***************************
        this->indexType=0;
        cout<<"PLL"<<endl;
        PLLIndexConstruct();
        //***********************************************************
    }else if(indexType == 2){
        this->indexType = 2;
        cout<<"Batch PSL"<<endl;
        vSm.reserve(nodenum);
        for(int i = 0; i < nodenum; i++)
        {
            Semaphore* s = new Semaphore(1);
            vSm.push_back(s);
        }
        BVCIndexConstructMThread();
    }



//    CorrectnessCheckCore();
}

void Graph::OverlayGraphProcess(){
    //core's adjacent complete through B*B
    /*for(int pid=0;pid<partiNum;pid++){//for each partition
        int ID1,ID2,dis;
        for(int i=0;i<BoundVertex[pid].size();i++){//for each boundary vertex
            ID1=BoundVertex[pid][i];
//            if(AdjaParti[pid][ID1].size()!=0){
            if(AdjaGraph[ID1].size()!=0){
                //if(PartiVertexInverted[pid].find(ID1)!=PartiVertexInverted[pid].end()){//only if ID1 has neighbors in partition pid
                for(int j=i+1;j<BoundVertex[pid].size();j++){
                    ID2=BoundVertex[pid][j];
                    if(AdjaGraph[ID2].size()!=0){
                        //if(PartiVertexInverted[pid].find(ID2)!=PartiVertexInverted[pid].end()){
                        dis=QueryH2HPartition(ID1,ID2,pid);//query between two boundary vertex should be answered by core labels

                        if(dis<INF){
                            SuppPartiID[ID1][ID2].insert(pid);
                            SuppPartiID[ID2][ID1].insert(pid);


                            if(AdjaCoreMap[ID1].find(ID2)!=AdjaCoreMap[ID1].end()){
                                if(AdjaCoreMap[ID1][ID2]>dis){
                                    cout<<"Changed. "<<AdjaCoreMap[ID1][ID2]<<" "<<dis<<endl;
                                    AdjaCoreMap[ID1][ID2]=dis;
                                    AdjaCoreMap[ID2][ID1]=dis;
                                    SuppPartiIDReal[ID1][ID2].clear();
                                    SuppPartiIDReal[ID2][ID1].clear();
                                    SuppPartiIDReal[ID1][ID2].insert(pid);
                                    SuppPartiIDReal[ID2][ID1].insert(pid);
                                }else if(AdjaCoreMap[ID1][ID2]==dis){
                                    SuppPartiIDReal[ID1][ID2].insert(pid);
                                    SuppPartiIDReal[ID2][ID1].insert(pid);
                                }
                            }else{
//                                AdjaCoreMap[ID1].insert(make_pair(ID2,dis));
//                                AdjaCoreMap[ID2].insert(make_pair(ID1,dis));
                                SuppPartiIDReal[ID1][ID2].insert(pid);
                                SuppPartiIDReal[ID2][ID1].insert(pid);
                            }
                        }

                    }
                }
            }
        }
    }*/

    /*vector<pair<int,int>> vecpair;
    vecpair.clear();
    AdjaCore.assign(nodenum,vecpair);

    for(int i=0;i<nodenum;++i){
        if(CoreTag[i]==-1){//if core vertex
            for(auto it=Emap[i].begin();it!=Emap[i].end();++it){
                if(CoreTag[it->first]==-1){
                    AdjaCore[i].push_back(*it);
                    AdjaCore[it->first].emplace_back(i,it->second);
                }
            }
        }
    }


    for(int i=0;i<AdjaCoreMap.size();i++){
        if(AdjaCoreMap[i].empty())
            continue;

        for(map<int,int>::iterator it=AdjaCoreMap[i].begin();it!=AdjaCoreMap[i].end();it++){
            AdjaCore[i].push_back(make_pair((*it).first, (*it).second));
        }
//        if(AdjaGraph[i].empty()){
//            AdjaGraph[i] = AdjaCore[i];
//        }else{
//            cout<<"Wrong!! AdjaGraph["<<i<<"] is not empty!"<<endl;
//            exit(1);
//        }
    }*/


//    Neighbors = AdjaCore;
//    cout<<"***************Finish Core's graph construction***************"<<endl;
//    CorrectnessCheckCore();
}

//Query processing for debug
int Graph::QueryDebug(int ID1, int ID2){
    int dis=INF;

    if(CoreTag[ID1]==-1 && CoreTag[ID2]==-1){//Case 1: both in core
        cout<<"Core-Core"<<endl;
        dis=QueryCoreDebug(ID1, ID2);

    }else if(CoreTag[ID1]==-1 && CoreTag[ID2]!=-1){//Case 2: ID2 in partition, ID1 in core
        cout<<"Core-Parti"<<endl;
        if(ifExtension){
            dis=QueryPartiCoreExtDebug(ID2, ID1);
        }else{
            dis=QueryPartiCoreDebug(ID2, ID1);
        }
    }else if(CoreTag[ID1]!=-1 && CoreTag[ID2]==-1){//Case 2: ID1 in partition, ID2 in core
        cout<<"Parti-Core"<<endl;
        if(ifExtension){
            dis=QueryPartiCoreExtDebug(ID1, ID2);
        }else{
            dis=QueryPartiCoreDebug(ID1, ID2);
        }
    }else if(CoreTag[ID1]!=-1 && CoreTag[ID2]!=-1){//both in partition

        if(ifExtension){
            dis=QueryPartiPartiExtDebug(ID1,ID2);
        }else{
//            dis=QueryPartiParti(ID1,ID2);
            if(CoreTag[ID1] != CoreTag[ID2]){//Case 3: in different peripheries
                cout<<"Parti-Parti"<<endl;
                int d=INF;
                int b1,b2,d1,d2;//final results
                int pid1=CoreTag[ID1];
                int pid2=CoreTag[ID2];

                vector<int> B1=BoundVertex[pid1];
                vector<int> B2=BoundVertex[pid2];

                map<int,int> m1,m2;
                m1.clear();
                m2.clear();
                int bID1, bID2, tempdis;
                for(int i=0;i<B1.size()-1;i++){
                    bID1=B1[i];
                    assert(Tree[rank[ID1]].disInf.find(bID1)!=Tree[rank[ID1]].disInf.end());
//            m1.insert(make_pair(bID1,Tree[rank[ID1]].disInf[i]));
                    m1.insert(make_pair(bID1,Tree[rank[ID1]].disInf[bID1]));
                }
                for(int j=0;j<B2.size()-1;j++){
                    bID2=B2[j];
                    assert(Tree[rank[ID2]].disInf.find(bID2)!=Tree[rank[ID2]].disInf.end());
//            m2.insert(make_pair(bID2,Tree[rank[ID2]].disInf[j]));
                    m2.insert(make_pair(bID2,Tree[rank[ID2]].disInf[bID2]));
                }

                for(int k=0;k<B1.size()-1;k++){
                    bID1=B1[k];

                    if(m1[bID1]>d)
                        continue;

                    for(int z=0;z<B2.size()-1;z++){
                        bID2=B2[z];

                        if(m2[bID2]>d)
                            continue;

                        tempdis=m1[bID1]+QueryCore(bID1,bID2)+m2[bID2];
                        if(tempdis<d){
                            d=tempdis;
                            b1=bID1; b2=bID2; d1=m1[bID1]; d2=m2[bID2];
                        }
                    }
                }
                dis=d;
                int d_12=QueryCore(b1,b2), dDijk_s=Dijkstra(ID1,b1,Neighbor), dDijk_12=Dijkstra(b1,b2,Neighbor), dDijk_t=Dijkstra(b2,ID2,Neighbor);
                cout<<ID1<<" "<<b1<<"("<<NodeOrder[b1]<<") "<<b2<<"("<<NodeOrder[b2]<<") "<<ID2<<" : "<<d1<<" "<<d_12<<" "<<d2<<" ; "<<dDijk_s<<" "<<dDijk_12<<"("<<DijkstraCore(b1,b2)<<") "<<dDijk_t<<endl;

//                if(d1!=dDijk_s){
//                    DijkstraPath(ID1,b1);
//                }
//                if(d_12!=dDijk_12){
//                    DijkstraPath(b1,b2);
//                }
//                if(d2!=dDijk_t){
//                    DijkstraPath(b2,ID2);
//                }

            }else{//Case 4: in the same periphery
                cout<<"Same-Parti"<<endl;
//                dis= QuerySameParti(ID1,ID2);
                int d=INF;
                int b1,b2,df1,df2;
                int pid1=CoreTag[ID1];
                int pid2=CoreTag[ID2];

                int temp_dis = QueryPeripheryTree(ID1, ID2, pid1);/// d2 may be wrong sometimes
                if (temp_dis < d){
                    d = temp_dis;//QueryH2HPartition(ID1,ID2,pid1);
                    b1=b2=-1;
                    df1=df2=-1;
                }

                vector<int> B = BoundVertex[pid1];
                map<int, int> m1, m2;
                m1.clear();
                m2.clear();
                vector<int> B1, B2;
                B1.clear();
                B2.clear();
                int bID, d1, d2;
                for (int i = 0; i < B.size() - 1; i++) {
                    bID = B[i];
                    assert(Tree[rank[ID1]].disInf.find(bID) != Tree[rank[ID1]].disInf.end());
                    assert(Tree[rank[ID2]].disInf.find(bID) != Tree[rank[ID2]].disInf.end());
                    d1 = Tree[rank[ID1]].disInf[bID];
                    d2 = Tree[rank[ID2]].disInf[bID];

                    if (d1 < d) {
                        B1.push_back(bID);
                        m1.insert(make_pair(bID, d1));
                    }
                    if (d2 < d) {
                        B2.push_back(bID);
                        m2.insert(make_pair(bID, d2));
                    }
                }

                int bID1, bID2, tempdis;
                if (!B1.empty() && !B2.empty()) {
                    for (int k = 0; k < B1.size(); k++) {
                        bID1 = B1[k];
                        if (m1[bID1] > d)
                            continue;
                        for (int z = 0; z < B2.size(); z++) {
                            bID2 = B2[z];
                            if (m2[bID2] > d)
                                continue;
                            tempdis = m1[bID1] + QueryCore(bID1, bID2) + m2[bID2];
                            if (tempdis < d){
                                d = tempdis;
                                b1=bID1;b2=bID2;
                                df1=m1[bID1];df2=m2[bID2];
                            }
                        }
                    }
                }

                if(b1!=-1){
                    cout<<"d4: "<<ID1<<" "<<b1<<" "<<b2<<" "<<ID2<<" : "<<df1<<" "<<QueryCore(b1,b2)<<" "<<df2<<" ; "<<Dijkstra(ID1,b1,Neighbor)<<" "<<Dijkstra(b1,b2,Neighbor)<<" "<<Dijkstra(b2,ID2,Neighbor)<<endl;
                }else{
                    int dDijk2 = Dijkstra(ID1,ID2,Neighbor);
                    cout<<"d2: "<<d<<"; "<<dDijk2<<endl;
                    if(d!=dDijk2){
//                        DijkstraPath(ID1,ID2);
                    }
                }

                dis = d;

            }
        }
    }
    return dis;
}

//Query processing
int Graph::Query(int ID1, int ID2){
    int dis=INF;

    if(CoreTag[ID1]==-1 && CoreTag[ID2]==-1){//Case 1: both in core
//        cout<<"Core-Core"<<endl;
        dis=QueryCore(ID1, ID2);
    }else if(CoreTag[ID1]==-1 && CoreTag[ID2]!=-1){//Case 2: ID2 in partition, ID1 in core
//        cout<<"Core-Parti"<<endl;
        if(ifExtension){
            dis=QueryPartiCoreExt(ID2, ID1);
        }else{
            dis=QueryPartiCore(ID2, ID1);
        }
    }else if(CoreTag[ID1]!=-1 && CoreTag[ID2]==-1){//Case 2: ID1 in partition, ID2 in core
//        cout<<"Parti-Core"<<endl;
        if(ifExtension){
            dis=QueryPartiCoreExt(ID1, ID2);
        }else{
            dis=QueryPartiCore(ID1, ID2);
        }
    }else if(CoreTag[ID1]!=-1 && CoreTag[ID2]!=-1){//both in partition

        if(ifExtension){
            dis=QueryPartiPartiExt(ID1,ID2);
        }else{
            if(CoreTag[ID1] != CoreTag[ID2]){//Case 3: in different peripheries
                dis= QueryPartiParti(ID1,ID2);
            }else{//Case 4: in the same periphery
                dis= QuerySameParti(ID1,ID2);
            }
        }
    }
    return dis;
}

int Graph::QueryPartiCore(int ID1, int ID2){//ID1 partition, ID2 core
	int d=INF;

	int pid=CoreTag[ID1];
	int bid;
	int dis1,dis2;

//    for(int k=0;k<Tree[rank[ID1]].disInf.size();++k){
//        bid=BoundVertex[pid][k];
//        dis1=Tree[rank[ID1]].disInf[k];
//        dis2= QueryCore(bid,ID2);
//        if(d>dis1+dis2)
//            d=dis1+dis2;
//    }
    for(auto it=Tree[rank[ID1]].disInf.begin();it!=Tree[rank[ID1]].disInf.end();++it){
        bid=it->first;
        dis1=it->second;
        dis2= QueryCore(bid,ID2);
        if(d>dis1+dis2)
            d=dis1+dis2;
    }
	return d;
}

int Graph::QueryPartiCoreDebug(int ID1, int ID2){//ID1 partition, ID2 core
    int d=INF;

    assert(CoreTag[ID1]>=0 && CoreTag[ID2]==-1);

    int pid=CoreTag[ID1];
    int bid;
    int dis1,dis2;
    int bfinal=-1,d1=INF,d2=INF;

    for(auto it=Tree[rank[ID1]].disInf.begin();it!=Tree[rank[ID1]].disInf.end();++it){
        bid=it->first;
        dis1=it->second;
        dis2= QueryCore(bid,ID2);
        if(d>dis1+dis2){
            d=dis1+dis2;
            bfinal = bid;
            d1 = dis1; d2 = dis2;
        }

    }

    int dDijk_s=Dijkstra(ID1,bfinal,Neighbor);
    int dDijk_t=Dijkstra(bfinal,ID2,Neighbor);
    cout<<ID1<<" "<<bfinal<<"("<<NodeOrder[bfinal]<<") "<<ID2<<"("<<NodeOrder[ID2]<<") : "<<d1<<" "<<d2<<" ; "<<dDijk_s<<" "<<dDijk_t<<"("<<DijkstraCore(bfinal,ID2)<<") "<<endl;

    if(d1!=dDijk_s){
        DijkstraPath(ID1,bfinal);
    }

    if(d2!=dDijk_t){
        DijkstraPath(bfinal,ID2);
    }

    return d;
}
//function of constructing the extended index for given pid
void Graph::ExtensionIndex(pair<int,int> pidRange, bool ifIncrease){
    int ID1;
    int b1, b2;
    int temp_dis=INF;

    for(int pid=pidRange.first;pid<pidRange.second;++pid){
        vector<int> B=BoundVertex[pid];
        for(auto it=PartiVertex[pid].begin();it!=PartiVertex[pid].end();++it){
            ID1 = *it;
            if(BoundVertexSet[pid].find(ID1) != BoundVertexSet[pid].end()){//if found, i.e., it is boundary vertex
                cout<<"Wrong! "<<ID1<<" of partition "<<pid<<" should not be a boundary vertex!"<<endl;
                exit(1);
            }
            if(ifIncrease){
                for(auto it2=IndexExt[ID1].begin();it2!=IndexExt[ID1].end();++it2){
                    IndexExt[ID1][it2->first] = INF;
                }
//                IndexExt[ID1].clear();
            }
            // for each non-boundary vertex
            for(int i=0;i<B.size()-1;++i){
                b1 = B[i];
                for(auto it2 = Label[b1].begin();it2!= Label[b1].end();++it2){//for each label of b1
                    b2 = it2->first;
                    if(IndexExt[ID1].find(b2) == IndexExt[ID1].end()){//if not found
                        IndexExt[ID1][b2] = Tree[rank[ID1]].disInf[b1] + it2->second;
                    }else{
                        temp_dis = Tree[rank[ID1]].disInf[b1] + it2->second;
                        if(IndexExt[ID1][b2] > temp_dis){
                            IndexExt[ID1][b2] = temp_dis;
                        }
                    }
                }

            }
        }
    }

}
//function of constructing the extended index for given pid
void Graph::ExtensionIndex2(pair<int,int> pidRange, bool ifIncrease, vector<int>& partiForUpdate){
    int ID1;
    int b1, b2;
    int temp_dis=INF;

    for(int id=pidRange.first;id<pidRange.second;++id){
        int pid = partiForUpdate[id];
        vector<int> B=BoundVertex[pid];
        for(auto it=PartiVertex[pid].begin();it!=PartiVertex[pid].end();++it){
            ID1 = *it;
            if(BoundVertexSet[pid].find(ID1) != BoundVertexSet[pid].end()){//if found, i.e., it is boundary vertex
                cout<<"Wrong! "<<ID1<<" of partition "<<pid<<" should not be a boundary vertex!"<<endl;
                exit(1);
            }
            if(ifIncrease){
                for(auto it2=IndexExt[ID1].begin();it2!=IndexExt[ID1].end();++it2){
                    IndexExt[ID1][it2->first] = INF;
                }
//                IndexExt[ID1].clear();
            }
            // for each non-boundary vertex
            for(int i=0;i<B.size()-1;++i){
                b1 = B[i];
                assert(Tree[rank[ID1]].disInf.find(b1)!=Tree[rank[ID1]].disInf.end());
                for(auto it2 = Label[b1].begin();it2!= Label[b1].end();++it2){//for each label of b1
                    b2 = it2->first;
                    if(IndexExt[ID1].find(b2) == IndexExt[ID1].end()){//if not found
                        IndexExt[ID1][b2] = Tree[rank[ID1]].disInf[b1] + it2->second;
                    }else{
                        temp_dis = Tree[rank[ID1]].disInf[b1] + it2->second;
                        if(IndexExt[ID1][b2] > temp_dis){
                            IndexExt[ID1][b2] = temp_dis;
                        }
                    }
                }

            }
        }
    }

}

//function for constructing the extension index for periphery vertex
void Graph::ExtensionIndexConstruct(bool ifParallel, bool ifIncrease){
    int ID1;
    int b1, b2;
    int temp_dis=INF;
    IndexExt.assign(nodenum,unordered_map<int,int>());
    if(ifParallel){//use multi-thread
    //multiple thread
        if(partiNum>threadnum){
            int step=partiNum/threadnum;
            boost::thread_group threadf;
            for(int i=0;i<threadnum;i++){
                pair<int,int> p;
                p.first=i*step;
                if(i==threadnum-1)
                    p.second=partiNum;
                else
                    p.second=(i+1)*step;
                threadf.add_thread(new boost::thread(&Graph::ExtensionIndex, this, p, ifIncrease));
            }
            threadf.join_all();
        }else{
            boost::thread_group threadf;
            for(int pid=0;pid<partiNum;++pid) {
                threadf.add_thread(new boost::thread(&Graph::ExtensionIndex, this, make_pair(pid,pid+1), ifIncrease));
            }
            threadf.join_all();
        }
    }
    else{//single thread
        for(int pid=0;pid<partiNum;++pid){//for each partition
            vector<int> B=BoundVertex[pid];
            for(auto it=PartiVertex[pid].begin();it!=PartiVertex[pid].end();++it){//for each periphery vertex (excluding interface vertex)
                ID1 = *it;
                if(BoundVertexSet[pid].find(ID1) != BoundVertexSet[pid].end()){//if found, i.e., it is boundary vertex
                    cout<<"Wrong! "<<ID1<<" of partition "<<pid<<" should not be a boundary vertex!"<<endl;
                    exit(1);
                }
                if(ifIncrease){
                    for(auto it2=IndexExt[ID1].begin();it2!=IndexExt[ID1].end();++it2){
                        IndexExt[ID1][it2->first] = INF;
                    }
//                IndexExt[ID1].clear();
                }
                // for each non-boundary vertex
                for(int i=0;i<B.size()-1;++i){
                    b1 = B[i];
                    for(auto it2 = Label[b1].begin();it2!= Label[b1].end();++it2){//for each label of b1
                        b2 = it2->first;//the hub
                        if(IndexExt[ID1].find(b2) == IndexExt[ID1].end()){//if not found
                            IndexExt[ID1][b2] = Tree[rank[ID1]].disInf[b1] + it2->second;
                        }else{
                            temp_dis = Tree[rank[ID1]].disInf[b1] + it2->second;
                            if(IndexExt[ID1][b2] > temp_dis){
                                IndexExt[ID1][b2] = temp_dis;
                            }
                        }
                    }

                }
            }
        }
    }

}
//function for updating the extension index for periphery vertex
void Graph::ExtensionIndexUpdate(bool ifParallel, bool ifIncrease, vector<int>& partiForUpdate){
    int ID1;
    int b1, b2;
    int temp_dis=INF;
    int pNum = partiForUpdate.size();

    if(ifParallel){//use multi-thread
        //multiple thread
        if(pNum>threadnum){
            int step=pNum/threadnum;
            boost::thread_group threadf;
            for(int i=0;i<threadnum;i++){
                pair<int,int> p;
                p.first=i*step;
                if(i==threadnum-1)
                    p.second=pNum;
                else
                    p.second=(i+1)*step;
                threadf.add_thread(new boost::thread(&Graph::ExtensionIndex2, this, p, ifIncrease,partiForUpdate));
            }
            threadf.join_all();
        }else{
            boost::thread_group threadf;
            for(int pid=0;pid<pNum;++pid) {
                threadf.add_thread(new boost::thread(&Graph::ExtensionIndex2, this, make_pair(pid,pid+1), ifIncrease,partiForUpdate));
            }
            threadf.join_all();
        }
    }
    else{//single thread
        for(int i=0;i<pNum;++i){//for each partition
            int pid = partiForUpdate[i];
            vector<int> B=BoundVertex[pid];
            bool flag=true;
            for(auto it=PartiVertex[pid].begin();it!=PartiVertex[pid].end();++it){//for each periphery vertex (excluding interface vertex)
                ID1 = *it;
                if(BoundVertexSet[pid].find(ID1) != BoundVertexSet[pid].end()){//if found, i.e., it is boundary vertex
                    cout<<"Wrong! "<<ID1<<" of partition "<<pid<<" should not be a boundary vertex!"<<endl;
                    exit(1);
                }
                if(ifIncrease){
//                    if(flag){
//                        cout<<"Increase partition "<<pid<<endl;
//                        flag = false;
//                    }

                    for(auto it2=IndexExt[ID1].begin();it2!=IndexExt[ID1].end();++it2){
                        IndexExt[ID1][it2->first] = INF;
                    }
//                IndexExt[ID1].clear();
                }
                // for each non-boundary vertex
                for(int i=0;i<B.size()-1;++i){
                    b1 = B[i];
                    assert(Tree[rank[ID1]].disInf.find(b1)!=Tree[rank[ID1]].disInf.end());
                    for(auto it2 = Label[b1].begin();it2!= Label[b1].end();++it2){//for each label of b1
                        b2 = it2->first;//the hub
                        if(IndexExt[ID1].find(b2) == IndexExt[ID1].end()){//if not found
                            IndexExt[ID1][b2] = Tree[rank[ID1]].disInf[b1] + it2->second;
                        }else{
                            temp_dis = Tree[rank[ID1]].disInf[b1] + it2->second;
                            if(IndexExt[ID1][b2] > temp_dis){
                                IndexExt[ID1][b2] = temp_dis;
                            }
                        }
                    }

                }
            }
        }
    }

}
//function of processing the query between in-partition vertex and core vertex by extended index
int Graph::QueryPartiCoreExt(int ID1, int ID2){//ID1 is in partition, ID2 is in core
    int d=INF;

    bool flag_union=false;

    unordered_map<int,int>::iterator it;
    int hub, dis1, dis2;
    int hubfinal=-1,dis1final=INF,dis2final=INF;
//    cout << "ID1 label Size:" << IndexExt[ID1].size() << "\tID2 label Size:" <<Label[ID2].size()<<endl;
    for(it=IndexExt[ID1].begin();it!=IndexExt[ID1].end();it++){
        hub=(*it).first;
        dis1=(*it).second;
        if(Label[ID2].find(hub)!=Label[ID2].end()){
            flag_union = true;
            dis2=Label[ID2][hub];
            if(dis1+dis2<d){
                d=dis1+dis2;
                hubfinal=hub;
                dis1final=dis1;
                dis2final=dis2;
                //cout<<"hub "<<hub<<",dis "<<d<<endl;
                //cout<<"check labeling "<<dis1<<" "<<DijkstraCore(ID1,hub)<<" "<<dis2<<" "<<DijkstraCore(ID1,hub)<<endl;
            }
        }
    }
    if(!flag_union){
        cout<<"There is no intersection of these two vertex label!!!"<<endl;
    }
    if(d == INF){
        cout<<"Wrong!!! "<<ID1<<" "<<ID2<<" "<<INF<<endl;
    }
    return d;
}
//function of processing the query between in-partition vertex and core vertex by extended index
int Graph::QueryPartiCoreExtDebug(int ID1, int ID2){//ID1 is in partition, ID2 is in core
    int d=INF;

    unordered_map<int,int>::iterator it;
    int hub, dis1, dis2;
    int hubfinal=-1,dis1final=INF,dis2final=INF;
//    cout << "ID1 label Size:" << IndexExt[ID1].size() << "\tID2 label Size:" <<Label[ID2].size()<<endl;
    for(it=IndexExt[ID1].begin();it!=IndexExt[ID1].end();it++){
        hub=(*it).first;
        dis1=(*it).second;
        if(Label[ID2].find(hub)!=Label[ID2].end()){
            dis2=Label[ID2][hub];
            if(dis1+dis2<d){
                d=dis1+dis2;
                hubfinal=hub;
                dis1final=dis1;
                dis2final=dis2;
                //cout<<"hub "<<hub<<",dis "<<d<<endl;
                //cout<<"check labeling "<<dis1<<" "<<DijkstraCore(ID1,hub)<<" "<<dis2<<" "<<DijkstraCore(ID1,hub)<<endl;
            }
        }
    }
    if(hubfinal != -1){
        int dDijk_s=Dijkstra(ID1,hubfinal,Neighbor), dDijk_t=Dijkstra(hubfinal,ID2,Neighbor);
        cout<<ID1<<" "<<hubfinal<<" "<<ID2<<" : "<<dis1final<<" "<<dis2final<<" ; "<<dDijk_s<<" "<<dDijk_t<<endl;
    }
    if(d == INF){
        cout<<"Wrong!!! "<<ID1<<" "<<ID2<<" "<<INF<<endl;
    }
    return d;
}

//function of processing the query between two in-partition vertices by extended index
int Graph::QueryPartiPartiExt(int ID1, int ID2){
    int d=INF;

    int pid1=CoreTag[ID1];
    int pid2=CoreTag[ID2];
    if(pid1==pid2){//if in the same partition
        int temp_dis = QueryPeripheryTree(ID1,ID2,pid1);//d2
        if(temp_dis<d)//QueryH2HPartition(ID1,ID2,pid1)
            d=temp_dis;//QueryH2HPartition(ID1,ID2,pid1);
        unordered_map<int,int>::iterator it;
        int hub, dis1, dis2;
        int hubfinal,dis1final,dis2final;
        for(it=IndexExt[ID1].begin();it!=IndexExt[ID1].end();it++){
            hub=(*it).first;
            dis1=(*it).second;
            if(IndexExt[ID2].find(hub)!=IndexExt[ID2].end()){
                dis2=IndexExt[ID2][hub];
                if(dis1+dis2<d){
                    d=dis1+dis2;
                    hubfinal=hub;
                    dis1final=dis1;
                    dis2final=dis2;
                    //cout<<"hub "<<hub<<",dis "<<d<<endl;
                    //cout<<"check labeling "<<dis1<<" "<<DijkstraCore(ID1,hub)<<" "<<dis2<<" "<<DijkstraCore(ID1,hub)<<endl;
                }
            }
        }
    }else{//if in different partitions
        unordered_map<int,int>::iterator it;
        int hub, dis1, dis2;
        int hubfinal,dis1final,dis2final;
        for(it=IndexExt[ID1].begin();it!=IndexExt[ID1].end();it++){
            hub=(*it).first;
            dis1=(*it).second;
            if(IndexExt[ID2].find(hub)!=IndexExt[ID2].end()){
                dis2=IndexExt[ID2][hub];
                if(dis1+dis2<d){
                    d=dis1+dis2;
                    hubfinal=hub;
                    dis1final=dis1;
                    dis2final=dis2;
                    //cout<<"hub "<<hub<<",dis "<<d<<endl;
                    //cout<<"check labeling "<<dis1<<" "<<DijkstraCore(ID1,hub)<<" "<<dis2<<" "<<DijkstraCore(ID1,hub)<<endl;
                }
            }
        }
    }
    if(d == INF){
        cout<<"Wrong!!! "<<ID1<<" "<<ID2<<" "<<INF<<endl;
    }
    return d;
}
//function of processing the query between two in-partition vertices by extended index
int Graph::QueryPartiPartiExtDebug(int ID1, int ID2){
    int d=INF;

    int pid1=CoreTag[ID1];
    int pid2=CoreTag[ID2];
    if(pid1==pid2){//if in the same partition
        cout<<"Same-Parti"<<endl;
        int temp_dis = QueryPeripheryTree(ID1,ID2,pid1);//d2
        if(temp_dis<d)//QueryH2HPartition(ID1,ID2,pid1)
            d=temp_dis;//QueryH2HPartition(ID1,ID2,pid1);
        unordered_map<int,int>::iterator it;
        int hub, dis1, dis2;
        int hubfinal=-1,dis1final=INF,dis2final=INF;
        for(it=IndexExt[ID1].begin();it!=IndexExt[ID1].end();it++){
            hub=(*it).first;
            dis1=(*it).second;
            if(IndexExt[ID2].find(hub)!=IndexExt[ID2].end()){
                dis2=IndexExt[ID2][hub];
                if(dis1+dis2<d){
                    d=dis1+dis2;//d4
                    hubfinal=hub;
                    dis1final=dis1;
                    dis2final=dis2;
                }
            }
        }
        if(hubfinal!=-1){
            cout<<"d4: "<<ID1<<" "<<hubfinal<<" "<<ID2<<" : "<<dis1final<<" "<<dis2final<<" ; "<<Dijkstra(ID1,hubfinal,Neighbor)<<" "<<Dijkstra(hubfinal,ID2,Neighbor)<<endl;
        }else{
            int dDijk2 = Dijkstra(ID1,ID2,Neighbor);
            cout<<"d2: "<<d<<"; "<<dDijk2<<endl;
//            if(d!=dDijk2){
//                DijkstraPath(ID1,ID2);
//            }
        }

    }else{//if in different partitions
        cout<<"Parti-Parti"<<endl;
        unordered_map<int,int>::iterator it;
        int hub, dis1, dis2;
        int hubfinal=-1,dis1final=INF,dis2final=INF;
        for(it=IndexExt[ID1].begin();it!=IndexExt[ID1].end();it++){
            hub=(*it).first;
            dis1=(*it).second;
            if(IndexExt[ID2].find(hub)!=IndexExt[ID2].end()){
                dis2=IndexExt[ID2][hub];
                if(dis1+dis2<d){
                    d=dis1+dis2;
                    hubfinal=hub;
                    dis1final=dis1;
                    dis2final=dis2;
                    //cout<<"hub "<<hub<<",dis "<<d<<endl;
                    //cout<<"check labeling "<<dis1<<" "<<DijkstraCore(ID1,hub)<<" "<<dis2<<" "<<DijkstraCore(ID1,hub)<<endl;
                }
            }
        }
        if(hubfinal != -1){
            int dDijk_s=Dijkstra(ID1,hubfinal,Neighbor), dDijk_t=Dijkstra(hubfinal,ID2,Neighbor);
            cout<<ID1<<" "<<hubfinal<<" "<<ID2<<" : "<<dis1final<<" "<<dis2final<<" ; "<<dDijk_s<<" "<<dDijk_t<<endl;

//                if(d1!=dDijk_s){
//                    DijkstraPath(ID1,b1);
//                }
//                if(d_12!=dDijk_12){
//                    DijkstraPath(b1,b2);
//                }
//                if(d2!=dDijk_t){
//                    DijkstraPath(b2,ID2);
//                }
        }else{
            cout<<"Wrong! There is no hub!"<<endl;
        }

    }
    if(d == INF){
        cout<<"Wrong!!! "<<ID1<<" "<<ID2<<" "<<INF<<endl;
    }
    return d;
}

int Graph::QuerySameParti(int ID1, int ID2){//both are within the partition
    int d=INF;

    int pid1=CoreTag[ID1];
    int pid2=CoreTag[ID2];
    if(pid1==pid2){//if in the same partition
//        cout<<"Same-Parti"<<endl;
        int temp_dis = QueryPeripheryTree(ID1,ID2,pid1);/// d2 may be wrong sometimes
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
        for(int i=0;i<B.size()-1;i++){
            bID=B[i];
            assert(Tree[rank[ID1]].disInf.find(bID)!=Tree[rank[ID1]].disInf.end());
            assert(Tree[rank[ID2]].disInf.find(bID)!=Tree[rank[ID2]].disInf.end());
//            d1=Tree[rank[ID1]].disInf[i];
//            d2=Tree[rank[ID2]].disInf[i];
            d1=Tree[rank[ID1]].disInf[bID];
            d2=Tree[rank[ID2]].disInf[bID];

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

int Graph::QueryPartiParti(int ID1, int ID2){//both are within the partition
    int d=INF;

    int pid1=CoreTag[ID1];
    int pid2=CoreTag[ID2];
    if(pid1==pid2){//if in the same partition
        cout<<"Wrong for partition-partition query!"<<endl;
        exit(1);

    }else{//if in different partitions
//        cout<<"Parti-Parti"<<endl;
        vector<int> B1=BoundVertex[pid1];
        vector<int> B2=BoundVertex[pid2];

        map<int,int> m1,m2;
        m1.clear();
        m2.clear();
        int bID1, bID2, tempdis;
        int b1,b2,d1,d2;
        for(int i=0;i<B1.size()-1;i++){
            bID1=B1[i];
            assert(Tree[rank[ID1]].disInf.find(bID1)!=Tree[rank[ID1]].disInf.end());
//            m1.insert(make_pair(bID1,Tree[rank[ID1]].disInf[i]));
            m1.insert(make_pair(bID1,Tree[rank[ID1]].disInf[bID1]));
        }
        for(int j=0;j<B2.size()-1;j++){
            bID2=B2[j];
            assert(Tree[rank[ID2]].disInf.find(bID2)!=Tree[rank[ID2]].disInf.end());
//            m2.insert(make_pair(bID2,Tree[rank[ID2]].disInf[j]));
            m2.insert(make_pair(bID2,Tree[rank[ID2]].disInf[bID2]));
        }

        for(int k=0;k<B1.size()-1;k++){
            bID1=B1[k];

            if(m1[bID1]>d)
                continue;

            for(int z=0;z<B2.size()-1;z++){
                bID2=B2[z];

                if(m2[bID2]>d)
                    continue;

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


void Graph::indexsizeCTH2H(){
	long long m=0,m1=0,m2=0,m3=0,m4=0,m5=0;

	//core index
	for(int k=0;k<Label.size();k++){
		m1+=Label[k].size()*2*sizeof(int);
	}


    for(int i=0;i<PruningPointList.size();i++){
        for(auto it=PruningPointList[i].begin();it!=PruningPointList[i].end();it++){
            m2+=(1+(*it).second.size())*sizeof(int);
        }
    }

	//Periphery index
    for(int i=0;i<Tree.size();i++){
        m3+=Tree[i].dis.size()*sizeof(int);
        m3+=Tree[i].disInf.size()*sizeof(int);
    }

    for(int i=0;i< SCconNodesMT.size();i++){
        for(auto it=SCconNodesMT[i].begin(); it!=SCconNodesMT[i].end(); it++){
            m4+=sizeof(int)+(*it).second.size()*sizeof(int);
        }
    }

    //extended index
    if(ifExtension){
        assert(!IndexExt.empty());
        for(int i=0;i<nodenum;++i){
            if(!IndexExt[i].empty()){
                m5+=sizeof(int)+IndexExt[i].size()*2*sizeof(int);
            }
        }
        cout<<"Extended label size: "<<(double)m5/1024/1024<<" MB"<<endl;
    }


	//cout<<"Index size "<<(double)m1/1024/1024<<", Pruning point size "<<(double)m2/1024/1024<<endl;
	m=m1+m2+m3+m4+m5;
    cout<<"Pruning point size: "<<(double)m2/1024/1024<<" MB"<<endl;
    cout<<"H2H update index size: "<<(double)m4/1024/1024<<" MB"<<endl;
	cout<<"Minimum index size "<<(double)(m1+m3)/1024/1024<<" MB"<<endl;
    cout<<"Overall index size "<<(double)m/1024/1024<<" MB"<<endl;
}

void Graph::Decrease(int a, int b, int oldW, int newW){
    for(int i=0;i<Neighbor[a].size();i++){
        if(Neighbor[a][i].first==b){
            Neighbor[a][i].second=newW;
            break;
        }
    }
    for(int i=0;i<Neighbor[b].size();i++){
        if(Neighbor[b][i].first==a){
            Neighbor[b][i].second=newW;
            break;
        }
    }

    int pid = -1;
    extUpdate = false;

    AdjaCoreMapOld = AdjaCoreMap;

    if((CoreTag[a]==-1 && BoundTag[a].first==0) || (CoreTag[b]==-1 && BoundTag[b].first==0)){// either endpoint is non-boundary core vertex
        //cout<<"******************change 1*******************"<<endl;
//        cout<<"Core-Core"<<endl;
//		DecreasePLL(a,b,oldW,newW,AdjaCore,Label);
        DecreasePSL(a,b,oldW,newW,AdjaCore,Label);
        extUpdate = true;
    }else if(CoreTag[a]!=-1 || CoreTag[b]!=-1){//either endpoint is in-partition vertex
        //cout<<"******************change 2*******************"<<endl;
        if(CoreTag[a]!=-1)//if a is periphery vertex
            pid=CoreTag[a];
        else
            pid=CoreTag[b];

        //cout<<"<<<<<<<<<<"<<endl;
        //cout<<CoreTag[a]<<" "<<BoundTag[a]<<" "<<CoreTag[b]<<" "<<BoundTag[b]<<" pid "<<pid<<endl;
        //cout<<"decrease "<<a<<" "<<b<<" "<<oldW<<" "<<newW<<endl;
        DecreaseH2HNew(a,b,newW,Neighbor,Tree,rank,heightMax,true);
//        DecreaseH2HNew(a,b,newW,Neighbor,Tree,rank,heightMax,false);

        //cout<<">>>>>>>>>>"<<endl;

        /// check whether the update of periphery has affected core
        int ID1,ID2,dis,olddis;
        for(int i=0;i<BoundVertex[pid].size()-1;i++){
            ID1=BoundVertex[pid][i];

            //if(PartiVertexInverted[pid].find(ID1)!=PartiVertexInverted[pid].end()){//only if ID1 has neighbors in partition pid
            for(int j=i+1;j<BoundVertex[pid].size()-1;j++){
                ID2=BoundVertex[pid][j];

                //if(PartiVertexInverted[pid].find(ID2)!=PartiVertexInverted[pid].end()){
//                dis = QueryPeripheryTree(ID1,ID2,pid);///
                dis = AdjaCoreMap[ID1][ID2];
//                if(dis != temp_dis)
//                    cout<<dis<<" "<<temp_dis<<endl;
                olddis=AdjaCoreMapOld[ID1][ID2];
                if(olddis>dis){
//                   DecreasePLL(ID1,ID2,olddis,dis,AdjaCore,Label);
                    DecreasePSL(ID1,ID2,olddis,dis,AdjaCore,Label);
                    extUpdate = true;
                }


            }

        }
    }else if(BoundTag[a].first==1 && BoundTag[b].first==1){//Both end points are interface vertex
        //cout<<"******************change 3*******************"<<endl;
        ///periphery update
//        cout<<"Both interface vertex!"<<endl;
        int olddis=AdjaCoreMap[a][b];
        if(olddis>newW){
            AdjaCoreMap[a][b]=newW;
            AdjaCoreMap[b][a]=newW;
//            DecreasePLL(a,b,olddis,newW,AdjaCore,Label);
            DecreasePSL(a,b,olddis,newW,AdjaCore,Label);
            extUpdate = true;
        }

    }
    if(ifExtension){
        /// extension index update
        vector<int> partiForUpdate;
        if(extUpdate){
            for(int i=0;i<partiNum;++i){
                partiForUpdate.push_back(i);
                if(PartiUpdateExt[i]){
                    PartiUpdateExt[i]=false;
                }
            }
            ExtensionIndexUpdate(ifParallel,false,partiForUpdate);
        }else{
            for(int i=0;i<partiNum;++i){
                if(PartiUpdateExt[i]){
                    partiForUpdate.push_back(i);
                    PartiUpdateExt[i]=false;
                }
            }
            ExtensionIndexUpdate(ifParallel,false,partiForUpdate);
        }
    }

}
//New version
void Graph::Increase(int a, int b, int oldW, int newW){
    for(int i=0;i<Neighbor[a].size();i++){
        if(Neighbor[a][i].first==b){
            if(oldW != Neighbor[a][i].second){
                cout<<"Old edge weight is incorrect! "<<a<<" "<<b<<": "<<oldW<<" "<<Neighbor[a][i].second<<endl;
                oldW = Neighbor[a][i].second;
            }
            Neighbor[a][i].second=newW;
            break;
        }
    }
    for(int i=0;i<Neighbor[b].size();i++){
        if(Neighbor[b][i].first==a){
            if(oldW != Neighbor[b][i].second){
                cout<<"Old edge weight is incorrect! "<<b<<" "<<a<<": "<<oldW<<" "<<Neighbor[b][i].second<<endl;
                oldW = Neighbor[b][i].second;
            }
            Neighbor[b][i].second=newW;
            break;
        }
    }

    AdjaCoreMapOld = AdjaCoreMap;

    extUpdate = false;
    int pid = -1;
	if((CoreTag[a]==-1 && !BoundTag[a].first) || (CoreTag[b]==-1 && !BoundTag[b].first)){//edges in the core
		//cout<<"edge in core"<<endl;
//		IncreasePLL(a,b,oldW,newW,AdjaCore,Label,PruningPointNew,NoSupportedPair);
//        IncreasePSL(a,b,oldW,newW,AdjaCore,Label,PruningPointNew,NoSupportedPair);
//        IncreasePSL(a,b,oldW,newW,AdjaCore,Label,PruningPointList);
        IncreasePSL2(a,b,oldW,newW,AdjaCore,Label,PruningPointList);
        extUpdate = true;
	}else if(CoreTag[a]!=-1 || CoreTag[b]!=-1){//edges in one partition
		//cout<<"edge in partition"<<endl;

		if(CoreTag[a]!=-1)
			pid=CoreTag[a];
		else
			pid=CoreTag[b];

        int ID1,ID2,dis,olddis;
//        ID1=142502;ID2=143695;olddis=45551;dis=45703;
//        cout<<ID1<<" "<<ID2<<" "<<olddis<<" "<<dis<<endl;
//        IncreasePSL(ID1,ID2,olddis,dis,AdjaCore,Label,PruningPointNew,NoSupportedPair);
//        QueryCoreDebug(142488,143850);
//        exit(0);


		//cout<<"<<<<<<<<<<"<<endl;
		//cout<<CoreTag[a]<<" "<<BoundTag[a]<<" "<<CoreTag[b]<<" "<<BoundTag[b]<<" pid "<<pid<<endl;
        IncreaseH2HNew(a,b,oldW,newW,Neighbor,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,true);
//        IncreaseH2HNew(a,b,oldW,newW,Neighbor,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,false);
		//cout<<">>>>>>>>>>"<<endl;

        ofstream OF("/Users/zhouxj/Documents/1-Research/Datasets/NY/NYCore4.update");
        for(int i=0;i<BoundVertex[pid].size()-1;i++){
            ID1=BoundVertex[pid][i];

            //if(PartiVertexInverted[pid].find(ID1)!=PartiVertexInverted[pid].end()){//only if ID1 has neighbors in partition pid
            for(int j=i+1;j<BoundVertex[pid].size()-1;j++){
                ID2=BoundVertex[pid][j];

                dis = AdjaCoreMap[ID1][ID2];
                olddis=AdjaCoreMapOld[ID1][ID2];
                if(olddis < dis){
//                    cout<<ID1<<" "<<ID2<<" "<<olddis<<" "<<dis<<endl;
//                    if(ID1==142502 && ID2==143695){
//                        cout<<ID1<<" "<<ID2<<" "<<olddis<<" "<<dis<<endl;
//                        QueryCoreDebug(142488,143850);
//                    }

                    OF<<ID1<<" "<<ID2<<" "<<olddis<<" "<<dis<<endl;
                    IncreasePSL(ID1,ID2,olddis,dis,AdjaCore,Label,PruningPointNew,NoSupportedPair);
//                    IncreasePSL(ID1,ID2,olddis,dis,AdjaCore,Label,PruningPointList);
//                    IncreasePSL2(ID1,ID2,olddis,dis,AdjaCore,Label,PruningPointList);//list version with NoSupportedPair
//                    IncreasePLL(ID1,ID2,olddis,dis,AdjaCore,Label,PruningPoint);
                    extUpdate = true;
//                    if(ID1==142502 && ID2==143695){
////                        cout<<ID1<<" "<<ID2<<" "<<olddis<<" "<<dis<<endl;
//                        QueryCoreDebug(142488,143850);
//                    }

                }
            }
        }
        OF.close();

	}else{//Both end points are boundary vertex
		//cout<<CoreTag[a]<<" "<<BoundTag[a]<<" "<<CoreTag[b]<<" "<<BoundTag[b]<<endl;
		//cout<<"edge between boundary vertex"<<endl;
//        cout<<"Both interface vertex!"<<endl;
        int olddis=AdjaCoreMap[a][b];
        assert(olddis <= oldW);
        if(olddis==oldW){//it indicates the update of e(a,b) may influence AdjaCoreMap[a][b]
            int newDis=INF;
            set<int> newSets;
            for(auto it=SuppPartiID[a][b].begin();it!=SuppPartiID[a][b].end();++it){
                if(it->second <= newW){
                    if(it->second < newDis){
                        newDis = it->second;
                        newSets.clear();
                        newSets.insert(it->first);
                    }else if(it->second == newDis){
                        newSets.insert(it->first);
                    }
                }
            }
            if(newDis == INF){//no supportive vertex can obtain super edge weight lower than newW
                AdjaCoreMap[a][b] = newW;
                AdjaCoreMap[b][a] = newW;
                SuppPartiIDReal[a][b].first = newW;
                SuppPartiIDReal[b][a].first = newW;
                SuppPartiIDReal[a][b].second.clear();
                SuppPartiIDReal[b][a].second.clear();
//                IncreasePSL(a,b,oldW,newW,AdjaCore,Label,PruningPointNew,NoSupportedPair);
//                IncreasePSL(a,b,oldW,newW,AdjaCore,Label,PruningPointList);
                IncreasePSL2(a,b,oldW,newW,AdjaCore,Label,PruningPointList);
                extUpdate = true;
            }else{
                AdjaCoreMap[a][b] = newDis;
                AdjaCoreMap[b][a] = newDis;
                SuppPartiIDReal[a][b].first = newDis;
                SuppPartiIDReal[b][a].first = newDis;
                SuppPartiIDReal[a][b].second = newSets;
                SuppPartiIDReal[b][a].second = newSets;
//                IncreasePSL(a,b,oldW,newDis,AdjaCore,Label,PruningPointNew,NoSupportedPair);
//                IncreasePSL(a,b,oldW,newDis,AdjaCore,Label,PruningPointList);
                IncreasePSL2(a,b,oldW,newDis,AdjaCore,Label,PruningPointList);
                extUpdate = true;
            }

        }

	}

    if(ifExtension){
        /// extension index update
        vector<int> partiForUpdate;
        if(extUpdate){
//            cout<<"extUpdate: "<<extUpdate<<endl;
            for(int i=0;i<partiNum;++i){
                partiForUpdate.push_back(i);
                if(PartiUpdateExt[i]){
                    PartiUpdateExt[i]=false;
                }
            }
            ExtensionIndexUpdate(true,true,partiForUpdate);
        }else{
            for(int i=0;i<partiNum;++i){
                if(PartiUpdateExt[i]){
                    partiForUpdate.push_back(i);
                    PartiUpdateExt[i]=false;
                }
            }
            ExtensionIndexUpdate(ifParallel,true,partiForUpdate);
        }
    }

}

