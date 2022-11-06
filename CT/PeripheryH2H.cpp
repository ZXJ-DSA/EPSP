/*
 * PeripheryH2H.cpp
 *
 *  Created on: 14 Oct 2022
 *      Author: zhangmengxuan
 */
#include "head.h"

extern bool ifExtensiion;

void Graph::PartitionPostProcess(){
	AdjaCoreMap.clear();
	map<int,int> mii;
	mii.clear();
	AdjaCoreMap.assign(nodenum,mii);

	vector<pair<int,int>> vecp;
	vecp.clear();
	vector<vector<pair<int,int>>> vecvecpair;
	vecvecpair.assign(nodenum,vecp);
	AdjaParti.assign(partiNum,vecvecpair);//hug memory consumption!!
    PartiVertex.assign(partiNum,unordered_set<int>());

	map<int,set<int>> mapvec;
	mapvec.clear();
	SuppPartiID.assign(nodenum, mapvec);
	map<int,set<int>> mapset;
	mapset.clear();
	SuppPartiIDReal.assign(nodenum, mapset);

	for(int NodeID=0;NodeID<nodenum;NodeID++){
		if(CoreTag[NodeID]==-1){//core vertex + boundary vertex
			if(BoundTag[NodeID]==1){//boundary vertex

				for(int nei=0;nei<Neighbor[NodeID].size();nei++){//for each neighbor
					int neiID=Neighbor[NodeID][nei].first;

					if(CoreTag[neiID]==-1){//**********neighbor is core vertex or boundary vertex
						AdjaCoreMap[NodeID].insert(make_pair(Neighbor[NodeID][nei].first,Neighbor[NodeID][nei].second));

						if(BoundTag[neiID]==1){//if the neighbor is boundary vertex
							SuppPartiIDReal[NodeID][neiID].insert(partiNum);
							SuppPartiIDReal[neiID][NodeID].insert(partiNum);
							SuppPartiID[NodeID][neiID].insert(partiNum);
							SuppPartiID[neiID][NodeID].insert(partiNum);
						}
					}else{//**********neighbor is partition vertex
						int pid=CoreTag[neiID];
						AdjaParti[pid][NodeID].push_back(Neighbor[NodeID][nei]);
					}
				}

			}else{//core vertex
				//CoreVertexInverted.insert(make_pair(NodeID, AdjaCoreMap.size()));
				map<int,int> m1;
				m1.clear();
				for(int k=0;k<Neighbor[NodeID].size();k++)
					AdjaCoreMap[NodeID].insert(make_pair(Neighbor[NodeID][k].first,Neighbor[NodeID][k].second));
				//AdjaCoreMap.push_back(m1);
			}
		}else{//partition vertex
			int pid=CoreTag[NodeID];
            PartiVertex[pid].insert(NodeID);
			for(int k=0;k<Neighbor[NodeID].size();k++)
				AdjaParti[pid][NodeID].push_back(Neighbor[NodeID][k]);

			//PartiVertexInverted[pid].insert(make_pair(NodeID,AdjaParti[pid].size()));
			//AdjaParti[pid].push_back(Neighbor[NodeID]);
		}
	}

	/*vector<int> PeripherySize;
	PeripherySize.clear();
	for(int pid=0;pid<partiNum;pid++){
		int psize=0;
		for(int id=0;id<nodenum;id++){
			if(AdjaParti[pid][id].size()>0)
				psize+=1;
		}
		PeripherySize.push_back(psize);
		//cout<<pid<<" "<<psize<<endl;
		if(psize==0)
			cout<<"Wrong: no vertex in periphery "<<pid<<endl;
	}*/

//	int totalperiNum=0;
//	for(int i=0;i<PeripherySize.size();i++)
//		totalperiNum+=PeripherySize[i];
//	cout<<"totalperiNum "<<totalperiNum<<endl;

}

void Graph::OverlayGraph(){
	//core's adjacent complete through B*B
	for(int pid=0;pid<partiNum;pid++){
		int ID1,ID2,dis;
		for(int i=0;i<BoundVertex[pid].size();i++){
			ID1=BoundVertex[pid][i];
			if(AdjaParti[pid][ID1].size()!=0){
			//if(PartiVertexInverted[pid].find(ID1)!=PartiVertexInverted[pid].end()){//only if ID1 has neighbors in partition pid
				for(int j=i+1;j<BoundVertex[pid].size();j++){
					ID2=BoundVertex[pid][j];
					if(AdjaParti[pid][ID2].size()!=0){
					//if(PartiVertexInverted[pid].find(ID2)!=PartiVertexInverted[pid].end()){
						dis=QueryH2HPartition(ID1,ID2,pid);

						if(dis<INF){
							SuppPartiID[ID1][ID2].insert(pid);
							SuppPartiID[ID2][ID1].insert(pid);


							if(AdjaCoreMap[ID1].find(ID2)!=AdjaCoreMap[ID1].end()){
								if(AdjaCoreMap[ID1][ID2]>dis){
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
								AdjaCoreMap[ID1].insert(make_pair(ID2,dis));
								AdjaCoreMap[ID2].insert(make_pair(ID1,dis));
								SuppPartiIDReal[ID1][ID2].insert(pid);
								SuppPartiIDReal[ID2][ID1].insert(pid);
							}
						}

					}
				}
			}
		}
	}

	vector<pair<int,int>> vecpair;
	vecpair.clear();
	AdjaCore.assign(nodenum,vecpair);
	for(int i=0;i<AdjaCoreMap.size();i++){
		for(map<int,int>::iterator it=AdjaCoreMap[i].begin();it!=AdjaCoreMap[i].end();it++){
			AdjaCore[i].push_back(make_pair((*it).first, (*it).second));
		}
	}
	cout<<"***************Finish Core's graph construction***************"<<endl;
	//CorrectnessCheckCore();
}

//Query processing
int Graph::Query(int ID1, int ID2){
    int dis=INF;

    if(CoreTag[ID1]==-1 && CoreTag[ID2]==-1){//both in core
        //cout<<"Core-Core"<<endl;
        dis=QueryCore(ID1, ID2);
    }else if(CoreTag[ID1]==-1 && CoreTag[ID2]!=-1){//ID2 in partition, ID1 in core
//        cout<<"Core-Parti"<<endl;
        if(ifExtensiion){
            dis=QueryPartiCoreExt(ID2, ID1);
        }else{
            dis=QueryPartiCore(ID2, ID1);
        }
    }else if(CoreTag[ID1]!=-1 && CoreTag[ID2]==-1){//ID1 in partition, ID2 in core
//        cout<<"Parti-Core"<<endl;
        if(ifExtensiion){
            dis=QueryPartiCoreExt(ID1, ID2);
        }else{
            dis=QueryPartiCore(ID1, ID2);
        }
    }else if(CoreTag[ID1]!=-1 && CoreTag[ID2]!=-1){//both in partition
//        cout<<"Parti-Parti"<<endl;
        if(ifExtensiion){
            dis=QueryPartiPartiExt(ID1,ID2);
        }else{
            dis=QueryPartiParti(ID1,ID2);
        }
    }
    return dis;
}

int Graph::QueryPartiCore(int ID1, int ID2){
	int d=INF;

	int pid=CoreTag[ID1];
	int bid;
	int dis1,dis2;
	for(int k=0;k<BoundVertex[pid].size();k++){
		bid=BoundVertex[pid][k];//For each boundary vertex
		dis1=QueryH2HPartition(ID1,bid,pid);
		dis2=QueryCore(bid,ID2);
		if(d>dis1+dis2)
			d=dis1+dis2;
	}

	return d;
}

void Graph::ExtensionIndex(pair<int,int> pidRange){
    int ID1, ID2;
    int b1, b2;
    int temp_dis=INF;

    for(int pid=pidRange.first;pid<pidRange.second;++pid){
        vector<int> B=BoundVertex[pid];
        for(auto it=PartiVertex[pid].begin();it!=PartiVertex[pid].end();++it){
            ID1 = *it;
            if(BoundVertexSet[pid].find(ID1) != BoundVertexSet[pid].end()){//if found, i.e., it is boundary vertex
                continue;
            }
            // for each non-boundary vertex
            for(int i=0;i<B.size();++i){
                b1 = B[i];
                for(auto it2 = IndexCore[b1].begin();it2!= IndexCore[b1].end();++it2){//for each label of b1
                    b2 = it2->first;
                    if(IndexExt[ID1].find(b2) == IndexExt[ID1].end()){//if not found
                        IndexExt[ID1][b2] = QueryH2HPartition(ID1,b1,pid) + it2->second;
                    }else{
                        temp_dis = QueryH2HPartition(ID1,b1,pid) + it2->second;
                        if(IndexExt[ID1][b2] > temp_dis){
                            IndexExt[ID1][b2] = temp_dis;
                        }
                    }
                }
                // for the root vertex
                if(IndexExt[ID1].find(b1) == IndexExt[ID1].end()){//if not found
                    IndexExt[ID1][b1] = QueryH2HPartition(ID1,b1,pid);
                }else{
                    temp_dis = QueryH2HPartition(ID1,b1,pid);
                    if(IndexExt[ID1][b1] > temp_dis){
                        IndexExt[ID1][b1] = temp_dis;
                    }
                }
            }
        }
    }

}
//function of constructing the extension index for periphery vertex
void Graph::ExtensionIndexConstruct(bool ifParallel){
    int ID1, ID2;
    int b1, b2;
    int temp_dis=INF;
    IndexExt.assign(nodenum,unordered_map<int,int>());
    if(ifParallel){
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
                threadf.add_thread(new boost::thread(&Graph::ExtensionIndex, this, p));
            }
            threadf.join_all();
        }else{
            boost::thread_group threadf;
            for(int pid=0;pid<partiNum;++pid) {
                threadf.add_thread(new boost::thread(&Graph::ExtensionIndex, this, make_pair(pid,pid+1)));
            }
            threadf.join_all();
        }
    }
    else{
        for(int pid=0;pid<partiNum;++pid){//for each partition
            vector<int> B=BoundVertex[pid];
            for(auto it=PartiVertex[pid].begin();it!=PartiVertex[pid].end();++it){
                ID1 = *it;
                if(BoundVertexSet[pid].find(ID1) != BoundVertexSet[pid].end()){//if found, i.e., it is boundary vertex
                    continue;
                }
                // for each non-boundary vertex
                for(int i=0;i<B.size();++i){
                    b1 = B[i];
                    for(auto it2 = IndexCore[b1].begin();it2!= IndexCore[b1].end();++it2){//for each label of b1
                        b2 = it2->first;
                        if(IndexExt[ID1].find(b2) == IndexExt[ID1].end()){//if not found
                            IndexExt[ID1][b2] = QueryH2HPartition(ID1,b1,pid) + it2->second;
                        }else{
                            temp_dis = QueryH2HPartition(ID1,b1,pid) + it2->second;
                            if(IndexExt[ID1][b2] > temp_dis){
                                IndexExt[ID1][b2] = temp_dis;
                            }
                        }
                    }
                    // for the root vertex
                    if(IndexExt[ID1].find(b1) == IndexExt[ID1].end()){//if not found
                        IndexExt[ID1][b1] = QueryH2HPartition(ID1,b1,pid);
                    }else{
                        temp_dis = QueryH2HPartition(ID1,b1,pid);
                        if(IndexExt[ID1][b1] > temp_dis){
                            IndexExt[ID1][b1] = temp_dis;
                        }
                    }
                }
            }
        }
    }

}

int Graph::QueryPartiCoreExt(int ID1, int ID2){//ID1 is in partition, ID2 is in core
    int d=INF;

    bool flag_union=false;

    unordered_map<int,int>::iterator it;
    int hub, dis1, dis2;
    int hubfinal,dis1final,dis2final;
//    cout << "ID1 label Size:" << IndexExt[ID1].size() << "\tID2 label Size:" <<IndexCore[ID2].size()<<endl;
    for(it=IndexExt[ID1].begin();it!=IndexExt[ID1].end();it++){
        hub=(*it).first;
        dis1=(*it).second;
        if(IndexCore[ID2].find(hub)!=IndexCore[ID2].end()){
            flag_union = true;
            dis2=IndexCore[ID2][hub];
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
        cout<<"There is no intersection of these two vertex label!"<<endl;
    }
    if(d == INF){
        cout<<"Wrong!!! "<<ID1<<" "<<ID2<<" "<<INF<<endl;
    }
    return d;
}

int Graph::QueryPartiPartiExt(int ID1, int ID2){
    int d=INF;

    int pid1=CoreTag[ID1];
    int pid2=CoreTag[ID2];
    if(pid1==pid2){//if in the same partition
        int temp_dis = QueryH2HPartition(ID1,ID2,pid1);//d2
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

int Graph::QueryPartiParti(int ID1, int ID2){
	int d=INF;

	int pid1=CoreTag[ID1];
	int pid2=CoreTag[ID2];
	if(pid1==pid2){//if in the same partition
        int temp_dis = QueryH2HPartition(ID1,ID2,pid1);
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
	}else{//if in different partitions
		vector<int> B1=BoundVertex[pid1];
		vector<int> B2=BoundVertex[pid2];

		map<int,int> m1,m2;
		m1.clear();
		m2.clear();
		int bID1, bID2, tempdis;
		for(int i=0;i<B1.size();i++){
			bID1=B1[i];
			m1.insert(make_pair(bID1,QueryH2HPartition(ID1,bID1,pid1)));
		}
		for(int j=0;j<B2.size();j++){
			bID2=B2[j];
			m2.insert(make_pair(bID2,QueryH2HPartition(ID2,bID2,pid2)));
		}

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

	return d;
}

void Graph::indexsizeCTH2H(){
	long long m=0,m1=0,m2=0,m3=0;

	//core index
	for(int k=0;k<IndexCore.size();k++){
		m1+=IndexCore[k].size()*2*sizeof(int);
	}

	for(int i=0;i<PruningPointCore.size();i++){
		for(unordered_map<int,vector<int>>::iterator it=PruningPointCore[i].begin();it!=PruningPointCore[i].end();it++){
			m2+=(1+(*it).second.size())*sizeof(int);
		}
	}

	//Periphery index
	for(int pid=0;pid<partiNum;pid++){
		for(int i=0;i<Trees[pid].size();i++){
			m3+=Trees[pid][i].dis.size()*sizeof(int);
		}
		for(int i=0;i< SCconNodesMTs[pid].size();i++){
			for(auto it=SCconNodesMTs[pid][i].begin(); it!=SCconNodesMTs[pid][i].end(); it++){
				m3+=sizeof(int)+(*it).second.size()*sizeof(int);
			}
		}
	}

	//cout<<"Index size "<<(double)m1/1024/1024<<", Pruning point size "<<(double)m2/1024/1024<<endl;
	m=m1+m2+m3;
	cout<<"Index size "<<(double)m/1024/1024<<" MB"<<endl;
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
	if((CoreTag[a]==-1 && BoundTag[a]==0) || (CoreTag[b]==-1 && BoundTag[b]==0)){//edges in the core
		//cout<<"******************change 1*******************"<<endl;
//		DecreasePLL(a,b,oldW,newW,AdjaCore,IndexCore);
        DecreasePSL(a,b,oldW,newW,AdjaCore,IndexCore);
	}else if(CoreTag[a]!=-1 || CoreTag[b]!=-1){//edges in one partition
		//cout<<"******************change 2*******************"<<endl;
		if(CoreTag[a]!=-1)
			pid=CoreTag[a];
		else
			pid=CoreTag[b];

		//cout<<"<<<<<<<<<<"<<endl;
		//cout<<CoreTag[a]<<" "<<BoundTag[a]<<" "<<CoreTag[b]<<" "<<BoundTag[b]<<" pid "<<pid<<endl;
		//cout<<"decrease "<<a<<" "<<b<<" "<<oldW<<" "<<newW<<endl;
		DecreaseH2H(a,b,newW,AdjaParti[pid],Trees[pid],ranks[pid],heightMaxs[pid]);
		//cout<<">>>>>>>>>>"<<endl;

		/*//change the partition graph
		for(int k=0;k<AdjaParti[pid][a].size();k++){
			if(AdjaParti[pid][a][k].first==b){
				AdjaParti[pid][a][k].second=newW;
				break;
			}
		}
		for(int k=0;k<AdjaParti[pid][b].size();k++){
			if(AdjaParti[pid][b][k].first==a){
				AdjaParti[pid][b][k].second=newW;
				break;
			}
		}*/

		int ID1,ID2,dis,olddis;
		for(int i=0;i<BoundVertex[pid].size();i++){
			ID1=BoundVertex[pid][i];
			if(AdjaParti[pid][ID1].size()!=0){
			//if(PartiVertexInverted[pid].find(ID1)!=PartiVertexInverted[pid].end()){//only if ID1 has neighbors in partition pid
				for(int j=i+1;j<BoundVertex[pid].size();j++){
					ID2=BoundVertex[pid][j];
					if(AdjaParti[pid][ID2].size()!=0){
					//if(PartiVertexInverted[pid].find(ID2)!=PartiVertexInverted[pid].end()){
						dis=QueryH2HPartition(ID1,ID2,pid);
						olddis=AdjaCoreMap[ID1][ID2];
						if(olddis>dis){
//                            DecreasePLL(ID1,ID2,olddis,dis,AdjaCore,IndexCore);
                            DecreasePSL(ID1,ID2,olddis,dis,AdjaCore,IndexCore);
                        }

					}
				}
			}
		}
	}else if(BoundTag[a]==1 && BoundTag[b]==1){//Both end points are boundary vertex
		//cout<<"******************change 3*******************"<<endl;
		int olddis=AdjaCoreMap[a][b];
		if(olddis>newW){
//            DecreasePLL(a,b,olddis,newW,AdjaCore,IndexCore);
            DecreasePSL(a,b,olddis,newW,AdjaCore,IndexCore);
        }

	}

    if(ifExtensiion){
        /// extension index update
        if(extUpdate){
            ExtensionIndexConstruct(ifParallel);
        }else if(pid != -1){
            ExtensionIndex(make_pair(pid,pid+1));
        }
    }

}

void Graph::Increase(int a, int b, int oldW, int newW){
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
    extUpdate = false;
    int pid = -1;
	if((CoreTag[a]==-1 && BoundTag[a]==0) || (CoreTag[b]==-1 && BoundTag[b]==0)){//edges in the core
		//cout<<"edge in core"<<endl;
//		IncreasePLL(a,b,oldW,newW,AdjaCore,IndexCore,PruningPointCore,NoSupportedPairCore);
        IncreasePSL(a,b,oldW,newW,AdjaCore,IndexCore,PruningPointCore,NoSupportedPairCore);
	}else if(CoreTag[a]!=-1 || CoreTag[b]!=-1){//edges in one partition
		//cout<<"edge in partition"<<endl;

		if(CoreTag[a]!=-1)
			pid=CoreTag[a];
		else
			pid=CoreTag[b];

		//cout<<"<<<<<<<<<<"<<endl;
		//cout<<CoreTag[a]<<" "<<BoundTag[a]<<" "<<CoreTag[b]<<" "<<BoundTag[b]<<" pid "<<pid<<endl;
		IncreaseH2H(a,b,oldW,newW,AdjaParti[pid],Trees[pid],ranks[pid],heightMaxs[pid],SCconNodesMTs[pid],VidtoTNids[pid]);
		//cout<<">>>>>>>>>>"<<endl;

		/*//change the partition graph
		for(int k=0;k<AdjaParti[pid][a].size();k++){
			if(AdjaParti[pid][a][k].first==b){
				AdjaParti[pid][a][k].second=newW;
				break;
			}
		}
		for(int k=0;k<AdjaParti[pid][b].size();k++){
			if(AdjaParti[pid][b][k].first==a){
				AdjaParti[pid][b][k].second=newW;
				break;
			}
		}*/

		int ID1,ID2,dis,olddis;
		for(int i=0;i<BoundVertex[pid].size();i++){
			ID1=BoundVertex[pid][i];
			if(AdjaParti[pid][ID1].size()!=0){
			//if(PartiVertexInverted[pid].find(ID1)!=PartiVertexInverted[pid].end()){//only if ID1 has neighbors in partition pid
				for(int j=i+1;j<BoundVertex[pid].size();j++){
					ID2=BoundVertex[pid][j];
					if(AdjaParti[pid][ID2].size()!=0){
					//if(PartiVertexInverted[pid].find(ID2)!=PartiVertexInverted[pid].end()){

						if(SuppPartiIDReal[ID1][ID2].find(pid)!=SuppPartiIDReal[ID1][ID2].end()){
							int olddis=AdjaCoreMap[ID1][ID2];
							int newdis=QueryH2HPartition(ID1,ID2,pid);
							if(newdis>olddis){
								if(SuppPartiIDReal[ID1][ID2].size()==1){//refine the real supported partition
									SuppPartiIDReal[ID1][ID2].clear();
									int p, partidis;
									int finaldis=newdis;
									for(set<int>::iterator its=SuppPartiID[a][b].begin();its!=SuppPartiID[a][b].end();its++){
										p=*its;
										if(p==partiNum){
											for(int n1=0;n1<Neighbor[ID1].size();n1++){
												if(Neighbor[ID1][n1].first==ID2){
													partidis=Neighbor[ID1][n1].second;
													break;
												}
											}
										}else{
											partidis=QueryH2HPartition(ID1,ID2,p);
										}

										if(partidis<finaldis){
											SuppPartiIDReal[ID1][ID2].clear();
											SuppPartiIDReal[ID1][ID2].insert(p);
											finaldis=partidis;
										}else if(partidis==finaldis){
											SuppPartiIDReal[ID1][ID2].insert(p);
										}
									}
									//synchronize the supp information
									SuppPartiIDReal[ID2][ID1]=SuppPartiIDReal[ID1][ID2];

//									IncreasePLL(ID1,ID2,olddis,finaldis,AdjaCore,IndexCore,PruningPointCore,NoSupportedPairCore);
                                    IncreasePSL(ID1,ID2,olddis,finaldis,AdjaCore,IndexCore,PruningPointCore,NoSupportedPairCore);
								}else if(SuppPartiIDReal[ID1][ID2].size()>1){//delete the current supported partition
									SuppPartiIDReal[ID1][ID2].erase(pid);
									SuppPartiIDReal[ID2][ID1].erase(pid);
								}
							}
						}

					}
				}
			}
		}


	}else{//Both end points are boundary vertex
		//cout<<CoreTag[a]<<" "<<BoundTag[a]<<" "<<CoreTag[b]<<" "<<BoundTag[b]<<endl;
		//cout<<"edge between boundary vertex"<<endl;

		//for(set<int>::iterator its=SuppPartiID[a][b].begin();its!=SuppPartiID[a][b].end();its++)
			//cout<<*its<<endl;

		//cout<<"sup number "<<SuppPartiID[a][b].size()<<endl;

		if(SuppPartiIDReal[a][b].find(partiNum)!=SuppPartiIDReal[a][b].end()){
			if(SuppPartiIDReal[a][b].size()==1){//refine the real supported partition
				SuppPartiIDReal[a][b].clear();
				int p, partidis;
				int finaldis=newW;
				for(set<int>::iterator its=SuppPartiID[a][b].begin();its!=SuppPartiID[a][b].end();its++){
					p=*its;
					if(p==partiNum){
						partidis=newW;
					}else{
						partidis=QueryH2HPartition(a,b,p);
					}
					//cout<<p<<" "<<partidis<<endl;
					if(partidis<finaldis){
						SuppPartiIDReal[a][b].clear();
						SuppPartiIDReal[a][b].insert(p);
						finaldis=partidis;
					}else if(partidis==finaldis){
						SuppPartiIDReal[a][b].insert(p);
					}
				}
				SuppPartiIDReal[b][a]=SuppPartiIDReal[a][b];
//				IncreasePLL(a,b,oldW,finaldis,AdjaCore,IndexCore,PruningPointCore,NoSupportedPairCore);
                IncreasePSL(a,b,oldW,finaldis,AdjaCore,IndexCore,PruningPointCore,NoSupportedPairCore);
			}else if(SuppPartiIDReal[a][b].size()>1){//delete the current supported partition
				SuppPartiIDReal[a][b].erase(partiNum);
				SuppPartiIDReal[b][a].erase(partiNum);
			}
		}
	}

    if(ifExtensiion){
        /// extension index update
        if(extUpdate){
            ExtensionIndexConstruct(ifParallel);
        }else if(pid != -1){
            ExtensionIndex(make_pair(pid,pid+1));
        }
    }

}
