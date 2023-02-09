/*
 * PSL.cpp
 *
 *  Created on: 13 Oct 2022
 *      Author: zhangmengxuan
 */
#include "head3.h"


void Graph::IndexConstructMThread2New(){
	unordered_map<int,int> map0; map0.clear();
	Label.assign(nodenum, map0);
	Dhop.assign(nodenum, map0);
//	PruningPointNew.clear();
//	unordered_map<int,vector<int>> map1;
//	PruningPointNew.assign(nodenum,map1);

    PruningPoint.clear();
    unordered_map<int,set<int>> map1; map1.clear();
    PruningPoint.assign(nodenum,map1);

//	PruningPointStepNew.clear();
//	unordered_map<int,unordered_map<int,int>> map2;
//	PruningPointStepNew.assign(nodenum, map2);
//	NoSupportedPair.clear();

	DvertexNew.assign(nodenum, true);

	for(int i=0;i<nodenum;i++){
		Label[i].insert(make_pair(i,0));
		Dhop[i].insert(make_pair(i,0));
		//Dvectex.push_back(i);
	}

	//LabelStep.push_back(Dhop);

	bool flag=true;
	int step=0;
	while(flag){
		cout<<"step "<<step<<" finish!"<<endl;
//		LabelStep.push_back(Dhop);
		flag=DhopLableRefreshMulti2New(step+1);
		step+=1;
	}
	cout<<"Index finish construction"<<endl;
}

bool Graph::DhopLableRefreshMulti2New(int step){
	bool flag=false;
	vector<unordered_map<int,int>> newDhop;
	unordered_map<int,int> m0; m0.clear();
	newDhop.assign(nodenum,m0);
	vector<int> newDvec;
	vector<vector<tri>> PP;
	vector<tri> prp; prp.clear();
	PP.assign(threadnum, prp);

	boost::thread_group thread;

	//int stepthr=Dvectex.size()/threadnum; //cout<<" Dvectex size "<<Dvectex.size()<<endl;
	vector<vector<int>> ProcessID;
	vector<int> vvv; ProcessID.assign(threadnum,vvv);
	threadDistribute(ProcessID);

	for(int i=0;i<ProcessID.size();i++){
		vector<int> p=ProcessID[i];
		thread.add_thread(new boost::thread(&Graph::labelMultiThread2New, this, boost::ref(newDhop), p,step));
	}
	thread.join_all();

	//cout<<"lllllllllllllllllllllllllllllll"<<endl;

	int zerolabel=0;
	Dhop.assign(newDhop.begin(), newDhop.end());
	Dvectex.clear(); set<int> Dset;
	DvertexNew.assign(nodenum,false);
	for(int nodeID=0;nodeID<nodenum;nodeID++){
		if(Dhop[nodeID].size()>0){
			flag=true;
			for(unordered_map<int,int>::iterator it=Dhop[nodeID].begin();it!=Dhop[nodeID].end();it++){
				if(Label[nodeID].find((*it).first)!=Dhop[nodeID].end()){//if found
					Label[nodeID][(*it).first]=(*it).second;
				}else{//if not found
					Label[nodeID].insert(*it);
				}
			}
			for(int i=0;i<AdjaCore[nodeID].size();i++){//Neighbors
				if(Dset.find(AdjaCore[nodeID][i].first)==Dset.end()){
					Dset.insert(AdjaCore[nodeID][i].first);
					//Dvectex.push_back(Neighbors[nodeID][i].first);
					DvertexNew[AdjaCore[nodeID][i].first]=true;
				}
			}
		}else if(Dhop[nodeID].size()==0)
			zerolabel+=1;
	}

	//cout<<"zero "<<zerolabel<<" Vertex to change "<<Dvectex.size()<<endl;
	return flag;
}

void Graph::threadDistribute(vector<vector<int>>& processID){
	int ID;
	int cnt=0;
	int threadOrder;
	int a;
	for(int r=0;r<nodenum;r++){
		ID=vNodeOrder[r];
		if(DvertexNew[ID]){
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

void Graph::labelMultiThread2New(vector<unordered_map<int,int>>& newDhop, vector<int>& p,int step){
//	sm->wait();

	for(int i=0;i<p.size();i++){
		int nodeID=p[i];
		unordered_map<int,int> Dhop0; Dhop0.clear();
		int neighID, neighW;
		for(int Index=0;Index<AdjaCore[nodeID].size();Index++){
			neighID=AdjaCore[nodeID][Index].first;
			//cout<<"neighID "<<neighID<<" Dhop.size "<<Dhop[neighID].size()<<endl;
			if(Dhop[neighID].size()>0){
				neighW=AdjaCore[nodeID][Index].second;

				unordered_map<int,int>::iterator it=Dhop[neighID].begin();
				int hub, dis, d;
				for(;it!=Dhop[neighID].end();it++){
					hub=(*it).first; dis=(*it).second;d=neighW+dis;
					if(NodeOrder[hub]>NodeOrder[nodeID]){//if r(hub) > r(nodeID), i.e., r(w) > r(v)
						int TempDis; vector<int> SupNode;
						ShortestDisQuery1(nodeID, hub,SupNode,TempDis); //cout<<"nodeID "<<nodeID<<" TempDis "<<TempDis<<" d "<<d<<endl;
						//ShortestDisQuery2(nodeID, hub,SupNode,TempDis,d);
						if(TempDis>d){
							if(Dhop0.find(hub)!=Dhop0.end()){
								if(Dhop0[hub]>d) Dhop0[hub]=d;
							}else{
								Dhop0.insert(make_pair(hub,d));
							}
						}
						else{//if Q(hub,nodeID,L<h) <= e(nodeID,neighID) + Lh-1(neighID)[hub], i.e., Q(w,v,L<h) <= e(v,u) + Lh-1(u)[w]
							for(int k=0;k<SupNode.size();k++){
								int supn=SupNode[k];

								if(supn!=nodeID && supn!=hub){
									vSm[nodeID]->wait();
//									PruningPointNew[nodeID][supn].push_back(hub);//((v,c),w)
//									PruningPointStepNew[nodeID][supn].insert(make_pair(hub,step));
                                    PruningPoint[nodeID][supn].insert(hub);

									vSm[nodeID]->notify();

									vSm[hub]->wait();
//									PruningPointNew[hub][supn].push_back(nodeID);//((w,c),v)
//									PruningPointStepNew[hub][supn].insert(make_pair(nodeID,step));
                                    PruningPoint[hub][supn].insert(nodeID);

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

//	sm->notify();
}

int Graph::ShortestDisQuery1(int ID1,int ID2,vector<int>& SupNode, int& d){
	d=INF;
	unordered_map<int,int>::iterator it;
	int hub, dis1, dis2;
	for(it=Label[ID1].begin();it!=Label[ID1].end();it++){
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
//old version
/*void Graph::IncreasePSL(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label,vector<unordered_map<int,vector<int>>> &PruningPointNew,set<pair<int,int>> &NoSupportedPair){
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

	int LID,HID;
	if(NodeOrder[a]>NodeOrder[b]){
		LID=b; HID=a;
	}else{
		LID=a; HID=b;
	}

	int dis,disvally,dispeak,peakhub;
	//activate or not
	dis=DisQueryLower1(LID,HID,Neighbors,Label);//d1', via a neighbor of LID (which has lower order than HID) to reach HID
	dispeak=DisQueryPeak(LID,HID,Label).first;//d2, via the common hub vertex of LID and HID
	if(dispeak<=oldW)//if d2 is lower than oldW, the increase of oldW will not affect the shortest distance between LID and HID
        return;
    if(ifDebug){
        if(dis == INF)
            cout<<DisQueryLower1(LID,HID,Neighbors,Label);
        if(dis==oldW){
            cout<<dis<<" "<<oldW<<endl;
        }
    }

    if(dis>oldW){//if d1' is larger than oldW, it indicates the shortest distance between LID and HID is equal to oldW, index update triggered
		vector<vector<pair<int,int>>> Change;//the label that is changed
		vector<pair<int,int>> vec;
		Change.assign(nodenum,vec);
		set<int> WaitPro;//the vertex waited for label propagation, i.e., AL1
		vector<vector<int>> ChangeP;
		vector<int> vecint;
		ChangeP.assign(nodenum,vecint);
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
                        if(ifDebug){
                            int tempD = DijkstraCore(LID,hubID);
                            if(tempD != Label[LID][hubID]){
//                                cout<<LID<<" "<<hubID<<": "<<Label[LID][hubID]<<" "<<tempD<<" "<<Dijkstra(LID,hubID,Neighbor)<<endl;////
                                cout<<LID<<" "<<hubID<<": "<<Label[LID][hubID]<<" "<<tempD<<endl;////
                            }
                        }

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
                        if(ifDebug){
                            int tempD = DijkstraCore(HID,hubID);
                            if(tempD != Label[HID][hubID]){
//                                cout<<HID<<" "<<hubID<<": "<<Label[HID][hubID]<<" "<<tempD<<" "<<Dijkstra(HID,hubID,Neighbor)<<endl;////
                                cout<<HID<<" "<<hubID<<": "<<Label[HID][hubID]<<" "<<tempD<<endl;////
                            }
                        }

                        //cout<<"weight "<<HID<<" "<<hubID<<" "<<disvally<<endl;
                    }
                }
			}
		}

		while(WaitProP.size()>0 || WaitPro.size()>0){
			set<int> WaitProTem;
			vector<vector<pair<int,int>>> ChangeTem;
			vector<pair<int,int>> vec;
			ChangeTem.assign(nodenum,vec);
			set<int> WaitProPTem;
			vector<int> vecint;
			vector<vector<int>> ChangePTem;
			ChangePTem.assign(nodenum,vecint);

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
                                disvally=DisQueryVally(neiID,hID,Neighbors,Label);
                                if(Label[neiID][hID]<disvally){
                                    WaitProTem.insert(neiID);
                                    ChangeTem[neiID].push_back(make_pair(hID, neiDis+hDis));
                                    Label[neiID][hID]=disvally;//may not be the final correct value at this moment

                                    //cout<<"1--1 "<<neiID<<" "<<hID<<" "<<disvally<<endl;
                                }
                            }

						}
					}

					//Change->ChangeP: AL1->AL2
					if(PruningPointNew[curID].find(hID)!=PruningPointNew[curID].end()){//check whether the pruned vertex should be inserted to curID again
						for(int snum=0;snum<PruningPointNew[curID][hID].size();snum++){
							int s=PruningPointNew[curID][hID][snum];

							if(NodeOrder[s]<NodeOrder[curID] && NoSupportedPair.find(make_pair(s,curID))==NoSupportedPair.end()){//it is not in NoSupportedPair
								disvally=DisQueryVally(s,curID,Neighbors,Label);
                                pair<int,int> peakPair;
								peakPair=DisQueryPeak(s,curID,Label);
                                dispeak=peakPair.first; peakhub=peakPair.second;
								if(dispeak>disvally){
									WaitProPTem.insert(s);
									ChangePTem[s].push_back(curID);
									Label[s][curID]=disvally;
									NoSupportedPair.insert(make_pair(s,curID));

                                    if(ifDebug){
                                        int tempD = DijkstraCore(s,curID);
                                        if(tempD != Label[s][curID]){
//                                            cout<<s<<" "<<curID<<": "<<Label[s][curID]<<" "<<tempD<<" "<<Dijkstra(s,curID,Neighbor)<<endl;//////
                                            cout<<s<<" "<<curID<<": "<<Label[s][curID]<<" "<<tempD<<endl;
                                        }
                                    }

									//cout<<"1--2 "<<s<<" "<<curID<<" "<<disvally<<endl;
								}
							}else if(NodeOrder[s]>NodeOrder[curID] && NoSupportedPair.find(make_pair(curID,s))==NoSupportedPair.end()){
								disvally=DisQueryVally(curID,s,Neighbors,Label);
                                pair<int,int> peakPair;
                                peakPair=DisQueryPeak(curID,s,Label);//the old P2 distance dominates P1 distance, but the increase of d(curID,hID) may violate the previous balance
                                dispeak=peakPair.first; peakhub=peakPair.second;
                                if(dispeak>disvally){
									WaitProPTem.insert(curID);
									ChangePTem[curID].push_back(s);
									Label[curID][s]=disvally;//should be the final correct value
									NoSupportedPair.insert(make_pair(curID,s));//insert to indicate that there is no pruning between curID and s

                                    if(ifDebug){
                                        int tempD = DijkstraCore(curID,s);
                                        if(tempD != Label[curID][s]){
//                                            cout<<curID<<" "<<s<<": "<<Label[curID][s]<<" "<<tempD<<" "<<Dijkstra(curID,s,Neighbor)<<endl;//////
                                            cout<<curID<<" "<<s<<": "<<Label[curID][s]<<" "<<tempD<<endl;
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
							disvally=DisQueryVally(neiid, u,Neighbors,Label);
							dispeak=DisQueryPeak(neiid, u,Label).first;
							if(disvally<dispeak){
								if(Label[neiid].find(u)==Label[neiid].end() || (Label[neiid].find(u)!=Label[neiid].end() && Label[neiid][u]>disvally)){//if not found or found but disvally is lower
                                    WaitProPTem.insert(neiid);
                                    ChangePTem[neiid].push_back(u);
                                    Label[neiid][u]=disvally;//should be the final correct value
                                    NoSupportedPair.insert(make_pair(neiid,u));

                                    if(ifDebug){
                                        int tempD = DijkstraCore(neiid,u);
                                        if(tempD != Label[neiid][u]){
//                                                cout<<neiid<<" "<<u<<": "<<Label[neiid][u]<<" "<<tempD<<" "<<Dijkstra(neiid,u,Neighbor)<<endl;
                                            cout<<neiid<<" "<<u<<": "<<Label[neiid][u]<<" "<<tempD<<endl;
                                        }
                                    }

                                    //cout<<"2--2 "<<neiid<<" "<<u<<" "<<disvally<<endl;
								}
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
}*/
//new version by set
void Graph::IncreasePSL(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label,vector<unordered_map<int,set<int>>> &PruningPointNew){
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

    int LID,HID;
    if(NodeOrder[a]>NodeOrder[b]){
        LID=b; HID=a;
    }else{
        LID=a; HID=b;
    }

    int dis,disvally,dispeak,peakhub;
    pair<int,int> peakPair;
    set<tuple<int,int,int>> outdatedPruning;
    //activate or not
    dis=DisQueryLower1(LID,HID,Neighbors,Label);//d1', via a neighbor of LID (which has lower order than HID) to reach HID
    dispeak=DisQueryPeak(LID,HID,Label).first;//d2, via the common hub vertex of LID and HID
    if(dispeak<=oldW)//if d2 is lower than oldW, the increase of oldW will not affect the shortest distance between LID and HID
        return;
    if(ifDebug){
        if(dis == INF)
            cout<<DisQueryLower1(LID,HID,Neighbors,Label);
        if(dis==oldW){
            cout<<dis<<" "<<oldW<<endl;
        }
    }

    if(dis>oldW){//if d1' is larger than oldW, it indicates the shortest distance between LID and HID is equal to oldW, index update triggered
        vector<vector<pair<int,int>>> Change;//the label that is changed
        vector<pair<int,int>> vec;
        Change.assign(nodenum,vec);
        set<int> WaitPro;//the vertex waited for label propagation, i.e., AL1
        vector<vector<int>> ChangeP;
        vector<int> vecint;
        ChangeP.assign(nodenum,vecint);
        set<int> WaitProP;//the vertex waited for label propagation, i.e., AL2

        //check pruning point
//        if(!PruningPointNew[HID].empty()){
//            for(auto it=PruningPointNew[HID].begin();it!=PruningPointNew[HID].end();++it){
//                int lowID = it->first;
//                if(NodeOrder[lowID] < NodeOrder[LID]){
//                    disvally=DisQueryVally(lowID,LID, Neighbors,Label);
//                    peakPair=DisQueryPeak(lowID,LID,Label);
//                    dispeak=peakPair.first; peakhub=peakPair.second;
//                    if(dispeak <= disvally){
//                        for(auto it2=it->second.begin();it2!=it->second.end();++it2){
//                            if(*it2 == LID){
//                                WaitPro.insert(lowID);
//                                Change[lowID].push_back(make_pair(LID, dispeak));
//                                if(Label[lowID].find(HID)!=Label[lowID].end()){//if found
//                                    Change[lowID].push_back(make_pair(HID, Label[lowID][HID]));
//                                }
//                                Label[lowID][LID] = dispeak;
//                            }
//                        }
//                    }
//
//                }
//
//            }
//        }


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
                        if(ifDebug){
                            int tempD = DijkstraCore(LID,hubID);
                            if(tempD != Label[LID][hubID]){
//                                cout<<LID<<" "<<hubID<<": "<<Label[LID][hubID]<<" "<<tempD<<" "<<Dijkstra(LID,hubID,Neighbor)<<endl;////
                                cout<<"Wrong for LID's label! "<<LID<<" "<<hubID<<": "<<Label[LID][hubID]<<" "<<tempD<<endl;////
                            }
                        }

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
                        if(ifDebug){
                            int tempD = DijkstraCore(HID,hubID);
                            if(tempD != Label[HID][hubID]){
//                                cout<<HID<<" "<<hubID<<": "<<Label[HID][hubID]<<" "<<tempD<<" "<<Dijkstra(HID,hubID,Neighbor)<<endl;////
                                cout<<"Wrong for HID's label! "<<HID<<" "<<hubID<<": "<<Label[HID][hubID]<<" "<<tempD<<endl;////
                            }
                        }

                        //cout<<"weight "<<HID<<" "<<hubID<<" "<<disvally<<endl;
                    }
                }
            }
        }



        while(WaitProP.size()>0 || WaitPro.size()>0){//the vertex set that corresponding labels have been changed
            set<int> WaitProTem;
            vector<vector<pair<int,int>>> ChangeTem;
            vector<pair<int,int>> vec;
            ChangeTem.assign(nodenum,vec);
            set<int> WaitProPTem;
            vector<int> vecint;
            vector<vector<int>> ChangePTem;
            ChangePTem.assign(nodenum,vecint);

            //Change->Change & ChangeP: AL1->AL1 and AL1->AL2
            for(auto it=WaitPro.begin();it!=WaitPro.end();it++){
                int curID=*it;
                vector<pair<int,int>> curChange=Change[curID];
                int neiID, neiDis, hID, hDis;

//                if(curID == 142488 || curID == 142515 || curID == 143850)
//                {
//                    cout<<curID<<endl;
//                }

                for(int j=0;j<curChange.size();j++){
                    hID=curChange[j].first; hDis=curChange[j].second;
                    //Change->Change: AL1->AL1
                    for(int k=0;k<Neighbors[curID].size();k++){
                        neiID=Neighbors[curID][k].first; neiDis=Neighbors[curID][k].second;
                        if(Label[neiID].find(hID)!=Label[neiID].end() ){//&& neiDis+hDis==Label[neiID][hID]
                            if(neiDis+hDis==Label[neiID][hID]){
                                disvally=DisQueryVally(neiID,hID,Neighbors,Label);
                                if(Label[neiID][hID]<disvally){
                                    WaitProTem.insert(neiID);
                                    ChangeTem[neiID].push_back(make_pair(hID, neiDis+hDis));
                                    Label[neiID][hID]=disvally;//may not be the final correct value at this moment

                                    //cout<<"1--1 "<<neiID<<" "<<hID<<" "<<disvally<<endl;
                                }
                            }

                        }
                    }

                    //Change->ChangeP: AL1->AL2
                    outdatedPruning.clear();
                    if(PruningPointNew[curID].find(hID)!=PruningPointNew[curID].end()){//check whether the pruned vertex should be inserted to curID again
//                        for(int snum=0;snum<PruningPointNew[curID][hID].size();snum++){
                        for(auto snum=PruningPointNew[curID][hID].begin();snum!=PruningPointNew[curID][hID].end();++snum){
                            int s=*snum;

                            /*if(PruningPoint.find(make_pair(curID,hID))!=PruningPoint.end()){
                                for(int snum=0;snum<PruningPoint[make_pair(curID,hID)].size();snum++){
                                    int s=PruningPoint[make_pair(curID,hID)][snum];*/
//                            if(NodeOrder[s]<NodeOrder[curID] && NoSupportedPair.find(make_pair(s,curID))==NoSupportedPair.end()){//it is not in NoSupportedPair
                            if(NodeOrder[s]<NodeOrder[curID]){//it is not in NoSupportedPair
                                disvally=DisQueryVally(s,curID,Neighbors,Label);
                                peakPair=DisQueryPeak(s,curID,Label);
                                dispeak=peakPair.first; peakhub=peakPair.second;
                                if(dispeak>disvally){
                                    WaitProPTem.insert(s);
                                    ChangePTem[s].push_back(curID);
                                    Label[s][curID]=disvally;
//                                    NoSupportedPair.insert(make_pair(s,curID));
                                    outdatedPruning.insert(make_tuple(curID,hID,s));
                                    outdatedPruning.insert(make_tuple(hID,curID,s));

                                    if(ifDebug){
                                        int tempD = DijkstraCore(s,curID);
                                        if(tempD != Label[s][curID]){
//                                            cout<<s<<" "<<curID<<": "<<Label[s][curID]<<" "<<tempD<<" "<<Dijkstra(s,curID,Neighbor)<<endl;//////
                                            cout<<"AL1 to AL2. "<<s<<" "<<curID<<": "<<Label[s][curID]<<" "<<tempD<<endl;
                                        }
                                    }

                                    //cout<<"1--2 "<<s<<" "<<curID<<" "<<disvally<<endl;
                                }
                                else {//if dispeak<=disvally
                                    if(Label[s].find(curID) != Label[s].end()) {
                                        if (dispeak != Label[s][curID]) {
                                            Label[s][curID] = dispeak;
                                        }
                                    }

                                    if(peakhub != -1 && peakhub != hID){
                                        outdatedPruning.insert(make_tuple(curID,hID,s));
                                        outdatedPruning.insert(make_tuple(hID,curID,s));
                                        if(PruningPointNew[curID].find(peakhub) == PruningPointNew[curID].end()){//if not found
                                            PruningPointNew[curID].insert({peakhub,set<int>()});
                                            PruningPointNew[peakhub].insert({curID,set<int>()});
                                        }
                                        PruningPointNew[curID][peakhub].insert(s);
                                        PruningPointNew[peakhub][curID].insert(s);

                                    }
                                }

                            }else if(NodeOrder[s]>NodeOrder[curID]){// && NoSupportedPair.find(make_pair(curID,s))==NoSupportedPair.end()
                                disvally=DisQueryVally(curID,s,Neighbors,Label);
                                peakPair=DisQueryPeak(curID,s,Label);//the old P2 distance dominates P1 distance, but the increase of d(curID,hID) may violate the previous balance
                                dispeak=peakPair.first; peakhub=peakPair.second;
                                if(dispeak>disvally){
                                    WaitProPTem.insert(curID);
                                    ChangePTem[curID].push_back(s);
                                    Label[curID][s]=disvally;//should be the final correct value
//                                    NoSupportedPair.insert(make_pair(curID,s));//insert to indicate that there is no pruning between curID and s
                                    outdatedPruning.insert(make_tuple(curID,hID,s));
                                    outdatedPruning.insert(make_tuple(hID,curID,s));

                                    if(ifDebug){
                                        int tempD = DijkstraCore(curID,s);
                                        if(tempD != Label[curID][s]){
//                                            cout<<curID<<" "<<s<<": "<<Label[curID][s]<<" "<<tempD<<" "<<Dijkstra(curID,s,Neighbor)<<endl;//////
                                            cout<<"AL1 to AL2. "<<curID<<" "<<s<<": "<<Label[curID][s]<<" "<<tempD<<endl;
                                        }
                                    }

                                    //cout<<"1--2 "<<curID<<" "<<s<<" "<<disvally<<endl;
                                }
                                else {//if dispeak<=disvally
                                    if(Label[curID].find(s) != Label[curID].end()) {
                                        if (dispeak != Label[curID][s]) {
                                            Label[curID][s] = dispeak;
                                        }
                                    }

                                    if(peakhub != -1 && peakhub != hID) {
                                        outdatedPruning.insert(make_tuple(curID, hID, s));
                                        outdatedPruning.insert(make_tuple(hID, curID, s));
                                        if (PruningPointNew[curID].find(peakhub) ==
                                            PruningPointNew[curID].end()) {//if not found
                                            PruningPointNew[curID].insert({peakhub, set<int>()});
                                            PruningPointNew[peakhub].insert({curID, set<int>()});
                                        }
                                        PruningPointNew[curID][peakhub].insert(s);
                                        PruningPointNew[peakhub][curID].insert(s);
                                    }
                                }

                            }
                        }
                    }
                    //check higher-order vertex of the affected distance label
                    if(!PruningPointNew[hID].empty()){
                        for(auto it=PruningPointNew[hID].begin();it!=PruningPointNew[hID].end();++it){
                            int lowID = it->first;
                            if(NodeOrder[lowID] < NodeOrder[curID]){
                                for(auto it2=it->second.begin();it2!=it->second.end();++it2){
                                    if(*it2 == curID){
                                        disvally=DisQueryVally(lowID,curID,Neighbors,Label);
                                        peakPair=DisQueryPeak(lowID,curID,Label);
                                        dispeak=peakPair.first; peakhub=peakPair.second;
                                        if(dispeak > disvally){
                                            WaitProPTem.insert(lowID);
                                            ChangePTem[lowID].push_back(curID);
                                            Label[lowID][curID] = disvally;
                                            outdatedPruning.insert(make_tuple(lowID,hID,curID));
                                            outdatedPruning.insert(make_tuple(hID,lowID,curID));
                                        }
                                        else {
                                            if(Label[lowID].find(curID) != Label[lowID].end()) {
                                                if (dispeak != Label[lowID][curID]) {
                                                    Label[lowID][curID] = dispeak;
                                                }
                                            }

                                            if(peakhub != -1 && peakhub != hID){
                                                outdatedPruning.insert(make_tuple(lowID,hID,curID));
                                                outdatedPruning.insert(make_tuple(hID,lowID,curID));
                                                if(PruningPointNew[lowID].find(peakhub) == PruningPointNew[lowID].end()){//if not found
                                                    PruningPointNew[lowID].insert({peakhub,set<int>()});
                                                    PruningPointNew[peakhub].insert({lowID,set<int>()});
                                                }
                                                PruningPointNew[lowID][peakhub].insert(curID);
                                                PruningPointNew[peakhub][lowID].insert(curID);

                                            }
                                        }

                                    }
                                }


                            }

                        }
                    }


                    for(auto it=outdatedPruning.begin();it!=outdatedPruning.end();++it){
                        PruningPointNew[get<0>(*it)][get<1>(*it)].erase(get<2>(*it));
//                        PruningPoint[get<1>(*it)][get<0>(*it)].erase(get<2>(*it));
                    }

                }
            }

            //ChangeP->CHangeP: AL2->AL2
            int v,u,neiid,neiw;
            outdatedPruning.clear();
            for(auto itp=WaitProP.begin();itp!=WaitProP.end();itp++){
                v=*itp;
                for(int k=0;k<ChangeP[v].size();k++){
                    u=ChangeP[v][k];
                    for(int l=0;l<Neighbors[v].size();l++){
                        neiid=Neighbors[v][l].first; neiw=Neighbors[v][l].second;
                        if(NodeOrder[neiid]<NodeOrder[u]){
                            disvally=DisQueryVally(neiid, u,Neighbors,Label);
                            peakPair=DisQueryPeak(neiid, u,Label);
                            dispeak=peakPair.first; peakhub=peakPair.second;
                            if(dispeak>disvally){
                                if(Label[neiid].find(u)==Label[neiid].end() || (Label[neiid].find(u)!=Label[neiid].end() && Label[neiid][u]>disvally)){//if not found or found but disvally is lower
                                    WaitProPTem.insert(neiid);
                                    ChangePTem[neiid].push_back(u);
                                    Label[neiid][u]=disvally;//should be the final correct value
//                                    NoSupportedPair.insert(make_pair(neiid,u));
//                                    outdatedPruning.insert(make_tuple(neiid,v,u));//
//                                    outdatedPruning.insert(make_tuple(v,neiid,u));//
                                    if(ifDebug){
                                        int tempD = DijkstraCore(neiid,u);
                                        if(tempD != Label[neiid][u]){
//                                                cout<<neiid<<" "<<u<<": "<<Label[neiid][u]<<" "<<tempD<<" "<<Dijkstra(neiid,u,Neighbor)<<endl;
                                            cout<<"AL2 to AL2. "<<neiid<<" "<<u<<": "<<Label[neiid][u]<<" "<<tempD<<endl;
                                        }
                                    }
                                    //cout<<"2--2 "<<neiid<<" "<<u<<" "<<disvally<<endl;
                                }
                            }
                            else {//if dispeak<=disvally
                                if(Label[neiid].find(u) != Label[neiid].end()) {
                                    if (dispeak != Label[neiid][u]) {
                                        Label[neiid][u] = dispeak;
                                    }
                                }

                                if(peakhub != -1 && peakhub != u){
//                                    outdatedPruning.insert(make_tuple(neiid,v,u));//
//                                    outdatedPruning.insert(make_tuple(v,neiid,u));//
                                    if(PruningPointNew[neiid].find(peakhub) == PruningPointNew[neiid].end()){//if not found
                                        PruningPointNew[neiid].insert({peakhub,set<int>()});
                                        PruningPointNew[peakhub].insert({neiid,set<int>()});
                                    }
                                    PruningPointNew[neiid][peakhub].insert(u);
                                    PruningPointNew[peakhub][neiid].insert(u);

                                }
                            }


                        }
                    }
                }
            }
            for(auto it=outdatedPruning.begin();it!=outdatedPruning.end();++it){
                if(PruningPointNew[get<0>(*it)].find(get<1>(*it)) != PruningPointNew[get<0>(*it)].end()){//if found
                    PruningPointNew[get<0>(*it)][get<1>(*it)].erase(get<2>(*it));
//                        PruningPoint[get<1>(*it)][get<0>(*it)].erase(get<2>(*it));
                }
            }

            WaitPro=WaitProTem;
            Change=ChangeTem;
            WaitProP=WaitProPTem;
            ChangeP=ChangePTem;
        }

    }
}


void Graph::DecreasePSL(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label){
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
	int Dab=ShortestDisQuery(a,b,Label);//query by PSL label

	int LID,HID;
	if(NodeOrder[a]>NodeOrder[b]){
		LID=b; HID=a;
	}else{
		LID=a; HID=b;
	}

	if(Dab>newW){//the index change is triggered
		vector<vector<pair<int,int>>> Change;
		vector<pair<int,int>> vec;
		Change.assign(nodenum,vec);
		set<int> WaitPro;

		Label[LID][HID]=newW;
		Change[LID].push_back(make_pair(HID,newW));
		WaitPro.insert(LID);

		//check the label of a,b
		int hubid, hubdis;
		unordered_map<int,int>::iterator it1=Label[LID].begin();
		for(;it1!=Label[LID].end();it1++){
			hubid=(*it1).first; hubdis=(*it1).second;
			if(NodeOrder[hubid]>NodeOrder[HID] && newW+hubdis<ShortestDisQuery(HID,hubid,Label)){
				Label[HID][hubid]=newW+hubdis;
				Change[HID].push_back(make_pair(hubid, newW+hubdis));
				WaitPro.insert(HID);
			}
		}
		unordered_map<int,int>::iterator it2=Label[HID].begin();
		for(;it2!=Label[HID].end();it2++){
			hubid=(*it2).first; hubdis=(*it2).second;
			if(newW+hubdis<ShortestDisQuery(LID, hubid,Label)){
				Label[LID][hubid]=newW+hubdis;
				Change[LID].push_back(make_pair(hubid, newW+hubdis));
				WaitPro.insert(LID);
			}
		}

		//check the label of their neighbors step by step
		while(WaitPro.size()>0){
			set<int> WaitProTem;
			vector<vector<pair<int,int>>> ChangeTem;
			vector<pair<int,int>> vec;
			ChangeTem.assign(nodenum,vec);

			for(set<int>::iterator it=WaitPro.begin();it!=WaitPro.end();it++){
				int curID=*it;
				vector<pair<int,int>> curChange=Change[curID];
				int neiID, neiDis, hID, hDis;
				for(int i=0;i<Neighbors[curID].size();i++){
					neiID=Neighbors[curID][i].first;
					neiDis=Neighbors[curID][i].second;

					for(int j=0;j<curChange.size();j++){
						hID=curChange[j].first; hDis=curChange[j].second;
						if(NodeOrder[hID]>NodeOrder[neiID] && ShortestDisQuery(neiID, hID,Label)>neiDis+hDis){
							Label[neiID][hID]=neiDis+hDis;
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

int Graph::ShortestDisQuery(int ID1,int ID2,vector<unordered_map<int,int>> &Label){
	int d=INF;
	unordered_map<int,int>::iterator it;
	int hub, dis1, dis2;
	for(it=Label[ID1].begin();it!=Label[ID1].end();it++){
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
//function of computing P1
int Graph::DisQueryVally2(int ID1, int ID2, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label){
	int neiID,neiDis;
	int d=INF;
	for(int i=0;i<Neighbors[ID1].size();i++){
		neiID=Neighbors[ID1][i].first;
		neiDis=Neighbors[ID1][i].second;
		if(NodeOrder[neiID]<=NodeOrder[ID2] && Label[neiID].find(ID2)!=Label[neiID].end()){//if found
            //new
            if(NodeOrder[neiID] < NodeOrder[ID1] && Label[ID1].find(ID2)!=Label[ID1].end()){
                bool flag_continue = false;
                if(Label[neiID][ID2] == neiDis + Label[ID1][ID2]){
                    for(auto it=Label[neiID].begin();it!=Label[neiID].end();++it){
                        if(NodeOrder[it->first]<NodeOrder[ID2] && Label[it->first].find(ID2)!=Label[it->first].end()){
                            if(it->second + Label[it->first][ID2] == Label[neiID][ID2]){
                                flag_continue = true;
                                break;
                            }
                        }else if(NodeOrder[it->first]>NodeOrder[ID2] && Label[ID2].find(it->first)!=Label[ID2].end()){
                            if(it->second + Label[ID2][it->first] == Label[neiID][ID2]){
                                flag_continue = true;
                                break;
                            }
                        }
                    }
                    if(flag_continue)
                        continue;
                }
            }

			if(neiDis+Label[neiID][ID2]<d){
				d=neiDis+Label[neiID][ID2];
			}
		}
	}
	return d;
}
//old version. function of computing the shortest distance through the label of ID1's neighbors, i.e., d1.
int Graph::DisQueryVally(int ID1, int ID2, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label){
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

//function of computing d2, new version
pair<int,int> Graph::DisQueryPeak(int ID1, int ID2,vector<unordered_map<int,int>> &Label){
	int d=INF;
	unordered_map<int,int>::iterator it;
	int hub, dis1, finalHub=-1;
	for(it=Label[ID1].begin();it!=Label[ID1].end();it++){
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
//old version: function of computing d2
//int Graph::DisQueryPeak(int ID1, int ID2,vector<unordered_map<int,int>> &Label){
//    int d=INF;
//    unordered_map<int,int>::iterator it;
//    int hub, dis1;
//    for(it=Label[ID1].begin();it!=Label[ID1].end();it++){
//        hub=(*it).first;
//        dis1=(*it).second;
//        if(hub!=ID2 && Label[ID2].find(hub)!=Label[ID2].end()){//hub is the hub of ID1 and ID2 and hub != ID2
//            if(dis1+Label[ID2][hub]<d){
//                d=dis1+Label[ID2][hub];
//            }
//        }
//    }
//    return d;
//}
//function of computing distance from ID1(LID) to ID2(HID), via a neighbor of ID1 (which has lower order than HID) on updated graph
int Graph::DisQueryLower1(int ID1, int ID2, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label){
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

