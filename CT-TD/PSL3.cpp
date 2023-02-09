/*
 * PSL.cpp
 *
 *  Created on: 13 Oct 2022
 *     Authors: zhangmengxuan, Xinjie ZHOU
 */
#include "head3.h"

void Graph::BVCIndexConstructMThread(){
    unordered_map<int,int> map0; map0.clear();
    Label.assign(nodenum, map0);
    Dhop.assign(nodenum, map0);

    PruningPointList.clear();
    unordered_map<int,list<int>> map3; map3.clear();
    PruningPointList.assign(nodenum,map3);

    DvertexNew.assign(nodenum, true);

    for(int i=0;i<nodenum;i++){
        Label[i].insert(make_pair(i,0));
        Dhop[i].insert(make_pair(i,0));
        //Dvectex.push_back(i);
    }

    for(int i=0;i<5;++i){
        cout<<"Order "<<nodenum-i-1<<" : "<<vNodeOrder[nodenum-i-1]<<" "<<AdjaCore[vNodeOrder[nodenum-i-1]].size()<<endl;
    }

    int batchSize = 512;
    int batchNum = ceil((double)nodenum/batchSize);
    int stepShow = batchNum / 100;
    for(int b1=0;b1<batchNum;++b1){


        vector<int> bNodes;
        int maxR = min((b1+1)*batchSize,nodenum);
        for(int i=b1*batchSize;i<maxR;++i){
            bNodes.push_back(vNodeOrder[nodenum-i-1]);
        }
        // process each batch
        bool flag=true;
        int step=0;
        while(flag){
            if(b1%stepShow==0 && step%5==0){
                cout<<"Batch "<<b1<<": step "<<step<<" finish!"<<endl;
            }

            flag=BVCDhopLableRefreshMulti(bNodes, step+1);
            step+=1;
        }
    }


    cout<<"Index finish construction"<<endl;
}

bool Graph::BVCDhopLableRefreshMulti(vector<int>& bNodes, int step){
    bool flag=false;
    vector<unordered_map<int,int>> newDhop;
    unordered_map<int,int> m0; m0.clear();
    newDhop.assign(nodenum,m0);
    vector<int> newDvec;
    vector<vector<tri>> PP;
    vector<tri> prp; prp.clear();
    PP.assign(threadnum, prp);

    if(true){//use multiple thread
//    if(false){
        boost::thread_group thread;
        vector<vector<int>> ProcessID;
        vector<int> vvv; ProcessID.assign(threadnum,vvv);
        threadDistribute2(bNodes,ProcessID);

        for(int i=0;i<ProcessID.size();i++){
            vector<int> p=ProcessID[i];
            thread.add_thread(new boost::thread(&Graph::labelMultiThread2New, this, boost::ref(newDhop), p,step));
        }
        thread.join_all();
    }else{//use single thread
        vector<vector<int>> ProcessID;
        vector<int> vvv; ProcessID.assign(threadnum,vvv);
        threadDistribute2(bNodes,ProcessID);

        for(int i=0;i<ProcessID.size();i++){
            vector<int> p=ProcessID[i];
            labelMultiThread2New(newDhop,p,step);
        }
    }


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

void Graph::threadDistribute2(vector<int>& bNodes, vector<vector<int>>& processID){
    int ID;
    int cnt=0;
    int threadOrder;
    int a;

    for(int r=0;r<bNodes.size();r++){
//        ID=vNodeOrder[r];
        ID=bNodes[r];
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

void Graph::IndexConstructMThread2New(){
	unordered_map<int,int> map0; map0.clear();
	Label.assign(nodenum, map0);
	Dhop.assign(nodenum, map0);
//	PruningPointNew.clear();
//	unordered_map<int,vector<int>> map1;
//	PruningPointNew.assign(nodenum,map1);

//    PruningPoint.clear();
//    unordered_map<int,set<int>> map1; map1.clear();
//    PruningPoint.assign(nodenum,map1);

    PruningPointList.clear();
    unordered_map<int,list<int>> map3; map3.clear();
    PruningPointList.assign(nodenum,map3);

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

    if(true){//use multiple thread
//    if(false){
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
    }else{//use single thread
        vector<vector<int>> ProcessID;
        vector<int> vvv; ProcessID.assign(threadnum,vvv);
        threadDistribute(ProcessID);

        for(int i=0;i<ProcessID.size();i++){
            vector<int> p=ProcessID[i];
            labelMultiThread2New(newDhop,p,step);
        }
    }


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
						ShortestDisQuery1(nodeID, hub,SupNode,TempDis); //query by label, TempDis is the result
                        //cout<<"nodeID "<<nodeID<<" TempDis "<<TempDis<<" d "<<d<<endl;
						//ShortestDisQuery2(nodeID, hub,SupNode,TempDis,d);
						if(TempDis>d){
							if(Dhop0.find(hub)!=Dhop0.end()){
								if(Dhop0[hub]>d)
                                    Dhop0[hub]=d;
							}else{
								Dhop0.insert(make_pair(hub,d));
							}
						}
						else{//although d is not the shortest distance, we still need to record it. if Q(hub,nodeID,L<h) <= e(nodeID,neighID) + Lh-1(neighID)[hub], i.e., Q(w,v,L<h) <= e(v,u) + Lh-1(u)[w]
							for(int k=0;k<SupNode.size();k++){
								int supn=SupNode[k];

								if(supn!=nodeID && supn!=hub){
									vSm[nodeID]->wait();
//									PruningPointNew[nodeID][supn].push_back(hub);//((v,c),w)
//									PruningPointStepNew[nodeID][supn].insert(make_pair(hub,step));
//                                    PruningPoint[nodeID][supn].insert(hub);
                                    PruningPointList[nodeID][supn].push_back(hub);

									vSm[nodeID]->notify();

									vSm[hub]->wait();
//									PruningPointNew[hub][supn].push_back(nodeID);//((w,c),v)
//									PruningPointStepNew[hub][supn].insert(make_pair(nodeID,step));
//                                    PruningPoint[hub][supn].insert(nodeID);
                                    PruningPointList[hub][supn].push_back(nodeID);

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
//original old version
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
    cout<<LID<<"("<<NodeOrder[LID]<<") "<<HID<<"("<<NodeOrder[HID]<<")"<<endl;
	int dis,disvally,dispeak,peakhub;
	//activate or not
	dis=DisQueryLower1(LID,HID,Neighbors,Label);//d1', via a neighbor of LID (which has lower order than HID) to reach HID
	dispeak=DisQueryPeak(LID,HID,Label);//d2, via the common hub vertex of LID and HID
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
                                    if(neiID == 207450 && hID == 211814){
                                        cout<<"AL1->AL1: "<<neiID<<" "<<hID<<" "<<disvally<<endl;
                                    }
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
								peakPair=DisQueryPeak2(s,curID,Label);
                                dispeak=peakPair.first; //peakhub=peakPair.second;
								if(dispeak>disvally){
									WaitProPTem.insert(s);
									ChangePTem[s].push_back(curID);
									Label[s][curID]=disvally;
									NoSupportedPair.insert(make_pair(s,curID));
                                    if(s == 207450 && curID == 211814){
                                        cout<<"AL1->AL2 1: "<<s<<" "<<curID<<" "<<disvally<<endl;
                                    }

									//cout<<"1--2 "<<s<<" "<<curID<<" "<<disvally<<endl;
								}
							}else if(NodeOrder[s]>NodeOrder[curID] && NoSupportedPair.find(make_pair(curID,s))==NoSupportedPair.end()){
								disvally=DisQueryVally(curID,s,Neighbors,Label);
                                pair<int,int> peakPair;
                                peakPair=DisQueryPeak2(curID,s,Label);//the old P2 distance dominates P1 distance, but the increase of d(curID,hID) may violate the previous balance
                                dispeak=peakPair.first; //peakhub=peakPair.second;
                                if(dispeak>disvally){
									WaitProPTem.insert(curID);
									ChangePTem[curID].push_back(s);
									Label[curID][s]=disvally;//should be the final correct value
									NoSupportedPair.insert(make_pair(curID,s));//insert to indicate that there is no pruning between curID and s
                                    if(curID == 207450 && s == 211814){
                                        cout<<"AL1->AL2 2: "<<curID<<" "<<s<<" "<<disvally<<endl;
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
							dispeak=DisQueryPeak2(neiid, u,Label).first;
							if(disvally<dispeak){
//                                if(Label[neiid].find(u)==Label[neiid].end() || (Label[neiid].find(u)!=Label[neiid].end())){/// if not found or found but disvally is lower
								if(Label[neiid].find(u)==Label[neiid].end() || (Label[neiid].find(u)!=Label[neiid].end() && Label[neiid][u]>disvally)){/// if not found or found but disvally is lower

//                                    cout<<neiid<<" "<<u<<" "<<disvally;
//                                    if(Label[neiid].find(u)!=Label[neiid].end()){//if found
//                                        cout<<" "<<Label[neiid][u];
//                                    }
//                                    cout<<endl;

                                    WaitProPTem.insert(neiid);
                                    ChangePTem[neiid].push_back(u);
                                    Label[neiid][u]=disvally;//should be the final correct value
                                    NoSupportedPair.insert(make_pair(neiid,u));///

                                    if(neiid == 207450 && u == 211814){
                                        cout<<"AL2->AL2: "<<neiid<<" "<<u<<" "<<disvally<<" "<<dispeak<<" "<<DijkstraCore(neiid,u)<<endl;
                                    }

//                                    cout<<neiid<<" "<<u<<" "<<Label[neiid][u]<<endl;
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

//original old version
void Graph::IncreasePSL(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label,vector<unordered_map<int,vector<int>>> &PruningPointNew,set<pair<int,int>> &NoSupportedPair){
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
    cout<<LID<<"("<<NodeOrder[LID]<<") "<<HID<<"("<<NodeOrder[HID]<<")"<<endl;
    int dis,disvally,dispeak,peakhub;
    //activate or not
    dis=DisQueryLower1(LID,HID,Neighbors,Label);//d1', via a neighbor of LID (which has lower order than HID) to reach HID
    dispeak=DisQueryPeak(LID,HID,Label);//d2, via the common hub vertex of LID and HID
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
        ChangeP.assign(nodenum,vector<int>());
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
                                    if(neiID == 207450 && hID == 211814){
                                        cout<<"AL1->AL1: "<<neiID<<" "<<hID<<" "<<disvally<<endl;
                                    }
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
                                peakPair=DisQueryPeak2(s,curID,Label);
                                dispeak=peakPair.first; //peakhub=peakPair.second;
                                if(dispeak>disvally){
                                    WaitProPTem.insert(s);
                                    ChangePTem[s].push_back(curID);
                                    Label[s][curID]=disvally;
                                    NoSupportedPair.insert(make_pair(s,curID));
                                    if(s == 207450 && curID == 211814){
                                        cout<<"AL1->AL2 1: "<<s<<" "<<curID<<" "<<disvally<<endl;
                                    }

                                    //cout<<"1--2 "<<s<<" "<<curID<<" "<<disvally<<endl;
                                }
                            }else if(NodeOrder[s]>NodeOrder[curID] && NoSupportedPair.find(make_pair(curID,s))==NoSupportedPair.end()){
                                disvally=DisQueryVally(curID,s,Neighbors,Label);
                                pair<int,int> peakPair;
                                peakPair=DisQueryPeak2(curID,s,Label);//the old P2 distance dominates P1 distance, but the increase of d(curID,hID) may violate the previous balance
                                dispeak=peakPair.first; //peakhub=peakPair.second;
                                if(dispeak>disvally){
                                    WaitProPTem.insert(curID);
                                    ChangePTem[curID].push_back(s);
                                    Label[curID][s]=disvally;//should be the final correct value
                                    NoSupportedPair.insert(make_pair(curID,s));//insert to indicate that there is no pruning between curID and s
                                    if(curID == 207450 && s == 211814){
                                        cout<<"AL1->AL2 2: "<<curID<<" "<<s<<" "<<disvally<<endl;
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
                            dispeak=DisQueryPeak2(neiid, u,Label).first;
                            if(disvally<dispeak){
                                if(Label[neiid].find(u)==Label[neiid].end()){/// if not found or found but disvally is lower
//                                if(Label[neiid].find(u)==Label[neiid].end() || (Label[neiid].find(u)!=Label[neiid].end() && Label[neiid][u]>disvally)){/// if not found or found but disvally is lower

//                                    cout<<neiid<<" "<<u<<" "<<disvally;
//                                    if(Label[neiid].find(u)!=Label[neiid].end()){//if found
//                                        cout<<" "<<Label[neiid][u];
//                                    }
//                                    cout<<endl;

                                    WaitProPTem.insert(neiid);
                                    ChangePTem[neiid].push_back(u);
                                    Label[neiid][u]=disvally;//should be the final correct value

                                    NoSupportedPair.insert(make_pair(neiid,u));///

                                    if(neiid == 207450 && u == 211814){
                                        cout<<"AL2->AL2: "<<neiid<<" "<<u<<" "<<disvally<<" "<<dispeak<<" "<<DijkstraCore(neiid,u)<<endl;
                                    }

//                                    cout<<neiid<<" "<<u<<" "<<Label[neiid][u]<<endl;
                                    //cout<<"2--2 "<<neiid<<" "<<u<<" "<<disvally<<endl;
                                }
                                else{//if found
                                    if(Label[neiid][u]>disvally){// and they are unequal
                                        WaitProPTem.insert(neiid);
                                        ChangePTem[neiid].push_back(u);
                                        Label[neiid][u]=disvally;//should be the final correct value
                                        NoSupportedPair.insert(make_pair(neiid,u));///
                                    }else if(Label[neiid][u]<disvally){

                                    }


                                    if(neiid == 207450 && u == 211814){
                                        cout<<"AL2->AL2 2: "<<neiid<<" "<<u<<" "<<disvally<<" "<<dispeak<<" "<<DijkstraCore(neiid,u)<<endl;
                                    }
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
}

//old version, vector version with NoSuportedPair
void Graph::IncreasePSL(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label,vector<unordered_map<int,vector<int>>> &PruningPointNew){
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
    dispeak=DisQueryPeak(LID,HID,Label);//d2, via the common hub vertex of LID and HID
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
                                peakPair=DisQueryPeak2(s,curID,Label);
                                dispeak=peakPair.first; //peakhub=peakPair.second;
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
                                peakPair=DisQueryPeak2(curID,s,Label);//the old P2 distance dominates P1 distance, but the increase of d(curID,hID) may violate the previous balance
                                dispeak=peakPair.first; //peakhub=peakPair.second;
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
                            dispeak=DisQueryPeak2(neiid, u,Label).first;
                            if(disvally<dispeak){
                                if(Label[neiid].find(u)==Label[neiid].end() || (Label[neiid].find(u)!=Label[neiid].end() && Label[neiid][u]>disvally)){/// if not found or found but disvally is lower
                                    WaitProPTem.insert(neiid);
                                    ChangePTem[neiid].push_back(u);
                                    Label[neiid][u]=disvally;//should be the final correct value
                                    NoSupportedPair.insert(make_pair(neiid,u));///

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
}

//new version: list version with NoSupportedPair
void Graph::IncreasePSL2(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label,vector<unordered_map<int,list<int>>> &PruningPointNew){
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

    int LID,HID;
    if(NodeOrder[a]>NodeOrder[b]){
        LID=b; HID=a;
    }else{
        LID=a; HID=b;
    }

    int dis,disvally,dispeak,peakhub;
    pair<int,int> peakPair;//distance, hubID
    set<tuple<int,int,int>> outdatedPruning;//<nodeID,supportNode,prunedID>
    map<pair<int,int>,int> newPruningPoints;
    //activate or not
    dis=DisQueryLower1(LID,HID,Neighbors,Label);

    //if(dispeak<=oldW) return;
    if(Label[LID].find(HID)!=Label[LID].end() && dis>oldW){//index update triggered
        vector<vector<pair<int,int>>> Change;
        vector<pair<int,int>> vec;
        Change.assign(nodenum,vec);
        set<int> WaitPro;
        vector<vector<int>> ChangeP;
        vector<int> vecint;
        ChangeP.assign(nodenum,vecint);
//        vector<vector<pair<int,int>>> ChangeP;
//        ChangeP.assign(nodenum,vector<pair<int,int>>());
        set<int> WaitProP;

        WaitPro.insert(LID);
        Change[LID].push_back(make_pair(HID, oldW));
        disvally=DisQueryVally(LID,HID,Neighbors,Label);
        Label[LID][HID]=disvally;//correct to the new value

        NoSupportedPair.clear();

        //affected by the w(a.b)
        int hubID, hDis;
        int dis, cnt;
        for(unordered_map<int,int>::iterator it=Label[HID].begin();it!=Label[HID].end();it++){
            hubID=(*it).first; hDis=(*it).second;
            if(hubID!=HID && Label[LID].find(hubID)!=Label[LID].end() && oldW+hDis==Label[LID][hubID]){
                disvally=DisQueryVally(LID,hubID,Neighbors,Label);
                if(Label[LID][hubID]<disvally){
                    WaitPro.insert(LID);
                    Change[LID].push_back(make_pair(hubID, oldW+hDis));
                    Label[LID][hubID]=disvally;

                }
            }
        }

        for(unordered_map<int,int>::iterator it=Label[LID].begin();it!=Label[LID].end();it++){
            hubID=(*it).first; hDis=(*it).second;
            if(Label[HID].find(hubID)!=Label[HID].end() && oldW+hDis==Label[HID][hubID]){
                disvally=DisQueryVally(HID,hubID,Neighbors,Label);
                if(Label[HID][hubID]<disvally){
                    WaitPro.insert(HID);
                    Change[HID].push_back(make_pair(hubID, oldW+hDis));
                    Label[HID][hubID]=disvally;

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
//            vector<vector<pair<int,int>>> ChangePTem;
//            ChangePTem.assign(nodenum,vector<pair<int,int>>());


            //Change->Change & ChangeP
            for(set<int>::iterator it=WaitPro.begin();it!=WaitPro.end();it++){
                int curID=*it;
                vector<pair<int,int>> curChange=Change[curID];
                int neiID, neiDis, hID, hDis;

                for(int j=0;j<curChange.size();j++){
                    hID=curChange[j].first; hDis=curChange[j].second;
                    //Change->Change
                    for(int k=0;k<Neighbors[curID].size();k++){
                        neiID=Neighbors[curID][k].first; neiDis=Neighbors[curID][k].second;

                        if(Label[neiID].find(hID)!=Label[neiID].end() && neiDis+hDis==Label[neiID][hID]){
                            disvally=DisQueryVally(neiID,hID,Neighbors,Label);
                            if(Label[neiID][hID]<disvally){
                                WaitProTem.insert(neiID);
                                ChangeTem[neiID].push_back(make_pair(hID, neiDis+hDis));
                                Label[neiID][hID]=disvally;

                            }
                        }
                    }

                    //Change->ChangeP
                    if(PruningPointNew[curID].find(hID)!=PruningPointNew[curID].end()){
//                        outdatedPruning.clear();
                        for(auto snum=PruningPointNew[curID][hID].begin();snum!=PruningPointNew[curID][hID].end();++snum){//for each pruned vertex
                            int s=*snum;

                            if(NodeOrder[s]<NodeOrder[curID] && NoSupportedPair.find(make_pair(s,curID))==NoSupportedPair.end()){
//                            if(NodeOrder[s]<NodeOrder[curID]){

                                disvally=DisQueryVally(s,curID,Neighbors,Label);
//                                dispeak=DisQueryPeak(s,curID,Label);
                                peakPair=DisQueryPeak2(s,curID,Label);//the original dispeak is smaller than disvally
                                dispeak=peakPair.first; peakhub=peakPair.second;
                                if(dispeak>disvally){
                                    WaitProPTem.insert(s);
                                    ChangePTem[s].push_back(curID);
//                                    ChangePTem[s].push_back(make_pair(curID,dispeak));
                                    Label[s][curID]=disvally;

                                    outdatedPruning.insert(make_tuple(curID,hID,s));
                                    outdatedPruning.insert(make_tuple(s,hID,curID));

                                    NoSupportedPair.insert(make_pair(s,curID));
                                }
                                else {//if dispeak<=disvally
                                    if(Label[s].find(curID) != Label[s].end()) {
                                        cout<<"Remove label1 "<<s<<" "<<curID<<": "<<Label[s][curID]<<endl;
                                        Label[s].erase(curID);
                                    }

                                    if(peakhub != -1 && peakhub != hID){
                                        newPruningPoints[make_pair(s,curID)] = peakhub;
//                                        outdatedPruning.insert(make_tuple(curID,hID,s));
//                                        outdatedPruning.insert(make_tuple(s,hID,curID));
//
//                                        PruningPointNew[curID][peakhub].push_back(s);
//                                        PruningPointNew[s][peakhub].push_back(curID);

                                    }
                                }
                            }else if(NodeOrder[s]>NodeOrder[curID] && NoSupportedPair.find(make_pair(curID,s))==NoSupportedPair.end()){
//                            } else if(NodeOrder[s]>NodeOrder[curID]){

                                disvally=DisQueryVally(curID,s,Neighbors,Label);
//                                dispeak=DisQueryPeak(curID,s,Label);
                                peakPair=DisQueryPeak2(curID,s,Label);//the old P2 distance dominates P1 distance, but the increase of d(curID,hID)
                                dispeak=peakPair.first; peakhub=peakPair.second;
                                if(dispeak>disvally){
                                    WaitProPTem.insert(curID);
                                    ChangePTem[curID].push_back(s);
//                                    ChangePTem[curID].push_back(make_pair(s,dispeak));
                                    Label[curID][s]=disvally;
                                    outdatedPruning.insert(make_tuple(curID,hID,s));
                                    outdatedPruning.insert(make_tuple(s,hID,curID));

                                    NoSupportedPair.insert(make_pair(curID,s));
                                }
                                else {//if dispeak<=disvally
                                    if(Label[curID].find(s) != Label[curID].end()) {
                                        cout<<"Remove label2 "<<curID<<" "<<s<<": "<<Label[curID][s]<<endl;
                                        Label[curID].erase(s);
                                    }

                                    if(peakhub != -1 && peakhub != hID) {
                                        newPruningPoints[make_pair(curID,s)] = peakhub;
//                                        outdatedPruning.insert(make_tuple(curID, hID, s));
//                                        outdatedPruning.insert(make_tuple(s, hID, curID));
//
//                                        PruningPointNew[curID][peakhub].push_back(s);
//                                        PruningPointNew[s][peakhub].push_back(curID);
                                    }
                                }
                            }
                        }
//                        for(auto it2=outdatedPruning.begin();it2!=outdatedPruning.end();++it2){
//                            if(PruningPointNew[get<0>(*it2)].find(get<1>(*it2)) != PruningPointNew[get<0>(*it2)].end()) {//if found
//                                vector<list<int>::iterator> temp_its;
//                                for(auto it=PruningPointNew[get<0>(*it2)][get<1>(*it2)].begin();it!=PruningPointNew[get<0>(*it2)][get<1>(*it2)].end();++it){
//                                    if(*it==get<2>(*it2)){
//                                        temp_its.push_back(it);
//                                    }
//                                }
//                                for(int i=0;i<temp_its.size();++i){
//                                    PruningPointNew[get<0>(*it2)][get<1>(*it2)].erase(temp_its[i]);
//                                }
//                            }
//
//                        }
                    }
                }
            }

            //ChangeP->CHangeP
            int v,u,neiid,neiw,cDis;
            for(set<int>::iterator itp=WaitProP.begin();itp!=WaitProP.end();itp++){
                v=*itp;
                for(int k=0;k<ChangeP[v].size();k++){
                    u=ChangeP[v][k];
//                    u=ChangeP[v][k].first, cDis=ChangeP[v][k].second;
//                    outdatedPruning.clear();
                    for(int l=0;l<Neighbors[v].size();l++){
                        neiid=Neighbors[v][l].first; neiw=Neighbors[v][l].second;

                        if(NodeOrder[neiid]<NodeOrder[u]){
                            disvally=DisQueryVally(neiid, u,Neighbors,Label);
//                            disvally=DisQueryVallyVert(neiid, u,Neighbors,Label);
//                            dispeak=DisQueryPeak(neiid, u,Label);
                            peakPair=DisQueryPeak2(neiid, u,Label);
                            dispeak=peakPair.first; peakhub=peakPair.second;

                            if(disvally<dispeak){
                                if(Label[neiid].find(u)==Label[neiid].end() || (Label[neiid].find(u)!=Label[neiid].end() && Label[neiid][u]!=disvally)){/// !=
                                    WaitProPTem.insert(neiid);
                                    ChangePTem[neiid].push_back(u);
//                                    ChangePTem[neiid].push_back(make_pair(u,dispeak));
                                    Label[neiid][u]=disvally;
//                                    outdatedPruning.insert(make_tuple(neiid,hubID,u));//
//                                    outdatedPruning.insert(make_tuple(u,hubID,neiid));//

                                    NoSupportedPair.insert(make_pair(neiid,u));
                                }
                            }
                            else {//if dispeak<=disvally
//                                if(Label[neiid].find(u) != Label[neiid].end()) {///wrong part
//                                    Label[neiid].erase(u);
//                                }

                                if(peakhub != -1 && peakhub != hubID){
                                    newPruningPoints[make_pair(neiid,u)] = peakhub;
////                                    outdatedPruning.insert(make_tuple(neiid,hubID,hID));//
////                                    outdatedPruning.insert(make_tuple(hID,hubID,neiid));//
//                                    PruningPointNew[neiid][peakhub].push_back(u);
//                                    PruningPointNew[u][peakhub].push_back(neiid);
                                }
                            }
                        }
                    }
//                    for(auto it2=outdatedPruning.begin();it2!=outdatedPruning.end();++it2){
//                        if(PruningPointNew[get<0>(*it2)].find(get<1>(*it2)) != PruningPointNew[get<0>(*it2)].end()){//if found
//                            vector<list<int>::iterator> temp_its;
//                            for(auto it=PruningPointNew[get<0>(*it2)][get<1>(*it2)].begin();it!=PruningPointNew[get<0>(*it2)][get<1>(*it2)].end();++it){
//                                if(*it==get<2>(*it2)){
//                                    temp_its.push_back(it);
//                                }
//                            }
//                            for(int i=0;i<temp_its.size();++i){
//                                PruningPointNew[get<0>(*it2)][get<1>(*it2)].erase(temp_its[i]);
//                            }
//                        }
//                    }
                }
            }

            WaitPro=WaitProTem;
            Change=ChangeTem;
            WaitProP=WaitProPTem;
            ChangeP=ChangePTem;
        }

        ///remove old pruning point
        for(auto it2=outdatedPruning.begin();it2!=outdatedPruning.end();++it2){
            if(PruningPointNew[get<0>(*it2)].find(get<1>(*it2)) != PruningPointNew[get<0>(*it2)].end()){//if found
                vector<list<int>::iterator> temp_its;
                for(auto it=PruningPointNew[get<0>(*it2)][get<1>(*it2)].begin();it!=PruningPointNew[get<0>(*it2)][get<1>(*it2)].end();++it){
                    if(*it==get<2>(*it2)){
                        temp_its.push_back(it);
                    }
                }
                for(int i=0;i<temp_its.size();++i){
                    PruningPointNew[get<0>(*it2)][get<1>(*it2)].erase(temp_its[i]);
                }
            }
        }
        ///add new pruning point
        for(auto it=newPruningPoints.begin();it!=newPruningPoints.end();++it){
            PruningPointNew[it->first.first][it->second].push_back(it->first.second);
            PruningPointNew[it->first.second][it->second].push_back(it->first.first);
        }
    }

}

//new version: list version without NoSupportedPair, correct
void Graph::IncreasePSL(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label,vector<unordered_map<int,list<int>>> &PruningPointNew){
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

    int LID,HID;
    if(NodeOrder[a]>NodeOrder[b]){
        LID=b; HID=a;
    }else{
        LID=a; HID=b;
    }

    int dis,disvally,dispeak,peakhub;
    pair<int,int> peakPair;//distance, hubID
    set<tuple<int,int,int>> outdatedPruning;//<nodeID,supportNode,prunedID>
    //activate or not
    dis=DisQueryLower1(LID,HID,Neighbors,Label);

    //if(dispeak<=oldW) return;
    if(Label[LID].find(HID)!=Label[LID].end() && dis>oldW){//index update triggered
        vector<vector<pair<int,int>>> Change;
        vector<pair<int,int>> vec;
        Change.assign(nodenum,vec);
        set<int> WaitPro;
//        vector<vector<int>> ChangeP;
        vector<vector<pair<int,int>>> ChangeP;
//        vector<int> vecint;
//        ChangeP.assign(nodenum,vecint);
        ChangeP.assign(nodenum,vector<pair<int,int>>());
        set<int> WaitProP;

        WaitPro.insert(LID);
        Change[LID].push_back(make_pair(HID, oldW));
        disvally=DisQueryVally(LID,HID,Neighbors,Label);
        Label[LID][HID]=disvally;//correct to the new value

        //affected by the w(a.b)
        int hubID, hDis;
        int dis, cnt;
        for(unordered_map<int,int>::iterator it=Label[HID].begin();it!=Label[HID].end();it++){
            hubID=(*it).first; hDis=(*it).second;
            if(hubID!=HID && Label[LID].find(hubID)!=Label[LID].end() && oldW+hDis==Label[LID][hubID]){
                disvally=DisQueryVally(LID,hubID,Neighbors,Label);
                if(Label[LID][hubID]<disvally){
                    WaitPro.insert(LID);
                    Change[LID].push_back(make_pair(hubID, oldW+hDis));
                    Label[LID][hubID]=disvally;
                    /// No need of dispeak
//                    peakPair=DisQueryPeak2(LID,hubID, Label);
//                    dispeak=peakPair.first; peakhub=peakPair.second;
//                    if(dispeak>disvally){
//                        Label[LID][hubID]=disvally;//may not be correct at this moment
//                    }
//                    else{
////                            Label[LID][hubID]=dispeak;//may not be correct at this moment
//                        Label[LID].erase(hubID);
//
//                        PruningPointNew[LID][peakhub].insert(hubID);
//                        PruningPointNew[hubID][peakhub].insert(LID);
//                    }
                }
            }
        }

        for(unordered_map<int,int>::iterator it=Label[LID].begin();it!=Label[LID].end();it++){
            hubID=(*it).first; hDis=(*it).second;
            if(Label[HID].find(hubID)!=Label[HID].end() && oldW+hDis==Label[HID][hubID]){
                disvally=DisQueryVally(HID,hubID,Neighbors,Label);
                if(Label[HID][hubID]<disvally){
                    WaitPro.insert(HID);
                    Change[HID].push_back(make_pair(hubID, oldW+hDis));
                    Label[HID][hubID]=disvally;
                    /// No need of dispeak
//                    peakPair=DisQueryPeak2(HID,hubID,Label);
//                    dispeak=peakPair.first; peakhub=peakPair.second;
//                    if(dispeak>disvally){
//                        Label[HID][hubID]=disvally;//may not be correct at this moment, e.g., larger than correct value
//
//                    }
//                    else{
////                            Label[HID][hubID]=dispeak;
//                        Label[HID].erase(hubID);
//
//                        PruningPointNew[HID][peakhub].insert(hubID);
//                        PruningPointNew[hubID][peakhub].insert(HID);
//                    }
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
//            vector<vector<int>> ChangePTem;
//            ChangePTem.assign(nodenum,vecint);
            vector<vector<pair<int,int>>> ChangePTem;
            ChangePTem.assign(nodenum,vector<pair<int,int>>());


            //Change->Change & ChangeP
            for(set<int>::iterator it=WaitPro.begin();it!=WaitPro.end();it++){
                int curID=*it;
                vector<pair<int,int>> curChange=Change[curID];
                int neiID, neiDis, hID, hDis;

                for(int j=0;j<curChange.size();j++){
                    hID=curChange[j].first; hDis=curChange[j].second;
                    //Change->Change
                    for(int k=0;k<Neighbors[curID].size();k++){
                        neiID=Neighbors[curID][k].first; neiDis=Neighbors[curID][k].second;

//                        if(neiID == 99110 && hID == 212434){
//                            cout<<"AL1->AL1: "<<neiID<<" "<<hID<<endl;
//                        }
                        if(Label[neiID].find(hID)!=Label[neiID].end() && neiDis+hDis==Label[neiID][hID]){
                            disvally=DisQueryVally(neiID,hID,Neighbors,Label);
                            if(Label[neiID][hID]<disvally){
                                WaitProTem.insert(neiID);
                                ChangeTem[neiID].push_back(make_pair(hID, neiDis+hDis));
                                Label[neiID][hID]=disvally;
                                /// No need of dispeak
//                                peakPair=DisQueryPeak2(neiID,hID,Label);
//                                dispeak=peakPair.first; peakhub=peakPair.second;
//                                if(dispeak>disvally){
//                                    Label[neiID][hID]=disvally;//may not be correct at this moment, e.g., larger than correct value
//
//                                }
//                                else{
////                                        Label[neiID][hID]=dispeak;
//                                    Label[neiID].erase(hID);
//
//                                    PruningPointNew[neiID][peakhub].insert(hID);
//                                    PruningPointNew[hID][peakhub].insert(neiID);
//
//                                }
                            }
                        }
                    }

                    //Change->ChangeP
                    if(PruningPointNew[curID].find(hID)!=PruningPointNew[curID].end()){
                        outdatedPruning.clear();
                        for(auto snum=PruningPointNew[curID][hID].begin();snum!=PruningPointNew[curID][hID].end();++snum){//for each pruned vertex
                            int s=*snum;

//                            if(NodeOrder[s]<NodeOrder[curID] && NoSupportedPair.find(make_pair(s,curID))==NoSupportedPair.end()){
                            if(NodeOrder[s]<NodeOrder[curID]){

                                disvally=DisQueryVally(s,curID,Neighbors,Label);
//                                dispeak=DisQueryPeak(s,curID,Label);
                                peakPair=DisQueryPeak2(s,curID,Label);//the original dispeak is smaller than disvally
                                dispeak=peakPair.first; peakhub=peakPair.second;
                                if(dispeak>disvally){
                                    WaitProPTem.insert(s);
//                                    ChangePTem[s].push_back(curID);
                                    ChangePTem[s].push_back(make_pair(curID,dispeak));
                                    Label[s][curID]=disvally;

                                    outdatedPruning.insert(make_tuple(curID,hID,s));
                                    outdatedPruning.insert(make_tuple(s,hID,curID));

//                                    NoSupportedPair.insert(make_pair(s,curID));
                                }
                                else {//if dispeak<=disvally
                                    if(Label[s].find(curID) != Label[s].end()) {
                                        Label[s].erase(curID);
//                                        if (dispeak != Label[s][curID]) {
//                                            Label[s][curID] = dispeak;
//                                        }
                                    }

                                    if(peakhub != -1 && peakhub != hID){
                                        outdatedPruning.insert(make_tuple(curID,hID,s));
                                        outdatedPruning.insert(make_tuple(s,hID,curID));

                                        PruningPointNew[curID][peakhub].push_back(s);
                                        PruningPointNew[s][peakhub].push_back(curID);

                                    }
                                }
//                            }else if(NodeOrder[s]>NodeOrder[curID] && NoSupportedPair.find(make_pair(curID,s))==NoSupportedPair.end()){
                            } else if(NodeOrder[s]>NodeOrder[curID]){

                                disvally=DisQueryVally(curID,s,Neighbors,Label);
//                                dispeak=DisQueryPeak(curID,s,Label);
                                peakPair=DisQueryPeak2(curID,s,Label);//the old P2 distance dominates P1 distance, but the increase of d(curID,hID)
                                dispeak=peakPair.first; peakhub=peakPair.second;
                                if(dispeak>disvally){
                                    WaitProPTem.insert(curID);
//                                    ChangePTem[curID].push_back(s);
                                    ChangePTem[curID].push_back(make_pair(s,dispeak));
                                    Label[curID][s]=disvally;
                                    outdatedPruning.insert(make_tuple(curID,hID,s));
                                    outdatedPruning.insert(make_tuple(s,hID,curID));

//                                    NoSupportedPair.insert(make_pair(curID,s));
                                }
                                else {//if dispeak<=disvally
                                    if(Label[curID].find(s) != Label[curID].end()) {
                                        Label[curID].erase(s);
//                                        if (dispeak != Label[curID][s]) {
//                                            Label[curID][s] = dispeak;//should be correct value
//                                        }
                                    }

                                    if(peakhub != -1 && peakhub != hID) {
                                        outdatedPruning.insert(make_tuple(curID, hID, s));
                                        outdatedPruning.insert(make_tuple(s, hID, curID));

                                        PruningPointNew[curID][peakhub].push_back(s);
                                        PruningPointNew[s][peakhub].push_back(curID);
                                    }
                                }
                            }
                        }
                        for(auto it2=outdatedPruning.begin();it2!=outdatedPruning.end();++it2){
                            if(PruningPointNew[get<0>(*it2)].find(get<1>(*it2)) != PruningPointNew[get<0>(*it2)].end()) {//if found
                                vector<list<int>::iterator> temp_its;
                                for(auto it=PruningPointNew[get<0>(*it2)][get<1>(*it2)].begin();it!=PruningPointNew[get<0>(*it2)][get<1>(*it2)].end();++it){
                                    if(*it==get<2>(*it2)){
                                        temp_its.push_back(it);
                                    }
                                }
                                for(int i=0;i<temp_its.size();++i){
                                    PruningPointNew[get<0>(*it2)][get<1>(*it2)].erase(temp_its[i]);
                                }
                            }

                        }
                    }
                }
            }

            //ChangeP->CHangeP
            int v,u,neiid,neiw,cDis;
            for(set<int>::iterator itp=WaitProP.begin();itp!=WaitProP.end();itp++){
                v=*itp;
                for(int k=0;k<ChangeP[v].size();k++){
//                    u=ChangeP[v][k];
                    u=ChangeP[v][k].first, cDis=ChangeP[v][k].second;
                    outdatedPruning.clear();
                    for(int l=0;l<Neighbors[v].size();l++){
                        neiid=Neighbors[v][l].first; neiw=Neighbors[v][l].second;
//                        if(Label[neiid].find(u)!=Label[neiid].end() && neiw+cDis==Label[neiid][u]){
//
//                        }
                        if(NodeOrder[neiid]<NodeOrder[u]){
                            disvally=DisQueryVally(neiid, u,Neighbors,Label);
//                            disvally=DisQueryVallyVert(neiid, u,Neighbors,Label);
//                            dispeak=DisQueryPeak(neiid, u,Label);
                            peakPair=DisQueryPeak2(neiid, u,Label);
                            dispeak=peakPair.first; peakhub=peakPair.second;

//                            if(neiid == 99110 && u == 212434){
//                                cout<<"AL2->AL2: "<<neiid<<" "<<u<<" "<<disvally<<" "<<dispeak;
//                                if(Label[neiid].find(u)!=Label[neiid].end()){
//                                    cout<<" "<<Label[neiid][u];
//                                }
//                                cout<<endl;
////                                DisQueryVallyDebug(neiid, u,Neighbors,Label);
//                            }
                            if(disvally<dispeak){
                                if(Label[neiid].find(u)==Label[neiid].end() || (Label[neiid].find(u)!=Label[neiid].end() && Label[neiid][u]!=disvally)){///
                                    WaitProPTem.insert(neiid);
//                                    ChangePTem[neiid].push_back(u);
                                    ChangePTem[neiid].push_back(make_pair(u,dispeak));
                                    Label[neiid][u]=disvally;
//                                    outdatedPruning.insert(make_tuple(neiid,hubID,u));//
//                                    outdatedPruning.insert(make_tuple(u,hubID,neiid));//

//                                    NoSupportedPair.insert(make_pair(neiid,u));
                                }
                            }
                            else {//if dispeak<=disvally
//                                if(Label[neiid].find(u) != Label[neiid].end()) {///wrong part
//                                    Label[neiid].erase(u);
//                                }

                                if(peakhub != -1 && peakhub != hubID){
//                                    outdatedPruning.insert(make_tuple(neiid,hubID,hID));//
//                                    outdatedPruning.insert(make_tuple(hID,hubID,neiid));//

                                    PruningPointNew[neiid][peakhub].push_back(u);
                                    PruningPointNew[u][peakhub].push_back(neiid);
                                }
                            }
                        }
                    }
                    for(auto it2=outdatedPruning.begin();it2!=outdatedPruning.end();++it2){
                        if(PruningPointNew[get<0>(*it2)].find(get<1>(*it2)) != PruningPointNew[get<0>(*it2)].end()){//if found
                            vector<list<int>::iterator> temp_its;
                            for(auto it=PruningPointNew[get<0>(*it2)][get<1>(*it2)].begin();it!=PruningPointNew[get<0>(*it2)][get<1>(*it2)].end();++it){
                                if(*it==get<2>(*it2)){
                                    temp_its.push_back(it);
                                }
                            }
                            for(int i=0;i<temp_its.size();++i){
                                PruningPointNew[get<0>(*it2)][get<1>(*it2)].erase(temp_its[i]);
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

//new version: set version without NoSupportedPair, correct
/*void Graph::IncreasePSL(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label,vector<unordered_map<int,set<int>>> &PruningPointNew){
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

    int LID,HID;
    if(NodeOrder[a]>NodeOrder[b]){
        LID=b; HID=a;
    }else{
        LID=a; HID=b;
    }

    int dis,disvally,dispeak,peakhub;
    pair<int,int> peakPair;//distance, hubID
    set<tuple<int,int,int>> outdatedPruning;//<nodeID,supportNode,prunedID>
    //activate or not
    dis=DisQueryLower1(LID,HID,Neighbors,Label);

    //if(dispeak<=oldW) return;
    if(Label[LID].find(HID)!=Label[LID].end() && dis>oldW){//index update triggered
        vector<vector<pair<int,int>>> Change;
        vector<pair<int,int>> vec;
        Change.assign(nodenum,vec);
        set<int> WaitPro;
//        vector<vector<int>> ChangeP;
        vector<vector<pair<int,int>>> ChangeP;
//        vector<int> vecint;
//        ChangeP.assign(nodenum,vecint);
        ChangeP.assign(nodenum,vector<pair<int,int>>());
        set<int> WaitProP;

        WaitPro.insert(LID);
        Change[LID].push_back(make_pair(HID, oldW));
        disvally=DisQueryVally(LID,HID,Neighbors,Label);
        Label[LID][HID]=disvally;//correct to the new value

        //affected by the w(a.b)
        int hubID, hDis;
        int dis, cnt;
        for(unordered_map<int,int>::iterator it=Label[HID].begin();it!=Label[HID].end();it++){
            hubID=(*it).first; hDis=(*it).second;
            if(hubID!=HID && Label[LID].find(hubID)!=Label[LID].end() && oldW+hDis==Label[LID][hubID]){
                disvally=DisQueryVally(LID,hubID,Neighbors,Label);
                if(Label[LID][hubID]<disvally){
                    WaitPro.insert(LID);
                    Change[LID].push_back(make_pair(hubID, oldW+hDis));
                    Label[LID][hubID]=disvally;
                    /// No need of dispeak
//                    peakPair=DisQueryPeak2(LID,hubID, Label);
//                    dispeak=peakPair.first; peakhub=peakPair.second;
//                    if(dispeak>disvally){
//                        Label[LID][hubID]=disvally;//may not be correct at this moment
//                    }
//                    else{
////                            Label[LID][hubID]=dispeak;//may not be correct at this moment
//                        Label[LID].erase(hubID);
//
//                        PruningPointNew[LID][peakhub].insert(hubID);
//                        PruningPointNew[hubID][peakhub].insert(LID);
//                    }
                }
            }
        }

        for(unordered_map<int,int>::iterator it=Label[LID].begin();it!=Label[LID].end();it++){
            hubID=(*it).first; hDis=(*it).second;
            if(Label[HID].find(hubID)!=Label[HID].end() && oldW+hDis==Label[HID][hubID]){
                disvally=DisQueryVally(HID,hubID,Neighbors,Label);
                if(Label[HID][hubID]<disvally){
                    WaitPro.insert(HID);
                    Change[HID].push_back(make_pair(hubID, oldW+hDis));
                    Label[HID][hubID]=disvally;
                    /// No need of dispeak
//                    peakPair=DisQueryPeak2(HID,hubID,Label);
//                    dispeak=peakPair.first; peakhub=peakPair.second;
//                    if(dispeak>disvally){
//                        Label[HID][hubID]=disvally;//may not be correct at this moment, e.g., larger than correct value
//
//                    }
//                    else{
////                            Label[HID][hubID]=dispeak;
//                        Label[HID].erase(hubID);
//
//                        PruningPointNew[HID][peakhub].insert(hubID);
//                        PruningPointNew[hubID][peakhub].insert(HID);
//                    }
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
//            vector<vector<int>> ChangePTem;
//            ChangePTem.assign(nodenum,vecint);
            vector<vector<pair<int,int>>> ChangePTem;
            ChangePTem.assign(nodenum,vector<pair<int,int>>());


            //Change->Change & ChangeP
            for(set<int>::iterator it=WaitPro.begin();it!=WaitPro.end();it++){
                int curID=*it;
                vector<pair<int,int>> curChange=Change[curID];
                int neiID, neiDis, hID, hDis;

                for(int j=0;j<curChange.size();j++){
                    hID=curChange[j].first; hDis=curChange[j].second;
                    //Change->Change
                    for(int k=0;k<Neighbors[curID].size();k++){
                        neiID=Neighbors[curID][k].first; neiDis=Neighbors[curID][k].second;

//                        if(neiID == 99110 && hID == 212434){
//                            cout<<"AL1->AL1: "<<neiID<<" "<<hID<<endl;
//                        }
                        if(Label[neiID].find(hID)!=Label[neiID].end() && neiDis+hDis==Label[neiID][hID]){
                            disvally=DisQueryVally(neiID,hID,Neighbors,Label);
                            if(Label[neiID][hID]<disvally){
                                WaitProTem.insert(neiID);
                                ChangeTem[neiID].push_back(make_pair(hID, neiDis+hDis));
                                Label[neiID][hID]=disvally;
                                /// No need of dispeak
//                                peakPair=DisQueryPeak2(neiID,hID,Label);
//                                dispeak=peakPair.first; peakhub=peakPair.second;
//                                if(dispeak>disvally){
//                                    Label[neiID][hID]=disvally;//may not be correct at this moment, e.g., larger than correct value
//
//                                }
//                                else{
////                                        Label[neiID][hID]=dispeak;
//                                    Label[neiID].erase(hID);
//
//                                    PruningPointNew[neiID][peakhub].insert(hID);
//                                    PruningPointNew[hID][peakhub].insert(neiID);
//
//                                }
                            }
                        }
                    }

                    //Change->ChangeP
                    if(PruningPointNew[curID].find(hID)!=PruningPointNew[curID].end()){
                        outdatedPruning.clear();
                        for(auto snum=PruningPointNew[curID][hID].begin();snum!=PruningPointNew[curID][hID].end();++snum){//for each pruned vertex
                            int s=*snum;
//                        for(int sum=0;sum<PruningPointNew[curID][hID].size();sum++){
//                            int s=PruningPointNew[curID][hID][sum];

//                            if(NodeOrder[s]<NodeOrder[curID] && NoSupportedPair.find(make_pair(s,curID))==NoSupportedPair.end()){
                            if(NodeOrder[s]<NodeOrder[curID]){
//                                if(s == 99110 && curID == 212434){
//                                    cout<<"AL1->AL2: "<<s<<" "<<curID<<endl;
//                                }

                                disvally=DisQueryVally(s,curID,Neighbors,Label);
//                                dispeak=DisQueryPeak(s,curID,Label);
                                peakPair=DisQueryPeak2(s,curID,Label);//the original dispeak is smaller than disvally
                                dispeak=peakPair.first; peakhub=peakPair.second;
                                if(dispeak>disvally){
                                    WaitProPTem.insert(s);
//                                    ChangePTem[s].push_back(curID);
                                    ChangePTem[s].push_back(make_pair(curID,dispeak));
                                    Label[s][curID]=disvally;

                                    outdatedPruning.insert(make_tuple(curID,hID,s));
                                    outdatedPruning.insert(make_tuple(s,hID,curID));

//                                    NoSupportedPair.insert(make_pair(s,curID));
                                }
                                else {//if dispeak<=disvally
                                    if(Label[s].find(curID) != Label[s].end()) {
                                        Label[s].erase(curID);
//                                        if (dispeak != Label[s][curID]) {
//                                            Label[s][curID] = dispeak;
//                                        }
                                    }

                                    if(peakhub != -1 && peakhub != hID){
                                        assert(PruningPointNew[curID].find(hID) != PruningPointNew[curID].end() && PruningPointNew[curID][hID].find(s) != PruningPointNew[curID][hID].end());
                                        assert(PruningPointNew[s].find(hID) != PruningPointNew[s].end() && PruningPointNew[s][hID].find(curID) != PruningPointNew[s][hID].end());
                                        outdatedPruning.insert(make_tuple(curID,hID,s));
                                        outdatedPruning.insert(make_tuple(s,hID,curID));

                                        PruningPointNew[curID][peakhub].insert(s);
                                        PruningPointNew[s][peakhub].insert(curID);

                                    }
                                }
//                            }else if(NodeOrder[s]>NodeOrder[curID] && NoSupportedPair.find(make_pair(curID,s))==NoSupportedPair.end()){
                            } else if(NodeOrder[s]>NodeOrder[curID]){
//                                if(curID == 99110 && s == 212434){
//                                    cout<<"AL1->AL2: "<<s<<" "<<curID<<endl;
//                                }

                                disvally=DisQueryVally(curID,s,Neighbors,Label);
//                                dispeak=DisQueryPeak(curID,s,Label);
                                peakPair=DisQueryPeak2(curID,s,Label);//the old P2 distance dominates P1 distance, but the increase of d(curID,hID)
                                dispeak=peakPair.first; peakhub=peakPair.second;
                                if(dispeak>disvally){
                                    WaitProPTem.insert(curID);
//                                    ChangePTem[curID].push_back(s);
                                    ChangePTem[curID].push_back(make_pair(s,dispeak));
                                    Label[curID][s]=disvally;
                                    outdatedPruning.insert(make_tuple(curID,hID,s));
                                    outdatedPruning.insert(make_tuple(s,hID,curID));

//                                    NoSupportedPair.insert(make_pair(curID,s));
                                }
                                else {//if dispeak<=disvally
                                    if(Label[curID].find(s) != Label[curID].end()) {
                                        Label[curID].erase(s);
//                                        if (dispeak != Label[curID][s]) {
//                                            Label[curID][s] = dispeak;//should be correct value
//                                        }
                                    }

                                    if(peakhub != -1 && peakhub != hID) {
                                        assert(PruningPointNew[curID].find(hID) != PruningPointNew[curID].end() && PruningPointNew[curID][hID].find(s) != PruningPointNew[curID][hID].end());
                                        assert(PruningPointNew[s].find(hID) != PruningPointNew[s].end() && PruningPointNew[s][hID].find(curID) != PruningPointNew[s][hID].end());
                                        outdatedPruning.insert(make_tuple(curID, hID, s));
                                        outdatedPruning.insert(make_tuple(s, hID, curID));

                                        PruningPointNew[curID][peakhub].insert(s);
                                        PruningPointNew[s][peakhub].insert(curID);
                                    }
                                }
                            }
                        }
                        for(auto it2=outdatedPruning.begin();it2!=outdatedPruning.end();++it2){
                            if(PruningPointNew[get<0>(*it2)].find(get<1>(*it2)) != PruningPointNew[get<0>(*it2)].end()) {//if found
                                PruningPointNew[get<0>(*it2)][get<1>(*it2)].erase(get<2>(*it2));
                            }

                        }
                    }
                }
            }

            //ChangeP->CHangeP
            int v,u,neiid,neiw,cDis;
            for(set<int>::iterator itp=WaitProP.begin();itp!=WaitProP.end();itp++){
                v=*itp;
                for(int k=0;k<ChangeP[v].size();k++){
//                    u=ChangeP[v][k];
                    u=ChangeP[v][k].first, cDis=ChangeP[v][k].second;
                    outdatedPruning.clear();
                    for(int l=0;l<Neighbors[v].size();l++){
                        neiid=Neighbors[v][l].first; neiw=Neighbors[v][l].second;
//                        if(Label[neiid].find(u)!=Label[neiid].end() && neiw+cDis==Label[neiid][u]){
//
//                        }
                        if(NodeOrder[neiid]<NodeOrder[u]){
                            disvally=DisQueryVally(neiid, u,Neighbors,Label);
//                            disvally=DisQueryVallyVert(neiid, u,Neighbors,Label);
//                            dispeak=DisQueryPeak(neiid, u,Label);
                            peakPair=DisQueryPeak2(neiid, u,Label);
                            dispeak=peakPair.first; peakhub=peakPair.second;

//                            if(neiid == 99110 && u == 212434){
//                                cout<<"AL2->AL2: "<<neiid<<" "<<u<<" "<<disvally<<" "<<dispeak;
//                                if(Label[neiid].find(u)!=Label[neiid].end()){
//                                    cout<<" "<<Label[neiid][u];
//                                }
//                                cout<<endl;
////                                DisQueryVallyDebug(neiid, u,Neighbors,Label);
//                            }
                            if(disvally<dispeak){
                                if(Label[neiid].find(u)==Label[neiid].end() || (Label[neiid].find(u)!=Label[neiid].end() && Label[neiid][u]>disvally)){///
                                    WaitProPTem.insert(neiid);
//                                    ChangePTem[neiid].push_back(u);
                                    ChangePTem[neiid].push_back(make_pair(u,dispeak));
                                    Label[neiid][u]=disvally;
//                                    outdatedPruning.insert(make_tuple(neiid,hubID,u));//
//                                    outdatedPruning.insert(make_tuple(u,hubID,neiid));//

//                                    NoSupportedPair.insert(make_pair(neiid,u));
                                }
                            }
                            else {//if dispeak<=disvally
//                                if(Label[neiid].find(u) != Label[neiid].end()) {///wrong part
//                                    Label[neiid].erase(u);
//                                }

                                if(peakhub != -1 && peakhub != hubID){
//                                    outdatedPruning.insert(make_tuple(neiid,hubID,hID));//
//                                    outdatedPruning.insert(make_tuple(hID,hubID,neiid));//

                                    PruningPointNew[neiid][peakhub].insert(u);
                                    PruningPointNew[u][peakhub].insert(neiid);
                                }
                            }
                        }
                    }
                    for(auto it2=outdatedPruning.begin();it2!=outdatedPruning.end();++it2){
                        if(PruningPointNew[get<0>(*it2)].find(get<1>(*it2)) != PruningPointNew[get<0>(*it2)].end()){//if found
                            PruningPointNew[get<0>(*it2)][get<1>(*it2)].erase(get<2>(*it2));
//                        PruningPoint[get<1>(*it)][get<0>(*it)].erase(get<2>(*it));
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

//original version: correct version
/*void Graph::IncreasePSL(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label,vector<unordered_map<int,set<int>>> &PruningPointNew){
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

    int LID,HID;
    if(NodeOrder[a]>NodeOrder[b]){
        LID=b; HID=a;
    }else{
        LID=a; HID=b;
    }

    int dis,disvally,dispeak,peakhub;
    pair<int,int> peakPair;//distance, hubID
    set<tuple<int,int,int>> outdatedPruning;//<nodeID,supportNode,prunedID>
    //activate or not
    dis=DisQueryLower1(LID,HID,Neighbors,Label);

    //if(dispeak<=oldW) return;
    if(Label[LID].find(HID)!=Label[LID].end() && dis>oldW){//index update triggered
        vector<vector<pair<int,int>>> Change;
        vector<pair<int,int>> vec;
        Change.assign(nodenum,vec);
        set<int> WaitPro;
//        vector<vector<int>> ChangeP;
        vector<vector<pair<int,int>>> ChangeP;
//        vector<int> vecint;
//        ChangeP.assign(nodenum,vecint);
        ChangeP.assign(nodenum,vector<pair<int,int>>());
        set<int> WaitProP;

        WaitPro.insert(LID);
        Change[LID].push_back(make_pair(HID, oldW));
        disvally=DisQueryVally(LID,HID,Neighbors,Label);
        Label[LID][HID]=disvally;//correct to the new value

        //affected by the w(a.b)
        int hubID, hDis;
        int dis, cnt;
        for(unordered_map<int,int>::iterator it=Label[HID].begin();it!=Label[HID].end();it++){
            hubID=(*it).first; hDis=(*it).second;
            if(hubID!=HID && Label[LID].find(hubID)!=Label[LID].end() && oldW+hDis==Label[LID][hubID]){
                disvally=DisQueryVally(LID,hubID,Neighbors,Label);
                if(Label[LID][hubID]<disvally){
                    WaitPro.insert(LID);
                    Change[LID].push_back(make_pair(hubID, oldW+hDis));
//                    Label[LID][hubID]=disvally;

                    peakPair=DisQueryPeak2(LID,hubID, Label);
                    dispeak=peakPair.first; peakhub=peakPair.second;
                    if(dispeak>disvally){
                        Label[LID][hubID]=disvally;//may not be correct at this moment
                    }
                    else{
//                            Label[LID][hubID]=dispeak;//may not be correct at this moment
                        Label[LID].erase(hubID);

                        PruningPointNew[LID][peakhub].insert(hubID);
                        PruningPointNew[hubID][peakhub].insert(LID);
                    }
                }
            }
        }

        for(unordered_map<int,int>::iterator it=Label[LID].begin();it!=Label[LID].end();it++){
            hubID=(*it).first; hDis=(*it).second;
            if(Label[HID].find(hubID)!=Label[HID].end() && oldW+hDis==Label[HID][hubID]){
                disvally=DisQueryVally(HID,hubID,Neighbors,Label);
                if(Label[HID][hubID]<disvally){
                    WaitPro.insert(HID);
                    Change[HID].push_back(make_pair(hubID, oldW+hDis));
//                    Label[HID][hubID]=disvally;

                    peakPair=DisQueryPeak2(HID,hubID,Label);
                    dispeak=peakPair.first; peakhub=peakPair.second;
                    if(dispeak>disvally){
                        Label[HID][hubID]=disvally;//may not be correct at this moment, e.g., larger than correct value

                    }
                    else{
//                            Label[HID][hubID]=dispeak;
                        Label[HID].erase(hubID);

                        PruningPointNew[HID][peakhub].insert(hubID);
                        PruningPointNew[hubID][peakhub].insert(HID);
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
//            vector<vector<int>> ChangePTem;
//            ChangePTem.assign(nodenum,vecint);
            vector<vector<pair<int,int>>> ChangePTem;
            ChangePTem.assign(nodenum,vector<pair<int,int>>());


            //Change->Change & ChangeP
            for(set<int>::iterator it=WaitPro.begin();it!=WaitPro.end();it++){
                int curID=*it;
                vector<pair<int,int>> curChange=Change[curID];
                int neiID, neiDis, hID, hDis;

                for(int j=0;j<curChange.size();j++){
                    hID=curChange[j].first; hDis=curChange[j].second;
                    //Change->Change
                    for(int k=0;k<Neighbors[curID].size();k++){
                        neiID=Neighbors[curID][k].first; neiDis=Neighbors[curID][k].second;

                        if(neiID == 99110 && hID == 212434){
                            cout<<"AL1->AL1: "<<neiID<<" "<<hID<<endl;
                        }
                        if(Label[neiID].find(hID)!=Label[neiID].end() && neiDis+hDis==Label[neiID][hID]){
                            disvally=DisQueryVally(neiID,hID,Neighbors,Label);
                            if(Label[neiID][hID]<disvally){
                                WaitProTem.insert(neiID);
                                ChangeTem[neiID].push_back(make_pair(hID, neiDis+hDis));
//                                Label[neiID][hID]=disvally;

                                peakPair=DisQueryPeak2(neiID,hID,Label);
                                dispeak=peakPair.first; peakhub=peakPair.second;
                                if(dispeak>disvally){
                                    Label[neiID][hID]=disvally;//may not be correct at this moment, e.g., larger than correct value

                                }
                                else{
//                                        Label[neiID][hID]=dispeak;
                                    Label[neiID].erase(hID);

                                    PruningPointNew[neiID][peakhub].insert(hID);
                                    PruningPointNew[hID][peakhub].insert(neiID);

                                }
                            }
                        }
                    }

                    //Change->ChangeP
                    if(PruningPointNew[curID].find(hID)!=PruningPointNew[curID].end()){
                        outdatedPruning.clear();
                        for(auto snum=PruningPointNew[curID][hID].begin();snum!=PruningPointNew[curID][hID].end();++snum){//for each pruned vertex
                            int s=*snum;
//                        for(int sum=0;sum<PruningPointNew[curID][hID].size();sum++){
//                            int s=PruningPointNew[curID][hID][sum];


//                            if(NodeOrder[s]<NodeOrder[curID] && NoSupportedPair.find(make_pair(s,curID))==NoSupportedPair.end()){
                            if(NodeOrder[s]<NodeOrder[curID]){
                                if(s == 99110 && curID == 212434){
                                    cout<<"AL1->AL2: "<<s<<" "<<curID<<endl;
                                }

                                disvally=DisQueryVally(s,curID,Neighbors,Label);
//                                dispeak=DisQueryPeak(s,curID,Label);
                                peakPair=DisQueryPeak2(s,curID,Label);//the original dispeak is smaller than disvally
                                dispeak=peakPair.first; peakhub=peakPair.second;
                                if(dispeak>disvally){
                                    WaitProPTem.insert(s);
//                                    ChangePTem[s].push_back(curID);
                                    ChangePTem[s].push_back(make_pair(curID,dispeak));
                                    Label[s][curID]=disvally;

                                    outdatedPruning.insert(make_tuple(curID,hID,s));
                                    outdatedPruning.insert(make_tuple(s,hID,curID));

//                                    NoSupportedPair.insert(make_pair(s,curID));
                                }
                                else {//if dispeak<=disvally
                                    if(Label[s].find(curID) != Label[s].end()) {
                                        Label[s].erase(curID);
//                                        if (dispeak != Label[s][curID]) {
//                                            Label[s][curID] = dispeak;
//                                        }
                                    }

                                    if(peakhub != -1 && peakhub != hID){
                                        assert(PruningPointNew[curID].find(hID) != PruningPointNew[curID].end() && PruningPointNew[curID][hID].find(s) != PruningPointNew[curID][hID].end());
                                        assert(PruningPointNew[s].find(hID) != PruningPointNew[s].end() && PruningPointNew[s][hID].find(curID) != PruningPointNew[s][hID].end());
                                        outdatedPruning.insert(make_tuple(curID,hID,s));
                                        outdatedPruning.insert(make_tuple(s,hID,curID));

                                        PruningPointNew[curID][peakhub].insert(s);
                                        PruningPointNew[s][peakhub].insert(curID);

                                    }
                                }
//                            }else if(NodeOrder[s]>NodeOrder[curID] && NoSupportedPair.find(make_pair(curID,s))==NoSupportedPair.end()){
                            }
                            else if(NodeOrder[s]>NodeOrder[curID]){
                                if(curID == 99110 && s == 212434){
                                    cout<<"AL1->AL2: "<<s<<" "<<curID<<endl;
                                }

                                disvally=DisQueryVally(curID,s,Neighbors,Label);
//                                dispeak=DisQueryPeak(curID,s,Label);
                                peakPair=DisQueryPeak2(curID,s,Label);//the old P2 distance dominates P1 distance, but the increase of d(curID,hID)
                                dispeak=peakPair.first; peakhub=peakPair.second;
                                if(dispeak>disvally){
                                    WaitProPTem.insert(curID);
//                                    ChangePTem[curID].push_back(s);
                                    ChangePTem[curID].push_back(make_pair(s,dispeak));
                                    Label[curID][s]=disvally;
                                    outdatedPruning.insert(make_tuple(curID,hID,s));
                                    outdatedPruning.insert(make_tuple(s,hID,curID));

//                                    NoSupportedPair.insert(make_pair(curID,s));
                                }
                                else {//if dispeak<=disvally
                                    if(Label[curID].find(s) != Label[curID].end()) {
                                        Label[curID].erase(s);
//                                        if (dispeak != Label[curID][s]) {
//                                            Label[curID][s] = dispeak;//should be correct value
//                                        }
                                    }

                                    if(peakhub != -1 && peakhub != hID) {
                                        assert(PruningPointNew[curID].find(hID) != PruningPointNew[curID].end() && PruningPointNew[curID][hID].find(s) != PruningPointNew[curID][hID].end());
                                        assert(PruningPointNew[s].find(hID) != PruningPointNew[s].end() && PruningPointNew[s][hID].find(curID) != PruningPointNew[s][hID].end());
                                        outdatedPruning.insert(make_tuple(curID, hID, s));
                                        outdatedPruning.insert(make_tuple(s, hID, curID));

                                        PruningPointNew[curID][peakhub].insert(s);
                                        PruningPointNew[s][peakhub].insert(curID);
                                    }
                                }
                            }
                        }
                        for(auto it2=outdatedPruning.begin();it2!=outdatedPruning.end();++it2){
                            if(PruningPointNew[get<0>(*it2)].find(get<1>(*it2)) != PruningPointNew[get<0>(*it2)].end()) {//if found
                                PruningPointNew[get<0>(*it2)][get<1>(*it2)].erase(get<2>(*it2));
                            }

                        }
                    }
                }
            }

            //ChangeP->CHangeP
            int v,u,neiid,neiw,cDis;
            for(set<int>::iterator itp=WaitProP.begin();itp!=WaitProP.end();itp++){
                v=*itp;
                for(int k=0;k<ChangeP[v].size();k++){
//                    u=ChangeP[v][k];
                    u=ChangeP[v][k].first, cDis=ChangeP[v][k].second;
                    outdatedPruning.clear();
                    for(int l=0;l<Neighbors[v].size();l++){
                        neiid=Neighbors[v][l].first; neiw=Neighbors[v][l].second;
//                        if(Label[neiid].find(u)!=Label[neiid].end() && neiw+cDis==Label[neiid][u]){
//
//                        }
                        if(NodeOrder[neiid]<NodeOrder[u]){
                            disvally=DisQueryVally(neiid, u,Neighbors,Label);
//                            disvally=DisQueryVallyVert(neiid, u,Neighbors,Label);
//                            dispeak=DisQueryPeak(neiid, u,Label);
                            peakPair=DisQueryPeak2(neiid, u,Label);
                            dispeak=peakPair.first; peakhub=peakPair.second;

                            if(neiid == 99110 && u == 212434){
                                cout<<"AL2->AL2: "<<neiid<<" "<<u<<" "<<disvally<<" "<<dispeak;
                                if(Label[neiid].find(u)!=Label[neiid].end()){
                                    cout<<" "<<Label[neiid][u];
                                }
                                cout<<endl;
//                                DisQueryVallyDebug(neiid, u,Neighbors,Label);
                            }
                            if(disvally<dispeak){
                                if(Label[neiid].find(u)==Label[neiid].end() || (Label[neiid].find(u)!=Label[neiid].end() && Label[neiid][u]>disvally)){///
                                    WaitProPTem.insert(neiid);
//                                    ChangePTem[neiid].push_back(u);
                                    ChangePTem[neiid].push_back(make_pair(u,dispeak));
                                    Label[neiid][u]=disvally;
//                                    outdatedPruning.insert(make_tuple(neiid,hubID,u));//
//                                    outdatedPruning.insert(make_tuple(u,hubID,neiid));//

//                                    NoSupportedPair.insert(make_pair(neiid,u));
                                }
                            }
                            else {//if dispeak<=disvally
//                                if(Label[neiid].find(u) != Label[neiid].end()) {///wrong part
//                                    Label[neiid].erase(u);
//                                }

                                if(peakhub != -1 && peakhub != hubID){
//                                    outdatedPruning.insert(make_tuple(neiid,hubID,hID));//
//                                    outdatedPruning.insert(make_tuple(hID,hubID,neiid));//

                                    PruningPointNew[neiid][peakhub].insert(u);
                                    PruningPointNew[u][peakhub].insert(neiid);
                                }
                            }
                        }
                    }
                    for(auto it2=outdatedPruning.begin();it2!=outdatedPruning.end();++it2){
                        if(PruningPointNew[get<0>(*it2)].find(get<1>(*it2)) != PruningPointNew[get<0>(*it2)].end()){//if found
                            PruningPointNew[get<0>(*it2)][get<1>(*it2)].erase(get<2>(*it2));
//                        PruningPoint[get<1>(*it)][get<0>(*it)].erase(get<2>(*it));
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


//new new version by set
//void Graph::IncreasePSL(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label,vector<unordered_map<int,set<int>>> &PruningPointNew){
//    for(int i=0;i<Neighbors[a].size();i++){
//        if(Neighbors[a][i].first==b){
////            cout<<i<<" "<<a<<" "<<b<<" "<<Neighbors[a][i].second<<" "<<newW<<endl;
//            if(oldW != Neighbors[a][i].second){
//                cout<<"Old edge weight is incorrect! "<<a<<" "<<b<<": "<<oldW<<" "<<Neighbors[a][i].second<<endl;
//                oldW = Neighbors[a][i].second;
//            }
//            Neighbors[a][i].second=newW;
//            break;
//        }
//    }
//    for(int i=0;i<Neighbors[b].size();i++){
//        if(Neighbors[b][i].first==a){
////            cout<<i<<" "<<a<<" "<<b<<" "<<Neighbors[b][i].second<<" "<<newW<<endl;
//            if(oldW != Neighbors[b][i].second){
//                cout<<"Old edge weight is incorrect! "<<b<<" "<<a<<": "<<oldW<<" "<<Neighbors[b][i].second<<endl;
//                oldW = Neighbors[b][i].second;
//            }
//            Neighbors[b][i].second=newW;
//            break;
//        }
//    }
//
//    bool ifDebug = false;//false;
////    ifDebug = true;
//
//    int LID,HID;
//    if(NodeOrder[a]>NodeOrder[b]){
//        LID=b; HID=a;
//    }else{
//        LID=a; HID=b;
//    }
//
//    int tempNDis;
//    int dis,disvally,dispeak,peakhub;
//    pair<int,int> peakPair;//distance, hubID
//    set<tuple<int,int,int>> outdatedPruning;//<nodeID,supportNode,prunedID>
//    //activate or not
//    dis=DisQueryLower1(LID,HID,Neighbors,Label);//d1', via a neighbor of LID (which has lower order than HID) to reach HID
//    dispeak=DisQueryPeak(LID,HID,Label);//d2, via the common hub vertex of LID and HID
//    if(dispeak<=oldW)//if d2 is lower than oldW, the increase of oldW will not affect the shortest distance between LID and HID
//        return;
//    if(ifDebug){
//        if(dis == INF)
//            cout<<DisQueryLower1(LID,HID,Neighbors,Label);
//        if(dis==oldW){
//            cout<<dis<<" "<<oldW<<endl;
//        }
//    }
//
////    if(Label[196778].find(195124) != Label[196778].end()){
////        cout<<Label[196778][195124]<<" "<<DijkstraCore(196778,195124)<<endl;
////    }
//
//    if(Label[LID].find(HID)!=Label[LID].end() && dis>oldW){//if d1' is larger than oldW, it indicates the shortest distance between LID and HID is equal to oldW, index update triggered
//        vector<vector<pair<int,int>>> Change;//the label that is changed, (curID,<hID,oldDis>)
//        vector<pair<int,int>> vec;
//        Change.assign(nodenum,vec);
//        set<int> WaitPro;//the vertex waited for label propagation, i.e., AL1
////        vector<vector<int>> ChangeP;
//        vector<vector<pair<int,int>>> ChangeP;//the label that is changed, (curID,<hID,hubID>)
//        ChangeP.assign(nodenum,vector<pair<int,int>>());
//        set<int> WaitProP;//the vertex waited for label propagation, i.e., AL2
//
//
//        WaitPro.insert(LID);
//        Change[LID].push_back(make_pair(HID, oldW));
//        disvally=DisQueryVally(LID,HID,Neighbors,Label);//shortest distance query through the neighbor's label
////        disvally=DisQueryVallyVert(LID,HID,Neighbors,Label);//shortest distance query through the neighbor's label, correct one
//        Label[LID][HID]=disvally;//correct
//
//        //cout<<"start "<<LID<<" "<<HID<<" "<<disvally<<endl;
//
//        //affected by the w(a,b)
//        int hubID, hDis;
//        int dis, cnt;
//        for(auto it=Label[HID].begin();it!=Label[HID].end();++it){//check the vertex has higher order than HID and update LID's label (except HID)
//            hubID=(*it).first; hDis=(*it).second;
//            if(hubID!=HID && Label[LID].find(hubID)!=Label[LID].end()){// && oldW+hDis==Label[LID][hubID]
//                if(oldW+hDis==Label[LID][hubID]){//check the affected label, if the shortest path between LID and hubID pass through e(LID,HID)
//                    disvally=DisQueryVally(LID,hubID,Neighbors,Label);
////                    disvally=DisQueryVallyVert(LID,hubID,Neighbors,Label);
//                    if(Label[LID][hubID]<disvally){//if the new distance of P1 is higher, update it
//                        WaitPro.insert(LID);
//                        Change[LID].push_back(make_pair(hubID, oldW+hDis));//record the changed label of LID with the old distance
//                        Label[LID][hubID]=disvally;//may not be correct at this moment
////                        cout<<LID<<" "<<hubID<<" "<<oldW+hDis<<endl;
//
//
////                        peakPair=DisQueryPeak2(LID,hubID, Label);
////                        dispeak=peakPair.first; peakhub=peakPair.second;
////                        if(dispeak>disvally){
////                            Label[LID][hubID]=disvally;//may not be correct at this moment
////                            tempNDis = disvally;
////                        }
////                        else{
//////                            Label[LID][hubID]=dispeak;//may not be correct at this moment
////                            Label[LID].erase(hubID);
////                            tempNDis = dispeak;
////
////                            PruningPointNew[LID][peakhub].insert(hubID);
////                            PruningPointNew[hubID][peakhub].insert(LID);
////                        }
//
//                        //cout<<"weight "<<LID<<" "<<hubID<<" "<<disvally<<endl;
//
//                    }
//                }
//
//            }
//        }
//
//        for(auto it=Label[LID].begin();it!=Label[LID].end();++it){//update HID's label
//            hubID=(*it).first; hDis=(*it).second;
//            if(Label[HID].find(hubID)!=Label[HID].end()){//&& oldW+hDis==Label[HID][hubID]
//                if(oldW+hDis==Label[HID][hubID]){//check the affected label, if the shortest path between HID and hubID pass through e(LID,HID)
//                    disvally=DisQueryVally(HID,hubID,Neighbors,Label);
////                    disvally=DisQueryVallyVert(HID,hubID,Neighbors,Label);
//                    if(Label[HID][hubID]<disvally){
//                        WaitPro.insert(HID);
//                        Change[HID].push_back(make_pair(hubID, oldW+hDis));//record the changed label
//                        Label[HID][hubID]=disvally;//may not be correct at this moment, e.g., larger than correct value
////                        cout<<HID<<" "<<hubID<<" "<<oldW+hDis<<endl;
////                        peakPair=DisQueryPeak2Vert(HID,hubID,Neighbors,Label);
//
////                        peakPair=DisQueryPeak2(HID,hubID,Label);
////                        dispeak=peakPair.first; peakhub=peakPair.second;
////                        if(dispeak>disvally){
////                            Label[HID][hubID]=disvally;//may not be correct at this moment, e.g., larger than correct value
////                            tempNDis = disvally;
////                        }else{
//////                            Label[HID][hubID]=dispeak;
////                            Label[HID].erase(hubID);
////                            tempNDis=dispeak;
////
////                            PruningPointNew[HID][peakhub].insert(hubID);
////                            PruningPointNew[hubID][peakhub].insert(HID);
////                        }
//
//                        //cout<<"weight "<<HID<<" "<<hubID<<" "<<disvally<<endl;
//                    }
//                }
//            }
//        }
//
//        while(WaitProP.size()>0 || WaitPro.size()>0){//the vertex set that corresponding labels have been changed
//            set<int> WaitProTem;
//            vector<vector<pair<int,int>>> ChangeTem;
//            vector<pair<int,int>> vec;
//            ChangeTem.assign(nodenum,vec);
//            set<int> WaitProPTem;
//
//            vector<vector<pair<int,int>>> ChangePTem;
//            ChangePTem.assign(nodenum,vector<pair<int,int>>());
//
//            //Change->Change & ChangeP: AL1->AL1 and AL1->AL2
//            for(auto it=WaitPro.begin();it!=WaitPro.end();it++){
//                int curID=*it;
//                vector<pair<int,int>> curChange=Change[curID];//the changed label
//                int neiID, neiDis, hID, hDis;
//
//                for(int j=0;j<curChange.size();j++){
//                    hID=curChange[j].first; hDis=curChange[j].second;//the hub of changed label
//                    if(curID == 146346 || hID == 207156){
////                        cout<<"AL1: "<<curID<<" "<<hID<<endl;
//                    }
////                    if(Label[196778].find(195124) != Label[196778].end()){
////                        cout<<Label[196778][195124]<<" "<<DijkstraCore(196778,195124)<<" "<<curID<<" "<<hID<<endl;
////                    }
//                    //Change->Change: AL1->AL1
//                    ///check neighbors of curID
//                    for(int k=0;k<Neighbors[curID].size();k++){//identify the affected neighbors of curID
//                        neiID=Neighbors[curID][k].first; neiDis=Neighbors[curID][k].second;
//                        if(hID == 207156 && neiID == 146346){
//                            cout<<"neiID: "<<neiID<<" "<<hID<<endl;
//                        }
//                        if(Label[neiID].find(hID)!=Label[neiID].end() ){//&& neiDis+hDis==Label[neiID][hID]
//                            if(neiDis+hDis==Label[neiID][hID]){//check the affected label, if the shortest path between neiID and hID pass e(curID,hID)
////                                if(curID == 196731 && hID == 195124 && neiID == 196755){
////                                    cout<<curID<<"("<<NodeOrder[curID]<<") "<<hID<<"("<<NodeOrder[hID]<<") "<<neiID<<"("<<NodeOrder[neiID]<<")"<<endl;
////                                    cout<<"196778 "<<NodeOrder[196778]<<endl;
////                                    if(Label[curID].find(hID) != Label[curID].end()){
////                                        cout<<Label[curID][hID]<<endl;
////                                    }
////                                }
//                                disvally=DisQueryVally(neiID,hID,Neighbors,Label);
////                                disvally=DisQueryVallyVert2(neiID,hID,Neighbors,Label, hDis);
////                                assert(Label[curID].find(hID) != Label[curID].end());
////                                disvally=neiDis+Label[curID][hID];
//                                if(Label[neiID][hID]<disvally){
//                                    WaitProTem.insert(neiID);
//                                    ChangeTem[neiID].push_back(make_pair(hID, neiDis+hDis));
////                                    Label[neiID][hID]=disvally;//may not be the final correct value at this moment
//
////                                    peakPair=DisQueryPeak2Vert(neiID,hID,Neighbors,Label);
//                                    peakPair=DisQueryPeak2(neiID,hID,Label);
//                                    dispeak=peakPair.first; peakhub=peakPair.second;
//                                    if(dispeak>disvally){
//                                        Label[neiID][hID]=disvally;//may not be correct at this moment, e.g., larger than correct value
//                                        tempNDis=disvally;
//                                    }else{
////                                        Label[neiID][hID]=dispeak;
//                                        Label[neiID].erase(hID);
//                                        tempNDis=dispeak;
//
//                                        PruningPointNew[neiID][peakhub].insert(hID);
//                                        PruningPointNew[hID][peakhub].insert(neiID);
//
//                                    }
//
//                                    //cout<<"1--1 "<<neiID<<" "<<hID<<" "<<disvally<<endl;
//
//                                }
//                            }
//
//                        }
//                    }
//
//
//                    //Change->ChangeP: AL1->AL2
//                    if(PruningPointNew[curID].find(hID)!=PruningPointNew[curID].end()){//if hID has pruned some vertices, we need to check whether the pruned label should be added again as d(curID,hID) has changed
////                        for(int snum=0;snum<PruningPointNew[curID][hID].size();snum++){
//                        outdatedPruning.clear();
//                        for(auto snum=PruningPointNew[curID][hID].begin();snum!=PruningPointNew[curID][hID].end();++snum){//for each pruned vertex
//                            int s=*snum;
//                            //if the pruned shortest path pass the increased edge
////                            if(NodeOrder[s]<NodeOrder[curID] && NoSupportedPair.find(make_pair(s,curID))==NoSupportedPair.end()){//it is not in NoSupportedPair
//                            if(NodeOrder[s]<NodeOrder[curID]){//if s has lower order than curID, check whether curID should be added as a hub of s
//                                disvally=DisQueryVally(s,curID,Neighbors,Label);
//                                peakPair=DisQueryPeak2(s,curID,Label);//the original dispeak is smaller than disvally
////                                disvally=DisQueryVallyVert(s,curID,Neighbors,Label);
////                                peakPair=DisQueryPeakVert(s,curID,Neighbors,Label);
//                                int dispeakOld = 0;
//                                if(Label[s].find(hID) != Label[s].end()){
//                                    dispeakOld = hDis + Label[s][hID];
//                                }
//                                assert(dispeakOld <= disvally);
//                                dispeak=peakPair.first; peakhub=peakPair.second;
//                                if(dispeak>disvally){
//                                    WaitProPTem.insert(s);
//                                    ChangePTem[s].push_back(make_pair(curID,hID));
//                                    Label[s][curID]=disvally;
////                                    NoSupportedPair.insert(make_pair(s,curID));
//                                    outdatedPruning.insert(make_tuple(curID,hID,s));
//                                    outdatedPruning.insert(make_tuple(s,hID,curID));
//
//                                    //cout<<"1--2 "<<s<<" "<<curID<<" "<<disvally<<endl;
//                                }
//                                else {//if dispeak<=disvally
//                                    if(Label[s].find(curID) != Label[s].end()) {
//                                        Label[s].erase(curID);
////                                        if (dispeak != Label[s][curID]) {
////                                            Label[s][curID] = dispeak;
////                                        }
//                                    }
//
//                                    if(peakhub != -1 && peakhub != hID){
//                                        assert(PruningPointNew[curID].find(hID) != PruningPointNew[curID].end() && PruningPointNew[curID][hID].find(s) != PruningPointNew[curID][hID].end());
//                                        assert(PruningPointNew[s].find(hID) != PruningPointNew[s].end() && PruningPointNew[s][hID].find(curID) != PruningPointNew[s][hID].end());
//                                        outdatedPruning.insert(make_tuple(curID,hID,s));
//                                        outdatedPruning.insert(make_tuple(s,hID,curID));
//
//                                        PruningPointNew[curID][peakhub].insert(s);
//                                        PruningPointNew[s][peakhub].insert(curID);
//
//                                    }
//                                }
//
//                            }
////                            else if(NodeOrder[s]>NodeOrder[curID] && NoSupportedPair.find(make_pair(curID,s))==NoSupportedPair.end()){//
//                            else if(NodeOrder[s]>NodeOrder[curID]){// && NoSupportedPair.find(make_pair(curID,s))==NoSupportedPair.end()
//                                disvally=DisQueryVally(curID,s,Neighbors,Label);
//                                peakPair=DisQueryPeak2(curID,s,Label);//the old P2 distance dominates P1 distance, but the increase of d(curID,hID) may violate the previous balance
////                                disvally=DisQueryVallyVert(curID,s,Neighbors,Label);
////                                peakPair=DisQueryPeakVert(curID,s,Neighbors,Label);
//                                int dispeakOld = 0;
//                                if(Label[s].find(hID) != Label[s].end()){
//                                    dispeakOld = hDis + Label[s][hID];
//                                }
//                                assert(dispeakOld <= disvally);
//                                dispeak=peakPair.first; peakhub=peakPair.second;
//                                if(dispeak>disvally){
//                                    WaitProPTem.insert(curID);
//                                    ChangePTem[curID].push_back(make_pair(s,hID));
//                                    Label[curID][s]=disvally;//should be the final correct value
////                                    NoSupportedPair.insert(make_pair(curID,s));//insert to indicate that there is no pruning between curID and s
//                                    outdatedPruning.insert(make_tuple(curID,hID,s));
//                                    outdatedPruning.insert(make_tuple(s,hID,curID));
//
//                                    //cout<<"1--2 "<<curID<<" "<<s<<" "<<disvally<<endl;
//                                }
//                                else {//if dispeak<=disvally
//                                    if(Label[curID].find(s) != Label[curID].end()) {
//                                        Label[curID].erase(s);
////                                        if (dispeak != Label[curID][s]) {
////                                            Label[curID][s] = dispeak;//should be correct value
////                                        }
//                                    }
//
//                                    if(peakhub != -1 && peakhub != hID) {
//                                        assert(PruningPointNew[curID].find(hID) != PruningPointNew[curID].end() && PruningPointNew[curID][hID].find(s) != PruningPointNew[curID][hID].end());
//                                        assert(PruningPointNew[s].find(hID) != PruningPointNew[s].end() && PruningPointNew[s][hID].find(curID) != PruningPointNew[s][hID].end());
//                                        outdatedPruning.insert(make_tuple(curID, hID, s));
//                                        outdatedPruning.insert(make_tuple(s, hID, curID));
//
//                                        PruningPointNew[curID][peakhub].insert(s);
//                                        PruningPointNew[s][peakhub].insert(curID);
//                                    }
//                                }
//
//                            }
//
//                        }
//                        for(auto it2=outdatedPruning.begin();it2!=outdatedPruning.end();++it2){
//                            if(PruningPointNew[get<0>(*it2)].find(get<1>(*it2)) != PruningPointNew[get<0>(*it2)].end()) {//if found
//                                PruningPointNew[get<0>(*it2)][get<1>(*it2)].erase(get<2>(*it2));
//                            }
////                        PruningPoint[get<1>(*it)][get<0>(*it)].erase(get<2>(*it));
////                        if(PruningPointNew[get<0>(*it)][get<1>(*it)].empty())
////                            PruningPointNew[get<0>(*it)].erase(get<1>(*it));
//                        }
//                    }
//
//                }
//            }
//
//            //ChangeP->CHangeP: AL2->AL2
////            int v,u,neiid,neiw;
//            for(auto itp=WaitProP.begin();itp!=WaitProP.end();itp++){
//                int curID=*itp;//v is curID
//                int neiid,neiw;
//                for(int k=0;k<ChangeP[curID].size();k++){
//                    int hID=ChangeP[curID][k].first;//u is hID
//                    int hubID=ChangeP[curID][k].second;//
//
//                    if(curID == 146346){
////                        cout<<"AL2: "<<curID<<" "<<hID<<endl;
//                    }
//                    outdatedPruning.clear();
//                    for(int l=0;l<Neighbors[curID].size();l++){
//                        neiid=Neighbors[curID][l].first; neiw=Neighbors[curID][l].second;
//                        if(neiw+hDis==Label[neiid][hID]) {//check the affected label, if the shortest path between neiID and hID pass e(curID,hID)
//                        }
//                        if(NodeOrder[neiid]<NodeOrder[hID]){
//                            disvally=DisQueryVally(neiid, hID,Neighbors,Label);
//                            peakPair=DisQueryPeak2(neiid, hID,Label);
////                            disvally=DisQueryVallyVert(neiid, hID,Neighbors,Label);
////                            peakPair=DisQueryPeakVert(neiid, u,Neighbors,Label);
//                            dispeak=peakPair.first; peakhub=peakPair.second;
//                            if(hID == 207156 && neiid == 146346){
//                                cout<<"neiid: "<<neiid<<" "<<hID<<" ; "<<disvally<<" "<<dispeak<<endl;///
//                            }
//                            if(dispeak>disvally){
//                                if(Label[neiid].find(hID)==Label[neiid].end() || (Label[neiid].find(hID)!=Label[neiid].end() && Label[neiid][hID]!=disvally)){//if not found or found but disvally is lower///
//                                    WaitProPTem.insert(neiid);
////                                    ChangePTem[neiid].push_back(hID);
//                                    ChangePTem[neiid].push_back(make_pair(hID,hubID));
//                                    Label[neiid][hID]=disvally;//should be the final correct value
//                                    NoSupportedPair.insert(make_pair(neiid,hID));
////                                    outdatedPruning.insert(make_tuple(neiid,hubID,hID));//
////                                    outdatedPruning.insert(make_tuple(hID,hubID,neiid));//
//
//                                    //cout<<"2--2 "<<neiid<<" "<<u<<" "<<disvally<<endl;
//                                }
//                            }
//                            else {//if dispeak<=disvally
//                                if(Label[neiid].find(hID) != Label[neiid].end()) {
//                                    Label[neiid].erase(hID);
////                                    if (dispeak != Label[neiid][u]) {
////                                        Label[neiid][u] = dispeak;
////                                    }
//                                }
//
//                                if(peakhub != -1 && peakhub != hubID){
////                                    outdatedPruning.insert(make_tuple(neiid,hubID,hID));//
////                                    outdatedPruning.insert(make_tuple(hID,hubID,neiid));//
//
//                                    PruningPointNew[neiid][peakhub].insert(hID);
//                                    PruningPointNew[hID][peakhub].insert(neiid);
//                                }
//                            }
//                        }
//                    }
//                    for(auto it2=outdatedPruning.begin();it2!=outdatedPruning.end();++it2){
//                        if(PruningPointNew[get<0>(*it2)].find(get<1>(*it2)) != PruningPointNew[get<0>(*it2)].end()){//if found
//                            PruningPointNew[get<0>(*it2)][get<1>(*it2)].erase(get<2>(*it2));
////                        PruningPoint[get<1>(*it)][get<0>(*it)].erase(get<2>(*it));
//                        }
//                    }
//                }
//            }
//
//
//            WaitPro=WaitProTem;
//            Change=ChangeTem;
//            WaitProP=WaitProPTem;
//            ChangeP=ChangePTem;
//        }
//
//    }
//}

//correct new version by set
//Correct for pruning point
/*void Graph::IncreasePSL(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label,vector<unordered_map<int,set<int>>> &PruningPointNew){
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

    int LID,HID;
    if(NodeOrder[a]>NodeOrder[b]){
        LID=b; HID=a;
    }else{
        LID=a; HID=b;
    }

    int tempNDis;
    int dis,disvally,dispeak,peakhub;
    pair<int,int> peakPair;//distance, hubID
    set<tuple<int,int,int>> outdatedPruning;//<nodeID,supportNode,prunedID>
    //activate or not
    dis=DisQueryLower1(LID,HID,Neighbors,Label);//d1', via a neighbor of LID (which has lower order than HID) to reach HID
    dispeak=DisQueryPeak(LID,HID,Label);//d2, via the common hub vertex of LID and HID
    if(dispeak<=oldW)//if d2 is lower than oldW, the increase of oldW will not affect the shortest distance between LID and HID
        return;
    if(ifDebug){
        if(dis == INF)
            cout<<DisQueryLower1(LID,HID,Neighbors,Label);
        if(dis==oldW){
            cout<<dis<<" "<<oldW<<endl;
        }
    }

//    if(Label[196778].find(195124) != Label[196778].end()){
//        cout<<Label[196778][195124]<<" "<<DijkstraCore(196778,195124)<<endl;
//    }

    if(Label[LID].find(HID)!=Label[LID].end() && dis>oldW){//if d1' is larger than oldW, it indicates the shortest distance between LID and HID is equal to oldW, index update triggered
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
//        disvally=DisQueryVallyVert(LID,HID,Neighbors,Label);//shortest distance query through the neighbor's label, correct one
        Label[LID][HID]=disvally;//correct

        //cout<<"start "<<LID<<" "<<HID<<" "<<disvally<<endl;

        //affected by the w(a,b)
        int hubID, hDis;
        int dis, cnt;
        for(auto it=Label[HID].begin();it!=Label[HID].end();++it){//check the vertex has higher order than HID and update LID's label (except HID)
            hubID=(*it).first; hDis=(*it).second;
            if(hubID!=HID && Label[LID].find(hubID)!=Label[LID].end()){// && oldW+hDis==Label[LID][hubID]
                if(oldW+hDis==Label[LID][hubID]){//if the shortest path between LID and hubID pass through e(LID,HID)
                    disvally=DisQueryVally(LID,hubID,Neighbors,Label);
//                    disvally=DisQueryVallyVert(LID,hubID,Neighbors,Label);
                    if(Label[LID][hubID]<disvally){//if the new distance of P1 is higher, update it
                        WaitPro.insert(LID);
                        Change[LID].push_back(make_pair(hubID, oldW+hDis));//record the changed label of LID with the old distance
//                        cout<<LID<<" "<<hubID<<" "<<oldW+hDis<<endl;

                        peakPair=DisQueryPeak2(LID,hubID, Label);
                        dispeak=peakPair.first; peakhub=peakPair.second;
                        if(dispeak>disvally){
                            Label[LID][hubID]=disvally;//may not be correct at this moment
                            tempNDis = disvally;
                        }else{
//                            Label[LID][hubID]=dispeak;//may not be correct at this moment
                            Label[LID].erase(hubID);
                            tempNDis = dispeak;

                            PruningPointNew[LID][peakhub].insert(hubID);
                            PruningPointNew[hubID][peakhub].insert(LID);
                        }
//                        Label[LID][hubID]=disvally;//may not be correct at this moment
                        //cout<<"weight "<<LID<<" "<<hubID<<" "<<disvally<<endl;

                    }
                }

            }
        }

        for(auto it=Label[LID].begin();it!=Label[LID].end();it++){//update HID's label
            hubID=(*it).first; hDis=(*it).second;
            if(Label[HID].find(hubID)!=Label[HID].end()){//&& oldW+hDis==Label[HID][hubID]
                if(oldW+hDis==Label[HID][hubID]){//if the shortest path between HID and hubID pass through e(LID,HID)
                    disvally=DisQueryVally(HID,hubID,Neighbors,Label);
//                    disvally=DisQueryVallyVert(HID,hubID,Neighbors,Label);
                    if(Label[HID][hubID]<disvally){
                        WaitPro.insert(HID);
                        Change[HID].push_back(make_pair(hubID, oldW+hDis));//record the changed label
//                        cout<<HID<<" "<<hubID<<" "<<oldW+hDis<<endl;
//                        peakPair=DisQueryPeak2Vert(HID,hubID,Neighbors,Label);
                        peakPair=DisQueryPeak2(HID,hubID,Label);
                        dispeak=peakPair.first; peakhub=peakPair.second;
                        if(dispeak>disvally){
                            Label[HID][hubID]=disvally;//may not be correct at this moment, e.g., larger than correct value
                            tempNDis = disvally;
                        }else{
//                            Label[HID][hubID]=dispeak;
                            Label[HID].erase(hubID);
                            tempNDis=dispeak;

                            PruningPointNew[HID][peakhub].insert(hubID);
                            PruningPointNew[hubID][peakhub].insert(HID);
                        }
//                        Label[HID][hubID]=disvally;//may not be correct at this moment, e.g., larger than correct value
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
                vector<pair<int,int>> curChange=Change[curID];//the changed label
                int neiID, neiDis, hID, hDis;

                for(int j=0;j<curChange.size();j++){
                    hID=curChange[j].first; hDis=curChange[j].second;//the hub of changed label
//                    if(Label[196778].find(195124) != Label[196778].end()){
//                        cout<<Label[196778][195124]<<" "<<DijkstraCore(196778,195124)<<" "<<curID<<" "<<hID<<endl;
//                    }
                    //Change->Change: AL1->AL1
                    //check neighbors of curID
                    for(int k=0;k<Neighbors[curID].size();k++){//identify the affected neighbors of curID
                        neiID=Neighbors[curID][k].first; neiDis=Neighbors[curID][k].second;
                        if(Label[neiID].find(hID)!=Label[neiID].end() ){//&& neiDis+hDis==Label[neiID][hID]
                            if(neiDis+hDis==Label[neiID][hID]){//if the shortest path between neiID and hID pass e(curID,hID)
//                                if(curID == 196731 && hID == 195124 && neiID == 196755){
//                                    cout<<curID<<"("<<NodeOrder[curID]<<") "<<hID<<"("<<NodeOrder[hID]<<") "<<neiID<<"("<<NodeOrder[neiID]<<")"<<endl;
//                                    cout<<"196778 "<<NodeOrder[196778]<<endl;
//                                    if(Label[curID].find(hID) != Label[curID].end()){
//                                        cout<<Label[curID][hID]<<endl;
//                                    }
//                                }
                                disvally=DisQueryVally(neiID,hID,Neighbors,Label);
//                                disvally=DisQueryVallyVert2(neiID,hID,Neighbors,Label, hDis);
//                                assert(Label[curID].find(hID) != Label[curID].end());
//                                disvally=neiDis+Label[curID][hID];
                                if(Label[neiID][hID]<disvally){
                                    WaitProTem.insert(neiID);
                                    ChangeTem[neiID].push_back(make_pair(hID, neiDis+hDis));

                                    peakPair=DisQueryPeak2Vert(neiID,hID,Neighbors,Label);
                                    dispeak=peakPair.first; peakhub=peakPair.second;
                                    if(dispeak>disvally){
                                        Label[neiID][hID]=disvally;//may not be correct at this moment, e.g., larger than correct value
                                        tempNDis=disvally;
                                    }else{
//                                        Label[neiID][hID]=dispeak;
                                        Label[neiID].erase(hID);
                                        tempNDis=dispeak;


                                        PruningPointNew[neiID][peakhub].insert(hID);
                                        PruningPointNew[hID][peakhub].insert(neiID);

                                    }
//                                    Label[neiID][hID]=disvally;//may not be the final correct value at this moment
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

//                            if(NodeOrder[s]<NodeOrder[curID] && NoSupportedPair.find(make_pair(s,curID))==NoSupportedPair.end()){//it is not in NoSupportedPair
                            if(NodeOrder[s]<NodeOrder[curID]){//it is not in NoSupportedPair
                                disvally=DisQueryVally(s,curID,Neighbors,Label);
                                peakPair=DisQueryPeak2(s,curID,Label);
//                                disvally=DisQueryVallyVert(s,curID,Neighbors,Label);
//                                peakPair=DisQueryPeakVert(s,curID,Neighbors,Label);
                                dispeak=peakPair.first; peakhub=peakPair.second;
                                if(dispeak>disvally){
                                    WaitProPTem.insert(s);
                                    ChangePTem[s].push_back(curID);
                                    Label[s][curID]=disvally;
//                                    NoSupportedPair.insert(make_pair(s,curID));
                                    outdatedPruning.insert(make_tuple(curID,hID,s));
                                    outdatedPruning.insert(make_tuple(s,hID,curID));


                                    //cout<<"1--2 "<<s<<" "<<curID<<" "<<disvally<<endl;
                                }
                                else {//if dispeak<=disvally
                                    if(Label[s].find(curID) != Label[s].end()) {
                                        Label[s].erase(curID);
//                                        if (dispeak != Label[s][curID]) {
//                                            Label[s][curID] = dispeak;
//                                        }

                                    }

                                    if(peakhub != -1 && peakhub != hID){
                                        outdatedPruning.insert(make_tuple(curID,hID,s));
                                        outdatedPruning.insert(make_tuple(s,hID,curID));

                                        PruningPointNew[curID][peakhub].insert(s);
                                        PruningPointNew[s][peakhub].insert(curID);


                                    }
                                }

                            }else if(NodeOrder[s]>NodeOrder[curID]){// && NoSupportedPair.find(make_pair(curID,s))==NoSupportedPair.end()
                                disvally=DisQueryVally(curID,s,Neighbors,Label);
                                peakPair=DisQueryPeak2(curID,s,Label);//the old P2 distance dominates P1 distance, but the increase of d(curID,hID) may violate the previous balance
//                                disvally=DisQueryVallyVert(curID,s,Neighbors,Label);
//                                peakPair=DisQueryPeakVert(curID,s,Neighbors,Label);
                                dispeak=peakPair.first; peakhub=peakPair.second;
                                if(dispeak>disvally){
                                    WaitProPTem.insert(curID);
                                    ChangePTem[curID].push_back(s);
                                    Label[curID][s]=disvally;//should be the final correct value
//                                    NoSupportedPair.insert(make_pair(curID,s));//insert to indicate that there is no pruning between curID and s
                                    outdatedPruning.insert(make_tuple(curID,hID,s));
                                    outdatedPruning.insert(make_tuple(s,hID,curID));


                                    //cout<<"1--2 "<<curID<<" "<<s<<" "<<disvally<<endl;
                                }
                                else {//if dispeak<=disvally
                                    if(Label[curID].find(s) != Label[curID].end()) {
                                        Label[curID].erase(s);
//                                        if (dispeak != Label[curID][s]) {
//                                            Label[curID][s] = dispeak;//should be correct value
//                                        }

                                    }

                                    if(peakhub != -1 && peakhub != hID) {
                                        outdatedPruning.insert(make_tuple(curID, hID, s));
                                        outdatedPruning.insert(make_tuple(s, hID, curID));

//                                        for(auto ii=peakhubs.begin();ii!=peakhubs.end();++ii) {
//                                            peakhub = *ii;

                                        PruningPointNew[curID][peakhub].insert(s);
                                        PruningPointNew[s][peakhub].insert(curID);
//                                        }
                                    }
                                }

                            }
                        }
                    }


                    for(auto it=outdatedPruning.begin();it!=outdatedPruning.end();++it){
                        if(PruningPointNew[get<0>(*it)].find(get<1>(*it)) != PruningPointNew[get<0>(*it)].end()) {//if found
                            PruningPointNew[get<0>(*it)][get<1>(*it)].erase(get<2>(*it));
                        }
//                        PruningPoint[get<1>(*it)][get<0>(*it)].erase(get<2>(*it));
//                        if(PruningPointNew[get<0>(*it)][get<1>(*it)].empty())
//                            PruningPointNew[get<0>(*it)].erase(get<1>(*it));
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
                            peakPair=DisQueryPeak2(neiid, u,Label);
//                            disvally=DisQueryVallyVert(neiid, u,Neighbors,Label);
//                            peakPair=DisQueryPeakVert(neiid, u,Neighbors,Label);
                            dispeak=peakPair.first; peakhub=peakPair.second;
                            if(dispeak>disvally){
                                if(Label[neiid].find(u)==Label[neiid].end() || (Label[neiid].find(u)!=Label[neiid].end() && Label[neiid][u]!=disvally)){//if not found or found but disvally is lower///
                                    WaitProPTem.insert(neiid);
                                    ChangePTem[neiid].push_back(u);
                                    Label[neiid][u]=disvally;//should be the final correct value
//                                    NoSupportedPair.insert(make_pair(neiid,u));
//                                    outdatedPruning.insert(make_tuple(neiid,v,u));//
//                                    outdatedPruning.insert(make_tuple(u,v,neiid));//

                                    //cout<<"2--2 "<<neiid<<" "<<u<<" "<<disvally<<endl;
                                }
                            }
                            else {//if dispeak<=disvally
                                if(Label[neiid].find(u) != Label[neiid].end()) {
                                    Label[neiid].erase(u);
//                                    if (dispeak != Label[neiid][u]) {
//                                        Label[neiid][u] = dispeak;
//                                    }

                                }

                                if(peakhub != -1 && peakhub != v){
//                                    outdatedPruning.insert(make_tuple(neiid,v,u));//
//                                    outdatedPruning.insert(make_tuple(u,v,neiid));//

//                                    for(auto ii=peakhubs.begin();ii!=peakhubs.end();++ii) {
//                                        peakhub = *ii;

                                    PruningPointNew[neiid][peakhub].insert(u);
                                    PruningPointNew[u][peakhub].insert(neiid);
//                                    }

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
}*/

//// new version by set
//void Graph::IncreasePSL(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label,vector<unordered_map<int,set<int>>> &PruningPointNew){
//    for(int i=0;i<Neighbors[a].size();i++){
//        if(Neighbors[a][i].first==b){
////            cout<<i<<" "<<a<<" "<<b<<" "<<Neighbors[a][i].second<<" "<<newW<<endl;
//            if(oldW != Neighbors[a][i].second){
//                cout<<"Old edge weight is incorrect! "<<a<<" "<<b<<": "<<oldW<<" "<<Neighbors[a][i].second<<endl;
//                oldW = Neighbors[a][i].second;
//            }
//            Neighbors[a][i].second=newW;
//            break;
//        }
//    }
//    for(int i=0;i<Neighbors[b].size();i++){
//        if(Neighbors[b][i].first==a){
////            cout<<i<<" "<<a<<" "<<b<<" "<<Neighbors[b][i].second<<" "<<newW<<endl;
//            if(oldW != Neighbors[b][i].second){
//                cout<<"Old edge weight is incorrect! "<<b<<" "<<a<<": "<<oldW<<" "<<Neighbors[b][i].second<<endl;
//                oldW = Neighbors[b][i].second;
//            }
//            Neighbors[b][i].second=newW;
//            break;
//        }
//    }
//
//    bool ifDebug = false;//false;
////    ifDebug = true;
//
//    int LID,HID;
//    if(NodeOrder[a]>NodeOrder[b]){
//        LID=b; HID=a;
//    }else{
//        LID=a; HID=b;
//    }
//
//    int tempNDis;
//    int dis,disvally,dispeak,peakhub;
//    pair<int,int> peakPair;//distance, hubID
//    set<tuple<int,int,int>> outdatedPruning;//<nodeID,supportNode,prunedID>
//    //activate or not
//    dis=DisQueryLower1(LID,HID,Neighbors,Label);//d1', via a neighbor of LID (which has lower order than HID) to reach HID
//    dispeak=DisQueryPeak2(LID,HID,Label).first;//d2, via the common hub vertex of LID and HID
//    if(dispeak<=oldW)//if d2 is lower than oldW, the increase of oldW will not affect the shortest distance between LID and HID
//        return;
//    if(ifDebug){
//        if(dis == INF)
//            cout<<DisQueryLower1(LID,HID,Neighbors,Label);
//        if(dis==oldW){
//            cout<<dis<<" "<<oldW<<endl;
//        }
//    }
//
////    if(Label[196778].find(195124) != Label[196778].end()){
////        cout<<Label[196778][195124]<<" "<<DijkstraCore(196778,195124)<<endl;
////    }
//
//    if(dis>oldW){//if d1' is larger than oldW, it indicates the shortest distance between LID and HID is equal to oldW, index update triggered
//        vector<vector<pair<int,int>>> Change;//the label that is changed
//        vector<pair<int,int>> vec;
//        Change.assign(nodenum,vec);
//        set<int> WaitPro;//the vertex waited for label propagation, i.e., AL1
//        vector<vector<int>> ChangeP;
//        vector<int> vecint;
//        ChangeP.assign(nodenum,vecint);
//        set<int> WaitProP;//the vertex waited for label propagation, i.e., AL2
//
//
//        WaitPro.insert(LID);
//        Change[LID].push_back(make_pair(HID, oldW));
//        disvally=DisQueryVally(LID,HID,Neighbors,Label);//shortest distance query through the neighbor's label
////        disvally=DisQueryVallyVert(LID,HID,Neighbors,Label);//shortest distance query through the neighbor's label, correct one
//        Label[LID][HID]=disvally;//correct
//
//        //cout<<"start "<<LID<<" "<<HID<<" "<<disvally<<endl;
//
//        //affected by the w(a,b)
//        int hubID, hDis;
//        int dis, cnt;
//        for(auto it=Label[HID].begin();it!=Label[HID].end();++it){//check the vertex has higher order than HID and update LID's label (except HID)
//            hubID=(*it).first; hDis=(*it).second;
//            if(hubID!=HID && Label[LID].find(hubID)!=Label[LID].end()){// && oldW+hDis==Label[LID][hubID]
//                if(oldW+hDis==Label[LID][hubID]){//if the shortest path between LID and hubID pass through e(LID,HID)
////                    disvally=DisQueryVally(LID,hubID,Neighbors,Label);
//                    disvally=DisQueryVallyVert(LID,hubID,Neighbors,Label);
//                    if(Label[LID][hubID]<disvally){//if the new distance of P1 is higher, update it
//                        WaitPro.insert(LID);
//                        Change[LID].push_back(make_pair(hubID, oldW+hDis));//record the changed label of LID with the old distance
////                        cout<<LID<<" "<<hubID<<" "<<oldW+hDis<<endl;
//
//                        peakPair=DisQueryPeak2Vert(LID,hubID,Neighbors, Label); dispeak=peakPair.first; peakhub=peakPair.second;
//                        if(dispeak>disvally){
//                            Label[LID][hubID]=disvally;//may not be correct at this moment
//                            tempNDis = disvally;
//                        }else{
////                            Label[LID][hubID]=dispeak;//may not be correct at this moment
//                            Label[LID].erase(hubID);
//                            tempNDis = dispeak;
//
////                            for(auto ii=peakhubs.begin();ii!=peakhubs.end();++ii){
////                                peakhub = *ii;
////                                if(PruningPointNew[LID].find(peakhub) == PruningPointNew[LID].end()){//if not found
////                                    PruningPointNew[LID].insert({peakhub,set<int>()});
////                                    PruningPointNew[hubID].insert({peakhub,set<int>()});
////                                }
//                                PruningPointNew[LID][peakhub].insert(hubID);
//                                PruningPointNew[hubID][peakhub].insert(LID);
////                            }
//
//                        }
////                        Label[LID][hubID]=disvally;//may not be correct at this moment
//                        //cout<<"weight "<<LID<<" "<<hubID<<" "<<disvally<<endl;
//                        if(ifDebug){
//                            int tempD = DijkstraCore(LID,hubID);
//                            if(tempD != tempNDis){
//                                cout<<"LID label error. "<<LID<<" "<<hubID<<": "<<tempNDis<<" "<<tempD<<endl;
//                            }
//                        }
//
//                    }
//                }
//
//            }
//        }
//
//        for(auto it=Label[LID].begin();it!=Label[LID].end();it++){//update HID's label
//            hubID=(*it).first; hDis=(*it).second;
//            if(Label[HID].find(hubID)!=Label[HID].end()){//&& oldW+hDis==Label[HID][hubID]
//                if(oldW+hDis==Label[HID][hubID]){//if the shortest path between HID and hubID pass through e(LID,HID)
////                    disvally=DisQueryVally(HID,hubID,Neighbors,Label);
//                    disvally=DisQueryVallyVert(HID,hubID,Neighbors,Label);
//                    if(Label[HID][hubID]<disvally){
//                        WaitPro.insert(HID);
//                        Change[HID].push_back(make_pair(hubID, oldW+hDis));//record the changed label
////                        cout<<HID<<" "<<hubID<<" "<<oldW+hDis<<endl;
//                        peakPair=DisQueryPeak2Vert(HID,hubID,Neighbors,Label); dispeak=peakPair.first; peakhub=peakPair.second;
//                        if(dispeak>disvally){
//                            Label[HID][hubID]=disvally;//may not be correct at this moment, e.g., larger than correct value
//                            tempNDis = disvally;
//                        }else{
////                            Label[HID][hubID]=dispeak;
//                            Label[HID].erase(hubID);
//                            tempNDis=dispeak;
//
////                            for(auto ii=peakhubs.begin();ii!=peakhubs.end();++ii) {
////                                peakhub = *ii;
////                                if(PruningPointNew[HID].find(peakhub) == PruningPointNew[HID].end()){//if not found
////                                    PruningPointNew[HID].insert({peakhub,set<int>()});
////                                    PruningPointNew[hubID].insert({peakhub,set<int>()});
////                                }
//                                PruningPointNew[HID][peakhub].insert(hubID);
//                                PruningPointNew[hubID][peakhub].insert(HID);
////                            }
//                        }
////                        Label[HID][hubID]=disvally;//may not be correct at this moment, e.g., larger than correct value
//                        //cout<<"weight "<<HID<<" "<<hubID<<" "<<disvally<<endl;
//                        if(ifDebug){
//                            int tempD = DijkstraCore(HID,hubID);
//                            if(tempD != tempNDis){
//                                cout<<"HID label error. "<<HID<<" "<<hubID<<": "<<tempNDis<<" "<<tempD<<endl;
//                            }
//                        }
//
//                    }
//                }
//            }
//        }
//
//        while(WaitProP.size()>0 || WaitPro.size()>0){//the vertex set that corresponding labels have been changed
//            set<int> WaitProTem;
//            vector<vector<pair<int,int>>> ChangeTem;
//            vector<pair<int,int>> vec;
//            ChangeTem.assign(nodenum,vec);
//            set<int> WaitProPTem;
//            vector<int> vecint;
//            vector<vector<int>> ChangePTem;
//            ChangePTem.assign(nodenum,vecint);
//
//            //Change->Change & ChangeP: AL1->AL1 and AL1->AL2
//            for(auto it=WaitPro.begin();it!=WaitPro.end();it++){
//                int curID=*it;
//                vector<pair<int,int>> curChange=Change[curID];//the changed label
//                int neiID, neiDis, hID, hDis;
//
//                for(int j=0;j<curChange.size();j++){
//                    hID=curChange[j].first; hDis=curChange[j].second;//the hub of changed label
////                    if(Label[196778].find(195124) != Label[196778].end()){
////                        cout<<Label[196778][195124]<<" "<<DijkstraCore(196778,195124)<<" "<<curID<<" "<<hID<<endl;
////                    }
//                    //Change->Change: AL1->AL1
//                    //check neighbors of curID
//                    for(int k=0;k<Neighbors[curID].size();k++){//identify the affected neighbors of curID
//                        neiID=Neighbors[curID][k].first; neiDis=Neighbors[curID][k].second;
//                        if(Label[neiID].find(hID)!=Label[neiID].end() ){//&& neiDis+hDis==Label[neiID][hID]
//                            if(neiDis+hDis==Label[neiID][hID]){//if the shortest path between neiID and hID pass e(curID,hID)
////                                if(curID == 196731 && hID == 195124 && neiID == 196755){
////                                    cout<<curID<<"("<<NodeOrder[curID]<<") "<<hID<<"("<<NodeOrder[hID]<<") "<<neiID<<"("<<NodeOrder[neiID]<<")"<<endl;
////                                    cout<<"196778 "<<NodeOrder[196778]<<endl;
////                                    if(Label[curID].find(hID) != Label[curID].end()){
////                                        cout<<Label[curID][hID]<<endl;
////                                    }
////                                }
////                                disvally=DisQueryVally(neiID,hID,Neighbors,Label);
//                                disvally=DisQueryVallyVert2(neiID,hID,Neighbors,Label, hDis);
////                                assert(Label[curID].find(hID) != Label[curID].end());
////                                disvally=neiDis+Label[curID][hID];
//                                if(Label[neiID][hID]<disvally){
//                                    WaitProTem.insert(neiID);
//                                    ChangeTem[neiID].push_back(make_pair(hID, neiDis+hDis));
//                                    peakPair=DisQueryPeak2Vert(neiID,hID,Neighbors,Label); dispeak=peakPair.first; peakhub=peakPair.second;
//                                    if(dispeak>disvally){
//                                        Label[neiID][hID]=disvally;//may not be correct at this moment, e.g., larger than correct value
//                                        tempNDis=disvally;
//                                    }else{
////                                        Label[neiID][hID]=dispeak;
//                                        Label[neiID].erase(hID);
//                                        tempNDis=dispeak;
//
////                                        for(auto ii=peakhubs.begin();ii!=peakhubs.end();++ii) {
////                                            peakhub = *ii;
////                                            if(PruningPointNew[neiID].find(peakhub) == PruningPointNew[neiID].end()){//if not found
////                                                PruningPointNew[neiID].insert({peakhub,set<int>()});
////                                                PruningPointNew[hID].insert({peakhub,set<int>()});
////                                            }
//                                            PruningPointNew[neiID][peakhub].insert(hID);
//                                            PruningPointNew[hID][peakhub].insert(neiID);
////                                        }
//                                    }
////                                    Label[neiID][hID]=disvally;//may not be the final correct value at this moment
//                                    //cout<<"1--1 "<<neiID<<" "<<hID<<" "<<disvally<<endl;
//                                    if(ifDebug){
//                                        int tempD = DijkstraCore(neiID,hID);
//                                        if(tempD != tempNDis){
//                                            cout<<"AL1 to AL1. "<<neiID<<" "<<hID<<": "<<tempNDis<<" "<<tempD<<endl;
////                                            DijkstraCorePath(neiID,hID);
//                                        }
//                                    }
//                                }
//                            }
//
//                        }
//                    }
//                    //check neighbors of hID
//                    /*for(int k=0;k<Neighbors[hID].size();k++){//identify the affected neighbors of curID
//                        neiID=Neighbors[hID][k].first; neiDis=Neighbors[hID][k].second;
//                        if(Label[curID].find(neiID)!=Label[curID].end() ){//&& neiDis+hDis==Label[neiID][hID]
//                            if(neiDis+hDis==Label[curID][neiID]){//if the shortest path between neiID and hID pass e(curID,hID)
//                                if(curID == 196731 && hID == 195124 && neiID == 196755){
//                                    cout<<curID<<"("<<NodeOrder[curID]<<") "<<hID<<"("<<NodeOrder[hID]<<") "<<neiID<<"("<<NodeOrder[neiID]<<")"<<endl;
//                                    cout<<"196778 "<<NodeOrder[196778]<<endl;
//                                    if(Label[curID].find(hID) != Label[curID].end()){
//                                        cout<<Label[curID][hID]<<endl;
//                                    }
//                                }
////                                disvally=DisQueryVally(neiID,hID,Neighbors,Label);
//                                disvally=DisQueryVallyVert(curID,neiID,Neighbors,Label);
//                                assert(Label[curID].find(neiID) != Label[curID].end());
////                                disvally=neiDis+Label[curID][hID];
//                                if(Label[curID][neiID]<disvally){
//                                    WaitProTem.insert(curID);
//                                    ChangeTem[curID].push_back(make_pair(neiID, neiDis+hDis));
//                                    peakPair=DisQueryPeakVert(curID,neiID,Neighbors,Label); dispeak=peakPair.first; peakhub=peakPair.second;
//                                    if(dispeak>disvally){
//                                        Label[curID][neiID]=disvally;//may not be correct at this moment, e.g., larger than correct value
//                                        tempNDis = disvally;
//                                    }else{
////                                        Label[curID][neiID]=dispeak;
//                                        Label[curID].erase(neiID);
//                                        tempNDis = dispeak;
//                                        if(PruningPointNew[neiID].find(peakhub) == PruningPointNew[neiID].end()){//if not found
//                                            PruningPointNew[neiID].insert({peakhub,set<int>()});
//                                            PruningPointNew[curID].insert({peakhub,set<int>()});
//                                        }
//                                        PruningPointNew[neiID][peakhub].insert(curID);
//                                        PruningPointNew[curID][peakhub].insert(neiID);
//                                    }
////                                    Label[neiID][hID]=disvally;//may not be the final correct value at this moment
//                                    //cout<<"1--1 "<<neiID<<" "<<hID<<" "<<disvally<<endl;
//                                    if(ifDebug){
//                                        int tempD = DijkstraCore(curID, neiID);
//                                        if(tempD !=tempNDis){
//                                            cout<<"AL1 to AL1. "<<curID<<" "<<neiID<<": "<<tempNDis<<" "<<tempD<<endl;
////                                            DijkstraCorePath(neiID,hID);
//                                        }
//                                    }
//                                }
//                            }
//
//                        }
//                    }*/
//
//                    //Change->ChangeP: AL1->AL2
//                    outdatedPruning.clear();
//                    if(PruningPointNew[curID].find(hID)!=PruningPointNew[curID].end()){//check whether the pruned vertex should be inserted to curID again
////                        for(int snum=0;snum<PruningPointNew[curID][hID].size();snum++){
//                        for(auto snum=PruningPointNew[curID][hID].begin();snum!=PruningPointNew[curID][hID].end();++snum){
//                            int s=*snum;
//
//                            /*if(PruningPoint.find(make_pair(curID,hID))!=PruningPoint.end()){
//                                for(int snum=0;snum<PruningPoint[make_pair(curID,hID)].size();snum++){
//                                    int s=PruningPoint[make_pair(curID,hID)][snum];*/
////                            if(NodeOrder[s]<NodeOrder[curID] && NoSupportedPair.find(make_pair(s,curID))==NoSupportedPair.end()){//it is not in NoSupportedPair
//                            if(NodeOrder[s]<NodeOrder[curID]){//it is not in NoSupportedPair
//                                disvally=DisQueryVally(s,curID,Neighbors,Label);
//                                peakPair=DisQueryPeak2(s,curID,Label);
////                                disvally=DisQueryVallyVert(s,curID,Neighbors,Label);
////                                peakPair=DisQueryPeakVert(s,curID,Neighbors,Label);
//                                dispeak=peakPair.first; peakhub=peakPair.second;
//                                if(dispeak>disvally){
//                                    WaitProPTem.insert(s);
//                                    ChangePTem[s].push_back(curID);
//                                    Label[s][curID]=disvally;
////                                    NoSupportedPair.insert(make_pair(s,curID));
//                                    outdatedPruning.insert(make_tuple(curID,hID,s));
//                                    outdatedPruning.insert(make_tuple(s,hID,curID));
//
//                                    if(ifDebug){
//                                        int tempD = DijkstraCore(s,curID);
//                                        if(tempD != Label[s][curID]){
////                                            cout<<s<<" "<<curID<<": "<<Label[s][curID]<<" "<<tempD<<" "<<Dijkstra(s,curID,Neighbor)<<endl;//////
//                                            cout<<"AL1 to AL2. "<<s<<" "<<curID<<": "<<Label[s][curID]<<" "<<tempD<<endl;
//                                        }
//                                    }
//
//                                    //cout<<"1--2 "<<s<<" "<<curID<<" "<<disvally<<endl;
//                                }
//                                else {//if dispeak<=disvally
//                                    if(Label[s].find(curID) != Label[s].end()) {
//                                        if (dispeak != Label[s][curID]) {
//                                            Label[s][curID] = dispeak;
//                                        }
//                                        if(ifDebug){
//                                            int tempD = DijkstraCore(s,curID);
//                                            if(tempD != Label[s][curID]){
//                                                cout<<"New prune. AL1 to AL2. "<<s<<" "<<curID<<": "<<Label[s][curID]<<" "<<tempD<<endl;
//                                            }
//                                        }
//                                    }
//
//                                    if(peakhub != -1 && peakhub != hID){
//                                        outdatedPruning.insert(make_tuple(curID,hID,s));
//                                        outdatedPruning.insert(make_tuple(s,hID,curID));
//
////                                        for(auto ii=peakhubs.begin();ii!=peakhubs.end();++ii) {
////                                            peakhub = *ii;
////                                            if(PruningPointNew[curID].find(peakhub) == PruningPointNew[curID].end()){//if not found
////                                                PruningPointNew[curID].insert({peakhub,set<int>()});
////                                                PruningPointNew[s].insert({peakhub,set<int>()});
////                                            }
//                                            PruningPointNew[curID][peakhub].insert(s);
//                                            PruningPointNew[s][peakhub].insert(curID);
////                                        }
//
//                                    }
//                                }
//
//                            }else if(NodeOrder[s]>NodeOrder[curID]){// && NoSupportedPair.find(make_pair(curID,s))==NoSupportedPair.end()
//                                disvally=DisQueryVally(curID,s,Neighbors,Label);
//                                peakPair=DisQueryPeak2(curID,s,Label);//the old P2 distance dominates P1 distance, but the increase of d(curID,hID) may violate the previous balance
////                                disvally=DisQueryVallyVert(curID,s,Neighbors,Label);
////                                peakPair=DisQueryPeakVert(curID,s,Neighbors,Label);
//                                dispeak=peakPair.first; peakhub=peakPair.second;
//                                if(dispeak>disvally){
//                                    WaitProPTem.insert(curID);
//                                    ChangePTem[curID].push_back(s);
//                                    Label[curID][s]=disvally;//should be the final correct value
////                                    NoSupportedPair.insert(make_pair(curID,s));//insert to indicate that there is no pruning between curID and s
//                                    outdatedPruning.insert(make_tuple(curID,hID,s));
//                                    outdatedPruning.insert(make_tuple(s,hID,curID));
//
//                                    if(ifDebug){
//                                        int tempD = DijkstraCore(curID,s);
//                                        if(tempD != Label[curID][s]){
//                                            cout<<"AL1 to AL2. "<<curID<<" "<<s<<": "<<Label[curID][s]<<" "<<tempD<<endl;
//                                        }
//                                    }
//
//                                    //cout<<"1--2 "<<curID<<" "<<s<<" "<<disvally<<endl;
//                                }
//                                else {//if dispeak<=disvally
//                                    if(Label[curID].find(s) != Label[curID].end()) {
//                                        if (dispeak != Label[curID][s]) {
//                                            Label[curID][s] = dispeak;//should be correct value
//                                        }
//                                        if(ifDebug){
//                                            int tempD = DijkstraCore(curID,s);
//                                            if(tempD != Label[curID][s]){
//                                                cout<<"New prune. AL1 to AL2. "<<curID<<" "<<s<<": "<<Label[curID][s]<<" "<<tempD<<endl;
//                                            }
//                                        }
//                                    }
//
//                                    if(peakhub != -1 && peakhub != hID) {
//                                        outdatedPruning.insert(make_tuple(curID, hID, s));
//                                        outdatedPruning.insert(make_tuple(s, hID, curID));
//
////                                        for(auto ii=peakhubs.begin();ii!=peakhubs.end();++ii) {
////                                            peakhub = *ii;
////                                            if (PruningPointNew[curID].find(peakhub)==PruningPointNew[curID].end()) {//if not found
////                                                PruningPointNew[curID].insert({peakhub, set<int>()});
////                                                PruningPointNew[s].insert({peakhub, set<int>()});
////                                            }
//                                            PruningPointNew[curID][peakhub].insert(s);
//                                            PruningPointNew[s][peakhub].insert(curID);
////                                        }
//                                    }
//                                }
//
//                            }
//                        }
//                    }
//                    //check higher-order vertex of the affected distance label
//                    /*if(!PruningPointNew[hID].empty()){
//                        for(auto it=PruningPointNew[hID].begin();it!=PruningPointNew[hID].end();++it){
//                            int lowID = it->first;
//                            if(NodeOrder[lowID] < NodeOrder[curID]){
//                                cout<<"Higher order works!!"<<endl;
//                                for(auto it2=it->second.begin();it2!=it->second.end();++it2){
//                                    if(*it2 == curID){
//                                        disvally=DisQueryVally(lowID,curID,Neighbors,Label);
//                                        peakPair=DisQueryPeak(lowID,curID,Label);
//                                        dispeak=peakPair.first; peakhub=peakPair.second;
//                                        if(dispeak > disvally){
//                                            WaitProPTem.insert(lowID);
//                                            ChangePTem[lowID].push_back(curID);
//                                            Label[lowID][curID] = disvally;
//                                            outdatedPruning.insert(make_tuple(lowID,hID,curID));
//                                            outdatedPruning.insert(make_tuple(curID,hID,lowID));
//                                        }
//                                        else {
//                                            if(Label[lowID].find(curID) != Label[lowID].end()) {
//                                                if (dispeak != Label[lowID][curID]) {
//                                                    Label[lowID][curID] = dispeak;
//                                                }
//                                            }
//
//                                            if(peakhub != -1 && peakhub != hID){
//                                                outdatedPruning.insert(make_tuple(lowID,hID,curID));
//                                                outdatedPruning.insert(make_tuple(curID,hID,lowID));
//                                                if(PruningPointNew[lowID].find(peakhub) == PruningPointNew[lowID].end()){//if not found
//                                                    PruningPointNew[lowID].insert({peakhub,set<int>()});
//                                                    PruningPointNew[curID].insert({peakhub,set<int>()});
//                                                }
//                                                PruningPointNew[lowID][peakhub].insert(curID);
//                                                PruningPointNew[curID][peakhub].insert(lowID);
//
//                                            }
//                                        }
//
//                                    }
//                                }
//
//
//                            }
//
//                        }
//                    }*/
//
//
//                    for(auto it=outdatedPruning.begin();it!=outdatedPruning.end();++it){
//                        if(PruningPointNew[get<0>(*it)].find(get<1>(*it)) != PruningPointNew[get<0>(*it)].end()) {//if found
//                            PruningPointNew[get<0>(*it)][get<1>(*it)].erase(get<2>(*it));
//                        }
////                        PruningPoint[get<1>(*it)][get<0>(*it)].erase(get<2>(*it));
////                        if(PruningPointNew[get<0>(*it)][get<1>(*it)].empty())
////                            PruningPointNew[get<0>(*it)].erase(get<1>(*it));
//                    }
//
//                }
//            }
//
//            //ChangeP->CHangeP: AL2->AL2
//            int v,u,neiid,neiw;
//            outdatedPruning.clear();
//            for(auto itp=WaitProP.begin();itp!=WaitProP.end();itp++){
//                v=*itp;
//                for(int k=0;k<ChangeP[v].size();k++){
//                    u=ChangeP[v][k];
//                    for(int l=0;l<Neighbors[v].size();l++){
//                        neiid=Neighbors[v][l].first; neiw=Neighbors[v][l].second;
//                        if(NodeOrder[neiid]<NodeOrder[u]){
//                            disvally=DisQueryVally(neiid, u,Neighbors,Label);
//                            peakPair=DisQueryPeak2(neiid, u,Label);
////                            disvally=DisQueryVallyVert(neiid, u,Neighbors,Label);
////                            peakPair=DisQueryPeakVert(neiid, u,Neighbors,Label);
//                            dispeak=peakPair.first; peakhub=peakPair.second;
//                            if(dispeak>disvally){
//                                if(Label[neiid].find(u)==Label[neiid].end() || (Label[neiid].find(u)!=Label[neiid].end() && Label[neiid][u]>disvally)){//if not found or found but disvally is lower
//                                    WaitProPTem.insert(neiid);
//                                    ChangePTem[neiid].push_back(u);
//                                    Label[neiid][u]=disvally;//should be the final correct value
////                                    NoSupportedPair.insert(make_pair(neiid,u));
////                                    outdatedPruning.insert(make_tuple(neiid,v,u));//
////                                    outdatedPruning.insert(make_tuple(u,v,neiid));//
//                                    if(ifDebug){
//                                        int tempD = DijkstraCore(neiid,u);
//                                        if(tempD != Label[neiid][u]){
//                                            cout<<"AL2 to AL2. "<<neiid<<" "<<u<<": "<<Label[neiid][u]<<" "<<tempD<<endl;
//                                        }
//                                    }
//                                    //cout<<"2--2 "<<neiid<<" "<<u<<" "<<disvally<<endl;
//                                }
//                            }
//                            else {//if dispeak<=disvally
//                                if(Label[neiid].find(u) != Label[neiid].end()) {
//                                    if (dispeak != Label[neiid][u]) {
//                                        Label[neiid][u] = dispeak;
//                                    }
//                                    if(ifDebug){
//                                        int tempD = DijkstraCore(neiid,u);
//                                        if(tempD != Label[neiid][u]){
//                                            cout<<"New prune. AL2 to AL2. "<<neiid<<" "<<u<<": "<<Label[neiid][u]<<" "<<tempD<<endl;
//                                        }
//                                    }
//                                }
//
//                                if(peakhub != -1 && peakhub != v){
////                                    outdatedPruning.insert(make_tuple(neiid,v,u));//
////                                    outdatedPruning.insert(make_tuple(u,v,neiid));//
//
////                                    for(auto ii=peakhubs.begin();ii!=peakhubs.end();++ii) {
////                                        peakhub = *ii;
////                                        if(PruningPointNew[neiid].find(peakhub) == PruningPointNew[neiid].end()){//if not found
////                                            PruningPointNew[neiid].insert({peakhub,set<int>()});
////                                            PruningPointNew[u].insert({peakhub,set<int>()});
////                                        }
//                                        PruningPointNew[neiid][peakhub].insert(u);
//                                        PruningPointNew[u][peakhub].insert(neiid);
////                                    }
//
//                                }
//                            }
//
//
//                        }
//                    }
//                }
//            }
//            for(auto it=outdatedPruning.begin();it!=outdatedPruning.end();++it){
//                if(PruningPointNew[get<0>(*it)].find(get<1>(*it)) != PruningPointNew[get<0>(*it)].end()){//if found
//                    PruningPointNew[get<0>(*it)][get<1>(*it)].erase(get<2>(*it));
////                        PruningPoint[get<1>(*it)][get<0>(*it)].erase(get<2>(*it));
//                }
//            }
//
//            WaitPro=WaitProTem;
//            Change=ChangeTem;
//            WaitProP=WaitProPTem;
//            ChangeP=ChangePTem;
//        }
//
//    }
//}

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
int Graph::DisQueryVally2(int ID1, int ID2, int ID3, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label){
	int neiID,neiDis;
	int d=INF;
	for(int i=0;i<Neighbors[ID1].size();i++){
		neiID=Neighbors[ID1][i].first;
		neiDis=Neighbors[ID1][i].second;
        if(neiID != ID3){
            if(NodeOrder[neiID]<=NodeOrder[ID2] && Label[neiID].find(ID2)!=Label[neiID].end()){//if found
                if(neiDis+Label[neiID][ID2]<d){
                    d=neiDis+Label[neiID][ID2];
                }
            }
        }

	}
	return d;
}
//function of updating the label with the changed edge
int Graph::DisQueryVallyVert2(int ID1, int ID2, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label, int oldW){
    int neiID,neiDis;
    int d=INF;
    int tempUpdate=INF;
    for(int i=0;i<Neighbors[ID1].size();i++){
        neiID=Neighbors[ID1][i].first;
        neiDis=Neighbors[ID1][i].second;
        if(NodeOrder[neiID]<=NodeOrder[ID2] && Label[neiID].find(ID2)!=Label[neiID].end()){//if found
            if(neiDis+oldW == Label[neiID][ID2]){//if the shortest path between neiID and ID2 pass ID1
                //case 1: there is another shortest vally path that does not pass ID1 which is equal to d(neiID,ID2), then Label[neiID][ID2] is correct
                for(int p=0;p<Neighbors[neiID].size();++p){
                    int tempID=Neighbors[neiID][p].first, tempDis=Neighbors[neiID][p].second;
                    if(tempID != ID1){
                        if(NodeOrder[tempID]<=NodeOrder[ID2] && Label[tempID].find(ID2)!=Label[tempID].end()){
                            if(tempDis+Label[tempID][ID2] == Label[neiID][ID2]){//if there is another shortest vally path that does not pass ID1
                                if(neiDis+Label[neiID][ID2]<d){
                                    d=neiDis+Label[neiID][ID2];
                                }
                            }else if(tempDis+Label[tempID][ID2] > Label[neiID][ID2]){
                                if(tempDis+Label[tempID][ID2] < tempUpdate){
                                    tempUpdate = tempDis+Label[tempID][ID2];//tempUpdate is the correct value of Label[neiID][ID2] if Label[neiID][ID2] does not pass Label[ID1][ID2]
                                }
                                if(neiDis+tempUpdate<d){
                                    d=neiDis+tempUpdate;
                                }
                            }
                        }
                    }
                }
                //case 2: the Label[neiID][ID2] has changed (increased)
                //subcase 1: the changed Label[neiID][ID2] still pass Label[ID1][ID2], then Label[ID1][ID2] does not pass neiID, we left the update of Label[neiID][ID2] for future
                //subcase 2: the changed Label[neiID][ID2] does not pass Label[ID1][ID2], we should further check Label[neiID][ID2]
//                int disVally = DisQueryVally2(neiID,ID2,ID1,Neighbors,Label);
            }else{
                if(neiDis+Label[neiID][ID2]<d){
                    d=neiDis+Label[neiID][ID2];
                }
            }
        }
//        if(tempUpdate < INF){
//            Label[neiID][ID2] = tempUpdate;
//        }
    }
    return d;
}
//function of updating the label with the changed edge
int Graph::DisQueryVallyVert(int ID1, int ID2, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label){
    int neiID,neiDis;
    int d=INF;
    int tempUpdate=INF;
    for(int i=0;i<Neighbors[ID1].size();i++){
        neiID=Neighbors[ID1][i].first;
        neiDis=Neighbors[ID1][i].second;
        if(NodeOrder[neiID]<=NodeOrder[ID2] && Label[neiID].find(ID2)!=Label[neiID].end()){//if found
            if(neiDis+Label[ID1][ID2] == Label[neiID][ID2]){//if the shortest path between neiID and ID2 pass ID1
                //case 1: there is another shortest vally path that does not pass ID1 which is equal to d(neiID,ID2), then Label[neiID][ID2] is correct
                for(int p=0;p<Neighbors[neiID].size();++p){
                    int tempID=Neighbors[neiID][p].first, tempDis=Neighbors[neiID][p].second;
                    if(tempID != ID1){
                        if(NodeOrder[tempID]<=NodeOrder[ID2] && Label[tempID].find(ID2)!=Label[tempID].end()){
                            if(tempDis+Label[tempID][ID2] == Label[neiID][ID2]){//if there is another shortest vally path that does not pass ID1
                                if(neiDis+Label[neiID][ID2]<d){
                                    d=neiDis+Label[neiID][ID2];
                                }
                            }else if(tempDis+Label[tempID][ID2] > Label[neiID][ID2]){
                                if(tempDis+Label[tempID][ID2] < tempUpdate){
                                    tempUpdate = tempDis+Label[tempID][ID2];//tempUpdate is the correct value of Label[neiID][ID2] if Label[neiID][ID2] does not pass Label[ID1][ID2]
                                }
                                if(neiDis+tempUpdate<d){
                                    d=neiDis+tempUpdate;
                                }
                            }
                        }
                    }
                }
                //case 2: the Label[neiID][ID2] has changed (increased)
                //subcase 1: the changed Label[neiID][ID2] still pass Label[ID1][ID2], then Label[ID1][ID2] does not pass neiID, we left the update of Label[neiID][ID2] for future
                //subcase 2: the changed Label[neiID][ID2] does not pass Label[ID1][ID2], we should further check Label[neiID][ID2]
//                int disVally = DisQueryVally2(neiID,ID2,ID1,Neighbors,Label);
            }else{
                if(neiDis+Label[neiID][ID2]<d){
                    d=neiDis+Label[neiID][ID2];
                }
            }
        }
//        if(tempUpdate < INF){
//            Label[neiID][ID2] = tempUpdate;
//        }
    }
    return d;
}
//new version. function of computing the shortest distance through the label of ID1's neighbors, i.e., d1.
int Graph::DisQueryVallyNew(int ID1, int ID2, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label){
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
//debug version. function of computing the shortest distance through the label of ID1's neighbors, i.e., d1.
int Graph::DisQueryVallyDebug(int ID1, int ID2, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label){
    int neiID,neiDis;
    int d=INF;
    for(int i=0;i<Neighbors[ID1].size();i++){
        neiID=Neighbors[ID1][i].first;
        neiDis=Neighbors[ID1][i].second;
        if(NodeOrder[neiID]<=NodeOrder[ID2] && Label[neiID].find(ID2)!=Label[neiID].end()){//if found
            if(neiDis+Label[neiID][ID2]<d){
                d=neiDis+Label[neiID][ID2];
                cout<<ID1<<" "<<ID2<<": "<<neiID<<" "<<d<<" "<<neiDis<<" "<<Label[neiID][ID2]<<endl;
            }
        }
    }
    return d;
}
//function of computing d2 through neighbor's label, new version
pair<int,set<int>> Graph::DisQueryPeak2Vert2(int ID1, int ID2, vector<vector<pair<int,int>>> &Neighbors, vector<unordered_map<int,int>> &Label){
    int d=INF;
    unordered_map<int,int>::iterator it;
    int hub, dis1, finalHub=-1;
    int neiID,neiDis;
    set<int> hubs;
    for(int i=0;i<Neighbors[ID1].size();i++){
        neiID=Neighbors[ID1][i].first;
        neiDis=Neighbors[ID1][i].second;

        for(it=Label[neiID].begin();it!=Label[neiID].end();++it){
            hub=(*it).first;
            dis1=(*it).second;
            if(hub!=ID2 && Label[ID2].find(hub)!=Label[ID2].end()){//hub is the hub of ID1 and ID2 and hub != ID2
                if(neiDis+dis1+Label[ID2][hub]<d){
                    d=neiDis+dis1+Label[ID2][hub];
                    finalHub = hub;
                    hubs.clear();
                    hubs.insert(finalHub);
                }
//                else if(neiDis+dis1+Label[ID2][hub] == d){
//                    hubs.insert(hub);
//                    cout<<"!!!!!!!!!!!!!!!!!!"<<endl;
//                }
            }
        }


    }
    return make_pair(d,hubs);
}
//function of computing d2 through neighbor's label, new version
pair<int,int> Graph::DisQueryPeak2Vert(int ID1, int ID2, vector<vector<pair<int,int>>> &Neighbors, vector<unordered_map<int,int>> &Label){
    int d=INF;
    unordered_map<int,int>::iterator it;
    int hub, dis1, finalHub=-1;
    int neiID,neiDis;

    for(int i=0;i<Neighbors[ID1].size();i++){
        neiID=Neighbors[ID1][i].first;
        neiDis=Neighbors[ID1][i].second;

        for(it=Label[neiID].begin();it!=Label[neiID].end();++it){
            hub=(*it).first;
            dis1=(*it).second;
            if(hub!=ID2 && Label[ID2].find(hub)!=Label[ID2].end()){//hub is the hub of ID1 and ID2 and hub != ID2
                if(neiDis+dis1+Label[ID2][hub]<d){
                    d=neiDis+dis1+Label[ID2][hub];
                    finalHub = hub;
                }
            }
        }


    }
    return make_pair(d,finalHub);
}
//function of computing d2, new version
pair<int,set<int>> Graph::DisQueryPeak22(int ID1, int ID2,vector<unordered_map<int,int>> &Label){
	int d=INF;
	unordered_map<int,int>::iterator it;
	int hub, dis1, finalHub=-1;
    set<int> hubs;
	for(it=Label[ID1].begin();it!=Label[ID1].end();it++){
		hub=(*it).first;
		dis1=(*it).second;
		if(hub!=ID2 && Label[ID2].find(hub)!=Label[ID2].end()){//hub is the hub of ID1 and ID2 and hub != ID2
			if(dis1+Label[ID2][hub]<d){
				d=dis1+Label[ID2][hub];
                finalHub = hub;
                hubs.clear();
                hubs.insert(finalHub);
			}
//            else if(dis1+Label[ID2][hub]==d){
//                hubs.insert(hub);
//                cout<<"!!!!!!!!!!!!!!!!!!"<<endl;
//            }
		}
	}
	return make_pair(d,hubs);
}
//function of computing d2, new version
pair<int,int> Graph::DisQueryPeak2(int ID1, int ID2,vector<unordered_map<int,int>> &Label){
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
int Graph::DisQueryPeak(int ID1, int ID2,vector<unordered_map<int,int>> &Label){
    int d=INF;
    unordered_map<int,int>::iterator it;
    int hub, dis1;
    for(it=Label[ID1].begin();it!=Label[ID1].end();it++){
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
