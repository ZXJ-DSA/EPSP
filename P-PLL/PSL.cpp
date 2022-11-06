/*
 * PSL.cpp
 *
 *  Created on: 24 Oct 2022
 *      Author: zhangmengxuan
 */
#include "head.h"

void Graph::IndexConstructMThread2New(){
	unordered_map<int,int> map0; map0.clear();
	Label.assign(nodenum, map0);
	Dhop.assign(nodenum, map0);
	PruningPointNew.clear();
	unordered_map<int,vector<int>> map1;
	PruningPointNew.assign(nodenum,map1);

	PruningPointStepNew.clear();
	unordered_map<int,unordered_map<int,int>> map2;
	PruningPointStepNew.assign(nodenum, map2);
	NoSupportedPair.clear();

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
		//cout<<"step "<<step<<" finish!"<<endl;
		LabelStep.push_back(Dhop);
		//cout<<"1111111111111"<<endl;
		flag=DhopLableRefreshMulti2New(step+1);
		//cout<<"2222222222222222 "<<flag<<endl;
		step+=1;
	}
	//cout<<"Index finish construction"<<endl;
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
//cout<<"33333333333333333333"<<endl;
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
				if(Label[nodeID].find((*it).first)!=Dhop[nodeID].end()){
					Label[nodeID][(*it).first]=(*it).second;
				}else{
					Label[nodeID].insert(*it);
				}
			}
			for(int i=0;i<Neighbors[nodeID].size();i++){
				if(Dset.find(Neighbors[nodeID][i].first)==Dset.end()){
					Dset.insert(Neighbors[nodeID][i].first);
					//Dvectex.push_back(Neighbors[nodeID][i].first);
					DvertexNew[Neighbors[nodeID][i].first]=true;
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
	//sm->wait();

	for(int i=0;i<p.size();i++){
		int nodeID=p[i];
		unordered_map<int,int> Dhop0; Dhop0.clear();
		int neighID, neighW;
		for(int Index=0;Index<Neighbors[nodeID].size();Index++){
			neighID=Neighbors[nodeID][Index].first;
			//cout<<"neighID "<<neighID<<" Dhop.size "<<Dhop[neighID].size()<<endl;
			if(Dhop[neighID].size()>0){
				neighW=Neighbors[nodeID][Index].second;

				unordered_map<int,int>::iterator it=Dhop[neighID].begin();
				int hub, dis, d;
				for(;it!=Dhop[neighID].end();it++){
					hub=(*it).first; dis=(*it).second;d=neighW+dis;
					if(NodeOrder[hub]>NodeOrder[nodeID]){
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
						else{
							for(int k=0;k<SupNode.size();k++){
								int supn=SupNode[k];

								if(supn!=nodeID && supn!=hub){
									vSm[nodeID]->wait();
									PruningPointNew[nodeID][supn].push_back(hub);
									PruningPointStepNew[nodeID][supn].insert(make_pair(hub,step));

									vSm[nodeID]->notify();

									vSm[hub]->wait();
									PruningPointNew[hub][supn].push_back(nodeID);
									PruningPointStepNew[hub][supn].insert(make_pair(nodeID,step));

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

	//sm->notify();
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

void Graph::IncreasePSL(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label,vector<unordered_map<int,vector<int>>> &PruningPointNew,set<pair<int,int>> &NoSupportedPair){
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

	int dis,disvally,dispeak;
	//activate or not
	dis=DisQueryLower1(LID,HID,Neighbors,Label);
	dispeak=DisQueryPeak(LID,HID,Label);
	if(dispeak<=oldW) return;
	if(dis>oldW){//index update triggered
		vector<vector<pair<int,int>>> Change;
		vector<pair<int,int>> vec;
		Change.assign(nodenum,vec);
		set<int> WaitPro;
		vector<vector<int>> ChangeP;
		vector<int> vecint;
		ChangeP.assign(nodenum,vecint);
		set<int> WaitProP;

		WaitPro.insert(LID);
		Change[LID].push_back(make_pair(HID, oldW));
		disvally=DisQueryVally(LID,HID,Neighbors,Label);
		Label[LID][HID]=disvally;//correct to the new value

		//cout<<"start "<<LID<<" "<<HID<<" "<<disvally<<endl;

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

					//cout<<"weight "<<LID<<" "<<hubID<<" "<<disvally<<endl;
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

					//cout<<"weight "<<HID<<" "<<hubID<<" "<<disvally<<endl;
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

								//cout<<"1--1 "<<neiID<<" "<<hID<<" "<<disvally<<endl;
							}
						}
					}

					//Change->ChangeP
					if(PruningPointNew[curID].find(hID)!=PruningPointNew[curID].end()){
						for(int snum=0;snum<PruningPointNew[curID][hID].size();snum++){
							int s=PruningPointNew[curID][hID][snum];

					/*if(PruningPoint.find(make_pair(curID,hID))!=PruningPoint.end()){
						for(int snum=0;snum<PruningPoint[make_pair(curID,hID)].size();snum++){
							int s=PruningPoint[make_pair(curID,hID)][snum];*/
							if(NodeOrder[s]<NodeOrder[curID] && NoSupportedPair.find(make_pair(s,curID))==NoSupportedPair.end()){
								disvally=DisQueryVally(s,curID,Neighbors,Label);
								dispeak=DisQueryPeak(s,curID,Label);
								if(dispeak>disvally){
									WaitProPTem.insert(s);
									ChangePTem[s].push_back(curID);
									Label[s][curID]=disvally;
									NoSupportedPair.insert(make_pair(s,curID));

									//cout<<"1--2 "<<s<<" "<<curID<<" "<<disvally<<endl;
								}
							}else if(NodeOrder[s]>NodeOrder[curID] && NoSupportedPair.find(make_pair(curID,s))==NoSupportedPair.end()){
								disvally=DisQueryVally(curID,s,Neighbors,Label);
								dispeak=DisQueryPeak(curID,s,Label);
								if(dispeak>disvally){
									WaitProPTem.insert(curID);
									ChangePTem[curID].push_back(s);
									Label[curID][s]=disvally;
									NoSupportedPair.insert(make_pair(curID,s));

									//cout<<"1--2 "<<curID<<" "<<s<<" "<<disvally<<endl;
								}
							}
						}
					}
				}
			}

			//ChangeP->CHangeP
			int v,u,neiid,neiw;
			for(set<int>::iterator itp=WaitProP.begin();itp!=WaitProP.end();itp++){
				v=*itp;
				for(int k=0;k<ChangeP[v].size();k++){
					u=ChangeP[v][k];
					for(int l=0;l<Neighbors[v].size();l++){
						neiid=Neighbors[v][l].first; neiw=Neighbors[v][l].second;
						if(NodeOrder[neiid]<NodeOrder[u]){
							disvally=DisQueryVally(neiid, u,Neighbors,Label);
							dispeak=DisQueryPeak(neiid, u,Label);
							if(disvally<dispeak){
								//if(Label[neiid].find(u)==Label[neiid].end() || (Label[neiid].find(u)!=Label[neiid].end() && Label[neiid][u]>disvally)){
								if(Label[neiid].find(u)==Label[neiid].end() || (Label[neiid].find(u)!=Label[neiid].end() && Label[neiid][u]!=disvally)){
									WaitProPTem.insert(neiid);
									ChangePTem[neiid].push_back(u);
									Label[neiid][u]=disvally;
									NoSupportedPair.insert(make_pair(neiid,u));

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

void Graph::DecreasePSL(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label){
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
	int Dab=ShortestDisQuery(a,b,Label);

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

int Graph::DisQueryVally(int ID1, int ID2, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label){
	int neiID,neiDis;
	int d=INF;
	for(int i=0;i<Neighbors[ID1].size();i++){
		neiID=Neighbors[ID1][i].first;
		neiDis=Neighbors[ID1][i].second;
		if(NodeOrder[neiID]<=NodeOrder[ID2] && Label[neiID].find(ID2)!=Label[neiID].end()){
			if(neiDis+Label[neiID][ID2]<d){
				d=neiDis+Label[neiID][ID2];
			}
		}
	}
	return d;
}

int Graph::DisQueryPeak(int ID1, int ID2,vector<unordered_map<int,int>> &Label){
	int d=INF;
	unordered_map<int,int>::iterator it;
	int hub, dis1;
	for(it=Label[ID1].begin();it!=Label[ID1].end();it++){
		hub=(*it).first;
		dis1=(*it).second;
		if(hub!=ID2 && Label[ID2].find(hub)!=Label[ID2].end()){
			if(dis1+Label[ID2][hub]<d){
				d=dis1+Label[ID2][hub];
			}
		}
	}
	return d;
}

int Graph::DisQueryLower1(int ID1, int ID2, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label){
	int neiID,neiDis;
	int d=INF;
	for(int i=0;i<Neighbors[ID1].size();i++){
		neiID=Neighbors[ID1][i].first;
		neiDis=Neighbors[ID1][i].second;
		if(NodeOrder[neiID]<NodeOrder[ID2] && Label[neiID].find(ID2)!=Label[neiID].end()){
			if(neiDis+Label[neiID][ID2]<d){
				d=neiDis+Label[neiID][ID2];
			}
		}
	}
	return d;
}



