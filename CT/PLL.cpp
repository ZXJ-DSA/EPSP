/*
 * PLL.cpp
 *
 *  Created on: 13 Oct 2022
 *      Author: zhangmengxuan
 */
#include "head.h"

void Graph::IndexConst1(){
	unordered_map<int,int> m;
	Label.assign(nodenum,m);
	PruningPointNew.clear();
	unordered_map<int,vector<int>> map1;
	PruningPointNew.assign(nodenum,map1);

	int ID;
	int cnt=0;
	for(int i=nodenum-1;i>=0;i--){
		ID=vNodeOrder[i];
		vector<pair<int,int>> vp;
		DijksPrune1(ID,vp);
		/*if(cnt%1000==0){
			cout<<"cnt "<<cnt<<"ID "<<ID<<" vp.size "<<vp.size()<<endl;
			//cout<<"ID "<<ID<<" vp.size "<<vp.size()<<endl;
		}*/
		cnt+=1;
		for(int j=0;j<vp.size();j++){
			Label[vp[j].first].insert(make_pair(ID, vp[j].second));
			//cout<<vp[j].first<<" "<<vp[j].second<<endl;
		}
	}
}

void Graph::DijksPrune1(int nodeID, vector<pair<int,int>>& vp){
	benchmark::heap<2, int, int> pqueue(nodenum);
	pqueue.update(nodeID,0);

	vector<bool> closed(nodenum, false);
	vector<int> distance(nodenum, INF);

	distance[nodeID]=0;
	int topNodeID, topNodeDis;
	vector<pair<int,int>>::iterator it;
	int NNodeID,NWeigh;

	while(!pqueue.empty()){
		pqueue.extract_min(topNodeID, topNodeDis);
		closed[topNodeID]=true;

		int TempDis; vector<int> SupNode;
		ShortestDisQuery1(nodeID, topNodeID,SupNode,TempDis);
		if(TempDis<=topNodeDis){

			if(NodeOrder[topNodeID]<NodeOrder[nodeID]){
				for(int k=0;k<SupNode.size();k++){
					int supn=SupNode[k];
					unordered_map<int,vector<int>> map1;
					PruningPointNew[topNodeID][supn].push_back(nodeID);
					PruningPointNew[nodeID][supn].push_back(topNodeID);
					//PruningPoint[make_pair(topNodeID,supn)].push_back(nodeID);
					//PruningPoint[make_pair(nodeID,supn)].push_back(topNodeID);
					//cout<<"pruning point "<<topNodeID<<" "<<nodeID<<" "<<supn<<endl;
				}
			}
			continue;
		}


		//Label[topNodeID].insert(nodeID, topNodeDis);
		vp.push_back(make_pair(topNodeID,topNodeDis));
		for(it=Neighbors[topNodeID].begin();it!=Neighbors[topNodeID].end();it++){
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

void Graph::IncreasePLL(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label,vector<unordered_map<int,vector<int>>> &PruningPointNew,set<pair<int,int>> &NoSupportedPair){
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

	//if(dispeak<=oldW) return;
	if(Label[LID].find(HID)!=Label[LID].end() && dis>oldW){//index update triggered
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
						for(int sum=0;sum<PruningPointNew[curID][hID].size();sum++){
							int s=PruningPointNew[curID][hID][sum];
					//if(PruningPoint.find(make_pair(curID,hID))!=PruningPoint.end()){
						//for(int snum=0;snum<PruningPoint[make_pair(curID,hID)].size();snum++){
							//int s=PruningPoint[make_pair(curID,hID)][snum];
							if(NodeOrder[s]<NodeOrder[curID] && NoSupportedPair.find(make_pair(s,curID))==NoSupportedPair.end()){
								disvally=DisQueryVally(s,curID,Neighbors,Label);
								dispeak=DisQueryPeak(s,curID,Label);
								if(dispeak>disvally){
									WaitProPTem.insert(s);
									ChangePTem[s].push_back(curID);
									Label[s][curID]=disvally;

									NoSupportedPair.insert(make_pair(s,curID));
								}
							}else if(NodeOrder[s]>NodeOrder[curID] && NoSupportedPair.find(make_pair(curID,s))==NoSupportedPair.end()){
								disvally=DisQueryVally(curID,s,Neighbors,Label);
								dispeak=DisQueryPeak(curID,s,Label);
								if(dispeak>disvally){
									WaitProPTem.insert(curID);
									ChangePTem[curID].push_back(s);
									Label[curID][s]=disvally;

									NoSupportedPair.insert(make_pair(curID,s));
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
								if(Label[neiid].find(u)==Label[neiid].end() || (Label[neiid].find(u)!=Label[neiid].end() && Label[neiid][u]!=disvally)){
									WaitProPTem.insert(neiid);
									ChangePTem[neiid].push_back(u);
									Label[neiid][u]=disvally;

									NoSupportedPair.insert(make_pair(neiid,u));
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
        extUpdate = true;
	}
}

void Graph::DecreasePLL(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors,vector<unordered_map<int,int>> &Label){
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
        /// update flag
        extUpdate = true;
	}
}
