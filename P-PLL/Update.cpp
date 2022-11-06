/*
 * UpdateNew.cpp
 *
 *  Created on: 5 Oct 2022
 *      Author: zhangmengxuan
 */
#include "head.h"

void Graph::DecreaseSingle(int a, int b, int oldW, int newW){
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

	int pid=EtoParID[a][b];
	if(pid<pnum){
		vector<pair<pair<int,int>,pair<int,int>>> weightOverlay;//collect the changed edges on overlay graph
		weightOverlay.clear();
		DecreasePLL(a,b,oldW,newW,NeighborsParti[pid],LabelParti[pid]);

		//weightOverlay collect the changed edges on overlay graph
		vector<int> Bid=BoundVer[pid];
		//check the boundary edge within partition
		int bid1,bid2,olddis,newdis;
		for(int i=0;i<Bid.size();i++){
			bid1=Bid[i];
			for(int j=i+1;j<Bid.size();j++){
				bid2=Bid[j];
				olddis=BedgePID[bid1][bid2].first;
				newdis=HopQueryLocal(bid1,bid2,LabelParti[pid]);
				if(newdis<olddis){
					BedgePID[bid1][bid2].first=newdis;
					BedgePID[bid2][bid1].first=newdis;
					weightOverlay.push_back(make_pair(make_pair(bid1,bid2),make_pair(olddis,newdis)));
				}
			}
		}

		//update the overlay graph index, after partition index update
		for(int l=0;l<weightOverlay.size();l++){
			DecreasePLL(weightOverlay[l].first.first,weightOverlay[l].first.second,weightOverlay[l].second.first,weightOverlay[l].second.second,NeighborsOverlay,LabelOverlay);
		}
	}else{
		DecreasePLL(a, b, oldW, newW, NeighborsOverlay, LabelOverlay);
	}
}

void Graph::IncreaseSingle(int a, int b, int oldW, int newW){
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

	int pid=EtoParID[a][b];
	//cout<<"increase edge in partition "<<pid<<endl;
	if(pid<pnum){
		vector<pair<pair<int,int>,pair<int,int>>> weightOverlay;//collect the changed edges on overlay graph
		weightOverlay.clear();
		IncreasePLL(a,b,oldW,newW,NeighborsParti[pid],LabelParti[pid],PruningPointParti[pid],NoSupportedParti[pid]);

		//boundary edges check
		vector<int> Bid=BoundVer[pid];
		int bid1,bid2,olddis,newdis;
		for(int i=0;i<Bid.size();i++){
			bid1=Bid[i];
			for(int j=i+1;j<Bid.size();j++){
				bid2=Bid[j];
				olddis=BedgePID[bid1][bid2].first;
				newdis=HopQueryLocal(bid1,bid2,LabelParti[pid]);
				BedgePID[bid1][bid2].first=newdis;
				BedgePID[bid2][bid1].first=newdis;
				if(newdis>olddis)//if '=', not problem; if '<', problem
					weightOverlay.push_back(make_pair(make_pair(bid1,bid2),make_pair(olddis,newdis)));
				else if(newdis<olddis)
					cout<<"Something wrong happens."<<endl;
			}
		}

		//update the overlay graph index, after partition index update
		for(int l=0;l<weightOverlay.size();l++){
			IncreasePLL(weightOverlay[l].first.first,weightOverlay[l].first.second,weightOverlay[l].second.first,weightOverlay[l].second.second,NeighborsOverlay,LabelOverlay,PruningPointOverlay,NoSupportedOverlay);
		}
	}else{
		IncreasePLL(a, b, oldW, newW,NeighborsOverlay,LabelOverlay,PruningPointOverlay,NoSupportedOverlay);
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
	}
}
