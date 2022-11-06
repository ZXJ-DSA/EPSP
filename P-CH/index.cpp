/*
 * index.cpp
 *
 *  Created on: 19 Oct 2022
 *      Author: zhangmengxuan
 */
#include "head.h"

vector<int> NodeOrders;
struct OrderComp{
	int x;
	int y;//order(x)<order(y)
	OrderComp(int _x, int _y){
		x=_x; y=_y;
	}
	bool operator< (const OrderComp& d) const{
		if(x==d.x && y==d.y){//avoid the redundant
			return false;
		}else{
			if(x!=d.x)
				return NodeOrders[x]<NodeOrders[d.x];
			if(y!=d.y)
				return NodeOrders[y]<NodeOrders[d.y];
		}
	}
};
/*void Graph::CHindex(){
	//initialize E
	map<int,pair<int,int>> m;
	E.assign(nodenum,m);
	for(int i=0;i<Neighbor.size();i++){
		for(int j=0;j<Neighbor[i].size();j++)
			E[i].insert(make_pair(Neighbor[i][j].first,make_pair(Neighbor[i][j].second,1)));
	}

	vector<bool> exist; exist.assign(nodenum,true);
	vector<bool> change; change.assign(nodenum,false);

	//cout<<"Begin to contract"<<endl;
	SCconNodes.clear();
	for(int nodeorder=0;nodeorder<nodenum;nodeorder++){
		int x=vNodeOrder[nodeorder];
		//if(x!=-1){//to identify and exclude the isolated vertices
		if(Neighbor[x].size()!=0){//to identify and exclude the isolated vertices
			exist[x]=false;

			vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();

			for(auto it=E[x].begin();it!=E[x].end();it++){
				if(exist[(*it).first]){
					Neigh.push_back(*it);
				}
			}
			NeighborCon[x].assign(Neigh.begin(),Neigh.end());

			for(int i=0;i<Neigh.size();i++){
				int y=Neigh[i].first;
				deleteEorder(x,y);
			}

			for(int i=0;i<Neigh.size();i++){
				for(int j=i+1;j<Neigh.size();j++){
					insertEorder(Neigh[i].first,Neigh[j].first,Neigh[i].second.first+Neigh[j].second.first);
					if(Neigh[i].first<Neigh[j].first)
						SCconNodes[make_pair(Neigh[i].first,Neigh[j].first)].push_back(x);//no direction
					else if(Neigh[j].first<Neigh[i].first)
						SCconNodes[make_pair(Neigh[j].first,Neigh[i].first)].push_back(x);
				}
			}
		}
	}
}*/

void Graph::CHindexMT(){
	map<int, vector<int>> mi;
	SCconNodesMT.assign(nodenum, mi);
	vector<pair<int,pair<int,int>>> vect;
	NeighborCon.assign(nodenum,vect);

	//initialize E
	map<int,pair<int,int>> m;
	E.assign(nodenum,m);
	for(int i=0;i<Neighbor.size();i++){
		for(int j=0;j<Neighbor[i].size();j++)
			E[i].insert(make_pair(Neighbor[i][j].first,make_pair(Neighbor[i][j].second,1)));
	}

	vector<bool> exist; exist.assign(nodenum,true);
	//vector<bool> change; change.assign(nodenum,false);

	//cout<<"Begin to contract"<<endl;
	for(int nodeorder=0;nodeorder<nodenum;nodeorder++){
		int x=vNodeOrder[nodeorder];
		if(Neighbor[x].size()>0){
		//if(x!=-1){//to identify and exclude the isolated vertices
			exist[x]=false;

			vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();

			for(auto it=E[x].begin();it!=E[x].end();it++){
				if(exist[(*it).first]){
					Neigh.push_back(*it);
				}
			}
			NeighborCon[x].assign(Neigh.begin(),Neigh.end());

			for(int i=0;i<Neigh.size();i++){
				int y=Neigh[i].first;
				deleteEorder(x,y);
				//change[y]=true;
			}

			if(Neigh.size()<=100){
				//single thread
				for(int i=0;i<Neigh.size();i++){
					for(int j=i+1;j<Neigh.size();j++){
						insertEorder(Neigh[i].first,Neigh[j].first,Neigh[i].second.first+Neigh[j].second.first);
						if(Neigh[i].first<Neigh[j].first)
							SCconNodesMT[Neigh[i].first][Neigh[j].first].push_back(x);//no direction
						else if(Neigh[j].first<Neigh[i].first)
							SCconNodesMT[Neigh[j].first][Neigh[i].first].push_back(x);
					}
				}
			}else{
				//multiple thread
				int step=Neigh.size()/threadnum;
				boost::thread_group thread;
				for(int i=0;i<threadnum;i++){
					pair<int,int> p;
					p.first=i*step;
					if(i==threadnum-1)
						p.second=Neigh.size();
					else
						p.second=(i+1)*step;
					thread.add_thread(new boost::thread(&Graph::NeighborComorder, this, boost::ref(Neigh), p, x));
				}
				thread.join_all();
			}

		}
	}
}

void Graph::NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x){
//	sm->wait();
	int ID1, w1;
	int ID2, w2;
	for(int k=p.first;k<p.second;k++){
		ID1=Neighvec[k].first;
		w1=Neighvec[k].second.first;
		for(int h=0;h<Neighvec.size();h++){
			ID2=Neighvec[h].first;
			w2=Neighvec[h].second.first;
			insertEMTorder(ID1, ID2, w1+w2);
			if(ID1<ID2)
				SCconNodesMT[ID1][ID2].push_back(x);
		}
	}
//	sm->notify();
}

void Graph::insertEMTorder(int u,int v,int w){
	if(E[u].find(v)==E[u].end()){
		E[u].insert(make_pair(v,make_pair(w,1)));
	}
	else{
		if(E[u][v].first>w)
			E[u][v]=make_pair(w,1);
		else if(E[u][v].first==w)
			E[u][v].second+=1;
	}
}

void Graph::insertEorder(int u,int v,int w){
	if(E[u].find(v)==E[u].end()){
		E[u].insert(make_pair(v,make_pair(w,1)));
	}
	else{
		if(E[u][v].first>w)
			E[u][v]=make_pair(w,1);
		else if(E[u][v].first==w)
			E[u][v].second+=1;
	}

	if(E[v].find(u)==E[v].end()){
		E[v].insert(make_pair(u,make_pair(w,1)));
	}
	else{
		if(E[v][u].first>w)
			E[v][u]=make_pair(w,1);
		else if(E[v][u].first==w)
			E[v][u].second+=1;
	}
}

void Graph::deleteEorder(int u,int v){
	if(E[u].find(v)!=E[u].end()){
		E[u].erase(E[u].find(v));
	}

	if(E[v].find(u)!=E[v].end()){
		E[v].erase(E[v].find(u));
	}
}

//Query processing
int	Graph::QueryCH(int ID1, int ID2, vector<vector<pair<int,pair<int,int>>>> &NeighborCon){
	if(ID1==ID2) return 0;
	if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
	int d=INF;
	benchmark::heap<2,int,int> fHeapForward(nodenum);
	benchmark::heap<2, int, int> fHeapBackward(nodenum);

	//closed or not
	vector<bool> vVisitedF(nodenum, false);
	vector<bool> vVisitedB(nodenum, false);
	//the existing shortest distance
	vector<int>	vDistanceForward(nodenum, INF);
	vector<int>	vDistanceBackward(nodenum, INF);
	//stop search or not
	bool bF = false;
	bool bB = false;
	vDistanceForward[ID1] = 0;
	vDistanceBackward[ID2] = 0;
	fHeapForward.update(ID1,0);
	fHeapBackward.update(ID2,0);

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

			for(auto out=NeighborCon[topNodeIDForward].begin();out!=NeighborCon[topNodeIDForward].end();out++){
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

			for(auto in=NeighborCon[topNodeIDBackward].begin();in!=NeighborCon[topNodeIDBackward].end();in++){
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

int	Graph::QueryCHInIn(int ID1, int ID2, int PID){
	if(ID1==ID2) return 0;
	if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
	int d=INF;
	benchmark::heap<2,int,int> fHeapForward(nodenum);
	benchmark::heap<2, int, int> fHeapBackward(nodenum);

	//closed or not
	vector<bool> vVisitedF(nodenum, false);
	vector<bool> vVisitedB(nodenum, false);
	//the existing shortest distance
	vector<int>	vDistanceForward(nodenum, INF);
	vector<int>	vDistanceBackward(nodenum, INF);
	//stop search or not
	bool bF = false;
	bool bB = false;
	vDistanceForward[ID1] = 0;
	vDistanceBackward[ID2] = 0;
	fHeapForward.update(ID1,0);
	fHeapBackward.update(ID2,0);

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

			for(auto out=NeighborCons[PID][topNodeIDForward].begin();out!=NeighborCons[PID][topNodeIDForward].end();out++){
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

			for(auto out=NeighborConOverlay[topNodeIDForward].begin();out!=NeighborConOverlay[topNodeIDForward].end();out++){
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

			for(auto in=NeighborCons[PID][topNodeIDBackward].begin();in!=NeighborCons[PID][topNodeIDBackward].end();in++){
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

			for(auto in=NeighborConOverlay[topNodeIDBackward].begin();in!=NeighborConOverlay[topNodeIDBackward].end();in++){
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

int	Graph::QueryCHInOut(int ID1, int ID2, int PID1){
	if(ID1==ID2) return 0;
	if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
	int d=INF;
	benchmark::heap<2,int,int> fHeapForward(nodenum);
	benchmark::heap<2, int, int> fHeapBackward(nodenum);

	//closed or not
	vector<bool> vVisitedF(nodenum, false);
	vector<bool> vVisitedB(nodenum, false);
	//the existing shortest distance
	vector<int>	vDistanceForward(nodenum, INF);
	vector<int>	vDistanceBackward(nodenum, INF);
	//stop search or not
	bool bF = false;
	bool bB = false;
	vDistanceForward[ID1] = 0;
	vDistanceBackward[ID2] = 0;
	fHeapForward.update(ID1,0);
	fHeapBackward.update(ID2,0);

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
				for(auto out=NeighborCons[PID1][topNodeIDForward].begin();out!=NeighborCons[PID1][topNodeIDForward].end();out++){
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
				for(auto out=NeighborConOverlay[topNodeIDForward].begin();out!=NeighborConOverlay[topNodeIDForward].end();out++){
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

			for(auto in=NeighborConOverlay[topNodeIDBackward].begin();in!=NeighborConOverlay[topNodeIDBackward].end();in++){
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

int	Graph::QueryCHOutOut(int ID1, int ID2, int PID1, int PID2){
	if(ID1==ID2) return 0;
	if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
	int d=INF;
	benchmark::heap<2,int,int> fHeapForward(nodenum);
	benchmark::heap<2, int, int> fHeapBackward(nodenum);

	//closed or not
	vector<bool> vVisitedF(nodenum, false);
	vector<bool> vVisitedB(nodenum, false);
	//the existing shortest distance
	vector<int>	vDistanceForward(nodenum, INF);
	vector<int>	vDistanceBackward(nodenum, INF);
	//stop search or not
	bool bF = false;
	bool bB = false;
	vDistanceForward[ID1] = 0;
	vDistanceBackward[ID2] = 0;
	fHeapForward.update(ID1,0);
	fHeapBackward.update(ID2,0);

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

			for(auto out=NeighborCons[PID1][topNodeIDForward].begin();out!=NeighborCons[PID1][topNodeIDForward].end();out++){
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

			for(auto out=NeighborConOverlay[topNodeIDForward].begin();out!=NeighborConOverlay[topNodeIDForward].end();out++){
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

			for(auto in=NeighborCons[PID2][topNodeIDBackward].begin();in!=NeighborCons[PID2][topNodeIDBackward].end();in++){
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

			for(auto in=NeighborConOverlay[topNodeIDBackward].begin();in!=NeighborConOverlay[topNodeIDBackward].end();in++){
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

//CH index update
void Graph::CHdec(int a,int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbor,vector<vector<pair<int,pair<int,int>>>> &NeighborCon){
	//modify the original graph information
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
	//NodeOrders.clear();
	NodeOrders.assign(NodeOrder.begin(),NodeOrder.end());
	set<OrderComp> OC;
	map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the fresh distance and avoid search in the adjacent list
	//OC.clear(); OCdis.clear();

	if(NodeOrder[a]<NodeOrder[b]){
		for(int i=0;i<NeighborCon[a].size();i++){
			if(NeighborCon[a][i].first==b){
				if(NeighborCon[a][i].second.first>newW){
					NeighborCon[a][i].second.first=newW;
					NeighborCon[a][i].second.second=1;

					OCdis[make_pair(a,b)]=newW;
					OC.insert(OrderComp(a,b));
				}else if(NeighborCon[a][i].second.first==newW)
					NeighborCon[a][i].second.second+=1;
				break;
			}
		}
	}else{
		for(int i=0;i<NeighborCon[b].size();i++){
			if(NeighborCon[b][i].first==a){
				if(NeighborCon[b][i].second.first>newW){
					NeighborCon[b][i].second.first=newW;
					NeighborCon[b][i].second.second=1;

					OCdis[make_pair(b,a)]=newW;
					OC.insert(OrderComp(b,a));
				}else if(NeighborCon[b][i].second.first==newW)
					NeighborCon[b][i].second.second+=1;
				break;
			}
		}
	}


	while(!OC.empty()){
		int s=(*OC.begin()).x; int t=(*OC.begin()).y;
		int wt;
		OC.erase(OC.begin());

		wt=OCdis[make_pair(s,t)];
		map<int,int> InM2t; //InM2t.clear();
		vector<pair<int,int>> InMLower; //InMLower.clear();
		for(int i=0;i<NeighborCon[s].size();i++){
			if(NodeOrder[NeighborCon[s][i].first]>NodeOrder[t])
				InM2t.insert(make_pair(NeighborCon[s][i].first,NeighborCon[s][i].second.first));
			else if(NodeOrder[NeighborCon[s][i].first]<NodeOrder[t])
				InMLower.push_back(make_pair(NeighborCon[s][i].first,NeighborCon[s][i].second.first));
		}
		int inID,inW,inWt;
		for(int i=0;i<NeighborCon[t].size();i++){
			inID=NeighborCon[t][i].first;
			if(InM2t.find(inID)!=InM2t.end()){
				inW=InM2t[inID];
				inWt=NeighborCon[t][i].second.first;
				if(inWt>inW+wt){
					NeighborCon[t][i].second.first=inW+wt;
					NeighborCon[t][i].second.second=1;
					OCdis[make_pair(t,inID)]=inW+wt;
					OrderComp oc={t,inID};
					OC.insert(oc);
				}else if(inWt==inW+wt){
					NeighborCon[t][i].second.second+=1;
				}
			}
		}

		for(int i=0;i<InMLower.size();i++){
			inID=InMLower[i].first; inW=InMLower[i].second;
			for(int j=0;j<NeighborCon[inID].size();j++){
				if(NeighborCon[inID][j].first==t){
					inWt=NeighborCon[inID][j].second.first;
					if(inWt>inW+wt){
						NeighborCon[inID][j].second.first=inW+wt;
						NeighborCon[inID][j].second.second=1;

						OCdis[make_pair(inID,t)]=inW+wt;
						OrderComp oc={inID,t};
						OC.insert(oc);
					}else if(inWt==inW+wt)
						NeighborCon[inID][j].second.second+=1;
					break;
				}
			}
		}
	}//finish change index
}

void Graph::CHinc(int a,int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbor,vector<vector<pair<int,pair<int,int>>>> &NeighborCon,vector<map<int, vector<int>>> &SCconNodesMT){
	//modify the original graph information
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

	//NodeOrders.clear();
		NodeOrders.assign(NodeOrder.begin(),NodeOrder.end());
		set<OrderComp> OC; //OC.clear();
		map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
		//OCdis.clear();

		if(NodeOrder[a]<NodeOrder[b]){
			for(int i=0;i<NeighborCon[a].size();i++){
				if(NeighborCon[a][i].first==b){
					if(NeighborCon[a][i].second.first==oldW){
						NeighborCon[a][i].second.second-=1;
						if(NeighborCon[a][i].second.second<1){
							OrderComp oc={a,b};
							OC.insert(oc);
							OCdis[make_pair(a,b)]=oldW;
						}
					}
					break;
				}
			}
		}else{
			for(int i=0;i<NeighborCon[b].size();i++){
				if(NeighborCon[b][i].first==a){
					if(NeighborCon[b][i].second.first==oldW){
						NeighborCon[b][i].second.second-=1;
						if(NeighborCon[b][i].second.second<1){
							OrderComp oc={b,a};
							OC.insert(oc);
							OCdis[make_pair(b,a)]=oldW;
						}
					}
					break;
				}
			}
		}


		while(!OC.empty()){
		int s=(*OC.begin()).x; int t=(*OC.begin()).y;
		int wt;
		OC.erase(OC.begin());

		wt=OCdis[make_pair(s,t)];//distance of s--->t before change
		int inID,inW;
		map<int,int> HigherIn; vector<pair<int,int>> LowerIn;
		//HigherIn.clear(); LowerIn.clear();
		//the shortcuts infected by s-->t
		for(int i=0;i<NeighborCon[s].size();i++){
			inID=NeighborCon[s][i].first;
			inW=NeighborCon[s][i].second.first;
			if(NodeOrder[inID]<NodeOrder[t]){
				LowerIn.push_back(make_pair(inID,inW));
			}else if(NodeOrder[inID]>NodeOrder[t]){
				HigherIn.insert(make_pair(inID,inW));
			}
		}
		for(int i=0;i<NeighborCon[t].size();i++){
			inID=NeighborCon[t][i].first;
			if(HigherIn.find(inID)!=HigherIn.end()){
				inW=HigherIn[inID];
				if(NeighborCon[t][i].second.first==wt+inW){
					NeighborCon[t][i].second.second-=1;
					if(NeighborCon[t][i].second.second<1){
						OrderComp oc={t,inID};
						OC.insert(oc);
						OCdis[make_pair(t,inID)]=wt+inW;
					}
				}
			}
		}
		for(int i=0;i<LowerIn.size();i++){
			inID=LowerIn[i].first; inW=LowerIn[i].second;
			for(int j=0;j<NeighborCon[inID].size();j++){
				if(NeighborCon[inID][j].first==t){
					if(NeighborCon[inID][j].second.first==inW+wt){
						NeighborCon[inID][j].second.second-=1;
						if(NeighborCon[inID][j].second.second<1){
							OrderComp oc={inID,t};
							OC.insert(oc);
							OCdis[make_pair(inID,t)]=wt+inW;
						}
					}
					break;
				}
			}
		}

		//get the new weight value of s-->t
		wt=INF; int countwt=0;
		for(int i=0;i<Neighbor[s].size();i++){
			if(Neighbor[s][i].first==t){
				wt=Neighbor[s][i].second;//the weight value in the original graph
				countwt=1;
				break;
			}
		}
		int ssw,wtt,wid;
		vector<int> Wnodes; //Wnodes.clear();
		if(s<t){
			//Wnodes=SCconNodes[make_pair(s,t)]; //cout<<"wid num "<<Wnodes.size()<<endl;
			Wnodes=SCconNodesMT[s][t];
		}else{
			//Wnodes=SCconNodes[make_pair(t,s)];
			Wnodes=SCconNodesMT[t][s];
		}

		for(int i=0;i<Wnodes.size();i++){
			wid=Wnodes[i];
			for(int j=0;j<NeighborCon[wid].size();j++){
				if(NeighborCon[wid][j].first==s){
					ssw=NeighborCon[wid][j].second.first;
				}
				if(NeighborCon[wid][j].first==t){
					wtt=NeighborCon[wid][j].second.first;
				}
			}

			if(ssw+wtt<wt){
				wt=ssw+wtt;
				countwt=1;
			}else if(ssw+wtt==wt){
				countwt+=1;
			}
		}

		//refresh the weight value of s--t in the index
		for(int i=0;i<NeighborCon[s].size();i++){
			if(NeighborCon[s][i].first==t){
				NeighborCon[s][i].second.first=wt;
				NeighborCon[s][i].second.second=countwt;
				break;
			}
		}
	}
}

void Graph::CHdecBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbor,vector<vector<pair<int,pair<int,int>>>> &NeighborCon){
	//maintain the index caused by the weight change
	//NodeOrders.clear();
	NodeOrders.assign(NodeOrder.begin(),NodeOrder.end());
	set<OrderComp> OC;
	map<pair<int,int>,int> OCdis;//{(s,t),d} maintain the fresh distance and avoid search in the adjacent list
	//OC.clear(); OCdis.clear();

	int a,b,newW;//the weight of (a,b) decrease to newW
	for(int k=0;k<wBatch.size();k++){
		a=wBatch[k].first.first;
		b=wBatch[k].first.second;
		newW=wBatch[k].second.second;

		//modify the information in original graph
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

		if(NodeOrder[a]<NodeOrder[b]){
			for(int i=0;i<NeighborCon[a].size();i++){
				if(NeighborCon[a][i].first==b){
					if(NeighborCon[a][i].second.first>newW){
						//cout<<OutNeighborCon[a][i].second.first<<"..........."<<newW<<endl;
						NeighborCon[a][i].second.first=newW;
						NeighborCon[a][i].second.second=1;

						OCdis[make_pair(a,b)]=newW;
						OC.insert(OrderComp(a,b));
					}else if(NeighborCon[a][i].second.first==newW)
						NeighborCon[a][i].second.second+=1;
					break;
				}
			}
		}else{
			for(int i=0;i<NeighborCon[b].size();i++){
				if(NeighborCon[b][i].first==a){
					if(NeighborCon[b][i].second.first>newW){
						NeighborCon[b][i].second.first=newW;
						NeighborCon[b][i].second.second=1;

						OCdis[make_pair(b,a)]=newW;
						OC.insert(OrderComp(b,a));
					}else if(NeighborCon[b][i].second.first==newW)
						NeighborCon[b][i].second.second+=1;
					break;
				}
			}
		}
	}


	while(!OC.empty()){
		int s=(*OC.begin()).x; int t=(*OC.begin()).y;
		int wt;
		OC.erase(OC.begin());
		wt=OCdis[make_pair(s,t)];
		map<int,int> InM2t; //InM2t.clear();
		vector<pair<int,int>> InMLower; //InMLower.clear();
		for(int i=0;i<NeighborCon[s].size();i++){
			if(NodeOrder[NeighborCon[s][i].first]>NodeOrder[t])
				InM2t.insert(make_pair(NeighborCon[s][i].first,NeighborCon[s][i].second.first));
			else if(NodeOrder[NeighborCon[s][i].first]<NodeOrder[t])
				InMLower.push_back(make_pair(NeighborCon[s][i].first,NeighborCon[s][i].second.first));
		}
		int inID,inW,inWt;
		for(int i=0;i<NeighborCon[t].size();i++){
			inID=NeighborCon[t][i].first;
			if(InM2t.find(inID)!=InM2t.end()){
				inW=InM2t[inID];
				inWt=NeighborCon[t][i].second.first;
				if(inWt>inW+wt){
					NeighborCon[t][i].second.first=inW+wt;
					NeighborCon[t][i].second.second=1;
					OCdis[make_pair(t,inID)]=inW+wt;
					OrderComp oc={t,inID};
					OC.insert(oc);
				}else if(inWt==inW+wt){
					NeighborCon[t][i].second.second+=1;
				}
			}
		}

		for(int i=0;i<InMLower.size();i++){
			inID=InMLower[i].first; inW=InMLower[i].second;
			for(int j=0;j<NeighborCon[inID].size();j++){
				if(NeighborCon[inID][j].first==t){
					inWt=NeighborCon[inID][j].second.first;
					if(inWt>inW+wt){
						NeighborCon[inID][j].second.first=inW+wt;
						NeighborCon[inID][j].second.second=1;

						OCdis[make_pair(inID,t)]=inW+wt;
						OrderComp oc={inID,t};
						OC.insert(oc);
					}else if(inWt==inW+wt)
						NeighborCon[inID][j].second.second+=1;
					break;
				}
			}
		}
	}//finish change index
}

void Graph::CHincBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbor,vector<vector<pair<int,pair<int,int>>>> &NeighborCon,vector<map<int, vector<int>>> &SCconNodesMT){
	//NodeOrders.clear();
	NodeOrders.assign(NodeOrder.begin(),NodeOrder.end());
	set<OrderComp> OC; //OC.clear();
	map<pair<int,int>,int> OCdis;//{(s,t),d} maintain the old distance before refreshed and avoid search in the adjacent list
	//OCdis.clear();

	for(int wb=0;wb<wBatch.size();wb++){
		int a=wBatch[wb].first.first;
		int b=wBatch[wb].first.second;
		int oldW=wBatch[wb].second.first;
		int newW=wBatch[wb].second.second;

		//modify the original graph information
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


		if(NodeOrder[a]<NodeOrder[b]){
			for(int i=0;i<NeighborCon[a].size();i++){
				if(NeighborCon[a][i].first==b){
					if(NeighborCon[a][i].second.first==oldW){
						NeighborCon[a][i].second.second-=1;
						if(NeighborCon[a][i].second.second<1){
							OrderComp oc={a,b};
							OC.insert(oc);
							OCdis[make_pair(a,b)]=oldW;
						}
					}
					break;
				}
			}
		}else{
			for(int i=0;i<NeighborCon[b].size();i++){
				if(NeighborCon[b][i].first==a){
					if(NeighborCon[b][i].second.first==oldW){
						NeighborCon[b][i].second.second-=1;
						if(NeighborCon[b][i].second.second<1){
							OrderComp oc={b,a};
							OC.insert(oc);
							OCdis[make_pair(b,a)]=oldW;
						}
					}
					break;
				}
			}
		}
	}

	while(!OC.empty()){
		int s=(*OC.begin()).x; int t=(*OC.begin()).y;
		int wt;
		OC.erase(OC.begin());
		wt=OCdis[make_pair(s,t)];//distance of s--->t before change
		int inID,inW;
		map<int,int> HigherIn; vector<pair<int,int>> LowerIn;
		//HigherIn.clear(); LowerIn.clear();
		//the shortcuts infected by s-->t
		for(int i=0;i<NeighborCon[s].size();i++){
			inID=NeighborCon[s][i].first;
			inW=NeighborCon[s][i].second.first;
			if(NodeOrder[inID]<NodeOrder[t]){
				LowerIn.push_back(make_pair(inID,inW));
			}else if(NodeOrder[inID]>NodeOrder[t]){
				HigherIn.insert(make_pair(inID,inW));
			}
		}
		for(int i=0;i<NeighborCon[t].size();i++){
			inID=NeighborCon[t][i].first;
			if(HigherIn.find(inID)!=HigherIn.end()){
				inW=HigherIn[inID];
				if(NeighborCon[t][i].second.first==wt+inW){
					NeighborCon[t][i].second.second-=1;
					if(NeighborCon[t][i].second.second<1){
						OrderComp oc={t,inID};
						OC.insert(oc);
						OCdis[make_pair(t,inID)]=wt+inW;
					}
				}
			}
		}
		for(int i=0;i<LowerIn.size();i++){
			inID=LowerIn[i].first; inW=LowerIn[i].second;
			for(int j=0;j<NeighborCon[inID].size();j++){
				if(NeighborCon[inID][j].first==t){
					if(NeighborCon[inID][j].second.first==inW+wt){
						NeighborCon[inID][j].second.second-=1;
						if(NeighborCon[inID][j].second.second<1){
							OrderComp oc={inID,t};
							OC.insert(oc);
							OCdis[make_pair(inID,t)]=wt+inW;
						}
					}
					break;
				}
			}
		}

		//get the new weight value of s-->t
		wt=INF; int countwt=0;
		for(int i=0;i<Neighbor[s].size();i++){
			if(Neighbor[s][i].first==t){
				wt=Neighbor[s][i].second;//the weight value in the original graph
				countwt=1;
				break;
			}
		}
		int ssw,wtt,wid;
		vector<int> Wnodes; //Wnodes.clear();
		if(s<t){
			//Wnodes=SCconNodes[make_pair(s,t)]; //cout<<"wid num "<<Wnodes.size()<<endl;
			Wnodes=SCconNodesMT[s][t];
		}else{
			//Wnodes=SCconNodes[make_pair(t,s)];
			Wnodes=SCconNodesMT[s][t];
		}

		for(int i=0;i<Wnodes.size();i++){
			wid=Wnodes[i];
			for(int j=0;j<NeighborCon[wid].size();j++){
				if(NeighborCon[wid][j].first==s){
					ssw=NeighborCon[wid][j].second.first;
				}
				if(NeighborCon[wid][j].first==t){
					wtt=NeighborCon[wid][j].second.first;
				}
			}

			if(ssw+wtt<wt){
				wt=ssw+wtt;
				countwt=1;
			}else if(ssw+wtt==wt){
				countwt+=1;
			}
		}

		//refresh the weight value of s--t in the index
		for(int i=0;i<NeighborCon[s].size();i++){
			if(NeighborCon[s][i].first==t){
				NeighborCon[s][i].second.first=wt;
				NeighborCon[s][i].second.second=countwt;
				break;
			}
		}
	}
}

//CH index update with boundary shortcuts recorded
vector<pair<pair<int,int>,pair<int,int>>> Graph::CHdecNew(int a,int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbor,vector<vector<pair<int,pair<int,int>>>> &NeighborCon){
	vector<pair<pair<int,int>,pair<int,int>>> totalupdate;//record the updated shortcuts between boundary vertex
	totalupdate.clear();
	map<pair<int,int>,int> PreW;
	map<pair<int,int>,int> NowW;
	PreW.clear();
	NowW.clear();

	//modify the original graph information
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
	//NodeOrders.clear();
	NodeOrders.assign(NodeOrder.begin(),NodeOrder.end());
	set<OrderComp> OC;
	map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the fresh distance and avoid search in the adjacent list
	//OC.clear(); OCdis.clear();

	if(NodeOrder[a]<NodeOrder[b]){
		for(int i=0;i<NeighborCon[a].size();i++){
			if(NeighborCon[a][i].first==b){
				if(NeighborCon[a][i].second.first>newW){

					if(TotalBoundSet.find(a)!=TotalBoundSet.end() && TotalBoundSet.find(b)!=TotalBoundSet.end()){
						if(PreW.find(make_pair(a,b))==PreW.end()){
							PreW.insert(make_pair(make_pair(a,b),NeighborCon[a][i].second.first));
						}
						NowW[make_pair(a,b)]=newW;
					}

					NeighborCon[a][i].second.first=newW;
					NeighborCon[a][i].second.second=1;

					OCdis[make_pair(a,b)]=newW;
					OC.insert(OrderComp(a,b));
				}else if(NeighborCon[a][i].second.first==newW)
					NeighborCon[a][i].second.second+=1;
				break;
			}
		}
	}else{
		for(int i=0;i<NeighborCon[b].size();i++){
			if(NeighborCon[b][i].first==a){
				if(NeighborCon[b][i].second.first>newW){

					if(TotalBoundSet.find(a)!=TotalBoundSet.end() && TotalBoundSet.find(b)!=TotalBoundSet.end()){
						if(PreW.find(make_pair(b,a))==PreW.end()){
							PreW.insert(make_pair(make_pair(b,a),NeighborCon[b][i].second.first));
						}
						NowW[make_pair(b,a)]=newW;
					}

					NeighborCon[b][i].second.first=newW;
					NeighborCon[b][i].second.second=1;

					OCdis[make_pair(b,a)]=newW;
					OC.insert(OrderComp(b,a));
				}else if(NeighborCon[b][i].second.first==newW)
					NeighborCon[b][i].second.second+=1;
				break;
			}
		}
	}


	while(!OC.empty()){
		int s=(*OC.begin()).x; int t=(*OC.begin()).y;
		int wt;
		OC.erase(OC.begin());

		wt=OCdis[make_pair(s,t)];
		map<int,int> InM2t; //InM2t.clear();
		vector<pair<int,int>> InMLower; //InMLower.clear();
		for(int i=0;i<NeighborCon[s].size();i++){
			if(NodeOrder[NeighborCon[s][i].first]>NodeOrder[t])
				InM2t.insert(make_pair(NeighborCon[s][i].first,NeighborCon[s][i].second.first));
			else if(NodeOrder[NeighborCon[s][i].first]<NodeOrder[t])
				InMLower.push_back(make_pair(NeighborCon[s][i].first,NeighborCon[s][i].second.first));
		}
		int inID,inW,inWt;
		for(int i=0;i<NeighborCon[t].size();i++){
			inID=NeighborCon[t][i].first;
			if(InM2t.find(inID)!=InM2t.end()){
				inW=InM2t[inID];
				inWt=NeighborCon[t][i].second.first;
				if(inWt>inW+wt){

					if(TotalBoundSet.find(t)!=TotalBoundSet.end() && TotalBoundSet.find(inID)!=TotalBoundSet.end()){
						if(PreW.find(make_pair(t,inID))==PreW.end()){
							PreW.insert(make_pair(make_pair(t,inID),inWt));
						}
						NowW[make_pair(t,inID)]=inW+wt;
					}

					NeighborCon[t][i].second.first=inW+wt;
					NeighborCon[t][i].second.second=1;
					OCdis[make_pair(t,inID)]=inW+wt;
					OrderComp oc={t,inID};
					OC.insert(oc);
				}else if(inWt==inW+wt){
					NeighborCon[t][i].second.second+=1;
				}
			}
		}

		for(int i=0;i<InMLower.size();i++){
			inID=InMLower[i].first; inW=InMLower[i].second;
			for(int j=0;j<NeighborCon[inID].size();j++){
				if(NeighborCon[inID][j].first==t){
					inWt=NeighborCon[inID][j].second.first;
					if(inWt>inW+wt){

						if(TotalBoundSet.find(inID)!=TotalBoundSet.end() && TotalBoundSet.find(t)!=TotalBoundSet.end()){
							if(PreW.find(make_pair(inID,t))==PreW.end()){
								PreW.insert(make_pair(make_pair(inID,t),inWt));
							}
							NowW[make_pair(inID,t)]=inW+wt;
						}

						NeighborCon[inID][j].second.first=inW+wt;
						NeighborCon[inID][j].second.second=1;

						OCdis[make_pair(inID,t)]=inW+wt;
						OrderComp oc={inID,t};
						OC.insert(oc);
					}else if(inWt==inW+wt)
						NeighborCon[inID][j].second.second+=1;
					break;
				}
			}
		}
	}//finish change index

	for(auto it=PreW.begin();it!=PreW.end();it++){
		totalupdate.push_back(make_pair((*it).first,make_pair((*it).second,NowW[(*it).first])));
	}
	return totalupdate;
}

vector<pair<pair<int,int>,pair<int,int>>> Graph::CHincNew(int a,int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbor,vector<vector<pair<int,pair<int,int>>>> &NeighborCon,vector<map<int, vector<int>>> &SCconNodesMT){
	vector<pair<pair<int,int>,pair<int,int>>> totalupdate;//record the updated shortcuts between boundary vertex
	totalupdate.clear();
	map<pair<int,int>,int> PreW;
	map<pair<int,int>,int> NowW;
	PreW.clear();
	NowW.clear();

	//modify the original graph information
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

	//NodeOrders.clear();
		NodeOrders.assign(NodeOrder.begin(),NodeOrder.end());
		set<OrderComp> OC; //OC.clear();
		map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
		//OCdis.clear();

		if(NodeOrder[a]<NodeOrder[b]){
			for(int i=0;i<NeighborCon[a].size();i++){
				if(NeighborCon[a][i].first==b){
					if(NeighborCon[a][i].second.first==oldW){
						NeighborCon[a][i].second.second-=1;
						if(NeighborCon[a][i].second.second<1){
							OrderComp oc={a,b};
							OC.insert(oc);
							OCdis[make_pair(a,b)]=oldW;
						}
					}
					break;
				}
			}
		}else{
			for(int i=0;i<NeighborCon[b].size();i++){
				if(NeighborCon[b][i].first==a){
					if(NeighborCon[b][i].second.first==oldW){
						NeighborCon[b][i].second.second-=1;
						if(NeighborCon[b][i].second.second<1){
							OrderComp oc={b,a};
							OC.insert(oc);
							OCdis[make_pair(b,a)]=oldW;
						}
					}
					break;
				}
			}
		}


		while(!OC.empty()){
		int s=(*OC.begin()).x; int t=(*OC.begin()).y;
		int wt;
		OC.erase(OC.begin());

		wt=OCdis[make_pair(s,t)];//distance of s--->t before change
		int inID,inW;
		map<int,int> HigherIn; vector<pair<int,int>> LowerIn;
		//HigherIn.clear(); LowerIn.clear();
		//the shortcuts infected by s-->t
		for(int i=0;i<NeighborCon[s].size();i++){
			inID=NeighborCon[s][i].first;
			inW=NeighborCon[s][i].second.first;
			if(NodeOrder[inID]<NodeOrder[t]){
				LowerIn.push_back(make_pair(inID,inW));
			}else if(NodeOrder[inID]>NodeOrder[t]){
				HigherIn.insert(make_pair(inID,inW));
			}
		}
		for(int i=0;i<NeighborCon[t].size();i++){
			inID=NeighborCon[t][i].first;
			if(HigherIn.find(inID)!=HigherIn.end()){
				inW=HigherIn[inID];
				if(NeighborCon[t][i].second.first==wt+inW){
					NeighborCon[t][i].second.second-=1;
					if(NeighborCon[t][i].second.second<1){
						OrderComp oc={t,inID};
						OC.insert(oc);
						OCdis[make_pair(t,inID)]=wt+inW;
					}
				}
			}
		}
		for(int i=0;i<LowerIn.size();i++){
			inID=LowerIn[i].first; inW=LowerIn[i].second;
			for(int j=0;j<NeighborCon[inID].size();j++){
				if(NeighborCon[inID][j].first==t){
					if(NeighborCon[inID][j].second.first==inW+wt){
						NeighborCon[inID][j].second.second-=1;
						if(NeighborCon[inID][j].second.second<1){
							OrderComp oc={inID,t};
							OC.insert(oc);
							OCdis[make_pair(inID,t)]=wt+inW;
						}
					}
					break;
				}
			}
		}

		//get the new weight value of s-->t
		wt=INF; int countwt=0;
		for(int i=0;i<Neighbor[s].size();i++){
			if(Neighbor[s][i].first==t){
				wt=Neighbor[s][i].second;//the weight value in the original graph
				countwt=1;
				break;
			}
		}
		int ssw,wtt,wid;
		vector<int> Wnodes; //Wnodes.clear();
		if(s<t){
			//Wnodes=SCconNodes[make_pair(s,t)]; //cout<<"wid num "<<Wnodes.size()<<endl;
			Wnodes=SCconNodesMT[s][t];
		}else{
			//Wnodes=SCconNodes[make_pair(t,s)];
			Wnodes=SCconNodesMT[t][s];
		}

		for(int i=0;i<Wnodes.size();i++){
			wid=Wnodes[i];
			for(int j=0;j<NeighborCon[wid].size();j++){
				if(NeighborCon[wid][j].first==s){
					ssw=NeighborCon[wid][j].second.first;
				}
				if(NeighborCon[wid][j].first==t){
					wtt=NeighborCon[wid][j].second.first;
				}
			}

			if(ssw+wtt<wt){
				wt=ssw+wtt;
				countwt=1;
			}else if(ssw+wtt==wt){
				countwt+=1;
			}
		}

		//refresh the weight value of s--t in the index
		for(int i=0;i<NeighborCon[s].size();i++){
			if(NeighborCon[s][i].first==t){

				if(TotalBoundSet.find(s)!=TotalBoundSet.end() && TotalBoundSet.find(t)!=TotalBoundSet.end()){
					if(PreW.find(make_pair(s,t))==PreW.end()){
						PreW.insert(make_pair(make_pair(s,t),NeighborCon[s][i].second.first));
					}
					NowW[make_pair(s,t)]=wt;
				}

				NeighborCon[s][i].second.first=wt;
				NeighborCon[s][i].second.second=countwt;
				break;
			}
		}
	}

	for(auto it=PreW.begin();it!=PreW.end();it++){
		totalupdate.push_back(make_pair((*it).first,make_pair((*it).second,NowW[(*it).first])));
	}
	return totalupdate;
}
