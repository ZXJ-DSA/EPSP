/*
 * update.cpp
 *
 *  Created on: 19 Oct 2022
 *      Author: zhangmengxuan
 */
#include "head.h"

void Graph::DecreaseSingle(int a, int b, int oldW, int newW){
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

	int pid=EtoParID[a][b];
	if(pid<pnum){
		vector<pair<pair<int,int>,pair<int,int>>> weightOverlay;//collect the changed edges on overlay graph
		weightOverlay.clear();
		weightOverlay=CHdecNew(a,b,oldW,newW,NeighborsParti[pid],NeighborCons[pid]);
		/*CHdec(a,b,oldW,newW,NeighborsParti[pid],NeighborCons[pid]);

		//weightOverlay collect the changed edges on overlay graph
		vector<int> Bid=BoundVer[pid];
		//check the boundary edge within partition
		int bid1,bid2,olddis,newdis;
		for(int i=0;i<Bid.size();i++){
			bid1=Bid[i];
			for(int j=i+1;j<Bid.size();j++){
				bid2=Bid[j];
				olddis=BedgePID[bid1][bid2];
				newdis=QueryCH(bid1,bid2,NeighborCons[pid]);
				if(newdis<olddis){
					BedgePID[bid1][bid2]=newdis;
					BedgePID[bid2][bid1]=newdis;
					weightOverlay.push_back(make_pair(make_pair(bid1,bid2),make_pair(olddis,newdis)));
				}
			}
		}*/

		CHdecBat(weightOverlay,NeighborsOverlay,NeighborConOverlay);
	}else{
		CHdec(a,b,oldW,newW,NeighborsOverlay,NeighborConOverlay);
	}
}

void Graph::IncreaseSingle(int a, int b, int oldW, int newW){
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

	int pid=EtoParID[a][b];
	//cout<<"increase edge in partition "<<pid<<endl;
	if(pid<pnum){
		//cout<<"zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz"<<endl;
		vector<pair<pair<int,int>,pair<int,int>>> weightOverlay;//collect the changed edges on overlay graph
		weightOverlay.clear();
		weightOverlay=CHincNew(a,b,oldW,newW,NeighborsParti[pid],NeighborCons[pid],SCconNodesMTs[pid]);
		/*CHinc(a,b,oldW,newW,NeighborsParti[pid],NeighborCons[pid],SCconNodesMTs[pid]);

		//boundary edges check
		vector<int> Bid=BoundVer[pid];
		int bid1,bid2,olddis,newdis;
		for(int i=0;i<Bid.size();i++){
			bid1=Bid[i];
			for(int j=i+1;j<Bid.size();j++){
				bid2=Bid[j];
				olddis=BedgePID[bid1][bid2];
				newdis=QueryCH(bid1,bid2,NeighborCons[pid]);
				if(newdis>olddis){
					BedgePID[bid1][bid2]=newdis;
					BedgePID[bid2][bid1]=newdis;
					weightOverlay.push_back(make_pair(make_pair(bid1,bid2),make_pair(olddis,newdis)));
				}
			}
		}*/

		CHincBat(weightOverlay,NeighborsOverlay,NeighborConOverlay,SCconNodesMTOverlay);

	}else{
		CHinc(a, b, oldW, newW,NeighborsOverlay,NeighborConOverlay,SCconNodesMTOverlay);
	}
}



