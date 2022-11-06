/*
 * update.cpp
 *
 *  Created on: 20 Sep 2022
 *      Author: zhangmengxuan
 */
#include "head.h"

void Graph::DecreaseSingle(int a, int b, int newW){
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
		vector<pair<pair<int,int>,int>> weightOverlay;//collect the changed edges on overlay graph
		weightOverlay.clear();
		Decrease(a,b,newW,NeighborsParti[pid],Trees[pid],ranks[pid],heightMaxs[pid]);

		//weightOverlay collect the changed edges on overlay graph
		vector<int> Bid=BoundVer[pid];
		//check the boundary edge within partition
		int bid1,bid2,olddis,newdis;
		for(int i=0;i<Bid.size();i++){
			bid1=Bid[i];
			for(int j=i+1;j<Bid.size();j++){
				bid2=Bid[j];
				olddis=BedgePID[bid1][bid2].first;
				newdis=QueryH2HPartition(bid1,bid2,pid);
				if(newdis<olddis){
					BedgePID[bid1][bid2].first=newdis;
					BedgePID[bid2][bid1].first=newdis;
					weightOverlay.push_back(make_pair(make_pair(bid1,bid2),newdis));
				}
			}
		}

		DecreaseBatch(weightOverlay,NeighborsOverlay,TreeOverlay,rankOverlay,heightMaxOverlay);
		//cout<<"Overlay update number "<<weightOverlay.size()<<endl;
		//update the overlay graph index, after partition index update
		/*for(int l=0;l<weightOverlay.size();l++){
			Decrease(weightOverlay[l].first.first,weightOverlay[l].first.second,weightOverlay[l].second,NeighborsOverlay,TreeOverlay,rankOverlay,heightMaxOverlay);
		}*/
	}else{
		Decrease(a, b, newW,NeighborsOverlay,TreeOverlay,rankOverlay,heightMaxOverlay);
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
		//cout<<"zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz"<<endl;
		vector<pair<pair<int,int>,pair<int,int>>> weightOverlay;//collect the changed edges on overlay graph
		weightOverlay.clear();
		Increase(a,b,oldW,newW,NeighborsParti[pid],Trees[pid],ranks[pid],heightMaxs[pid],SCconNodesMTs[pid],VidtoTNids[pid]);
		//cout<<"/////////////////////////////////////////"<<endl;

		//cout<<"boundary edge checkkkkkkkkkkkkkkkkkkkkkkkkkkk"<<endl;
		//boundary edges check
		vector<int> Bid=BoundVer[pid];
		int bid1,bid2,olddis,newdis;
		for(int i=0;i<Bid.size();i++){
			bid1=Bid[i];
			for(int j=i+1;j<Bid.size();j++){
				bid2=Bid[j];
				olddis=BedgePID[bid1][bid2].first;
				newdis=QueryH2HPartition(bid1,bid2,pid);
				BedgePID[bid1][bid2].first=newdis;
				BedgePID[bid2][bid1].first=newdis;
				int overlaydis=HopQueryOverlay(bid1,bid2);
				if(newdis>olddis)//if '=', not problem; if '<', problem
					weightOverlay.push_back(make_pair(make_pair(bid1,bid2),make_pair(olddis,newdis)));
				else if(newdis<olddis)
					cout<<"Something wrong happens."<<endl;

			}
		}

		IncreaseBatch(weightOverlay,NeighborsOverlay,TreeOverlay,rankOverlay,heightMaxOverlay,SCconNodesOverlayMT,VidtoTNidOverlay);
		//update the overlay graph index, after partition index update
		/*cout<<"Overlay update number "<<weightOverlay.size()<<endl;
		for(int l=0;l<weightOverlay.size();l++){
			Increase(weightOverlay[l].first.first,weightOverlay[l].first.second,weightOverlay[l].second.first,weightOverlay[l].second.second,NeighborsOverlay,TreeOverlay,rankOverlay,heightMaxOverlay,SCconNodesOverlayMT,VidtoTNidOverlay);
		}*/
		//cout<<"''''''''''''''''''''''''''''''''''''''''''"<<endl;
	}else{
		Increase(a, b, oldW, newW,NeighborsOverlay,TreeOverlay,rankOverlay,heightMaxOverlay,SCconNodesOverlayMT,VidtoTNidOverlay);
		//cout<<",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,"<<endl;
	}
}
