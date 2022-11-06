/*
 * PostBoundary.cpp
 *
 *  Created on: 26 Oct 2022
 *      Author: zhangmengxuan
 */
#include "head.h"

//query processing
void Graph::CorrectnessCheckPost(){
	srand (time(NULL));
	int s, t, d1, d2, d3;
	for(int i=0;i<1000;i++){
		s=rand()%nodenum;
		t=rand()%nodenum;
		d1=Dijkstra(s,t,Neighbors);
		d2=HopQueryPost(s,t);
		if(d1!=d2){
			cout<<"InCorrect!"<<s<<" "<<t<<" "<<d1<<" "<<d2<<endl;
		}
	}
}

int Graph::HopQueryPost(int ID1, int ID2){
	int d=INF;

	if(TotalBoundSet.find(ID1)!=TotalBoundSet.end() && TotalBoundSet.find(ID2)!=TotalBoundSet.end()){
		d=HopQueryLocal(ID1,ID2, LabelOverlay);
		//cout<<"out-out"<<endl;
	}else if(TotalBoundSet.find(ID1)!=TotalBoundSet.end()){
		d=HopQueryInOut(ID2,ID1);
		//cout<<"in-out"<<endl;
	}else if(TotalBoundSet.find(ID2)!=TotalBoundSet.end()){
		d=HopQueryInOut(ID1,ID2);
		//cout<<"in-out"<<endl;
	}else{
		if(VtoParID[ID1]==VtoParID[ID2]){//the same partition
			d=HopQueryLocal(ID1,ID2,LabelPartiPost[VtoParID[ID1]]);
			//d=HopQueryInIn(ID1,ID2);
			//cout<<"in-in same"<<endl;
		}else{//different partitions
			d=HopQueryOutOut(ID1,ID2);
			//cout<<"in-in different"<<endl;
		}
	}

	return d;
}

void Graph::EffiCheckPost1(){
	int num, ID1, ID2;
	vector<pair<int,int>> ODpair1,ODpair2,ODpair3,ODpair4;

	vector<vector<int>> VertexCategory;
	vector<int> vec;
	VertexCategory.assign(pnum+1,vec);//put vertices within the same partition together

	for(int ID=0;ID<VtoParID.size();ID++){
		if(TotalBoundSet.find(ID)==TotalBoundSet.end())
			VertexCategory[VtoParID[ID]].push_back(ID);
	}
	for(auto it=TotalBoundSet.begin();it!=TotalBoundSet.end();it++)
		VertexCategory[pnum].push_back(*it);

	//for(int k=0;k<VertexCategory.size();k++)
		//cout<<k<<" "<<VertexCategory[k].size()<<endl;

//cout<<"/////////////////start query///////////////////"<<endl;

	srand (time(NULL));
	int s, t, k, l;
	//ODpair1: both boundary vertices
	for(int i=0;i<1000;i++){
		s=rand()%VertexCategory[pnum].size();
		t=rand()%VertexCategory[pnum].size();
		ODpair1.push_back(make_pair(VertexCategory[pnum][s],VertexCategory[pnum][t]));
	}
//cout<<"111111111111111111"<<endl;
	//ODpair2: vertices within the same partition
	for(int i=0;i<1000;i++){
		k=rand()%pnum;
		s=rand()%VertexCategory[k].size();
		t=rand()%VertexCategory[k].size();
		ODpair2.push_back(make_pair(VertexCategory[k][s],VertexCategory[k][t]));
	}
//cout<<"22222222222222222222222"<<endl;
	//ODpair3: vertices in different partitions
	for(int i=0;i<1000;i++){
		k=rand()%pnum;
		l=rand()%pnum;
		if(k!=l){
			s=rand()%VertexCategory[k].size();
			t=rand()%VertexCategory[l].size();
			ODpair3.push_back(make_pair(VertexCategory[k][s],VertexCategory[l][t]));
		}else
			i--;
	}
//cout<<"3333333333333333333333"<<endl;
	//ODpair4: one within partition, one boundary vertex
	for(int i=0;i<1000;i++){
		k=rand()%pnum;
		s=rand()%VertexCategory[k].size();
		t=rand()%VertexCategory[pnum].size();
		ODpair4.push_back(make_pair(VertexCategory[k][s],VertexCategory[pnum][t]));
	}
//cout<<"44444444444444444444444444"<<endl;
	std::chrono::high_resolution_clock::time_point t1, t2;
	std::chrono::duration<double> time_span;
	double runT;

	t1=std::chrono::high_resolution_clock::now();
	for(int i=0;i<ODpair1.size();i++){
		HopQueryLocal(ODpair1[i].first,ODpair1[i].second,LabelOverlay);
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"1: both OD are boundary: "<<1000*runT/ODpair1.size()<<" ms."<<endl;

	t1=std::chrono::high_resolution_clock::now();
	for(int i=0;i<ODpair2.size();i++){
		HopQueryLocal(ODpair2[i].first,ODpair2[i].second,LabelPartiPost[VtoParID[ODpair2[i].first]]);
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"2: OD within same partition: "<<1000*runT/ODpair2.size()<<" ms."<<endl;

	t1=std::chrono::high_resolution_clock::now();
	for(int i=0;i<ODpair3.size();i++){
		HopQueryOutOut(ODpair3[i].first,ODpair3[i].second);
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"3: OD in different partition: "<<1000*runT/ODpair3.size()<<" ms."<<endl;

	t1=std::chrono::high_resolution_clock::now();
	for(int i=0;i<ODpair4.size();i++){
		HopQueryInOut(ODpair4[i].first,ODpair4[i].second);
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"4: one boundary, one within partition: "<<1000*runT/ODpair4.size()<<" ms."<<endl;
}


void Graph::EffiCheckPost(string filename,int runtimes){
	ifstream IF(filename);
	if(!IF){
		cout<<"Cannot open Map "<<filename<<endl;
	}
	int num, ID1, ID2;
	vector<pair<int,int>> ODpair;
	IF>>num;
	for(int k=0;k<num;k++){
		IF>>ID1>>ID2;
		ODpair.push_back(make_pair(ID1, ID2));
	}
    if(runtimes > num){
        runtimes = num;
    }
    cout<<"Run times: "<<runtimes<<endl;
	int s, t;
	std::chrono::high_resolution_clock::time_point t1, t2;
	std::chrono::duration<double> time_span;
	double runT;

	t1=std::chrono::high_resolution_clock::now();
	for(int i=0;i<runtimes;i++){
		HopQueryPost(ODpair[i].first,ODpair[i].second);
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"PLL query time: "<<1000*runT/runtimes<<" ms."<<endl;
}

void Graph::IndexsizePost(){
	long long m=0;

	for(int pid=0;pid<pnum;pid++){
		for(int k=0;k<LabelParti[pid].size();k++){
			m+=LabelParti[pid][k].size()*2*sizeof(int);
		}

		for(int i=0;i<PruningPointParti[pid].size();i++){
			for(unordered_map<int,vector<int>>::iterator it=PruningPointParti[pid][i].begin();it!=PruningPointParti[pid][i].end();it++){
				m+=(1+(*it).second.size())*sizeof(int);
			}
		}
	}

	for(int pid=0;pid<pnum;pid++){
		for(int k=0;k<LabelPartiPost[pid].size();k++){
			m+=LabelPartiPost[pid][k].size()*2*sizeof(int);
		}

		for(int i=0;i<PruningPointPartiPost[pid].size();i++){
			for(unordered_map<int,vector<int>>::iterator it=PruningPointPartiPost[pid][i].begin();it!=PruningPointPartiPost[pid][i].end();it++){
				m+=(1+(*it).second.size())*sizeof(int);
			}
		}
	}

	for(int k=0;k<LabelOverlay.size();k++){
		m+=LabelOverlay[k].size()*2*sizeof(int);
	}

	for(int i=0;i<PruningPointOverlay.size();i++){
		for(unordered_map<int,vector<int>>::iterator it=PruningPointOverlay[i].begin();it!=PruningPointOverlay[i].end();it++){
			m+=(1+(*it).second.size())*sizeof(int);
		}
	}

	cout<<"Index size "<<(double)m/1024/1024<<" MB"<<endl;
}

//insert edge into subgraphs
void Graph::insertion(int a, int b, int newW, vector<vector<pair<int,int>>> &Neighbors, vector<unordered_map<int,int>> &Label){
	//(a,b) not exist beforehand
	Neighbors[a].push_back(make_pair(b,newW));
	Neighbors[b].push_back(make_pair(a,newW));


	//Chnum=0;
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
		//Chnum+=1;
		Change[LID].push_back(make_pair(HID,newW));
		WaitPro.insert(LID);

		//check the label of a,b
		int hubid, hubdis;
		unordered_map<int,int>::iterator it1=Label[LID].begin();
		for(;it1!=Label[LID].end();it1++){
			hubid=(*it1).first; hubdis=(*it1).second;
			if(NodeOrder[hubid]>NodeOrder[HID] && newW+hubdis<ShortestDisQuery(HID,hubid,Label)){
				Label[HID][hubid]=newW+hubdis;
				//Chnum+=1;
				Change[HID].push_back(make_pair(hubid, newW+hubdis));
				WaitPro.insert(HID);
			}
		}
		unordered_map<int,int>::iterator it2=Label[HID].begin();
		for(;it2!=Label[HID].end();it2++){
			hubid=(*it2).first; hubdis=(*it2).second;
			if(newW+hubdis<ShortestDisQuery(LID, hubid,Label)){
				Label[LID][hubid]=newW+hubdis;
				//Chnum+=1;
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
							//Chnum+=1;
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

void Graph::PartitionUpdate(){
	int ID1,ID2,wlocal,woverlay;
	for(int k=0;k<pnum;k++){
		//boundary edges
		for(int i=0;i<BoundVer[k].size();i++){
			for(int j=i+1;j<BoundVer[k].size();j++){
				ID1=BoundVer[k][i];
				ID2=BoundVer[k][j];
				wlocal=HopQueryLocal(ID1,ID2,LabelParti[k]);
				woverlay=HopQueryLocal(ID1,ID2,LabelOverlay);
				bool found=false;//whether the boundary edge exist or not
				int wei;
				if(woverlay<wlocal){
					for(int h=0;h<NeighborsPartiPost[k][ID1].size();h++){
						if(NeighborsPartiPost[k][ID1][h].first==ID2){
							found=true;
							wei=NeighborsPartiPost[k][ID1][h].second;
							break;
						}
					}
					if(found)
						DecreasePLL(ID1,ID2,wei,woverlay,NeighborsPartiPost[k],LabelPartiPost[k]);
					else
						insertion(ID1,ID2,woverlay,NeighborsPartiPost[k],LabelPartiPost[k]);
				}else if(woverlay>wlocal)
					cout<<"Something wrong: shortest path in the overlay graph rather than in the subgraph."<<endl;
			}
		}
	}
}

/*void Graph::DecreasePost(int a, int b, int oldW, int newW){



}


void Graph::IncreasePost(int a, int b, int oldW, int newW){



}*/
