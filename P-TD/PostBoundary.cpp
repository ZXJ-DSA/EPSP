/*
 * PostBoundary.cpp
 *
 *  Created on: 27 Oct 2022
 *      Author: zhangmengxuan
 */
#include "head.h"

void Graph::PartitionUpdate(){
	int ID1,ID2,wlocal,woverlay;
	for(int k=0;k<pnum;k++){
		//boundary edges
		for(int i=0;i<BoundVer[k].size();i++){
			for(int j=i+1;j<BoundVer[k].size();j++){
				ID1=BoundVer[k][i];
				ID2=BoundVer[k][j];
				wlocal=QueryH2HPartition(ID1,ID2,k);
				woverlay=HopQueryOverlay(ID1,ID2);
				bool found=false;//whether the boundary edge exist or not
				int wei;
				if(woverlay<wlocal){
					for(int h=0;h<NeighborsPartiPost[k][ID1].size();h++){
						if(NeighborsPartiPost[k][ID1][h].first==ID2){
							found=true;
							NeighborsPartiPost[k][ID1][h].second=woverlay;
							break;
						}
					}
					if(found){
						for(int h=0;h<NeighborsPartiPost[k][ID2].size();h++){
							if(NeighborsPartiPost[k][ID2][h].first==ID1){
								NeighborsPartiPost[k][ID2][h].second=woverlay;
								break;
							}
						}
					}else{
						NeighborsPartiPost[k][ID1].push_back(make_pair(ID2,woverlay));
						NeighborsPartiPost[k][ID2].push_back(make_pair(ID1,woverlay));
					}

				}else if(woverlay>wlocal)
					cout<<"Something wrong: shortest path in the overlay graph rather than in the subgraph."<<endl;
			}
		}
	}
}

int Graph::HopQueryPost(int ID1, int ID2){
	int d=INF;

	if(TotalBoundSet.find(ID1)!=TotalBoundSet.end() && TotalBoundSet.find(ID2)!=TotalBoundSet.end()){
		d=HopQueryOverlay(ID1,ID2);
	}else if(TotalBoundSet.find(ID1)!=TotalBoundSet.end()){
		d=HopQueryInOut(ID2,ID1);
	}else if(TotalBoundSet.find(ID2)!=TotalBoundSet.end()){
		d=HopQueryInOut(ID1,ID2);
	}else{
		if(VtoParID[ID1]==VtoParID[ID2]){//the same partition
			d=QueryH2HPartition(ID1,ID2,VtoParID[ID2]);
		}else{//different partitions
			d=HopQueryOutOut(ID1,ID2);
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
        while(s==t){
            t=rand()%VertexCategory[pnum].size();
        }
		ODpair1.push_back(make_pair(VertexCategory[pnum][s],VertexCategory[pnum][t]));
	}
//cout<<"111111111111111111"<<endl;
	//ODpair2: vertices within the same partition
	for(int i=0;i<1000;i++){
		k=rand()%pnum;
		s=rand()%VertexCategory[k].size();
		t=rand()%VertexCategory[k].size();
        while(s==t){
            t=rand()%VertexCategory[k].size();
        }
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
		}

	}
//cout<<"3333333333333333333333"<<endl;
	//ODpair4: one within partition, one boundary vertex
	for(int i=0;i<1000;i++){
		k=rand()%pnum;
		s=rand()%VertexCategory[k].size();
		t=rand()%VertexCategory[pnum].size();
        while(s==t){
            t=rand()%VertexCategory[pnum].size();
        }
		ODpair4.push_back(make_pair(VertexCategory[k][s],VertexCategory[pnum][t]));
	}
//cout<<"44444444444444444444444444"<<endl;
	std::chrono::high_resolution_clock::time_point t1, t2;
	std::chrono::duration<double> time_span;
	double runT;

	t1=std::chrono::high_resolution_clock::now();
	for(int i=0;i<ODpair1.size();i++){
		HopQueryOverlay(ODpair1[i].first,ODpair1[i].second);
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"1: both OD are boundary: "<<1000*runT/ODpair1.size()<<" ms."<<endl;

	t1=std::chrono::high_resolution_clock::now();
	for(int i=0;i<ODpair2.size();i++){
		QueryH2HPartition(ODpair2[i].first,ODpair2[i].second,VtoParID[ODpair2[i].first]);
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
	cout<<"average query time: "<<1000*runT/runtimes<<" ms."<<endl;
}

void Graph::IndexsizePost(){
	long long m=0;

	for(int pid=0;pid<pnum;pid++){
		for(int i=0;i<Trees[pid].size();i++){
			m+=Trees[pid][i].dis.size()*sizeof(int);
			m+=Trees[pid][i].vert.size()*3*sizeof(int);
		}
		for(int i=0;i< SCconNodesMTs[pid].size();i++){
			for(auto it=SCconNodesMTs[pid][i].begin(); it!=SCconNodesMTs[pid][i].end(); it++){
				m+=sizeof(int)+(*it).second.size()*sizeof(int);
			}
		}
		for(int i=0;i<VidtoTNids[pid].size();i++){
			m+=VidtoTNids[pid][i].size()*sizeof(int);
		}
	}
	m=m*2;

	for(int i=0;i<TreeOverlay.size();i++){
		m+=TreeOverlay[i].dis.size()*sizeof(int);
		m+=TreeOverlay[i].vert.size()*3*sizeof(int);
	}
	for(int i=0;i< SCconNodesOverlayMT.size();i++){
		for(auto it=SCconNodesOverlayMT[i].begin(); it!=SCconNodesOverlayMT[i].end(); it++){
			m+=sizeof(int)+(*it).second.size()*sizeof(int);
		}
	}
	for(int i=0;i<VidtoTNidOverlay.size();i++)
		m+=VidtoTNidOverlay[i].size()*sizeof(int);
	cout<<"Index size "<<(double)m/1024/1024<<" MB"<<endl;
}
