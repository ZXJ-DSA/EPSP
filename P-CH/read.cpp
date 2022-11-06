/*
 * read.cpp
 *
 *  Created on: 19 Oct 2022
 *      Author: zhangmengxuan
 */
#include "head.h"

void Graph::UpdateRead(string filename,vector<pair<pair<int,int>,int>>& TestData){
    TestData.clear();

    int num, ID1, ID2, neww;
    ifstream IF(filename);
    if(!IF){
        cout<<"Cannot open "<<filename<<endl;
        exit(1);
    }
    IF>>num;
    for(int i=0;i<num;i++){
        IF>>ID1>>ID2>>neww;
        TestData.push_back(make_pair(make_pair(ID1, ID2), neww));
    }
    IF.close();
}

void Graph::GraphRead(string filename){
	ifstream inGraph(filename);
	if(!inGraph){
		cout<<"Cannot open Map "<<filename<<endl;
        exit(1);
	}

	int num1,num2;
	int ID1,ID2,weight;
	inGraph>>num1>>num2;
	nodenum=num1;
	edgenum=num2;

	vector<pair<int,int>> vecp; vecp.clear();
	Neighbor.assign(nodenum, vecp);
	EdgeWei.clear();
	for(int i=0;i<num2;i++){
		inGraph>>ID1>>ID2>>weight;
		Neighbor[ID1].push_back(make_pair(ID2,weight));
		//Neighbors[ID2].push_back(make_pair(ID1,weight));//use this line when one edge is written once in the graph file
		if(ID1<ID2){
			EdgeWei[make_pair(ID1,ID2)]=weight;
		}else{
			EdgeWei[make_pair(ID2,ID1)]=weight;
		}
	}
	cout<<"Finish Graph Reading! Nnum "<<nodenum<<" Enum "<<edgenum<<endl;
}

void Graph::OrderRead(string filename){
	NodeOrder.assign(nodenum,-1);
	vNodeOrder.assign(nodenum,-1);

	ifstream IF(filename);
	if(!IF){
        cout<<"Cannot open file "<<filename<<endl;
        exit(1);
    }

	int num,ID,order;
	IF>>num;
	for(int i=0;i<num;i++){
		IF>>ID>>order;
		NodeOrder[ID]=order;
		vNodeOrder[order]=ID;
	}
    for(int i=0;i<3;++i){
        cout<<i<<" "<<vNodeOrder[i]<<" "<<Neighbor[vNodeOrder[i]].size()<<endl;
    }
}

void Graph::GraphPartitionRead(string filename){
	//read the partitioned graphs
	vector<pair<int,int>> vec; vec.clear();
	vector<vector<pair<int,int>>> vecvec;
	vecvec.assign(nodenum,vec);

	//initialize edges' PID
	unordered_map<int,int> unomap;
	unomap.clear();
	EtoParID.assign(nodenum,unomap);

	ifstream IF(filename+"/subgraph_edge");
	if(!IF){
        cout<<"Cannot open file "<<"subgraph_edge"<<endl;
        exit(1);
    }

	int pnum1;
	IF>>pnum1;
    if(partiName == "NC"){
        pnum = pnum1;
    }
    cout<<"Real Partition number: "<<pnum1<<endl;
    NeighborsParti.assign(pnum, vecvec);
	for(int k=0;k<pnum1;k++){
		int edgenum0,ID1,ID2,weight;
		IF>>edgenum0;
		for(int i=0;i<edgenum0;i++){
			IF>>ID1>>ID2>>weight;
			NeighborsParti[k][ID1].push_back(make_pair(ID2,weight));

			if(EtoParID[ID1].find(ID2)!=EtoParID[ID1].end())
				cout<<"something wrong: edge ("<<ID1<<", "<<ID2<<") alrealy in other partitions"<<endl;
			EtoParID[ID1].insert(make_pair(ID2,k));
		}
	}

	//read the vertex in each partition
	VtoParID.assign(nodenum,-1);
	ifstream IF1(filename+"/subgraph_vertex");
	if(!IF1){
        cout<<"Cannot open file "<<"subgraph_vertex"<<endl;
        exit(1);
    }

	int pnum2;
	IF1>>pnum2;
	for(int k=0;k<pnum2;k++){
		int vernum,ID;
		IF1>>vernum;
		for(int i=0;i<vernum;i++){
			IF1>>ID;

			/*if(ID>=nodenum)
				cout<<"ID "<<ID<<" ,partition ID "<<k<<endl;*/

			if(VtoParID[ID]==-1)
				VtoParID[ID]=k;
			else
				cout<<"vertex already in one partition!"<<ID<<" "<<VtoParID[ID]<<" "<<k<<endl;
		}
	}
	//further check that each vertex is in one and only one partition
	for(int vid=0;vid<nodenum;vid++){
		if(VtoParID[vid]==-1)
			cout<<"vertex "<<vid<<" not within any partition"<<endl;
	}

	//read the cut edges
	CutEdges.clear();
	set<int> ss; ss.clear();
	BoundVerSet.assign(pnum,ss);
	ifstream IF2(filename+"/cut_edges");
	if(!IF2){
        cout<<"Cannot open file "<<"cut_edges"<<endl;
        exit(1);
    }

	int ednum,ID1,ID2,weight;
	IF2>>ednum;
	for(int i=0;i<ednum;i++){
		IF2>>ID1>>ID2>>weight;
		CutEdges.push_back(make_pair(make_pair(ID1,ID2),weight));

		if(VtoParID[ID1]==VtoParID[ID2])
			cout<<"two end points of cut edge are in the same partition"<<endl;

		BoundVerSet[VtoParID[ID1]].insert(ID1);

		if(EtoParID[ID1].find(ID2)!=EtoParID[ID1].end())
			cout<<"something wrong: cut edge ("<<ID1<<", "<<ID2<<") alrealy in partitions"<<endl;

		EtoParID[ID1].insert(make_pair(ID2,pnum));
	}

	//further check the edges
	int calEdgeNum=0;
	for(int i=0;i<EtoParID.size();i++){
		calEdgeNum+=EtoParID[i].size();
	}
	cout<<"calculated edge number "<<calEdgeNum<<", graph edge number "<<edgenum<<endl;

	vector<int> vv; vv.clear();
	BoundVer.assign(pnum,vv);
	TotalBoundSet.clear();
	for(int k=0;k<pnum;k++){
		for(set<int>::iterator it=BoundVerSet[k].begin();it!=BoundVerSet[k].end();it++){
			BoundVer[k].push_back(*it);
			TotalBoundSet.insert(*it);
		}
	}

	cout<<"Partition data finish reading!"<<endl;
}

void Graph::EffiCheck(string filename, int runtimes){
	ifstream IF(filename);
	if(!IF){
		cout<<"Cannot open Map "<<filename<<endl;
	}
	int num, ID1, ID2;
	vector<pair<int,int>> ODpair;
	IF>>num;
    if(runtimes>num){
        runtimes = num;
    }
    cout<<"Run times: "<<runtimes<<endl;
	for(int k=0;k<num;k++){
		IF>>ID1>>ID2;
		ODpair.push_back(make_pair(ID1, ID2));
	}

	int s, t;
	std::chrono::high_resolution_clock::time_point t1, t2;
	std::chrono::duration<double> time_span;
	double runT;

	t1=std::chrono::high_resolution_clock::now();
	for(int i=0;i<runtimes;i++){//ODpair.size()
		Query(ODpair[i].first,ODpair[i].second);
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"average query time: "<<1000*runT/runtimes<<" ms."<<endl;//ODpair.size()
}

void Graph::EffiCheck(){
    int num, ID1, ID2;
    vector<pair<int,int>> ODpair1;
    vector<pair<int,int>> ODpair2;
    vector<pair<int,int>> ODpair3;
    vector<pair<int,int>> ODpair4;

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

    //ODpair3: vertices in different partitions
    for(int i=0;i<1000;i++){
        k=rand()%pnum;
        l=rand()%pnum;
        if(k!=l){
            s=rand()%VertexCategory[k].size();
            t=rand()%VertexCategory[l].size();
            ODpair3.push_back(make_pair(VertexCategory[k][s],VertexCategory[l][t]));
        }else{
            i--;
        }
    }

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

    std::chrono::high_resolution_clock::time_point t1, t2;
    std::chrono::duration<double> time_span;
    double runT;

    runT=0;
    double dijkT=0;
    int d1,d2;
    for(int i=0;i<ODpair1.size();i++){
        t1=std::chrono::high_resolution_clock::now();
        d1=QueryCH(ODpair1[i].first,ODpair1[i].second,NeighborConOverlay);
        t2=std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
        runT+= time_span.count();
    }
    cout<<"both OD are boundary: "<<1000*runT/ODpair1.size()<<" ms."<<endl;

    t1=std::chrono::high_resolution_clock::now();
    for(int i=0;i<ODpair2.size();i++){
        QueryCHInIn(ODpair2[i].first,ODpair2[i].second,VtoParID[ODpair2[i].first]);
    }
    t2=std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
    runT= time_span.count();
    cout<<"OD within same partition: "<<1000*runT/ODpair2.size()<<" ms."<<endl;

    t1=std::chrono::high_resolution_clock::now();
    for(int i=0;i<ODpair3.size();i++){
        QueryCHOutOut(ODpair3[i].first,ODpair3[i].second,VtoParID[ODpair3[i].first],VtoParID[ODpair3[i].second]);
    }
    t2=std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
    runT= time_span.count();
    cout<<"OD in different partition: "<<1000*runT/ODpair3.size()<<" ms."<<endl;

    runT=0; dijkT=0;
    for(int i=0;i<ODpair4.size();i++){
        t1=std::chrono::high_resolution_clock::now();
        d1=QueryCHInOut(ODpair4[i].first,ODpair4[i].second,VtoParID[ODpair4[i].first]);
        t2=std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
        runT += time_span.count();
    }

    cout<<"one boundary, one within partition: "<<1000*runT/ODpair4.size()<<" ms."<<endl;
}

//overlay graph construction
/*void Graph::OverlayGraphConstructPost(){
	vector<unordered_map<int,int>> OverlayGraph1;
	//initialize map of Neighbors
	//used to get minimum value
	unordered_map<int,int> map0; map0.clear();
	OverlayGraph1.assign(nodenum, map0);

	//initialize Neighbors
	vector<pair<int,int>> vecp; vecp.clear();
	NeighborsOverlay.assign(nodenum, vecp);

	//initialize Boundary edges' supported PID
	unordered_map<int,int> mp;
	mp.clear();
	BedgePID.assign(nodenum,mp);

	int ID1,ID2,weight;
	for(int k=0;k<pnum;k++){
		//boundary edges
		for(int i=0;i<BoundVer[k].size();i++){
			for(int j=i+1;j<BoundVer[k].size();j++){
				ID1=BoundVer[k][i];
				ID2=BoundVer[k][j];

				weight=QueryCH(ID1,ID2,NeighborCons[k]);

				BedgePID[ID1][ID2]=weight;
				BedgePID[ID2][ID1]=weight;
				OverlayGraph1[ID1].insert(make_pair(ID2,weight));
				OverlayGraph1[ID2].insert(make_pair(ID1,weight));
			}
		}
	}

	for(int h=0;h<CutEdges.size();h++){
		ID1=CutEdges[h].first.first;
		ID2=CutEdges[h].first.second;
		weight=CutEdges[h].second;

		if(OverlayGraph1[ID1].find(ID2)==OverlayGraph1[ID1].end()){
			OverlayGraph1[ID1].insert(make_pair(ID2,weight));
		}else{
			cout<<"abnormal since cut edge is also boundary edge"<<endl;
		}
	}

	for(int i=0;i<nodenum;i++){
		for(unordered_map<int,int>::iterator it=OverlayGraph1[i].begin();it!=OverlayGraph1[i].end();it++){
			NeighborsOverlay[i].push_back(make_pair((*it).first, (*it).second));
		}
	}
}*/

//only use the shortcut between boundary vertices to construct the overlay graph
void Graph::OverlayGraphConstructPost(){
	//initialize Neighbors
	vector<pair<int,int>> vecp; vecp.clear();
	NeighborsOverlay.assign(nodenum, vecp);

	int ID1,ID2,weight;
	for(int k=0;k<pnum;k++){
		//boundary edges
		for(int i=0;i<BoundVer[k].size();i++){
			ID1=BoundVer[k][i];
			for(int j=0;j<NeighborCons[k][ID1].size();j++){
				ID2=NeighborCons[k][ID1][j].first;
				weight=NeighborCons[k][ID1][j].second.first;
				if(BoundVerSet[k].find(ID2)!=BoundVerSet[k].end()){
					NeighborsOverlay[ID1].push_back(make_pair(ID2,weight));
					NeighborsOverlay[ID2].push_back(make_pair(ID1,weight));
				}
			}
			/*for(int j=i+1;j<BoundVer[k].size();j++){
				ID1=BoundVer[k][i];
				ID2=BoundVer[k][j];

				weight=QueryCH(ID1,ID2,NeighborCons[k]);

				BedgePID[ID1][ID2]=weight;
				BedgePID[ID2][ID1]=weight;
				OverlayGraph1[ID1].insert(make_pair(ID2,weight));
				OverlayGraph1[ID2].insert(make_pair(ID1,weight));
			}*/
		}
	}

	for(int h=0;h<CutEdges.size();h++){
		ID1=CutEdges[h].first.first;
		ID2=CutEdges[h].first.second;
		weight=CutEdges[h].second;

		NeighborsOverlay[ID1].push_back(make_pair(ID2,weight));
	}
}

void Graph::CorrectnessCheckOverlay(){
	srand (time(NULL));
	int s, t, d1, d2, d3;

	set<int> Bvertex;
	Bvertex.clear();
	for(int pid=0;pid<BoundVer.size();pid++){
		for(int i=0;i<BoundVer[pid].size();i++)
			Bvertex.insert(BoundVer[pid][i]);
	}
	vector<int> Bver;
	Bver.clear();
	for(set<int>::iterator it=Bvertex.begin();it!=Bvertex.end();it++){
		Bver.push_back(*it);
	}
	//cout<<"Boundary vertices size "<<Bver.size()<<endl;
	for(int i=0;i<1000;i++){
		s=rand()%Bver.size();
		t=rand()%Bver.size();
		d1=Dijkstra(Bver[s],Bver[t]);
		d2=DijkstraOverlay(Bver[s],Bver[t]);

		if(d1!=d2){
			cout<<"InCorrect!"<<Bver[s]<<" "<<Bver[t]<<" "<<d1<<" "<<d2<<endl;
		}
	}
	//cout<<"Correctness Check of Hierarchical graph is finished!"<<endl;
}

int Graph::Dijkstra(int ID1, int ID2){
	benchmark::heap<2, int, int> pqueue(nodenum);
	pqueue.update(ID1,0);

	vector<bool> closed(nodenum, false);
	vector<int> distance(nodenum, INF);
//	vector<int> prece(nodenum, 0);
	distance[ID1]=0;
	int topNodeID, topNodeDis;
	vector<pair<int,int>>::iterator it;
	int NNodeID,NWeigh;

	int d=INF;//initialize d to infinite for the unreachable case

	while(!pqueue.empty()){
		pqueue.extract_min(topNodeID, topNodeDis);
		if(topNodeID==ID2){
			d=distance[ID2];
			break;
		}
		closed[topNodeID]=true;

		for(it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
			NNodeID=(*it).first;
			NWeigh=(*it).second+topNodeDis;
			if(!closed[NNodeID]){
				if(distance[NNodeID]>NWeigh){
					distance[NNodeID]=NWeigh;
					pqueue.update(NNodeID, NWeigh);
				//	prece[NNodeID]=topNodeID;
				}
			}
		}
	}

	return d;
}

int Graph::DijkstraOverlay(int ID1, int ID2){
	benchmark::heap<2, int, int> pqueue(nodenum);
	pqueue.update(ID1,0);

	vector<bool> closed(nodenum, false);
	vector<int> distance(nodenum, INF);
//	vector<int> prece(nodenum, 0);
	distance[ID1]=0;
	int topNodeID, topNodeDis;
	vector<pair<int,int>>::iterator it;
	int NNodeID,NWeigh;

	int d=INF;//initialize d to infinite for the unreachable case

	while(!pqueue.empty()){
		pqueue.extract_min(topNodeID, topNodeDis);
		if(topNodeID==ID2){
			d=distance[ID2];
			break;
		}
		closed[topNodeID]=true;

		for(it=NeighborsOverlay[topNodeID].begin();it!=NeighborsOverlay[topNodeID].end();it++){
			NNodeID=(*it).first;
			NWeigh=(*it).second+topNodeDis;
			if(!closed[NNodeID]){
				if(distance[NNodeID]>NWeigh){
					distance[NNodeID]=NWeigh;
					pqueue.update(NNodeID, NWeigh);
				//	prece[NNodeID]=topNodeID;
				}
			}
		}
	}

	return d;
}

void Graph::CorrectnessCheck(){
	srand (time(NULL));
	int s, t, d1, d2, d3;
	for(int i=0;i<1000;i++){
		//cout<<i<<endl;
		s=rand()%nodenum;
		t=rand()%nodenum;
		d1=Dijkstra(s,t);
		d2=Query(s,t);
		if(d1!=d2){
			cout<<"InCorrect!"<<s<<" "<<t<<" "<<d1<<" "<<d2<<endl;
		}
	}
}

int Graph::Query(int ID1, int ID2){
	int d=INF;

	if(TotalBoundSet.find(ID1)!=TotalBoundSet.end() && TotalBoundSet.find(ID2)!=TotalBoundSet.end()){
		d=QueryCH(ID1,ID2,NeighborConOverlay);
	}else if(TotalBoundSet.find(ID1)!=TotalBoundSet.end()){
		d=QueryCHInOut(ID2,ID1,VtoParID[ID2]);
	}else if(TotalBoundSet.find(ID2)!=TotalBoundSet.end()){
		d=QueryCHInOut(ID1,ID2,VtoParID[ID1]);
	}else{
		if(VtoParID[ID1]==VtoParID[ID2]){//the same partition
			d=QueryCHInIn(ID1,ID2,VtoParID[ID1]);
		}else{//different partitions
			d=QueryCHOutOut(ID1,ID2,VtoParID[ID1],VtoParID[ID2]);
		}
	}

	return d;
}

void Graph::Indexsize(){
	long long m=0;

	for(int pid=0;pid<pnum;pid++){
		for(int i=0;i<NeighborCons[pid].size();i++){
			m+=NeighborCons[pid][i].size()*3*sizeof(int);
		}
		for(int i=0;i< SCconNodesMTs[pid].size();i++){
			for(auto it=SCconNodesMTs[pid][i].begin(); it!=SCconNodesMTs[pid][i].end(); it++){
				m+=sizeof(int)+(*it).second.size()*sizeof(int);
			}
		}
	}

	for(int i=0;i<NeighborConOverlay.size();i++){
		m+=NeighborConOverlay[i].size()*3*sizeof(int);
	}
	for(int i=0;i< SCconNodesMTOverlay.size();i++){
		for(auto it=SCconNodesMTOverlay[i].begin(); it!=SCconNodesMTOverlay[i].end(); it++){
			m+=sizeof(int)+(*it).second.size()*sizeof(int);
		}
	}
	cout<<"Index size "<<(double)m/1024/1024<<" MB"<<endl;
}
