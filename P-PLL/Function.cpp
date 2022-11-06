/*
 * Function.cpp
 *
 *  Created on: 6 Jul 2022
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
	}

	int num1,num2;
	int ID1,ID2,weight;
	inGraph>>num1>>num2;
	nodenum=num1;
	edgenum=num2;

	vector<pair<int,int>> vecp; vecp.clear();
	Neighbors.assign(nodenum, vecp);
	EdgeWei.clear();
	for(int i=0;i<num2;i++){
		inGraph>>ID1>>ID2>>weight;
		Neighbors[ID1].push_back(make_pair(ID2,weight));
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
	if(!IF)
		cout<<"Cannot open file "<<filename<<endl;
	int num,ID,order;
	IF>>num;
	for(int i=0;i<num;i++){
		IF>>ID>>order;
		NodeOrder[ID]=order;
		vNodeOrder[order]=ID;
	}
    for(int i=0;i<3;++i){
        cout<<i<<" "<<vNodeOrder[i]<<" "<<Neighbors[vNodeOrder[i]].size()<<endl;
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
	if(!IF)
		cout<<"Cannot open file "<<"subgraph_edge"<<endl;
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
	if(!IF1)
		cout<<"Cannot open file "<<"subgraph_vertex"<<endl;
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
	if(!IF2)
		cout<<"Cannot open file "<<"cut_edges"<<endl;
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
        exit(1);
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
		HopQuery(ODpair[i].first,ODpair[i].second);
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"PLL query time: "<<1000*runT/runtimes<<" ms."<<endl;//ODpair.size()
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
        ODpair1.push_back(make_pair(VertexCategory[pnum][s],VertexCategory[pnum][t]));
    }

    //ODpair2: vertices within the same partition
    for(int i=0;i<1000;i++){
        k=rand()%pnum;
        s=rand()%VertexCategory[k].size();
        t=rand()%VertexCategory[k].size();
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
        }else
            i--;
    }

    //ODpair4: one within partition, one boundary vertex
    for(int i=0;i<1000;i++){
        k=rand()%pnum;
        s=rand()%VertexCategory[k].size();
        t=rand()%VertexCategory[pnum].size();
        ODpair4.push_back(make_pair(VertexCategory[k][s],VertexCategory[pnum][t]));
    }

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
    cout<<"both OD are boundary: "<<1000*runT/ODpair1.size()<<" ms."<<endl;

    t1=std::chrono::high_resolution_clock::now();
    for(int i=0;i<ODpair2.size();i++){
        HopQueryInIn(ODpair2[i].first,ODpair2[i].second);
    }
    t2=std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
    runT= time_span.count();
    cout<<"OD within same partition: "<<1000*runT/ODpair2.size()<<" ms."<<endl;

    t1=std::chrono::high_resolution_clock::now();
    for(int i=0;i<ODpair3.size();i++){
        HopQueryOutOut(ODpair3[i].first,ODpair3[i].second);
    }
    t2=std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
    runT= time_span.count();
    cout<<"OD in different partition: "<<1000*runT/ODpair3.size()<<" ms."<<endl;

    t1=std::chrono::high_resolution_clock::now();
    for(int i=0;i<ODpair4.size();i++){
        HopQueryInOut(ODpair4[i].first,ODpair4[i].second);
    }
    t2=std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
    runT= time_span.count();
    cout<<"one boundary, one within partition: "<<1000*runT/ODpair4.size()<<" ms."<<endl;
}

void Graph::Indexsize(){
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
	map<int,vector<int>> mapvec;
	mapvec.clear();
	BedgeAllPID.assign(nodenum, mapvec);
	map<int,pair<int,set<int>>> mp;
	mp.clear();
	BedgePID.assign(nodenum,mp);

	int ID1,ID2,weight;
	for(int k=0;k<pnum;k++){
		//boundary edges
		for(int i=0;i<BoundVer[k].size();i++){
			for(int j=i+1;j<BoundVer[k].size();j++){
				ID1=BoundVer[k][i];
				ID2=BoundVer[k][j];
				weight=HopQueryLocal(ID1,ID2,LabelParti[k]);
				BedgePID[ID1][ID2].first=weight;
				BedgePID[ID2][ID1].first=weight;
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
	cout<<"Overlay Graph Finish Construction!"<<endl;
}*/

void Graph::OverlayGraphConstructPost(){
	//initialize Neighbors
	vector<pair<int,int>> vecp;
	vecp.clear();
	NeighborsOverlay.assign(nodenum, vecp);

	//initialize Boundary edges' supported PID
	map<int,vector<int>> mapvec;
	mapvec.clear();
	BedgeAllPID.assign(nodenum, mapvec);
	map<int,pair<int,set<int>>> mp;
	mp.clear();
	BedgePID.assign(nodenum,mp);

	int ID1,ID2,weight;
	for(int k=0;k<pnum;k++){
		//boundary edges
		for(int i=0;i<BoundVer[k].size();i++){
			for(int j=i+1;j<BoundVer[k].size();j++){
				ID1=BoundVer[k][i];
				ID2=BoundVer[k][j];
				weight=HopQueryLocal(ID1,ID2,LabelParti[k]);
				BedgePID[ID1][ID2].first=weight;
				BedgePID[ID2][ID1].first=weight;
				NeighborsOverlay[ID1].push_back(make_pair(ID2,weight));
				NeighborsOverlay[ID2].push_back(make_pair(ID1,weight));
			}
		}
	}

	for(int h=0;h<CutEdges.size();h++){
		ID1=CutEdges[h].first.first;
		ID2=CutEdges[h].first.second;
		weight=CutEdges[h].second;

		NeighborsOverlay[ID1].push_back(make_pair(ID2,weight));
	}
	cout<<"Overlay Graph Finish Construction!"<<endl;

	vector<pair<int,pair<int,int>>> Gedge;
	Gedge.clear();
	for(int i=0;i<nodenum;i++){
		for(int k=0;k<NeighborsOverlay[i].size();k++)
			Gedge.push_back(make_pair(i,NeighborsOverlay[i][k]));
	}

	/*ofstream OF("/home/s4451682/MengxuanZhang/partiPLL/OverlayGraph");
	OF<<nodenum<<" "<<Gedge.size()<<endl;
	for(int j=0;j<Gedge.size();j++)
		OF<<Gedge[j].first<<" "<<Gedge[j].second.first<<" "<<Gedge[j].second.second<<endl;*/
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
		d1=Dijkstra(Bver[s],Bver[t],Neighbors);
		d2=Dijkstra(Bver[s],Bver[t],NeighborsOverlay);
		//cout<<i<<endl;
		if(d1!=d2){
			cout<<"InCorrect!"<<Bver[s]<<" "<<Bver[t]<<" "<<d1<<" "<<d2<<endl;
		}
	}
	//cout<<"Correctness Check of Hierarchical graph is finished!"<<endl;
}

int Graph::Dijkstra(int ID1, int ID2, vector<vector<pair<int,int>>> &Neighbors){
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

		for(it=Neighbors[topNodeID].begin();it!=Neighbors[topNodeID].end();it++){
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
		s=rand()%nodenum;
		t=rand()%nodenum;
		d1=Dijkstra(s,t,Neighbors);
		d2=HopQuery(s,t);
		if(d1!=d2){
			cout<<"InCorrect!"<<s<<" "<<t<<" "<<d1<<" "<<d2<<endl;
		}
	}
}

int Graph::HopQuery(int ID1, int ID2){
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
			d=HopQueryInIn(ID1,ID2);
			//cout<<"in-in same"<<endl;
		}else{//different partitions
			d=HopQueryOutOut(ID1,ID2);
			//cout<<"in-in different"<<endl;
		}
	}

	return d;
}

int Graph::HopQueryLocal(int ID1, int ID2, vector<unordered_map<int,int>> &Label){
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

int Graph::HopQueryInOut(int in, int out){
	int d=INF;

	int boundid, tempdis;
	for(int i=0;i<BoundVer[VtoParID[in]].size();i++){
		boundid=BoundVer[VtoParID[in]][i];
		tempdis=HopQueryLocal(in, boundid, LabelParti[VtoParID[in]])+HopQueryLocal(boundid,out, LabelOverlay);
		if(tempdis<d)
			d=tempdis;
	}

	return d;
}

int Graph::HopQueryOutOut(int ID1, int ID2){
	int d=INF;

	vector<int> B1=BoundVer[VtoParID[ID1]];
	vector<int> B2=BoundVer[VtoParID[ID2]];

	map<int,int> m1,m2;
	m1.clear();
	m2.clear();
	int bID1, bID2, tempdis;
	for(int i=0;i<B1.size();i++){
		bID1=B1[i];
		m1.insert(make_pair(bID1,HopQueryLocal(ID1,bID1,LabelParti[VtoParID[ID1]])));
	}
	for(int j=0;j<B2.size();j++){
		bID2=B2[j];
		m2.insert(make_pair(bID2,HopQueryLocal(ID2,bID2,LabelParti[VtoParID[ID2]])));
	}

	for(int k=0;k<B1.size();k++){
		bID1=B1[k];

		if(m1[bID1]>d)
			continue;

		for(int z=0;z<B2.size();z++){
			bID2=B2[z];

			if(m2[bID2]>d)
				continue;

			tempdis=m1[bID1]+HopQueryLocal(bID1,bID2,LabelOverlay)+m2[bID2];
			if(tempdis<d)
				d=tempdis;
		}
	}

	return d;
}

int Graph::HopQueryInIn(int ID1, int ID2){
	int d=INF;

	int pid=VtoParID[ID1];
	int dis0=HopQueryLocal(ID1,ID2,LabelParti[pid]);
	if(dis0<d)
		d=dis0;

	vector<int> B=BoundVer[pid];
	map<int,int> m1,m2;
	m1.clear();
	m2.clear();
	vector<int> B1,B2;
	B1.clear();
	B2.clear();
	int bID,d1,d2;
	for(int i=0;i<B.size();i++){
		bID=B[i];
		d1=HopQueryLocal(ID1,bID,LabelParti[pid]);
		d2=HopQueryLocal(ID2,bID,LabelParti[pid]);
		if(d1<d){
			B1.push_back(bID);
			m1.insert(make_pair(bID,d1));
		}
		if(d2<d){
			B2.push_back(bID);
			m2.insert(make_pair(bID,d2));
		}
	}

	int bID1, bID2, tempdis;
	for(int k=0;k<B1.size();k++){
		bID1=B1[k];
		if(m1[bID1]>d)
			continue;
		for(int z=0;z<B2.size();z++){
			bID2=B2[z];
			if(m2[bID2]>d)
				continue;
			tempdis=m1[bID1]+HopQueryLocal(bID1,bID2,LabelOverlay)+m2[bID2];
			if(tempdis<d)
				d=tempdis;
		}
	}

	return d;
}

//PLL index construction
void Graph::IndexConst(){
	unordered_map<int,int> m;
	m.clear();
	Label.assign(nodenum,m);
	unordered_map<int,vector<int>> map1;
	map1.clear();
	PruningPointNew.assign(nodenum,map1);

	int ID;
	int cnt=0;
	for(int i=nodenum-1;i>=0;i--){
		ID=vNodeOrder[i];
		vector<pair<int,int>> vp;
		DijksPrune(ID,vp);
		cnt+=1;
		for(int j=0;j<vp.size();j++){
			Label[vp[j].first].insert(make_pair(ID, vp[j].second));
		}
	}
}
//function of DijksPrun for parallel implementation
/*void Graph::DijksPruneP(int nodeID, vector<pair<int,int>>& vp, vector<unordered_map<int,vector<int>>> & PruningPointNew){
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
        ShortestDisQuery(nodeID, topNodeID,SupNode,TempDis);
        if(TempDis<=topNodeDis){

            if(NodeOrder[topNodeID]<NodeOrder[nodeID]){
                for(int k=0;k<SupNode.size();k++){
                    int supn=SupNode[k];
                    PruningPointNew[topNodeID][supn].push_back(nodeID);
                    PruningPointNew[nodeID][supn].push_back(topNodeID);
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
}*/

void Graph::DijksPrune(int nodeID, vector<pair<int,int>>& vp){
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
		ShortestDisQuery(nodeID, topNodeID,SupNode,TempDis);
		if(TempDis<=topNodeDis){

			if(NodeOrder[topNodeID]<NodeOrder[nodeID]){
				for(int k=0;k<SupNode.size();k++){
					int supn=SupNode[k];
					PruningPointNew[topNodeID][supn].push_back(nodeID);
					PruningPointNew[nodeID][supn].push_back(topNodeID);
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

int Graph::ShortestDisQuery(int ID1,int ID2,vector<int>& SupNode, int& d){
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

//PLL index construction in multiple threads
/*void Graph::IndexConst1(){
	unordered_map<int,int> m; m.clear();
	vector<unordered_map<int,int>> vecmap;
	vecmap.assign(nodenum,m);
	LabelParti.assign(pnum,vecmap);
	unordered_map<int,vector<int>> map1; map1.clear();
	vector<unordered_map<int,vector<int>>> vecmapvec;
	vecmapvec.assign(nodenum,map1);
	PruningPointParti.assign(pnum,vecmapvec);

	boost::thread_group thread;
	int step=pnum/threadnum;
	for(int i=0;i<threadnum;i++){
		pair<int,int> p;
		p.first=i*step;
		if(i==threadnum-1)
			p.second=pnum;
		else
			p.second=(i+1)*step;
		cout<<"p's first and second value: "<<p.first<<" "<<p.second<<endl;
		thread.add_thread(new boost::thread(&Graph::IndexConstMT, this, p));
	}
	thread.join_all();
}

void Graph::IndexConstMT(pair<int,int> p){
	for(int k=p.first;k<p.second;k++){
		for(int i=nodenum-1;i>=0;i--){
			int ID=vNodeOrder[i];
			vector<pair<int,int>> vp;
			DijksPruneMT(ID,vp,k);
			for(int j=0;j<vp.size();j++){
				LabelParti[k][vp[j].first].insert(make_pair(ID, vp[j].second));
			}
		}
	}
}

void Graph::DijksPruneMT(int nodeID, vector<pair<int,int>>& vp,int partiID){
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
		ShortestDisQueryMT(nodeID, topNodeID,SupNode,TempDis,partiID);
		if(TempDis<=topNodeDis){

			if(NodeOrder[topNodeID]<NodeOrder[nodeID]){
				for(int k=0;k<SupNode.size();k++){
					int supn=SupNode[k];
					PruningPointParti[partiID][topNodeID][supn].push_back(nodeID);
					PruningPointParti[partiID][nodeID][supn].push_back(topNodeID);
					//cout<<"pruning point "<<topNodeID<<" "<<nodeID<<" "<<supn<<endl;
				}
			}
			continue;
		}


		//Label[topNodeID].insert(nodeID, topNodeDis);
		vp.push_back(make_pair(topNodeID,topNodeDis));
		for(it=NeighborsParti[partiID][topNodeID].begin();it!=NeighborsParti[partiID][topNodeID].end();it++){
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

int Graph::ShortestDisQueryMT(int ID1,int ID2,vector<int>& SupNode, int& d,int partiID){
	d=INF;
	unordered_map<int,int>::iterator it;
	int hub, dis1, dis2;
	for(it=LabelParti[partiID][ID1].begin();it!=LabelParti[partiID][ID1].end();it++){
		hub=(*it).first;
		dis1=(*it).second;
		if(LabelParti[partiID][ID2].find(hub)!=LabelParti[partiID][ID2].end()){
			dis2=LabelParti[partiID][ID2][hub];
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
}*/

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
	map<int,vector<int>> mapvec;
	mapvec.clear();
	BedgeAllPID.assign(nodenum, mapvec);
	map<int,pair<int,set<int>>> mp;
	mp.clear();
	BedgePID.assign(nodenum,mp);

	int ID1,ID2,weight;
	for(int k=0;k<pnum;k++){
		//boundary edges
		cout<<"k "<<k<<"////////////////////////////"<<endl;
		for(int i=0;i<BoundVer[k].size();i++){
			for(int j=i+1;j<BoundVer[k].size();j++){
				ID1=BoundVer[k][i];
				ID2=BoundVer[k][j];
				weight=HopQueryLocal(ID1,ID2,LabelParti[k]);
				BedgeAllPID[ID1][ID2].push_back(k);
				BedgeAllPID[ID2][ID1].push_back(k);

				if(OverlayGraph1[ID1].find(ID2)==OverlayGraph1[ID1].end()){
					OverlayGraph1[ID1].insert(make_pair(ID2,weight));
					OverlayGraph1[ID2].insert(make_pair(ID1,weight));
					BedgePID[ID1][ID2].first=weight;
					BedgePID[ID1][ID2].second.insert(k);
					BedgePID[ID2][ID1].first=weight;
					BedgePID[ID2][ID1].second.insert(k);
				}else if(OverlayGraph1[ID1][ID2]>weight){
					cout<<"situation b//////////////////"<<endl;
					OverlayGraph1[ID1][ID2]=weight;
					OverlayGraph1[ID2][ID1]=weight;
					BedgePID[ID1][ID2].first=weight;
					BedgePID[ID1][ID2].second.clear();
					BedgePID[ID1][ID2].second.insert(k);
					BedgePID[ID2][ID1].first=weight;
					BedgePID[ID2][ID1].second.clear();
					BedgePID[ID2][ID1].second.insert(k);
				}else if(OverlayGraph1[ID1][ID2]==weight){
					cout<<"situation c//////////////////"<<endl;
					BedgePID[ID1][ID2].second.insert(k);
					BedgePID[ID2][ID1].second.insert(k);
				}

			}
		}
	}

	for(int h=0;h<CutEdges.size();h++){
		ID1=CutEdges[h].first.first;
		ID2=CutEdges[h].first.second;
		weight=CutEdges[h].second;

		if(OverlayGraph1[ID1].find(ID2)!=OverlayGraph1[ID1].end()){
			cout<<"boundary edge is also cut edge "<<ID1<<" "<<ID2<<" ,partition ID "<<VtoParID[ID1]<<" "<<VtoParID[ID2]<<endl;
			if(OverlayGraph1[ID1][ID2]>weight){
				OverlayGraph1[ID1][ID2]=weight;
			}
		}else{
			OverlayGraph1[ID1].insert(make_pair(ID2,weight));
		}
	}

	for(int i=0;i<nodenum;i++){
		for(unordered_map<int,int>::iterator it=OverlayGraph1[i].begin();it!=OverlayGraph1[i].end();it++){
			NeighborsOverlay[i].push_back(make_pair((*it).first, (*it).second));
		}
	}
}*/

/*int Graph::DisQueryVallyParti(int ID1, int ID2, int PID){
	int neiID,neiDis;
	int d=INF;
	for(int i=0;i<NeighborsParti[PID][ID1].size();i++){
		neiID=NeighborsParti[PID][ID1][i].first;
		neiDis=NeighborsParti[PID][ID1][i].second;
		if(NodeOrder[neiID]<=NodeOrder[ID2] && LabelParti[PID][neiID].find(ID2)!=LabelParti[PID][neiID].end()){
			if(neiDis+LabelParti[PID][neiID][ID2]<d){
				d=neiDis+LabelParti[PID][neiID][ID2];
			}
		}
	}
	return d;
}*/

/*int Graph::DisQueryVallyOverlay(int ID1, int ID2){
	int neiID,neiDis;
	int d=INF;
	for(int i=0;i<NeighborsOverlay[ID1].size();i++){
		neiID=NeighborsOverlay[ID1][i].first;
		neiDis=NeighborsOverlay[ID1][i].second;
		if(NodeOrder[neiID]<=NodeOrder[ID2] && LabelOverlay[neiID].find(ID2)!=LabelOverlay[neiID].end()){
			if(neiDis+LabelOverlay[neiID][ID2]<d){
				d=neiDis+LabelOverlay[neiID][ID2];
			}
		}
	}
	return d;
}*/

/*int Graph::DisQueryVally(int ID1, int ID2, vector<vector<pair<int,int>>> &Neighbors, vector<unordered_map<int,int>> &Label){
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
}*/


/*int Graph::DisQueryPeakParti(int ID1, int ID2, int PID){
	int d=INF;
	unordered_map<int,int>::iterator it;
	int hub, dis1;
	for(it=LabelParti[PID][ID1].begin();it!=LabelParti[PID][ID1].end();it++){
		hub=(*it).first;
		dis1=(*it).second;
		if(hub!=ID2 && LabelParti[PID][ID2].find(hub)!=LabelParti[PID][ID2].end()){
			if(dis1+LabelParti[PID][ID2][hub]<d){
				d=dis1+LabelParti[PID][ID2][hub];
			}
		}
	}
	return d;
}*/

/*int Graph::DisQueryPeakOverlay(int ID1, int ID2){
	int d=INF;
	unordered_map<int,int>::iterator it;
	int hub, dis1;
	for(it=LabelOverlay[ID1].begin();it!=LabelOverlay[ID1].end();it++){
		hub=(*it).first;
		dis1=(*it).second;
		if(hub!=ID2 && LabelOverlay[ID2].find(hub)!=LabelOverlay[ID2].end()){
			if(dis1+LabelOverlay[ID2][hub]<d){
				d=dis1+LabelOverlay[ID2][hub];
			}
		}
	}
	return d;
}*/

/*int Graph::DisQueryPeak(int ID1, int ID2, vector<unordered_map<int,int>> &Label){
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
}*/

/*int Graph::DisQueryLower1Parti(int ID1, int ID2, int PID){
	int neiID,neiDis;
	int d=INF;
	for(int i=0;i<NeighborsParti[PID][ID1].size();i++){
		neiID=NeighborsParti[PID][ID1][i].first;
		neiDis=NeighborsParti[PID][ID1][i].second;
		if(NodeOrder[neiID]<NodeOrder[ID2] && LabelParti[PID][neiID].find(ID2)!=LabelParti[PID][neiID].end()){
			if(neiDis+LabelParti[PID][neiID][ID2]<d){
				d=neiDis+LabelParti[PID][neiID][ID2];
			}
		}
	}
	return d;
}*/

/*int Graph::DisQueryLower1Overlay(int ID1, int ID2){
	int neiID,neiDis;
	int d=INF;
	for(int i=0;i<NeighborsOverlay[ID1].size();i++){
		neiID=NeighborsOverlay[ID1][i].first;
		neiDis=NeighborsOverlay[ID1][i].second;
		if(NodeOrder[neiID]<NodeOrder[ID2] && LabelOverlay[neiID].find(ID2)!=LabelOverlay[neiID].end()){
			if(neiDis+LabelOverlay[neiID][ID2]<d){
				d=neiDis+LabelOverlay[neiID][ID2];
			}
		}
	}
	return d;
}*/

/*int Graph::DisQueryLower1(int ID1, int ID2, vector<vector<pair<int,int>>> &Neighbors, vector<unordered_map<int,int>> &Label){
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
}*/
