/*
 * read.cpp
 *
 *  Created on: 20 Sep 2022
 *      Author: zhangmengxuan
 */

#include "head.h"

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
	Neighbors.assign(nodenum, vecp);
	EdgeWei.clear();
	for(int i=0;i<num2;i++){
		inGraph>>ID1>>ID2>>weight;
		Neighbors[ID1].push_back(make_pair(ID2,weight));
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
        cout<<i<<" "<<vNodeOrder[i]<<" "<<Neighbors[vNodeOrder[i]].size()<<endl;
    }
}

void Graph::GraphPartitionRead(string filename){
	//read the partitioned graphs
	vector<pair<int,int>> vec; vec.clear();
	vector<vector<pair<int,int>>> vecvec;
	vecvec.assign(nodenum,vec);
//	NeighborsParti.assign(pnum, vecvec);

	//initialize edges' PID
	unordered_map<int,int> unomap;
	unomap.clear();
	EtoParID.assign(nodenum,unomap);

    VtoParID.assign(nodenum,-1);
    bool flag_minus = false;


    ifstream IF1(filename+"/subgraph_vertex");
    if(!IF1){
        cout<<"Cannot open file "<<"subgraph_vertex"<<endl;
        exit(1);
    }

    int pnum2;
    IF1>>pnum2;
    if(algoName == "NC"){
//        flag_minus = true;
        pnum = pnum2;
    }else if(algoName == "SC" || algoName == "MT"){
//        flag_minus = true;
//        pnum2 = pnum;
    }
    cout<<"Partition number: "<<pnum2<<endl;
    NeighborsParti.assign(pnum, vecvec);
    for(int k=0;k<pnum2;k++){
        int vernum,ID;
        IF1>>vernum;
        for(int i=0;i<vernum;i++){
            IF1>>ID;
//            if(flag_minus){
//                ID = ID-1;
//            }
            assert(ID>=0);
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

	ifstream IF(filename+"/subgraph_edge");
	if(!IF){
        cout<<"Cannot open file "<<"subgraph_edge"<<endl;
        exit(1);
    }

	int pnum1;
	IF>>pnum1;
//    if(algoName == "NC"){
//        pnum1 = pnum2;
//    }else if(algoName == "SC" || algoName == "MT"){
//        pnum1 = pnum;
//    }
	for(int k=0;k<pnum1;k++){
		int edgenum0,ID1,ID2,weight;
		IF>>edgenum0;
		for(int i=0;i<edgenum0;i++){
			IF>>ID1>>ID2>>weight;
//            if(flag_minus){
//                ID1 = ID1-1; ID2 = ID2-1;
//            }
            assert(ID1>=0 && ID1 <nodenum);
            assert(ID2>=0 && ID2 <nodenum);
			NeighborsParti[k][ID1].push_back(make_pair(ID2,weight));

			if(EtoParID[ID1].find(ID2)!=EtoParID[ID1].end())
				cout<<"something wrong: edge ("<<ID1<<", "<<ID2<<") alrealy in other partitions"<<endl;
			EtoParID[ID1].insert(make_pair(ID2,k));
		}
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
//        if(flag_minus){
//            ID1 = ID1-1; ID2 = ID2-1;
//        }

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


//function of reading update edges
void Graph::ReadUpdates(string filename){
    int ID1, ID2, weight;
    ifstream inFile(filename, ios::in);
    if (!inFile) { // if not exist
        cout << "Fail to open file" << filename << endl;
        exit(1);
    }
    int num;
    inFile >> num;
    for(int i=0;i<num;i++){
        inFile>>ID1>>ID2>>weight;
        updateEdges.emplace_back(make_pair(ID1, ID2), weight);
    }
    inFile.close();
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
	for(int k=0;k<num;k++){
		IF>>ID1>>ID2;
		ODpair.push_back(make_pair(ID1, ID2));
	}

	int s, t;
	std::chrono::high_resolution_clock::time_point t1, t2;
	std::chrono::duration<double> time_span;
	double runT;
    if(runtimes > ODpair.size()){
        runtimes = ODpair.size();
    }
    cout<<"Run times: "<<runtimes<<endl;
	t1=std::chrono::high_resolution_clock::now();
	for(int i=0;i<ODpair.size();i++){
		HopQuery(ODpair[i].first,ODpair[i].second);
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"average query time: "<<runT*1000/ODpair.size()<<" ms."<<endl;
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
        d1=HopQueryOverlay(ODpair1[i].first,ODpair1[i].second);
        t2=std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
        runT+= time_span.count();
        //Timer tt;
        //tt.start();
        //d2=Dijkstra(ODpair1[i].first,ODpair1[i].second);
        //tt.stop();
        //dijkT+=tt.GetRuntime();
        //if(d1 != d2){
          //  cout<<"Incorrect! "<<ODpair1[i].first<<" "<<ODpair1[i].second<<" "<<d1<<" "<<d2<<endl;
        //}
    }
    cout<<"both OD are boundary: "<<1000*runT/ODpair1.size()<<" ms."<<endl;
    //cout<<"both OD are boundary (Dijkstra): "<<1000*dijkT/ODpair1.size()<<" ms."<<endl;

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

    runT=0; dijkT=0;
    for(int i=0;i<ODpair4.size();i++){
        t1=std::chrono::high_resolution_clock::now();
        d1=HopQueryInOut(ODpair4[i].first,ODpair4[i].second);
        t2=std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
        runT += time_span.count();
        /*Timer tt;
        tt.start();
        d2=Dijkstra(ODpair4[i].first,ODpair4[i].second);
        tt.stop();
        dijkT+=tt.GetRuntime();
        if(d1 != d2){
            cout<<"Incorrect! "<<ODpair4[i].first<<" "<<ODpair4[i].second<<" "<<d1<<" "<<d2<<endl;
        }*/
    }

    cout<<"one boundary, one within partition: "<<1000*runT/ODpair4.size()<<" ms."<<endl;
    //cout<<"one boundary, one within partition (Dijkstra): "<<1000*dijkT/ODpair4.size()<<" ms."<<endl;
}

//overlay graph construction
void Graph::OverlayGraphConstructPost(){
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
				weight=QueryH2HPartition(ID1,ID2,k);
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

void Graph::CorrectnessCheck(int runtimes){
	srand (time(NULL));
	int s, t, d1, d2, d3;
	for(int i=0;i<runtimes;i++){
		//cout<<i<<endl;
		s=rand()%nodenum;
		t=rand()%nodenum;
		d1=Dijkstra(s,t);
		d2=HopQuery(s,t);
		if(d1!=d2){
			cout<<"InCorrect!"<<s<<" "<<t<<" "<<d1<<" "<<d2<<endl;
		}
	}
}

int Graph::HopQuery(int ID1, int ID2){
	int d=INF;

	if(TotalBoundSet.find(ID1)!=TotalBoundSet.end() && TotalBoundSet.find(ID2)!=TotalBoundSet.end()){
		d=HopQueryOverlay(ID1,ID2);
	}else if(TotalBoundSet.find(ID1)!=TotalBoundSet.end()){
		d=HopQueryInOut(ID2,ID1);
	}else if(TotalBoundSet.find(ID2)!=TotalBoundSet.end()){
		d=HopQueryInOut(ID1,ID2);
	}else{
		if(VtoParID[ID1]==VtoParID[ID2]){//the same partition
			d=HopQueryInIn(ID1,ID2);
		}else{//different partitions
			d=HopQueryOutOut(ID1,ID2);
		}
	}

	return d;
}

int Graph::HopQueryOverlay(int ID1, int ID2){
	if(ID1==ID2) return 0;
	if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
	int r1=rankOverlay[ID1], r2=rankOverlay[ID2];
	int LCA=LCAQueryOverlay(r1,r2);

	if(LCA==r1)
		return TreeOverlay[r2].dis[TreeOverlay[r1].pos.back()];
	else if(LCA==r2)
		return TreeOverlay[r1].dis[TreeOverlay[r2].pos.back()];
	else{
		int tmp=INF;
		for(int i=0;i<TreeOverlay[LCA].pos.size();i++){
			if(tmp>TreeOverlay[r1].dis[TreeOverlay[LCA].pos[i]]+TreeOverlay[r2].dis[TreeOverlay[LCA].pos[i]])
				tmp=TreeOverlay[r1].dis[TreeOverlay[LCA].pos[i]]+TreeOverlay[r2].dis[TreeOverlay[LCA].pos[i]];
		}
		return tmp;
	}
}

int Graph::LCAQueryOverlay(int _p, int _q){
	int p = toRMQOverlay[_p], q = toRMQOverlay[_q];
	if (p > q){
		int x = p;
		p = q;
		q = x;
	}
	int len = q - p + 1;
	int i = 1, k = 0;
	while (i * 2 < len){
		i *= 2;
		k++;
	}
	q = q - i + 1;
	if (TreeOverlay[RMQIndexOverlay[k][p]].height < TreeOverlay[RMQIndexOverlay[k][q]].height)
		return RMQIndexOverlay[k][p];
	else return RMQIndexOverlay[k][q];
}

//Query within one partition
int Graph::QueryH2HPartition(int ID1, int ID2, int PID){
	if(ID1==ID2) return 0;
	int r1=ranks[PID][ID1], r2=ranks[PID][ID2];
	int LCA=LCAQueryPartition(r1,r2,PID);

	if(LCA==r1)
		return Trees[PID][r2].dis[Trees[PID][r1].pos.back()];
	else if(LCA==r2)
		return Trees[PID][r1].dis[Trees[PID][r2].pos.back()];
	else{
		int tmp=INF;
		for(int i=0;i<Trees[PID][LCA].pos.size();i++){
			if(tmp>Trees[PID][r1].dis[Trees[PID][LCA].pos[i]]+Trees[PID][r2].dis[Trees[PID][LCA].pos[i]])
				tmp=Trees[PID][r1].dis[Trees[PID][LCA].pos[i]]+Trees[PID][r2].dis[Trees[PID][LCA].pos[i]];
		}
		return tmp;
	}
}

int Graph::LCAQueryPartition(int _p, int _q, int PID){
	int p = toRMQs[PID][_p], q = toRMQs[PID][_q];
	if (p > q){
		int x = p;
		p = q;
		q = x;
	}
	int len = q - p + 1;
	int i = 1, k = 0;
	while (i * 2 < len){
		i *= 2;
		k++;
	}
	q = q - i + 1;
	if (Trees[PID][RMQIndexs[PID][k][p]].height < Trees[PID][RMQIndexs[PID][k][q]].height)
		return RMQIndexs[PID][k][p];
	else return RMQIndexs[PID][k][q];
}


int Graph::HopQueryInOut(int in, int out){
	int d=INF;

	int boundid, tempdis;
	for(int i=0;i<BoundVer[VtoParID[in]].size();i++){
		boundid=BoundVer[VtoParID[in]][i];
		tempdis=QueryH2HPartition(in, boundid, VtoParID[in])+HopQueryOverlay(boundid,out);
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
		m1.insert(make_pair(bID1,QueryH2HPartition(ID1,bID1,VtoParID[ID1])));
	}
	for(int j=0;j<B2.size();j++){
		bID2=B2[j];
		m2.insert(make_pair(bID2,QueryH2HPartition(ID2,bID2,VtoParID[ID2])));
	}

	for(int k=0;k<B1.size();k++){
		bID1=B1[k];

		if(m1[bID1]>d)
			continue;

		for(int z=0;z<B2.size();z++){
			bID2=B2[z];

			if(m2[bID2]>d)
				continue;

			tempdis=m1[bID1]+HopQueryOverlay(bID1,bID2)+m2[bID2];
			if(tempdis<d)
				d=tempdis;
		}
	}

	return d;
}

int Graph::HopQueryInIn(int ID1, int ID2){
	int d=INF;

	int pid=VtoParID[ID1];
	int d0=QueryH2HPartition(ID1,ID2,pid);
	if(d0<d)
		d=d0;

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
		d1=QueryH2HPartition(ID1,bID,pid);
		d2=QueryH2HPartition(ID2,bID,pid);
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
			tempdis=m1[bID1]+HopQueryOverlay(bID1,bID2)+m2[bID2];
			if(tempdis<d)
				d=tempdis;
		}
	}

	return d;
}

void Graph::Indexsize(){
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
    cout<<"Index size: "<<(double)m/1024/1024<<" MB"<<endl;
}
