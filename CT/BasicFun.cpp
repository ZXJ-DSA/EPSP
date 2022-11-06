/*
 * BasicFun.cpp
 *
 *  Created on: 13 Oct 2022
 *      Author: zhangmengxuan
 */
#include "head.h"

void Graph::ReadGraph(string filename){
	ifstream inGraph(filename);
	if(!inGraph){
		cout<<"Cannot open Map "<<filename<<endl;
        exit(1);
	}

	string line;
//	do{
//		getline(inGraph,line);
//		if(line[0]=='p'){
//			vector<string> vs;
//			boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
//			nodenum=stoi(vs[2]); edgenum=0;
//		}
//	}while(line[0]=='c'||line[0]=='p');

    getline(inGraph,line);
    vector<string> vs;
    boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
    nodenum=stoi(vs[0]); edgenum=0;
    getline(inGraph,line);

	//graph g initialize
	vector<pair<int,int>> vecp; vecp.clear();
	Neighbor.assign(nodenum, vecp);
	set<int> m; m.clear();
	vector<set<int>> E;
	E.assign(nodenum,m);

	int ID1,ID2, weight;
	while(!line.empty()){
		vector<string> vs;
		boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
		ID1=stoi(vs[0]); ID2=stoi(vs[1]); weight=stoi(vs[2]);

		if(E[ID1].find(ID2)==E[ID1].end()){
			edgenum+=2;
			Neighbor[ID1].push_back(make_pair(ID2,weight));
			Neighbor[ID2].push_back(make_pair(ID1,weight));
			E[ID1].insert(ID2);
			E[ID2].insert(ID1);
		}
		if(inGraph.eof()) break;
		getline(inGraph,line);
	}
	cout<<"Finish Reading! nodenum "<<nodenum<<endl;
}

void Graph::UpdateGene(int num, string filename){
	vector<pair<pair<int,int>, pair<int,int>>> UpdateData;

	set<pair<int,int>> Edges;
	vector<pair<pair<int,int>,int>> ENodeID;
	int ID1,ID2,wei;
	for(int i=0;i<Neighbor.size();i++){
		ID1=i;
		for(int j=0;j<Neighbor[i].size();j++){
			ID2=Neighbor[i][j].first;
			wei=Neighbor[i][j].second;
			if(ID1<ID2 && Edges.find(make_pair(ID1,ID2))==Edges.end()){
				Edges.insert(make_pair(ID1,ID2));
				ENodeID.push_back(make_pair(make_pair(ID1,ID2),wei));
			}
			else if(ID2<ID1 && Edges.find(make_pair(ID2,ID1))==Edges.end()){
				Edges.insert(make_pair(ID2,ID1));
				ENodeID.push_back(make_pair(make_pair(ID2,ID1),wei));
			}
		}
	}

	ofstream OF(filename);
	OF<<num<<endl;
	set<int> eid;
	for(int k=0;k<num;k++){
		int edgeid=rand()%ENodeID.size();
		if(eid.find(edgeid)==eid.end()){
			OF<<ENodeID[edgeid].first.first<<" "<<ENodeID[edgeid].first.second<<" "<<ENodeID[edgeid].second<<endl;
			eid.insert(edgeid);
		}else{
			k--;
		}
	}
	OF.close();
}

void Graph::ODGene(int num, string filename){
	set<pair<int,int>> ODpair;
	vector<pair<int,int>> ODpairVec;

	srand (time(NULL));
	int s, t;
	for(int i=0;i<num;i++){
		s=rand()%nodenum;
		t=rand()%nodenum;
		if(ODpair.find(make_pair(s,t))==ODpair.end()){
			ODpairVec.push_back(make_pair(s,t));
			ODpair.insert(make_pair(s,t));
			ODpair.insert(make_pair(t,s));
		}else{
			i--;
		}
	}
	cout<<"generated OD pair number "<<ODpairVec.size()<<endl;

	ofstream OF(filename);
	OF<<ODpairVec.size()<<endl;
	for(int k=0;k<ODpairVec.size();k++){
		OF<<ODpairVec[k].first<<" "<<ODpairVec[k].second<<endl;
	}
	OF.close();
}

void Graph::StainingMethod(int ID){
	queue<int> Q;

	vector<bool> Stained;
	Stained.assign(nodenum, false);

	Q.push(ID);
	Stained[ID]=true;
	int frontid, neiid;
	while(!Q.empty()){
		frontid=Q.front();
		Q.pop();
		for(int k=0;k<Neighbor[frontid].size();k++){
			neiid=Neighbor[frontid][k].first;
			if(!Stained[neiid]){
				Q.push(neiid);
				Stained[neiid]=true;
			}
		}
	}

	int stainNum=0;
	for(int i=0;i<nodenum;i++){
		if(Stained[i])
			stainNum+=1;
	}
	//cout<<"Stained Number "<<stainNum<<endl;

	vector<int> VertexInverted;
	VertexInverted.assign(nodenum, -1);
	int j=0;
	for(int i=0;i<nodenum;i++){
		if(Stained[i]){
			VertexInverted[i]=j;
			j+=1;
		}
	}
	//cout<<"Check j= "<<j<<", stainNum= "<<stainNum<<endl;

	int OriNodeNum=nodenum;
	nodenum=stainNum;
	vector<vector<pair<int,int>>> Neighbor1=Neighbor;
	Neighbor.clear();
	vector<pair<int,int>> vecpair;
	vecpair.clear();
	Neighbor.assign(nodenum, vecpair);
	int InvertedID, nei, Invertednei, wei;
	for(int ID=0;ID<OriNodeNum;ID++){
		if(VertexInverted[ID]!=-1){
			InvertedID=VertexInverted[ID];
			for(int k=0;k<Neighbor1[ID].size();k++){
				nei=Neighbor1[ID][k].first;
				wei=Neighbor1[ID][k].second;
				if(VertexInverted[nei]!=-1){
					Invertednei=VertexInverted[nei];
					Neighbor[InvertedID].push_back(make_pair(Invertednei,wei));
				}
			}
		}
	}

}

int Graph::Dijkstra(int ID1, int ID2,vector<vector<pair<int,int>>> &Neighbor){
	benchmark::heap<2, int, int> pqueue(nodenum);
	pqueue.update(ID1,0);

	vector<bool> closed(nodenum, false);
	vector<int> distance(nodenum, INF);
	vector<int> prece(nodenum, 0);
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
					prece[NNodeID]=topNodeID;
				}
			}
		}
	}

	return d;
}

int Graph::DijkstraPath(int ID1, int ID2){
	benchmark::heap<2, int, int> pqueue(nodenum);
	pqueue.update(ID1,0);

	vector<bool> closed(nodenum, false);
	vector<int> distance(nodenum, INF);
	vector<int> prece(nodenum, 0);
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
					prece[NNodeID]=topNodeID;
				}
			}
		}
	}

	//path retrieval
	vector<int> path;
	path.clear();
	path.push_back(ID2);
	int preID=prece[ID2];
	while(preID!=ID1){
		path.push_back(preID);
		preID=prece[preID];
	}
	path.push_back(ID1);

	cout<<"path: ";
	for(int i=path.size()-1;i>-1;i--){
		cout<<" "<<path[i]<<" "<<CoreTag[path[i]]<<" "<<BoundTag[path[i]]<<" "<<NodeOrder[path[i]]<<endl;
	}
	cout<<endl;

	return d;
}

void Graph::PartitionPreProcess(){//consider the case when there is no children for boundary tree node
	H2HconCore();//construct core

	//identify the partition from tree-periphery
	vector<int> cr;//the vertex rank of root's children
	cr.clear();
    int temp_ch=0;
    int temp_vert=0;
	for(int k=0;k<TreeCore[0].ch.size();k++){//TreeCore[0] is the core
		int childrank=TreeCore[0].ch[k];

		if(TreeCore[childrank].ch.size()>0){
			cr.push_back(childrank);
            temp_ch += TreeCore[childrank].ch.size();
            temp_vert += TreeCore[childrank].vert.size();
		}
	}
//    cout<<"TreeCore[0].ch.size(): "<<TreeCore[0].ch.size()<<endl;
//    cout<<"Accumulated TreeCore[childrank].ch.size(): "<<temp_ch<<endl;
//    cout<<"Accumulated TreeCore[childrank].vert.size(): "<<temp_vert<<endl;
	partiNum=cr.size();

	cout<<"Periphery Number: "<<partiNum<<endl;

	vector<int> vec;
	vec.clear();
	BoundVertex.assign(partiNum,vec);
	set<int> sset;
	sset.clear();
	BoundVertexSet.assign(partiNum,sset);
	BoundTag.assign(nodenum,0);
	map<int,int> PartiRoot;//partition root & partition ID
	PartiRoot.clear();
	for(int PID=0;PID<partiNum;PID++){
		int childrank=cr[PID];
		for(int i=0;i<TreeCore[childrank].vert.size();i++){
			BoundVertex[PID].push_back(TreeCore[childrank].vert[i].first);
			BoundVertexSet[PID].insert(TreeCore[childrank].vert[i].first);
			BoundTag[TreeCore[childrank].vert[i].first]=1;
		}
		BoundVertex[PID].push_back(TreeCore[childrank].uniqueVertex);
		BoundVertexSet[PID].insert(TreeCore[childrank].uniqueVertex);
		BoundTag[TreeCore[childrank].uniqueVertex]=1;

		PartiRoot.insert(make_pair(childrank,PID));
	}

	//if(PartiRoot.size()!=cr.size())
		//cout<<"something wrong with the boundary node"<<endl;

	CoreTag.assign(nodenum,-1);
	int NodeID,RootNode,parentNode;
	for(int len=HighestOrder-1;len>=0;len--){
		NodeID=vNodeOrder[len];
		RootNode=TreeCore[rankCore[NodeID]].treeroot;
		parentNode=TreeCore[rankCore[NodeID]].pa;
		if(parentNode!=0){
			CoreTag[NodeID]=PartiRoot[RootNode];
		}
	}
}

/*void Graph::PartitionPreProcess(){
	H2HconCore();

	//identify the partition from tree-periphery
	partiNum=TreeCore[0].ch.size();

	cout<<"Periphery Number "<<partiNum<<endl;

	vector<int> vec;
	vec.clear();
	BoundVertex.assign(partiNum,vec);
	set<int> sset;
	sset.clear();
	BoundVertexSet.assign(partiNum,sset);
	BoundTag.assign(nodenum,0);
	map<int,int> PartiRoot;//partition root & partition ID
	PartiRoot.clear();
	for(int PID=0;PID<partiNum;PID++){
		int childrank=TreeCore[0].ch[PID];

		//if(TreeCore[childrank].ch.size()==0)
			//cout<<"there is no point in partition "<<PID<<endl;

		if(TreeCore[childrank].ch.size()>0){
			for(int i=0;i<TreeCore[childrank].vert.size();i++){
				BoundVertex[PID].push_back(TreeCore[childrank].vert[i].first);
				BoundVertexSet[PID].insert(TreeCore[childrank].vert[i].first);
				BoundTag[TreeCore[childrank].vert[i].first]=1;
			}
			BoundVertex[PID].push_back(TreeCore[childrank].uniqueVertex);
			BoundVertexSet[PID].insert(TreeCore[childrank].uniqueVertex);
			BoundTag[TreeCore[childrank].uniqueVertex]=1;

			PartiRoot.insert(make_pair(childrank,PID));
		}

	}

	CoreTag.assign(nodenum,-1);
	int NodeID,RootNode,parentNode;
	for(int len=HighestOrder-1;len>=0;len--){
		NodeID=vNodeOrder[len];
		RootNode=TreeCore[rankCore[NodeID]].treeroot;
		parentNode=TreeCore[rankCore[NodeID]].pa;
		if(parentNode!=0){
			CoreTag[NodeID]=PartiRoot[RootNode];
		}
	}
}*/

void Graph::WritePartition(string filename){
	//********Boundary vertex********//
	ofstream OF1(filename+"Boundary");
	OF1<<partiNum<<endl;
	for(int pid=0;pid<partiNum;pid++){
		OF1<<pid<<" "<<BoundVertex[pid].size();
		for(int i=0;i<BoundVertex[pid].size();i++)
			OF1<<" "<<BoundVertex[pid][i];
		OF1<<endl;
	}
	OF1.close();

	//*********vertex tag*********//
	ofstream OF2(filename+"Tag");
	OF2<<nodenum<<endl;
	for(int i=0;i<nodenum;i++){
		OF2<<i<<" "<<CoreTag[i]<<" "<<BoundTag[i]<<" "<<NodeOrder[i]<<endl;
	}
	OF2.close();
}

void Graph::ReadPartition(string filename){
	//********Boundary vertex********//
	ifstream IF1(filename+"Boundary");
	int pnum, pid, bnum, bid;
	IF1>>pnum;
	vector<int> vec;
	vec.clear();
	partiNum=pnum;
	BoundVertex.assign(pnum,vec);
	IF1>>pid>>bnum;
	for(int k=0;k<bnum;k++){
		IF1>>bid;
		BoundVertex[pid].push_back(bid);
	}
	IF1.close();

	//*********vertex tag*********//
	ifstream IF2(filename+"Tag");
	int num,id,core,bound,order;
	IF2>>num;
	CoreTag.assign(num,0);
	BoundTag.assign(num,0);
	NodeOrder.assign(num,0);
	for(int k=0;k<num;k++){
		IF2>>id>>core>>bound>>order;
		CoreTag[id]=core;
		BoundTag[id]=bound;
		NodeOrder[id]=order;
	}
	IF2.close();
}

int Graph::DijkstraCore(int ID1, int ID2){
	benchmark::heap<2, int, int> pqueue(nodenum);
	pqueue.update(ID1,0);

	vector<bool> closed(nodenum, false);
	vector<int> distance(nodenum, INF);
	vector<int> prece(nodenum, 0);
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

		for(it=AdjaCore[topNodeID].begin();it!=AdjaCore[topNodeID].end();it++){
			NNodeID=(*it).first;
			NWeigh=(*it).second+topNodeDis;
			if(!closed[NNodeID]){
				if(distance[NNodeID]>NWeigh){
					distance[NNodeID]=NWeigh;
					pqueue.update(NNodeID, NWeigh);
					prece[NNodeID]=topNodeID;
				}
			}
		}
	}

	return d;
}

int Graph::DijkstraCorePath(int ID1, int ID2){
	benchmark::heap<2, int, int> pqueue(nodenum);
	pqueue.update(ID1,0);

	vector<bool> closed(nodenum, false);
	vector<int> distance(nodenum, INF);
	vector<int> prece(nodenum, 0);
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

		for(it=AdjaCore[topNodeID].begin();it!=AdjaCore[topNodeID].end();it++){
			NNodeID=(*it).first;
			NWeigh=(*it).second+topNodeDis;
			if(!closed[NNodeID]){
				if(distance[NNodeID]>NWeigh){
					distance[NNodeID]=NWeigh;
					pqueue.update(NNodeID, NWeigh);
					prece[NNodeID]=topNodeID;
				}
			}
		}
	}

	//path retrieval
	vector<int> path;
	path.clear();
	path.push_back(ID2);
	int preID=prece[ID2];
	while(preID!=ID1){
		path.push_back(preID);
		preID=prece[preID];
	}
	path.push_back(ID1);

	cout<<"Core path: ";
	for(int i=path.size()-1;i>-1;i--){
		cout<<" "<<path[i]<<" "<<CoreTag[path[i]]<<" "<<BoundTag[path[i]]<<endl;
	}
	cout<<endl;

	return d;
}

////Query processing
//int Graph::Query(int ID1, int ID2){
//	int dis=INF;
//
//	if(CoreTag[ID1]==-1 && CoreTag[ID2]==-1){//both in core
//		//cout<<"1"<<endl;
//		dis=QueryCore(ID1, ID2);
//	}else if(CoreTag[ID1]==-1 && CoreTag[ID2]!=-1){//ID2 in partition, ID1 in core
//		//cout<<"2"<<endl;
//		dis=QueryPartiCore(ID2, ID1);
//	}else if(CoreTag[ID1]!=-1 && CoreTag[ID2]==-1){//ID1 in partition, ID2 in core
//		//cout<<"3"<<endl;
//		dis=QueryPartiCore(ID1, ID2);
//	}else if(CoreTag[ID1]!=-1 && CoreTag[ID2]!=-1){//both in partition
//		//cout<<"4"<<endl;
//		dis=QueryPartiParti(ID1,ID2);
//	}
//	return dis;
//}

int Graph::QueryCore(int ID1, int ID2){
	int hubfinal,dis1final,dis2final;
	int d=INF;
	unordered_map<int,int>::iterator it;
	int hub, dis1, dis2;
	for(it=IndexCore[ID1].begin();it!=IndexCore[ID1].end();it++){
		hub=(*it).first;
		dis1=(*it).second;
		if(IndexCore[ID2].find(hub)!=IndexCore[ID2].end()){
			dis2=IndexCore[ID2][hub];
			if(dis1+dis2<d){
				d=dis1+dis2;
				hubfinal=hub;
				dis1final=dis1;
				dis2final=dis2;
				//cout<<"hub "<<hub<<",dis "<<d<<endl;
				//cout<<"check labeling "<<dis1<<" "<<DijkstraCore(ID1,hub)<<" "<<dis2<<" "<<DijkstraCore(ID1,hub)<<endl;
			}
		}
	}

	//cout<<"hub "<<hubfinal<<",dis "<<d<<endl;
	//cout<<"check labeling "<<dis1final<<" "<<DijkstraCore(ID1,hubfinal)<<" "<<dis2final<<" "<<DijkstraCore(ID2,hubfinal)<<endl;

	return d;
}

void Graph::CorrectnessCheck(){
	srand (time(NULL));
	int s, t, d1, d2, d3;
    cout<<"Correctness check..."<<endl;
	for(int i=0;i<100;i++){
		s=rand()%nodenum;
		t=rand()%nodenum;
		d1=Dijkstra(s,t,Neighbor);
		d2=Query(s,t);
		if(d1!=d2){
			cout<<"InCorrect!"<<s<<" "<<t<<" "<<d1<<" "<<d2<<endl;
		}
	}
}

void Graph::CorrectnessCheckCore(){
	srand (time(NULL));
	int s, t, d1, d2, d3;
	for(int i=0;i<1000;i++){
		s=rand()%nodenum;
		t=rand()%nodenum;
		if(CoreTag[s]==-1 && CoreTag[t]==-1){
			d1=Dijkstra(s,t,Neighbor);
			d2=DijkstraCore(s,t);
			if(d1!=d2){
				cout<<"InCorrect!"<<s<<" "<<t<<" "<<d1<<" "<<d2<<endl;
				DijkstraPath(s,t);
				DijkstraCorePath(s,t);
			}
		}else
			i--;
	}
}

void Graph::TestDataRead(string filename, double Ratio, vector<pair<pair<int,int>,pair<int,int>>>& Data){
	ifstream IF(filename);
	int num;
	int a,b,w,newW;
	IF>>num;
	for(int cnum=0;cnum<num;cnum++){
		IF>>a>>b>>w;
		newW=Ratio*w;
		Data.push_back(make_pair(make_pair(a,b), make_pair(w, newW)));
	}
}

void Graph::ReadUpdate(string filename,vector<pair<pair<int,int>,int>>& TestData){
    TestData.clear();

    int num, ID1, ID2, oldw;
    ifstream IF(filename);
    if(!IF){
        cout<<"Cannot open file "<<filename<<endl;
        exit(1);
    }
    IF>>num;
    for(int i=0;i<num;i++){
        IF>>ID1>>ID2>>oldw;
        TestData.push_back(make_pair(make_pair(ID1, ID2), oldw));
    }
    IF.close();
}

void Graph::EffiCheck(string filename,int runtimes){
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
		Query(ODpair[i].first,ODpair[i].second);
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"Average Query Time: "<<(double)runT*1000/runtimes<<" ms."<<endl;
}

/*void Graph::WriteCoreGraph(string graphfile){
	ofstream OF(graphfile);
	OF<<"p"<<" "<<0<<" "<<nodenum<<endl;
	for(int ID1=0;ID1<nodenum;ID1++){
		for(int j=0;j<AdjaCore[ID1].size();j++){
			int ID2=AdjaCore[ID1][j].first;
			int wei=AdjaCore[ID1][j].second;
			OF<<ID1<<" "<<ID2<<" "<<wei<<endl;
		}
	}
}

void Graph::ReadCoreGraph(string filename){
	ifstream inGraph(filename);
	if(!inGraph){
		cout<<"Cannot open Map "<<filename<<endl;
	}

	string line;
	do{
		getline(inGraph,line);
		if(line[0]=='p'){
			vector<string> vs;
			boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
			nodenum=stoi(vs[2]); edgenum=0;
		}
	}while(line[0]=='c'||line[0]=='p');

	//graph g initialize
	vector<pair<int,int>> vecp; vecp.clear();
	Neighbors.assign(nodenum, vecp);
	set<int> m; m.clear();
	E.assign(nodenum,m);

	int ID1,ID2, weight;
	while(!line.empty()){
		vector<string> vs;
		boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
		ID1=stoi(vs[0]); ID2=stoi(vs[1]); weight=stoi(vs[2]);

		if(E[ID1].find(ID2)==E[ID1].end()){
			edgenum+=2;
			Neighbors[ID1].push_back(make_pair(ID2,weight));
			Neighbors[ID2].push_back(make_pair(ID1,weight));
			E[ID1].insert(ID2);
			E[ID2].insert(ID1);
			if(ID1<ID2) EdgeEnum.push_back(make_pair(make_pair(ID1,ID2), weight));
			else EdgeEnum.push_back(make_pair(make_pair(ID2,ID1), weight));
		}
		if(inGraph.eof()) break;
		getline(inGraph,line);
	}
	cout<<"Finish Graph Reading! Nnum "<<nodenum<<" Enum "<<edgenum<<endl;

	benchmark::heap<2, int, int> Q(nodenum);
	for(int i=0;i<nodenum;i++){
		Q.update(i,Neighbors[i].size());
		if(Neighbors.size()==0) cout<<"Isolated node "<<i<<endl;
	}
	int topID, topdegree;
	int cnt=0;
	NodeOrder.assign(nodenum, 0);
	vNodeOrder.assign(nodenum, 0);
	while(!Q.empty()){
		Q.extract_min(topID, topdegree);
		NodeOrder[topID]=cnt;
		vNodeOrder[cnt]=topID;
		cnt+=1;
	}
	cout<<"Node finish ordering!"<<endl;
}*/
