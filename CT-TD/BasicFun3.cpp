/*
 * BasicFun.cpp
 *
 *  Created on: 13 Oct 2022
 *      Author: zhangmengxuan
 */
#include "head3.h"

void Graph::ReadGraph(string filename){
	ifstream inGraph(filename);
	if(!inGraph){
		cout<<"Cannot open Map "<<filename<<endl;
        exit(1);
	}

	string line;

    getline(inGraph,line);
    vector<string> vs;
    boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
    nodenum=stoi(vs[0]); edgenum=0;
    int tempENum=stoi(vs[1]);
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

		if(E[ID1].find(ID2)==E[ID1].end()){//if not found
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
    assert(edgenum == tempENum);
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
    if(stainNum != nodenum){
        cout<<"Incorrect!!! stain number: "<<stainNum<<" ; node number: "<<nodenum<<endl;
    }

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

	cout<<"path from "<<ID1<<" to "<<ID2<<": "<<endl;
	for(int i=path.size()-1;i>-1;i--){
		cout<<" "<<path[i]<<"("<<CoreTag[path[i]]<<","<<BoundTag[path[i]].first<<","<<NodeOrder[path[i]]<<") ";//<<endl;
        if(i>0){
            for(int j=0;j<Neighbor[path[i]].size();++j){
                if(Neighbor[path[i]][j].first == path[i-1]){
                    cout<<Neighbor[path[i]][j].second<<endl;
                    break;
                }
            }
        }
	}
	cout<<endl;

	return d;
}


int Graph::DijkstraPathCore(int ID1, int ID2){
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

    pair<int,int> highestVertex(-1,0);//ID, order
    cout<<"path between "<<ID1<<" and "<<ID2<<": "<<endl;
    for(int i=path.size()-1;i>-1;i--){
        cout<<path[i]<<" "<<CoreTag[path[i]]<<" "<<NodeOrder[path[i]];
        if(NodeOrder[path[i]] > highestVertex.second){
            highestVertex.second = NodeOrder[path[i]];
            highestVertex.first = path[i];
        }
        if(i>0){
            for(int j=0;j<AdjaCore[path[i]].size();++j){
                if(AdjaCore[path[i]][j].first == path[i-1]){
                    cout<<"; "<<AdjaCore[path[i]][j].second;
                    break;
                }
            }
        }
        cout<<endl;
    }
    cout<<"Highest-order vertex: "<<highestVertex.first<<" ("<<highestVertex.second<<")"<<endl;

    return d;
}

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
		OF2<<i<<" "<<CoreTag[i]<<" "<<BoundTag[i].first<<" "<<NodeOrder[i]<<endl;
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
	BoundTag.assign(num,make_pair(false,set<int>()));
	NodeOrder.assign(num,0);
	for(int k=0;k<num;k++){
		IF2>>id>>core>>bound>>order;
		CoreTag[id]=core;
		BoundTag[id].first=bound;
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

	cout<<"Core path: "<<endl;
	for(int i=path.size()-1;i>-1;i--){
		cout<<" "<<path[i]<<" "<<CoreTag[path[i]]<<" "<<NodeOrder[path[i]]<<endl;
	}
	cout<<endl;

	return d;
}

int Graph::QueryCoreDebug(int ID1, int ID2){
    int hubfinal,dis1final,dis2final;
    int d=INF;
    unordered_map<int,int>::iterator it;
    int hub, dis1, dis2, temp;
    if(CoreTag[ID1]!=-1 ||CoreTag[ID2]!=-1){
        cout<<"Core Tag wrong! "<<CoreTag[ID1]<< " "<<CoreTag[ID2]<<endl;
        exit(1);
    }
//    cout<<NodeOrder[ID1]<<" "<<NodeOrder[ID2]<<endl;
    for(it=Label[ID1].begin();it!=Label[ID1].end();it++){
        hub=(*it).first;
        dis1=(*it).second;
        if(hub == 207156)//146447
            cout<<"hub: "<<hub<<endl;
        if(Label[ID2].find(hub)!=Label[ID2].end()){
            dis2=Label[ID2][hub];
            temp=dis1+dis2;
            cout<<"temp hub "<<hub<<",dis "<<temp<<" "<<d<<endl;
            if(temp<d){
                d=temp;
                hubfinal=hub;
                dis1final=dis1;
                dis2final=dis2;

                //cout<<"check labeling "<<dis1<<" "<<DijkstraCore(ID1,hub)<<" "<<dis2<<" "<<DijkstraCore(ID1,hub)<<endl;
            }
        }
    }

    cout<<ID1<<" "<<ID2<<": hub "<<hubfinal<<"("<<NodeOrder[hubfinal]<<"), dis "<<d<<endl;
    cout<<"check labeling "<<dis1final<<" "<<DijkstraCore(ID1,hubfinal)<<"("<<DijkstraCore(hubfinal,ID1)<<") "<<dis2final<<" "<<DijkstraCore(ID2,hubfinal)<<"("<<DijkstraCore(hubfinal,ID2)<<")"<<endl;
    DijkstraPathCore(ID1,ID2);
    return d;
}

int Graph::QueryCore(int ID1, int ID2){
	int hubfinal,dis1final,dis2final;
	int d=INF;
	unordered_map<int,int>::iterator it;
	int hub, dis1, dis2;
    if(CoreTag[ID1]!=-1 ||CoreTag[ID2]!=-1){
        cout<<"Core Tag wrong! "<<CoreTag[ID1]<< " "<<CoreTag[ID2]<<endl;
        exit(1);
    }
//    cout<<NodeOrder[ID1]<<" "<<NodeOrder[ID2]<<endl;
	for(it=Label[ID1].begin();it!=Label[ID1].end();it++){
		hub=(*it).first;
		dis1=(*it).second;
		if(Label[ID2].find(hub)!=Label[ID2].end()){
			dis2=Label[ID2][hub];
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

void Graph::CorrectnessCheck(int runtimes){
	srand (time(NULL));
	int s, t, d1, d2, d3;
    runtimes = 1;
    cout<<"Correctness check ("<<runtimes<<" rounds) ..."<<endl;
	for(int i=0;i<runtimes;i++){
//        if(i%100==0) cout<<i<<endl;
		s=rand()%nodenum;
		t=rand()%nodenum;
//        if(i==0){
//            s=197277;t=165625;//parti-parti
//        }else if(i==1){
//            s=159503;t=188356;//parti-parti
//        }
        s=164495;t=207358;//parti-parti
//        s=142488;t=143850;//core-core
//        s=208312;t=210695;//same-parti
//        s=327243,t=752625;//parti-core
//        s=92495,t=401841;//core-parti
		d1=Dijkstra(s,t,Neighbor);
		d2=Query(s,t);
//        d1=DijkstraCore(s,t);
//        d2=QueryCore(s,t);
//        cout<<s<<"("<<CoreTag[s]<<") "<<t<<"("<<CoreTag[t]<<") "<<d2<<" "<<d1<<endl;
		if(d1!=d2){
			cout<<"InCorrect! "<<i<<": "<<s<<" "<<t<<" "<<d2<<" "<<d1;
            cout<<" ; CoreTag: "<<CoreTag[s]<<" "<<CoreTag[t]<<endl;
            QueryDebug(s,t);
		}
	}
}

void Graph::CorrectnessCheckCore(){
	srand (time(NULL));
	int s, t, d1, d2, d3;
    vector<int> coreVertex;
    for(int i=0;i<nodenum;++i){
        if(CoreTag[i] == -1){
            coreVertex.emplace_back(i);
        }
    }
    int corenum=coreVertex.size();
	for(int i=0;i<100;i++){
		s=coreVertex[rand()%corenum];
		t=coreVertex[rand()%corenum];
		if(CoreTag[s]==-1 && CoreTag[t]==-1){//for core vertex
			d1=QueryCore(s,t);
            d2=DijkstraCore(s,t);


			if(d1!=d2){
				cout<<"InCorrect!"<<s<<" "<<t<<" "<<d1<<" "<<d2<<endl;
//				DijkstraPath(s,t);
//				DijkstraCorePath(s,t);
			}
		}else
			i--;
	}
}

void Graph::CorrectnessCheckCoreAll(){
    int s,t,d1,d2;
    cout<<"Correctness checking..."<<endl;
    for(auto it1=CoreVertex.begin();it1!=CoreVertex.end();++it1){
        s=*it1;
        for(auto it2=it1;it2!=CoreVertex.end();++it2){
            t=*it2;
            if(s==t)
                continue;
            d1=DijkstraCore(s,t);
            d2=QueryCore(s,t);
            if(d1!=d2){
                cout<<"InCorrect! "<<s<<" "<<t<<" "<<d2<<" "<<d1<<endl;
                QueryCoreDebug(s,t);
            }
        }
    }
//    cout<<"Done."<<endl;
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


void Graph::ReadUpdate2(string filename,vector<pair<pair<int,int>,pair<int,int>>>& TestData){
    TestData.clear();

    int num, ID1, ID2, oldW,newW;
    ifstream IF(filename);
    if(!IF){
        cout<<"Cannot open file "<<filename<<endl;
        exit(1);
    }
    string line;
    getline(IF,line);

    while(!line.empty()){
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
        ID1=stoi(vs[0]); ID2=stoi(vs[1]); oldW=stoi(vs[2]); newW=stoi(vs[3]);
        TestData.push_back(make_pair(make_pair(ID1, ID2), make_pair(oldW,newW)));
        if(IF.eof())
            break;
        getline(IF,line);
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
//function of reading core label to disk
void Graph::ReadLabels(string file){
    unordered_map<int,int> map0; map0.clear();
    Label.assign(nodenum, map0);


//    PruningPoint.clear();
//    unordered_map<int,set<int>> map1; map1.clear();
//    PruningPoint.assign(nodenum,map1);

    PruningPointNew.clear();
    unordered_map<int,vector<int>> map2; map2.clear();
    PruningPointNew.assign(nodenum,map2);

    PruningPointList.clear();
    unordered_map<int,list<int>> map3; map3.clear();
    PruningPointList.assign(nodenum,map3);

    cout<<"Reading labels..."<<endl;
    /// read label
    ifstream IF(file+".label");
    if(!IF){
        cout<<"Cannot open file "<<file+".label"<<endl;
        exit(1);
    }
    int nodeNum;
    int tempE=0;
    int ID1,ID2,weight;
    string line;

    getline(IF,line);
    vector<string> vs;
    boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
    nodeNum = stoi(vs[0]);
    assert(nodenum == nodeNum);
    getline(IF,line);

    while(!line.empty()){
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
        ID1=stoi(vs[0]); ID2=stoi(vs[1]); weight=stoi(vs[2]);

        Label[ID1].insert({ID2,weight});

        if(IF.eof())
            break;
        getline(IF,line);
    }

    IF.close();
    /// read pruning point
    ifstream IF2(file+".prune");
    if(!IF2){
        cout<<"Cannot open file "<<file+".prune"<<endl;
        exit(1);
    }


//    set<int> vSet;
    getline(IF2,line);
    while(!line.empty()){
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
        ID1=stoi(vs[0]); ID2=stoi(vs[1]);
        set<int> vSet; vSet.clear();
        for(int i=2;i<vs.size();++i){
            vSet.insert(stoi(vs[i]));
            PruningPointNew[ID1][ID2].push_back(stoi(vs[i]));//((v,c),w)
            PruningPointList[ID1][ID2].push_back(stoi(vs[i]));
        }

//        PruningPoint[ID1].insert({ID2,vSet});


        if(IF2.eof())
            break;
        getline(IF2,line);
    }

    IF2.close();
}
//function of writing core label to disk
void Graph::WriteLabels(string file){
    ofstream OF(file+".label");
    OF<<nodenum<<endl;
    int tempE=0;
    int ID2,wei;
    for(int ID1=0;ID1<nodenum;ID1++){
        if(Label[ID1].size()>1){
            for(auto it=Label[ID1].begin();it!=Label[ID1].end();++it){
                ID2=it->first;
                wei=it->second;
                OF<<ID1<<" "<<ID2<<" "<<wei<<endl;
                tempE++;
            }
        }
    }
    OF.close();
    /// Write pruning point
    ofstream OF2(file+".prune");
    for(int ID1=0;ID1<nodenum;ID1++){
        if(!PruningPointList[ID1].empty()){
            for(auto it=PruningPointList[ID1].begin();it!=PruningPointList[ID2].end();++it){
                ID2=it->first;
                if(!it->second.empty()){
                    OF2<<ID1<<" "<<ID2;
                    for(auto it2=it->second.begin();it2!=it->second.end();++it2){
                        wei = *it2;
                        OF2<<" "<<wei;
                    }
                    OF2<<endl;
                }
            }
        }

    }
    OF2.close();
}

void Graph::WriteCoreGraph(string graphfile){
    cout<<"Writing core graph..."<<endl;
	ofstream OF(graphfile);
	OF<<nodenum<<endl;
    int tempE=0;
	for(int ID1=0;ID1<nodenum;ID1++){
		for(int j=0;j<AdjaCore[ID1].size();j++){
			int ID2=AdjaCore[ID1][j].first;
			int wei=AdjaCore[ID1][j].second;
			OF<<ID1<<" "<<ID2<<" "<<wei<<endl;
            tempE++;
		}
	}
    OF.close();
    ofstream OF2(graphfile+".order");
    OF2<<nodenum<<" "<<tempE<<endl;
    for(int ID1=0;ID1<nodenum;ID1++){
            OF2<<ID1<<" "<<NodeOrder[ID1]<<endl;
    }
    OF2.close();
    cout<<"Done."<<endl;
}

void Graph::ReadCoreGraph(string filename){
	ifstream inGraph(filename);
	if(!inGraph){
		cout<<"Cannot open Map "<<filename<<endl;
        exit(1);
	}

	string line;

	//graph g initialize


    getline(inGraph,line);
    vector<string> vs;
    boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
    nodenum=stoi(vs[0]); edgenum=0;

    vector<pair<int,int>> vecp; vecp.clear();
    AdjaCore.assign(nodenum, vecp);
    AdjaCoreMap.assign(nodenum,map<int,int>());
    getline(inGraph,line);

	int ID1,ID2, weight;
	while(!line.empty()){
		vector<string> vs;
		boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
		ID1=stoi(vs[0]); ID2=stoi(vs[1]); weight=stoi(vs[2]);

        CoreVertex.insert(ID1);
		AdjaCore[ID1].push_back(make_pair(ID2,weight));
        AdjaCoreMap[ID1].insert({ID2,weight});
        edgenum++;
		if(inGraph.eof())
            break;
		getline(inGraph,line);
	}

	cout<<"Finish Graph Reading! Nnum "<<CoreVertex.size()<<" Enum "<<edgenum<<endl;

    for(ID1=0;ID1<nodenum;++ID1){
        for(auto it=AdjaCore[ID1].begin();it!=AdjaCore[ID1].end();++it){
            ID2=it->first; weight=it->second;
            if(AdjaCoreMap[ID2].find(ID1) == AdjaCoreMap[ID2].end()){//if not found
                cout<<"Wrong! "<<ID1<<" "<<ID2<<endl;
            }else if(weight != AdjaCoreMap[ID2][ID1])
            {
                cout<<"Wrong! "<<ID1<<" "<<ID2<<weight<<AdjaCoreMap[ID2][ID1]<<endl;
            }

        }
    }

    CoreTag.assign(nodenum,0);
    for(auto it=CoreVertex.begin();it!=CoreVertex.end();++it){
        CoreTag[*it] = -1;
    }


    ifstream inGraph2(filename+".order");
    if(!inGraph2){
        cout<<"Cannot open Map "<<filename+".order"<<endl;
    }
    getline(inGraph2,line);
    vs.clear();
    boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
    int nodeNum=stoi(vs[0]);
    int edgeNum=stoi(vs[1]);
    assert(nodeNum == nodenum);
    assert(edgeNum == edgenum);
    getline(inGraph2,line);
    NodeOrder.assign(nodenum,-1);
    vNodeOrder.assign(nodenum,-1);

    int order;
    while(!line.empty()){
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
        ID1=stoi(vs[0]); order=stoi(vs[1]);
        NodeOrder[ID1] = order;
        vNodeOrder[order] = ID1;
        if(inGraph2.eof())
            break;
        getline(inGraph2,line);
    }
	cout<<"Node finish ordering!"<<endl;
}

void Graph::WriteOrder(string filename){
    //Write order file to disk
    ofstream OF(filename);
    if(!OF){
        cout<<"Cannot open Map "<<filename<<endl;
        exit(1);
    }
    OF<<nodenum<<endl;
    int tempE=0;
    int ID2,wei;
    for(int ID1=0;ID1<nodenum;ID1++){

         OF<<ID1<<" "<<NodeOrder[ID1]<<endl;

    }
    OF.close();
}
//function of checking the connectivity, set_A: the vertex set
vector<int> Graph::DFS_CC(vector<map<int,int>> & Edges, set<int> set_A, set<int> & set_LCC, int nodenum) {
    /// DFS for connected component
    stack<int> stack_A;
//    set<int> set_A;//nodes waiting for visited
    unordered_set<int> set_B;//nodes visited for current component
    set_B.clear();
    int item_id,temp_id;
    vector<bool> flag_visited(nodenum,false);
    bool flag_finish = false;
    int temp_num = 0;
    int component_i = 0;
    pair<unordered_set<int>,int> LCC;
    vector<int> CCs;//the vertex size of each connected component

//    for(int i=0;i<nodenum;++i){
//        set_A.insert(i);
//    }
    int seed = *set_A.begin();
    stack_A.push(seed);
    set_A.erase(seed);
    set_B.insert(seed);
    flag_visited[seed] = true;
    //Get the connected components by DFS
    while(!set_A.empty()) {//if not finish
        temp_num = 0;
        while (!stack_A.empty()) {
            item_id = stack_A.top();
            stack_A.pop();
            for (auto it = Edges[item_id].begin(); it != Edges[item_id].end(); ++it) {
                temp_id = it->first;
                temp_num += 1;
                if (!flag_visited[temp_id]) {//if not visited
                    stack_A.push(temp_id);
                    set_A.erase(temp_id);
                    set_B.insert(temp_id);
                    flag_visited[temp_id] = true;
                }
            }
        }
        if (set_B.size() > LCC.first.size()) {
            LCC.first.clear();
            LCC.first = set_B;
            LCC.second = temp_num;// /2
        }
        assert(!set_B.empty());
        CCs.push_back(set_B.size());
//        if(!set_B.empty() && set_B.size() < mcc.first.size()){
//            cout<<"Smaller connected component with vertex size "<<set_B.size()<<": ";
//            for(auto it=set_B.begin();it!=set_B.end();++it){
//                cout<<*it<<" ";
//            }
//            cout<<"; degree: ";
//            for(auto it=set_B.begin();it!=set_B.end();++it){
//                cout<<Edges[*it].size()<<" ";
//            }
//            cout<<endl;
//        }
        ++component_i;
        set_B.clear();
        if (!set_A.empty()) {
            stack_A.push(*set_A.begin());
            set_B.insert(*set_A.begin());
            flag_visited[*set_A.begin()] = true;
            set_A.erase(*set_A.begin());
        } else {
            break;
        }
    }
    if(component_i==1){
        cout<<"This graph has only one connected component. ";
        cout<<"Nodes size of graph: "<< LCC.first.size() << " ; ";
        cout<<"Edges size of graph: "<< LCC.second << endl;
    }else{
        cout<<"!!! This graph has "<< component_i <<" connected component!"<<endl;
        cout<<"Nodes size of the largest connected component is: "<<LCC.first.size()<<endl;
        cout<<"Edges size of the largest connected component is: "<<LCC.second<<endl;
    }
    for(auto it=LCC.first.begin();it!=LCC.first.end();++it){
        set_LCC.insert(*it);
    }
    std::sort(CCs.begin(), CCs.end());
    return CCs;
//    return component_i;
}

struct PQEdgePairCompareLess
{
    bool operator () (const pair<uint,pair<int,int>>& a, const pair<uint,pair<int,int>>& b) const
    {
//        return (a.first > b.first);
        return (a.first > b.first) || (b.first <= a.first && (a.second.first > b.second.first));
    }

    pair<uint,pair<int,int>> min_value() const
    {   pair<uint,pair<int,int>> a;
        a.first=std::numeric_limits<uint>::max();
//        a.second=pair<int,int>(std::numeric_limits<int>::max(),std::numeric_limits<int>::max());
        a.second=make_pair(std::numeric_limits<int>::max(),std::numeric_limits<int>::max());
        return a; }
};

//Function of counting the minimal spanning tree
int Graph::MinSpanTree(vector<vector<pair<int,int>>> & Neighbors)//Lazy Prim's algorithm, complexity O(|E|log(|E|)). Eager Prim's algorithm use an indexed priority queue which can efficiently update and poll key-value pairs, further reduce the complexity to O(|E|log|V|).
{
//    benchmark::heap<2, int, int> pqueue(nodenum);
    priority_queue<pair<int,pair<int,int>>,vector<pair<int,pair<int,int>>>,PQEdgePairCompareLess> pqueue;
    unsigned long long mstCost=0;
    pair<int,pair<int,int>> item_;
    int item_id,item_dis;
    int temp_id,temp_dis;
    long long int temp_w;
    int nodeNum = Neighbors.size();
    vector<bool> closed(nodeNum, false); //flag vector of whether closed
//    vector<map<int,int>> MSTGraph(nodeNum,map<int,int>());
    vector<vector<pair<int,int>>> MSTGraph(nodeNum,vector<pair<int,int>>());

    uint mst_edgenum = nodeNum-1;//the number of undirected edges
    uint edgeCount = 0;
    double ave_w=0;

    //Initiation of start node
    item_id = 0;
    closed[item_id] = true;
    for (auto it = Neighbors[item_id].begin(); it != Neighbors[item_id].end(); ++it) {
        temp_id = it->first;
        if(closed[temp_id])
            continue;
        temp_w = it->second;
        pqueue.push(make_pair(temp_w,make_pair(item_id,temp_id)));
    }

    //Iteration
    while (!pqueue.empty() && edgeCount!=mst_edgenum) {//for every node in pqueue
        item_ = pqueue.top();// top min item
        pqueue.pop();
        item_id = item_.second.second;
        if(closed[item_id])
            continue;
        //push the minimal edges into mst
//        MSTGraph[item_.second.first].insert({item_.second.second,item_.first});
//        MSTGraph[item_.second.second].insert({item_.second.first,item_.first});
        MSTGraph[item_.second.first].emplace_back(item_.second.second,item_.first);
        MSTGraph[item_.second.second].emplace_back(item_.second.first,item_.first);

        ++edgeCount;
        mstCost += temp_w;
        closed[item_id] = true;
        for (auto it = Neighbor[item_id].begin(); it != Neighbor[item_id].end(); ++it) {
            temp_id = it->first;
            if(closed[temp_id])
                continue;
            temp_w = it->second;
            pqueue.push(make_pair(temp_w,make_pair(item_id,temp_id)));
        }
    }
    cout<<"Size of spanning tree: "<<edgeCount<<endl;
//        cout<<"MST cost is "<<mstCost<<endl;

    if(edgeCount!=mst_edgenum){
        cout<<"Minimum spanning tree does not exist! edgeCount: "<<edgeCount<<", mst_edgenum: "<<mst_edgenum<<endl;
        return -1;
    } else{
//        int s,t,d1,d2;
//        int runtimes = 100;
//        for(int i=0;i<runtimes;i++){
////        if(i%100==0) cout<<i<<endl;
//            s=rand()%nodeNum;
//            t=rand()%nodeNum;
//            d1=Dijkstra(s,t,MSTGraph);
//            d2=Dijkstra(s,t,Neighbors);
//            if(d1!=d2){
//                cout<<"InCorrect! "<<i<<": "<<s<<" "<<t<<" "<<d1<<" "<<d2<<endl;
//            }
//        }
        return 0;
    }
}

void Graph::CoreGraphDebug(string graphfile) {
    std::chrono::high_resolution_clock::time_point t10;
    std::chrono::high_resolution_clock::time_point t11;
    std::chrono::duration<double> time_span1;
    double runT1=0;
    int s,t,d1,d2;
    int ID1,ID2,oldW,newW;

    ReadCoreGraph(graphfile);
//    set<int> setB;
//    DFS_CC(AdjaCoreMap,CoreVertex,setB,nodenum);

    vector<int> coreVertex;
    for(auto it=CoreVertex.begin();it!=CoreVertex.end();++it){
        coreVertex.push_back(*it);
    }
    int num=coreVertex.size();


//    Construct_core(2);//Batch PSL
//    Construct_core(1);//PSL
    Construct_core(0);//PLL
//    WriteLabels(graphfile+"2");


    if(dataset == "Test"){
        Construct_core(0);// PLL index construction
//    Construct_core(1);// PSL index construction
//    WriteLabels(graphfile+".labelPLL");
    }else{
//        ReadLabels(graphfile);
    }


    CorrectnessCheckCore();


    //update
    string updateFile2 = graphfile+".update";
    vector<pair<pair<int,int>,pair<int,int>>> testdata2;
    ReadUpdate2(updateFile2, testdata2);
    cout<<"Batch number: "<<testdata2.size()<<endl;
//    for(int u=0;u<testdata2.size();++u){
    for(int u=0;u<1;++u){
        ID1=testdata2[u].first.first;
        ID2=testdata2[u].first.second;
        oldW=testdata2[u].second.first;
        newW=testdata2[u].second.second;
        cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;

        cout<<"Before update: "<<QueryCore(164496,207450)<<" "<<DijkstraCore(164496,207450)<<endl;
        cout<<QueryCore(207450,211814)<<" "<<DijkstraCore(207450,211814)<<endl;
        if(Label[207450].find(211814) != Label[207450].end()){//if found
            cout<<Label[207450][211814]<<" "<<QueryCore(207450,211814)<<" "<<DijkstraCore(207450,211814)<<endl;
        }
//
//        DijkstraPathCore(212434,99110);
//        QueryCoreDebug(164496,207450);
        t10=std::chrono::high_resolution_clock::now();
//        IncreasePSL(ID1,ID2,oldW,newW,AdjaCore,Label,PruningPointNew,NoSupportedPair);//original vector version with NoSuportedPair
//        IncreasePSL2(ID1,ID2,oldW,newW,AdjaCore,Label,PruningPointList);//list version with NoSuportedPair
        IncreasePSL(ID1,ID2,oldW,newW,AdjaCore,Label,PruningPointList);//list version without NoSuportedPair
//        IncreasePLL(ID1,ID2,oldW,newW,AdjaCore,Label,PruningPointNew,NoSupportedPair);
//        IncreasePLL(ID1,ID2,oldW,newW,AdjaCore,Label,PruningPointList);
        t11=std::chrono::high_resolution_clock::now();

        time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
        runT1 += time_span1.count();
//                g.QueryCoreDebug(142488,143850);
//        g.CorrectnessCheck(200);
        cout<<"After update: "<<QueryCore(164496,207450)<<" "<<DijkstraCore(164496,207450)<<endl;
        if(Label[207450].find(211814) != Label[207450].end()){//if found
            cout<<Label[207450][211814]<<" "<<QueryCore(207450,211814)<<" "<<DijkstraCore(207450,211814)<<endl;
        }
        if(u==0){// || u==279
//            WriteCoreGraph(graphfile+"2");
//            WriteLabels(graphfile+"2");
//            cout<<"Done."<<endl;
//            QueryCoreDebug(142488,143850);
        }

        if(dataset == "Test"){
            CorrectnessCheckCoreAll();
        }
        else{
//            CorrectnessCheckCoreAll();
            cout<<"Correctness checking..."<<endl;
            for(int i=0;i<1;i++){
                s=coreVertex[rand()%num];
                t=coreVertex[rand()%num];
                s=164496;t=207450;
//                s=150791;t=196252;//core-core
//                s=212434;t=99110;//core-core, for core3
//                s=212434;t=99509;//core-core, for core3
                assert(CoreVertex.find(s)!=CoreVertex.end());
                assert(CoreVertex.find(t)!=CoreVertex.end());
                d1=DijkstraCore(s,t);
                d2=QueryCore(s,t);

//            cout<<s<<"("<<g.CoreTag[s]<<") "<<t<<"("<<g.CoreTag[t]<<") "<<d2<<" "<<d1<<" "<<d3<<endl;
                if(d1!=d2){
                    int d3=DijkstraCore(t,s);
                    cout<<"InCorrect! "<<i<<": "<<s<<" "<<t<<" : "<<d2<<" "<<d1<<" "<<d3<<endl;
                    QueryDebug(s,t);
                }
            }
        }

    }
    exit(0);
}
