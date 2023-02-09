/*
 * H2H.cpp
 *
 *  Created on: 14 Oct 2022
 *      Author: zhangmengxuan
 */
#include "head3.h"

vector<int> NodeOrderss;
struct OrderCompp{//prior to reture the vertex with smaller order
	int x;
	OrderCompp(int _x){
		x=_x;
	}
	bool operator< (const OrderCompp& d) const{
		if(x==d.x){//avoid the redundant
			return false;
		}else{
			if(x!=d.x)
				return NodeOrderss[x]<NodeOrderss[d.x];
		}
	}
};

//void Graph::H2Hindex(){
//	vector<pair<int,pair<int,int>>> vect;
//	NeighborCon.assign(nodenum,vect);
//
//	//CHindex();//single threads
//	CHindexMT();//multiple threads
//	makeTree();
//	makeIndex();
//}
/*
void Graph::CHindex(){
	//initialize E
	map<int,pair<int,int>> m;
	E.assign(nodenum,m);
	for(int i=0;i<Neighbors.size();i++){
		for(int j=0;j<Neighbors[i].size();j++)
			E[i].insert(make_pair(Neighbors[i][j].first,make_pair(Neighbors[i][j].second,1)));
	}

	vector<bool> exist; exist.assign(nodenum,true);
	vector<bool> change; change.assign(nodenum,false);

	//cout<<"Begin to contract"<<endl;
	SCconNodes.clear();
	for(int nodeorder=0;nodeorder<nodenum;nodeorder++){
		int x=vNodeOrder[nodeorder];
		//if(x!=-1){//to identify and exclude the isolated vertices
		if(Neighbors[x].size()!=0){//to identify and exclude the isolated vertices
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
}

void Graph::CHindexMT(){
	map<int, vector<int>> mi;
	SCconNodesMT.assign(nodenum, mi);

	//initialize E
	map<int,pair<int,int>> m;
	E.assign(nodenum,m);
	for(int i=0;i<Neighbors.size();i++){
		for(int j=0;j<Neighbors[i].size();j++)
			E[i].insert(make_pair(Neighbors[i][j].first,make_pair(Neighbors[i][j].second,1)));
	}

	vector<bool> exist; exist.assign(nodenum,true);
	//vector<bool> change; change.assign(nodenum,false);

	//cout<<"Begin to contract"<<endl;
	for(int nodeorder=0;nodeorder<nodenum;nodeorder++){
		int x=vNodeOrder[nodeorder];
		//if(x!=-1){//to identify and exclude the isolated vertices
		if(Neighbors[x].size()!=0){//to identify and exclude the isolated vertices
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
}*/

//void Graph::NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x){
//	sm->wait();
//	int ID1, w1;
//	int ID2, w2;
//	for(int k=p.first;k<p.second;k++){
//		ID1=Neighvec[k].first;
//		w1=Neighvec[k].second.first;
//		for(int h=0;h<Neighvec.size();h++){
//			ID2=Neighvec[h].first;
//			w2=Neighvec[h].second.first;
//			insertEMTorder(ID1, ID2, w1+w2);
//			if(ID1<ID2)
//				SCconNodesMT[ID1][ID2].push_back(x);
//		}
//	}
//	sm->notify();
//}

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
//function of inserting edge (u,v)
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
//function of erasing edge (u,v), i.e., erase u from v's adjacency list.
void Graph::deleteEorder(int u,int v){
//	if(E[u].find(v)!=E[u].end()){
//		E[u].erase(E[u].find(v));
//	}

	if(E[v].find(u)!=E[v].end()){
		E[v].erase(E[v].find(u));
	}
}

/*void Graph::makeTree(){
	vector<int> vecemp; //vecemp.clear();
	VidtoTNid.assign(nodenum,vecemp);

	//rank.assign(nodenum,0);
	rank.assign(nodenum,-1);
	//Tree.clear();
	int len=vNodeOrder.size()-1;
	heightMax=0;

	Node rootn;
	int x=vNodeOrder[len];
	//cout<<"len "<<len<<" , ID "<<x<<endl;
	//while(x==-1){//to skip those vertices whose ID is -1
	while(Neighbors[x].size()==0){//to skip those vertices whose ID is -1
		len--;
		x=vNodeOrder[len];
		//cout<<"len "<<len<<" , ID "<<x<<endl;
	}
	rootn.vert=NeighborCon[x];
	rootn.uniqueVertex=x;
	rootn.pa=-1;
	rootn.height=1;
	rank[x]=0;
	Tree.push_back(rootn);
	len--;

	int nn;
	for(;len>=0;len--){
		int x=vNodeOrder[len];
		if(Neighbors[x].size()!=0){
			Node nod;
			nod.vert=NeighborCon[x];
			nod.uniqueVertex=x;
			int pa=match(x,NeighborCon[x]);
			Tree[pa].ch.push_back(Tree.size());
			nod.pa=pa;
			nod.height=Tree[pa].height+1;

			nod.hdepth=Tree[pa].height+1;
			for(int i=0;i<NeighborCon[x].size();i++){
				nn=NeighborCon[x][i].first;
				VidtoTNid[nn].push_back(Tree.size());
				if(Tree[rank[nn]].hdepth<Tree[pa].height+1)
					Tree[rank[nn]].hdepth=Tree[pa].height+1;
			}
			if(nod.height>heightMax) heightMax=nod.height;
			rank[x]=Tree.size();
			Tree.push_back(nod);
		}
		//cout<<"len "<<len<<" , ID "<<x<<endl;
	}
}*/

int Graph::match(int x,vector<pair<int,pair<int,int>>> &vert){
	int nearest=vert[0].first;
	for(int i=1;i<vert.size();i++){
		if(rank[vert[i].first]>rank[nearest])//get the vertex with highest rank
			nearest=vert[i].first;
	}
	int p=rank[nearest];
	return p;
}


void Graph::makeRMQCore(){
    //EulerSeq.clear();
    toRMQ.assign(nodenum,0);
    //RMQIndex.clear();
    makeRMQDFSCore(0, 1);
    RMQIndex.push_back(EulerSeq);

    int m = EulerSeq.size();
    for (int i = 2, k = 1; i < m; i = i * 2, k++){
        vector<int> tmp;
        //tmp.clear();
        tmp.assign(m,0);
        for (int j = 0; j < m - i; j++){
            int x = RMQIndex[k - 1][j], y = RMQIndex[k - 1][j + i / 2];
            if (Tree[x].height < Tree[y].height)
                tmp[j] = x;
            else tmp[j] = y;
        }
        RMQIndex.push_back(tmp);
    }
}

void Graph::makeRMQDFSCore(int p, int height){
    toRMQ[p] = EulerSeq.size();//toRMQ[p] records the first occurrence of node p in the Euler tour, i.e., R[p]
    EulerSeq.push_back(p);//EulerSeq is the Euler tour, i.e., E[1,...,2n-1]
    for (int i = 0; i < Tree[p].ch.size(); i++){
        makeRMQDFSCore(Tree[p].ch[i], height + 1);
        EulerSeq.push_back(p);
    }
}

void Graph::makeRMQDFS(int p, int height){
	toRMQ[p] = EulerSeq.size();//toRMQ[p] records the first occurrence of node p in the Euler tour, i.e., R[p]
	EulerSeq.push_back(p);//EulerSeq is the Euler tour, i.e., E[1,...,2n-1]
	for (int i = 0; i < Tree[p].ch.size(); i++){
		makeRMQDFS(Tree[p].ch[i], height + 1);
		EulerSeq.push_back(p);
	}
}
//function of computing the H2H label of peripheries: version 2, without the correct info from core
void Graph::makeIndexDFSCore(int p, vector<int>& list, vector<int> & interface, map<int,unordered_map<int,int>> & disInfs){
//    if(p == rank[210695]){
//        cout<<210695<<" "<<CoreTag[210695]<<" "<<Tree[Tree[p].treeroot].uniqueVertex<<" "<<Tree[p].treeroot<<endl;
//    }
    //initialize
    int NeiNum=Tree[p].vert.size();
    Tree[p].pos.assign(NeiNum+1,-1);
    Tree[p].dis.assign(list.size(),INF);
    Tree[p].cnt.assign(list.size(),0);
    Tree[p].FN.assign(list.size(),true);
//    Tree[p].FNInf.assign(interface.size(),true);
    for(int i=0;i<interface.size();++i){
        Tree[p].FNInf.insert({interface[i],true});
    }
    Tree[p].vAncestor.assign(list.begin(), list.end());
    Tree[p].vAncestor.push_back(Tree[p].uniqueVertex);
    Tree[p].dis.push_back(0);
    //for interface
    int InfNum = interface.size();
//    Tree[p].disInf.assign(InfNum,INF);
    int ID1=Tree[p].uniqueVertex;
    int ID2;

    /// interface
    for(int j=0;j<interface.size();j++){//for each interface vertex
        ID2=interface[j];
//        if(p == rank[210695] && ID2 == 208312)
//            cout<<ID2<<endl;
        Tree[p].disInf.insert({ID2,INF});
        for(int i=0;i<NeiNum;i++){//for each neighbor
            int x=Tree[p].vert[i].first;
            int disvb=Tree[p].vert[i].second.first;//
            if(x == interface[j]){//if it is current interface vertex
                if(Tree[p].disInf[ID2] > disvb){
                    Tree[p].disInf[ID2]=disvb;//
                    Tree[p].FNInf[ID2]=true;//obtain from vert
                }
                continue;
            }
            int z;
            if(rank[x] > 0 ){//if x is not boundary vertex, i.e., if it is ancestor
                z = Tree[rank[x]].disInf[ID2];
                if(Tree[p].disInf[ID2]>z+disvb){
                    Tree[p].disInf[ID2]=z+disvb;
//                    if(!Tree[rank[x]].FNInf[j]){//if the ancestor's interface distance is not from vert
//                        Tree[p].FNInf[ID2]=false;
//                    }
                    Tree[p].FNInf[ID2]=false;
                }
            }else{//if it is other interface vertex
                if(x<=ID2){
//                    assert(disInfs[x].find(ID2) != disInfs[x].end());
                    z = disInfs[x][ID2];
                }else{
//                    assert(disInfs[ID2].find(x) != disInfs[ID2].end());
                    z = disInfs[ID2][x];
                }
                assert(z>0);
                if(Tree[p].disInf[ID2]>z+disvb){
                    Tree[p].disInf[ID2]=z+disvb;
                    Tree[p].FNInf[ID2]=false;
                }
            }
        }
//        int d1= Tree[p].disInf[j];
//        int d2= Dijkstra(ID1,ID2,Neighbor);
//        if(d1!=d2){
//            cout<<"Incorrect! "<<ID1<<" "<<ID2<<": "<<d1<<" "<<d2<<endl;
//        }
    }

    /// ancestor
    //map<int,Nei> Nmap; Nmap.clear();//shortcut infor ordered by the pos ID
    for(int i=0;i<NeiNum;i++){
        for(int j=1;j<list.size();j++){
            if(Tree[p].vert[i].first==list[j]){
                Tree[p].pos[i]=j;
                Tree[p].dis[j]=Tree[p].vert[i].second.first;//
                Tree[p].cnt[j]=1;
                break;
            }
        }
    }
    Tree[p].pos[NeiNum]=list.size();

    bool flag= false;
    for(int j=1;j<list.size();j++){
        int y=list[j];//the jth ancestor is y
        int z;//the distance from x to y
//        if(p == rank[210695] && y == 208312)
//            cout<<y<<endl;
        for(int i=0;i<NeiNum;i++){
            int x=Tree[p].vert[i].first;
            int disvb=Tree[p].vert[i].second.first;//shortcut
            int k=Tree[p].pos[i];//the kth ancestor is x

            if(k==-1){//if the neighbor is interface vertex
                flag = false;
                for(int id=0;id<interface.size();++id){
                    if(interface[id] == x){
                        flag = true;
                        z=Tree[rank[y]].disInf[interface[id]];//from the ancestor to interface
                        if(Tree[p].dis[j]>z+disvb){
                            Tree[p].dis[j]=z+disvb;
                            Tree[p].FN[j]=false;
                            Tree[p].cnt[j]=1;
                        }else if(Tree[p].dis[j]==z+disvb){
                            Tree[p].cnt[j]+=1;//record how many path has the shortest distance
                        }
                        break;
                    }
                }
                if(!flag){
                    cout<<"There is non-partition vertex in the adjacency list! "<<x<<endl;
                    exit(1);
                }
            }
            else{//if the neighbor is ancestor
                if(k!=j){
                    if(k<j)
                        z=Tree[rank[y]].dis[k];//get the distance to the lower-order vertex
                    else if(k>j)
                        z=Tree[rank[x]].dis[j];

                    if(Tree[p].dis[j]>z+disvb){
                        Tree[p].dis[j]=z+disvb;
                        Tree[p].FN[j]=false;
                        Tree[p].cnt[j]=1;
                    }else if(Tree[p].dis[j]==z+disvb){
                        Tree[p].cnt[j]+=1;
                    }
                }
            }

        }
//        int d1= Tree[p].dis[j];
//        int d2= Dijkstra(ID1,y,Neighbor);
//        if(d1!=d2){
//            cout<<"Incorrect! "<<ID1<<" "<<y<<": "<<d1<<" "<<d2<<endl;
//        }
    }
    assert(Tree[p].disInf.size()==InfNum);
    //nested loop
    list.push_back(Tree[p].uniqueVertex);
    for(int i=0;i<Tree[p].ch.size();i++){
        makeIndexDFSCore(Tree[p].ch[i],list,interface,disInfs);
    }
    list.pop_back();
}

//function of computing the H2H label of peripheries: version 1, leveraging the correct info from core
/*void Graph::makeIndexDFSCore(int p, vector<int>& list, vector<int> & interface, map<int,unordered_map<int,int>> & disInfs){
    //initialize
    int NeiNum=Tree[p].vert.size();
    Tree[p].pos.assign(NeiNum+1,-1);
    Tree[p].dis.assign(list.size(),INF);
    Tree[p].cnt.assign(list.size(),0);
    Tree[p].FN.assign(list.size(),true);
    Tree[p].vAncestor.assign(list.begin(), list.end());
    Tree[p].vAncestor.push_back(Tree[p].uniqueVertex);
    Tree[p].dis.push_back(0);
    //for interface
    int InfNum = interface.size();
//    Tree[p].disInf.assign(InfNum,INF);
    int ID1=Tree[p].uniqueVertex;
    int ID2;

    /// interface
    for(int j=0;j<interface.size();j++){
        ID2=interface[j];
        Tree[p].disInf.insert({ID2,INF});
        for(int i=0;i<NeiNum;i++){
            int x=Tree[p].vert[i].first;
            int disvb=Tree[p].vert[i].second.first;//
            if(x == interface[j]){
                if(Tree[p].disInf[ID2] > disvb){
                    Tree[p].disInf[ID2]=disvb;//
                }
                continue;
            }
            int z;
            if(rank[x] > 0 ){//if not boundary vertex, i.e., if it is ancestor
                z = Tree[rank[x]].disInf[ID2];
                if(Tree[p].disInf[ID2]>z+disvb){
                    Tree[p].disInf[ID2]=z+disvb;
                }
            }else{//if it is interface vertex
                if(x<=ID2){
//                    assert(disInfs[x].find(ID2) != disInfs[x].end());
                    z = disInfs[x][ID2];
                }else{
//                    assert(disInfs[ID2].find(x) != disInfs[ID2].end());
                    z = disInfs[ID2][x];
                }
                assert(z>0);
                if(Tree[p].disInf[ID2]>z+disvb){
                    Tree[p].disInf[ID2]=z+disvb;
                }
            }
        }
//        int d1= Tree[p].disInf[j];
//        int d2= Dijkstra(ID1,ID2,Neighbor);
//        if(d1!=d2){
//            cout<<"Incorrect! "<<ID1<<" "<<ID2<<": "<<d1<<" "<<d2<<endl;
//        }
    }

    /// ancestor
    //map<int,Nei> Nmap; Nmap.clear();//shortcut infor ordered by the pos ID
    for(int i=0;i<NeiNum;i++){
        for(int j=1;j<list.size();j++){
            if(Tree[p].vert[i].first==list[j]){
                Tree[p].pos[i]=j;
                Tree[p].dis[j]=Tree[p].vert[i].second.first;//
                Tree[p].cnt[j]=1;
                break;
            }
        }
    }
    Tree[p].pos[NeiNum]=list.size();

    bool flag= false;
    for(int j=1;j<list.size();j++){
        int y=list[j];//the jth ancestor is y
        int z;//the distance from x to y
        for(int i=0;i<NeiNum;i++){
            int x=Tree[p].vert[i].first;
            int disvb=Tree[p].vert[i].second.first;//shortcut
            int k=Tree[p].pos[i];//the kth ancestor is x
            if(k==-1){//if the neighbor is interface vertex
                flag = false;
                for(int id=0;id<interface.size();++id){
                    if(interface[id] == x){
                        flag = true;
                        z=Tree[rank[y]].disInf[interface[id]];//from the ancestor to interface
                        if(Tree[p].dis[j]>z+disvb){
                            Tree[p].dis[j]=z+disvb;
                            Tree[p].FN[j]=false;
                            Tree[p].cnt[j]=1;
                        }else if(Tree[p].dis[j]==z+disvb){
                            Tree[p].cnt[j]+=1;
                        }
                        break;
                    }
                }
                if(!flag){
                    cout<<"There is non-partition vertex in the adjacency list! "<<x<<endl;
                    exit(1);
                }
            }
            else{//if the neighbor is ancestor
                if(k!=j){
                    if(k<j)
                        z=Tree[rank[y]].dis[k];//get the distance to the lower-order vertex
                    else if(k>j)
                        z=Tree[rank[x]].dis[j];

                    if(Tree[p].dis[j]>z+disvb){
                        Tree[p].dis[j]=z+disvb;
                        Tree[p].FN[j]=false;
                        Tree[p].cnt[j]=1;
                    }else if(Tree[p].dis[j]==z+disvb){
                        Tree[p].cnt[j]+=1;
                    }
                }
            }

        }
//        int d1= Tree[p].dis[j];
//        int d2= Dijkstra(ID1,y,Neighbor);
//        if(d1!=d2){
//            cout<<"Incorrect! "<<ID1<<" "<<y<<": "<<d1<<" "<<d2<<endl;
//        }
    }
    assert(Tree[p].disInf.size()==InfNum);
    //nested loop
    list.push_back(Tree[p].uniqueVertex);
    for(int i=0;i<Tree[p].ch.size();i++){
        makeIndexDFSCore(Tree[p].ch[i],list,interface,disInfs);
    }
    list.pop_back();
}*/

/*void Graph::makeIndexDFS(int p, vector<int>& list){
	//initialize
	int NeiNum=Tree[p].vert.size();
	Tree[p].pos.assign(NeiNum+1,0);
	Tree[p].dis.assign(list.size(),INF);
	Tree[p].cnt.assign(list.size(),0);
	Tree[p].FN.assign(list.size(),true);

	//pos
	//map<int,Nei> Nmap; Nmap.clear();//shortcut infor ordered by the pos ID
	for(int i=0;i<NeiNum;i++){
		for(int j=0;j<list.size();j++){
			if(Tree[p].vert[i].first==list[j]){
				Tree[p].pos[i]=j;
				Tree[p].dis[j]=Tree[p].vert[i].second.first;
				Tree[p].cnt[j]=1;
				break;
			}
		}
	}
	Tree[p].pos[NeiNum]=list.size();


	//dis
	for(int i=0;i<NeiNum;i++){
		int x=Tree[p].vert[i].first;
		int disvb=Tree[p].vert[i].second.first;
		int k=Tree[p].pos[i];//the kth ancestor is x

		for(int j=0;j<list.size();j++){
			int y=list[j];//the jth ancestor is y

			int z;//the distance from x to y
			if(k!=j){
				if(k<j)
					z=Tree[rank[y]].dis[k];
				else if(k>j)
					z=Tree[rank[x]].dis[j];

				if(Tree[p].dis[j]>z+disvb){
					Tree[p].dis[j]=z+disvb;
					Tree[p].FN[j]=false;
					Tree[p].cnt[j]=1;
				}else if(Tree[p].dis[j]==z+disvb){
					Tree[p].cnt[j]+=1;
				}
			}
		}
	}

	//nested loop
	list.push_back(Tree[p].uniqueVertex);
	for(int i=0;i<Tree[p].ch.size();i++){
		makeIndexDFS(Tree[p].ch[i],list);
	}
	list.pop_back();
}*/

int Graph::QueryH2H(int ID1,int ID2){
	if(ID1==ID2) return 0;
	if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
	int r1=rank[ID1], r2=rank[ID2];
	int LCA=LCAQuery(r1,r2);

	if(LCA==r1)
		return Tree[r2].dis[Tree[r1].pos.back()];
	else if(LCA==r2)
		return Tree[r1].dis[Tree[r2].pos.back()];
	else{
		int tmp=INF;
		for(int i=0;i<Tree[LCA].pos.size();i++){
			if(tmp>Tree[r1].dis[Tree[LCA].pos[i]]+Tree[r2].dis[Tree[LCA].pos[i]])
				tmp=Tree[r1].dis[Tree[LCA].pos[i]]+Tree[r2].dis[Tree[LCA].pos[i]];
		}
		return tmp;
	}
}

int Graph::LCAQuery(int _p, int _q){
	int p = toRMQ[_p], q = toRMQ[_q];
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
	if (Tree[RMQIndex[k][p]].height < Tree[RMQIndex[k][q]].height)
		return RMQIndex[k][p];
	else return RMQIndex[k][q];
}

//H2H index update for either vertex is periphery vertex
void Graph::DecreaseH2HNew(int a,int b, int newW, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifParallel){// Neighbors is useless in edge decrease update
    map<int,int> checkedDis;//map<tree node ID, distance index>

    for(int i=0;i<Tree.size();i++){
        Tree[i].DisRe.clear();
    }

    int lid,hid;
    if(NodeOrder[a]<NodeOrder[b]){
        lid=a;hid=b;
    }else{
        lid=b;hid=a;
    }
//    cout<<"CoreTag[a]: "<<CoreTag[a]<<" ; CoreTag[b]: "<<CoreTag[b]<<endl;
    int IniH=Tree[rank[lid]].height;//the height where weight change begins
    int ProH=Tree[rank[lid]].height; int ProID=lid;
    vector<set<int>> SCre;//record the shortcut change in each height
    set<int> ss;//ss.clear();
    SCre.assign(ProH+1,ss);

    //map<int,set<int>> DisRe;//rankid; record the distance change caused by the shortcut in each height
    //DisRe.clear();

    int pid = CoreTag[lid];
    int rootVertex = BoundVertex[pid][BoundVertex[pid].size()-1];
    assert(Tree[rank[rootVertex]].pa == 0);

    int MinH;

    bool tri=false;
    for(int i=0;i<Tree[rank[lid]].vert.size();i++){
        if(Tree[rank[lid]].vert[i].first==hid){//if the ancestor is b
            if(Tree[rank[lid]].vert[i].second.first>newW){
                Tree[rank[lid]].vert[i].second.first=newW;
                Tree[rank[lid]].vert[i].second.second=1;
                tri=true;
                SCre[ProH].insert(hid);
                MinH=IniH;
            }else if(Tree[rank[lid]].vert[i].second.first==newW){
                Tree[rank[lid]].vert[i].second.second+=1;
            }
            break;
        }
    }

    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distance labeling has changed

    vector<int> ProIDRecord; ProIDRecord.assign(ProH+1,0);
//cout<<"33333333333"<<endl;
    //int ProBeginH;
    int ProBeginID;
    bool ProBeginIDValue=false;
    set<pair<int,int>> changedCorePairs; changedCorePairs.clear();
    bool flagHigher = false;
    set<int> ProBeginIDInf;
    int MinHInf=INF;
    set<int> affectedParti;

    if(tri){
        //cout<<"Bottom-up ;;;;;;;;;;;;;;;;;; "<<endl;
        /// Bottom-up SE update
        while(ProH>=MinH){
            if(CoreTag[ProID] == -1)//if ProID is interface vertex
                break;
            ProIDRecord[ProH]=ProID;
            vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;//the accessory vertex of lower-order tree node
            bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
            for(auto it=SCre[ProH].begin();it!=SCre[ProH].end();it++){
                int Cid=*it; //higher-order tree node which is the endpoint of a changed shortcut
                int Cw=INF;//=OCdis[make_pair(ProID,Cid)]; //the weight of changed shortcut
                int cidH=Tree[rank[Cid]].height-1;

                map<int,int> Hnei; //Hnei.clear(); // higher-order neighbors
                vector<pair<int,int>> Lnei; //Lnei.clear(); // lower-order neighbors
                for(int j=0;j<Vert.size();j++){
                    if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                        Hnei[Vert[j].first]=Vert[j].second.first;//Hnei may contain the interface vertex of ProID
                    }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                        Lnei.emplace_back(make_pair(Vert[j].first,Vert[j].second.first));//Lnei contains the neighbors of ProID which has lower order than Cid
                    }else{
                        Cw=Vert[j].second.first;//vert distance from ProID to Cid
                    }
                }

                ///deal with the interface vertex
                if(CoreTag[Cid] == -1){
                    if(Tree[rank[ProID]].height<MinHInf){
                        MinHInf = Tree[rank[ProID]].height;
                        ProBeginIDInf.clear();
                        ProBeginIDInf.insert(ProID);
                    }else if(Tree[rank[ProID]].height==MinHInf){
                        ProBeginIDInf.insert(ProID);
                    }
                    int wsum,lidHeight;
                    /// For Hnei, update the super edges between Cid and other interface vertex, i.e., the shortcuts in AdjaCore
                    for(int i=0;i<AdjaCore[Cid].size();++i){
                        hid=AdjaCore[Cid][i].first;
                        if(Hnei.find(hid)!=Hnei.end()){//if hid is also the higher-order neighbor of ProID (a)
                            wsum=Cw+Hnei[hid];
                            if(AdjaCoreMap[Cid].find(hid) == AdjaCoreMap[Cid].end()){
                                cout<<"Cannot find hid in Cid's adjacency list!!"<<endl;
                                exit(1);
                            }else if(AdjaCoreMap[hid].find(Cid) == AdjaCoreMap[hid].end()){
                                cout<<"Cannot find Cid in hid's adjacency list!!"<<endl;
                                exit(1);
                            }
//                            if(wsum<AdjaCore[Cid][i].second){//update neighbor distance of Cid
                            if(wsum<AdjaCoreMap[Cid][hid]){//update neighbor distance of Cid

                                assert(AdjaCoreMapOld[Cid][hid] == AdjaCore[Cid][i].second);
                                assert(AdjaCoreMapOld[hid][Cid] == AdjaCore[Cid][i].second);
                                AdjaCoreMap[Cid][hid]=wsum;
                                AdjaCoreMap[hid][Cid]=wsum;
//                                AdjaCore[Cid][i].second=wsum;
//                                for(int j=0;j<AdjaCore[hid].size();++j){
//                                    if(AdjaCore[hid][j].first == Cid){
//                                        AdjaCore[hid][j].second = wsum;
//                                    }
//                                }
                                ///identify the affected partitions
                                if(SuppPartiID[Cid][hid].size()>1){
                                    for(auto it=SuppPartiID[Cid][hid].begin();it!=SuppPartiID[Cid][hid].end();++it){
                                        if(it->first != pid){
                                            affectedParti.insert(it->first);
                                        }
                                    }
                                }

                            }

                        }
                    }
                    /// For Lnei, update the super edges of Lnei vertex
                    for(int j=0;j<Lnei.size();j++){
                        lid=Lnei[j].first;
                        if(CoreTag[lid] == -1){//if lid is interface vertex
                            for(int i=0;i<AdjaCore[lid].size();++i) {
                                int vertid = AdjaCore[lid][i].first;
                                if (vertid == Cid) {//Only deal with Cid for Lnei vertex; Tree[rank[lid]].vert[k].first
                                    wsum=Cw+Lnei[j].second;
                                    if(AdjaCoreMap[Cid].find(lid) == AdjaCoreMap[Cid].end()){
                                        cout<<"Cannot find lid in Cid's adjacency list!!"<<endl;
                                        exit(1);
                                    }else if(AdjaCoreMap[lid].find(Cid) == AdjaCoreMap[lid].end()){
                                        cout<<"Cannot find Cid in lid's adjacency list!!"<<endl;
                                        exit(1);
                                    }
                                    if(AdjaCoreMap[lid][Cid]>wsum){//update neighbor distance of Lnei
                                        assert(AdjaCoreMapOld[Cid][lid] == AdjaCore[lid][i].second);
                                        assert(AdjaCoreMapOld[lid][Cid] == AdjaCore[lid][i].second);
                                        AdjaCoreMap[lid][Cid]=wsum;
                                        AdjaCoreMap[Cid][lid]=wsum;
//                                        AdjaCore[lid][i].second=wsum;
//                                        for(int k=0;k<AdjaCore[Cid].size();++k){
//                                            if(AdjaCore[Cid][k].first == lid){
//                                                AdjaCore[Cid][k].second = wsum;
//                                            }
//                                        }
                                        ///identify the affected partitions
                                        if(SuppPartiID[lid][Cid].size()>1){
                                            for(auto it=SuppPartiID[lid][Cid].begin();it!=SuppPartiID[lid][Cid].end();++it){
                                                if(it->first != pid){
                                                    affectedParti.insert(it->first);
                                                }
                                            }
                                        }
                                    }
                                    break;
                                }
                            }
                        }else{//if lid is periphery vertex
                            lidHeight=Tree[rank[lid]].height-1;
                            for(int k=0;k<Tree[rank[lid]].vert.size();k++){
                                int vertid=Tree[rank[lid]].vert[k].first;
                                if(vertid==Cid){//Only deal with Cid for Lnei vertex; Tree[rank[lid]].vert[k].first
                                    wsum=Cw+Lnei[j].second;
                                    if(Tree[rank[lid]].vert[k].second.first>wsum){//update neighbor distance of Lnei
                                        Tree[rank[lid]].vert[k].second.first=wsum;
                                        Tree[rank[lid]].vert[k].second.second=1;
                                        SCre[Tree[rank[lid]].height].insert(Cid);//if lid has lower order than ProID, the inserted Cid will not be processed
                                        if(Tree[rank[lid]].height<MinH)
                                            MinH=Tree[rank[lid]].height;

                                        if(Tree[rank[lid]].height<MinHInf){
                                            MinHInf = Tree[rank[lid]].height;
                                            ProBeginIDInf.clear();
                                            ProBeginIDInf.insert(lid);
                                        }else if(Tree[rank[lid]].height==MinHInf){
                                            ProBeginIDInf.insert(lid);
                                        }
                                        /// interface
                                        if(CoreTag[vertid] == -1){//if the neighbor is interface vertex, never enter this branch
                                            assert(Tree[rank[lid]].disInf.find(vertid) != Tree[rank[lid]].disInf.end());
                                            if(Tree[rank[lid]].disInf[vertid] > wsum){
                                                Tree[rank[lid]].DisReInf.insert(vertid);//record the vertex id that the interface label should be updated
//                                                Tree[rank[lid]].disInf[vertid] = wsum;
                                            }
                                        }
                                    }else if(Tree[rank[lid]].vert[k].second.first==wsum){
                                        Tree[rank[lid]].vert[k].second.second+=1;
                                    }

                                    break;
                                }
                            }
                        }
                    }
                    continue;
                }

//                if(Tree[rank[ProID]].dis[cidH]>=Cw){//if the distance to ancestor Cid has changed
                if(Tree[rank[ProID]].dis[cidH]>Cw){//if the distance to ancestor Cid has changed
                    Tree[rank[ProID]].dis[cidH]=Cw;
                    Tree[rank[ProID]].FN[cidH]=true;
                    ProIDdisCha=true;
                    Tree[rank[ProID]].DisRe.insert(Cid);//record the vertex id that the distance label should be updated
                    //DisRe[rank[ProID]].insert(Cid); //cout<<"dischange Cid "<<Cid<<endl;
                }else if(Tree[rank[ProID]].dis[cidH]==Cw){
                    cout<<"Equal distance............."<<endl;
                }
                /// For Hnei, update the super edges of Cid (b)
                int hid,hidHeight,lid,lidHeight,wsum;
                for(int j=0;j<Tree[rank[Cid]].vert.size();j++){
                    hid=Tree[rank[Cid]].vert[j].first;hidHeight=Tree[rank[hid]].height-1;
                    if(Hnei.find(hid)!=Hnei.end()){//if hid is also the higher-order neighbor of ProID (a)
                        wsum=Cw+Hnei[hid];
                        if(wsum<Tree[rank[Cid]].vert[j].second.first){//update neighbor distance of Cid
                            Tree[rank[Cid]].vert[j].second.first=wsum;
                            Tree[rank[Cid]].vert[j].second.second=1;
                            SCre[Tree[rank[Cid]].height].insert(hid);
                            if(Tree[rank[Cid]].height<MinH)
                                MinH=Tree[rank[Cid]].height;

                            /// interface
                            if(CoreTag[hid] == -1){//if the neighbor is interface vertex
                                assert(Tree[rank[Cid]].disInf.find(hid) != Tree[rank[Cid]].disInf.end());
                                if(Tree[rank[Cid]].disInf[hid] > wsum){
                                    Tree[rank[Cid]].DisReInf.insert(hid);//record the vertex id that the interface label should be updated
                                    if(Tree[rank[Cid]].height<MinHInf){
                                        MinHInf = Tree[rank[Cid]].height;
                                        ProBeginIDInf.clear();
                                        ProBeginIDInf.insert(Cid);
                                    }else if(Tree[rank[Cid]].height==MinHInf){
                                        ProBeginIDInf.insert(Cid);
                                    }
//                                    Tree[rank[Cid]].disInf[hid] = wsum;
                                }
                            }
                        }else if(wsum==Tree[rank[Cid]].vert[j].second.first){
                            Tree[rank[Cid]].vert[j].second.second+=1;
                        }

                    }
                }
                /// For Lnei, update the super edges of Lnei vertex
                for(int j=0;j<Lnei.size();j++){
                    lid=Lnei[j].first;lidHeight=Tree[rank[lid]].height-1;
                    for(int k=0;k<Tree[rank[lid]].vert.size();k++){
                        int vertid=Tree[rank[lid]].vert[k].first;
                        if(vertid==Cid){//Only deal with Cid for Lnei vertex; Tree[rank[lid]].vert[k].first
                            wsum=Cw+Lnei[j].second;
                            if(Tree[rank[lid]].vert[k].second.first>wsum){//update neighbor distance of Lnei
                                Tree[rank[lid]].vert[k].second.first=wsum;
                                Tree[rank[lid]].vert[k].second.second=1;
                                SCre[Tree[rank[lid]].height].insert(Cid);//if lid has lower order than ProID, the inserted Cid will not be processed
                                if(Tree[rank[lid]].height<MinH)
                                    MinH=Tree[rank[lid]].height;

                            }else if(Tree[rank[lid]].vert[k].second.first==wsum){
                                Tree[rank[lid]].vert[k].second.second+=1;
                            }

                            break;
                        }
                    }
                }
            }

            if(ProIDdisCha){//if the distance labeling is detected changed
                vertexIDChL.insert(ProID);//the vertex set that the distance label has changed
                //ProBeginH=ProH;
                ProBeginID=ProID;
                ProBeginIDValue=true;
                //cout<<"PorBeginID "<<ProBeginID<<endl;
            }

            ProH-=1;
            ProID=Tree[Tree[rank[ProID]].pa].uniqueVertex;//propagate to parent vertex
        }
        assert(MinH > 1);

        /// Top-down label update

        //top-down process
//        flag_fromRoot = true;
        vector<int> interface;
        for(auto it=Tree[rank[rootVertex]].vert.begin();it!=Tree[rank[rootVertex]].vert.end();++it){
            interface.emplace_back(it->first);
        }
//        if(!ProBeginIDInf.empty()){
//            int minH = MinHInf;
//            if(minH <= MinH){
//                for(auto it=ProBeginIDInf.begin();it!=ProBeginIDInf.end();++it){
//                    int rid=*it;
//                    InterfacePropagate( rank[rid],interface,Tree,false);
//                }
//                flagHigher = true;
//            }else{
//                InterfacePropagate( rank[ProBeginID],interface,Tree,false);
//            }
//        }else{
//            InterfacePropagate( rank[ProBeginID],interface,Tree,false);
//        }

        PartiUpdateExt[pid] = true;
        if(affectedParti.size()>=1){
//            cout<<"Flag 3: affectedParti number: "<<affectedParti.size()<<endl;
            if(ifParallel){//use multi-thread
                boost::thread_group thread;
                vector<int> vParti;
                vParti.push_back(pid);
                for(auto it=affectedParti.begin();it!=affectedParti.end();++it){
                    vParti.push_back(*it);
                    PartiUpdateExt[*it] = true;//if the interface entry needs to be updated, we may also need to update the extended labels
                }
                cout<<"vParti size: "<<vParti.size()<<endl;
                if(vParti.size()>threadnum){
                    int step=vParti.size()/threadnum;
                    boost::thread_group threadf;
                    for(int i=0;i<threadnum;i++){
                        pair<int,int> p;
                        p.first=i*step;
                        if(i==threadnum-1)
                            p.second=vParti.size();
                        else
                            p.second=(i+1)*step;
                        threadf.add_thread(new boost::thread(&Graph::InterfacePropagateParallel, this, p, boost::ref(vParti), false));
                    }
                    threadf.join_all();
                }else{
                    boost::thread_group threadf;
                    for(int i=0;i<vParti.size();i++){
                        threadf.add_thread(new boost::thread(&Graph::InterfacePropagateParallel, this, make_pair(i,i+1), boost::ref(vParti), false));
                    }
                    threadf.join_all();
                }

            }else{//use single-thread
                InterfacePropagate( rank[rootVertex], interface,Tree,false);//propagate from root vertex

                if(!affectedParti.empty()){//if we need to update the shortcut of other partition
                    for(auto it=affectedParti.begin();it!=affectedParti.end();++it){
                        int Pid = *it;
                        int rootID = BoundVertex[Pid][BoundVertex[Pid].size()-1];
                        vector<int> interfaceP;
                        PartiUpdateExt[*it] = true;//if the interface entry needs to be updated, we may also need to update the extended labels
                        for(auto it2=Tree[rank[rootID]].vert.begin();it2!=Tree[rank[rootID]].vert.end();++it2){
                            interfaceP.emplace_back(it2->first);
                        }
                        InterfacePropagate( rank[rootID], interfaceP,Tree,false);
                    }
                }
            }
        }else{//use single thread
            InterfacePropagate( rank[rootVertex], interface,Tree,false);//propagate from root vertex

            if(!affectedParti.empty()){//if we need to update the shortcut of other partition
                for(auto it=affectedParti.begin();it!=affectedParti.end();++it){
                    int Pid = *it;
                    int rootID = BoundVertex[Pid][BoundVertex[Pid].size()-1];
                    vector<int> interfaceP;
                    PartiUpdateExt[*it] = true;//if the interface entry needs to be updated, we may also need to update the extended labels
                    for(auto it2=Tree[rank[rootID]].vert.begin();it2!=Tree[rank[rootID]].vert.end();++it2){
                        interfaceP.emplace_back(it2->first);
                    }
                    InterfacePropagate( rank[rootID], interfaceP,Tree,false);
                }
            }
        }


//        if(flag_fromRoot){//update interface label
//            InterfacePropagate( rank[rootVertex], interface,Tree,false);
//        }

        if(ProBeginIDValue){//update ancestor label
            vector<int> ancestors; //linee.clear(); // ancestors
            ancestors.reserve(heightMax);

            int pachidd=Tree[Tree[rank[ProBeginID]].pa].uniqueVertex;
            //cout<<"aaaaaaaaaaaa "<<pachidd<<" "<<Tree[rank[pachidd]].height<<endl;
//            while(Tree[rank[pachidd]].height>1){
            while(pachidd != -1){
                ancestors.insert(ancestors.begin(),pachidd);
                pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
                //cout<<"height "<<Tree[rank[pachidd]].height<<endl;

            }

            ancestors.insert(ancestors.begin(),pachidd);//begin from virtual root

            EachNodeProBDis5New(rank[ProBeginID], ancestors, interface, vertexIDChL, checkedDis, Tree, rank);
//            EachNodeProBDis5(rank[ProBeginID], ancestors, vertexIDChL, checkedDis, Tree, rank);
        }
    }
    //return checkedDis.size();
}

//function for the interface entry propagation
void Graph::InterfacePropagate(int child, vector<int>& interfaces, vector<Node> &Tree, bool ifIncrease){
//    if(child == rank[210695]){
//        cout<<210695<<" "<<Tree[Tree[child].treeroot].uniqueVertex<<endl;
//    }
    // Solution 1: check all interface entries
    for(int j=0;j<interfaces.size();j++){//for each interface vertex
        int ID2=interfaces[j];
//        if(child == rank[197979] && ID2 == 144762)
//            cout<<ID2<<endl;
        if(ifIncrease){
            Tree[child].disInf[ID2] = INF;
        }
        for(int k=0;k<Tree[child].vert.size();k++) {
            int b = Tree[child].vert[k].first;
            int vbW = Tree[child].vert[k].second.first;
            if(b == interfaces[j]){
                if(Tree[child].disInf[ID2] > vbW){
                    Tree[child].disInf[ID2]=vbW;//
                }
                continue;
            }
            int z = INF;
            if(rank[b] > 0 ){//if x is not boundary vertex, i.e., if it is ancestor
                z = Tree[rank[b]].disInf[ID2];
                if(Tree[child].disInf[ID2]>z+vbW){
                    Tree[child].disInf[ID2]=z+vbW;
                }
            }else{//if it is interface vertex
                if(AdjaCoreMap[b].find(ID2) != AdjaCoreMap[b].end()){
                    z = AdjaCoreMap[b][ID2];
                }
                assert(z>0);
                if(Tree[child].disInf[ID2]>z+vbW){// never enter this branch
                    Tree[child].disInf[ID2]=z+vbW;
                }
            }
        }
    }
    // Solution 2
    /*if(!Tree[child].DisReInf.empty()){
        for(int j=0;j<interfaces.size();j++) {//for each interface vertex
            int ID2=interfaces[j];
            if((Tree[child].DisReInf.find(ID2) != Tree[child].DisReInf.end()) || Tree[child].FNInf[j]){//if we find ID2 in DisReInf or the interface distance is obtained from ancestors
                for(int k=0;k<Tree[child].vert.size();k++) {
                    int b = Tree[child].vert[k].first;
                    int vbW = Tree[child].vert[k].second.first;
                    if(b == ID2){
                        if(Tree[child].disInf[ID2] > vbW){
                            Tree[child].disInf[ID2]=vbW;//
                            Tree[child].FNInf[j]=false;
                        }
                        continue;
                    }
                    int z = INF;
                    if(rank[b] > 0 ){//if x is not boundary vertex, i.e., if it is ancestor
                        z = Tree[rank[b]].disInf[ID2];
                        if(Tree[child].disInf[ID2]>z+vbW){
                            Tree[child].disInf[ID2]=z+vbW;
                            Tree[child].FNInf[j]=true;
                        }
                    }else{//if it is interface vertex
                        if(AdjaCoreMap[b].find(ID2) != AdjaCoreMap[b].end()){
                            z = AdjaCoreMap[b][ID2];
                        }
                        assert(z>0);
                        if(Tree[child].disInf[ID2]>z+vbW){// never enter this branch
                            Tree[child].disInf[ID2]=z+vbW;
                            Tree[child].FNInf[j]=false;
                        }
                    }
                }
            }
        }
    }
    else{
        for(int j=0;j<interfaces.size();j++){//for each interface vertex
            int ID2=interfaces[j];
//            if(Tree[child].FNInf[j]){//if from ancestors
                for(int k=0;k<Tree[child].vert.size();k++) {
                    int b = Tree[child].vert[k].first;
                    int vbW = Tree[child].vert[k].second.first;
                    if(b == interfaces[j]){
                        if(Tree[child].disInf[ID2] > vbW){
                            Tree[child].disInf[ID2]=vbW;//
                            Tree[child].FNInf[j]=false;
                        }
                        continue;
                    }
                    int z = INF;
                    if(rank[b] > 0 ){//if x is not boundary vertex, i.e., if it is ancestor
                        z = Tree[rank[b]].disInf[ID2];
                        if(Tree[child].disInf[ID2]>z+vbW){
                            Tree[child].disInf[ID2]=z+vbW;
                            Tree[child].FNInf[j]=true;
                        }
                    }else{//if it is interface vertex
                        if(AdjaCoreMap[b].find(ID2) != AdjaCoreMap[b].end()){
                            z = AdjaCoreMap[b][ID2];
                        }
                        assert(z>0);
                        if(Tree[child].disInf[ID2]>z+vbW){// never enter this branch
                            Tree[child].disInf[ID2]=z+vbW;
                            Tree[child].FNInf[j]=false;
                        }
                    }
                }
//            }
        }
    }*/

    for(int i=0;i<Tree[child].ch.size();i++){
        InterfacePropagate(Tree[child].ch[i], interfaces,Tree,ifIncrease);
    }
}

//function for the interface entry propagation: parallel version
void Graph::InterfacePropagateParallel(pair<int,int> pRange, vector<int>& pids, bool ifIncrease){
    for(int i=pRange.first;i<pRange.second;++i){
        int Pid = pids[i];
        int rootID = BoundVertex[Pid][BoundVertex[Pid].size()-1];
        vector<int> interfaceP;
        for(auto it=Tree[rank[rootID]].vert.begin();it!=Tree[rank[rootID]].vert.end();++it){
            interfaceP.emplace_back(it->first);
        }
        InterfacePropagate( rank[rootID], interfaceP,Tree,ifIncrease);
    }
}
// top-down label update: interface version
void Graph::EachNodeProBDis5New(int child,vector<int>& line, vector<int>& interfaces, set<int>& vertexIDChL, map<int,int>& checkedDis, vector<Node> &Tree, vector<int> &rank){
    bool ProIDdisCha=false;
    /// update ancestor entry
    if(!Tree[child].DisRe.empty()){//Tree[child].DisRe.size()!=0
        /// update ancestor entries
        for(int k=0;k<Tree[child].vert.size();k++){
            int b=Tree[child].vert[k].first;
            if(CoreTag[b] == -1){//if it is interface vertex, check whether the distance label may be affected by interface
                int vbW=Tree[child].vert[k].second.first;
                for(int i=1;i<line.size();i++){
                    assert(Tree[rank[line[i]]].disInf.find(b) != Tree[rank[line[i]]].disInf.end());
                    if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].disInf[b]){
                        Tree[child].dis[i]=vbW+Tree[rank[line[i]]].disInf[b];
                        Tree[child].FN[i]=false;
                        ProIDdisCha=true;
                    }
                }
            }else{//if not interface vertex
                int bH=Tree[rank[b]].height-1,vbW=Tree[child].vert[k].second.first;
                if(Tree[child].FN[bH]){//if distance label from child to b is directly sourced from their shortcut
                    if(Tree[child].DisRe.find(b)!=Tree[child].DisRe.end()){//if found b, all ancestor check
                        for(int i=1;i<bH;i++){
                            checkedDis.insert(make_pair(child,i));
                            if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){//update ancestor distance label
                                Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
                                Tree[child].FN[i]=false;
                                ProIDdisCha=true;
                            }
                        }
                        for(int i=bH+1;i<line.size();i++){
                            checkedDis.insert(make_pair(child,i));
                            if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
                                Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
                                Tree[child].FN[i]=false;
                                ProIDdisCha=true;
                            }
                        }

                    }else{//partial ancestor check if we cannot find b

                        if(vertexIDChL.find(b)!=vertexIDChL.end()){
                            for(int i=1;i<bH;i++){
                                checkedDis.insert(make_pair(child,i));
                                if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
                                    Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
                                    Tree[child].FN[i]=false;
                                    ProIDdisCha=true;
                                }
                            }
                        }
                        for(int i=bH+1;i<line.size();i++){
                            checkedDis.insert(make_pair(child,i));
                            if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
                                Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
                                Tree[child].FN[i]=false;
                                ProIDdisCha=true;
                            }
                        }

                    }
                }
            }


        }
    }
    else{// if there is no label for checking
        for(int k=0;k<Tree[child].vert.size();k++){
            int b=Tree[child].vert[k].first;
            if(CoreTag[b] == -1){//if it is interface vertex, check whether the distance label may be affected by interface
                int vbW=Tree[child].vert[k].second.first;
                for(int i=1;i<line.size();i++){
                    assert(Tree[rank[line[i]]].disInf.find(b) != Tree[rank[line[i]]].disInf.end());
                    if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].disInf[b]){
                        Tree[child].dis[i]=vbW+Tree[rank[line[i]]].disInf[b];
                        Tree[child].FN[i]=false;
                        ProIDdisCha=true;
                    }
                }
            }else {//if not interface vertex
                int bH = Tree[rank[b]].height - 1, vbW = Tree[child].vert[k].second.first;
                if (Tree[child].FN[bH]) {//Property 5
                    if (vertexIDChL.find(b) != vertexIDChL.end()) {//if the distance label of b is changed
                        for (int i = 1; i < bH; i++) {//check ancestor from 0 to bH
                            checkedDis.insert(make_pair(child, i));
                            if (Tree[child].dis[i] > vbW + Tree[rank[b]].dis[i]) {
                                Tree[child].dis[i] = vbW + Tree[rank[b]].dis[i];
                                Tree[child].FN[i] = false;
                                ProIDdisCha = true;
                            }
                        }
                    }
                    for (int i = bH + 1; i < line.size(); i++) {
                        checkedDis.insert(make_pair(child, i));
                        if (Tree[child].dis[i] > vbW + Tree[rank[line[i]]].dis[bH]) {
                            Tree[child].dis[i] = vbW + Tree[rank[line[i]]].dis[bH];
                            Tree[child].FN[i] = false;
                            ProIDdisCha = true;
                        }
                    }
                }
            }
        }
    }

    if(ProIDdisCha){
        vertexIDChL.insert(Tree[child].uniqueVertex);
    }

    line.push_back(Tree[child].uniqueVertex);
    for(int i=0;i<Tree[child].ch.size();i++){
        EachNodeProBDis5New(Tree[child].ch[i], line, interfaces,vertexIDChL,checkedDis,Tree, rank);
    }
    line.pop_back();

}

//New function for H2H increase update: interface version
void Graph::IncreaseH2HNew(int a,int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, bool ifParallel){
    int ChangeNum=0;
    map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
    //OCdis.clear();

    int lid,hid;
    if(NodeOrder[a]<NodeOrder[b]){
        lid=a;hid=b;
    }else{
        lid=b;hid=a;
    }
    cout<<"CoreTag[a]: "<<CoreTag[a]<<" ; CoreTag[b]: "<<CoreTag[b]<<endl;

    int lowestH = INF;
    if(CoreTag[a]!=-1){
//        cout<<"Height of a: "<<Tree[rank[a]].height;
        if(lowestH>Tree[rank[a]].height)
            lowestH = Tree[rank[a]].height;
    }
    if(CoreTag[b]!=-1){
//        cout<<" ; Height of b: "<<Tree[rank[b]].height;
        if(lowestH>Tree[rank[b]].height)
            lowestH = Tree[rank[b]].height;
    }
//    cout<<endl;

//    cout<<"lid: "<<lid<<" ; hid: "<<hid<<endl;
    int IniH=Tree[rank[lid]].height;//the height where weight change begins
    int ProH=Tree[rank[lid]].height; int ProID=lid;
    vector<set<int>> SCre;//record the shortcut change in each height
    set<int> vec; //vec.clear();
    SCre.assign(IniH+1,vec);

    int pid = CoreTag[lid];
    int rootVertex = BoundVertex[pid][BoundVertex[pid].size()-1];

    if(Tree[Tree[rank[lid]].treeroot].uniqueVertex != rootVertex){
        cout<<"Incorrect root vertex! "<<rootVertex<<" "<<Tree[Tree[rank[lid]].treeroot].uniqueVertex<<" "<<Tree[rank[lid]].treeroot<<endl;
        rootVertex = Tree[Tree[rank[lid]].treeroot].uniqueVertex;
    }
    assert(Tree[rank[rootVertex]].pa == 0);

    int MinH;
    set<int> ProBeginIDInf;
    int MinHInf = INF;//the minimum tree height for interface entry update

    vector<int> line; //line.clear();//ancestors of lid
    line.reserve(heightMax);
    int pachid=ProID;
    while(pachid!=-1){
        line.insert(line.begin(),pachid);
        pachid=Tree[Tree[rank[pachid]].pa].uniqueVertex;
    }
    line.insert(line.begin(),pachid);

    bool tri=false;
    for(int i=0;i<Tree[rank[lid]].vert.size();i++){
        if(Tree[rank[lid]].vert[i].first==hid){//hid must be the vert neighbor of lid
            if(Tree[rank[lid]].vert[i].second.first==oldW){
                Tree[rank[lid]].vert[i].second.second-=1;
                if(Tree[rank[lid]].vert[i].second.second<1){
                    OCdis[make_pair(lid,hid)]=oldW;
                    SCre[ProH].insert(hid);
                    MinH=IniH;
                    tri=true;//cout<<"Trigger the Shortcut Change or not? "<<tri<<endl;
                }
            }
            break;
        }
    }

    bool influence; int ProBeginID;
    bool ProBeginIDValue=false;
    bool flagHigher = false;
    set<int> affectedParti;

//    cout<<"Height: "<<Tree[rank[210695]].height<<" "<<Tree[rank[208312]].height<<" "<<Tree[rank[210694]].height<<endl;

    if(tri){
        set<pair<int,int>> SECoreForUpdate;//the set of shortcuts among interface vertices for updating, (s,t) (order[s]<order[t])
        SECoreForUpdate.clear();
        //shortcut update
        while(ProH>=MinH){
            if(CoreTag[ProID] == -1)//if ProID is interface vertex
                break;
            influence=false;
            vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;

//            if(ProID == 210695 || ProID == 208312 || ProID == 210694){
//                cout<<ProID<<endl;
//            }
            /// deal with each affected super edge
            for(auto it=SCre[ProH].begin();it!=SCre[ProH].end();it++){//shortcut e(ProID,Cid) is the affected super edge
                ///Step 1: identify the affected super edges caused by the update of e(ProID,Cid)
                int Cid=*it;
                int Cw=OCdis[make_pair(ProID,Cid)];//old weight of the affected super edge
                int cidH;

                map<int,int> Hnei; //Hnei.clear();
                vector<pair<int,int>> Lnei; //Lnei.clear();
                for(int j=0;j<Vert.size();j++){//classify the vert neighbors of ProID into Hnei and Lnei
                    if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                        Hnei[Vert[j].first]=Vert[j].second.first;
                    }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                        Lnei.emplace_back(make_pair(Vert[j].first,Vert[j].second.first));
                    }
                }

                ///identify the related super edges for updating
                if(CoreTag[Cid] == -1)///deal with the interface vertex
                {
//                    if(Cid == 48300 || Cid == 27115){
//                        cout<<"!!! "<<Cid<<endl;
//                    }

                    if(Tree[rank[ProID]].height<MinHInf){
                        MinHInf = Tree[rank[ProID]].height;
                        ProBeginIDInf.clear();
                        ProBeginIDInf.insert(ProID);
                    }else if(Tree[rank[ProID]].height==MinHInf){
                        ProBeginIDInf.insert(ProID);
                    }
                    int wsum,lidHeight;
                    /// For Hnei, update the super edges between Cid and other interface vertex, i.e., the shortcuts in AdjaCore
                    if(!Hnei.empty()){
                        for(int i=0;i<AdjaCore[Cid].size();++i){
                            hid=AdjaCore[Cid][i].first;
                            if(Hnei.find(hid)!=Hnei.end()){
                                wsum=Cw+Hnei[hid];
//                                if(hid == 48300 || hid == 27115)
//                                    cout<<"1 "<<Cid<<" "<<hid<<endl;
                                if(wsum==AdjaCore[Cid][i].second){
                                    //since we do not know whether e(Cid,hid) is only supported by wsum, we will check it after all in-periphery shortcuts are processed
                                    SECoreForUpdate.insert(make_pair(Cid,hid));
//                                SCre[Tree[rank[Cid]].height].insert(hid);
//                                if(Tree[rank[Cid]].height<MinH)
//                                    MinH=Tree[rank[Cid]].height;
//                                    OCdis[make_pair(Cid,hid)]=wsum;
                                }
                            }
                        }
                    }
                    /// For Lnei, update the super edges of Lnei vertex
                    for(int j=0;j<Lnei.size();j++){
                        lid=Lnei[j].first;
                        if(CoreTag[lid] == -1){//if lid is interface vertex
                            for(int i=0;i<AdjaCore[lid].size();++i) {
                                int vertid = AdjaCore[lid][i].first;
                                if (vertid == Cid) {//Only deal with Cid for Lnei vertex; Tree[rank[lid]].vert[k].first
                                    wsum=Cw+Lnei[j].second;
//                                    if(lid == 48300 || lid == 27115)
//                                        cout<<"2 "<<lid<<" "<<Cid<<endl;
                                    if(AdjaCore[lid][i].second==wsum){//update neighbor distance of Lnei
                                        SECoreForUpdate.insert(make_pair(lid,Cid));
                                    }
                                    break;
                                }
                            }
                        }else{//if lid is periphery vertex
                            lidHeight=Tree[rank[lid]].height-1;
                            for(int k=0;k<Tree[rank[lid]].vert.size();k++){
                                int vertid=Tree[rank[lid]].vert[k].first;
                                if(vertid==Cid){//Only deal with Cid for Lnei vertex; Tree[rank[lid]].vert[k].first
                                    wsum=Cw+Lnei[j].second;
                                    if(Tree[rank[lid]].vert[k].second.first==wsum){///check 3: the interface entry of lid
                                        Tree[rank[lid]].vert[k].second.second-=1;
                                        if(Tree[rank[lid]].vert[k].second.second<1){
                                            SCre[Tree[rank[lid]].height].insert(Cid);
                                            if(Tree[rank[lid]].height<MinH)
                                                MinH=Tree[rank[lid]].height;
                                            OCdis[make_pair(lid,Cid)]=Cw+Lnei[j].second;//old shortcut distance
                                            if(Tree[rank[lid]].height<MinHInf){
                                                MinHInf = Tree[rank[lid]].height;
                                                ProBeginIDInf.clear();
                                                ProBeginIDInf.insert(lid);
                                            }else if(Tree[rank[lid]].height==MinHInf){
                                                ProBeginIDInf.insert(lid);
                                            }
                                            /// interface
                                            if(CoreTag[vertid] == -1){//if the neighbor is interface vertex, never enter this branch
                                                assert(Tree[rank[lid]].disInf.find(vertid) != Tree[rank[lid]].disInf.end());
                                                if(Tree[rank[lid]].disInf[vertid] == Cw+Lnei[j].second){
                                                    Tree[rank[lid]].DisReInf.insert(vertid);//record the vertex id that the interface label should be updated
//                                        Tree[rank[lid]].disInf[vertid] = wsum;
                                                }
                                            }
                                        }
                                    }
                                    break;
                                }
                            }
                        }
                    }

                    /// before Cw=d(ProID,Cid) gets its new value, we first check which dis it will be influenced
//                    assert(Tree[rank[ProID]].FNInf.find(Cid) != Tree[rank[ProID]].FNInf.end());
//                    if(Tree[rank[ProID]].FNInf[Cid]){//the distance is directly obtained from this shortcut (vert)
//
//                    }
                }
                else/// if the Cid is not interface vertex, check the affected label related to Cid
                {
                    cidH=Tree[rank[Cid]].height-1;

                    //check the affected shortcuts
                    int hid,lid;
                    //check the label of Cid
                    for(int j=0;j<Tree[rank[Cid]].vert.size();j++){
                        hid=Tree[rank[Cid]].vert[j].first;
                        if(Hnei.find(hid)!=Hnei.end()){//for the neighbor of both ProID and Cid who has higher order than Cid
                            if(Cw+Hnei[hid]==Tree[rank[Cid]].vert[j].second.first){///check 1: check the ancestor entry of Cid
                                Tree[rank[Cid]].vert[j].second.second-=1;
                                if(Tree[rank[Cid]].vert[j].second.second<1){
                                    SCre[Tree[rank[Cid]].height].insert(hid);
                                    if(Tree[rank[Cid]].height<MinH)
                                        MinH=Tree[rank[Cid]].height;
                                    OCdis[make_pair(Cid,hid)]=Cw+Hnei[hid];

                                    /// interface
                                    if(CoreTag[hid] == -1){//if the neighbor is interface vertex
                                        assert(Tree[rank[Cid]].disInf.find(hid) != Tree[rank[Cid]].disInf.end());
                                        if(Tree[rank[Cid]].height<MinHInf){
                                            MinHInf = Tree[rank[Cid]].height;
                                            ProBeginIDInf.clear();
                                            ProBeginIDInf.insert(Cid);
                                        }else if(Tree[rank[Cid]].height==MinHInf){
                                            ProBeginIDInf.insert(Cid);
                                        }
                                        if(Tree[rank[Cid]].disInf[hid] == Cw+Hnei[hid]){
                                            Tree[rank[Cid]].DisReInf.insert(hid);//record the vertex id that the interface label should be updated
//                                    Tree[rank[Cid]].disInf[hid] = wsum;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    //check the label of Lnei (i.e., the vertex has lower order than Cid)
                    for(int j=0;j<Lnei.size();j++){
                        lid=Lnei[j].first;
                        for(int k=0;k<Tree[rank[lid]].vert.size();k++){
                            int vertid=Tree[rank[lid]].vert[k].first;
                            if(vertid==Cid){
                                if(Tree[rank[lid]].vert[k].second.first==Cw+Lnei[j].second){///check 2: check the ancestor entry (Cid) of lid
                                    Tree[rank[lid]].vert[k].second.second-=1;
                                    if(Tree[rank[lid]].vert[k].second.second<1){
                                        SCre[Tree[rank[lid]].height].insert(Cid);
                                        if(Tree[rank[lid]].height<MinH)
                                            MinH=Tree[rank[lid]].height;
                                        OCdis[make_pair(lid,Cid)]=Cw+Lnei[j].second;//old shortcut distance
//                                        cout<<"LidH: "<<Tree[rank[lid]].height<<" ; Lid: "<<lid<<" ; Cid: "<<Cid<<" ; OldW: "<<Tree[rank[lid]].vert[k].second.first<<endl;
                                    }
                                }
                                break;
                            }
                        }
                    }

                    /// before Cw=d(ProID,Cid) gets its new value, we first check which dis it will be influenced
                    if(Tree[rank[ProID]].FN[cidH]){//the distance is directly obtained from this shortcut (vert)
                        influence=true;//as the distance from ProID to Cid is equal to the affected super edge e(ProID,Cid), the update of the super edge will affect the ancestor label
                        //higher than Cid
                        for(int i=0;i<cidH;i++){
                            if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[Cid]].dis[i]){
                                Tree[rank[ProID]].cnt[i]-=1;
                            }
                        }

                        //equal to Cid
//                        Tree[rank[ProID]].FN[cidH]=false;//FN[cidH] should be assigned false only then the distance is shorter than corresponding vert
                        Tree[rank[ProID]].cnt[cidH]-=1;

                        //lower than Cid
                        for(int i=cidH+1;i<Tree[rank[ProID]].dis.size();i++){
                            if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[line[i]]].dis[cidH]){
                                Tree[rank[ProID]].cnt[i]-=1;
                            }
                        }
                    }

                }


                /// Step 2: compute the new value of super edges (shortcuts), i.e., e(ProID,Cid)
                //	cout<<Cw<<" increase to ";
                Cw=INF; int countwt=0;

                for(int i=0;i<Neighbors[ProID].size();i++){
                    if(Neighbors[ProID][i].first==Cid){
                        Cw=Neighbors[ProID][i].second;//the weight value in the original graph
                        countwt=1;
                        break;
                    }
                }

                int ssw=INF,wtt=INF,wid=-1;
                vector<pair<int,int>> Wnodes; //Wnodes.clear(); //Supportive vertices

                /// calculate the correct value by traversing supportive vertices
                assert(SCconNodesMT[ProID][Cid].size() == SCconNodesMT[Cid][ProID].size());
                Wnodes=SCconNodesMT[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;

                for(int i=0;i<Wnodes.size();i++){
                    wid=Wnodes[i].first;//wid must be periphery vertex
                    assert(CoreTag[wid] != -1);
                    if(CoreTag[wid] != pid)
                        continue;

                    for(int j=0;j<Tree[rank[wid]].vert.size();j++){
                        if(Tree[rank[wid]].vert[j].first==ProID){
                            ssw=Tree[rank[wid]].vert[j].second.first;
                        }
                        if(Tree[rank[wid]].vert[j].first==Cid){
                            wtt=Tree[rank[wid]].vert[j].second.first;
                        }
                    }
                    assert(ssw != INF);
                    assert(wtt != INF);
                    if(ssw+wtt<Cw){
                        Cw=ssw+wtt;
                        countwt=1;
                    }else if(ssw+wtt==Cw){
                        countwt+=1;
                    }
                }
                assert(Cw < INF);
                //cout<<Cw<<endl;
                //refresh the shortcut to the new value
                for(int i=0;i<Tree[rank[ProID]].vert.size();i++){
                    if(Tree[rank[ProID]].vert[i].first==Cid){
                        Tree[rank[ProID]].vert[i].second.first=Cw;
                        Tree[rank[ProID]].vert[i].second.second=countwt;
                        break;
                    }
                }
//                cout<<"ProH: "<<ProH<<" ; ProID: "<<ProID<<" ; Cid: "<<Cid<<" ; OldW: "<<OCdis[make_pair(ProID,Cid)]<<" ; NewW: "<<Cw<<endl;

            }

            influence = true;///
            if(influence){
                ProBeginID=ProID;
                ProBeginIDValue=true;
            }

            ProH-=1;
            ProID=Tree[Tree[rank[ProID]].pa].uniqueVertex;
        }
//        cout<<"MinH: "<<MinH<<" ; ProBeginID: "<<ProBeginID<<" ; ProBeginIDH: "<<Tree[rank[ProBeginID]].height<<endl;
        /// Step 3: update the shortcuts among interface vertices
        for(auto it=SECoreForUpdate.begin();it!=SECoreForUpdate.end();++it){
            lid=it->first; hid=it->second;
//            cout<<"Super edge: "<<lid<<" "<<hid<<endl;
            assert(SuppPartiID[lid].find(hid)!=SuppPartiID[lid].end());
            assert(!SuppPartiID[lid][hid].empty());

//            if(lid == 142488 || hid == 143850){
//                cout<<lid<<" "<<hid<<endl;
//            }
            int Cw=INF; int countwt=0;
            for(int i=0;i<Neighbors[lid].size();i++){
                if(Neighbors[lid][i].first==hid){
                    Cw=Neighbors[lid][i].second;//the weight value in the original graph
                    countwt=1;
                    break;
                }
            }
            int ssw=INF,wtt=INF,wid=-1;
            vector<pair<int,int>> Wnodes; //Wnodes.clear(); //Supportive vertices

            /// calculate the correct value by traversing supportive vertices
            assert(SCconNodesMT[lid][hid].size() == SCconNodesMT[hid][lid].size());
            Wnodes=SCconNodesMT[lid][hid]; //cout<<"wid num "<<Wnodes.size()<<endl;

            for(int i=0;i<Wnodes.size();i++){
                wid=Wnodes[i].first;//wid must be periphery vertex
                assert(CoreTag[wid] != -1);

                if(CoreTag[wid] != pid)
                    continue;

                for(int j=0;j<Tree[rank[wid]].vert.size();j++){
                    if(Tree[rank[wid]].vert[j].first==lid){//shortcut weight from wid to lid
                        ssw=Tree[rank[wid]].vert[j].second.first;
                    }
                    if(Tree[rank[wid]].vert[j].first==hid){//shortcut weight from wid to hid
                        wtt=Tree[rank[wid]].vert[j].second.first;
                    }
                }
                assert(ssw != INF);
                assert(wtt != INF);
                if(ssw+wtt<Cw){
                    Cw=ssw+wtt;
                    countwt=1;
                }else if(ssw+wtt==Cw){
                    countwt+=1;
                }
            }

//            cout<<"Super edge: "<<lid<<" "<<hid<<" "<<AdjaCoreMap[lid][hid]<<" "<<Cw<<endl;

            //assign new super edge weight
            assert(Cw < INF);
            if(SuppPartiIDReal[lid][hid].second.size() == 1){//if this super edge is only supported by this partition
                assert(*SuppPartiIDReal[lid][hid].second.begin()==pid);
                int tempDis=INF;
                set<int> tempPids;
                for(auto it=SuppPartiID[lid][hid].begin();it!=SuppPartiID[lid][hid].end();++it){
                    if(it->first == pid)
                        continue;
                    if(tempDis > it->second){
                        tempDis = it->second;
                        tempPids.clear();
                        tempPids.insert(it->first);
                    }else if(tempDis == it->second){
                        tempPids.insert(it->first);
                    }
                }
                //refresh the shortcut to the new value
                if(tempDis <= Cw){
                    AdjaCoreMap[lid][hid]=tempDis;
                    AdjaCoreMap[hid][lid]=tempDis;

                    if(tempDis == Cw){
                        tempPids.insert(pid);
                    }
                    SuppPartiIDReal[lid][hid].first = tempDis;
                    SuppPartiIDReal[hid][lid].first = tempDis;
                    SuppPartiIDReal[lid][hid].second = tempPids;
                    SuppPartiIDReal[hid][lid].second = tempPids;
                }else{//if tempDis > Cw
                    if(Cw > SuppPartiID[lid][hid][pid]){
                        AdjaCoreMap[lid][hid]=Cw;
                        AdjaCoreMap[hid][lid]=Cw;

//                        tempPids.clear();
//                        tempPids.insert(pid);
                        SuppPartiIDReal[lid][hid].first = Cw;
                        SuppPartiIDReal[hid][lid].first = Cw;
//                        SuppPartiIDReal[lid][hid].second = tempPids;
//                        SuppPartiIDReal[hid][lid].second = tempPids;
                    }
                }
            }
            else if(SuppPartiIDReal[lid][hid].second.size() > 1){//if this super edge is supported by more than one partition
                if(SuppPartiIDReal[lid][hid].second.find(pid) != SuppPartiIDReal[lid][hid].second.end()){//if found
                    if(Cw > SuppPartiIDReal[lid][hid].first){
                        SuppPartiIDReal[lid][hid].first = Cw;
                        SuppPartiIDReal[hid][lid].first = Cw;
                        SuppPartiIDReal[lid][hid].second.erase(pid);
                        SuppPartiIDReal[hid][lid].second.erase(pid);
                    }
                }
            }
            SuppPartiID[lid][hid][pid] = Cw;
            SuppPartiID[hid][lid][pid] = Cw;

            ///identify the affected partitions
            if(SuppPartiID[lid][hid].size()>1){
                for(auto it=SuppPartiID[lid][hid].begin();it!=SuppPartiID[lid][hid].end();++it){
                    if(it->first != pid){
                        affectedParti.insert(it->first);
//                        cout<<"affected parti: "<<it->first<<endl;
                    }
                }
            }

        }

        /// Top-down label update
        vector<int> interface;
        for(auto it=Tree[rank[rootVertex]].vert.begin();it!=Tree[rank[rootVertex]].vert.end();++it){
            interface.emplace_back(it->first);
        }
        //we cannot update interface entries simply by ProBeginIDInf as we do not know if the vertex of other branches are affected or not
//        if(!ProBeginIDInf.empty()){
//            int minH = MinHInf;
//            if(minH <= MinH){
//                for(auto it=ProBeginIDInf.begin();it!=ProBeginIDInf.end();++it){
//                    int rid=*it;
//                    InterfacePropagate( rank[rid],interface,Tree,true);
//                }
//                flagHigher = true;
//            }else{
//                InterfacePropagate( rank[ProBeginID],interface,Tree,true);
//            }
//        }else{
//            InterfacePropagate( rank[ProBeginID],interface,Tree,true);
//        }

        PartiUpdateExt[pid] = true;
        if(affectedParti.size()>=1){
            if(ifParallel){//use multi-thread ifParallel
//                cout<<"Update the interface entry with multiple threads."<<endl;
                boost::thread_group thread;
                vector<int> vParti;
                vParti.push_back(pid);
                for(auto it=affectedParti.begin();it!=affectedParti.end();++it){
                    vParti.push_back(*it);
                    PartiUpdateExt[*it] = true;
                }
                cout<<"Affected partition number: "<<vParti.size()<<endl;
                if(vParti.size()>threadnum){
//                    cout<<"entry 1"<<endl;
                    int step=vParti.size()/threadnum;
                    boost::thread_group threadf;
                    for(int i=0;i<threadnum;i++){
                        pair<int,int> p;
                        p.first=i*step;
                        if(i==threadnum-1)
                            p.second=vParti.size();
                        else
                            p.second=(i+1)*step;
                        threadf.add_thread(new boost::thread(&Graph::InterfacePropagateParallel, this, p, boost::ref(vParti), true));
                    }
                    threadf.join_all();
                }else{
//                    cout<<"entry 2"<<endl;
                    boost::thread_group threadf;
                    for(int i=0;i<vParti.size();i++){
                        threadf.add_thread(new boost::thread(&Graph::InterfacePropagateParallel, this, make_pair(i,i+1), boost::ref(vParti), true));
                    }
                    threadf.join_all();
                }


            }else{//use single-thread
//                cout<<"Update the interface entry with single thread."<<endl;
                InterfacePropagate( rank[rootVertex], interface,Tree,true);//propagate from root vertex

                if(!affectedParti.empty()){//if we need to update the shortcut of other partition
                    for(auto it=affectedParti.begin();it!=affectedParti.end();++it){
                        int Pid = *it;
                        int rootID = BoundVertex[Pid][BoundVertex[Pid].size()-1];
                        PartiUpdateExt[Pid] = true;
                        vector<int> interfaceP;
                        for(auto it2=Tree[rank[rootID]].vert.begin();it2!=Tree[rank[rootID]].vert.end();++it2){
                            interfaceP.emplace_back(it2->first);
                        }
                        InterfacePropagate( rank[rootID], interfaceP,Tree,true);
                    }
                }
            }
        }else{//use single-thread
//            cout<<"Update the interface entry with single thread."<<endl;
            InterfacePropagate( rank[rootVertex], interface,Tree,true);//propagate from root vertex

            if(!affectedParti.empty()){//if we need to update the shortcut of other partition
                for(auto it=affectedParti.begin();it!=affectedParti.end();++it){
                    int Pid = *it;
                    int rootID = BoundVertex[Pid][BoundVertex[Pid].size()-1];
                    PartiUpdateExt[Pid] = true;
                    vector<int> interfaceP;
                    for(auto it2=Tree[rank[rootID]].vert.begin();it2!=Tree[rank[rootID]].vert.end();++it2){
                        interfaceP.emplace_back(it2->first);
                    }
                    InterfacePropagate( rank[rootID], interfaceP,Tree,true);
                }
            }
        }



//        cout<<"Ancestors: ";
//        for(int i=MinH-1;i<Tree[rank[210695]].height;++i){
//            cout<<Tree[rank[210695]].vAncestor[i]<<" ";
//        }
//        cout<<endl;

        if(ProBeginIDValue){
            vector<int> ancestors; //line1.clear();
            ancestors.reserve(heightMax);
            pachid=Tree[Tree[rank[ProBeginID]].pa].uniqueVertex;
            while(pachid != -1){//Tree[rank[pachid]].height>1
                ancestors.insert(ancestors.begin(),pachid);
                pachid=Tree[Tree[rank[pachid]].pa].uniqueVertex;
            }
            ancestors.insert(ancestors.begin(),pachid);

//            for(int i=0;i<Tree[rank[210695]].vAncestor.size();++i){
//                if(Tree[rank[210695]].vAncestor[i] == ProBeginID){
//                    cout<<"Found! "<<i<<" "<<ProBeginID<<endl;
//                }
//            }

            eachNodeProcessIncrease1New(rank[ProBeginID],ancestors,interface,ChangeNum,Tree,rank,VidtoTNid,lowestH);
            //cout<<"9999999999999"<<endl;
        }
    }
    //return ChangeNum;
}

void Graph::eachNodeProcessIncrease1New(int child, vector<int>& line, vector<int>& interfaces, int& changelabel, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid, int lowestH){
    int childID=Tree[child].uniqueVertex;
    int childH=Tree[child].height-1;

//    if(child == rank[210695])
//        cout<<210695<<endl;

    /// update ancestor entries
    for(int di=1;di<Tree[child].dis.size()-1;di++){//omit the last distance, i.e., the distance to itself
//        if(child == rank[210695] && line[di] == 208312){
//            cout<<210695<<" "<<208312<<endl;
//        }

        if(di>=line.size())
            cout<<"Problem here! "<<di<<" "<<line.size()<<endl;
        if(rank[line[di]]==0)
            cout<<"Wrong! "<<line[di]<<" "<<rank[line[di]]<<" "<<CoreTag[line[di]]<<" "<<int(BoundTag[line[di]].first)<<endl;
//        if(Tree[child].cnt[di]==0 || Tree[child].height>=lowestH){//check the changed distance, if cnt[di]==0 which indicates the ancestor entry from child to di has changed
        if(true){//check the changed distance
            changelabel+=1;
            //firstly, check which dis can be infected
            int disBF=Tree[child].dis[di];//the changed ancestor distance
            int PID;

            //secondly, calculate the actual distance
            int dis=INF; int count=0;
            int Dvb; int b,bH; int DDvb=INF;
            for(int j=0;j<Tree[child].vert.size();j++){
                Dvb=Tree[child].vert[j].second.first;
                b=Tree[child].vert[j].first;

                if(CoreTag[b] == -1){//if it is interface vertex, check whether the distance label may be affected by interface
                    int vbW=Tree[child].vert[j].second.first;

                    if(dis>vbW+Tree[rank[line[di]]].disInf[b]){
                        dis=vbW+Tree[rank[line[di]]].disInf[b];
                        count=1;
                    }else if(dis==vbW+Tree[rank[line[di]]].disInf[b]){
                        count+=1;
                    }
//                    continue;
                }else{//if it is ancestor vertex
                    bH=Tree[rank[b]].height-1;
                    if(bH<di){//if the neighbor b has higher order than line[di]
                        if(Dvb+Tree[rank[line[di]]].dis[bH]<dis){
                            dis=Dvb+Tree[rank[line[di]]].dis[bH];
                            count=1;
                        }else if(Dvb+Tree[rank[line[di]]].dis[bH]==dis){
                            count+=1;
                        }
                    }else if(bH==di){//if the neighbor b is line[di]
                        DDvb=Dvb;
                        if(Dvb<dis){
                            dis=Dvb;
                            count=1;
                        }else if(Dvb==dis){
                            count+=1;
                        }
                    }else{//if the neighbor b has lower order than line[di]
                        if(Dvb+Tree[rank[b]].dis[di]<dis){
                            dis=Dvb+Tree[rank[b]].dis[di];
                            count=1;
                        }else if(Dvb+Tree[rank[b]].dis[di]==dis){
                            count+=1;
                        }
                    }
                }

            }
            if(DDvb==dis)
                Tree[child].FN[di]=true;
            else{
                Tree[child].FN[di]=false;
            }
            bool flagCheck=false;
            if(Tree[child].dis[di]<dis){
                Tree[child].dis[di]=dis;
                Tree[child].cnt[di]=count;
                flagCheck = true;
            }else if(Tree[child].dis[di]==dis && Tree[child].cnt[di]!=count){
                Tree[child].cnt[di]=count;
//                cout<<"Special case!!"<<endl;
            }

            if(flagCheck){
                //chidlID, check the children of childID,
                for(int k=0;k<VidtoTNid[childID].size();k++){//check the children who have direct vert shortcut to childID
                    PID=VidtoTNid[childID][k];
//                if(Tree[PID].FN[childH] && Tree[PID].dis[di]==disBF+Tree[PID].dis[childH]){//?? if the ancestor distance from PID to childID sources from vert
                    if(Tree[PID].dis[di]==disBF+Tree[PID].dis[childH]){//
                        Tree[PID].cnt[di]-=1;
                    }
                }

                //line[i], check the children of line[di]
                for(int k=0;k<VidtoTNid[line[di]].size();k++){//check the ancestor's child tree nodes
                    PID=VidtoTNid[line[di]][k];
                    if(PID>child){//? if child has higher order, i.e., the children have lower order
                        assert(NodeOrder[Tree[PID].uniqueVertex] < NodeOrder[childID]);
//                    if(Tree[PID].FN[di] && Tree[PID].dis[childH]==disBF+Tree[PID].dis[di]){//?? if the ancestor entry to line[di] is equal to vert and the ancestor distance to child is equal to
                        if(Tree[PID].dis[childH]==disBF+Tree[PID].dis[di]){//
                            Tree[PID].cnt[childH]-=1;
                        }
                    }
                }
            }

        }
    }

    line.push_back(childID);
    if(rank[childID]==0)
        cout<<"Wrong!2 "<<childID<<" "<<childH<<" "<<rank[childID]<<" "<<CoreTag[childID]<<" "<<int(BoundTag[childID].first)<<endl;
    for(int i=0;i<Tree[child].ch.size();i++){
        eachNodeProcessIncrease1New(Tree[child].ch[i],line,interfaces, changelabel,Tree,rank,VidtoTNid,lowestH);
    }
    line.pop_back();
}

//Query within one partition: No LCA version
/*int Graph::QueryPeripheryTree(int ID1, int ID2, int PID){
    int d=INF;

    int r1=rank[ID1], r2=rank[ID2];

    if(AdjaGraph[ID1].size()!=0 && AdjaGraph[ID2].size()!=0){//consider the vertex connectivity in subgraph
        if(ID1==ID2) return 0;

        unordered_map<int,int>::iterator it;
        int hub, dis1, dis2;
        int hubfinal,dis1final,dis2final;

        int height1 = Tree[r1].dis.size();
        int height2 = Tree[r2].dis.size();
//        cout<<"Tree depth: "<<height1<< " "<<height2<<endl;
        int temp_dis;
        if(height1 > height2){
            for(int i=1;i<height2;++i){
                if(Tree[r1].vAncestor[i] == Tree[r2].vAncestor[i]){//if common ancestor
                    temp_dis = Tree[r1].dis[i] + Tree[r2].dis[i];
                    if(d > temp_dis){
                        d = temp_dis;
                    }
                }
                else{
                    break;
                }
            }
        }else{
            for(int i=1;i<height1;++i){
                if(Tree[r1].vAncestor[i] == Tree[r2].vAncestor[i]){//if common ancestor
                    temp_dis = Tree[r1].dis[i] + Tree[r2].dis[i];
                    if(d > temp_dis){
                        d = temp_dis;
                    }
                }else{
                    break;
                }
            }
        }
        if(d == INF){
            cout<<"Wrong! d = "<<INF<<" "<<AdjaGraph[ID1].size()<<" "<<AdjaGraph[ID2].size()<<endl;
            for(int i=0;i<Tree[r1].vAncestor.size();++i){
                cout<<Tree[r1].vAncestor[i]<<"\t";
            }
            cout<<endl;
            for(int i=0;i<Tree[r2].vAncestor.size();++i){
                cout<<Tree[r2].vAncestor[i]<<"\t";
            }
            cout<<endl;
        }
    }

    return d;
}*/

//Query within one partition: LCA version
int Graph::QueryPeripheryTree(int ID1, int ID2, int PID){
//    assert(AdjaGraph[ID1].size()!=0 && AdjaGraph[ID2].size()!=0);
	if((CoreTag[ID1]==-1 && BoundTag[ID1].first) || (CoreTag[ID2]==-1 && BoundTag[ID2].first)){//if either endpoint is interface vertex
		if((CoreTag[ID1]==-1 && BoundTag[ID1].first) && (CoreTag[ID1]==-1 && BoundTag[ID1].first)){//if both endpoints are interface vertex
			int dis = INF;
			for(auto it=AdjaCore[ID1].begin();it!=AdjaCore[ID1].end();++it){//get the distance from contracted graph, AdjaCore
				if(it->first == ID2){
					dis = it->second;
					break;
				}
			}
			return dis;
		}else{//if only one endpoint is interface vertex
			if((CoreTag[ID1]==-1 && BoundTag[ID1].first)){//if ID1 is interface vertex
				if(Tree[rank[ID2]].disInf.find(ID1) != Tree[rank[ID2]].disInf.end()){
					return Tree[rank[ID2]].disInf[ID1];
				}else{
					cout<<"Wrong interface entry!"<<endl; exit(1);
				}
			}else if((CoreTag[ID2]==-1 && BoundTag[ID2].first)){//if ID2 is interface vertex
				if(Tree[rank[ID1]].disInf.find(ID2) != Tree[rank[ID1]].disInf.end()){
					return Tree[rank[ID1]].disInf[ID2];
				}else{
					cout<<"Wrong interface entry!"<<endl; exit(1);
				}
			}
		}
	}else if(CoreTag[ID1] != -1 && CoreTag[ID2] != -1){//if both endpoint are periphery vertex, use the ancestor entry
		if(ID1==ID2)
			return 0;
		int r1=rank[ID1], r2=rank[ID2];
		if(r1<=0 ||r2<=0){
			cout<<"rank error: "<<r1<<" "<<r2<<endl;
			exit(1);
		}
		int LCA=LCAQueryPartition(r1,r2,PID);
        int d1,d2;
		if(LCA==r1)
			return Tree[r2].dis[Tree[r1].pos.back()];
		else if(LCA==r2)
			return Tree[r1].dis[Tree[r2].pos.back()];
		else{
			int tmp=INF;
			for(int i=0;i<Tree[LCA].pos.size();i++){
				if(Tree[LCA].pos[i] == -1)
					continue;
				if(tmp>Tree[r1].dis[Tree[LCA].pos[i]]+Tree[r2].dis[Tree[LCA].pos[i]]){
                    tmp=Tree[r1].dis[Tree[LCA].pos[i]]+Tree[r2].dis[Tree[LCA].pos[i]];
                    d1=Tree[r1].dis[Tree[LCA].pos[i]];
                    d2=Tree[r2].dis[Tree[LCA].pos[i]];
                }

			}
//            cout<<"LCA: "<<LCA<<" ; d1: "<<d1<<" ; d2: "<<d2<<endl;
			return tmp;
		}
	}else{
		cout<<"Wrong query!"<<endl; exit(1);
	}



    return INF;
}

int Graph::LCAQueryPartition(int _p, int _q, int PID){//space complexity and preprocessing complexity: O(nlog(n))
	int p = toRMQ[_p], q = toRMQ[_q];
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
	if (Tree[RMQIndex[k][p]].height < Tree[RMQIndex[k][q]].height)
		return RMQIndex[k][p];
	else return RMQIndex[k][q];
}

/*int Graph::LCAQueryPartition(int _p, int _q, int PID){//space complexity and preprocessing complexity: O(nlog(n))
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
}*/
