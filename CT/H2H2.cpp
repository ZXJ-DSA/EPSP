/*
 * H2H.cpp
 *
 *  Created on: 14 Oct 2022
 *      Author: zhangmengxuan
 */
#include "head2.h"

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
//parallel version of H2Hindex
void Graph::H2HindexParallel(bool ifParallel){

    if(ifParallel){//if parallel in H2H index construction
        cout<<"Multi-thread in H2Hindex construction while single-thread in CHindex construction."<<endl;
        //multiple thread
        if(partiNum>threadnum){
            int step=partiNum/threadnum;
            cout<<"step: "<<step<<endl;
            boost::thread_group threadf;
            for(int i=0;i<threadnum;i++){
//                cout<<"threadnum "<<i<<endl;
                pair<int,int> p;
                p.first=i*step;
                if(i==threadnum-1)
                    p.second=partiNum;
                else
                    p.second=(i+1)*step;
                threadf.add_thread(new boost::thread(&Graph::H2HindexP, this, p, false));
            }
            threadf.join_all();
        }else{
            boost::thread_group threadf;
            for(int pid=0;pid<partiNum;++pid) {
                threadf.add_thread(new boost::thread(&Graph::H2HindexP, this, make_pair(pid,pid+1), false));
            }
            threadf.join_all();
        }
    }else{
        cout<<"Single-thread in H2Hindex construction while multi-thread in CHindex construction."<<endl;
        //single thread
        for(int pid=0;pid<partiNum;++pid) {
            H2HindexP(make_pair(pid,pid+1),true);
        }
    }

}

void Graph::H2HindexP(pair<int, int> p, bool ifParallel) {
    int ID;
    for(int partiID=p.first;partiID<p.second;++partiID){
        Graph h;
        h.nodenum=nodenum;
        h.NodeOrder=NodeOrder;
        h.vNodeOrder=vNodeOrder;
//        h.Neighbors=g.AdjaParti[partiID];
        h.Neighbors.assign(nodenum, vector<pair<int,int>>());
        for(auto it=AdjaPartiM[partiID].begin();it!=AdjaPartiM[partiID].end();++it){
            ID = it->first;
            if(it->second.empty()){
                cout<<"Wrong! Empty neighbor!"<<endl;
                exit(1);
            }
            h.Neighbors[ID] = it->second;
            Semaphore* s = new Semaphore(1);
            h.mSm[ID]=s;
        }

        h.H2Hindex(ifParallel);
        if(partiID >= Trees.size()){
            cout<<"Overflow!"<<endl;
        }
        Trees[partiID] = h.Tree;
        ranks[partiID] = h.rank;
        heightMaxs[partiID]=h.heightMax;
        //SCconNodesMTs[partiID]=h.SCconNodesMT;
        //VidtoTNids[partiID]=h.VidtoTNid;
    }
}

void Graph::H2Hindex(bool ifParallel){
	vector<pair<int,pair<int,int>>> vect;
	NeighborCon.assign(nodenum,vect);

    if(ifParallel){
        CHindexMT();//multiple threads
    }else{
        CHindex();//single threads
    }

	makeTree();
	makeIndex();
}


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
                cout<<"Multi thread step: "<<step<<endl;
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
}

void Graph::NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x){
	//sm->wait();
	int ID1, w1;
	int ID2, w2;
	for(int k=p.first;k<p.second;k++){
		ID1=Neighvec[k].first;
		w1=Neighvec[k].second.first;
		for(int h=0;h<Neighvec.size();h++){
            mSm[ID1]->wait();
			ID2=Neighvec[h].first;
			w2=Neighvec[h].second.first;
			insertEMTorder(ID1, ID2, w1+w2);
			if(ID1<ID2)
				SCconNodesMT[ID1][ID2].push_back(x);
            mSm[ID1]->notify();
		}
	}
	//sm->notify();
}

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

void Graph::deleteEorder(int u,int v){
	if(E[u].find(v)!=E[u].end()){
		E[u].erase(E[u].find(v));
	}

	if(E[v].find(u)!=E[v].end()){
		E[v].erase(E[v].find(u));
	}
}

void Graph::makeTree(){
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
//        if(len==0)
//            break;
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
}

int Graph::match(int x,vector<pair<int,pair<int,int>>> &vert){
	int nearest=vert[0].first;
	for(int i=1;i<vert.size();i++){
		if(rank[vert[i].first]>rank[nearest])//get the vertex with highest rank
			nearest=vert[i].first;
	}
	int p=rank[nearest];
	return p;
}

void Graph::makeIndex(){
	//makeRMQ();

	//initialize
	vector<int> list; //list.clear();
	list.push_back(Tree[0].uniqueVertex);
	Tree[0].pos.clear();
	Tree[0].pos.push_back(0);
    Tree[0].dis.push_back(0);
    Tree[0].vAncestor.push_back(Tree[0].uniqueVertex);

	for(int i=0;i<Tree[0].ch.size();i++){
		makeIndexDFS(Tree[0].ch[i],list);
	}
}
///?
void Graph::makeRMQ(){
	//EulerSeq.clear();
	toRMQ.assign(nodenum,0);
	//RMQIndex.clear();
	makeRMQDFS(0, 1);
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

void Graph::makeRMQDFS(int p, int height){
	toRMQ[p] = EulerSeq.size();
	EulerSeq.push_back(p);
	for (int i = 0; i < Tree[p].ch.size(); i++){
		makeRMQDFS(Tree[p].ch[i], height + 1);
		EulerSeq.push_back(p);
	}
}

void Graph::makeIndexDFS(int p, vector<int>& list){
	//initialize
	int NeiNum=Tree[p].vert.size();
	Tree[p].pos.assign(NeiNum+1,0);
	Tree[p].dis.assign(list.size(),INF);
	Tree[p].cnt.assign(list.size(),0);
	Tree[p].FN.assign(list.size(),true);
    Tree[p].vAncestor.assign(list.begin(), list.end());
    Tree[p].vAncestor.push_back(Tree[p].uniqueVertex);
    Tree[p].dis.push_back(0);

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
}

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

//H2H index update
/*void Graph::DecreaseH2H(int a,int b, int newW, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax){
	map<int,int> checkedDis;//map<tree node ID, distance index>
//cout<<"1111111111"<<endl;
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
//cout<<"2222222222"<<endl;
	for(int i=0;i<Tree.size();i++){
		Tree[i].DisRe.clear();
	}

	int lid,hid;
	if(NodeOrder[a]<NodeOrder[b]){
		lid=a;hid=b;
	}else{
		lid=b;hid=a;
	}

	int IniH=Tree[rank[lid]].height;//the height where weight change begins
	int ProH=Tree[rank[lid]].height; int ProID=lid;
	vector<set<int>> SCre;//record the shortcut change in each height
	set<int> ss;//ss.clear();
	SCre.assign(ProH+1,ss);

	//map<int,set<int>> DisRe;//rankid; record the distance change caused by the shortcut in each height
	//DisRe.clear();

	int MinH;

	bool tri=false;
	for(int i=0;i<Tree[rank[lid]].vert.size();i++){
		if(Tree[rank[lid]].vert[i].first==hid){
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

	set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed

	vector<int> ProIDRecord; ProIDRecord.assign(ProH+1,0);
//cout<<"33333333333"<<endl;
	//int ProBeginH;
	int ProBeginID;
	bool ProBeginIDValue=false;
	if(tri){
	//cout<<"Bottom-up ;;;;;;;;;;;;;;;;;; "<<endl;
//cout<<"4444444444444"<<endl;
		while(ProH>=MinH){

			ProIDRecord[ProH]=ProID;
			vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
			bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
			for(auto it=SCre[ProH].begin();it!=SCre[ProH].end();it++){
				int Cid=*it; int Cw;//=OCdis[make_pair(ProID,Cid)];
				int cidH=Tree[rank[Cid]].height-1;

				map<int,int> Hnei; //Hnei.clear();
				vector<pair<int,int>> Lnei; //Lnei.clear();
				for(int j=0;j<Vert.size();j++){
					if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
						Hnei[Vert[j].first]=Vert[j].second.first;
					}else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
						Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
					}else{
						Cw=Vert[j].second.first;
					}
				}

				if(Tree[rank[ProID]].dis[cidH]>=Cw){
					Tree[rank[ProID]].dis[cidH]=Cw;
					Tree[rank[ProID]].FN[cidH]=true;
					ProIDdisCha=true;
					Tree[rank[ProID]].DisRe.insert(Cid);
					//DisRe[rank[ProID]].insert(Cid); //cout<<"dischange Cid "<<Cid<<endl;
				}

				int hid,hidHeight,lid,lidHeight,wsum;
				for(int j=0;j<Tree[rank[Cid]].vert.size();j++){
					hid=Tree[rank[Cid]].vert[j].first;hidHeight=Tree[rank[hid]].height-1;
					if(Hnei.find(hid)!=Hnei.end()){
						wsum=Cw+Hnei[hid];
						if(wsum<Tree[rank[Cid]].vert[j].second.first){
							Tree[rank[Cid]].vert[j].second.first=wsum;
							Tree[rank[Cid]].vert[j].second.second=1;
							SCre[Tree[rank[Cid]].height].insert(hid);
							if(Tree[rank[Cid]].height<MinH) MinH=Tree[rank[Cid]].height;

						}else if(wsum==Tree[rank[Cid]].vert[j].second.first){
							Tree[rank[Cid]].vert[j].second.second+=1;
						}

					}
				}
				for(int j=0;j<Lnei.size();j++){
					lid=Lnei[j].first;lidHeight=Tree[rank[lid]].height-1;
					for(int k=0;k<Tree[rank[lid]].vert.size();k++){
						if(Tree[rank[lid]].vert[k].first==Cid){
							wsum=Cw+Lnei[j].second;
							if(Tree[rank[lid]].vert[k].second.first>wsum){
								Tree[rank[lid]].vert[k].second.first=wsum;
								Tree[rank[lid]].vert[k].second.second=1;
								SCre[Tree[rank[lid]].height].insert(Cid);
								if(Tree[rank[lid]].height<MinH) MinH=Tree[rank[lid]].height;

							}else if(Tree[rank[lid]].vert[k].second.first==wsum){
								Tree[rank[lid]].vert[k].second.second+=1;
							}

							break;
						}
					}
				}
			}

			if(ProIDdisCha){//if the distance labeling is detected changed
				vertexIDChL.insert(ProID);
				//ProBeginH=ProH;
				ProBeginID=ProID;
				ProBeginIDValue=true;
				//cout<<"PorBeginID "<<ProBeginID<<endl;
			}

			ProH-=1;
			ProID=Tree[Tree[rank[ProID]].pa].uniqueVertex;
		}

		if(ProBeginIDValue){
			vector<int> linee; //linee.clear();
			linee.reserve(heightMax);

			int pachidd=Tree[Tree[rank[ProBeginID]].pa].uniqueVertex;
		//cout<<"aaaaaaaaaaaa "<<pachidd<<" "<<Tree[rank[pachidd]].height<<endl;
			while(Tree[rank[pachidd]].height>1){
				linee.insert(linee.begin(),pachidd);
				pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
				//cout<<"height "<<Tree[rank[pachidd]].height<<endl;
			}

			linee.insert(linee.begin(),pachidd);

			//top-down process

			EachNodeProBDis5(rank[ProBeginID], linee, vertexIDChL, checkedDis, Tree, rank);
		}
	}
	//return checkedDis.size();
}*/

//H2H index update: Map version
void Graph::DecreaseH2H(int a,int b, int newW, unordered_map<int,vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax){
    map<int,int> checkedDis;//map<tree node ID, distance index>
//cout<<"1111111111"<<endl;
    if(Neighbors.find(a) == Neighbors.end() || Neighbors.find(b) == Neighbors.end()){//if not found a, b
        cout<<"Wrong! There is no a and b in Neighbors!!"<<endl;
        exit(1);
    }
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
//cout<<"2222222222"<<endl;
    for(int i=0;i<Tree.size();i++){
        Tree[i].DisRe.clear();
    }

    int lid,hid;
    if(NodeOrder[a]<NodeOrder[b]){
        lid=a;hid=b;
    }else{
        lid=b;hid=a;
    }

    int IniH=Tree[rank[lid]].height;//the height where weight change begins
    int ProH=Tree[rank[lid]].height; int ProID=lid;
    vector<set<int>> SCre;//record the shortcut change in each height
    set<int> ss;//ss.clear();
    SCre.assign(ProH+1,ss);

    //map<int,set<int>> DisRe;//rankid; record the distance change caused by the shortcut in each height
    //DisRe.clear();

    int MinH;

    bool tri=false;
    for(int i=0;i<Tree[rank[lid]].vert.size();i++){
        if(Tree[rank[lid]].vert[i].first==hid){
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

    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed

    vector<int> ProIDRecord; ProIDRecord.assign(ProH+1,0);
//cout<<"33333333333"<<endl;
    //int ProBeginH;
    int ProBeginID;
    bool ProBeginIDValue=false;
    if(tri){
        //cout<<"Bottom-up ;;;;;;;;;;;;;;;;;; "<<endl;
//cout<<"4444444444444"<<endl;
        while(ProH>=MinH){

            ProIDRecord[ProH]=ProID;
            vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
            bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
            for(auto it=SCre[ProH].begin();it!=SCre[ProH].end();it++){
                int Cid=*it; int Cw;//=OCdis[make_pair(ProID,Cid)];
                int cidH=Tree[rank[Cid]].height-1;

                map<int,int> Hnei; //Hnei.clear();
                vector<pair<int,int>> Lnei; //Lnei.clear();
                for(int j=0;j<Vert.size();j++){
                    if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                        Hnei[Vert[j].first]=Vert[j].second.first;
                    }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                        Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                    }else{
                        Cw=Vert[j].second.first;
                    }
                }

                if(Tree[rank[ProID]].dis[cidH]>=Cw){
                    Tree[rank[ProID]].dis[cidH]=Cw;
                    Tree[rank[ProID]].FN[cidH]=true;
                    ProIDdisCha=true;
                    Tree[rank[ProID]].DisRe.insert(Cid);
                    //DisRe[rank[ProID]].insert(Cid); //cout<<"dischange Cid "<<Cid<<endl;
                }

                int hid,hidHeight,lid,lidHeight,wsum;
                for(int j=0;j<Tree[rank[Cid]].vert.size();j++){
                    hid=Tree[rank[Cid]].vert[j].first;hidHeight=Tree[rank[hid]].height-1;
                    if(Hnei.find(hid)!=Hnei.end()){
                        wsum=Cw+Hnei[hid];
                        if(wsum<Tree[rank[Cid]].vert[j].second.first){
                            Tree[rank[Cid]].vert[j].second.first=wsum;
                            Tree[rank[Cid]].vert[j].second.second=1;
                            SCre[Tree[rank[Cid]].height].insert(hid);
                            if(Tree[rank[Cid]].height<MinH) MinH=Tree[rank[Cid]].height;

                        }else if(wsum==Tree[rank[Cid]].vert[j].second.first){
                            Tree[rank[Cid]].vert[j].second.second+=1;
                        }

                    }
                }
                for(int j=0;j<Lnei.size();j++){
                    lid=Lnei[j].first;lidHeight=Tree[rank[lid]].height-1;
                    for(int k=0;k<Tree[rank[lid]].vert.size();k++){
                        if(Tree[rank[lid]].vert[k].first==Cid){
                            wsum=Cw+Lnei[j].second;
                            if(Tree[rank[lid]].vert[k].second.first>wsum){
                                Tree[rank[lid]].vert[k].second.first=wsum;
                                Tree[rank[lid]].vert[k].second.second=1;
                                SCre[Tree[rank[lid]].height].insert(Cid);
                                if(Tree[rank[lid]].height<MinH) MinH=Tree[rank[lid]].height;

                            }else if(Tree[rank[lid]].vert[k].second.first==wsum){
                                Tree[rank[lid]].vert[k].second.second+=1;
                            }

                            break;
                        }
                    }
                }
            }

            if(ProIDdisCha){//if the distance labeling is detected changed
                vertexIDChL.insert(ProID);
                //ProBeginH=ProH;
                ProBeginID=ProID;
                ProBeginIDValue=true;
                //cout<<"PorBeginID "<<ProBeginID<<endl;
            }

            ProH-=1;
            ProID=Tree[Tree[rank[ProID]].pa].uniqueVertex;
        }

        if(ProBeginIDValue){
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);

            int pachidd=Tree[Tree[rank[ProBeginID]].pa].uniqueVertex;
            //cout<<"aaaaaaaaaaaa "<<pachidd<<" "<<Tree[rank[pachidd]].height<<endl;
            while(Tree[rank[pachidd]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
                //cout<<"height "<<Tree[rank[pachidd]].height<<endl;
            }

            linee.insert(linee.begin(),pachidd);

            //top-down process

            EachNodeProBDis5(rank[ProBeginID], linee, vertexIDChL, checkedDis, Tree, rank);
        }
    }
    //return checkedDis.size();
}

void Graph::EachNodeProBDis5(int child,vector<int>& line,set<int>& vertexIDChL, map<int,int>& checkedDis, vector<Node> &Tree, vector<int> &rank){
	bool ProIDdisCha=false;

	if(Tree[child].DisRe.size()!=0){
		for(int k=0;k<Tree[child].vert.size();k++){
			int b=Tree[child].vert[k].first, bH=Tree[rank[b]].height-1,vbW=Tree[child].vert[k].second.first;
			if(Tree[child].FN[bH]){
				if(Tree[child].DisRe.find(b)!=Tree[child].DisRe.end()){//all ancestor check
					for(int i=0;i<bH;i++){
						checkedDis.insert(make_pair(child,i));
						if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
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

				}else{//partial ancestor check

					if(vertexIDChL.find(b)!=vertexIDChL.end()){
						for(int i=0;i<bH;i++){
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
	}else{
		for(int k=0;k<Tree[child].vert.size();k++){
			int b=Tree[child].vert[k].first, bH=Tree[rank[b]].height-1,vbW=Tree[child].vert[k].second.first;
			if(Tree[child].FN[bH]){
				if(vertexIDChL.find(b)!=vertexIDChL.end()){
					for(int i=0;i<bH;i++){
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

	if(ProIDdisCha){
		vertexIDChL.insert(Tree[child].uniqueVertex);
	}

	line.push_back(Tree[child].uniqueVertex);
	for(int i=0;i<Tree[child].ch.size();i++){
		EachNodeProBDis5(Tree[child].ch[i], line, vertexIDChL,checkedDis,Tree, rank);
	}
	line.pop_back();

}

/*void Graph::IncreaseH2H(int a,int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<int>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid){
	int ChangeNum=0;
	map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
		//OCdis.clear();
//cout<<"11111111111111"<<endl;
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
//cout<<"22222222222222"<<endl;
		int lid,hid;
		if(NodeOrder[a]<NodeOrder[b]){
			lid=a;hid=b;
		}else{
			lid=b;hid=a;
		}

		int IniH=Tree[rank[lid]].height;//the height where weight change begins
		int ProH=Tree[rank[lid]].height; int ProID=lid;
		vector<set<int>> SCre;//record the shortcut change in each height
		set<int> vec; //vec.clear();
		SCre.assign(IniH+1,vec);
		int MinH;

		vector<int> line; //line.clear();
		line.reserve(heightMax);
		int pachid=ProID;
		while(Tree[rank[pachid]].height>1){
			line.insert(line.begin(),pachid);
			pachid=Tree[Tree[rank[pachid]].pa].uniqueVertex;
		}
		line.insert(line.begin(),pachid);

		bool tri=false;
		for(int i=0;i<Tree[rank[lid]].vert.size();i++){
			if(Tree[rank[lid]].vert[i].first==hid){
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
		if(tri){
			//shortcut update
			while(ProH>=MinH){
				influence=false;
				vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
				for(auto it=SCre[ProH].begin();it!=SCre[ProH].end();it++){
					int Cid=*it; int Cw=OCdis[make_pair(ProID,Cid)];
					int cidH=Tree[rank[Cid]].height-1;

					map<int,int> Hnei; //Hnei.clear();
					vector<pair<int,int>> Lnei; //Lnei.clear();
					for(int j=0;j<Vert.size();j++){
						if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
							Hnei[Vert[j].first]=Vert[j].second.first;
						}else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
							Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
						}
					}
					//check the affected shortcuts
					int hid,lid;
					for(int j=0;j<Tree[rank[Cid]].vert.size();j++){
						hid=Tree[rank[Cid]].vert[j].first;
						if(Hnei.find(hid)!=Hnei.end()){
							if(Cw+Hnei[hid]==Tree[rank[Cid]].vert[j].second.first){
								Tree[rank[Cid]].vert[j].second.second-=1;
								if(Tree[rank[Cid]].vert[j].second.second<1){
									SCre[Tree[rank[Cid]].height].insert(hid);
									if(Tree[rank[Cid]].height<MinH) MinH=Tree[rank[Cid]].height;
									OCdis[make_pair(Cid,hid)]=Cw+Hnei[hid];
								}
							}
						}
					}
					for(int j=0;j<Lnei.size();j++){
						lid=Lnei[j].first;
						for(int k=0;k<Tree[rank[lid]].vert.size();k++){
							if(Tree[rank[lid]].vert[k].first==Cid){
								if(Tree[rank[lid]].vert[k].second.first==Cw+Lnei[j].second){
									Tree[rank[lid]].vert[k].second.second-=1;
									if(Tree[rank[lid]].vert[k].second.second<1){
										SCre[Tree[rank[lid]].height].insert(Cid);
										if(Tree[rank[lid]].height<MinH) MinH=Tree[rank[lid]].height;
										OCdis[make_pair(lid,Cid)]=Cw+Lnei[j].second;
									}
								}
								break;
							}
						}
					}

					//before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
					if(Tree[rank[ProID]].FN[cidH]){
						influence=true;
						//higher than Cid
						for(int i=0;i<cidH;i++){
							if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[Cid]].dis[i]){
								Tree[rank[ProID]].cnt[i]-=1;
							}
						}

						//equal to Cid
						Tree[rank[ProID]].FN[cidH]=false;
						Tree[rank[ProID]].cnt[cidH]-=1;

						//lower than Cid
						for(int i=cidH+1;i<Tree[rank[ProID]].dis.size();i++){
							if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[line[i]]].dis[cidH]){
								Tree[rank[ProID]].cnt[i]-=1;
							}
						}
					}

					//get the new value of shortcut
				//	cout<<Cw<<" increase to ";
					Cw=INF; int countwt=0;

					for(int i=0;i<Neighbors[ProID].size();i++){
						if(Neighbors[ProID][i].first==Cid){
							Cw=Neighbors[ProID][i].second;//the weight value in the original graph
							countwt=1;
							break;
						}
					}

					int ssw,wtt,wid;
					vector<int> Wnodes; //Wnodes.clear();
					*//*if(ProID<Cid)
						Wnodes=SCconNodes[make_pair(ProID,Cid)]; //cout<<"wid num "<<Wnodes.size()<<endl;
					else
						Wnodes=SCconNodes[make_pair(Cid,ProID)];*//*

					if(ProID<Cid)
						Wnodes=SCconNodesMT[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
					else
						Wnodes=SCconNodesMT[Cid][ProID];
					for(int i=0;i<Wnodes.size();i++){
						wid=Wnodes[i];
						for(int j=0;j<Tree[rank[wid]].vert.size();j++){
							if(Tree[rank[wid]].vert[j].first==ProID){
								ssw=Tree[rank[wid]].vert[j].second.first;
							}
							if(Tree[rank[wid]].vert[j].first==Cid){
								wtt=Tree[rank[wid]].vert[j].second.first;
							}
						}

						if(ssw+wtt<Cw){
							Cw=ssw+wtt;
							countwt=1;
						}else if(ssw+wtt==Cw){
							countwt+=1;
						}
					}

					//cout<<Cw<<endl;
					//refresh the shortcut to the new value
					for(int i=0;i<Tree[rank[ProID]].vert.size();i++){
						if(Tree[rank[ProID]].vert[i].first==Cid){
							Tree[rank[ProID]].vert[i].second.first=Cw;
							Tree[rank[ProID]].vert[i].second.second=countwt;
							break;
						}
					}
				}

				if(influence){
					ProBeginID=ProID;
					ProBeginIDValue=true;
				}

				ProH-=1;
				ProID=Tree[Tree[rank[ProID]].pa].uniqueVertex;
			}


			if(ProBeginIDValue){
				vector<int> line1; //line1.clear();
				line1.reserve(heightMax);
				pachid=Tree[Tree[rank[ProBeginID]].pa].uniqueVertex;
				while(Tree[rank[pachid]].height>1){
					line1.insert(line1.begin(),pachid);
					pachid=Tree[Tree[rank[pachid]].pa].uniqueVertex;
				}
				line1.insert(line1.begin(),pachid);
		//cout<<"888888888888"<<endl;
				eachNodeProcessIncrease1(rank[ProBeginID],line1,ChangeNum,Tree,rank,VidtoTNid);
		//cout<<"9999999999999"<<endl;
			}
		}
		//return ChangeNum;
}*/

//Map version
void Graph::IncreaseH2H(int a,int b, int oldW, int newW, unordered_map<int,vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<int>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid){
    int ChangeNum=0;
    map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
    //OCdis.clear();
//cout<<"11111111111111"<<endl;
    if(Neighbors.find(a) == Neighbors.end() || Neighbors.find(b) == Neighbors.end()){//if not found a, b
        cout<<"Wrong! There is no a and b in Neighbors!!"<<endl;
        exit(1);
    }
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
//cout<<"22222222222222"<<endl;
    int lid,hid;
    if(NodeOrder[a]<NodeOrder[b]){
        lid=a;hid=b;
    }else{
        lid=b;hid=a;
    }

    int IniH=Tree[rank[lid]].height;//the height where weight change begins
    int ProH=Tree[rank[lid]].height; int ProID=lid;
    vector<set<int>> SCre;//record the shortcut change in each height
    set<int> vec; //vec.clear();
    SCre.assign(IniH+1,vec);
    int MinH;

    vector<int> line; //line.clear();
    line.reserve(heightMax);
    int pachid=ProID;
    while(Tree[rank[pachid]].height>1){
        line.insert(line.begin(),pachid);
        pachid=Tree[Tree[rank[pachid]].pa].uniqueVertex;
    }
    line.insert(line.begin(),pachid);

    bool tri=false;
    for(int i=0;i<Tree[rank[lid]].vert.size();i++){
        if(Tree[rank[lid]].vert[i].first==hid){
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
    if(tri){
        //shortcut update
        while(ProH>=MinH){
            influence=false;
            vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
            for(auto it=SCre[ProH].begin();it!=SCre[ProH].end();it++){
                int Cid=*it; int Cw=OCdis[make_pair(ProID,Cid)];
                int cidH=Tree[rank[Cid]].height-1;

                map<int,int> Hnei; //Hnei.clear();
                vector<pair<int,int>> Lnei; //Lnei.clear();
                for(int j=0;j<Vert.size();j++){
                    if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                        Hnei[Vert[j].first]=Vert[j].second.first;
                    }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                        Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                    }
                }
                //check the affected shortcuts
                int hid,lid;
                for(int j=0;j<Tree[rank[Cid]].vert.size();j++){
                    hid=Tree[rank[Cid]].vert[j].first;
                    if(Hnei.find(hid)!=Hnei.end()){
                        if(Cw+Hnei[hid]==Tree[rank[Cid]].vert[j].second.first){
                            Tree[rank[Cid]].vert[j].second.second-=1;
                            if(Tree[rank[Cid]].vert[j].second.second<1){
                                SCre[Tree[rank[Cid]].height].insert(hid);
                                if(Tree[rank[Cid]].height<MinH) MinH=Tree[rank[Cid]].height;
                                OCdis[make_pair(Cid,hid)]=Cw+Hnei[hid];
                            }
                        }
                    }
                }
                for(int j=0;j<Lnei.size();j++){
                    lid=Lnei[j].first;
                    for(int k=0;k<Tree[rank[lid]].vert.size();k++){
                        if(Tree[rank[lid]].vert[k].first==Cid){
                            if(Tree[rank[lid]].vert[k].second.first==Cw+Lnei[j].second){
                                Tree[rank[lid]].vert[k].second.second-=1;
                                if(Tree[rank[lid]].vert[k].second.second<1){
                                    SCre[Tree[rank[lid]].height].insert(Cid);
                                    if(Tree[rank[lid]].height<MinH) MinH=Tree[rank[lid]].height;
                                    OCdis[make_pair(lid,Cid)]=Cw+Lnei[j].second;
                                }
                            }
                            break;
                        }
                    }
                }

                //before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
                if(Tree[rank[ProID]].FN[cidH]){
                    influence=true;
                    //higher than Cid
                    for(int i=0;i<cidH;i++){
                        if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[Cid]].dis[i]){
                            Tree[rank[ProID]].cnt[i]-=1;
                        }
                    }

                    //equal to Cid
                    Tree[rank[ProID]].FN[cidH]=false;
                    Tree[rank[ProID]].cnt[cidH]-=1;

                    //lower than Cid
                    for(int i=cidH+1;i<Tree[rank[ProID]].dis.size();i++){
                        if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[line[i]]].dis[cidH]){
                            Tree[rank[ProID]].cnt[i]-=1;
                        }
                    }
                }

                //get the new value of shortcut
                //	cout<<Cw<<" increase to ";
                Cw=INF; int countwt=0;

                if(Neighbors.find(ProID) != Neighbors.end()){//if found
                    for(int i=0;i<Neighbors[ProID].size();i++){
                        if(Neighbors[ProID][i].first==Cid){
                            Cw=Neighbors[ProID][i].second;//the weight value in the original graph
                            countwt=1;
                            break;
                        }
                    }
                }

                int ssw,wtt,wid;
                vector<int> Wnodes; //Wnodes.clear();
                /*if(ProID<Cid)
                    Wnodes=SCconNodes[make_pair(ProID,Cid)]; //cout<<"wid num "<<Wnodes.size()<<endl;
                else
                    Wnodes=SCconNodes[make_pair(Cid,ProID)];*/

                if(ProID<Cid)
                    Wnodes=SCconNodesMT[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
                else
                    Wnodes=SCconNodesMT[Cid][ProID];
                for(int i=0;i<Wnodes.size();i++){
                    wid=Wnodes[i];
                    for(int j=0;j<Tree[rank[wid]].vert.size();j++){
                        if(Tree[rank[wid]].vert[j].first==ProID){
                            ssw=Tree[rank[wid]].vert[j].second.first;
                        }
                        if(Tree[rank[wid]].vert[j].first==Cid){
                            wtt=Tree[rank[wid]].vert[j].second.first;
                        }
                    }

                    if(ssw+wtt<Cw){
                        Cw=ssw+wtt;
                        countwt=1;
                    }else if(ssw+wtt==Cw){
                        countwt+=1;
                    }
                }

                //cout<<Cw<<endl;
                //refresh the shortcut to the new value
                for(int i=0;i<Tree[rank[ProID]].vert.size();i++){
                    if(Tree[rank[ProID]].vert[i].first==Cid){
                        Tree[rank[ProID]].vert[i].second.first=Cw;
                        Tree[rank[ProID]].vert[i].second.second=countwt;
                        break;
                    }
                }
            }

            if(influence){
                ProBeginID=ProID;
                ProBeginIDValue=true;
            }

            ProH-=1;
            ProID=Tree[Tree[rank[ProID]].pa].uniqueVertex;
        }


        if(ProBeginIDValue){
            vector<int> line1; //line1.clear();
            line1.reserve(heightMax);
            pachid=Tree[Tree[rank[ProBeginID]].pa].uniqueVertex;
            while(Tree[rank[pachid]].height>1){
                line1.insert(line1.begin(),pachid);
                pachid=Tree[Tree[rank[pachid]].pa].uniqueVertex;
            }
            line1.insert(line1.begin(),pachid);
            //cout<<"888888888888"<<endl;
            eachNodeProcessIncrease1(rank[ProBeginID],line1,ChangeNum,Tree,rank,VidtoTNid);
            //cout<<"9999999999999"<<endl;
        }
    }
    //return ChangeNum;
}

void Graph::eachNodeProcessIncrease1(int children, vector<int>& line, int& changelabel, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid){
	int childID=Tree[children].uniqueVertex;
	int childH=Tree[children].height-1;
	for(int i=0;i<Tree[children].dis.size();i++){
		if(Tree[children].cnt[i]==0){
			changelabel+=1;
			//firstly, check which dis can be infected
			int disBF=Tree[children].dis[i];
			int PID;
			//chidlID
			for(int k=0;k<VidtoTNid[childID].size();k++){
				PID=VidtoTNid[childID][k];
				if(Tree[PID].FN[childH] && Tree[PID].dis[i]==disBF+Tree[PID].dis[childH]){
					Tree[PID].cnt[i]-=1;
				}
			}

			//line[i]
			for(int k=0;k<VidtoTNid[line[i]].size();k++){
				PID=VidtoTNid[line[i]][k];
				if(PID>children){
					if(Tree[PID].FN[i] && Tree[PID].dis[childH]==disBF+Tree[PID].dis[i]){
						Tree[PID].cnt[childH]-=1;
					}
				}
			}

			//secondly, calculate the actual distance
			int dis=INF; int count=0;
			int Dvb; int b,bH; int DDvb=INF;
			for(int j=0;j<Tree[children].vert.size();j++){
				Dvb=Tree[children].vert[j].second.first;
				b=Tree[children].vert[j].first;
				bH=Tree[rank[b]].height-1;
				if(bH<i){
					if(Dvb+Tree[rank[line[i]]].dis[bH]<dis){
						dis=Dvb+Tree[rank[line[i]]].dis[bH];
						count=1;
					}else if(Dvb+Tree[rank[line[i]]].dis[bH]==dis){
						count+=1;
					}
				}else if(bH==i){
					DDvb=Dvb;
					if(Dvb<dis){
						dis=Dvb;
						count=1;
					}else if(Dvb==dis){
						count+=1;
					}
				}else{
					if(Dvb+Tree[rank[b]].dis[i]<dis){
						dis=Dvb+Tree[rank[b]].dis[i];
						count=1;
					}else if(Dvb+Tree[rank[b]].dis[i]==dis){
						count+=1;
					}
				}
			}
			if(DDvb==dis) Tree[children].FN[i]=true;
			Tree[children].dis[i]=dis;
			Tree[children].cnt[i]=count;
		}
	}

	line.push_back(childID);
	for(int i=0;i<Tree[children].ch.size();i++){
		eachNodeProcessIncrease1(Tree[children].ch[i],line,changelabel,Tree,rank,VidtoTNid);
	}
	line.pop_back();
}

//Query within one partition
/*int Graph::QueryH2HPartition(int ID1, int ID2, int PID){
	if(AdjaParti[PID][ID1].size()!=0 && AdjaParti[PID][ID2].size()!=0){//consider the vertex connectivity in subgraph
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
	}else
		return INF;
}*/

//Query within one partition: No LCA version
int Graph::QueryH2HPartition(int ID1, int ID2, int PID){
    int d=INF;

    int r1=ranks[PID][ID1], r2=ranks[PID][ID2];
    if(AdjaPartiM[PID].find(ID1) == AdjaPartiM[PID].end()) {
//        cout<<"QueryH2HPartition! "<<ID1<<" of "<<PID<<" is empty!"<<endl;
    }else if(AdjaPartiM[PID].find(ID2) == AdjaPartiM[PID].end()){
//        cout<<"QueryH2HPartition! "<<ID2<<" of "<<PID<<" is empty!"<<endl;
    }
    else{
        if(AdjaPartiM[PID][ID1].size()!=0 && AdjaPartiM[PID][ID2].size()!=0){//consider the vertex connectivity in subgraph
            if(ID1==ID2) return 0;

            unordered_map<int,int>::iterator it;
            int hub, dis1, dis2;
            int hubfinal,dis1final,dis2final;

            int height1 = Trees[PID][r1].dis.size();
            int height2 = Trees[PID][r2].dis.size();
//        cout<<"Tree depth: "<<height1<< " "<<height2<<endl;
            int temp_dis;
            if(height1 > height2){
                for(int i=0;i<height2;++i){
                    if(Trees[PID][r1].vAncestor[i] == Trees[PID][r2].vAncestor[i]){//if common ancestor
                        temp_dis = Trees[PID][r1].dis[i] + Trees[PID][r2].dis[i];
                        if(d > temp_dis){
                            d = temp_dis;
                        }
                    }
                    else{
                        break;
                    }
                }
            }else{
                for(int i=0;i<height1;++i){
                    if(Trees[PID][r1].vAncestor[i] == Trees[PID][r2].vAncestor[i]){//if common ancestor
                        temp_dis = Trees[PID][r1].dis[i] + Trees[PID][r2].dis[i];
                        if(d > temp_dis){
                            d = temp_dis;
                        }
                    }else{
                        break;
                    }
                }
            }
            if(d == INF){
                cout<<"Wrong! d = "<<INF<<" "<<AdjaPartiM[PID][ID1].size()<<" "<<AdjaPartiM[PID][ID2].size()<<endl;
                for(int i=0;i<Trees[PID][r1].vAncestor.size();++i){
                    cout<<Trees[PID][r1].vAncestor[i]<<"\t";
                }
                cout<<endl;
                for(int i=0;i<Trees[PID][r2].vAncestor.size();++i){
                    cout<<Trees[PID][r2].vAncestor[i]<<"\t";
                }
                cout<<endl;
            }
        }
    }



    return d;
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
