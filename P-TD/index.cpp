/*
 * func.cpp
 *
 *  Created on: 19 Sep 2022
 *      Author: zhangmengxuan
 */
#include "head.h"

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

void Graph::H2Hindex(){
	vector<pair<int,pair<int,int>>> vect;
	NeighborCon.assign(nodenum,vect);

	//CHindex();//single threads
	CHindexMT();//multiple threads
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
		if(x!=-1){//to identify and exclude the isolated vertices
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
}

void Graph::NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x){
	sm->wait();
	int ID1, w1;
	int ID2, w2;
	for(int k=p.first;k<p.second;k++){
		ID1=Neighvec[k].first;
		w1=Neighvec[k].second.first;
		for(int h=0;h<Neighvec.size();h++){
			ID2=Neighvec[h].first;
			w2=Neighvec[h].second.first;
			insertEMTorder(ID1, ID2, w1+w2);
			if(ID1<ID2)
				SCconNodesMT[ID1][ID2].push_back(x);
		}
	}
	sm->notify();
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
		if(rank[vert[i].first]>rank[nearest])
			nearest=vert[i].first;
	}
	int p=rank[nearest];
	return p;
}

void Graph::makeIndex(){
	makeRMQ();

	//initialize
	vector<int> list; //list.clear();
	list.push_back(Tree[0].uniqueVertex);
	Tree[0].pos.clear();
	Tree[0].pos.push_back(0);

	for(int i=0;i<Tree[0].ch.size();i++){
		makeIndexDFS(Tree[0].ch[i],list);
	}
}

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
void Graph::Decrease(int a,int b, int newW, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax){
	map<int,int> checkedDis;//map<tree node ID, distance index>

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

	//int ProBeginH;
	int ProBeginID;
	if(tri){
	//cout<<"Bottom-up ;;;;;;;;;;;;;;;;;; "<<endl;
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
			}

			ProH-=1;
			ProID=Tree[Tree[rank[ProID]].pa].uniqueVertex;
		}

		vector<int> linee; //linee.clear();
		linee.reserve(heightMax);
		int pachidd=Tree[Tree[rank[ProBeginID]].pa].uniqueVertex;
		while(Tree[rank[pachidd]].height>1){
			linee.insert(linee.begin(),pachidd);
			pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
		}
		linee.insert(linee.begin(),pachidd);

		//top-down process
		EachNodeProBDis5(rank[ProBeginID], linee, vertexIDChL, checkedDis, Tree, rank);
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

void Graph::Increase(int a,int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<int>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid){
	int ChangeNum=0;
	map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
		//OCdis.clear();

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
				}

				ProH-=1;
				ProID=Tree[Tree[rank[ProID]].pa].uniqueVertex;
			}
		}

		vector<int> line1; //line1.clear();
		line1.reserve(heightMax);
		pachid=Tree[Tree[rank[ProBeginID]].pa].uniqueVertex;
		while(Tree[rank[pachid]].height>1){
			line1.insert(line1.begin(),pachid);
			pachid=Tree[Tree[rank[pachid]].pa].uniqueVertex;
		}
		line1.insert(line1.begin(),pachid);

		eachNodeProcessIncrease1(rank[ProBeginID],line1,ChangeNum,Tree,rank,VidtoTNid);

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

void Graph::DecreaseBatch(vector<pair<pair<int,int>,int>>& wBatch, vector<vector<pair<int,int>>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax){
	map<int,int> checkedDis;

	for(int i=0;i<Tree.size();i++){
		Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
	}

	//NodeOrderss.clear();
	NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
	vector<set<int>> SCre; //SCre.clear();
	set<int> ss; //ss.clear();
	SCre.assign(nodenum,ss);//{vertexID, set<int>}
	set<OrderCompp> OC; //OC.clear();//vertexID in decreasing node order

	set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed

	int a,b,newW,lid,hid;
	for(int k=0;k<wBatch.size();k++){
		a=wBatch[k].first.first; b=wBatch[k].first.second ;newW=wBatch[k].second;
		if(NodeOrder[a]<NodeOrder[b]){
			lid=a;hid=b;
		}else{
			lid=b;hid=a;
		}

		for(int i=0;i<Neighbor[a].size();i++){
			if(Neighbor[a][i].first==b){
				Neighbor[a][i].second=newW;
				break;
			}
		}
		for(int i=0;i<Neighbor[b].size();i++){
			if(Neighbor[b][i].first==a){
				Neighbor[b][i].second=newW;
				break;
			}
		}

		for(int i=0;i<Tree[rank[lid]].vert.size();i++){
			if(Tree[rank[lid]].vert[i].first==hid){
				if(Tree[rank[lid]].vert[i].second.first>newW){
					Tree[rank[lid]].vert[i].second.first=newW;
					Tree[rank[lid]].vert[i].second.second=1;
					SCre[lid].insert(hid);
					OC.insert(OrderCompp(lid));
				}else if(Tree[rank[lid]].vert[i].second.first==newW){
					Tree[rank[lid]].vert[i].second.second+=1;
				}
				break;
			}
		}

	}

	vector<int> ProBeginVertexSet; //ProBeginVertexSet.clear();
	vector<int> ProBeginVertexSetNew;
	int ProBeginVertexID;
	int ProID;
	//processing the stars
	while(!OC.empty()){
		ProID=(*OC.begin()).x;
		OC.erase(OC.begin());
		vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
		bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
		for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
			int Cid=*it; int Cw;
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

			if(Tree[rank[ProID]].dis[cidH]>Cw){
				Tree[rank[ProID]].dis[cidH]=Cw;
				Tree[rank[ProID]].FN[cidH]=true;
				ProIDdisCha=true;
				Tree[rank[ProID]].DisRe.insert(Cid);
			}else if(Tree[rank[ProID]].dis[cidH]==Cw){
				Tree[rank[ProID]].FN[cidH]=true;
			}

			int hid,hidHeight,lid,lidHeight,wsum;
			for(int j=0;j<Tree[rank[Cid]].vert.size();j++){
				hid=Tree[rank[Cid]].vert[j].first;hidHeight=Tree[rank[hid]].height-1;
				if(Hnei.find(hid)!=Hnei.end()){
					wsum=Cw+Hnei[hid];
					if(wsum<Tree[rank[Cid]].vert[j].second.first){
						Tree[rank[Cid]].vert[j].second.first=wsum;
						Tree[rank[Cid]].vert[j].second.second=1;
						SCre[Cid].insert(hid);
						OC.insert(OrderCompp(Cid));
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
							SCre[lid].insert(Cid);
							OC.insert(OrderCompp(lid));
						}else if(Tree[rank[lid]].vert[k].second.first==wsum){
							Tree[rank[lid]].vert[k].second.second+=1;
						}

						break;
					}
				}
			}
		}

		if(ProIDdisCha){//if the distance labeling is dectected changed
			vertexIDChL.insert(ProID);
			ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSet.size()+1);
			ProBeginVertexSetNew.push_back(ProID);
			int rnew=rank[ProID],r;
			for(int i=0;i<ProBeginVertexSet.size();i++){
				r=rank[ProBeginVertexSet[i]];
				if(LCAQueryOverlay(rnew,r)!=rnew){
					ProBeginVertexSetNew.push_back(ProBeginVertexSet[i]);
				}
			}
			ProBeginVertexSet=ProBeginVertexSetNew;
		}
	}

	//cout<<"Finish bottom-up refresh"<<endl;
	for(int i=0;i<ProBeginVertexSet.size();i++){
		ProBeginVertexID=ProBeginVertexSet[i];
		vector<int> linee; //linee.clear();
		linee.reserve(heightMax);
		int pachidd=Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
		while(Tree[rank[pachidd]].height>1){
			linee.insert(linee.begin(),pachidd);
			pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
		}
		linee.insert(linee.begin(),pachidd);
		EachNodeProBDis5(rank[ProBeginVertexID], linee, vertexIDChL,checkedDis,Tree,rank);
	}
	//return checkedDis.size();
}

void Graph::IncreaseBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<int>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid){
	int checknum=0;
	map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
	//OCdis.clear();

	//NodeOrderss.clear();
	NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
	vector<set<int>> SCre; //SCre.clear();
	set<int> ss; //ss.clear();
	SCre.assign(nodenum,ss);//{vertexID, set<int>}
	set<OrderCompp> OC; OC.clear();//vertexID in decreasing node order

	for(int k=0;k<wBatch.size();k++){
		int a=wBatch[k].first.first;
		int b=wBatch[k].first.second;
		int oldW=wBatch[k].second.first;
		int newW=wBatch[k].second.second;

		if(oldW!=newW){
		for(int i=0;i<Neighbor[a].size();i++){
			if(Neighbor[a][i].first==b){
				Neighbor[a][i].second=newW;
				break;
			}
		}
		for(int i=0;i<Neighbor[b].size();i++){
			if(Neighbor[b][i].first==a){
				Neighbor[b][i].second=newW;
				break;
			}
		}

		int lid,hid;
		if(NodeOrder[a]<NodeOrder[b]){
			lid=a;hid=b;
		}else{
			lid=b;hid=a;
		}

		for(int i=0;i<Tree[rank[lid]].vert.size();i++){
			if(Tree[rank[lid]].vert[i].first==hid){
				if(Tree[rank[lid]].vert[i].second.first==oldW){
					Tree[rank[lid]].vert[i].second.second-=1;
					if(Tree[rank[lid]].vert[i].second.second<1){
						OCdis[make_pair(lid,hid)]=oldW;
						SCre[lid].insert(hid);
						OC.insert(OrderCompp(lid));
					}
				}
				break;
			}
		}
	}
	}

	vector<int> ProBeginVertexSet; ProBeginVertexSet.clear();
	vector<int> ProBeginVertexSetNew;
	bool influence;
	int ProID; vector<int> line;
	while(!OC.empty()){
		ProID=(*OC.begin()).x;
		OC.erase(OC.begin());
		vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
		influence=false;

		//each ProID corresponds to a line
		line.clear(); line.reserve(heightMax);
		int pachid=ProID;
		while(Tree[rank[pachid]].height>1){
			line.insert(line.begin(),pachid);
			pachid=Tree[Tree[rank[pachid]].pa].uniqueVertex;
		}
		line.insert(line.begin(),pachid);

		for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
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
							SCre[Cid].insert(hid);
							OC.insert(OrderCompp(Cid));
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
								SCre[lid].insert(Cid);
								OC.insert(OrderCompp(lid));
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

			for(int i=0;i<Neighbor[ProID].size();i++){
				if(Neighbor[ProID][i].first==Cid){
					Cw=Neighbor[ProID][i].second;//the weight value in the original graph
					countwt=1;
					break;
				}
			}

			int ssw,wtt,wid;
			vector<int> Wnodes;
			Wnodes.clear();
			/*if(SCconNodes.find(make_pair(ProID,Cid))!=SCconNodes.end())
				Wnodes=SCconNodes[make_pair(ProID,Cid)];
			else if(SCconNodes.find(make_pair(Cid,ProID))!=SCconNodes.end())
				Wnodes=SCconNodes[make_pair(Cid,ProID)];*/
			if(ProID<Cid)
				Wnodes=SCconNodesMT[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
			else
				Wnodes=SCconNodesMT[Cid][ProID];
			if(Wnodes.size()>0){
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
			ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSet.size()+1);
			ProBeginVertexSetNew.push_back(ProID);
			int rnew=rank[ProID],r;
			for(int i=0;i<ProBeginVertexSet.size();i++){
				r=rank[ProBeginVertexSet[i]];
				if(LCAQueryOverlay(rnew,r)!=rnew){
					ProBeginVertexSetNew.push_back(ProBeginVertexSet[i]);
				}
			}
			ProBeginVertexSet=ProBeginVertexSetNew;
		}

	}

	int ProBeginVertexID;
	for(int i=0;i<ProBeginVertexSet.size();i++){
		ProBeginVertexID=ProBeginVertexSet[i];
		vector<int> linee; //linee.clear();
		linee.reserve(heightMax);
		int pachidd=Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
		while(Tree[rank[pachidd]].height>1){
			linee.insert(linee.begin(),pachidd);
			pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
		}
		linee.insert(linee.begin(),pachidd);

		eachNodeProcessIncrease1(rank[ProBeginVertexID], linee,checknum,Tree,rank,VidtoTNid);
	}
	//return checknum;
}
