/*
 * Decomp.cpp
 *
 *  Created on: 13 Oct 2022
 *      Author: zhangmengxuan
 */
#include "head.h"

void Graph::H2HconCore(){
	CHconsCore();
	//cout<<"CH contraction"<<endl;
	makeTreeCore();
	//cout<<"Make Tree"<<endl;
}

vector<int> _DD_,_DD2_;//
struct DegComp1{
	int x;
	DegComp1(int _x){
		x=_x;
	}
	bool operator< (const DegComp1 d) const{
		if(_DD_[x]!=_DD_[d.x])
			return _DD_[x]<_DD_[d.x];
		if(_DD2_[x]!=_DD2_[x])
			return _DD2_[x]<_DD2_[d.x];
		return x<d.x;
	}
};

void Graph::CHconsCore(){
	//initialize E
	map<int,int> m;
	Emap.assign(nodenum,m);//edge weight map
	for(int i=0;i<Neighbor.size();i++){
		for(int j=0;j<Neighbor[i].size();j++)
			Emap[i].insert(make_pair(Neighbor[i][j].first,Neighbor[i][j].second));
	}

	_DD_.assign(nodenum,0);_DD2_.assign(nodenum,0);
	DD.assign(nodenum,0); DD2.assign(nodenum,0);

	set<DegComp1> Deg;//min first
	int degree;
	for(int i=0;i<nodenum;i++){
		degree=Neighbor[i].size();
		if(degree!=0){//get degree
			_DD_[i]=degree;
			_DD2_[i]=degree;
			DD[i]=degree;
			DD2[i]=degree;
			Deg.insert(DegComp1(i));
		}
	}

	vNodeOrder.clear();
	//vector<bool> exist;
	existCore.assign(nodenum,true);//if in the core, all vertices is originally in core
	vector<bool> change;
	change.assign(nodenum,false);//whether the neighbors have changed

	vector<pair<int,int>> vect;
	NeighborConCore.assign(nodenum,vect);

	bool CutLable=false;
	int count=0;
	while(!Deg.empty()){
		count+=1;
		int x=(*Deg.begin()).x;//minimum degree first

		while(true){
			if(change[x]){
				Deg.erase(DegComp1(x));
				_DD_[x]=DD[x];
				_DD2_[x]=DD2[x];
				Deg.insert(DegComp1(x));
				change[x]=false;
				x=(*Deg.begin()).x;
			}else
				break;
		}

		vNodeOrder.push_back(x);//least important vertex first
		Deg.erase(Deg.begin());

		vector<pair<int,int>> Neigh; //Neigh.clear();

		for(auto it=Emap[x].begin();it!=Emap[x].end();it++){
			if(existCore[(*it).first]){//if in the core
				Neigh.push_back(*it);
			}
		}
		NeighborConCore[x].assign(Neigh.begin(),Neigh.end());
		//in NeighborConCore, store the higher neighbors for vertices in tree;
		//and store both lower and higher neighbors for vertices in core

		/*if(NeighborConCore[x].size()==0){
			cout<<"Isolated vertex "<<x<<" Further check "<<Neighbor[x].size()<<endl;
			for(int k=0;k<10;k++){
				srand (time(NULL));
				int s=rand()%nodenum;
				cout<<Dijkstra(x,s)<<" ";
			}
			cout<<endl;
		}*/
        /// if still need to contract
		if(!CutLable){
			if(Neigh.size()<Width){//if the neighbor is smaller than tree width threshold, the vertex will be a part of tree
			//if(DD[x]<Width){
				existCore[x]=false;
				//delete the star
				for(int i=0;i<Neigh.size();i++){
					int y=Neigh[i].first;
					deleteECore(x,y);
					change[y]=true;
				}
				//add all-pair neighbors
				for(int i=0;i<Neigh.size();i++){
					for(int j=i+1;j<Neigh.size();j++){
						insertECore(Neigh[i].first,Neigh[j].first,Neigh[i].second+Neigh[j].second);
						change[Neigh[i].first]=true;
						change[Neigh[j].first]=true;
					}
				}
			}else{//else vertices in the core
				//if(!CutLable){
					HighestOrder=vNodeOrder.size()-1;
					//cout<<"HighesetOrder for vertex in periphery "<<HighestOrder<<",x "<<x<<endl;
					CutLable=true;
				//}
			}
		}


	}

	NodeOrder.assign(nodenum,-1);
	for(int k=0;k<vNodeOrder.size();k++){
		NodeOrder[vNodeOrder[k]]=k;
	}
}
//function of erasing edge (u,v), i.e., erase v from u's adjacency list, vice verse.
void Graph::deleteECore(int u,int v){
	if(Emap[u].find(v)!=Emap[u].end()){
		Emap[u].erase(Emap[u].find(v));
		DD[u]--;
	}

	if(Emap[v].find(u)!=Emap[v].end()){
		Emap[v].erase(Emap[v].find(u));
		DD[v]--;
	}
}
//function of inserting edge (u,v)
void Graph::insertECore(int u,int v,int w){
	if(Emap[u].find(v)==Emap[u].end()){//if not found
		Emap[u].insert(make_pair(v,w));
		DD[u]++;
		DD2[u]++;
	}
	else{//if found
		if(Emap[u][v]>w)
			Emap[u][v]=w;
	}

	if(Emap[v].find(u)==Emap[v].end()){
		Emap[v].insert(make_pair(u,w));
		DD[v]++;
		DD2[v]++;
	}
	else{
		if(Emap[v][u]>w)
			Emap[v][u]=w;
	}
}

int Graph::matchCore(int x,vector<pair<int,int>> &vert){
	int nearest=vert[0].first;
	for(int i=1;i<vert.size();i++){
		if(NodeOrder[vert[i].first]<NodeOrder[nearest])//get the least node order
			nearest=vert[i].first;
	}
	if(existCore[nearest])//if it exists in core
		return 0;
	else{
		//cout<<nearest<<" "<<rankCore[nearest]<<endl;
		return rankCore[nearest];
	}

}

void Graph::makeTreeCore(){
	//rank.assign(HighestOrder+2,0);
	rankCore.clear();
	rankCore.assign(nodenum,-1);
	int len=HighestOrder-1;

	NodeCore root;
	root.uniqueVertex=-1;
	root.height=1;
	TreeCore.push_back(root);

	for(;len>=0;len--){//check the vertices with order lower than HighestOrder
		int x=vNodeOrder[len];
		if(existCore[x]){
			cout<<"Wrong: should be out of core"<<endl;
		}
		NodeCore nod;
		nod.vert=NeighborConCore[x];//the higher order neighbors
		nod.uniqueVertex=x;
		int pa=matchCore(x,NeighborConCore[x]);

		//cout<<"pa "<<pa<<endl;

		TreeCore[pa].ch.push_back(TreeCore.size());
		nod.pa=pa;
		nod.height=TreeCore[pa].height+1;
		rankCore[x]=TreeCore.size();//the position of tree

		if(pa==0)
			nod.treeroot=rankCore[x];
		else
			nod.treeroot=TreeCore[pa].treeroot;

		TreeCore.push_back(nod);
	}
}




