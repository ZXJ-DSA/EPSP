/*
 * Decomp.cpp
 *
 *  Created on: 13 Oct 2022
 *      Author: zhangmengxuan
 */
#include "head3.h"

//void Graph::H2HconCore(){
//	CHconsCore();
//	//cout<<"CH contraction"<<endl;
//	makeTree();
//	//cout<<"Make Tree"<<endl;
//}

vector<int> _DD_,_DD2_;//true degree, temporal degree
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
//Function of contracting vertices for H2H
void Graph::MDEContract(){
    HighestOrder = nodenum;
    //for H2H update
    map<int, vector<pair<int,int>>> mi;
    SCconNodesMT.assign(nodenum, mi);

    //initialize E
    map<int,pair<int,int>> m;
    E.assign(nodenum,m);
    for(int i=0;i<Neighbor.size();i++){
        for(int j=0;j<Neighbor[i].size();j++)
            E[i].insert(make_pair(Neighbor[i][j].first,make_pair(Neighbor[i][j].second,1)));
    }

    _DD_.assign(nodenum,0);_DD2_.assign(nodenum,0);
    DD.assign(nodenum,0); DD2.assign(nodenum,0);

    set<DegComp1> Deg;//min first
//    vector<bool> active(nodenum,false);//flag that indicate whether a vertex is active for contraction
    int degree;
    for(int i=0;i<nodenum;i++){
        degree=Neighbor[i].size();

        if(degree > 0){//get degree
            _DD_[i]=degree;
            _DD2_[i]=degree;
            DD[i]=degree;
            DD2[i]=degree;
            Deg.insert(DegComp1(i));
//            active[i] = true;
        }else{
            cout<<"Wrong!! Degree of "<<i<<" is "<<degree<<endl;
            exit(1);
        }
    }

    vNodeOrder.clear();
    //vector<bool> exist;
    existCore.assign(nodenum,true);//if in the core, all vertices is originally in core
    vector<bool> change;
    change.assign(nodenum,false);//whether the neighbor (degree) has changed

    vector<pair<int,pair<int,int>>> vect;
    NeighborCon.assign(nodenum,vect);//temporal graph to store Neighbors in the core, for graph contraction
    NeighborConCH.assign(nodenum,vect);

    bool CutLabel=false;
    int count=0;
    int ID1,ID2;

    //Get the order of all vertices by MDE
    while(!Deg.empty()){
        count+=1;
        int x=(*Deg.begin()).x;//minimum degree first

        while(change[x]){//update the degree if it is changed
            Deg.erase(DegComp1(x));
            _DD_[x]=DD[x];
            _DD2_[x]=DD2[x];
            Deg.insert(DegComp1(x));
            change[x]=false;
            x=(*Deg.begin()).x;
        }

//        if(_DD_[x] >= Width)
//            break;

        vNodeOrder.push_back(x);//least important vertex first
        Deg.erase(Deg.begin());

//        vector<pair<int,int>> Neigh; Neigh.clear();
//        for(auto it=E[x].begin();it!=E[x].end();it++){//for the neighbors of x
//            if(existCore[(*it).first]){//if in the core
//                Neigh.push_back(*it);
//            }
//        }
//        NeighborConCore[x].assign(Neigh.begin(),Neigh.end());//for core vertex, the neighbors are all truly-core vertices
        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();
        for(auto it=E[x].begin();it!=E[x].end();it++){
            if(existCore[(*it).first]){
                Neigh.push_back(*it);
            }
        }
        NeighborCon[x].assign(Neigh.begin(),Neigh.end());

        /// if still need to contract
        if(!CutLabel){
            if(Neigh.size()<Width){//if the neighbor is smaller than tree width threshold, the vertex will be a part of tree
                NeighborConCH[x].assign(Neigh.begin(),Neigh.end());//assign the periphery vertex

                //if(DD[x]<Width){
                existCore[x]=false;
                //delete the star
                for(int i=0;i<Neigh.size();i++){
                    int y=Neigh[i].first;
                    deleteECore(x,y);//delete x from y's adjacency list
                    change[y]=true;
                }
                //add all-pair neighbors
                for(int i=0;i<Neigh.size();i++){
                    ID1=Neigh[i].first;
                    for(int j=i+1;j<Neigh.size();j++){
                        ID2=Neigh[j].first;
                        insertECore(ID1,ID2,Neigh[i].second.first+Neigh[j].second.first);
                        /// For TD update
//                        if(ID1<ID2)
//                            SCconNodesMT[ID1][ID2].push_back(x);//supportive vertex, no direction, may contain the supportive vertices for shortcuts between interface vertices
//                        else if(ID1>ID2)
//                            SCconNodesMT[ID2][ID1].push_back(x);
                        SCconNodesMT[ID1][ID2].emplace_back(x,Neigh[i].second.first+Neigh[j].second.first);//supportive vertex, no direction, may contain the supportive vertices for shortcuts between interface vertices
                        SCconNodesMT[ID2][ID1].emplace_back(x,Neigh[i].second.first+Neigh[j].second.first);

                        change[ID1]=true;
                        change[ID2]=true;
                    }
                }
            }else{//else vertices in the core
                //if(!CutLable){
                HighestOrder=vNodeOrder.size()-1;
                //cout<<"HighesetOrder for vertex in periphery "<<HighestOrder<<",x "<<x<<endl;
                CutLabel=true;
                //}
            }
        }

    }

    NodeOrder.assign(nodenum,-1);
    for(int k=0;k<vNodeOrder.size();k++){
        NodeOrder[vNodeOrder[k]]=k;
    }



}

/*void Graph::H2HContract_New(){
    HighestOrder = nodenum;
    //for H2H update
    //map<int, vector<pair<int,int>>> mi;
    //SCconNodesMT.assign(nodenum, mi);

    //initialize E
    map<int,pair<int,int>> m;
    E.assign(nodenum,m);
    for(int i=0;i<Neighbor.size();i++){
        for(int j=0;j<Neighbor[i].size();j++)
            E[i].insert(make_pair(Neighbor[i][j].first,make_pair(Neighbor[i][j].second,1)));
    }

    _DD_.assign(nodenum,0);_DD2_.assign(nodenum,0);
    DD.assign(nodenum,0); DD2.assign(nodenum,0);

    set<DegComp1> Deg;//min first
//    vector<bool> active(nodenum,false);//flag that indicate whether a vertex is active for contraction
    int degree;
    for(int i=0;i<nodenum;i++){
        degree=Neighbor[i].size();

        if(degree > 0){//get degree
            _DD_[i]=degree;
            _DD2_[i]=degree;
            DD[i]=degree;
            DD2[i]=degree;
            Deg.insert(DegComp1(i));
//            active[i] = true;
        }else{
            cout<<"Wrong!! Degree of "<<i<<" is "<<degree<<endl;
            exit(1);
        }
    }

    vNodeOrder.clear();
    //vector<bool> exist;
    existCore.assign(nodenum,true);//if in the core, all vertices is originally in core
    vector<bool> change;
    change.assign(nodenum,false);//whether the neighbor (degree) has changed

    vector<pair<int,pair<int,int>>> vect;
    NeighborCon.assign(nodenum,vect);//temporal graph to store Neighbors in the core, for graph contraction

    bool CutLabel=false;
    int count=0;
    int ID1,ID2;

    //Get the order of all vertices by MDE
    while(!Deg.empty()){
        count+=1;
        int x=(*Deg.begin()).x;//minimum degree first

        while(change[x]){//update the degree if it is changed
            Deg.erase(DegComp1(x));
            _DD_[x]=DD[x];
            _DD2_[x]=DD2[x];
            Deg.insert(DegComp1(x));
            change[x]=false;
            x=(*Deg.begin()).x;
        }

//        if(_DD_[x] >= Width)
//            break;

        vNodeOrder.push_back(x);//least important vertex first
        Deg.erase(Deg.begin());

//        vector<pair<int,int>> Neigh; Neigh.clear();
//        for(auto it=E[x].begin();it!=E[x].end();it++){//for the neighbors of x
//            if(existCore[(*it).first]){//if in the core
//                Neigh.push_back(*it);
//            }
//        }
//        NeighborConCore[x].assign(Neigh.begin(),Neigh.end());//for core vertex, the neighbors are all truly-core vertices
        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();
        for(auto it=E[x].begin();it!=E[x].end();it++){
            if(existCore[(*it).first]){
                Neigh.push_back(*it);
            }
        }
        NeighborCon[x].assign(Neigh.begin(),Neigh.end());

        /// if still need to contract
        if(!CutLabel){
            if(Neigh.size()<Width){//if the neighbor is smaller than tree width threshold, the vertex will be a part of tree
                //if(DD[x]<Width){
                existCore[x]=false;
                //delete the star
                for(int i=0;i<Neigh.size();i++){
                    int y=Neigh[i].first;
                    deleteECore(x,y);//delete x from y's adjacency list
                    change[y]=true;
                }
                //add all-pair neighbors
                for(int i=0;i<Neigh.size();i++){
                    ID1=Neigh[i].first;
                    for(int j=i+1;j<Neigh.size();j++){
                        ID2=Neigh[j].first;
                        insertECore(ID1,ID2,Neigh[i].second.first+Neigh[j].second.first);
                        /// For TD update
//                        if(ID1<ID2)
//                            SCconNodesMT[ID1][ID2].push_back(x);//supportive vertex, no direction, may contain the supportive vertices for shortcuts between interface vertices
//                        else if(ID1>ID2)
//                            SCconNodesMT[ID2][ID1].push_back(x);
                        //SCconNodesMT[ID1][ID2].emplace_back(x,Neigh[i].second.first+Neigh[j].second.first);//supportive vertex, no direction, may contain the supportive vertices for shortcuts between interface vertices
                        //SCconNodesMT[ID2][ID1].emplace_back(x,Neigh[i].second.first+Neigh[j].second.first);

                        change[ID1]=true;
                        change[ID2]=true;
                    }
                }
            }else{//else vertices in the core
                //if(!CutLable){
                HighestOrder=vNodeOrder.size()-1;
                //cout<<"HighesetOrder for vertex in periphery "<<HighestOrder<<",x "<<x<<endl;
                CutLabel=true;
                //}
            }
        }

    }

    NodeOrder.assign(nodenum,-1);
    for(int k=0;k<vNodeOrder.size();k++){
        NodeOrder[vNodeOrder[k]]=k;
    }



}*/

//Function of creating tree
void Graph::Create_tree(){
    cout<<"Creating tree..."<<endl;
    //// Get tree
    vector<int> vecemp; //vecemp.clear();
    VidtoTNid.assign(nodenum,vecemp);

    //rank.assign(HighestOrder+2,0);
    rank.clear();
    rank.assign(nodenum,0);//the vector index of tree nodes, map from vertex to tree node
    int len=HighestOrder-1;
    heightMax=0;

    Node root;//virtual root node
    int x = vNodeOrder[len];
    if(NeighborCon[x].empty()){
        cout<<"There exist non-virtual root!"<<endl;
        root.uniqueVertex=x;
        len--;
    }else{
        root.uniqueVertex=-1;
    }
    root.height=1;
    Tree.push_back(root);

    int nn;
    for(;len>=0;len--){//check the vertices with order lower than HighestOrder
        x=vNodeOrder[len];
        if(existCore[x]){
            cout<<"Wrong: should be out of core"<<endl; exit(1);
        }
        Node nod;
        nod.vert=NeighborCon[x];
        nod.uniqueVertex=x;
        int pa=matchCore(x,NeighborCon[x]);

        //cout<<"pa "<<pa<<endl;

        Tree[pa].ch.push_back(Tree.size());
        nod.pa=pa;
        nod.height=Tree[pa].height+1;
        /// for update
        nod.hdepth=Tree[pa].height+1;
        for(int i=0;i<NeighborCon[x].size();i++){//for the neighbors which have higher order
            nn=NeighborCon[x][i].first;
            if(existCore[nn])//skip core vertex
                continue;
            VidtoTNid[nn].push_back(Tree.size());//record the child tree node rank who has direct super edge to nn
            if(Tree[rank[nn]].hdepth<Tree[pa].height+1)
                Tree[rank[nn]].hdepth=Tree[pa].height+1;
        }
        if(nod.height>heightMax)
            heightMax=nod.height;

        rank[x]=Tree.size();//the position of tree, higher-order vertex has lower rank
        if(pa==0){//root node of this tree
            nod.treeroot=rank[x];//tree root is itself
//            Tree[0].ch.push_back(Tree.size());//get the children of virtual root node
        }
        else{
            nod.treeroot=Tree[pa].treeroot;
        }
        Tree.push_back(nod);
    }

    //// Get partitions
    vector<int> cr;//the vertex rank of root's children
    cr.clear();
    int temp_ch=0;
    int temp_vert=0;
    int sum_singular=0;
    for(int k=0;k<Tree[0].ch.size();k++){//Tree[0] is the core
        int childrank=Tree[0].ch[k];
        if(Tree[childrank].ch.size()>0){///only the root which has children is regarded as a tree
            cr.push_back(childrank);
            temp_ch += Tree[childrank].ch.size();
            temp_vert += Tree[childrank].vert.size();
        }
        else{
            ++sum_singular;
//            cout<<"Single vertex in tree "<< childrank<<"!!!"<<endl;
        }
    }
//    cout<<"Tree[0].ch.size(): "<<Tree[0].ch.size()<<endl;
//    cout<<"Accumulated Tree[childrank].ch.size(): "<<temp_ch<<endl;
//    cout<<"Accumulated Tree[childrank].vert.size(): "<<temp_vert<<endl;
    partiNum=cr.size();

    cout<<"Periphery Number: "<<partiNum<<endl;
    cout<<"Nominal periphery root number: "<<sum_singular+partiNum<<endl;
    cout<<"Vertex number in core: "<<nodenum - HighestOrder + sum_singular <<endl;

    ///Get boundary vertex
    PartiUpdateExt.assign(partiNum, false);//initiation of partition update vector
    vector<int> vec;
    vec.clear();
    BoundVertex.assign(partiNum,vec);
    set<int> sset;
    sset.clear();
    BoundVertexSet.assign(partiNum,sset);
    BoundTag.assign(nodenum, make_pair(false,set<int>()));//the boundary vertex is consisted by interface vertex and root vertex
    map<int,int> PartiRoot;//partition root & partition ID
    PartiRoot.clear();

    int ID1,ID2;

    for(int PID=0;PID<partiNum;PID++){
        int childrank=cr[PID];
        for(int i=0;i<Tree[childrank].vert.size();i++){
            ID1 = Tree[childrank].vert[i].first;
            BoundVertex[PID].push_back(ID1);
            BoundVertexSet[PID].insert(ID1);
            BoundTag[ID1].first=true;//interface vertex
            BoundTag[ID1].second.insert(PID);
        }
        BoundVertex[PID].push_back(Tree[childrank].uniqueVertex);
//        BoundVertexSet[PID].insert(Tree[childrank].uniqueVertex);
//        BoundTag[Tree[childrank].uniqueVertex]=true;//root vertex

        PartiRoot.insert(make_pair(childrank,PID));//map from tree id to partition id
    }

    //if(PartiRoot.size()!=cr.size())
    //cout<<"something wrong with the boundary node"<<endl;

    CoreTag.assign(nodenum,-1);//-1 indicates core vertex or root vertex, i>=0 indicates non-core vertex (i.e., the inner-partition vertex) and which partition it belongs to
    int NodeID,RootNode,parentNode;
    int count=0;
    for(int len=HighestOrder-1;len>=0;len--){
        NodeID=vNodeOrder[len];
        RootNode=Tree[rank[NodeID]].treeroot;
        parentNode=Tree[rank[NodeID]].pa;
        if(parentNode!=0){//if it is not root vertex
            CoreTag[NodeID]=PartiRoot[RootNode];
        }
        else if(!Tree[rank[NodeID]].ch.empty()){//if the root vertex has children
            CoreTag[NodeID]=PartiRoot[RootNode];
        }else{
            count++;
        }
    }
    cout<<"Single vertex periphery: "<<sum_singular<<" "<<count<<endl;

    /// Get partition info: AdjaGraph (only for CT-DS) and AdjaCore
    AdjaCoreMap.clear();
    map<int,int> mii;
    mii.clear();
    AdjaCoreMap.assign(nodenum,mii);

    /// for AdjaCore
    for(int NodeID=0;NodeID<nodenum;NodeID++){
        if(CoreTag[NodeID]==-1){//core vertex
            for(int nei=0;nei<NeighborCon[NodeID].size();nei++) {//for each neighbor
                assert(CoreTag[NodeID]==-1);
                int neiID=NeighborCon[NodeID][nei].first;
                if(CoreTag[neiID] != -1){
                    cout<<"Wrong! The contracted neighbor "<<neiID<<" of "<<NodeID<<" is not core vertex!!!"<<endl;
                    exit(1);
                }
                AdjaCoreMap[NodeID][neiID]=NeighborCon[NodeID][nei].second.first;//the root vertex is regarded as core vertex
                AdjaCoreMap[neiID][NodeID]=NeighborCon[NodeID][nei].second.first;//the root vertex is regarded as core vertex
            }
        }
    }
    vector<pair<int,int>> vecpair;
    vecpair.clear();
    AdjaCore.assign(nodenum,vecpair);
    for(int i=0;i<AdjaCoreMap.size();i++){
        if(AdjaCoreMap[i].empty())
            continue;

        for(map<int,int>::iterator it=AdjaCoreMap[i].begin();it!=AdjaCoreMap[i].end();it++){
            AdjaCore[i].push_back(make_pair((*it).first, (*it).second));
        }
//        if(AdjaGraph[i].empty()){
//            AdjaGraph[i] = AdjaCore[i];
//        }else{
//            cout<<"Wrong!! AdjaGraph["<<i<<"] is not empty!"<<endl;
//            exit(1);
//        }
    }

    /// for AdjaGraph: current version of AdjaGraph is not truly correct
//    AdjaGraph.clear();
//    AdjaGraph.resize(nodenum);
    map<int,map<int,int>> mapvec;
    mapvec.clear();
    SuppPartiID.assign(nodenum, mapvec);
    map<int,pair<int,set<int>>> mapset;
    mapset.clear();
    SuppPartiIDReal.assign(nodenum, mapset);
    PartiVertex.assign(partiNum,unordered_set<int>());//only insert in-partition vertex in this case

    for(int id=0;id<nodenum;++id){
        if(CoreTag[id]!=-1){//for periphery vertex
            PartiVertex[CoreTag[id]].insert(id);
        }
    }

//    cout<<"Vertex number of each periphery:";
//    vector<int> PartiVertexNum(partiNum,0);
//    double aveVnum = 0.0;
//    for(int i=0;i<partiNum;++i){
//        PartiVertexNum[i] = PartiVertex[i].size();
//        cout<<" "<<PartiVertexNum[i];
//        aveVnum += PartiVertexNum[i];
//    }
//    cout<<" ; Average vertex number: "<<aveVnum / partiNum<<endl;

    for(int PID=0;PID<partiNum;PID++){
        int childrank=cr[PID];
        for(int i=0;i<Tree[childrank].vert.size();i++){
            ID1 = Tree[childrank].vert[i].first;
//            PartiVertex[PID].insert(ID1);//insert interface vertex

            for(int j=i+1;j<Tree[childrank].vert.size();++j){
                ID2 = Tree[childrank].vert[j].first;
                if(SuppPartiID[ID1].find(ID2)==SuppPartiID[ID1].end()){//if we cannot find ID2
                    SuppPartiID[ID1][ID2]=map<int,int>();
                    SuppPartiID[ID2][ID1]=map<int,int>();
                }
                SuppPartiID[ID1][ID2].insert({PID,INF});
                SuppPartiID[ID2][ID1].insert({PID,INF});
                for(auto it=SCconNodesMT[ID1][ID2].begin();it!=SCconNodesMT[ID1][ID2].end();++it){
                    assert(AdjaCoreMap[ID1].find(ID2)!=AdjaCoreMap[ID1].end());//must exist
                    assert(AdjaCoreMap[ID2].find(ID1)!=AdjaCoreMap[ID2].end());//must exist
                    if(CoreTag[it->first] == PID){
                        if(SuppPartiID[ID1][ID2][PID] > it->second){
                            SuppPartiID[ID1][ID2][PID] = it->second;
                            SuppPartiID[ID2][ID1][PID] = it->second;
                        }

                        if(SuppPartiIDReal[ID1].find(ID2) == SuppPartiIDReal[ID1].end()){
                            SuppPartiIDReal[ID1][ID2]=make_pair(AdjaCoreMap[ID1][ID2],set<int>());
                            SuppPartiIDReal[ID2][ID1]=make_pair(AdjaCoreMap[ID1][ID2],set<int>());
                        }
                        if(AdjaCoreMap[ID1][ID2] == it->second){
                            SuppPartiIDReal[ID1][ID2].second.insert(CoreTag[it->first]);
                            SuppPartiIDReal[ID2][ID1].second.insert(CoreTag[it->first]);
                        }
                    }else{
                        CoreTag[it->first];
                    }

                }
            }
        }
    }
    //clear useless variables
    existCore.clear();
}

void Graph::TreeLabelCompute(pair<int,int> pidRange, vector<int> & pidRanks)
{
    int rankRoot;
    int ID1,ID2,temp_dis;
    for(int pid=pidRange.first;pid<pidRange.second;++pid){
        rankRoot = pidRanks[pid];

        vector<int> ancestors; //ancestors, the virtual root node is omited
        ancestors.push_back(Tree[0].uniqueVertex);
        ancestors.push_back(Tree[rankRoot].uniqueVertex);
        vector<int> interfaces;//interfaces
        Tree[rankRoot].pos.clear();
        Tree[rankRoot].pos.push_back(1);
        Tree[rankRoot].dis.push_back(INF);
        Tree[rankRoot].dis.push_back(0);
        ID1 = Tree[rankRoot].uniqueVertex;
        /// interface distance
        for(int i=0;i<Tree[rankRoot].vert.size();++i){
            ID2 = Tree[rankRoot].vert[i].first;
//            temp_dis = QueryCore(ID1, ID2);
            temp_dis = Tree[rankRoot].vert[i].second.first;
//            Tree[rankRoot].disInf.push_back(temp_dis);
            Tree[rankRoot].disInf.insert({ID2,temp_dis});
            interfaces.push_back(ID2);
            Tree[rankRoot].FNInf.insert({ID2,true});
            //correctness check
//                    int d1=Tree[rankRoot].disInf[i];
//                    int d2= Dijkstra(ID1,ID2,Neighbor);
//                    if(d1!=d2){
//                        cout<<"Incorrect! "<<ID1<<" "<<ID2<<": "<<d1<<" "<<d2<<endl;
//                    }
        }
        //get the all-pair distance among the interface vertices
        map<int,unordered_map<int,int>> disInfs;//distance from interface vertex u to another interface vertex v
        for(int i=0;i<Tree[rankRoot].vert.size();++i){
            ID1 = Tree[rankRoot].vert[i].first;
            for(int j=i+1;j<Tree[rankRoot].vert.size();++j){
                ID2 = Tree[rankRoot].vert[j].first;
//                int dis = QueryCore(ID1,ID2);
                int dis = INF;// dis might be INF
                if(AdjaCoreMap[ID1].find(ID2) != AdjaCoreMap[ID1].end()){//if found
                    dis = AdjaCoreMap[ID1][ID2];
                }
//                for(auto it=NeighborConCore[ID1].begin();it!=NeighborConCore[ID1].end();++it){
//                    if(it->first == ID2){
//                        dis = it->second;
//                        break;
//                    }
//                }
                if(dis == 0 ){
                    cout<<"Wrong! "<<ID1<<" "<<ID2<<" "<<dis<<endl; exit(1);
                }
                if(ID1 <= ID2){
                    disInfs[ID1].insert({ID2, dis});
                }else{
                    disInfs[ID2].insert({ID1, dis});
                }
            }
        }

        for(int i=0;i<Tree[rankRoot].ch.size();i++){
            makeIndexDFSCore(Tree[rankRoot].ch[i],ancestors,interfaces,disInfs);
        }

    }
}

//Function of tree-label index construction
void Graph::Compute_tree_label(bool ifParallel){
    printf( "Computing Tree Label...\n" );
    std::chrono::high_resolution_clock::time_point t10;
    std::chrono::high_resolution_clock::time_point t11;
    std::chrono::duration<double> time_span1;

    t10=std::chrono::high_resolution_clock::now();

    makeRMQCore();//all the root vertices is connected to the virtual root.

    int sum_r = 0;
    int ID1,ID2;
    int temp_dis;
    if(ifParallel){
        vector<int> pidRanks;
        for(int pid=0;pid<Tree[0].ch.size();++pid) {
            int rankRoot = Tree[0].ch[pid];
            if (Tree[rankRoot].ch.size() > 0) {
                pidRanks.emplace_back(rankRoot);
            }
        }
        assert(partiNum == pidRanks.size());
        if(partiNum>threadnum){
            int step=partiNum/threadnum;
            boost::thread_group threadf;
            for(int i=0;i<threadnum;i++){
                pair<int,int> p;
                p.first=i*step;
                if(i==threadnum-1)
                    p.second=partiNum;
                else
                    p.second=(i+1)*step;
                threadf.add_thread(new boost::thread(&Graph::TreeLabelCompute, this, p, pidRanks));
            }
            threadf.join_all();
        }else{
            boost::thread_group threadf;
            for(int pid=0;pid<partiNum;++pid) {
                threadf.add_thread(new boost::thread(&Graph::TreeLabelCompute, this, make_pair(pid,pid+1), pidRanks));
            }
            threadf.join_all();
        }
    }
    else{
        for(int pid=0;pid<Tree[0].ch.size();++pid){
            int rankRoot = Tree[0].ch[pid];
            if(Tree[rankRoot].ch.size() > 0){
                //initialize
                vector<int> ancestors; //ancestors, the virtual root node is omited
                ancestors.push_back(Tree[0].uniqueVertex);
                ancestors.push_back(Tree[rankRoot].uniqueVertex);
                vector<int> interfaces;//interfaces
                Tree[rankRoot].pos.clear();
                Tree[rankRoot].pos.push_back(1);
                Tree[rankRoot].dis.push_back(INF);
                Tree[rankRoot].dis.push_back(0);
                ID1 = Tree[rankRoot].uniqueVertex;
                /// interface distance
                for(int i=0;i<Tree[rankRoot].vert.size();++i){
                    ID2 = Tree[rankRoot].vert[i].first;
//                    temp_dis = QueryCore(ID1, ID2);//query on core label
                    temp_dis = Tree[rankRoot].vert[i].second.first;
                    Tree[rankRoot].disInf.insert({ID2,temp_dis});
                    interfaces.push_back(ID2);
                    Tree[rankRoot].FNInf.insert({ID2,true});
                    //correctness check
//                    int d1=Tree[rankRoot].disInf[i];
//                    int d2= Dijkstra(ID1,ID2,Neighbor);
//                    if(d1!=d2){
//                        cout<<"Incorrect! "<<ID1<<" "<<ID2<<": "<<d1<<" "<<d2<<endl;
//                    }
                }
                //get the all-pair distances among the interface vertices, which is equal to the shortcuts among them
                map<int,unordered_map<int,int>> disInfs;//distance from interface vertex u to another interface vertex v
                for(int i=0;i<Tree[rankRoot].vert.size();++i){
                    ID1 = Tree[rankRoot].vert[i].first;
                    for(int j=i+1;j<Tree[rankRoot].vert.size();++j){
                        ID2 = Tree[rankRoot].vert[j].first;
//                        int dis = QueryCore(ID1,ID2);
                        int dis = INF;// dis might be INF
                        if(AdjaCoreMap[ID1].find(ID2) != AdjaCoreMap[ID1].end()){//if found
                            dis = AdjaCoreMap[ID1][ID2];
                        }

                        if(dis == 0 ){
                            cout<<"Wrong! "<<ID1<<" "<<ID2<<" "<<dis<<endl; exit(1);
                        }
                        if(ID1 <= ID2){
                            disInfs[ID1].insert({ID2, dis});
                        }else{
                            disInfs[ID2].insert({ID1, dis});
                        }
                    }
                }
                assert(Tree[rankRoot].disInf.size()==interfaces.size());
                for(int i=0;i<Tree[rankRoot].ch.size();i++){
                    makeIndexDFSCore(Tree[rankRoot].ch[i],ancestors,interfaces,disInfs);
                }
            }
        }
    }

    t11=std::chrono::high_resolution_clock::now();
    time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
    double runT1 = time_span1.count();
    cout<<"Partition's label construction time: "<<runT1<<" s"<<endl;
}


//function of erasing edge (u,v), i.e., erase u from v's adjacency list.
void Graph::deleteECore(int u,int v){
//	if(Emap[u].find(v)!=Emap[u].end()){
//		Emap[u].erase(Emap[u].find(v));
//		DD[u]--;
//	}

	if(E[v].find(u)!=E[v].end()){
		E[v].erase(E[v].find(u));
		DD[v]--;
	}
}

//function of inserting edge (u,v)
void Graph::insertECore(int u,int v,int w){
	if(E[u].find(v)==E[u].end()){//if not found
		E[u].insert(make_pair(v,make_pair(w,1)));
		DD[u]++;
		DD2[u]++;
	}
	else{//if found
		if(E[u][v].first>w)
			E[u][v]= make_pair(w,1);
        else if(E[u][v].first==w)
            E[u][v].second+=1;
	}

	if(E[v].find(u)==E[v].end()){
        E[v].insert(make_pair(u,make_pair(w,1)));
		DD[v]++;
		DD2[v]++;
	}
	else{
        if(E[v][u].first>w)
            E[v][u]=make_pair(w,1);
        else if(E[v][u].first==w)
            E[v][u].second+=1;
	}
}

//
int Graph::matchCore(int x,vector<pair<int,pair<int,int>>> &vert){
	int nearest=vert[0].first;
	for(int i=1;i<vert.size();i++){
		if(NodeOrder[vert[i].first]<NodeOrder[nearest])//get the least node order
			nearest=vert[i].first;
	}
	if(existCore[nearest])//if it exists in core, i.e., the pa of root vertex is 0
		return 0;
	else{
		//cout<<nearest<<" "<<rankCore[nearest]<<endl;
		return rank[nearest];
	}

}

/*void Graph::makeTree(){
	//rank.assign(HighestOrder+2,0);
	rankCore.clear();
	rankCore.assign(nodenum,-1);
	int len=HighestOrder-1;

	NodeCore root;
	root.uniqueVertex=-1;
	root.height=1;
	Tree.push_back(root);

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

		Tree[pa].ch.push_back(Tree.size());
		nod.pa=pa;
		nod.height=Tree[pa].height+1;
		rankCore[x]=Tree.size();//the position of tree

		if(pa==0)
			nod.treeroot=rankCore[x];
		else
			nod.treeroot=Tree[pa].treeroot;

		Tree.push_back(nod);
	}
}*/




