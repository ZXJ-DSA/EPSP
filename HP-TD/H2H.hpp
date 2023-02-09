/*
 * Labelcon.cpp
 *
 *  Created on: 22 Dec 2020
 *      Author: zhangmengxuan, Xinjie ZHOU
 */
#include "gtreeTD.h"

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

//struct DegComp2{
//    int x;
//    DegComp2(int _x){
//        x=_x;
//    }
//    bool operator< (const DegComp2 d) const{
//        if(_DD_[x]!=_DD_[d.x])
//            return _DD_[x]<_DD_[d.x];
//        if(_DD2_[x]!=_DD2_[x])
//            return _DD2_[x]<_DD2_[d.x];
//        return x<d.x;
//    }
//};

void Gtree::H2HconOrderMT(){
//	CHconsorderMT(orderfile);
//    CHconsMTOrderGenerate(orderfile);
//    exit(0);
    MDEContract();//MDE-based contraction
//    LayeredMDEContract();
//    Contraction();
//    exit(0);
    cout<<"Done."<<endl;
    Timer tt;
    tt.start();
    makeTree();
    makeIndex();
    tt.stop();
    cout<<"The time used for index construction: "<<tt.GetRuntime()<<" s."<<endl;
}

void Gtree::getPartiGraph(){
    map<int,unordered_set<int>> partiVSet;
    vector<int> leafNodes;
    for(int i=0;i<GTree.size();++i){
        if(GTree[i].isleaf){
            leafNodes.push_back(i);

            unordered_set<int> tempSet; tempSet.clear();
            tempSet.insert(GTree[i].leafnodes.begin(), GTree[i].leafnodes.end());//get the partition vertex set
            partiVSet.insert({i,tempSet});
        }
    }
    int pID,ID;
    int nid,wei;
    for(int i=0;i<leafNodes.size();++i){
        pID = leafNodes[i];
        PartiGraph.insert({pID,unordered_map<int,vector<pair<int,int>>>()});
        for(int j=0;j<GTree[pID].leafnodes.size();++j){
            ID = GTree[pID].leafnodes[j];
            vector<pair<int,int>> adjs;
            for(int k=0;k<Nodes[ID].adjnodes.size();++k){
                nid = Nodes[ID].adjnodes[k];
                if(partiVSet[pID].find(nid) != partiVSet[pID].end()){//if found
                    wei = Nodes[ID].adjweight[k];
                    adjs.emplace_back(nid,wei);
                }
            }
            PartiGraph[pID].insert({ID,adjs});
        }
    }
}

void Gtree::getOverlayOrder(){
    OverlayGraph.clear();
    NodeOrder.assign(node_num,-1);
    vNodeOrder.assign(node_num,-1);
    vector<int> HighToLow;
    BoundaryTag.assign(node_num, false);

    int curID=0;
    int curLevel=0;
    //get boundary order by BFS
//    getBOrder(curID,HighToLow);
    queue<pair<int,int>> queue1;
    queue1.push(make_pair(curID,0));
    boundaryLevel.clear();
    pair<int,int> topPair;
    vector<int> currentB;

    while(!queue1.empty()){
        topPair = queue1.front();
        curID = topPair.first; int tempLevel = topPair.second;
        queue1.pop();
        if(curLevel != tempLevel){
            boundaryLevel.push_back(currentB);
            currentB.clear();
            curLevel = tempLevel;
        }
        if(!GTree[curID].children.empty()){

            for(int i=0;i<GTree[curID].children.size();++i){
                int pID = GTree[curID].children[i];
                queue1.push(make_pair(pID,tempLevel+1));

                for(int j=0;j<GTree[pID].borders.size();++j){
                    int ID = GTree[pID].borders[j];
                    if(Nodes[ID].isborder && NodeOrder[ID] == -1){//if ID is boundary vertex and it has never been assigned order
                        BoundaryTag[ID] = true;

                        HighToLow.push_back(ID);
                        NodeOrder[ID] = node_num - HighToLow.size();
                        vNodeOrder[NodeOrder[ID]] = ID;
                        OverlayGraph.insert({ID,vector<pair<int,int>>()});
                        currentB.push_back(ID);

                    }
                }
            }

        }
    }

    cout<<"Number of layers: "<<boundaryLevel.size()<<endl;
    cout<<"Node number of overlay graph: "<<OverlayGraph.size()<<endl;

//    if()

//    for(int i=0;i<node_num;++i){
//        if(NodeOrder[i] == -1){
//            HighToLow.push_back(i);
//            NodeOrder[i] = node_num - HighToLow.size();
//            vNodeOrder[NodeOrder[i]] = i;
//        }
//    }
//    assert(HighToLow.size() == node_num);

}

/*void Gtree::getBOrder(int curID,vector<int> & HighToLow){
    int pID,ID;
    if(!GTree[curID].children.empty()){
        for(int i=0;i<GTree[curID].children.size();++i){
            pID = GTree[curID].children[i];

            for(int j=0;j<GTree[pID].borders.size();++j){
                ID = GTree[pID].borders[j];
                if(Nodes[ID].isborder && NodeOrder[ID] == -1){//if ID is boundary vertex and it has never been assigned order
                    HighToLow.push_back(ID);
                    NodeOrder[ID] = node_num - HighToLow.size();
                    vNodeOrder[NodeOrder[ID]] = ID;
                    OverlayGraph.insert({ID,vector<pair<int,int>>()});
                }
            }

            getBOrder(pID,HighToLow);//DFS search
        }
    }
}*/

//function of deleting u from v's neighbor
void Gtree::deleteE(int u,int v){
//    if(E[u].find(v)!=E[u].end()){
//        E[u].erase(E[u].find(v));
//        //DD[u]--;
//    }

    if(E[v].find(u)!=E[v].end()){
        E[v].erase(E[v].find(u));
        //DD[v]--;
    }
}
//function of inserting u to v's neighbor and verse vice.
void Gtree::insertE(int u,int v,int w){
    if(E[u].find(v)==E[u].end()){
        E[u].insert(make_pair(v,make_pair(w,1)));
        //DD[u]++;
        //DD2[u]++;
    }
    else{
        if(E[u][v].first>w)
            E[u][v]=make_pair(w,1);
        else if(E[u][v].first==w)
            E[u][v].second+=1;
    }

    if(E[v].find(u)==E[v].end()){
        E[v].insert(make_pair(u,make_pair(w,1)));
        //DD[v]++;
        //DD2[v]++;
    }
    else{
        if(E[v][u].first>w)
            E[v][u]=make_pair(w,1);
        else if(E[v][u].first==w)
            E[v][u].second+=1;
    }
}

void Gtree::insertEMTorder(int u,int v,int w){
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



//Function of contracting vertices by MDE
void Gtree::MDEContract(){
    Timer tt;
    tt.start();
    cout<<"Start contraction..."<<endl;

    //for H2H update
    map<int, vector<int>> mi;
    SCconNodesMT.assign(node_num, mi);

    _DD_.assign(node_num,0);_DD2_.assign(node_num,0);
//    DD.assign(node_num,0); DD2.assign(node_num,0);

    set<DegComp1> Deg;//min first
//    vector<bool> active(node_num,false);//flag that indicate whether a vertex is active for contraction
    int degree;
    vector<bool> exist;
    exist.assign(node_num,false);//if in the core, all vertices is originally in core

    for(auto it=OverlayGraph.begin();it!=OverlayGraph.end();++it){
        int ID=it->first;
        degree=E[ID].size();
        exist[ID] = true;
        if(degree > 0){//get degree
            _DD_[ID]=degree;
            _DD2_[ID]=degree;
//            DD[i]=degree;
//            DD2[i]=degree;
            Deg.insert(DegComp1(ID));
//            active[i] = true;
        }else{
            cout<<"Wrong!! Degree of "<<ID<<" is "<<degree<<endl;
            exit(1);
        }
    }

    vNodeOrder.clear();
    while(vNodeOrder.size()<node_num-OverlayGraph.size()){
        vNodeOrder.push_back(-1);
    }
    vector<bool> change;
    change.assign(node_num,false);//whether the neighbor (degree) has changed

    vector<pair<int,pair<int,int>>> vect;
    NeighborCon.assign(node_num,vect);//temporal graph to store Neighbors in the core, for graph contraction

    bool CutLabel=false;
    int count=0;
    int ID1,ID2;

//    cout<<"192402: "<<E[192402].size()<<" "<<Nodes[192402].inleaf<<" "<<Nodes[192402].isborder<<endl;
//    for(auto i1=E[192402].begin();i1!=E[192402].end();++i1){
//        cout<<"E "<<i1->first<<" "<<Nodes[i1->first].inleaf<<" "<<Nodes[i1->first].isborder<<endl;;
//    }
//    for(auto i1=Nodes[192402].adjnodes.begin();i1!=Nodes[192402].adjnodes.end();++i1){
//        cout<<"Adj "<<*i1<<" "<<Nodes[*i1].inleaf<<" "<<Nodes[*i1].isborder<<endl;
//    }

    //Get the order of all vertices by MDE
    while(!Deg.empty()){
        count+=1;
        int x=(*Deg.begin()).x;//minimum degree first

        while(change[x]){//update the degree if it is changed
            Deg.erase(DegComp1(x));
            _DD_[x]=E[x].size();
            _DD2_[x]=E[x].size();
            Deg.insert(DegComp1(x));
            change[x]=false;
            x=(*Deg.begin()).x;
        }
//        cout<<x<<" "<<E[x].size()<<endl;
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
            if(exist[(*it).first]){
                Neigh.push_back(*it);
            }
        }
        if(Neigh.empty() && !Deg.empty()){
            cout<<"!!! Neigh is empty for "<<x<<" : "<<E[x].size()<<" "<<Nodes[x].adjnodes.size()<<" "<<count<<endl;
//            exit(1);
        }

        NeighborCon[x].assign(Neigh.begin(),Neigh.end());


        exist[x]=false;
        //delete the star
        for(int i=0;i<Neigh.size();i++){
            int y=Neigh[i].first;
            deleteE(x,y);//delete x from y's adjacency list
            change[y]=true;
        }
        //add all-pair neighbors
        for(int i=0;i<Neigh.size();i++) {
            ID1 = Neigh[i].first;
            for (int j = i + 1; j < Neigh.size(); j++) {
                ID2 = Neigh[j].first;
                insertE(ID1, ID2, Neigh[i].second.first + Neigh[j].second.first);
                /// For TD update

                if(Neigh[i].first<Neigh[j].first)
                    SCconNodesMT[Neigh[i].first][Neigh[j].first].push_back(x);//no direction
                else if(Neigh[j].first<Neigh[i].first)
                    SCconNodesMT[Neigh[j].first][Neigh[i].first].push_back(x);

                change[ID1] = true;
                change[ID2] = true;
            }
        }

    }
    tt.stop();
    cout<<"The time used for contraction: "<<tt.GetRuntime()<<" s."<<endl;

    NodeOrder.assign(node_num,-1);
    for(int k=0;k<vNodeOrder.size();k++){
        NodeOrder[vNodeOrder[k]]=k;
    }
    cout << "???" << endl;
//    while(vNodeOrder.size()<node_num){
//        vNodeOrder.push_back(-1);
//    }

}

//Function of contracting vertices by MDE
void Gtree::LayeredMDEContract(){
    Timer tt;
    tt.start();
    cout<<"Start contraction..."<<endl;

    //for H2H update
    map<int, vector<int>> mi;
    SCconNodesMT.assign(node_num, mi);

    _DD_.assign(node_num,0);_DD2_.assign(node_num,0);
//    DD.assign(node_num,0); DD2.assign(node_num,0);

    set<DegComp1> Deg;//min first
//    vector<bool> active(node_num,false);//flag that indicate whether a vertex is active for contraction
    int degree;
    vector<bool> exist;
    exist.assign(node_num,false);//if in the core, all vertices is originally in core

    vNodeOrder.clear();

    vector<bool> change;
    change.assign(node_num,false);//whether the neighbor (degree) has changed

    vector<pair<int,pair<int,int>>> vect;
    NeighborCon.assign(node_num,vect);//temporal graph to store Neighbors in the core, for graph contraction

    int count=0;
    int ID1,ID2;

    while(vNodeOrder.size()<node_num-OverlayGraph.size()){
        vNodeOrder.push_back(-1);
    }
    exist = BoundaryTag;
//    for(auto it=OverlayGraph.begin();it!=OverlayGraph.end();++it){
//        int ID=it->first;
//        exist[ID] = true;
//    }

    for ( int i = boundaryLevel.size() - 1; i >= 0; i-- ) {//start from the lowest level
        Deg.clear();

        cout<<"Layer "<<i<<" , size: "<<boundaryLevel[i].size()<<endl;


        for (int j = 0; j < boundaryLevel[i].size(); j++) {//for each partition in this level
            int ID = boundaryLevel[i][j];
            degree=E[ID].size();
            if(degree > 0) {//get degree
                _DD_[ID] = degree;
                _DD2_[ID] = degree;
                Deg.insert(DegComp1(ID));
//                exist[ID] = true;
            }else{
                cout<<"Wrong!! Degree of "<<ID<<" is "<<degree<<endl;
                exit(1);
            }
        }
        //Get the order of all vertices by MDE
        while(!Deg.empty()){
            count+=1;
            int x=(*Deg.begin()).x;//minimum degree first

            while(change[x]){//update the degree if it is changed
                Deg.erase(DegComp1(x));
                _DD_[x]=E[x].size();
                _DD2_[x]=E[x].size();
                Deg.insert(DegComp1(x));
                change[x]=false;
                x=(*Deg.begin()).x;
            }

            vNodeOrder.push_back(x);//least important vertex first
            Deg.erase(Deg.begin());

            vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();
            for(auto it=E[x].begin();it!=E[x].end();it++){
                if(exist[(*it).first]){
                    Neigh.push_back(*it);
                }
            }
            NeighborCon[x].assign(Neigh.begin(),Neigh.end());


            exist[x]=false;
            //delete the star
            for(int i=0;i<Neigh.size();i++){
                int y=Neigh[i].first;
                deleteE(x,y);//delete x from y's adjacency list
                change[y]=true;
            }
            //add all-pair neighbors
            for(int i=0;i<Neigh.size();i++) {
                ID1 = Neigh[i].first;
                for (int j = i + 1; j < Neigh.size(); j++) {
                    ID2 = Neigh[j].first;
                    insertE(ID1, ID2, Neigh[i].second.first + Neigh[j].second.first);
                    /// For TD update

                    if(Neigh[i].first<Neigh[j].first)
                        SCconNodesMT[Neigh[i].first][Neigh[j].first].push_back(x);//no direction
                    else if(Neigh[j].first<Neigh[i].first)
                        SCconNodesMT[Neigh[j].first][Neigh[i].first].push_back(x);

                    change[ID1] = true;
                    change[ID2] = true;
                }
            }

        }
    }




    tt.stop();
    cout<<"The time used for contraction: "<<tt.GetRuntime()<<" s."<<endl;

    NodeOrder.assign(node_num,-1);
    for(int k=node_num-OverlayGraph.size();k<vNodeOrder.size();k++){
        NodeOrder[vNodeOrder[k]]=k;
    }



}

//Function of contraction vertices by given order
void Gtree::Contraction(){
    Timer tt;
    tt.start();
    cout<<"Start contraction..."<<endl;
    vector<pair<int,pair<int,int>>> vect;
    NeighborCon.assign(node_num,vect);

    int stepShow = ceil(node_num/10000)*100;
    cout<<"Step for show: "<<stepShow<<endl;

    map<int, vector<int>> mi;
    SCconNodesMT.assign(node_num, mi);
    vector<bool> exist; exist.assign(node_num,false);
    for(auto it=OverlayGraph.begin();it!=OverlayGraph.end();++it){
        exist[it->first] = true;
    }
    //cout<<"Begin to contract"<<endl;
    for(int nodeorder=0;nodeorder<node_num;nodeorder++){//start from the most important vertex
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

            if(nodeorder % stepShow == 0){
                cout<<"node order: "<<nodeorder<<" ; neighbor size: "<<Neigh.size()<<endl;
            }

            for(int i=0;i<Neigh.size();i++){
                int y=Neigh[i].first;
                deleteE(x,y);//delete x from y's neighbor
                //change[y]=true;
            }

            if(Neigh.size()<=100){
                //single thread
                for(int i=0;i<Neigh.size();i++){
                    for(int j=i+1;j<Neigh.size();j++){
                        insertE(Neigh[i].first,Neigh[j].first,Neigh[i].second.first+Neigh[j].second.first);
                        if(Neigh[i].first<Neigh[j].first)
                            SCconNodesMT[Neigh[i].first][Neigh[j].first].push_back(x);//no direction
                        else if(Neigh[j].first<Neigh[i].first)
                            SCconNodesMT[Neigh[j].first][Neigh[i].first].push_back(x);
                    }
                }
            }else{
                //multiple thread
                int step=Neigh.size()/thread_num;
                boost::thread_group thread;
                for(int i=0;i<thread_num;i++){
                    pair<int,int> p;
                    p.first=i*step;
                    if(i==thread_num-1)
                        p.second=Neigh.size();
                    else
                        p.second=(i+1)*step;
                    thread.add_thread(new boost::thread(&Gtree::NeighborComorder, this, boost::ref(Neigh), p, x));
                }
                thread.join_all();
            }

        }
    }
    tt.stop();
    cout<<"The time used for contraction: "<<tt.GetRuntime()<<" s."<<endl;
}

void Gtree::NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x){
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

void Gtree::makeTree(){
    vector<int> vecemp; //vecemp.clear();
    VidtoTNid.assign(node_num,vecemp);

    rank.assign(node_num,0);
    //Tree.clear();
    int len=vNodeOrder.size()-1;
    heightMax=0;

    TDNode rootn;
    int x=vNodeOrder[len];
    //cout<<"len "<<len<<" , ID "<<x<<endl;
    while(x==-1){//to skip those vertices whose ID is -1
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
//        if(len%100==0){
//            cout<<"len: "<<len<<endl;
//        }
        x=vNodeOrder[len];
        if(x==-1){
            continue;
        }
        TDNode nod;
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
        if(nod.height>heightMax)
            heightMax=nod.height;
        rank[x]=Tree.size();
        Tree.push_back(nod);
//        cout<<"len "<<len<<" , ID "<<x<<" ; "<<NeighborCon[x].size()<<endl;
//        if(len == 263982){
//            cout<<len<<" "<<NeighborCon[x].size()<<endl;
//        }
    }
    cout<<"Finish make tree."<<endl;
}

void Gtree::makeIndex(){
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

int Gtree::match(int x,vector<pair<int,pair<int,int>>> &vert) {
    if(vert.empty()){
        cout<<"Empty vert! "<<x<<endl;
        exit(1);
    }
    int nearest = vert[0].first;
    for (int i = 1; i < vert.size(); i++) {
        if (rank[vert[i].first] > rank[nearest])
            nearest = vert[i].first;
    }
    int p = rank[nearest];
    return p;
}

void Gtree::makeRMQDFS(int p, int height){
    toRMQ[p] = EulerSeq.size();
    EulerSeq.push_back(p);
    for (int i = 0; i < Tree[p].ch.size(); i++){
        makeRMQDFS(Tree[p].ch[i], height + 1);
        EulerSeq.push_back(p);
    }
}

void Gtree::makeRMQ(){
    //EulerSeq.clear();
    toRMQ.assign(node_num,0);
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

int Gtree::LCAQuery(int _p, int _q){
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

void Gtree::makeIndexDFS(int p, vector<int>& list){
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

void Gtree::H2HdecBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch) {
    map<int, int> checkedDis;

    for (int i = 0; i < Tree.size(); i++) {
        Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
    }

//NodeOrderss.clear();
    NodeOrderss.assign(NodeOrder.begin(), NodeOrder.end());
    vector<set<int>> SCre; //SCre.clear();
    set<int> ss; //ss.clear();
    SCre.assign(node_num, ss);//{vertexID, set<int>}
    set<OrderCompp> OC; //OC.clear();//vertexID in decreasing node order

    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed

    int a, b, oldW, newW, lid, hid;
    for (int k = 0; k < wBatch.size(); k++) {
        a = wBatch[k].first.first;
        b = wBatch[k].first.second;
        oldW = wBatch[k].second.first;
        newW = wBatch[k].second.second;
        if (NodeOrder[a] < NodeOrder[b]) {
            lid = a;
            hid = b;
        } else {
            lid = b;
            hid = a;
        }

        for (int i = 0; i < Nodes[a].adjnodes.size(); i++) {
            if (Nodes[a].adjnodes[i] == b) {
                Nodes[a].adjweight[i] = newW;
                break;
            }
        }
        for (int i = 0; i < Nodes[b].adjnodes.size(); i++) {
            if (Nodes[b].adjnodes[i] == a) {
                Nodes[b].adjweight[i] = newW;
                break;
            }
        }

        for (int i = 0; i < Tree[rank[lid]].vert.size(); i++) {
            if (Tree[rank[lid]].vert[i].first == hid) {
                if (Tree[rank[lid]].vert[i].second.first > newW) {
                    Tree[rank[lid]].vert[i].second.first = newW;
                    Tree[rank[lid]].vert[i].second.second = 1;
                    SCre[lid].insert(hid);
                    OC.insert(OrderCompp(lid));
                } else if (Tree[rank[lid]].vert[i].second.first == newW) {
                    Tree[rank[lid]].vert[i].second.second += 1;
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
    while (!OC.empty()) {
        ProID = (*OC.begin()).x;
        OC.erase(OC.begin());
        vector<pair<int, pair<int, int>>> Vert = Tree[rank[ProID]].vert;
        bool ProIDdisCha = false;//to see if the distance labeling of proID change or not
        for (auto it = SCre[ProID].begin(); it != SCre[ProID].end(); it++) {
            int Cid = *it;
            int Cw;
            int cidH = Tree[rank[Cid]].height - 1;

            map<int, int> Hnei; //Hnei.clear();
            vector<pair<int, int>> Lnei; //Lnei.clear();
            for (int j = 0; j < Vert.size(); j++) {
                if (NodeOrder[Vert[j].first] > NodeOrder[Cid]) {
                    Hnei[Vert[j].first] = Vert[j].second.first;
                } else if (NodeOrder[Vert[j].first] < NodeOrder[Cid]) {
                    Lnei.push_back(make_pair(Vert[j].first, Vert[j].second.first));
                } else {
                    Cw = Vert[j].second.first;
                }
            }

            if (Tree[rank[ProID]].dis[cidH] > Cw) {
                Tree[rank[ProID]].dis[cidH] = Cw;
                Tree[rank[ProID]].FN[cidH] = true;
                ProIDdisCha = true;
                Tree[rank[ProID]].DisRe.insert(Cid);
            } else if (Tree[rank[ProID]].dis[cidH] == Cw) {
                Tree[rank[ProID]].FN[cidH] = true;
            }

            int hid, hidHeight, lid, lidHeight, wsum;
            for (int j = 0; j < Tree[rank[Cid]].vert.size(); j++) {
                hid = Tree[rank[Cid]].vert[j].first;
                hidHeight = Tree[rank[hid]].height - 1;
                if (Hnei.find(hid) != Hnei.end()) {
                    wsum = Cw + Hnei[hid];
                    if (wsum < Tree[rank[Cid]].vert[j].second.first) {
                        Tree[rank[Cid]].vert[j].second.first = wsum;
                        Tree[rank[Cid]].vert[j].second.second = 1;
                        SCre[Cid].insert(hid);
                        OC.insert(OrderCompp(Cid));
                    } else if (wsum == Tree[rank[Cid]].vert[j].second.first) {
                        Tree[rank[Cid]].vert[j].second.second += 1;
                    }

                }
            }
            for (int j = 0; j < Lnei.size(); j++) {
                lid = Lnei[j].first;
                lidHeight = Tree[rank[lid]].height - 1;
                for (int k = 0; k < Tree[rank[lid]].vert.size(); k++) {
                    if (Tree[rank[lid]].vert[k].first == Cid) {
                        wsum = Cw + Lnei[j].second;
                        if (Tree[rank[lid]].vert[k].second.first > wsum) {
                            Tree[rank[lid]].vert[k].second.first = wsum;
                            Tree[rank[lid]].vert[k].second.second = 1;
                            SCre[lid].insert(Cid);
                            OC.insert(OrderCompp(lid));
                        } else if (Tree[rank[lid]].vert[k].second.first == wsum) {
                            Tree[rank[lid]].vert[k].second.second += 1;
                        }

                        break;
                    }
                }
            }
        }

        if (ProIDdisCha) {//if the distance labeling is dectected changed
            vertexIDChL.insert(ProID);
            ProBeginVertexSetNew.clear();
            ProBeginVertexSetNew.reserve(ProBeginVertexSet.size() + 1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew = rank[ProID], r;
            for (int i = 0; i < ProBeginVertexSet.size(); i++) {
                r = rank[ProBeginVertexSet[i]];
                if (LCAQuery(rnew, r) != rnew) {
                    ProBeginVertexSetNew.push_back(ProBeginVertexSet[i]);
                }
            }
            ProBeginVertexSet = ProBeginVertexSetNew;
        }
    }

//cout<<"Finish bottom-up refresh"<<endl;
    for (int i = 0; i < ProBeginVertexSet.size(); i++) {
        ProBeginVertexID = ProBeginVertexSet[i];
        vector<int> linee; //linee.clear();
        linee.reserve(heightMax);
        int pachidd = Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
        while (Tree[rank[pachidd]].height > 1) {
            linee.insert(linee.begin(), pachidd);
            pachidd = Tree[Tree[rank[pachidd]].pa].uniqueVertex;
        }
        linee.insert(linee.begin(), pachidd);
        EachNodeProBDis5(rank[ProBeginVertexID], linee, vertexIDChL, checkedDis);
    }
//return checkedDis.size();
}

void Gtree::EachNodeProBDis5(int child,vector<int>& line,set<int>& vertexIDChL, map<int,int>& checkedDis){
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
        EachNodeProBDis5(Tree[child].ch[i], line, vertexIDChL,checkedDis);
    }
    line.pop_back();

}

void Gtree::H2HincBatMT(vector<pair<pair<int,int>,pair<int,int>>>& wBatch) {
    int checknum = 0;
    map<pair<int, int>, int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
//OCdis.clear();

//NodeOrderss.clear();
    NodeOrderss.assign(NodeOrder.begin(), NodeOrder.end());
    vector<set<int>> SCre; //SCre.clear();
    set<int> ss; //ss.clear();
    SCre.assign(node_num, ss);//{vertexID, set<int>}
    set<OrderCompp> OC;
    OC.clear();//vertexID in decreasing node order

    for (int k = 0; k < wBatch.size(); k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int oldW = wBatch[k].second.first;
        int newW = wBatch[k].second.second;

        if (oldW != newW) {
            for (int i = 0; i < Nodes[a].adjnodes.size(); i++) {
                if (Nodes[a].adjnodes[i] == b) {
                    Nodes[a].adjweight[i] = newW;
                    break;
                }
            }
            for (int i = 0; i < Nodes[b].adjnodes.size(); i++) {
                if (Nodes[b].adjnodes[i] == a) {
                    Nodes[b].adjweight[i] = newW;
                    break;
                }
            }

            int lid, hid;
            if (NodeOrder[a] < NodeOrder[b]) {
                lid = a;
                hid = b;
            } else {
                lid = b;
                hid = a;
            }

            for (int i = 0; i < Tree[rank[lid]].vert.size(); i++) {
                if (Tree[rank[lid]].vert[i].first == hid) {
                    if (Tree[rank[lid]].vert[i].second.first == oldW) {
                        Tree[rank[lid]].vert[i].second.second -= 1;
                        if (Tree[rank[lid]].vert[i].second.second < 1) {
                            OCdis[make_pair(lid, hid)] = oldW;
                            SCre[lid].insert(hid);
                            OC.insert(OrderCompp(lid));
                        }
                    }
                    break;
                }
            }
        }
    }

    vector<int> ProBeginVertexSet;
    ProBeginVertexSet.clear();
    vector<int> ProBeginVertexSetNew;
    bool influence;
    int ProID;
    vector<int> line;
    while (!OC.empty()) {
        ProID = (*OC.begin()).x;
        OC.erase(OC.begin());
        vector<pair<int, pair<int, int>>> Vert = Tree[rank[ProID]].vert;
        influence = false;

//each ProID corresponds to a line
        line.clear();
        line.reserve(heightMax);
        int pachid = ProID;
        while (Tree[rank[pachid]].height > 1) {
            line.insert(line.begin(), pachid);
            pachid = Tree[Tree[rank[pachid]].pa].uniqueVertex;
        }
        line.insert(line.begin(), pachid);

        for (auto it = SCre[ProID].begin(); it != SCre[ProID].end(); it++) {
            int Cid = *it;
            int Cw = OCdis[make_pair(ProID, Cid)];
            int cidH = Tree[rank[Cid]].height - 1;

            map<int, int> Hnei; //Hnei.clear();
            vector<pair<int, int>> Lnei; //Lnei.clear();
            for (int j = 0; j < Vert.size(); j++) {
                if (NodeOrder[Vert[j].first] > NodeOrder[Cid]) {
                    Hnei[Vert[j].first] = Vert[j].second.first;
                } else if (NodeOrder[Vert[j].first] < NodeOrder[Cid]) {
                    Lnei.push_back(make_pair(Vert[j].first, Vert[j].second.first));
                }
            }
//check the affected shortcuts
            int hid, lid;
            for (int j = 0; j < Tree[rank[Cid]].vert.size(); j++) {
                hid = Tree[rank[Cid]].vert[j].first;
                if (Hnei.find(hid) != Hnei.end()) {
                    if (Cw + Hnei[hid] == Tree[rank[Cid]].vert[j].second.first) {
                        Tree[rank[Cid]].vert[j].second.second -= 1;
                        if (Tree[rank[Cid]].vert[j].second.second < 1) {
                            SCre[Cid].insert(hid);
                            OC.insert(OrderCompp(Cid));
                            OCdis[make_pair(Cid, hid)] = Cw + Hnei[hid];
                        }
                    }
                }
            }
            for (int j = 0; j < Lnei.size(); j++) {
                lid = Lnei[j].first;
                for (int k = 0; k < Tree[rank[lid]].vert.size(); k++) {
                    if (Tree[rank[lid]].vert[k].first == Cid) {
                        if (Tree[rank[lid]].vert[k].second.first == Cw + Lnei[j].second) {
                            Tree[rank[lid]].vert[k].second.second -= 1;
                            if (Tree[rank[lid]].vert[k].second.second < 1) {
                                SCre[lid].insert(Cid);
                                OC.insert(OrderCompp(lid));
                                OCdis[make_pair(lid, Cid)] = Cw + Lnei[j].second;
                            }
                        }
                        break;
                    }
                }
            }


//before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
            if (Tree[rank[ProID]].FN[cidH]) {
                influence = true;
//higher than Cid
                for (int i = 0; i < cidH; i++) {
                    if (Tree[rank[ProID]].dis[i] == Cw + Tree[rank[Cid]].dis[i]) {
                        Tree[rank[ProID]].cnt[i] -= 1;
                    }
                }

//equal to Cid
                Tree[rank[ProID]].FN[cidH] = false;
                Tree[rank[ProID]].cnt[cidH] -= 1;

//lower than Cid
                for (int i = cidH + 1; i < Tree[rank[ProID]].dis.size(); i++) {
                    if (Tree[rank[ProID]].dis[i] == Cw + Tree[rank[line[i]]].dis[cidH]) {
                        Tree[rank[ProID]].cnt[i] -= 1;
                    }
                }
            }

//get the new value of shortcut
//	cout<<Cw<<" increase to ";
            Cw = INF;
            int countwt = 0;

            for (int i = 0; i < Nodes[ProID].adjnodes.size(); i++) {///Neighbor
                if (Nodes[ProID].adjnodes[i] == Cid) {
                    Cw = Nodes[ProID].adjweight[i];//the weight value in the original graph
                    countwt = 1;
                    break;
                }
            }

            int ssw, wtt, wid;
            vector<int> Wnodes;
            Wnodes.clear();

            if (ProID < Cid)
                Wnodes = SCconNodesMT[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
            else
                Wnodes = SCconNodesMT[Cid][ProID];
            if (Wnodes.size() > 0) {
                for (int i = 0; i < Wnodes.size(); i++) {
                    wid = Wnodes[i];
                    for (int j = 0; j < Tree[rank[wid]].vert.size(); j++) {
                        if (Tree[rank[wid]].vert[j].first == ProID) {
                            ssw = Tree[rank[wid]].vert[j].second.first;
                        }
                        if (Tree[rank[wid]].vert[j].first == Cid) {
                            wtt = Tree[rank[wid]].vert[j].second.first;
                        }
                    }

                    if (ssw + wtt < Cw) {
                        Cw = ssw + wtt;
                        countwt = 1;
                    } else if (ssw + wtt == Cw) {
                        countwt += 1;
                    }
                }
            }

//cout<<Cw<<endl;
//refresh the shortcut to the new value
            for (int i = 0; i < Tree[rank[ProID]].vert.size(); i++) {
                if (Tree[rank[ProID]].vert[i].first == Cid) {
                    Tree[rank[ProID]].vert[i].second.first = Cw;
                    Tree[rank[ProID]].vert[i].second.second = countwt;
                    break;
                }
            }
        }

        if (influence) {
            ProBeginVertexSetNew.clear();
            ProBeginVertexSetNew.reserve(ProBeginVertexSet.size() + 1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew = rank[ProID], r;
            for (int i = 0; i < ProBeginVertexSet.size(); i++) {
                r = rank[ProBeginVertexSet[i]];
                if (LCAQuery(rnew, r) != rnew) {
                    ProBeginVertexSetNew.push_back(ProBeginVertexSet[i]);
                }
            }
            ProBeginVertexSet = ProBeginVertexSetNew;
        }

    }

    int ProBeginVertexID;
    for (int i = 0; i < ProBeginVertexSet.size(); i++) {
        ProBeginVertexID = ProBeginVertexSet[i];
        vector<int> linee; //linee.clear();
        linee.reserve(heightMax);
        int pachidd = Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
        while (Tree[rank[pachidd]].height > 1) {
            linee.insert(linee.begin(), pachidd);
            pachidd = Tree[Tree[rank[pachidd]].pa].uniqueVertex;
        }
        linee.insert(linee.begin(), pachidd);

        eachNodeProcessIncrease1(rank[ProBeginVertexID], linee, checknum);
    }
//return checknum;
}

void Gtree::eachNodeProcessIncrease1(int children, vector<int>& line, int& changelabel){
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
        eachNodeProcessIncrease1(Tree[children].ch[i],line,changelabel);
    }
    line.pop_back();
}

int Gtree::QueryH2H(int ID1,int ID2){
    if(ID1==ID2) return 0;
    if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1)
        return INF;
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

int Gtree::QueryPartiBoundary(int ID1,int ID2){
    int d=INF;
    assert(BoundaryTag[ID1]==false && BoundaryTag[ID2]==true);
    int PID1=Nodes[ID1].inleaf;
    int PosID1=Nodes[ID1].inleafpos;
    vector<int> B1;
    int bID1,d1;
    int tempdis;
    B1 = GTree[PID1].borders;

    map<int,int> m1  ;
    m1.clear();

    int lnNum1 = GTree[PID1].leafnodes.size();

    //get distances to corresponding boundary vertex
    for(int i=0;i<B1.size();++i){
        bID1 = B1[i];
        d1 = GTree[PID1].mind[i*lnNum1 + PosID1];
        m1.insert({bID1,d1});
    }

    int b1f,b2f,d1f,d2f,dbb,dbbf;
    //get distances between boundary vertex
    for(int i=0;i<B1.size();++i){
        bID1 = B1[i];

        dbb = QueryH2H(bID1,ID2);
        tempdis = m1[bID1] + dbb;
        if(tempdis < d){
            d = tempdis;
            b1f = bID1;
            d1f = m1[bID1];
            dbbf = dbb;
        }

    }
    return d;
}

int Gtree::QueryPartiParti(int ID1,int ID2){
    int d=INF;
    assert(BoundaryTag[ID1]==false && BoundaryTag[ID2]==false);
    int PID1=Nodes[ID1].inleaf, PID2=Nodes[ID2].inleaf;
    int PosID1=Nodes[ID1].inleafpos, PosID2=Nodes[ID2].inleafpos;
    vector<int> B1,B2;
    int bID1,bID2,d1,d2;
    int tempdis;
    B1 = GTree[PID1].borders; B2 = GTree[PID2].borders;

    map<int,int> m1,m2;
    m1.clear();
    m2.clear();
    int lnNum1 = GTree[PID1].leafnodes.size();
    int lnNum2 = GTree[PID2].leafnodes.size();
    //get distances to corresponding boundary vertex
    for(int i=0;i<B1.size();++i){
        bID1 = B1[i];
        d1 = GTree[PID1].mind[i*lnNum1 + PosID1];
        m1.insert({bID1,d1});
    }
    for(int i=0;i<B2.size();++i){
        bID2 = B2[i];
        d2 = GTree[PID2].mind[i*lnNum2 + PosID2];
        m2.insert({bID2,d2});
    }
    int b1f,b2f,d1f,d2f,dbb,dbbf;
   //get distances between boundary vertex
   for(int i=0;i<B1.size();++i){
       bID1 = B1[i];
       for(int j=0;j<B2.size();++j){
           bID2 = B2[j];
           dbb = QueryH2H(bID1,bID2);
           tempdis = m1[bID1] + dbb + m2[bID2];
           if(tempdis < d){
               d = tempdis;
               b1f = bID1; b2f = bID2;
               d1f = m1[bID1]; d2f =  m2[bID2];
               dbbf = dbb;
           }
       }
   }
    return d;
}

//stage-based query
int Gtree::Query(int ID1, int ID2) {
    int d=INF;
    if(Nodes[ID1].inleaf == Nodes[ID2].inleaf){//in the same partition
        if(BoundaryTag[ID1] && BoundaryTag[ID2]){//if ID1 and ID2 are boundary vertices
            d = QueryH2H(ID1,ID2);
        }else if(BoundaryTag[ID1] && !BoundaryTag[ID2] ){//if only ID1 is boundary vertex
            d = QueryPartiBoundary(ID2,ID1);
        }else if(BoundaryTag[ID2] && !BoundaryTag[ID1]){
            d = QueryPartiBoundary(ID1,ID2);
        }else if(!BoundaryTag[ID1] && !BoundaryTag[ID2]){//if both are not boundary vertex

        }
        d = Dijkstra(ID1,ID2,Nodes);
    }else{//if in different partitions
        d = QueryPartiParti(ID1,ID2);
    }

    if(BoundaryTag[ID1] && BoundaryTag[ID2]){//if ID1 and ID2 are boundary vertices
        d = QueryH2H(ID1,ID2);
    }else if(BoundaryTag[ID1] && !BoundaryTag[ID2] ){//if only ID1 is boundary vertex
        d = QueryPartiBoundary(ID2,ID1);
    }else if(BoundaryTag[ID2] && !BoundaryTag[ID1]){
        d = QueryPartiBoundary(ID1,ID2);
    }else if(!BoundaryTag[ID1] && !BoundaryTag[ID2]){//if both are not boundary vertex

    }
    return d;
}

//stage-based query for No-Boundary
int Gtree::Query_No(int ID1, int ID2) {
    int d=INF;
    if(BoundaryTag[ID1] && BoundaryTag[ID2]){//if ID1 and ID2 are boundary vertices
        d = QueryH2H(ID1,ID2);
    }else if(BoundaryTag[ID1] && !BoundaryTag[ID2] ){//if only ID1 is boundary vertex
        d = QueryPartiBoundary(ID2,ID1);
    }else if(BoundaryTag[ID2] && !BoundaryTag[ID1]){
        d = QueryPartiBoundary(ID1,ID2);
    }else if(!BoundaryTag[ID1] && !BoundaryTag[ID2]){//if both are not boundary vertex
        if(Nodes[ID1].inleaf == Nodes[ID2].inleaf){//in the same partition
            d = Dijkstra(ID1,ID2,Nodes);
        }else{//if in different partitions
            d = QueryPartiParti(ID1,ID2);
        }
    }
    return d;
}

//stage-based query
int Gtree::QueryDebug(int ID1, int ID2) {
    int d=INF;
    if(BoundaryTag[ID1] && BoundaryTag[ID2]){//if ID1 and ID2 are boundary vertices
        cout<<"Boundary-Boundary"<<endl;
        d = QueryH2H(ID1,ID2);
    }else if(BoundaryTag[ID1] && !BoundaryTag[ID2] ){//if only ID1 is boundary vertex
        cout<<"Boundary-Parti"<<endl;
        d = QueryPartiBoundary(ID2,ID1);
    }else if(BoundaryTag[ID2] && !BoundaryTag[ID1]){
        cout<<"Parti-Boundary"<<endl;
        d = QueryPartiBoundary(ID1,ID2);
    }else if(!BoundaryTag[ID1] && !BoundaryTag[ID2]){//if both are not boundary vertex
        if(Nodes[ID1].inleaf == Nodes[ID2].inleaf){//in the same partition
            cout<<"Same-Parti"<<endl;
            d = Distance_query(ID1,ID2);
        }else{//if in different partitions
            cout<<"Parti-Parti"<<endl;
            d = QueryPartiParti(ID1,ID2);
        }
    }
    return d;
}
