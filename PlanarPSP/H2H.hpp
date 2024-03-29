/*
 * H2H.cpp
 *
 *  Created on: 22 Dec 2022
 *      Author: Xinjie Zhou
 */
#include "headSP.h"


vector<int> NodeOrder_;//nodeID order
vector<int> _DD_;//true degree, temporal degree ,_DD2_
vector<int> NodeOrders;

//// For CH index construction
void Graph::CHIndexConstruct(){
    string orderfile=orderPath;
//    orderfile="/home/data/xzhouby/datasets/"+dataset+"/"+dataset+"_NC_16/vertex_orderMDE2";
//    orderfile=sourcePath+dataset+"_NC_32/vertex_orderMDE2";
    double runT1=0, runT2=0, runT3=0;
    Timer tt;

    tt.start();
    MDEContraction(orderfile);
    tt.stop();
    runT1=tt.GetRuntime();
    cout<<"Time for MDE contraction: "<<runT1<<" s."<<endl;

    IndexsizeCHWP();
}

//function for computing the index size
void Graph::IndexsizeCHWP(){
    unsigned long long m=0,m1=0,m2=0,m3=0,m4=0;
    //Overlay index
    for(int i=0;i<NeighborCon.size();i++){
        m1+=NeighborCon[i].size()*2*sizeof(int);//dis
    }

    for(int i=0;i< SCconNodesMT.size();i++){
        for(auto it=SCconNodesMT[i].begin(); it!=SCconNodesMT[i].end(); it++){
            m4+=sizeof(int)+(*it).second.size()*2*sizeof(int);
        }
    }

    //cout<<"Index size "<<(double)m1/1024/1024<<", Pruning point size "<<(double)m2/1024/1024<<endl;
    m=m1+m2+m3+m4;
    cout<<"CH label size: "<<(double)m1/1024/1024<<" MB"<<endl;
    cout<<"CH Update information size: "<<(double)m4/1024/1024<<" MB"<<endl;
    cout<<"Overall index size "<<(double)m/1024/1024<<" MB"<<endl;
}

void Graph::H2HIndexConstruct() {
    string orderfile=orderPath;
//    orderfile="/home/data/xzhouby/datasets/"+dataset+"/"+dataset+"_NC_16/vertex_orderMDE2";
//    orderfile=sourcePath+dataset+"_NC_32/vertex_orderMDE2";
    double runT1=0, runT2=0, runT3=0;
    Timer tt;

    tt.start();
    MDEContraction(orderfile);
    tt.stop();
    runT1=tt.GetRuntime();
    cout<<"Time for MDE contraction: "<<runT1<<" s."<<endl;

    tt.start();
	makeTree();
    tt.stop();
    runT2=tt.GetRuntime();
    cout<<"Time for Tree construction: "<<runT2<<" s."<<endl;

    tt.start();
	makeIndex();
    tt.stop();
    runT3=tt.GetRuntime();
    cout<<"Time for Index building: "<<runT3<<" s."<<endl;

    cout<<"Overall index construction time: "<<runT1+runT2+runT3<<" s."<<endl;

    IndexsizeH2H();
}
//function for MDE contraction
void Graph::MDEContraction(string orderfile){
    cout<<"MDE contraction..."<<endl;
    vSm.reserve(node_num);
    for(int i = 0; i < node_num; i++)
    {
        Semaphore* s = new Semaphore(1);
        vSm.push_back(s);
    }
    ifstream IF(orderfile);
//    if(true){
    if(!IF.is_open()){/// if no order file, use MDE to generate order
        cout<<"Cannot open vertex ordering file "<<orderfile<<endl;
        int Twidth=0;//tree width
        //initialize SCconNodesMT
        map<int, vector<int>> mi;
        SCconNodesMT.assign(node_num, mi);

        //initialize E
        map<int,pair<int,int>> m;
        E.assign(node_num,m);
        for(int i=0;i<Neighbor.size();i++){
            for(int j=0;j<Neighbor[i].size();j++)
                E[i].insert(make_pair(Neighbor[i][j].first,make_pair(0,1)));
        }

        _DD_.assign(node_num,0);
        DD.assign(node_num,0);

        set<DegComp> Deg;
        int degree;
        for(int i=0;i<node_num;i++){
            degree=Neighbor[i].size();
            if(degree!=0){
                _DD_[i]=degree;
                DD[i]=degree;
                Deg.insert(DegComp(i));
            }
        }

        vector<bool> exist; exist.assign(node_num,true);
        vector<bool> change; change.assign(node_num,false);

        vector<pair<int,pair<int,int>>> vect;
        NeighborCon.assign(node_num,vect); //NeighborCon.clear();
        //SCconNodes.clear();

        //cout<<"Begin to contract"<<endl;
        int count=0;

        while(!Deg.empty()){
            if(count%10000==0)
                cout<<"count "<<count<<" , treewidth "<<Twidth<<endl;
            count+=1;
            int x=(*Deg.begin()).x;

            while(true){
                if(change[x]){
                    Deg.erase(DegComp(x));
                    _DD_[x]=DD[x];
                    Deg.insert(DegComp(x));
                    change[x]=false;
                    x=(*Deg.begin()).x;
                }else
                    break;
            }

            vNodeOrder.push_back(x);
            Deg.erase(Deg.begin());
            exist[x]=false;

            vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();

            for(auto it=E[x].begin();it!=E[x].end();it++){
                if(exist[(*it).first]){
                    Neigh.push_back(*it);
                }
            }

            if(Neigh.size()>Twidth)
                Twidth=Neigh.size();

            NeighborCon[x].assign(Neigh.begin(),Neigh.end());

            //multi threads for n^2 combination
            for(int i=0;i<Neigh.size();i++){
                int y=Neigh[i].first;
                deleteEOrderGenerate(x,y);
                change[y]=true;
            }

            int stepf=Neigh.size()/threadnum;
            boost::thread_group threadf;
            for(int i=0;i<threadnum;i++){
                pair<int,int> p;
                p.first=i*stepf;
                if(i==threadnum-1)
                    p.second=Neigh.size();
                else
                    p.second=(i+1)*stepf;
                threadf.add_thread(new boost::thread(&Graph::NeighborComOrderGenerate, this, boost::ref(Neigh), p));
            }
            threadf.join_all();
        }

        NodeOrder.assign(node_num,-1);
        for(int k=0;k<vNodeOrder.size();k++){
            NodeOrder[vNodeOrder[k]]=k;
        }
        ofstream ofile(orderfile);
        ofile << node_num << endl;
        for(int i = 0; i < NodeOrder.size(); i++)
            ofile << i << " " << NodeOrder[i] << endl;
        ofile.close();
        cout<<"Finish Contract"<<" , treewidth "<<Twidth<<endl;
        exit(0);
    }
    else{///if there is an order file
        cout<<"Reading vertex ordering..."<<endl;
        cout<<"Order file: "<<orderfile<<endl;
        NodeOrder.assign(node_num, -1);
        vNodeOrder.assign(node_num, -1);
        int num, nodeID, nodeorder;
        string line;
        getline(IF,line);
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
        num=stoi(vs[0]);
        for(int i=0;i<num;i++){
            vs.clear();
            getline(IF,line);
            boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
            if(vs.size()!=2){
                cout<<"Wrong syntax for ordering. "<<line<<endl; exit(1);
            }
            nodeID=stoi(vs[0]);nodeorder=stoi(vs[1]);

            NodeOrder[nodeID]=nodeorder;
            if(nodeorder!=-1){
                vNodeOrder[nodeorder]=nodeID;
            }else{
                cout<<"Wrong order! "<<nodeID<<" "<<nodeorder<<endl; exit(1);
            }
        }
        IF.close();
        unordered_set<int> vertices; vertices.clear();
        for(int i=0;i<node_num;++i){
            if(vertices.find(vNodeOrder[i])==vertices.end()){//if not found
                vertices.insert(vNodeOrder[i]);
            }
        }
        if(vertices.size()!=node_num){
            cout<<"Order wrong! "<<vertices.size()<<" "<<node_num<<endl; exit(1);
        }
        vertices.clear();

        for(int i=0;i<2;++i){
            int id=vNodeOrder[node_num-1-i];
            cout<<"Order "<<node_num-1-i<<": "<<id<<" "<<Neighbor[id].size()<<endl;
        }
        for(int i=0;i<2;++i){
            int id=vNodeOrder[i];
            cout<<"Order "<<i<<": "<<id<<" "<<Neighbor[id].size()<<endl;
        }

        vector<pair<int,pair<int,int>>> vect;
        NeighborCon.assign(node_num,vect);

        map<int, vector<int>> mi;
        SCconNodesMT.assign(node_num, mi);//record the supportive vertices of a shortcut, only record edge once by leveraging the ID positions of endpoints

        //initialize E
        map<int,pair<int,int>> m;
        E.assign(node_num,m);
        for(int i=0;i<Neighbor.size();i++){
            for(int j=0;j<Neighbor[i].size();j++)
                E[i].insert(make_pair(Neighbor[i][j].first,make_pair(Neighbor[i][j].second,1)));
        }

        vector<bool> exist; exist.assign(node_num,true);
        //vector<bool> change; change.assign(nodenum,false);

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

                for(int i=0;i<Neigh.size();i++){
                    int y=Neigh[i].first;
                    deleteEorder(x,y);
                    //change[y]=true;
                }

                int ID1,ID2;
                if(Neigh.size()<=100){
//                if(true){
                    //single thread
                    for(int i=0;i<Neigh.size();i++){
                        ID1=Neigh[i].first;
                        for(int j=i+1;j<Neigh.size();j++){
                            ID2=Neigh[j].first;
                            insertEorder(ID1,ID2,Neigh[i].second.first+Neigh[j].second.first);
                            if(ID1<ID2){
                                if(SCconNodesMT[ID1].find(ID2)==SCconNodesMT[ID1].end()){//if not found
                                    SCconNodesMT[ID1].insert({ID2,vector<int>()});
                                }
                                SCconNodesMT[ID1][ID2].push_back(x);//only record onece
                            }
                            else if(ID2<ID1){
                                if(SCconNodesMT[ID2].find(ID1)==SCconNodesMT[ID2].end()){//if not found
                                    SCconNodesMT[ID2].insert({ID1,vector<int>()});
                                }
                                SCconNodesMT[ID2][ID1].push_back(x);
                            }

                        }
                    }
                }else{
                    if(Neigh.size()>threadnum){
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
                    }else{
                        boost::thread_group thread;
                        for(int i=0;i<Neigh.size();i++){
                            pair<int,int> p;
                            p.first=i; p.second=(i+1);
                            thread.add_thread(new boost::thread(&Graph::NeighborComorder, this, boost::ref(Neigh), p, x));
                        }
                        thread.join_all();
                    }

                }

            }
            else{
                cout<<"Wrong order! "<<x<<" "<<nodeorder<<endl; exit(1);
            }
        }
    }
    NodeOrder_ = NodeOrder;
}
//function for H2H tree construction
void Graph::makeTree(){
    cout<<"Building H2H tree..."<<endl;
	vector<int> vecemp; //vecemp.clear();
	VidtoTNid.assign(node_num,vecemp);//record the tree node id whose unique vertex involves this vertex as neighbor

	rank.assign(node_num,0);
	//Tree.clear();
	int len=vNodeOrder.size()-1;
	heightMax=0;

	Node rootn;
	int x=vNodeOrder[len];//from the highest vertex
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

    int treewidth=0;
	int nn;
	for(;len>=0;len--){
		int x=vNodeOrder[len];
		Node nod;
		nod.vert=NeighborCon[x];
		nod.uniqueVertex=x;
        if(treewidth<NeighborCon[x].size()+1){
            treewidth=NeighborCon[x].size()+1;
        }
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
		if(nod.height>heightMax) {
            heightMax=nod.height;
        }
        rank[x]=Tree.size();

		Tree.push_back(nod);
		//cout<<"len "<<len<<" , ID "<<x<<endl;
	}

    cout<<"Tree height: "<<heightMax<<" ; treewidth: "<<treewidth<<endl;
}
//function of H2H index construction
void Graph::makeIndex(){
    cout<<"Building H2H index..."<<endl;
	makeRMQ();

	//initialize
	vector<int> list; //list.clear();
	list.push_back(Tree[0].uniqueVertex);
	Tree[0].pos.clear();
	Tree[0].pos.push_back(0);
    Tree[0].vAncestor=list;

	for(int i=0;i<Tree[0].ch.size();i++){
		makeIndexDFS(Tree[0].ch[i],list);
	}

}
/// Functions for MDE contraction
void Graph::deleteEOrderGenerate(int u,int v){
    if(E[u].find(v)!=E[u].end()){
        E[u].erase(E[u].find(v));
        DD[u]--;
    }

    if(E[v].find(u)!=E[v].end()){
        E[v].erase(E[v].find(u));
        DD[v]--;
    }
}

void Graph::NeighborComOrderGenerate(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p){
//    sm->wait();
    int ID1, w1;
    int ID2, w2;
    for(int k=p.first;k<p.second;k++){
        ID1=Neighvec[k].first;
        for(int h=0;h<Neighvec.size();h++){
            ID2=Neighvec[h].first;
//            vSm[ID1]->wait();
//            vSm[ID2]->wait();
            insertEMTOrderGenerate(ID1, ID2, 1);//for order generation, only need unweighted graph
//            vSm[ID1]->notify();
//            vSm[ID2]->notify();
        }
    }
//    sm->notify();
}

void Graph::insertEMTOrderGenerate(int u,int v,int w){
    if(E[u].find(v)==E[u].end()){
        E[u].insert(make_pair(v,make_pair(w,1)));
        DD[u]++;
//        DD2[u]++;
    }
}

void Graph::deleteEorder(int u,int v){
    if(E[u].find(v)!=E[u].end()){
        E[u].erase(E[u].find(v));
        //DD[u]--;
    }

    if(E[v].find(u)!=E[v].end()){
        E[v].erase(E[v].find(u));
        //DD[v]--;
    }
}

void Graph::insertEorder(int u,int v,int w){
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

void Graph::NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x){
//    sm->wait();
    int ID1, w1;
    int ID2, w2;
    for(int k=p.first;k<p.second;k++){
        ID1=Neighvec[k].first;
        w1=Neighvec[k].second.first;
        for(int h=0;h<Neighvec.size();h++){
            ID2=Neighvec[h].first;
            w2=Neighvec[h].second.first;
            if(ID1==ID2){
                continue;
            }
            insertEMTorder(ID1, ID2, w1+w2);
            if(ID1<ID2){
                if(SCconNodesMT[ID1].find(ID2)==SCconNodesMT[ID1].end()){//if not found
                    SCconNodesMT[ID1].insert({ID2,vector<int>()});
                }
                SCconNodesMT[ID1][ID2].push_back(x);
            }

        }
    }
//    sm->notify();
}
/// Functions for Tree contraction
int Graph::match(int x,vector<pair<int,pair<int,int>>> &vert){
    int nearest=vert[0].first;
    for(int i=1;i<vert.size();i++){
        if(rank[vert[i].first]>rank[nearest])
            nearest=vert[i].first;
    }
    int p=rank[nearest];
    return p;
}

/// Functions for Tree index contraction
void Graph::makeRMQ(){
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
            if(Tree[p].vert[i].first==list[j]){//get the original distance information by shortcuts
                Tree[p].pos[i]=j;
                Tree[p].dis[j]=Tree[p].vert[i].second.first;
                Tree[p].cnt[j]=1;
                break;
            }
        }
    }
    Tree[p].pos[NeiNum]=list.size();
//    Tree[p].dis.push_back(0);//distance to itself
    Tree[p].vAncestor=list;
    Tree[p].vAncestor.push_back(Tree[p].uniqueVertex);//the last vertex is the tree node

    //dis
    for(int i=0;i<NeiNum;i++){
        int x=Tree[p].vert[i].first;
        int disvb=Tree[p].vert[i].second.first;
        int k=Tree[p].pos[i];//the k-th ancestor is x, the i-th neighbor is also the k-th ancesotr

        for(int j=0;j<list.size();j++){
            int y=list[j];//the j-th ancestor is y

            int z;//the distance from x to y
            if(k!=j){
                if(k<j)//if x is the ancestor of y
                    z=Tree[rank[y]].dis[k];
                else if(k>j)//if y is the ancestor of x
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

//function for computing the index size
void Graph::IndexsizeH2H(){
    unsigned long long m=0,m1=0,m2=0,m3=0,m4=0;
    //Overlay index
    for(int i=0;i<Tree.size();i++){
        m1+=Tree[i].dis.size()*2*sizeof(int);//dis
        m3+=Tree[i].pos.size()*sizeof(int);//pos
        m2+=Tree[i].cnt.size()*sizeof(int);//cnt
        m2+=Tree[i].vert.size()*3*sizeof(int);//neighID/weight/count

    }

    for(int i=0;i< SCconNodesMT.size();i++){
        for(auto it=SCconNodesMT[i].begin(); it!=SCconNodesMT[i].end(); it++){
            m4+=sizeof(int)+(*it).second.size()*2*sizeof(int);
        }
    }

    //cout<<"Index size "<<(double)m1/1024/1024<<", Pruning point size "<<(double)m2/1024/1024<<endl;
    m=m1+m2+m3+m4;
    cout<<"Distance label size: "<<(double)m1/1024/1024<<" MB"<<endl;
    cout<<"H2H label size: "<<(double)(m1+m3)/1024/1024<<" MB"<<endl;
    cout<<"CH Update information size: "<<(double)m4/1024/1024<<" MB"<<endl;
    cout<<"H2H Update information size: "<<(double)(m2+m4)/1024/1024<<" MB"<<endl;
    cout<<"Overall index size "<<(double)m/1024/1024<<" MB"<<endl;
}


int	Graph::QueryCHWP(int ID1, int ID2){
    if(ID1==ID2) return 0;
    if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
    int d=INF;
    benchmark::heap<2,int,int> fHeapForward(node_num);
    benchmark::heap<2, int, int> fHeapBackward(node_num);

    //closed or not
    vector<bool> vVisitedF(node_num, false);
    vector<bool> vVisitedB(node_num, false);
    //the existing shortest distance
    vector<int>	vDistanceForward(node_num, INF);
    vector<int>	vDistanceBackward(node_num, INF);
    //stop search or not
    bool bF = false;
    bool bB = false;
    vDistanceForward[ID1] = 0;
    vDistanceBackward[ID2] = 0;
    fHeapForward.update(ID1,0);
    fHeapBackward.update(ID2,0);

    int topNodeIDForward, topNodeIDBackward,topDisForward,topDisBackward, neighborNodeID, neighborLength;

    while(!fHeapForward.empty() || !fHeapBackward.empty() )
    {
        if(bF && bB)
            break;
        if(bF && fHeapBackward.empty())
            break;
        if(bB && fHeapForward.empty())
            break;
        //Forward Search
        if(!fHeapForward.empty() && !bF)
        {
            fHeapForward.extract_min(topNodeIDForward, topDisForward);
            //cout<<topNodeIDForward<<" "<<topDisForward<<" "<<(NodeOrder[topNodeIDForward]>NodeOrder[137099])<<endl;

            if(vDistanceForward[topNodeIDForward] > d)
                bF = true;

            vVisitedF[topNodeIDForward] = true;

            if(vVisitedB[topNodeIDForward]){
                int distTmp=topDisForward+vDistanceBackward[topNodeIDForward];
                if(distTmp<d){
                    d=distTmp;
                    //cout<<"forwardtopvalue "<<topDisForward<<" "<<vDistanceBackward[topNodeIDForward]<<" "<<d<<"meet "<<topNodeIDForward<<endl;
                }
            }

            for(auto out=NeighborCon[topNodeIDForward].begin();out!=NeighborCon[topNodeIDForward].end();out++){
                neighborNodeID = (*out).first;
                neighborLength = (*out).second.first;

                int df = vDistanceForward[topNodeIDForward] + neighborLength;
                if(!vVisitedF[neighborNodeID]){
                    if(vDistanceForward[neighborNodeID] > df){
                        //if(neighborNodeID==37783) cout<<"neighdis "<<vDistanceForward[neighborNodeID]<<" adddis "<<df<<endl;
                        vDistanceForward[neighborNodeID] = df;
                        fHeapForward.update(neighborNodeID, df);
                    }
                }
            }
        }

        //Backward Search
        if(!fHeapBackward.empty() && !bB)
        {
            fHeapBackward.extract_min(topNodeIDBackward, topDisBackward);

            if(vDistanceBackward[topNodeIDBackward] > d)
                bB = true;

            vVisitedB[topNodeIDBackward] = true;

            if(vVisitedF[topNodeIDBackward]){
                int distTmp=topDisBackward+vDistanceForward[topNodeIDBackward];
                if(distTmp<d){
                    d=distTmp;
                    //cout<<"backtopvalue "<<topDisBackward<<" "<<vDistanceForward[topNodeIDBackward]<<" "<<d<<"meet "<<topNodeIDBackward<<endl;
                }
            }

            for(auto in=NeighborCon[topNodeIDBackward].begin();in!=NeighborCon[topNodeIDBackward].end();in++){
                neighborNodeID = (*in).first;
                neighborLength = (*in).second.first;

                int db = vDistanceBackward[topNodeIDBackward] + neighborLength;
                if(!vVisitedB[neighborNodeID]){
                    if(vDistanceBackward[neighborNodeID]>db){
                        vDistanceBackward[neighborNodeID] = db;
                        fHeapBackward.update(neighborNodeID, db);
                    }
                }
            }
        }
    }
    return d;
}

int Graph::QueryH2H(int ID1,int ID2){
    if(ID1==ID2) return 0;
    if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
    int r1=rank[ID1], r2=rank[ID2];
    int LCA=LCAQuery(r1,r2);
//    cout<<r1<<" "<<r2<<" "<<LCA<<" "<<Tree.size()<<endl;
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

//PLL


///////////////////////////// Index Maintenance ////////////////////////////////

/// CHWP algorithm

void Graph::CHdecBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch){
    //maintain the index caused by the weight change
    //NodeOrders.clear();
    NodeOrders.assign(NodeOrder.begin(),NodeOrder.end());
    set<OrderCompCH> OC;
    map<pair<int,int>,int> OCdis;//{(s,t),d} maintain the fresh distance and avoid search in the adjacent list
    //OC.clear(); OCdis.clear();

    int a,b,newW;//the weight of (a,b) decrease to newW
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first;
        b=wBatch[k].first.second;
        newW=wBatch[k].second.second;

        //modify the information in original graph
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

        if(NodeOrder[a]<NodeOrder[b]){
            for(int i=0;i<NeighborCon[a].size();i++){
                if(NeighborCon[a][i].first==b){
                    if(NeighborCon[a][i].second.first>newW){
                        //cout<<OutNeighborCon[a][i].second.first<<"..........."<<newW<<endl;
                        NeighborCon[a][i].second.first=newW;
                        NeighborCon[a][i].second.second=1;

                        OCdis[make_pair(a,b)]=newW;
                        OC.insert(OrderCompCH(a,b));
                    }else if(NeighborCon[a][i].second.first==newW)
                        NeighborCon[a][i].second.second+=1;
                    break;
                }
            }
        }else{
            for(int i=0;i<NeighborCon[b].size();i++){
                if(NeighborCon[b][i].first==a){
                    if(NeighborCon[b][i].second.first>newW){
                        NeighborCon[b][i].second.first=newW;
                        NeighborCon[b][i].second.second=1;

                        OCdis[make_pair(b,a)]=newW;
                        OC.insert(OrderCompCH(b,a));
                    }else if(NeighborCon[b][i].second.first==newW)
                        NeighborCon[b][i].second.second+=1;
                    break;
                }
            }
        }
    }


    while(!OC.empty()){
        int s=(*OC.begin()).x; int t=(*OC.begin()).y;
        int wt;
        OC.erase(OC.begin());
        wt=OCdis[make_pair(s,t)];
        map<int,int> InM2t; //InM2t.clear();
        vector<pair<int,int>> InMLower; //InMLower.clear();
        for(int i=0;i<NeighborCon[s].size();i++){
            if(NodeOrder[NeighborCon[s][i].first]>NodeOrder[t])
                InM2t.insert(make_pair(NeighborCon[s][i].first,NeighborCon[s][i].second.first));
            else if(NodeOrder[NeighborCon[s][i].first]<NodeOrder[t])
                InMLower.push_back(make_pair(NeighborCon[s][i].first,NeighborCon[s][i].second.first));
        }
        int inID,inW,inWt;
        for(int i=0;i<NeighborCon[t].size();i++){
            inID=NeighborCon[t][i].first;
            if(InM2t.find(inID)!=InM2t.end()){
                inW=InM2t[inID];
                inWt=NeighborCon[t][i].second.first;
                if(inWt>inW+wt){
                    NeighborCon[t][i].second.first=inW+wt;
                    NeighborCon[t][i].second.second=1;
                    OCdis[make_pair(t,inID)]=inW+wt;
                    OrderCompCH oc={t,inID};
                    OC.insert(oc);
                }else if(inWt==inW+wt){
                    NeighborCon[t][i].second.second+=1;
                }
            }
        }

        for(int i=0;i<InMLower.size();i++){
            inID=InMLower[i].first; inW=InMLower[i].second;
            for(int j=0;j<NeighborCon[inID].size();j++){
                if(NeighborCon[inID][j].first==t){
                    inWt=NeighborCon[inID][j].second.first;
                    if(inWt>inW+wt){
                        NeighborCon[inID][j].second.first=inW+wt;
                        NeighborCon[inID][j].second.second=1;

                        OCdis[make_pair(inID,t)]=inW+wt;
                        OrderCompCH oc={inID,t};
                        OC.insert(oc);
                    }else if(inWt==inW+wt)
                        NeighborCon[inID][j].second.second+=1;
                    break;
                }
            }
        }
    }//finish change index
}

void Graph::CHincBatMT(vector<pair<pair<int,int>,pair<int,int>>>& wBatch){
    //NodeOrders.clear();
    NodeOrders.assign(NodeOrder.begin(),NodeOrder.end());
    set<OrderCompCH> OC; //OC.clear();
    map<pair<int,int>,int> OCdis;//{(s,t),d} maintain the old distance before refreshed and avoid search in the adjacent list
    //OCdis.clear();

    for(int wb=0;wb<wBatch.size();wb++){
        int a=wBatch[wb].first.first;
        int b=wBatch[wb].first.second;
        int oldW=wBatch[wb].second.first;
        int newW=wBatch[wb].second.second;

        //modify the original graph information
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


        if(NodeOrder[a]<NodeOrder[b]){
            for(int i=0;i<NeighborCon[a].size();i++){
                if(NeighborCon[a][i].first==b){
                    if(NeighborCon[a][i].second.first==oldW){
                        NeighborCon[a][i].second.second-=1;
                        if(NeighborCon[a][i].second.second<1){
                            OrderCompCH oc={a,b};
                            OC.insert(oc);
                            OCdis[make_pair(a,b)]=oldW;
                        }
                    }
                    break;
                }
            }
        }else{
            for(int i=0;i<NeighborCon[b].size();i++){
                if(NeighborCon[b][i].first==a){
                    if(NeighborCon[b][i].second.first==oldW){
                        NeighborCon[b][i].second.second-=1;
                        if(NeighborCon[b][i].second.second<1){
                            OrderCompCH oc={b,a};
                            OC.insert(oc);
                            OCdis[make_pair(b,a)]=oldW;
                        }
                    }
                    break;
                }
            }
        }
    }

    while(!OC.empty()){
        int s=(*OC.begin()).x; int t=(*OC.begin()).y;
        int wt;
        OC.erase(OC.begin());
        wt=OCdis[make_pair(s,t)];//distance of s--->t before change
        int inID,inW;
        map<int,int> HigherIn; vector<pair<int,int>> LowerIn;
        //HigherIn.clear(); LowerIn.clear();
        //the shortcuts infected by s-->t
        for(int i=0;i<NeighborCon[s].size();i++){
            inID=NeighborCon[s][i].first;
            inW=NeighborCon[s][i].second.first;
            if(NodeOrder[inID]<NodeOrder[t]){
                LowerIn.emplace_back(inID,inW);
            }else if(NodeOrder[inID]>NodeOrder[t]){
                HigherIn.insert(make_pair(inID,inW));
            }
        }
        for(int i=0;i<NeighborCon[t].size();i++){
            inID=NeighborCon[t][i].first;
            if(HigherIn.find(inID)!=HigherIn.end()){
                inW=HigherIn[inID];
                if(NeighborCon[t][i].second.first==wt+inW){
                    NeighborCon[t][i].second.second-=1;
                    if(NeighborCon[t][i].second.second<1){
                        OrderCompCH oc={t,inID};
                        OC.insert(oc);
                        OCdis[make_pair(t,inID)]=wt+inW;
                    }
                }
            }
        }
        for(int i=0;i<LowerIn.size();i++){
            inID=LowerIn[i].first; inW=LowerIn[i].second;
            for(int j=0;j<NeighborCon[inID].size();j++){
                if(NeighborCon[inID][j].first==t){
                    if(NeighborCon[inID][j].second.first==inW+wt){
                        NeighborCon[inID][j].second.second-=1;
                        if(NeighborCon[inID][j].second.second<1){
                            OrderCompCH oc={inID,t};
                            OC.insert(oc);
                            OCdis[make_pair(inID,t)]=wt+inW;
                        }
                    }
                    break;
                }
            }
        }

        //get the new weight value of s-->t
        wt=INF; int countwt=0;
        for(int i=0;i<Neighbor[s].size();i++){
            if(Neighbor[s][i].first==t){
                wt=Neighbor[s][i].second;//the weight value in the original graph
                countwt=1;
                break;
            }
        }
        int ssw,wtt,wid;
        vector<int> Wnodes; //Wnodes.clear();
        if(s<t){
            //Wnodes=SCconNodes[make_pair(s,t)]; //cout<<"wid num "<<Wnodes.size()<<endl;
            Wnodes=SCconNodesMT[s][t];
        }else{
            //Wnodes=SCconNodes[make_pair(t,s)];
            Wnodes=SCconNodesMT[t][s];
        }

        for(int i=0;i<Wnodes.size();i++){
            wid=Wnodes[i];
            for(int j=0;j<NeighborCon[wid].size();j++){
                if(NeighborCon[wid][j].first==s){
                    ssw=NeighborCon[wid][j].second.first;
                }
                if(NeighborCon[wid][j].first==t){
                    wtt=NeighborCon[wid][j].second.first;
                }
            }

            if(ssw+wtt<wt){
                wt=ssw+wtt;
                countwt=1;
            }else if(ssw+wtt==wt){
                countwt+=1;
            }
        }

        //refresh the weight value of s--t in the index
        for(int i=0;i<NeighborCon[s].size();i++){
            if(NeighborCon[s][i].first==t){
//                cout<<"Refresh shortcut: "<<s<<" "<<t<<" "<<NeighborCon[s][i].second.first<<" "<<wt<<endl;
                NeighborCon[s][i].second.first=wt;
                NeighborCon[s][i].second.second=countwt;
                break;
            }
        }
    }
}

void Graph::H2HdecBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch){
	map<int,int> checkedDis;

	for(int i=0;i<Tree.size();i++){
		Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
	}

	//NodeOrderss.clear();
//	NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
	vector<set<int>> SCre; //SCre.clear();
	set<int> ss; //ss.clear();
	SCre.assign(node_num,ss);//{vertexID, set<int>}
	set<OrderCompp> OC; //OC.clear();//vertexID in decreasing node order

	set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed

	int a,b,oldW,newW,lid,hid;
	for(int k=0;k<wBatch.size();k++){
		a=wBatch[k].first.first; b=wBatch[k].first.second; oldW=wBatch[k].second.first;newW=wBatch[k].second.second;
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
				if(LCAQuery(rnew,r)!=rnew){
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
		EachNodeProBDis5(rank[ProBeginVertexID], linee, vertexIDChL,checkedDis);
	}
	//return checkedDis.size();
}

void Graph::EachNodeProBDis5(int child,vector<int>& line,set<int>& vertexIDChL, map<int,int>& checkedDis){
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

void Graph::H2HincBatMT(vector<pair<pair<int,int>,pair<int,int>>>& wBatch){
	int checknum=0;
	map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]), original weight of the affected shortcut, maintain the old distance before refreshed and avoid search in the adjacent list
	//OCdis.clear();

	//NodeOrderss.clear();
//	NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
	vector<set<int>> SCre; //SCre.clear(); the affected shortcut pair
	set<int> ss; ss.clear();
	SCre.assign(node_num,ss);//{vertexID, set<int>}
	set<OrderCompp> OC; OC.clear();//the lower-order vertex of the affected shortcut, vertexID in decreasing node order

	for(int k=0;k<wBatch.size();k++){
		int a=wBatch[k].first.first;
		int b=wBatch[k].first.second;
		int oldW=wBatch[k].second.first;
		int newW=wBatch[k].second.second;

		if(oldW<newW){
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
                        if(Tree[rank[lid]].vert[i].second.second<1){//the shortcut needs update, should be increased
                            OCdis[make_pair(lid,hid)]=oldW;//original weight of the affected shortcut
                            SCre[lid].insert(hid);//the affected shortcut pair
                            OC.insert(OrderCompp(lid));//the lower-order vertex of the affected shortcut
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
    /// Shortcut update
	while(!OC.empty()){
		ProID=(*OC.begin()).x;//from the lowest-order vertex
		OC.erase(OC.begin());
		vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
		influence=false;

		// get the ancestors of ProID, each ProID corresponds to a line
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
					Lnei.emplace_back(Vert[j].first,Vert[j].second.first);
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
			if(Tree[rank[ProID]].FN[cidH]){//if the distance label is from shortcut, then the label may be affected.
				influence=true;
				//higher than Cid
				for(int i=0;i<cidH;i++){
					if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[Cid]].dis[i]){
						Tree[rank[ProID]].cnt[i]-=1;
					}
				}

				//equal to Cid
				Tree[rank[ProID]].FN[cidH]=false;//? may still be the source
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

			if(ProID<Cid)
				Wnodes=SCconNodesMT[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
			else
				Wnodes=SCconNodesMT[Cid][ProID];
			if(!Wnodes.empty()){
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
				if(LCAQuery(rnew,r)!=rnew){//if they are in different branches
					ProBeginVertexSetNew.push_back(ProBeginVertexSet[i]);
				}
			}
			ProBeginVertexSet=ProBeginVertexSetNew;//identify the roots
		}

	}

	int ProBeginVertexID;
//    cout<<"Root number: "<<ProBeginVertexSet.size()<<endl;
	for(int i=0;i<ProBeginVertexSet.size();i++){//for each root
		ProBeginVertexID=ProBeginVertexSet[i];
//        cout<<i<<" "<<ProBeginVertexID<<" "<<Tree[rank[ProBeginVertexID]].height<<endl;
		vector<int> linee; //linee.clear();
		linee.reserve(heightMax);
		int pachidd=Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
		while(Tree[rank[pachidd]].height>1){
			linee.insert(linee.begin(),pachidd);
			pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
		}
		linee.insert(linee.begin(),pachidd);

		eachNodeProcessIncrease1(rank[ProBeginVertexID], linee,checknum);
	}
	//return checknum;
}

void Graph::eachNodeProcessIncrease1(int children, vector<int>& line, int& changelabel){
	int childID=Tree[children].uniqueVertex;
	int childH=Tree[children].height-1;
	for(int i=0;i<Tree[children].dis.size();i++){
		if(Tree[children].cnt[i]==0){//if the distance label to i-th ancestor should be maintained
//        if(true){
			changelabel+=1;
			//firstly, check which dis can be infected
			int disBF=Tree[children].dis[i];
			int PID;
			//chidlID
			for(int k=0;k<VidtoTNid[childID].size();k++){//for the tree node that contains childID as vert element
				PID=VidtoTNid[childID][k];
				if(Tree[PID].FN[childH] && Tree[PID].dis[i]==disBF+Tree[PID].dis[childH]){//if label is from shortcut
					Tree[PID].cnt[i]-=1;
				}
			}

			//line[i]
			for(int k=0;k<VidtoTNid[line[i]].size();k++){
				PID=VidtoTNid[line[i]][k];
//				if(Tree[PID].height>Tree[children].height){///modified for correctness
                if(Tree[PID].height>Tree[children].height && Tree[PID].vAncestor[childH] == childID){
                    if(PID>Tree.size()){
                        cout<<"PID error! "<<PID<<" "<<Tree.size()<<endl; exit(1);
                    }
                    if(childH>Tree[PID].dis.size()){
                        cout<<"childH error! "<<childH<<" "<<Tree[PID].dis.size()<<": "<<children<<"("<<Tree[children].height<<") "<<PID<<"("<<Tree[PID].height<<")"<<endl; exit(1);
                    }
					if(Tree[PID].FN[i] && Tree[PID].dis[childH]==disBF+Tree[PID].dis[i]){///
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
			if(DDvb==dis) {
                Tree[children].FN[i]=true;
            }
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
