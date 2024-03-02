/*
 * Construction.cpp
 *
 *  Created on: 24 August 2023
 *      Author: Xinjie ZHOU
 */
#include "headPSP.h"
#include "PH2H.hpp"
#include "PPLL.hpp"

/// Index Construction
void Graph::IndexConstruction(){
    if(algoChoice==CH){
        cout<<"Partitioned CH."<<endl;
        PCHIndexConstruct(PSPStrategy);
    }else if(algoChoice==H2H){
        cout<<"Partitioned H2H."<<endl;
        PH2HIndexConstruct(PSPStrategy);
    }else if(algoChoice==PLL){
        cout<<"Partitioned PLL."<<endl;
        PPLLIndexConstruct(PSPStrategy);
    }else if(algoChoice==Dijk){
        cout<<"Index-free search."<<endl;
        ReadGraph(graphfile);//
    }
}
//function of hybrid multi-stage SP index construction
void Graph::HybridSPIndexConstruct(){
    string orderfile=graphfile+".order";
    orderfile=graphfile+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_orderMDE2";

    double runT1=0, runT2=0, runT3=0;
    Timer tt;

    ReadGraph(graphfile);//

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
        SCconNodesMT.assign(node_num, map<int, vector<pair<int,int>>>());

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

        SCconNodesMT.assign(node_num, map<int, vector<pair<int,int>>>());//record the supportive vertices of a shortcut, only record edge once by leveraging the ID positions of endpoints

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
                            int temp=Neigh[i].second.first+Neigh[j].second.first;
                            insertEorder(ID1,ID2,temp);
                            if(ID1<ID2){
                                if(SCconNodesMT[ID1].find(ID2)==SCconNodesMT[ID1].end()){//if not found
                                    SCconNodesMT[ID1].insert({ID2,vector<pair<int,int>>()});
                                }
                                SCconNodesMT[ID1][ID2].emplace_back(x,temp);//only record onece
                            }
                            else if(ID2<ID1){
                                if(SCconNodesMT[ID2].find(ID1)==SCconNodesMT[ID2].end()){//if not found
                                    SCconNodesMT[ID2].insert({ID1,vector<pair<int,int>>()});
                                }
                                SCconNodesMT[ID2][ID1].emplace_back(x,temp);
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
                            thread.add_thread(new boost::thread(&Graph::NeighborComorderH2H, this, boost::ref(Neigh), p, x));
                        }
                        thread.join_all();
                    }else{
                        boost::thread_group thread;
                        for(int i=0;i<Neigh.size();i++){
                            pair<int,int> p;
                            p.first=i; p.second=(i+1);
                            thread.add_thread(new boost::thread(&Graph::NeighborComorderH2H, this, boost::ref(Neigh), p, x));
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

    int nn;
    for(;len>=0;len--){
        int x=vNodeOrder[len];
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
        if(nod.height>heightMax) {
            heightMax=nod.height;
        }
        rank[x]=Tree.size();

        Tree.push_back(nod);
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
//function of H2H index construction
void Graph::makeIndex(){
    cout<<"Building H2H index..."<<endl;
    makeRMQ(toRMQ, RMQIndex, Tree);//build LCA index

    //initialize
    vector<int> list; //list.clear();
    list.push_back(Tree[0].uniqueVertex);
    Tree[0].pos.clear();
    Tree[0].pos.push_back(0);

    for(int i=0;i<Tree[0].ch.size();i++){
        makeIndexDFS(Tree[0].ch[i],list);
    }

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

void Graph::PCHIndexConstruct(int strategy) {
    double runT1, runT2, runT3, runT4, runT5;
    runT1=0, runT2=0, runT3=0, runT4=0, runT5=0;

    /// Read order and partitions
    string orderfile=graphfile+".orderP";
    orderfile=graphfile+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_order";
    orderfile=graphfile+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_orderMDE";
    orderfile=graphfile+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_orderMDE2";
//    orderfile=graphfile+".order";
//#ifdef __APPLE__
////    cout<<"The platform is macOS."<<endl;
//#else
//    orderfile="/home/data/xzhouby/datasets/"+dataset+"/"+dataset+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_orderMDE2";
//#endif
    ReadOrder(orderfile);

    string partitionfile=graphfile+"_"+algoParti+"_"+to_string(partiNum);
    GraphPartitionRead(partitionfile);//read partitions

//    vSm.reserve(node_num);
//    for(int i = 0; i < node_num; i++)
//    {
//        Semaphore* s = new Semaphore(1);
//        vSm.push_back(s);
//    }

    Timer tt;
    if(strategy==PreBoundary){
        /// All-pair boundary shortcut computation
        tt.start();
        AllPairBoundaryDisCompute(true);
//        AllPairBoundaryDisCompute(false);
        tt.stop();
        runT1 = tt.GetRuntime();
        cout<<"All-pair boundary shortcut construction time: "<<runT1<<" s"<<endl;

        /// Partition index and construction
        tt.start();
        Construct_PartiIndex(true, false);
        tt.stop();
        runT2 = tt.GetRuntime();
        cout<<"Partition index construction time: "<<runT2<<" s"<<endl;

        /// Overlay index construction
        tt.start();
//    Construct_core(algoCoreC);
//    WriteCoreIndex(graphfile);
        Construct_OverlayIndex(false);
        tt.stop();
        runT3 = tt.GetRuntime();
        cout<<"Overlay index construction time: "<<runT3<<" s."<< endl;

    }
    else if(strategy>=NoBoundary){
        tt.start();
        /// Partition index and construction
        Construct_PartiIndex(true, false);
//        Construct_PartiIndex(false, false);//no parallel
        tt.stop();
        runT1 = tt.GetRuntime();
        cout<<"Partition index construction time: "<<runT1<<" s"<<endl;

        /// Overlay graph construction
        tt.start();
//        Construct_OverlayGraph(true,true);//all-pair boundary shortcut
        Construct_OverlayGraphNoAllPair(true,true);//all-pair boundary shortcut
        tt.stop();
        runT2 = tt.GetRuntime();
        cout<<"Overlay graph construction time: "<<runT2<<" s."<< endl;

        /// Overlay index construction
        tt.start();
        Construct_OverlayIndex(false);
        tt.stop();
        runT3 = tt.GetRuntime();
        cout<<"Overlay index construction time: "<<runT3<<" s."<< endl;


        /// Partition index repair
        if(strategy==PostBoundary){
            tt.start();
            ConstructPartitionPost(true, true);
//            ConstructPartitionPost(false, true);
            tt.stop();
            runT4 += tt.GetRuntime();
            cout<<"Time for post partition construction: "<<tt.GetRuntime()<<" s."<<endl;

            tt.start();
            ConstructPartitionPostIndex(true, false);
            tt.stop();
            runT4 += tt.GetRuntime();
            cout<<"Time for post partition index construction: "<<tt.GetRuntime()<<" s."<<endl;

//        tt.start();
//        Repair_PartiIndex(true);
////        Repair_PartiIndex(false);
//        tt.stop();
//        runT3 = tt.GetRuntime();
//        cout<<"Partition index repair time: "<<runT4<<" s"<<endl;
        }
//    Compute_tree_label(ifParallel, ifOpt);//Construct periphery index (H2H label + interface label)
    }

    cout<<"Overall Construction Time: "<<runT1+runT2+runT3+runT4+runT5<<" s."<<endl;

    IndexSizePH2H();//index (labeling+pruning point) size computation
}

void Graph::PH2HIndexConstruct(int strategy) {
    double runT1, runT2, runT3, runT4, runT5;
    runT1=0, runT2=0, runT3=0, runT4=0, runT5=0;

    /// Read order and partitions
    string orderfile=graphfile+".orderP";
    orderfile=graphfile+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_order";
    orderfile=graphfile+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_orderMDE";
    orderfile=graphfile+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_orderMDE2";
//#ifdef __APPLE__
////    cout<<"The platform is macOS."<<endl;
//#else
//    orderfile="/home/data/xzhouby/datasets/"+dataset+"/"+dataset+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_orderMDE2";
//#endif
    ReadOrder(orderfile);

    string partitionfile=graphfile+"_"+algoParti+"_"+to_string(partiNum);
    GraphPartitionRead(partitionfile);//read partitions

//    vSm.reserve(node_num);
//    for(int i = 0; i < node_num; i++)
//    {
//        Semaphore* s = new Semaphore(1);
//        vSm.push_back(s);
//    }

    Timer tt;


    if(strategy==PreBoundary){
        /// All-pair boundary shortcut computation
        tt.start();
        AllPairBoundaryDisCompute(true);
        tt.stop();
        runT1 = tt.GetRuntime();
        cout<<"All-pair boundary shortcut construction time: "<<runT1<<" s"<<endl;

        /// Partition index and construction
        tt.start();
        Construct_PartiIndex(true, true);
        tt.stop();
        runT2 = tt.GetRuntime();
        cout<<"Partition index construction time: "<<runT2<<" s"<<endl;

        /// Overlay index construction
        tt.start();
        Construct_OverlayIndex(true);
        tt.stop();
        runT3 = tt.GetRuntime();
        cout<<"Overlay index construction time: "<<runT3<<" s."<< endl;
    }else if(strategy>=NoBoundary){
        tt.start();
        /// Partition index and Overlay graph construction
//    Construct_PartiIndex(false);
        Construct_PartiIndex(true, true);
//    PMHLConstructPartiIndexCH(true);
        tt.stop();
        runT1 = tt.GetRuntime();
        cout<<"Partition index construction time: "<<runT1<<" s"<<endl;


        /// Overlay graph construction
        tt.start();
//        Construct_OverlayGraph(true, false);//all-pair boundary shortcut
        Construct_OverlayGraphNoAllPair(true,false);
        tt.stop();
        runT2 = tt.GetRuntime();
        cout<<"Overlay graph construction time: "<<runT2<<" s."<< endl;

        /// Overlay index construction
        tt.start();
        Construct_OverlayIndex(true);
        tt.stop();
        runT3 = tt.GetRuntime();
        cout<<"Overlay index construction time: "<<runT3<<" s."<< endl;

        /// Partition index repair
        if(strategy==PostBoundary){
            tt.start();
            ConstructPartitionPost(true, false);
            tt.stop();
            runT4 += tt.GetRuntime();
            cout<<"Time for post partition construction: "<<tt.GetRuntime()<<" s."<<endl;

            tt.start();
            ConstructPartitionPostIndex(true, true);
            tt.stop();
            runT4 += tt.GetRuntime();
            cout<<"Time for post partition index construction: "<<tt.GetRuntime()<<" s."<<endl;

//        tt.start();
//        Repair_PartiIndex(true);
////        Repair_PartiIndex(false);
//        tt.stop();
//        runT3 = tt.GetRuntime();
//        cout<<"Partition index repair time: "<<runT4<<" s"<<endl;
        }
//    Compute_tree_label(ifParallel, ifOpt);//Construct periphery index (H2H label + interface label)
    }



    cout<<"Overall Construction Time: "<<runT1+runT2+runT3+runT4+runT5<<" s."<<endl;

    IndexSizePH2H();//index (labeling+pruning point) size computation
}

void Graph::PPLLIndexConstruct(int strategy) {
    double runT1, runT2, runT3, runT4, runT5;
    runT1=0, runT2=0, runT3=0, runT4=0, runT5=0;

    /// Read order and partitions
    string orderfile=graphfile+".order";
//    orderfile=graphfile+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_order";
//    orderfile=graphfile+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_orderMDE";
//    orderfile=graphfile+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_orderMDE2";
//#ifdef __APPLE__
////    cout<<"The platform is macOS."<<endl;
//#else
//    orderfile="/home/data/xzhouby/datasets/"+dataset+"/"+dataset+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_orderMDE2";
//#endif
    ReadOrder(orderfile);

    string partitionfile=graphfile+"_"+algoParti+"_"+to_string(partiNum);
    GraphPartitionRead(partitionfile);//read partitions

    vSm.reserve(node_num);
    for(int i = 0; i < node_num; i++)
    {
        Semaphore* s = new Semaphore(1);
        vSm.push_back(s);
    }
    BPCLBSize=threadnum*10;
    Timer tt;

    if(strategy==PreBoundary){
        /// All-pair boundary shortcut computation
        tt.start();
        AllPairBoundaryDisCompute(true);
        tt.stop();
        runT1 = tt.GetRuntime();
        cout<<"All-pair boundary shortcut construction time: "<<runT1<<" s"<<endl;

        /// Overlay index construction
        tt.start();
        ConstructPLL_OverlayIndex(NeighborsOverlayV);
        tt.stop();
        runT3 = tt.GetRuntime();
        cout<<"Overlay index construction time: "<<runT3<<" s."<< endl;

        /// Partition index and construction
        tt.start();
        ConstructPLL_PartiIndex(true);
        tt.stop();
        runT2 = tt.GetRuntime();
        cout<<"Partition index construction time: "<<runT2<<" s"<<endl;


    }
    else if(strategy>=NoBoundary){
        tt.start();
        /// Partition index and Overlay graph construction
        ConstructPLL_PartiIndex(true);
//        ConstructPLL_PartiIndex(false);
        tt.stop();
        runT1 = tt.GetRuntime();
        cout<<"Partition index construction time: "<<runT1<<" s"<<endl;

        /// Overlay graph construction
        tt.start();
        ConstructPLL_OverlayGraph(true);
//        ConstructPLL_OverlayGraph(false);
        tt.stop();
        runT2 = tt.GetRuntime();
        cout<<"Overlay graph construction time: "<<runT2<<" s."<< endl;

        /// Overlay index construction
        tt.start();
        ConstructPLL_OverlayIndex(NeighborsOverlayV);
        tt.stop();
        runT3 = tt.GetRuntime();
        cout<<"Overlay index construction time: "<<runT3<<" s."<< endl;


        /// Partition index repair
        if(strategy==PostBoundary){
            NeighborsPartiPostV=NeighborsParti;
            tt.start();
            ConstructPartitionPostPLL(true);
            tt.stop();
            runT4 += tt.GetRuntime();
            cout<<"Time for post partition construction: "<<tt.GetRuntime()<<" s."<<endl;

            tt.start();
            ConstructPLL_PartiIndexPost(true);
            tt.stop();
            runT4 += tt.GetRuntime();
            cout<<"Time for post partition index construction: "<<tt.GetRuntime()<<" s."<<endl;

//        tt.start();
//        Repair_PartiIndex(true);
////        Repair_PartiIndex(false);
//        tt.stop();
//        runT3 = tt.GetRuntime();
//        cout<<"Partition index repair time: "<<runT4<<" s"<<endl;
        }
//    Compute_tree_label(ifParallel, ifOpt);//Construct periphery index (H2H label + interface label)
    }



    cout<<"Overall Construction Time: "<<runT1+runT2+runT3+runT4+runT5<<" s."<<endl;

    IndexSizePPLL();//index (labeling+pruning point) size computation
}

//function for computing the index size
void Graph::IndexSizePH2H(){
    unsigned long long m=0,m1=0,m2=0,m3=0,m4=0,m5=0;

    //Overlay index
    for(int i=0;i<Tree.size();i++){
        m1+=Tree[i].dis.size()*sizeof(int);//dis
        m1+=Tree[i].pos.size()*sizeof(int);//pos
        m2+=Tree[i].dis.size()*sizeof(int);//cnt
        m2+=Tree[i].vert.size()*3*sizeof(int);//neighID/weight/count
    }
    //overlay
    for(int i=0;i< SCconNodesMT.size();i++){
        for(auto it=SCconNodesMT[i].begin(); it!=SCconNodesMT[i].end(); it++){
            m2+=sizeof(int)+(*it).second.size()*2*sizeof(int);
        }
    }

    if(PSPStrategy>=PostBoundary){
        if(algoChoice==H2H){
            for(int pid=0;pid<TreesPost.size();++pid){
                for(int i=0;i<TreesPost[pid].size();i++){
                    m3+=TreesPost[pid][i].dis.size()*sizeof(int);//dis
                    m3+=TreesPost[pid][i].pos.size()*sizeof(int);//pos
                    m4+=TreesPost[pid][i].dis.size()*sizeof(int);//cnt
                    m4+=TreesPost[pid][i].vert.size()*3*sizeof(int);//neighID/weight/count
                }
            }

        }else if(algoChoice==CH){
            for(int pid=0;pid<TreesPost.size();++pid){
                for(int i=0;i<TreesPost[pid].size();i++){
//                    m3+=TreesPost[pid][i].dis.size()*sizeof(int);//dis
//                    m3+=TreesPost[pid][i].pos.size()*sizeof(int);//pos
//                    m4+=TreesPost[pid][i].dis.size()*sizeof(int);//cnt
                    m4+=TreesPost[pid][i].vert.size()*3*sizeof(int);//neighID/weight/count
                }
            }
        }

        for(int i=0;i< SCconNodesMTPost.size();i++){
            for(auto it=SCconNodesMTPost[i].begin(); it!=SCconNodesMTPost[i].end(); it++){
                m4+=sizeof(int)+(*it).second.size()*2*sizeof(int);
            }
        }

    }
    else if(PSPStrategy==NoBoundary || PSPStrategy==PreBoundary){
        if(algoChoice==H2H){
            //partitions
            for(int pid=0;pid<Trees.size();++pid){
                for(int i=0;i<Trees[pid].size();i++){
                    m3+=Trees[pid][i].dis.size()*sizeof(int);//dis
                    m3+=Trees[pid][i].pos.size()*sizeof(int);//pos
                    m4+=Trees[pid][i].dis.size()*sizeof(int);//cnt
                    m4+=Trees[pid][i].vert.size()*3*sizeof(int);//neighID/weight/count
                }
            }
        }else if(algoChoice==CH){
            for(int pid=0;pid<Trees.size();++pid){
                for(int i=0;i<Trees[pid].size();i++){
//                    m3+=Trees[pid][i].dis.size()*sizeof(int);//dis
//                    m3+=Trees[pid][i].pos.size()*sizeof(int);//pos
//                    m4+=Trees[pid][i].dis.size()*sizeof(int);//cnt
                    m4+=Trees[pid][i].vert.size()*3*sizeof(int);//neighID/weight/count
                }
            }
        }
        //partitions
        for(int i=0;i< SCconNodesMTP.size();i++){
            for(auto it=SCconNodesMTP[i].begin(); it!=SCconNodesMTP[i].end(); it++){
                m4+=sizeof(int)+(*it).second.size()*2*sizeof(int);
            }
        }
    }



    //cout<<"Index size "<<(double)m1/1024/1024<<", Pruning point size "<<(double)m2/1024/1024<<endl;
    m=m1+m2+m3+m4;
    cout<<"Distance labeling size: "<<(double)(m1+m3)/1024/1024<<" MB"<<endl;
    cout<<"Overlay graph index size: "<<(double)(m1+m2)/1024/1024<<" MB"<<endl;
    cout<<"Partition graphs index size "<<(double)(m3+m4)/1024/1024<<" MB"<<endl;
    cout<<"Overall index size "<<(double)m/1024/1024<<" MB"<<endl;
}

void Graph::IndexSizePPLL(){
    unsigned long long m=0,m1=0,m2=0,m3=0,m4=0,m5=0;

    //Overlay index
    for(int i=0;i<Label.size();i++){
        m1+=Label[i].size()*sizeof(int)*2;//dis
    }
    //overlay
    for(int i=0;i< PruningPointSet.size();i++){
        for(auto it=PruningPointSet[i].begin(); it!=PruningPointSet[i].end(); it++){
            m2+=sizeof(int)+(*it).second.size()*2*sizeof(int);
        }
    }

    if(PSPStrategy>=PostBoundary){
        for(int i=0;i<LabelsPost.size();++i){
            m3+=LabelsPost[i].size()*sizeof(int)*2;//dis
        }
        for(int i=0;i< PruningPointSetPost.size();i++){
            for(auto it=PruningPointSetPost[i].begin(); it!=PruningPointSetPost[i].end(); it++){
                m4+=sizeof(int)+(*it).second.size()*2*sizeof(int);
            }
        }

    }
    else if(PSPStrategy==NoBoundary || PSPStrategy==PreBoundary){
        for(int i=0;i<Labels.size();++i){
            m3+=Labels[i].size()*sizeof(int)*2;//dis
        }
        for(int i=0;i< PruningPointSetP.size();i++){
            for(auto it=PruningPointSetP[i].begin(); it!=PruningPointSetP[i].end(); it++){
                m4+=sizeof(int)+(*it).second.size()*2*sizeof(int);
            }
        }
    }



    //cout<<"Index size "<<(double)m1/1024/1024<<", Pruning point size "<<(double)m2/1024/1024<<endl;
    m=m1+m2+m3+m4;
    cout<<"Distance labeling size: "<<(double)(m1+m3)/1024/1024<<" MB"<<endl;
    cout<<"PruningPoint size: "<<(double)(m2+m4)/1024/1024<<" MB"<<endl;
    cout<<"Overlay graph index size: "<<(double)(m1+m2)/1024/1024<<" MB"<<endl;
    cout<<"Partition graphs index size "<<(double)(m3+m4)/1024/1024<<" MB"<<endl;
    cout<<"Overall index size "<<(double)m/1024/1024<<" MB"<<endl;
}

void Graph::PH2HVertexOrdering(int type){
    ReadGraph(graphfile);
    int pNum=partiNum;

//#ifdef __APPLE__
//    cout<<"The platform is macOS."<<endl;
//#else
//    string dirname="/home/data/xzhouby/datasets/"+dataset;
//    struct stat sb;
//    if (stat(dirname.c_str(), &sb) == 0){
//        dirname="/home/data/xzhouby/datasets/"+dataset+"/"+dataset+"_"+algoParti+"_"+ to_string(pNum);
//        mkdir(dirname.c_str(), 0777);
//    }
//    else{
//        cout << "The source path is invalid! "<<dirname<<endl;
//        exit(1);
//    }
//#endif

    switch (type) {
        case 0:{//MDE partition + distant MDE overlay
            cout<<"MDE ordering."<<endl;
            SketchGraphBuild();
            OverlayOrderingBuild();
            PartitionOrderingBuildMDE(true);
//            PartitionOrderingBuildMDE(false);
            OrderingAssemblyMDE(pNum);
            break;
        }
        case 1:{
            cout<<"Boundary-first ordering."<<endl;
            SketchGraphBuild();
            OverlayOrderingBuildBoundaryFirst(node_num, NeighborSketch);
            PartitionOrderingBuildMDE(true);
            OrderingAssemblyBoundaryFirst(pNum);
            break;
        }
        case 2:{
            cout<<"Boundary-first MDE ordering."<<endl;
            SketchGraphBuild();
            OverlayOrderingBuild();
            PartitionOrderingBuildMDE(true);
            OrderingAssemblyMDEBoundaryFirst(pNum);
            break;
        }
        default:{
            cout<<"Wrong ordering type! "<<type<<endl; exit(1);
        }
    }
    exit(0);
}
void Graph::OrderingAssemblyMDEBoundaryFirst(int pNum){
    string filename=graphfile+"_"+algoParti+"_"+to_string(pNum)+"/vertex_orderMDE2";
//#ifdef __APPLE__
////    cout<<"The platform is macOS."<<endl;
//#else
//    filename="/home/data/xzhouby/datasets/"+dataset+"/"+dataset+"_"+algoParti+"_"+ to_string(pNum)+"/vertex_orderMDE2";
//#endif

    ofstream OF(filename,ios::out);
    if(!OF.is_open()){
        cout<<"Cannot open file "<<filename<<endl;
        exit(1);
    }
    int PID,ID;
    NodeOrder.assign(node_num,-1);
    vNodeOrder.clear();
    int order_i=0;
    set<int> vertices;
    /// For non-boundary vertex
    for(int pid=0;pid<partiNum;++pid){
        PID=vNodeOrderOverlay[pid];
//        cout<<pid<<": "<<PID<<endl;
        for(auto it=vNodeOrderParti[PID].begin();it!=vNodeOrderParti[PID].end();++it){
//            cout<<PID<<": "<<it->first<<" "<<it->second<<endl;
//            ID=it->second;
            ID=*it;
            if(!PartiTag[ID].second){// if not boundary
                vNodeOrder.emplace_back(ID);
                vertices.insert(ID);
                ++order_i;
            }
        }
    }
    /// For boundary vertex
    for(int pid=0;pid<partiNum;++pid){
        PID=vNodeOrderOverlay[pid];
//        cout<<pid<<": "<<PID<<endl;
        for(auto it=vNodeOrderParti[PID].begin();it!=vNodeOrderParti[PID].end();++it){
//            cout<<PID<<": "<<it->first<<" "<<it->second<<endl;
//            ID=it->second;
            ID=*it;
            if(PartiTag[ID].second){// if boundary
                vNodeOrder.emplace_back(ID);
                vertices.insert(ID);
                ++order_i;
            }
        }
    }

//    cout<<verticesV.size()<<" "<<vertices.size()<<endl;
    if(order_i!=node_num || vertices.size()!=node_num || vNodeOrder.size()!=node_num){
        cout<<"Wrong order number! "<<order_i<<" "<<vertices.size()<<" "<<vNodeOrder.size()<<" "<<node_num<<endl; exit(1);
    }



    for(int i=0;i<node_num;++i){
        NodeOrder[vNodeOrder[i]]=i;
    }

    for(int i=0;i<nodeNumOverlay;++i){
        ID=vNodeOrder[node_num-i-1];
        if(!PartiTag[ID].second){
            cout<<nodeNumOverlay<<", "<<node_num-i<<": "<<ID<<"("<<NodeOrder[ID]<<") is not boundary vertex!"<<endl; exit(1);
        }
    }

    OF<<node_num<<endl;
    for(int i = 0; i < NodeOrder.size(); i++){
        if(NodeOrder[i]==-1){
            cout<<"Wrong order! "<<i<<"("<<PartiTag[i].first<<") "<<NodeOrder[i]<<endl; exit(1);
        }
        OF << i << " " << NodeOrder[i] << endl;
    }
    OF.close();
    cout<<"Finish "<<endl;
}
//function of MDE ordering assemblying
void Graph::OrderingAssemblyMDE(int pNum){
    string filename=graphfile+"_"+algoParti+"_"+to_string(pNum)+"/vertex_orderMDE";

    ofstream OF(filename,ios::out);
    if(!OF.is_open()){
        cout<<"Cannot open file "<<filename<<endl;
        exit(1);
    }
    int PID,ID;
    NodeOrder.assign(node_num,-1);
    int order_i=0;
    set<int> vertices;
    for(int pid=0;pid<partiNum;++pid){
        PID=vNodeOrderOverlay[pid];
//        cout<<pid<<": "<<PID<<endl;
        for(auto it=vNodeOrderParti[PID].begin();it!=vNodeOrderParti[PID].end();++it){
//            cout<<PID<<": "<<it->first<<" "<<it->second<<endl;
//            ID=it->second;
            ID=*it;
            vertices.insert(ID);
            NodeOrder[ID]=order_i;
            ++order_i;
        }
    }
//    cout<<verticesV.size()<<" "<<vertices.size()<<endl;
    if(order_i!=node_num){
        cout<<"Wrong order number! "<<order_i<<" "<<node_num<<endl; exit(1);
    }
    OF<<node_num<<endl;
    for(int i = 0; i < NodeOrder.size(); i++){
        if(NodeOrder[i]==-1){
            cout<<"Wrong order! "<<i<<"("<<PartiTag[i].first<<") "<<NodeOrder[i]<<endl; exit(1);
        }
        OF << i << " " << NodeOrder[i] << endl;
    }
    OF.close();
    cout<<"Finish "<<endl;
}
//function of boundary-first assemblying
void Graph::OrderingAssemblyBoundaryFirst(int pNum){
    string orderfile=graphfile+"_"+algoParti+"_"+to_string(pNum)+"/vertex_order2";
    set<int> vcheck;//use to check the redundant ordered vertex
    vcheck.clear();
    vNodeOrder.clear();
    int pid,ID;
    //for vertex within partition
    for(int k=0;k<vNodeOrderOverlay.size();k++){
        pid=vNodeOrderOverlay[k];

        for(auto it=vNodeOrderParti[k].begin();it!=vNodeOrderParti[k].end();++it){
//            ID=it->second;
            ID=*it;
            if(!PartiTag[ID].second){//if not boundary vertex
                vNodeOrder.push_back(ID);
                if(vcheck.find(ID)!=vcheck.end())
                    cout<<"wrong: redundant vertex ordered"<<endl;
                vcheck.insert(ID);
            }
        }
    }

    //for boundary vertex
    for(int k=0;k<vNodeOrderOverlay.size();k++){
        pid=vNodeOrderOverlay[k];
        for(int i=0;i<BoundVertex[pid].size();i++){
            ID=BoundVertex[pid][i];
            vNodeOrder.push_back(ID);
            if(vcheck.find(ID)!=vcheck.end())
                cout<<"wrong: redundant vertex ordered"<<endl;
            vcheck.insert(ID);
        }
    }

    //cout<<"total number of ordered vertex "<<vNodeOrder.size()<<endl;
    if(vNodeOrder.size()!=node_num)
        cout<<"Something wrong happened: some vertices do not have the vertex order!"<<endl;

    NodeOrder.assign(node_num,-1);
    for(int k=0;k<vNodeOrder.size();k++){
        NodeOrder[vNodeOrder[k]]=k;
    }

    ofstream OF(orderfile);
    if(!OF){
        cout<<"Cannot open file "<<orderfile<<endl;
        exit(1);
    }
    OF<<NodeOrder.size()<<endl;
    for(int i=0;i<NodeOrder.size();i++){
        OF<<i<<" "<<NodeOrder[i]<<endl;
    }
    OF.close();
    cout<<"Finished."<<endl;
}

void Graph::SketchOrder(vector<vector<pair<int,int>>> Neighbor, vector<int> &vNodeOrderSketch){
    map<int,pair<int,int>> m;
    E.assign(node_num,m);
    for(int i=0;i<Neighbor.size();i++){
        for(int j=0;j<Neighbor[i].size();j++)
            E[i].insert(make_pair(Neighbor[i][j].first,make_pair(1,1)));
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

    vector<vector<pair<int,pair<int,int>>>> NeighborCon;
    vector<pair<int,pair<int,int>>> vect;
    NeighborCon.assign(node_num,vect); //NeighborCon.clear();
    //SCconNodes.clear();

    //cout<<"Begin to contract"<<endl;
    int count=0;

    while(!Deg.empty()){
        //if(count%10==0) cout<<"count "<<count<<endl;
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

        vNodeOrderSketch.push_back(x);
        Deg.erase(Deg.begin());
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
            deleteEOrderGenerate(x,y);
            change[y]=true;
        }

        for(int i=0;i<Neigh.size();i++){
            for(int j=i+1;j<Neigh.size();j++){
                insertEOrderGenerate(Neigh[i].first,Neigh[j].first,Neigh[i].second.first+Neigh[j].second.first);
                change[Neigh[i].first]=true;
                change[Neigh[j].first]=true;
            }
        }
    }

    /*NodeOrderSketch.assign(nodenum,-1);
    for(int k=0;k<vNodeOrderSketch.size();k++){
        NodeOrderSketch[vNodeOrderSketch[k]]=k;
        cout<<"order "<<k<<", vertex "<<vNodeOrderSketch[k]<<", degree "<<NeighborSketch[vNodeOrderSketch[k]].size()<<endl;
    }*/
    //cout<<"Finish Contract"<<endl;
}

void Graph::SketchGraphBuild(){
    NeighborsParti.assign(node_num, vector<pair<vertex,int>>());
    NeighborsOverlay.assign(node_num,unordered_map<vertex,int>());
//    NeighborsOverlay.assign(node_num,vector<pair<vertex,int>>());
    PartiTag.assign(node_num, make_pair(-1,false));

    bool flag_minus = false;

    string filename=graphfile+"_"+algoParti+"_"+to_string(partiNum);
    ifstream IF1(filename+"/subgraph_vertex");
    if(!IF1){
        cout<<"Cannot open file "<<"subgraph_vertex"<<endl;
        exit(1);
    }

    int pnum2;
    IF1>>pnum2;
    if(algoParti == "NC"){
//        flag_minus = true;
        partiNum = pnum2;
    }else if(algoParti == "SC" || algoParti == "MT"){
//        flag_minus = true;
//        pnum2 = pnum;
    }
    cout<<"Partition number: "<<pnum2<<endl;

    PartiVertex.assign(partiNum,vector<vertex>());
    for(int k=0;k<pnum2;k++){
        int vernum,ID;
        IF1>>vernum;
        for(int i=0;i<vernum;i++){
            IF1>>ID;
//            if(flag_minus){
//                ID = ID-1;
//            }

            if(ID>=0 && ID<node_num){
                if(PartiTag[ID].first==-1){
                    PartiTag[ID].first=k;
                    PartiVertex[k].emplace_back(ID);
                }else{
                    cout<<"vertex already in one partition!"<<ID<<" "<<PartiTag[ID].first<<" "<<k<<endl;
                }
            }else{
                cout<<"Wrong vertex ID! "<<ID<<endl; exit(1);
            }

        }
    }
    //further check that each vertex is in one and only one partition
    for(int vid=0;vid<node_num;vid++){
        if(PartiTag[vid].first==-1){
            cout<<"vertex "<<vid<<" not within any partition"<<endl; exit(1);
        }
    }
    int nNum=0;
    for(int pid=0;pid<partiNum;++pid){
        nNum+=PartiVertex[pid].size();
    }
    if(nNum!=node_num){
        cout<<"Inconsistent node number! "<<nNum<<" "<<node_num<<endl; exit(1);
    }
    //record the vertex to PartiVertex in vertex order: from lower-rank vertex to higher-rank vertex

    ifstream IF(filename+"/subgraph_edge");
    if(!IF){
        cout<<"Cannot open file "<<"subgraph_edge"<<endl;
        exit(1);
    }

    int pnum1;
    IF>>pnum1;
    for(int k=0;k<pnum1;k++){
        int edgenum0,ID1,ID2,weight;
        IF>>edgenum0;
        for(int i=0;i<edgenum0;i++){
            IF>>ID1>>ID2>>weight;
//            if(flag_minus){
//                ID1 = ID1-1; ID2 = ID2-1;
//            }
            if(ID1>=0 && ID1 <node_num && ID2>=0 && ID2 <node_num && weight>0){
                NeighborsParti[ID1].emplace_back(ID2,weight);
            }else{
                cout<<"Wrong for subgraph_edge! "<<ID1<<" "<<ID2<<" "<<weight<<endl; exit(1);
            }


        }
    }

    vector<int> vecint;
    vecint.clear();
    NeighborSketch.assign(partiNum, vecint);
//    vector<set<int>> NeighborSketchS;
    NeighborSketchS.assign(partiNum, set<int>());

    BoundVertex.assign(partiNum,vector<vertex>());
    //read the cut edges
    ifstream IF2(filename+"/cut_edges");
    if(!IF2){
        cout<<"Cannot open file "<<"cut_edges"<<endl;
        exit(1);
    }

    int ednum,ID1,ID2,weight;
    int boundaryNum=0;
    int PID1, PID2;

    IF2>>ednum;
    for(int i=0;i<ednum;i++){
        IF2>>ID1>>ID2>>weight;
//        if(flag_minus){
//            ID1 = ID1-1; ID2 = ID2-1;
//        }

        if(ID1>=0 && ID1 <node_num && ID2>=0 && ID2 <node_num && weight>0){
            PID1=PartiTag[ID1].first, PID2=PartiTag[ID2].first;
            if(PartiTag[ID1].first==PartiTag[ID2].first){
                cout<<"two end points of cut edge are in the same partition"<<endl; exit(1);
            }
            if(!PartiTag[ID1].second){
                PartiTag[ID1].second=true;
                boundaryNum++;
                BoundVertex[PartiTag[ID1].first].emplace_back(ID1);
            }
            if(!PartiTag[ID2].second){
                PartiTag[ID2].second=true;
                boundaryNum++;
                BoundVertex[PartiTag[ID2].first].emplace_back(ID2);
            }

            NeighborsOverlay[ID1].insert({ID2,weight});
//            NeighborsOverlay[ID1].emplace_back(ID2,weight);

            if(NeighborSketchS[PID1].find(PID2)==NeighborSketchS[PID1].end()){//if not found, i.e., PID2 is not in the NeighborSketchS of PID1
                NeighborSketch[PID1].push_back(PID2);
                NeighborSketch[PID2].push_back(PID1);

                NeighborSketchS[PID1].insert(PID2);
                NeighborSketchS[PID2].insert(PID1);
            }
        }else{
            cout<<"Wrong for cut_edge! "<<ID1<<" "<<ID2<<" "<<weight<<endl; exit(1);
        }
    }

    nodeNumOverlay=boundaryNum;


    /*for(int k=0;k<pnum;k++){
        cout<<k<<" "<<NeighborSketch[k].size()<<endl;
    }*/
}

void Graph::OverlayOrderingBuild(){

    int lastParti = -1;
    map<int,pair<int,int>> m;
    E.assign(partiNum,m);
    for(int i=0;i<NeighborSketch.size();i++){
        for(auto it=NeighborSketch[i].begin();it!=NeighborSketch[i].end();++it){
            E[i].insert(make_pair(*it,make_pair(1,1)));
        }
    }
    _DD_.assign(partiNum,0);
    DD.assign(partiNum,0);
    _DD2_.assign(partiNum,0);

    set<DegComp> Deg;
    set<DegComp2> Deg2;
    int ID,degree;
    for(ID=0;ID<partiNum;ID++){
        degree=NeighborSketch[ID].size();
        if(degree!=0){
            _DD_[ID]=degree;
            DD[ID]=degree;
            Deg.insert(DegComp(ID));
        }else{
            _DD_[ID]=degree;
            DD[ID]=degree;
            Deg.insert(DegComp(ID));
            cout<<"Sketch graph build. Not single CC! degree is zero. "<<ID<<" "<<degree<<endl;
//            exit(1);
        }
    }

    vector<bool> exist; exist.assign(partiNum,true);
    vector<bool> change; change.assign(partiNum,false);

    vector<set<int>> NeighborCon(partiNum,set<int>());
    vector<int> neix;
    neix.clear();

    int count=0;
    int Twidth=0;
    int order_i=0;
    int ID1, ID2;
    int x;

    set<int> pSet;
    while(!Deg.empty()){

        x=(*Deg.begin()).x;
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
        Deg.erase(Deg.begin());



        if(lastParti!=-1 && NeighborCon[lastParti].find(x) != NeighborCon[lastParti].end()){//if x is lastParti's neighbor
            neix.clear();
            while(NeighborCon[lastParti].find(x) != NeighborCon[lastParti].end()){//if found
                _DD2_[x]++;
                neix.emplace_back(x);

                if(Deg.empty()){
                    break;
                }else{
                    x=Deg.begin()->x;
                    Deg.erase(Deg.begin());
                }
            }

            if(NeighborCon[lastParti].find(x) != NeighborCon[lastParti].end()){//if found, i.e., x is the neighbor of lastParti
                if(neix.size()>1){//if neix has more than one element
                    if(Deg.empty()){
                        Deg2.clear();
                        for(int i=0;i<neix.size();++i){
                            Deg2.insert(DegComp2(neix[i]));
                        }
                        x=Deg2.begin()->x;
                        Deg2.erase(Deg2.begin());
                        if(!Deg2.empty()){
                            for(auto it=Deg2.begin();it!=Deg2.end();++it){
                                ID=it->x;
                                Deg.insert(DegComp(ID));
                            }
                            Deg2.clear();
                        }
                    }else{
                        cout<<"Wrong! "<<endl; exit(1);
                    }
                }
            }//if not the neighbor
            else{
                if(!neix.empty()){
                    for(int i=0;i<neix.size();++i){
                        Deg.insert(DegComp(neix[i]));
                    }
                }
            }
        }

//        cout<<x<<" "<<Deg.size()<<endl;
        vNodeOrderOverlay.emplace_back(x);
        pSet.insert(x);
        lastParti = x;
        exist[x]=false;

        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();

        for(auto it=E[x].begin();it!=E[x].end();it++){
            if(exist[(*it).first]){
                Neigh.emplace_back(*it);
                NeighborCon[x].insert(it->first);
            }
        }

        if(Neigh.size()>Twidth)
            Twidth=Neigh.size();

        //multi threads for n^2 combination
        for(int i=0;i<Neigh.size();i++){
            int y=Neigh[i].first;
            deleteEOrderGenerate(x,y);//delete x from y's neighbor
            change[y]=true;
        }

        for(int i=0;i<Neigh.size();++i){
            ID1=Neigh[i].first;
            for(int j=i+1;j<Neigh.size();++j){
                ID2=Neigh[j].first;
                insertEOrderGenerate(ID1, ID2, 1);//for order generation, only need unweighted graph
            }
        }
    }
    if(vNodeOrderOverlay.size() != partiNum || pSet.size()!=partiNum){
        cout<<"Inconsistent size for sketch graph! "<<vNodeOrderOverlay.size()<<" "<<pSet.size() <<" "<< partiNum<<endl; exit(1);
    }

//    exit(0);
}

//based on SketchOrder function, to guarantee the order of neighboring vertices are not contiguous
void Graph::OverlayOrderingBuildBoundaryFirst(int nodenum, vector<vector<int>> Neighbor){
    map<int,pair<int,int>> m;
    E.assign(partiNum,m);
    for(int i=0;i<Neighbor.size();i++){
        for(int j=0;j<Neighbor[i].size();j++)
            E[i].insert(make_pair(Neighbor[i][j],make_pair(1,1)));
    }

    _DD_.assign(partiNum,0);
    DD.assign(partiNum,0);

    set<DegComp> Deg;
    int degree;
    for(int i=0;i<Neighbor.size();i++){
        degree=Neighbor[i].size();
        if(degree!=0){
            _DD_[i]=degree;
            DD[i]=degree;
            Deg.insert(DegComp(i));
        }
    }

    vector<bool> exist; exist.assign(partiNum,true);
    vector<bool> change; change.assign(partiNum,false);

    vector<vector<pair<int,pair<int,int>>>> NeighborCon;
    vector<pair<int,pair<int,int>>> vect;
    NeighborCon.assign(partiNum,vect); //NeighborCon.clear();
    //SCconNodes.clear();

    //cout<<"Begin to contract"<<endl;
    int count=0;
    int x;
    int lastx;
    vector<int> neix;
    neix.clear();
    while(!Deg.empty()){
        //if(count%10==0) cout<<"count "<<count<<endl;
        count+=1;

        while(true){
            if(Deg.empty()){//in case that all the remaining vertices are the neighbors of current vertex x
                x=neix[0];
                for(int j=1;j<neix.size();j++){
                    Deg.insert(DegComp(neix[j]));
//					cout<<"insert back/// "<<neix[j]<<endl;
                }
                neix.clear();
                break;///
            }
            else
                x=(*Deg.begin()).x;

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

            if(count==1){
                lastx=x;
                break;
            }else if(NeighborSketchS[lastx].find(x)==NeighborSketchS[lastx].end()){//if not found
                lastx=x;
                break;
            }else{
                Deg.erase(DegComp(x));
                neix.push_back(x);
//				cout<<"erase "<<x<<endl;
            }

//            if(Deg.empty())////
//                break;
        }

        if(neix.size()!=0){
            for(int j=0;j<neix.size();j++){
                Deg.insert(DegComp(neix[j]));
//				cout<<"insert back "<<neix[j]<<endl;
            }
        }
        neix.clear();

        vNodeOrderOverlay.push_back(x);
        Deg.erase(DegComp(x));
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
            deleteEOrderGenerate(x,y);
            change[y]=true;
        }

        for(int i=0;i<Neigh.size();i++){
            for(int j=i+1;j<Neigh.size();j++){
                insertEOrderGenerate(Neigh[i].first,Neigh[j].first,Neigh[i].second.first+Neigh[j].second.first);
                change[Neigh[i].first]=true;
                change[Neigh[j].first]=true;
            }
        }
    }

    /*NodeOrderSketch.assign(nodenum,-1);
    for(int k=0;k<vNodeOrderSketch.size();k++){
        NodeOrderSketch[vNodeOrderSketch[k]]=k;
        cout<<"order "<<k<<", vertex "<<vNodeOrderSketch[k]<<", degree "<<NeighborSketch[vNodeOrderSketch[k]].size()<<endl;
    }*/
    //cout<<"Finish Contract"<<endl;
}
//MDE-based vertex ordering for partition
void Graph::PartitionOrderingBuildMDE(bool ifParallel){
    //initialize E
    map<int,pair<int,int>> m;
    E.assign(node_num,m);
    for(int i=0;i<NeighborsParti.size();i++){
        for(auto it=NeighborsParti[i].begin();it!=NeighborsParti[i].end();++it){
            E[i].insert(make_pair(it->first,make_pair(it->second,1)));
        }
    }
    _DD_.assign(node_num,0);
    DD.assign(node_num,0);
//    vNodeOrderParti.assign(partiNum,map<int,int>());
    vNodeOrderParti.assign(partiNum,vector<int>());

    if(ifParallel){
        if(threadnum<partiNum){
            vector<vector<int>> processID;
            processID.assign(threadnum, vector<int>());
            vector<int> vertices;
            for(int pid=0;pid<partiNum;++pid){
                vertices.emplace_back(pid);
            }
            ThreadDistribute(vertices, processID);
            cout<<"Batch number: "<<processID[0].size()<<endl;
            boost::thread_group thread;
            for(auto j=0;j<processID.size();++j){
                thread.add_thread(new boost::thread(&Graph::PartitionOrderingV, this, boost::ref(processID[j]) ));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                thread.add_thread(new boost::thread(&Graph::PartitionOrdering, this, j));
            }
            thread.join_all();
        }
    }
    else{
        for(int i=0;i<partiNum;++i){
            PartitionOrdering(i);
        }
    }

}

void Graph::PartitionOrderingV(vector<int>& p){
    for(int i=0;i<p.size();++i){
        PartitionOrdering(p[i]);
    }
}

void Graph::PartitionOrdering(int pid){

    set<DegComp> Deg;
    int ID,degree;
    for(int i=0;i<PartiVertex[pid].size();i++){
        ID = PartiVertex[pid][i];
        degree=NeighborsParti[ID].size();
        if(degree!=0){
            _DD_[ID]=degree;
            DD[ID]=degree;
            Deg.insert(DegComp(ID));
        }else{
            _DD_[ID]=degree;
            DD[ID]=degree;
            Deg.insert(DegComp(ID));
            cout<<"PID "<<pid<<" . Not single CC! degree is zero. "<<ID<<" "<<degree<<endl;
//            exit(1);
        }
    }

    vector<bool> exist; exist.assign(node_num,true);
    vector<bool> change; change.assign(node_num,false);

    int count=0;
    int Twidth=0;
    int order_i=0;
    int ID1, ID2;
    while(!Deg.empty()){
//        if(count%10000==0)
//            cout<<"count "<<count<<" , treewidth "<<Twidth<<endl;
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

//        vNodeOrderParti[pid].insert({order_i,x});
        vNodeOrderParti[pid].push_back(x);
        order_i++;
        Deg.erase(Deg.begin());
        exist[x]=false;

        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();

        for(auto it=E[x].begin();it!=E[x].end();it++){
            if(exist[(*it).first]){
                Neigh.emplace_back(*it);
            }
        }

        if(Neigh.size()>Twidth)
            Twidth=Neigh.size();

        //multi threads for n^2 combination
        for(int i=0;i<Neigh.size();i++){
            int y=Neigh[i].first;
            deleteEOrderGenerate(x,y);//delete x from y's neighbor
            change[y]=true;
        }

        for(int i=0;i<Neigh.size();++i){
            ID1=Neigh[i].first;
            for(int j=i+1;j<Neigh.size();++j){
                ID2=Neigh[j].first;
                insertEOrderGenerate(ID1, ID2, 1);//for order generation, only need unweighted graph
            }
        }
    }
    if(vNodeOrderParti[pid].size() != PartiVertex[pid].size()){
        cout<<"Inconsistent size! "<< pid <<" "<<vNodeOrderParti[pid].size() <<" "<< PartiVertex[pid].size()<<endl; exit(1);
    }
}



/// Functions for MDE contraction
void Graph::deleteEOrderGenerate(int u,int v){//delete u from v's neighbor
    if(E[u].find(v)!=E[u].end()){
        E[u].erase(E[u].find(v));
        DD[u]--;
    }

    if(E[v].find(u)!=E[v].end()){
        E[v].erase(E[v].find(u));
        DD[v]--;
    }
}

void Graph::insertEOrderGenerate(int u,int v,int w){
    if(E[u].find(v)==E[u].end()){
        E[u].insert(make_pair(v,make_pair(w,1)));
        DD[u]++;
//        DD2[u]++;
    }
    if(E[v].find(u)==E[v].end()){
        E[v].insert(make_pair(u,make_pair(w,1)));
        DD[v]++;
//        DD2[u]++;
    }
}


/// Query Processing
//function for correctness check
void Graph::CorrectnessCheck(int runtimes){
    Timer tt;
    double runT=0;
    srand (time(NULL));
    int s, t, d1=INF, d2=INF, d3=INF;
//    runtimes = 1;
    cout<<"Correctness check ("<<runtimes<<" rounds) ...";
    for(int i=0;i<runtimes;i++){
//        if(i%100==0) cout<<i<<endl;
        s=rand()%node_num;
        t=rand()%node_num;
//        s=43465, t=46390;

        if(algoChoice==CH){
            if(PSPStrategy==PreBoundary || PSPStrategy==NoBoundary){
                tt.start();
                d2=QueryPCH(s,t,Trees);
                tt.stop();
            }else if(PSPStrategy==PostBoundary){
                tt.start();
                d2=QueryPCH(s,t,TreesPost);
                tt.stop();
            }

        }else if(algoChoice==H2H){
            tt.start();
            d2=QueryPH2H(s,t);
            tt.stop();
        }else if(algoChoice==PLL){
            tt.start();
            d2=QueryPPLL(s,t);
            tt.stop();
        }
        else if(algoChoice==0){
            tt.start();
            d2= Astar(s,t,Neighbor);
            tt.stop();
        }
        runT+=tt.GetRuntime();

        d1=Dijkstra(s,t,Neighbor);
//        cout<<"Algorithm "<<algoQuery<<", "<<i<<": "<<s<<"("<<NodeOrder[s]<<") "<<t<<"("<<NodeOrder[t]<<") "<<d2<<" "<<d1<<" ; Partition Tag: "<< PartiTag[s].first<<" "<<PartiTag[t].first<<"; Boundary Tag: "<<PartiTag[s].second<<" "<<PartiTag[t].second<<endl;
        if(d1!=d2){
            cout<<"Correct Test. InCorrect! Algorithm "<<algoChoice<<"+"<<PSPStrategy<<", "<<i<<": "<<s<<"("<<NodeOrder[s]<<") "<<t<<"("<<NodeOrder[t]<<") "<<d2<<" "<<d1<<" ; Partition Tag: "<< PartiTag[s].first<<" "<<PartiTag[t].first<<"; Boundary Tag: "<<PartiTag[s].second<<" "<<PartiTag[t].second<<endl;
            if(algoChoice==CH){
                if(PSPStrategy==PreBoundary || PSPStrategy==NoBoundary){
                    QueryPCHDebug(s,t,Trees);
                }else if(PSPStrategy==PostBoundary){
                    QueryPCHDebug(s,t,TreesPost);
                }
            }else if(algoChoice==H2H){
                QueryDebug(s,t);
            }


            exit(1);
        }
    }
    cout<<" Average Query Time: "<<1000*runT/runtimes<<" ms."<<endl;
}

//function for Query processing, debug version
int Graph::QueryDebug(int ID1, int ID2){
    int dis=INF;

    if(PartiTag[ID1].second && PartiTag[ID2].second){//Case 1: both in core
        cout<<"Core-Core"<<endl;
//        dis=QueryCoreDebug(ID1, ID2);

    }else if(PartiTag[ID1].second && !PartiTag[ID2].second){//Case 2: ID2 in partition, ID1 in core
        cout<<"Core-Parti"<<endl;
        dis=QueryPartiCoreDebug(ID2, ID1);

    }else if(!PartiTag[ID1].second && PartiTag[ID2].second){//Case 2: ID1 in partition, ID2 in core
        cout<<"Parti-Core"<<endl;
        dis=QueryPartiCoreDebug(ID1, ID2);

    }else if(!PartiTag[ID1].second && !PartiTag[ID2].second){//both in partition

        if(PartiTag[ID1].first != PartiTag[ID2].first){//Case 3: in different peripheries
            cout<<"Parti-Parti"<<endl;
            int d=INF;
            int b1,b2,d1,d2;//final results
            int pid1=PartiTag[ID1].first;
            int pid2=PartiTag[ID2].first;

            vector<int> B1=BoundVertex[pid1];
            vector<int> B2=BoundVertex[pid2];

            map<int,int> m1,m2;
            m1.clear();
            m2.clear();
            int bID1, bID2, tempdis;
            for(int i=0;i<B1.size();i++){
                bID1=B1[i];

//            m1.insert(make_pair(bID1,Tree[rank[ID1]].disInf[i]));
                m1.insert(make_pair(bID1, QueryH2HPartition(ID1,bID1,pid1)));
            }
            for(int j=0;j<B2.size();j++){
                bID2=B2[j];
//            m2.insert(make_pair(bID2,Tree[rank[ID2]].disInf[j]));
                m2.insert(make_pair(bID2, QueryH2HPartition(ID2,bID2,pid2)));
            }

            for(int k=0;k<B1.size();k++){
                bID1=B1[k];

                if(m1[bID1]>d)
                    continue;

                for(int z=0;z<B2.size();z++){
                    bID2=B2[z];

                    if(m2[bID2]>d)
                        continue;

                    tempdis=m1[bID1]+QueryCore(bID1,bID2)+m2[bID2];
                    if(tempdis<d){
                        d=tempdis;
                        b1=bID1; b2=bID2; d1=m1[bID1]; d2=m2[bID2];
//                        cout<<b1<<" "<<b2<<" "<<d<<endl;
                    }
                }
            }
            dis=d;
            int d_12=QueryCore(b1,b2), dDijk_s=Dijkstra(ID1,b1,Neighbor), dDijk_12=Dijkstra(b1,b2,Neighbor), dDijk_t=Dijkstra(b2,ID2,Neighbor);
            cout<<ID1<<" "<<b1<<"("<<NodeOrder[b1]<<","<<PartiTag[b1].first<<") "<<b2<<"("<<NodeOrder[b2]<<","<<PartiTag[b2].first<<") "<<ID2<<" : "<<d1<<" "<<d_12<<" "<<d2<<" ; "<<dDijk_s<<" "<<dDijk_12<<"("<<DijkstraCore(b1,b2)<<") "<<dDijk_t<<endl;

//            QueryCoreDebug(b1,b2);
//            QueryPartiPartiExtLCADebug(ID1,ID2);

//                if(d1!=dDijk_s){
//                    DijkstraPath(ID1,b1);
//                }
//                if(d_12!=dDijk_12){
//                    DijkstraPath(b1,b2);
//                }
//                if(d2!=dDijk_t){
//                    DijkstraPath(b2,ID2);
//                }

        }
        else{//Case 4: in the same periphery
            cout<<"Same-Parti"<<endl;
//                dis= QuerySameParti(ID1,ID2);
            int d=INF;
            int b1,b2,df1,df2;
            int pid1=PartiTag[ID1].first;
            int pid2=PartiTag[ID2].first;

            int temp_dis = QueryH2HPartition(ID1, ID2, pid1);/// d2 may be wrong sometimes
            if (temp_dis < d){
                d = temp_dis;//QueryH2HPartition(ID1,ID2,pid1);
                b1=b2=-1;
                df1=df2=-1;
            }

            vector<int> B = BoundVertex[pid1];
            map<int, int> m1, m2;
            m1.clear();
            m2.clear();
            vector<int> B1, B2;
            B1.clear();
            B2.clear();
            int bID, d1, d2;
            for (int i = 0; i < B.size(); i++) {
                bID = B[i];

                d1 = QueryH2HPartition(ID1,bID,pid1);
                d2 = QueryH2HPartition(ID2,bID,pid1);

                if (d1 < d) {
                    B1.push_back(bID);
                    m1.insert(make_pair(bID, d1));
                }
                if (d2 < d) {
                    B2.push_back(bID);
                    m2.insert(make_pair(bID, d2));
                }
            }

            int bID1, bID2, tempdis;
            if (!B1.empty() && !B2.empty()) {
                for (int k = 0; k < B1.size(); k++) {
                    bID1 = B1[k];
                    if (m1[bID1] > d)
                        continue;
                    for (int z = 0; z < B2.size(); z++) {
                        bID2 = B2[z];
                        if (m2[bID2] > d)
                            continue;
                        tempdis = m1[bID1] + QueryCore(bID1, bID2) + m2[bID2];
                        if (tempdis < d){
                            d = tempdis;
                            b1=bID1;b2=bID2;
                            df1=m1[bID1];df2=m2[bID2];
                        }
                    }
                }
            }

            if(b1!=-1){
                cout<<"d4: "<<ID1<<" "<<b1<<" "<<b2<<" "<<ID2<<" : "<<df1<<" "<<QueryCore(b1,b2)<<" "<<df2<<" ; "<<Dijkstra(ID1,b1,Neighbor)<<" "<<Dijkstra(b1,b2,Neighbor)<<" "<<Dijkstra(b2,ID2,Neighbor)<<endl;
            }else{
                int dDijk2 = Dijkstra(ID1,ID2,Neighbor);
                cout<<"d2: "<<d<<"; "<<dDijk2<<endl;
                if(d!=dDijk2){
//                        DijkstraPath(ID1,ID2);
                }
            }

            dis = d;

        }

    }
    return dis;
}

//function for core index correctness check
void Graph::CorrectnessCheckCore(int runtimes){
    srand (time(NULL));
    int s, t, d1, d2, d3;
    vector<int> coreVertex;
    for(int i=0;i<node_num;++i){
        if(PartiTag[i].second){
            coreVertex.emplace_back(i);
        }
    }
    int corenum=coreVertex.size();
    cout<<"Core graph correctness check ("<<runtimes<<" rounds)..."<<endl;
    for(int i=0;i<runtimes;i++){
        s=coreVertex[rand()%corenum];
        t=coreVertex[rand()%corenum];
        if(PartiTag[s].second && PartiTag[t].second){//for core vertex
            d1=QueryCore(s,t);
            d2=DijkstraCore(s,t);

            if(d1!=d2){
                cout<<"InCorrect! "<<s<<"("<<NodeOrder[s]<<") "<<t<<"("<<NodeOrder[t]<<"): "<<d1<<" "<<d2<<endl;
//				DijkstraPath(s,t);
//				DijkstraCorePath(s,t);
                exit(1);
            }
        }else
            i--;
    }
}

//function for efficiency test
void Graph::EffiCheck(string filename,int runtimes){
    ifstream IF(filename);
    if(!IF){
        cout<<"Cannot open Map "<<filename<<endl;
        exit(1);
    }
    cout<<"Query file: "<<filename<<endl;
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
    cout<<"Efficiency test. Run times: "<<runtimes<<endl;
    int s, t;
    Timer tt;

    double runT=0;
    int d1,d2;
    bool ifDebug=false;
//    ifDebug=true;

    if(ifDebug){
        cout<<"With correctness check."<<endl;
    }
    clock_t start = clock();
    vector<int> results(runtimes,-1);
    for(int i=0;i<runtimes;i++){
        ID1=ODpair[i].first, ID2=ODpair[i].second;
//        if(PartiTag[ID1].first!=PartiTag[ID2].first){
//            cout<<"Different Partition: "<<PartiTag[ID1].first<<" "<<PartiTag[ID2].first<<endl;
//        }

        if(algoChoice==CH){
            if(PSPStrategy==PreBoundary || PSPStrategy==NoBoundary){
                tt.start();
                d1=QueryPCH(ID1,ID2,Trees);
                tt.stop();
            }else if(PSPStrategy==PostBoundary){
                tt.start();
                d1=QueryPCH(ID1,ID2,TreesPost);
                tt.stop();
            }

        }else if(algoChoice==H2H){
            tt.start();
            d1=QueryPH2H(ID1,ID2);
            tt.stop();
        }else if(algoChoice==PLL){
            tt.start();
            d1=QueryPPLL(ID1,ID2);
            tt.stop();
        }else if(algoChoice==0){
            tt.start();
            d1=Astar(ID1,ID2,Neighbor);
//            d1=Dijkstra(ID1,ID2,Neighbor);
            tt.stop();
        }

        runT+=tt.GetRuntime();
        results[i]=d1;
        if(ifDebug){
            d2= Dijkstra(ID1,ID2,Neighbor);
            if(d1!=d2){
                if(algoChoice==H2H){
                    cout<<"Incorrect! "<<i<<": "<<ID1<<"("<<PartiTag[ID1].first<<") "<<ID2<<"("<<PartiTag[ID2].first<<"): "<<d1<<" "<<d2<<endl;
                    exit(1);
                }else if(algoChoice==CH){
                    cout<<"Incorrect! "<<i<<": "<<ID1<<"("<<PartiTag[ID1].first<<","<<PartiTag[ID1].second<<") "<<ID2<<"("<<PartiTag[ID2].first<<","<<PartiTag[ID2].second<<"): "<<d1<<" "<<d2<<endl; exit(1);
                }

            }
        }

    }


    cout<<"Average Query Time: "<<(double)runT*1000/runtimes<<" ms. "<<1000*(double)(clock() - start) / (CLOCKS_PER_SEC*runtimes)<<" ms."<<endl;
}


void Graph::DFSTree(vector<int>& tNodes, int id){
    tNodes.push_back(Tree[id].uniqueVertex);
    for(int i=0;i<Tree[id].ch.size();++i){
        DFSTree(tNodes,Tree[id].ch[i]);
    }
}





/// Index Maintenance


//function of testing the throughput of path-finding system, batchInterval is the time interval between two adjacent update batch (in seconds)
void Graph::IndexMaintenance(int updateType, int updateSize) {
    cout<<"Shortest path query throughput test..."<<endl;
    // read updates
    string file = graphfile + ".update";
    bool ifDebug=false;
//    ifDebug=true;
    vector<pair<pair<int,int>,pair<int,int>>> wBatch;
    int ID1, ID2, oldW, newW;
    srand(0);
    vector<pair<pair<int,int>,int>> updateData;
    ReadUpdate(file, updateData);

    cout<<"Update Number: "<<updateSize<<endl;

    string queryF = graphfile + ".query";
//    filename = graphfile + ".queryParti";
    ifstream IF(queryF);
    if(!IF){
        cout<<"Cannot open file "<<queryF<<endl;
        exit(1);
    }
    cout<<"Query file: "<<queryF<<endl;
    int num;
    vector<pair<int,int>> ODpair;
    IF>>num;
    for(int k=0;k<num;k++){
        IF>>ID1>>ID2;
        ODpair.emplace_back(ID1, ID2);
    }
    IF.close();

    Timer tt;
    double runT1=0, runT2 = 0;
    switch (updateType) {
        case 0:{
            break;
        }
        case 1:{
            //Decrease update
            cout<<"Update type: Decrease"<<endl;
            Graph g2=*this;
            for(int u=0;u<updateSize;u++){
                wBatch.clear();
                ID1 = updateData[u].first.first;
                ID2 = updateData[u].first.second;
                oldW = updateData[u].second;
                newW=oldW*0.5;
                if(newW < 1) {
                    cout<<"New edge weight is not positive! "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                    exit(1);
                }
                wBatch.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
                if(ifDebug){
                    cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                }

                tt.start();
                if(algoChoice==CH){
                    if(PSPStrategy==PreBoundary){
//                        PCHBatchUpdateDecPre(wBatch, u, runT1);
                        g2.PCHBatchUpdateDecPre(wBatch, u, runT1);
                    }else if(PSPStrategy==NoBoundary || PSPStrategy==PostBoundary){
//                        PCHBatchUpdateDec(wBatch, u, runT1);
                        g2.PCHBatchUpdateDec(wBatch, u, runT1);
                    }
                }else if(algoChoice==H2H){
                    if(PSPStrategy==PreBoundary){
//                        PH2HBatchUpdateDecPre(wBatch, u, runT1);
                        g2.PH2HBatchUpdateDecPre(wBatch, u, runT1);
                    }else if(PSPStrategy==NoBoundary || PSPStrategy==PostBoundary){
//                        PH2HBatchUpdateDec(wBatch, u, runT1);
                        g2.PH2HBatchUpdateDec(wBatch, u, runT1);
                    }
                }else if(algoChoice==PLL){
                    if(PSPStrategy==PreBoundary){
//                        PPLLBatchUpdateDecPre(wBatch, u, runT1);
                        g2.PPLLBatchUpdateDecPre(wBatch, u, runT1);
                    }else if(PSPStrategy==NoBoundary || PSPStrategy==PostBoundary){
//                        PPLLBatchUpdateDec(wBatch, u, runT1);
                        g2.PPLLBatchUpdateDec(wBatch, u, runT1);
                    }
                }
                tt.stop();
                runT1+=tt.GetRuntime();
                if(ifDebug){
//                        g2.CorrectnessCheckH2H(100);
                    CorrectnessCheck(100);
                }

            }
            cout<<"Average update time: "<<runT1/updateSize<<" s."<<endl;
//            break;
        }
        case 2:{
            //Increase update
            cout<<"Update type: Increase"<<endl;

            for(int u=0;u<updateSize;++u){
                wBatch.clear();
                ID1 = updateData[u].first.first;
                ID2 = updateData[u].first.second;
                oldW = updateData[u].second;
                newW = oldW * 1.5;
                wBatch.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
                if (ifDebug) {
                    cout << "Batch " << u << ": " << ID1 << " " << ID2 << " " << oldW << " " << newW << endl;
                }
                tt.start();
                if(algoChoice==CH){
                    if(PSPStrategy==PreBoundary){
                        PCHBatchUpdateIncPre(wBatch, u, runT2);
                    }else if(PSPStrategy==NoBoundary || PSPStrategy==PostBoundary){
                        PCHBatchUpdateInc(wBatch, u, runT2);
                    }
                }else if(algoChoice==H2H){
                    if(PSPStrategy==PreBoundary){
                        PH2HBatchUpdateIncPre(wBatch, u, runT2);
                    }else if(PSPStrategy==NoBoundary || PSPStrategy==PostBoundary){
                        PH2HBatchUpdateInc(wBatch, u, runT2);
                    }
                }else if(algoChoice==PLL){
                    if(PSPStrategy==PreBoundary){
                        PPLLBatchUpdateIncPre(wBatch, u, runT2);
                    }else if(PSPStrategy==NoBoundary || PSPStrategy==PostBoundary){
                        PPLLBatchUpdateInc(wBatch, u, runT2);
                    }
                }
                tt.stop();
                runT2+=tt.GetRuntime();
                if(ifDebug){
                    CorrectnessCheck(100);
                }
            }
            cout<<"Average update time: "<<runT2/updateSize<<" s."<<endl;
            break;
        }
        default:{
            cout<<"Wrong update type!"<<endl;
            break;
        }

    }
}

//function for throughput test of decrease update
void Graph::PCHBatchUpdateDec(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1) {
    Timer tt;
    Timer tt2;
    tt.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }

    tt2.start();
    vUpdated.assign(node_num, false);
    vector<vector<pair<pair<int,int>,int>>> updatedSCs;
    updatedSCs.assign(partiNum,vector<pair<pair<int,int>,int>>());
    if(!partiBatch.empty()){
        if(partiBatch.size()>threadnum){
            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
        }
//        cout<<"Update Partition number: "<<partiBatch.size()<<endl;
        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::DecreasePartiBatchUpdateCheckCH, this, pid, boost::ref(it->second), boost::ref(overlayBatch), false, boost::ref(updatedSCs[pid]) ));
        }
        thread.join_all();
    }
    map<pair<int,int>,pair<int,int>> updateSCTrue;
    int updateSCSize=0;
    for(auto it1=partiBatch.begin();it1!=partiBatch.end();++it1){
        int pid=it1->first;
        vector<int> Bid=BoundVertex[pid];
        //check the boundary edge within partition
        int bid1,bid2,olddis,newdis;
        for(auto it=updatedSCs[pid].begin();it!=updatedSCs[pid].end();++it){
            bid1=it->first.first, bid2=it->first.second; newdis=it->second;
            updateSCSize++;
            if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
                olddis=NeighborsOverlay[bid1][bid2];
            }else{//if not found
                cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
            }

            if(newdis<olddis){
//            cout<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
//                NeighborsOverlay[bid1][bid2]=newdis;
//                NeighborsOverlay[bid2][bid1]=newdis;
//                sm->wait();
                if(updateSCTrue.find(make_pair(bid1,bid2))==updateSCTrue.end()){
                    updateSCTrue.insert({make_pair(bid1,bid2), make_pair(olddis,newdis)});
                }else if(updateSCTrue[make_pair(bid1,bid2)].second>newdis){//if found and newdis is smaller
                    cout<<"More than one supportive vertices. "<<bid1<<" "<<bid2<<" "<<updateSCTrue[make_pair(bid1,bid2)].second<<" "<<newdis<<endl;
                    updateSCTrue[make_pair(bid1,bid2)].second=newdis;
                }
//                sm->notify();
            }
        }
    }
//    cout<<"updateSCTrue size: "<<updateSCTrue.size()<<" "<<updateSCSize<<endl;
    for(auto it=updateSCTrue.begin();it!=updateSCTrue.end();++it){
        overlayBatch.emplace_back(it->first,it->second);//weightOverlay collect the changed edges on overlay graph
    }
    if(!overlayBatch.empty()){
//        cout<<"OverlayBatch size: "<<overlayBatch.size()<<endl;
        DecreaseOverlayBatch(overlayBatch,NeighborsOverlay,Tree,rank,heightMax,false);
    }


    // repair the partition index
    if(PSPStrategy>=PostBoundary){
        tt2.start();
        Repair_PartiIndex(true, false, partiBatch);//post
        tt2.stop();

    }

    tt.stop();

}

//function for throughput test of decrease update
void Graph::PCHBatchUpdateDecPre(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1) {
    Timer tt;
    Timer tt2;
    tt.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }

    AllPairBoundaryDisUpdate(true, false, partiBatch, overlayBatch);

    vUpdated.assign(node_num, false);
    if(!partiBatch.empty()){
        if(partiBatch.size()>threadnum){
            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
        }
//        cout<<"Update Partition number: "<<partiBatch.size()<<endl;
        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::DecreasePartiBatchPre, this, pid, boost::ref(it->second), boost::ref(NeighborsParti), boost::ref(Trees[pid]), boost::ref(ranks[pid]), heightMaxs[pid], false));
        }
        thread.join_all();
    }

    if(!overlayBatch.empty()){
//        cout<<"OverlayBatch size: "<<overlayBatch.size()<<endl;
        DecreaseOverlayBatch(overlayBatch,NeighborsOverlay,Tree,rank,heightMax,false);
    }

    tt.stop();
}

//function for throughput test of decrease update
void Graph::PH2HBatchUpdateDec(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1) {
    Timer tt;
    Timer tt2;
    tt.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }

    tt2.start();
    vUpdated.assign(node_num, false);
    vector<vector<pair<pair<int,int>,int>>> updatedSCs;
    updatedSCs.assign(partiNum,vector<pair<pair<int,int>,int>>());
    if(!partiBatch.empty()){
        if(partiBatch.size()>threadnum){
            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
        }
//        cout<<"Update Partition number: "<<partiBatch.size()<<endl;
        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::DecreasePartiBatchUpdateCheckCH, this, pid, boost::ref(it->second), boost::ref(overlayBatch), false, boost::ref(updatedSCs[pid]) ));
        }
        thread.join_all();
    }
    map<pair<int,int>,pair<int,int>> updateSCTrue;
    int updateSCSize=0;
    for(auto it1=partiBatch.begin();it1!=partiBatch.end();++it1){
        int pid=it1->first;
        vector<int> Bid=BoundVertex[pid];
        //check the boundary edge within partition
        int bid1,bid2,olddis,newdis;
        for(auto it=updatedSCs[pid].begin();it!=updatedSCs[pid].end();++it){
            bid1=it->first.first, bid2=it->first.second; newdis=it->second;
            updateSCSize++;
            if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
                olddis=NeighborsOverlay[bid1][bid2];
            }else{//if not found
                cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
            }

            if(newdis<olddis){
//            cout<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
//                NeighborsOverlay[bid1][bid2]=newdis;
//                NeighborsOverlay[bid2][bid1]=newdis;
//                sm->wait();
                if(updateSCTrue.find(make_pair(bid1,bid2))==updateSCTrue.end()){
                    updateSCTrue.insert({make_pair(bid1,bid2), make_pair(olddis,newdis)});
                }else if(updateSCTrue[make_pair(bid1,bid2)].second>newdis){//if found and newdis is smaller
                    cout<<"More than one supportive vertices. "<<bid1<<" "<<bid2<<" "<<updateSCTrue[make_pair(bid1,bid2)].second<<" "<<newdis<<endl;
                    updateSCTrue[make_pair(bid1,bid2)].second=newdis;
                }
//                sm->notify();
            }
        }
    }
//    cout<<"updateSCTrue size: "<<updateSCTrue.size()<<" "<<updateSCSize<<endl;
    for(auto it=updateSCTrue.begin();it!=updateSCTrue.end();++it){
        overlayBatch.emplace_back(it->first,it->second);//weightOverlay collect the changed edges on overlay graph
    }
    if(!overlayBatch.empty()){
//        cout<<"OverlayBatch size: "<<overlayBatch.size()<<endl;
        DecreaseOverlayBatch(overlayBatch,NeighborsOverlay,Tree,rank,heightMax,false);
    }


//    cout<<"algoQuery: PCH-No"<<endl;

    DecreaseOverlayBatchLabel(Tree,rank,heightMax,ProBeginVertexSetOverlay,vertexIDChLOverlay);//vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int> &ProBeginVertexSet, set<int> &vertexIDChL

    if(!partiBatch.empty()){
        if(partiBatch.size()>threadnum){
            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
        }
//        cout<<"Update Partition number: "<<partiBatch.size()<<endl;
        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::DecreasePartiBatchLabel, this, boost::ref(Trees[pid]), boost::ref(ranks[pid]), heightMaxs[pid], boost::ref(ProBeginVertexSetParti[pid]), boost::ref(vertexIDChLParti[pid]) ));
        }
        thread.join_all();
    }

    tt2.stop();
    // repair the partition index
    if(PSPStrategy>=PostBoundary){
        tt2.start();
        Repair_PartiIndex(true, false, partiBatch);//post
        tt2.stop();

    }

    tt.stop();

}

//function for throughput test of decrease update
void Graph::PH2HBatchUpdateDecPre(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1) {
    Timer tt;
    Timer tt2;
    tt.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }

    AllPairBoundaryDisUpdate(true, false, partiBatch, overlayBatch);

    vUpdated.assign(node_num, false);
    if(!partiBatch.empty()){
        if(partiBatch.size()>threadnum){
            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
        }
//        cout<<"Update Partition number: "<<partiBatch.size()<<endl;
        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::DecreasePartiBatchPre, this, pid, boost::ref(it->second), boost::ref(NeighborsParti), boost::ref(Trees[pid]), boost::ref(ranks[pid]), heightMaxs[pid], true));
        }
        thread.join_all();
    }

    if(!overlayBatch.empty()){
//        cout<<"OverlayBatch size: "<<overlayBatch.size()<<endl;
        DecreaseOverlayBatch(overlayBatch,NeighborsOverlay,Tree,rank,heightMax,true);
    }

    tt.stop();
}

//function for throughput test of increase update
void Graph::PCHBatchUpdateInc(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1) {
    Timer tt;
    Timer tt2;
    tt.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }
    tt2.start();
    vUpdated.assign(node_num, false);
    vector<vector<pair<pair<int,int>,int>>> updatedSCs;
    updatedSCs.assign(partiNum,vector<pair<pair<int,int>,int>>());

    if(!partiBatch.empty()){
        if(partiBatch.size()>threadnum){
            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
        }
//        cout<<"Update Partition number: "<<partiBatch.size()<<endl;

        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::IncreasePartiBatchUpdateCheckCH, this, pid, boost::ref(it->second), boost::ref(overlayBatch), false, boost::ref(updatedSCs[pid]) ));
        }
        thread.join_all();

    }
    map<pair<int,int>,pair<int,int>> updateSCTrue;
    int updateSCSize=0;
    for(auto it1=partiBatch.begin();it1!=partiBatch.end();++it1){
        int pid=it1->first;
        vector<int> Bid=BoundVertex[pid];
        //check the boundary edge within partition
        int bid1,bid2,olddis,newdis;

        for(auto it=updatedSCs[pid].begin();it!=updatedSCs[pid].end();++it){
            bid1=it->first.first, bid2=it->first.second; newdis=it->second;
            updateSCSize++;
            if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
                olddis=NeighborsOverlay[bid1][bid2];
            }else{//if not found
                cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
            }
//        cout<<pid<<": "<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
            if(newdis>olddis){//if '=', not problem; if '<', problem
//                sm->wait();
                if(updateSCTrue.find(make_pair(bid1,bid2))==updateSCTrue.end()){
                    updateSCTrue.insert({make_pair(bid1,bid2), make_pair(olddis,newdis)});
                }else if(updateSCTrue[make_pair(bid1,bid2)].second>newdis){//if found and newdis is smaller
                    cout<<"More than one supportive vertices. "<<bid1<<" "<<bid2<<" "<<updateSCTrue[make_pair(bid1,bid2)].second<<" "<<newdis<<endl;
                    updateSCTrue[make_pair(bid1,bid2)].second=newdis;
                }
//                overlayBatch.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));
//                sm->notify();
            } else if(newdis<olddis){
                cout<<"Something wrong happens. "<<bid1<<"("<<PartiTag[bid1].first<<") "<<bid2<<"("<<PartiTag[bid2].first<<") : "<<newdis<<" "<<olddis<< endl;
                exit(1);
            }
        }
    }
//    cout<<"updateSCTrue size: "<<updateSCTrue.size()<<" "<< updateSCSize<<endl;
    for(auto it=updateSCTrue.begin();it!=updateSCTrue.end();++it){
        overlayBatch.emplace_back(it->first,it->second);
    }
    if(!overlayBatch.empty()){
//        cout<<"OverlayBatch size: "<<overlayBatch.size()<<endl;
        IncreaseOverlayBatch(overlayBatch,NeighborsOverlay,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,false);
    }



    // repair the partition index
    if(PSPStrategy>=PostBoundary){
        Repair_PartiIndex(true,true, partiBatch);
//        Repair_PartiIndex(false,true, partiBatch);
    }

    tt.stop();
//    CorrectnessCheck(100);
}

//function for throughput test of increase update
void Graph::PCHBatchUpdateIncPre(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1) {
    Timer tt;
    Timer tt2;
    tt.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }

    vUpdated.assign(node_num, false);
    AllPairBoundaryDisUpdate(true, true, partiBatch, overlayBatch);

    vUpdated.assign(node_num, false);
    if(!partiBatch.empty()){
        if(partiBatch.size()>threadnum){
            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
        }
//        cout<<"Update Partition number: "<<partiBatch.size()<<endl;
        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::IncreasePartiBatchPre, this, pid, boost::ref(it->second), boost::ref(NeighborsParti), boost::ref(Trees[pid]), boost::ref(ranks[pid]), heightMaxs[pid], true));
        }
        thread.join_all();
    }

    if(!overlayBatch.empty()){
//        cout<<"OverlayBatch size: "<<overlayBatch.size()<<endl;
        IncreaseOverlayBatch(overlayBatch,NeighborsOverlay,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,false);
    }
    tt.stop();
//    CorrectnessCheck(100);
}

//function for throughput test of increase update
void Graph::PH2HBatchUpdateInc(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1) {
    Timer tt;
    Timer tt2;
    tt.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }
    tt2.start();
    vUpdated.assign(node_num, false);
    vector<vector<pair<pair<int,int>,int>>> updatedSCs;
    updatedSCs.assign(partiNum,vector<pair<pair<int,int>,int>>());

    if(!partiBatch.empty()){
        if(partiBatch.size()>threadnum){
            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
        }
//        cout<<"Update Partition number: "<<partiBatch.size()<<endl;

        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::IncreasePartiBatchUpdateCheckCH, this, pid, boost::ref(it->second), boost::ref(overlayBatch), false, boost::ref(updatedSCs[pid]) ));
        }
        thread.join_all();

    }
    map<pair<int,int>,pair<int,int>> updateSCTrue;
    int updateSCSize=0;
    for(auto it1=partiBatch.begin();it1!=partiBatch.end();++it1){
        int pid=it1->first;
        vector<int> Bid=BoundVertex[pid];
        //check the boundary edge within partition
        int bid1,bid2,olddis,newdis;

        for(auto it=updatedSCs[pid].begin();it!=updatedSCs[pid].end();++it){
            bid1=it->first.first, bid2=it->first.second; newdis=it->second;
            updateSCSize++;
            if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
                olddis=NeighborsOverlay[bid1][bid2];
            }else{//if not found
                cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
            }
//        cout<<pid<<": "<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
            if(newdis>olddis){//if '=', not problem; if '<', problem
//                sm->wait();
                if(updateSCTrue.find(make_pair(bid1,bid2))==updateSCTrue.end()){
                    updateSCTrue.insert({make_pair(bid1,bid2), make_pair(olddis,newdis)});
                }else if(updateSCTrue[make_pair(bid1,bid2)].second>newdis){//if found and newdis is smaller
                    cout<<"More than one supportive vertices. "<<bid1<<" "<<bid2<<" "<<updateSCTrue[make_pair(bid1,bid2)].second<<" "<<newdis<<endl;
                    updateSCTrue[make_pair(bid1,bid2)].second=newdis;
                }
//                overlayBatch.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));
//                sm->notify();
            } else if(newdis<olddis){
                cout<<"Something wrong happens. "<<bid1<<"("<<PartiTag[bid1].first<<") "<<bid2<<"("<<PartiTag[bid2].first<<") : "<<newdis<<" "<<olddis<< endl;
                exit(1);
            }
        }
    }
//    cout<<"updateSCTrue size: "<<updateSCTrue.size()<<" "<< updateSCSize<<endl;
    for(auto it=updateSCTrue.begin();it!=updateSCTrue.end();++it){
        overlayBatch.emplace_back(it->first,it->second);
    }
    if(!overlayBatch.empty()){
//        cout<<"OverlayBatch size: "<<overlayBatch.size()<<endl;
        IncreaseOverlayBatch(overlayBatch,NeighborsOverlay,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,false);
    }

    IncreaseOverlayBatchLabel(Tree, rank, heightMax, ProBeginVertexSetOverlay, VidtoTNid);

    if(!partiBatch.empty()){
        if(partiBatch.size()>threadnum){
            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
        }
//        cout<<"Update Partition number: "<<partiBatch.size()<<endl;
        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::IncreasePartiBatchLabel, this, boost::ref(Trees[pid]), boost::ref(ranks[pid]), heightMaxs[pid], boost::ref(ProBeginVertexSetParti[pid]), boost::ref(VidtoTNidP) ));//vector<Node> &Tree, vector<int> &rank, int heightMax, int& checknum, vector<int>& ProBeginVertexSet, vector<vector<int>> &VidtoTNid
        }
        thread.join_all();
    }



    // repair the partition index
    if(PSPStrategy>=PostBoundary){
        Repair_PartiIndex(true,true, partiBatch);
//        Repair_PartiIndex(false,true, partiBatch);

    }

    tt.stop();
//    CorrectnessCheck(100);
}

//function for throughput test of increase update
void Graph::PH2HBatchUpdateIncPre(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1) {
    Timer tt;
    Timer tt2;
    tt.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }

    vUpdated.assign(node_num, false);
    AllPairBoundaryDisUpdate(true, true, partiBatch, overlayBatch);

    vUpdated.assign(node_num, false);
    if(!partiBatch.empty()){
        if(partiBatch.size()>threadnum){
            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
        }
//        cout<<"Update Partition number: "<<partiBatch.size()<<endl;
        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::IncreasePartiBatchPre, this, pid, boost::ref(it->second), boost::ref(NeighborsParti), boost::ref(Trees[pid]), boost::ref(ranks[pid]), heightMaxs[pid], true));
        }
        thread.join_all();
    }

    if(!overlayBatch.empty()){
//        cout<<"OverlayBatch size: "<<overlayBatch.size()<<endl;
        IncreaseOverlayBatch(overlayBatch,NeighborsOverlay,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,false);
    }
    tt.stop();
//    CorrectnessCheck(100);
}

void Graph::EachNodeProBDis5H2H(int child,vector<int>& line,set<int>& vertexIDChL, map<int,int>& checkedDis){
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
    }
    else{
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
        EachNodeProBDis5H2H(Tree[child].ch[i], line, vertexIDChL,checkedDis);
    }
    line.pop_back();

}


void Graph::eachNodeProcessIncrease1H2H(int children, vector<int>& line, int& changelabel){
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
//                if(Tree[PID].height>Tree[children].height){///modified for correctness, PID may not be the descendant of children
                if(Tree[PID].height>Tree[children].height && Tree[PID].vAncestor[childH] == childID){///modified for correctness, PID may not be the descendant of children
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
        eachNodeProcessIncrease1H2H(Tree[children].ch[i],line,changelabel);
    }
    line.pop_back();
}

//partition update of PH2H
void Graph::DecreasePartiBatchUpdateCheck(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay){
    //partition batch decrease update
    vector<pair<pair<int,int>,int>> updatedSC;

    DecreasePartiBatch(pid, wBatch, NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid], updatedSC, true);

    vector<int> Bid=BoundVertex[pid];
    //check the boundary edge within partition
    int bid1,bid2,olddis,newdis;
    for(int i=0;i<Bid.size();i++){
        bid1=Bid[i];
        for(int j=i+1;j<Bid.size();j++){
            bid2=Bid[j];
            if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
                olddis=NeighborsOverlay[bid1][bid2];
            }else{//if not found
//                    cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl;
                continue;//exit(1);
            }

            newdis=QueryH2HPartition(bid1,bid2,pid);
            if(newdis<olddis){
//                    cout<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
//                NeighborsOverlay[bid1][bid2]=newdis;
//                NeighborsOverlay[bid2][bid1]=newdis;
                sm->wait();
                weightOverlay.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));//weightOverlay collect the changed edges on overlay graph
                sm->notify();
            }
        }
    }
}
//partition update PCH
void Graph::DecreasePartiBatchUpdateCheckCH(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay, bool ifOpt, vector<pair<pair<int,int>,int>>& updatedSC){
    //partition batch decrease update
//    vector<pair<pair<int,int>,int>> updatedSC;


    if(ifOpt){
        DecreasePartiBatchForOpt(pid, wBatch, NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid], updatedSC, false, false);
    }else{
        DecreasePartiBatch(pid, wBatch, NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid], updatedSC, false);
    }


//    vector<int> Bid=BoundVertex[pid];
//    //check the boundary edge within partition
//    int bid1,bid2,olddis,newdis;
//    for(auto it=updatedSC.begin();it!=updatedSC.end();++it){
//        bid1=it->first.first, bid2=it->first.second; newdis=it->second;
//        if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
//            olddis=NeighborsOverlay[bid1][bid2];
//        }else{//if not found
//            cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
//        }
//
//        if(newdis<olddis){
////            cout<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
////                NeighborsOverlay[bid1][bid2]=newdis;
////                NeighborsOverlay[bid2][bid1]=newdis;
//            sm->wait();
//            weightOverlay.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));//weightOverlay collect the changed edges on overlay graph
//            sm->notify();
//        }
//    }

}


//for PH2H
void Graph::IncreasePartiBatchUpdateCheck(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay){
    //partition batch Increase update
    vector<pair<pair<int,int>,int>> updatedSC;

    IncreasePartiBatch(pid, wBatch,NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid],SCconNodesMTP,VidtoTNidP, updatedSC,true);

    vector<int> Bid=BoundVertex[pid];
    //check the boundary edge within partition
    int bid1,bid2,olddis,newdis;
    for(int i=0;i<Bid.size();i++){
        bid1=Bid[i];
        for(int j=i+1;j<Bid.size();j++){
            bid2=Bid[j];
            if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
                olddis=NeighborsOverlay[bid1][bid2];
            }else{//if not found
//                    cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl;
                continue;//exit(1);
            }

            newdis=QueryH2HPartition(bid1,bid2,pid);
            if(newdis>olddis){//if '=', not problem; if '<', problem
                sm->wait();
                weightOverlay.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));
                sm->notify();
            } else if(newdis<olddis){
                cout<<"Something wrong happens. "<<bid1<<"("<<PartiTag[bid1].first<<") "<<bid2<<"("<<PartiTag[bid2].first<<") : "<<newdis<<" "<<olddis<< endl;
                exit(1);
            }
        }
    }
}
//for PCH
void Graph::IncreasePartiBatchUpdateCheckCH(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay, bool ifOpt, vector<pair<pair<int,int>,int>>& updatedSC){
    //partition batch Increase update
//    vector<pair<pair<int,int>,int>> updatedSC;

    if(ifOpt){
        IncreasePartiBatchForOpt(pid, wBatch,NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid],SCconNodesMTP,VidtoTNidP, updatedSC,false);
    }else{
        IncreasePartiBatch(pid, wBatch,NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid],SCconNodesMTP,VidtoTNidP, updatedSC,false);
    }

//    cout<<pid<<". Size of updatedSC: "<<updatedSC.size()<<endl;
//    vector<int> Bid=BoundVertex[pid];
//    //check the boundary edge within partition
//    int bid1,bid2,olddis,newdis;
//
//    for(auto it=updatedSC.begin();it!=updatedSC.end();++it){
//        bid1=it->first.first, bid2=it->first.second; newdis=it->second;
//        if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
//            olddis=NeighborsOverlay[bid1][bid2];
//        }else{//if not found
//            cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
//        }
////        cout<<pid<<": "<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
//        if(newdis>olddis){//if '=', not problem; if '<', problem
//            sm->wait();
//            weightOverlay.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));
//            sm->notify();
//        } else if(newdis<olddis){
//            cout<<"Something wrong happens. "<<bid1<<"("<<PartiTag[bid1].first<<") "<<bid2<<"("<<PartiTag[bid2].first<<") : "<<newdis<<" "<<olddis<< endl;
//            exit(1);
//        }
//    }

}

