/*
 * CoreTree.cpp
 *
 *  Created on: 16 June 2023
 *      Author: Xinjie ZHOU
 */
#include "head.h"

/// Index Construction
void Graph::IndexConstruction() {
    if(algoTree==0){
        if(strategy!=NoBoundary){
            cout<<"Wrong PSP strategy! "<< strategy<<endl; exit(1);
        }
        CTIndexConstruct(true);
    }else if(algoTree==1){
        CTIndexConstruct(false);
    }
}

void Graph::CTIndexConstruct(bool ifCH){
    double runT1, runT2, runT3, runT4;
    runT1=0, runT2=0, runT3=0, runT4=0;
    Timer tt;
    tt.start();
    ///graph contraction and create trees
    Construct_tree(true);
    tt.stop();
    runT1 = tt.GetRuntime();
    cout<<"Graph contraction time: "<<runT1<<" s"<<endl;


//    g.WriteOrder(graphfile+".order");
//    WriteCoreGraph(graphfile+"C");
//    exit(0);

    ///core index construction
    tt.start();
    double T_core=Construct_core(algoCoreC);
//    WriteCoreIndex(graphfile);
//    ReadCoreIndex(graphfile+"-"+ to_string(bandWidth));
    tt.stop();
    runT2 = tt.GetRuntime();
    cout<<"Core's index construction time (with postprocessing): "<<runT2<<" s."<< endl;


    ///tree index construction
    if(!ifCH){
        tt.start();
        Compute_tree_label(ifParallel);//Construct periphery index (H2H label + interface label)
        tt.stop();
        runT3 = tt.GetRuntime();
        cout<<"Partition's index construction time: "<<runT3<<" s"<<endl;
    }
    cout<<"Overall Tree's construction time: "<<runT1+runT3<<" s."<<endl;

    cout<<"Overall Index Construction Time: "<<runT1+T_core+runT3+runT4<<" s."<<endl;
    cout<<"Overall Time: "<<runT1+runT2+runT3+runT4<<" s."<<endl;

    IndexsizeCTH2H();//index (labeling+pruning point) size computation
}

void Graph::QueryGenerationSameParti(){

    Construct_tree(true);
//#ifdef __APPLE__
////    cout<<"The platform is macOS."<<endl;
//    ODGeneSameParti(10000,graphfile+".querySamePartiCoreTree");//same partition
//#else
//    ODGeneSameParti(10000,"/home/data/xzhouby/datasets/"+dataset+"/"+dataset+".querySamePartiCoreTree");//same partition
//#endif

    ODGeneSameParti(10000,sourcePath+"tmp/"+dataset+"_sameParti_CT"+ to_string(bandWidth)+".query");//same partition

    exit(0);
}

//Function of generating realistic query
void Graph::ODGeneSameParti(int num, string filename){
    set<pair<int,int>> ODpair;
    vector<pair<int,int>> ODpairVec;

    srand (0);
    int s, t;
    cout<<"Generating same-partition queries..."<<endl;

    vector<vector<int>> PartiVertex;
    PartiVertex.assign(partiNum,vector<int>());
    int len=HighestOrder-1;
    int ID,PID;
    set<int> pSet;
    for(int i=0;i<node_num;i++){
        PID=CoreTag[i];
        if(CoreTag[i]!=-1){
            PartiVertex[PID].push_back(i);
            pSet.insert(PID);
        }
    }
    vector<int> pVector;
    for(auto it=pSet.begin();it!=pSet.end();++it){
        PID=*it;
        if(PartiVertex[PID].size()>=2){
            pVector.push_back(PID);
        }
    }
    cout<<"Partition number: "<<pSet.size()<<" ; "<<pVector.size()<<endl;

    for(int i=0;i<num;i++){
        int pid=pVector[rand()%pVector.size()];
        s=PartiVertex[pid][rand()%PartiVertex[pid].size()];
        t=PartiVertex[pid][rand()%PartiVertex[pid].size()];

        while(s==t){
            t=PartiVertex[pid][rand()%PartiVertex[pid].size()];
        }
        if(CoreTag[s]==-1 || CoreTag[t]==-1){
            cout<<"Not periphery vertex. "<<s<<"("<<CoreTag[s]<<") "<<t<<"("<<CoreTag[t]<<")"<<endl; exit(1);
        }
        if(ODpair.find(make_pair(s,t))==ODpair.end()){
            ODpairVec.emplace_back(s,t);
            ODpair.insert(make_pair(s,t));
            ODpair.insert(make_pair(t,s));
        }else{
            i--;
        }
    }

    cout<<"generated OD pair number "<<ODpairVec.size()<<endl;

    ofstream OF(filename);
    if (!OF) { // if not exist
        cout << "Fail to open file" << filename << endl;
        exit(1);
    }
    OF<<ODpairVec.size()<<endl;
    for(int k=0;k<ODpairVec.size();k++){
        OF<<ODpairVec[k].first<<" "<<ODpairVec[k].second<<endl;
    }
    OF.close();
    cout<<"Finish."<<endl;
}

void Graph::IndexMaintenanceTypeTest(int updateType, int updateBatch) {
    cout<<"Same core update"<<endl;
    IndexMaintenanceType(sourcePath+dataset+"sameCore_CT"+ to_string(bandWidth)+".query",updateType,updateBatch,0);
    cout<<"Same tree update"<<endl;
    IndexMaintenanceType(sourcePath+dataset+"sameTree_CT"+ to_string(bandWidth)+".query",updateType,updateBatch,1);
}

void Graph::IndexMaintenanceType(string updateFile, int updateType, int updateBatch, int type){
    cout<<"Update Batch: "<<updateBatch<<endl;
    vector<pair<pair<int,int>,int>> testdata;
    ReadUpdate(updateFile, testdata);
//    vector<pair<pair<int,int>,tuple<int,int,int>>> testdata2;
//    ReadUpdate3(updateFile,testdata2);

    int ID1,ID2, oldW,newW;
    double runT1=0;
    Timer tt1;
    bool ifDebug=false;
//    ifDebug=true;
    auto NeighborTemp = Neighbor; auto AdjaCoreTemp = AdjaCore; auto AdjaCoreMapTemp = AdjaCoreMap;
    auto LabelTemp = Label; auto TreeTemp = Tree; auto LabelVTemp=LabelV.Labels;

    switch (updateType) {
        case 1:{
            cout<<"Update type: Decrease"<<endl;
//            Graph g1=*this;
            runT1 = 0;
            for(int u=0;u<updateBatch;u++){
                ID1=testdata[u].first.first;
                ID2=testdata[u].first.second;
                oldW=testdata[u].second;
                newW=oldW*0.75;
                if(type==0){//same-core
                    if(CoreTag[ID1]!=-1 || CoreTag[ID2]!=-1){
                        cout<<"Wrong Same-Core Update pair! "<<ID1<<"("<<CoreTag[ID1]<<"), "<<ID2<<"("<<CoreTag[ID2]<<")."<<endl; exit(1);
                    }
                }
                else if(type==1){//same-tree
                    if(CoreTag[ID1]==-1 && CoreTag[ID2]==-1){
                        cout<<"Wrong Same-Tree Update pair! "<<ID1<<"("<<CoreTag[ID1]<<"), "<<ID2<<"("<<CoreTag[ID2]<<")."<<endl; exit(1);
                    }
                }
//                ID1=testdata2[u].first.first;
//                ID2=testdata2[u].first.second;
//                oldW=get<0>(testdata2[u].second);
//                newW=oldW-(oldW-get<1>(testdata2[u].second))/5;
//                if(newW==oldW){
//                    newW = get<1>(testdata2[u].second);
//                }
//                cout<<oldW<<" "<<get<1>(testdata2[u].second)<< " "<<newW<<endl;
                if(newW == oldW){
                    newW = oldW-1;
                }
                if(newW < 1) {
                    cout<<"New edge weight is not positive! "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                    exit(1);
                }
                cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW;//<<endl;

                tt1.start();
//                g1.Decrease(ID1,ID2,oldW,newW);
                Decrease(ID1,ID2,oldW,newW);
                tt1.stop();

                runT1 += tt1.GetRuntime();

                cout<<" Time: "<<tt1.GetRuntime()<< " s."<<endl;
                if(ifDebug){
                    CorrectnessCheck(100);
                }

//                testdata[u].second = newW;
            }
            cout<<"Average Decrease update Time: "<<runT1/updateBatch<<" s.\n"<<endl;
//            break;
            Neighbor = NeighborTemp; AdjaCore = AdjaCoreTemp; AdjaCoreMap = AdjaCoreMapTemp;
            Label = LabelTemp; Tree = TreeTemp; LabelV.Labels = LabelVTemp;
        }
        case 2:{
            cout<<"Update type: Increase"<<endl;
            runT1 = 0;
            for(int u=0;u<updateBatch;u++){
//            for(int u=406;u<updateBatch;u++){//core 30
                ID1=testdata[u].first.first;
                ID2=testdata[u].first.second;
                oldW=testdata[u].second;
                newW=oldW*1.25;
                if(type==0){//same-core
                    if(CoreTag[ID1]!=-1 || CoreTag[ID2]!=-1){
                        cout<<"Wrong Same-Core Update pair! "<<ID1<<"("<<CoreTag[ID1]<<"), "<<ID2<<"("<<CoreTag[ID2]<<")."<<endl; exit(1);
                    }
                }
                else if(type==1){//same-tree
                    if(CoreTag[ID1]==-1 && CoreTag[ID2]==-1){
                        cout<<"Wrong Same-Tree Update pair! "<<ID1<<"("<<CoreTag[ID1]<<"), "<<ID2<<"("<<CoreTag[ID2]<<")."<<endl; exit(1);
                    }
                }
//                ID1=testdata2[u].first.first;
//                ID2=testdata2[u].first.second;
//                oldW=get<0>(testdata2[u].second);
//                newW=get<2>(testdata2[u].second);
                if(newW == oldW){
                    newW = oldW+1;
                }
                cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW;//<<endl;

                tt1.start();
                Increase(ID1,ID2,oldW,newW);
                tt1.stop();

                runT1 += tt1.GetRuntime();
                cout<<" Time: "<<tt1.GetRuntime()<< " s."<<endl;

                if(u==43){//core 15
//                    g.WriteCoreGraph(graphfile+"Core153");
//                    g.WriteLabels(graphfile+"Core153");
                }

//                g.DijkstraPath(63055,233384);
//                g.QueryCoreDebug(142488,143850);
                if(ifDebug){
                    CorrectnessCheck(100);
                }

            }
            cout<<"Average Increase update Time: "<<runT1/updateBatch<<" s."<<endl;
            Neighbor = NeighborTemp; AdjaCore = AdjaCoreTemp; AdjaCoreMap = AdjaCoreMapTemp;
            Label = LabelTemp; Tree = TreeTemp; LabelV.Labels = LabelVTemp;
            break;
        }
        default:
            break;
    }
}

void Graph::IndexMaintenance(string updateFile, int updateType, int updateBatch){
    cout<<"Update Batch: "<<updateBatch<<endl;
    vector<pair<pair<int,int>,int>> testdata;
    ReadUpdate(updateFile, testdata);
//    vector<pair<pair<int,int>,tuple<int,int,int>>> testdata2;
//    ReadUpdate3(updateFile,testdata2);

    int ID1,ID2, oldW,newW;
    double runT1=0;
    Timer tt1;
    bool ifDebug=false;
//    ifDebug=true;
    auto NeighborTemp = Neighbor; auto AdjaCoreTemp = AdjaCore; auto AdjaCoreMapTemp = AdjaCoreMap;
    auto LabelTemp = Label; auto TreeTemp = Tree; //auto LabelVTemp=LabelV.Labels;

    switch (updateType) {
        case 1:{
            cout<<"Update type: Decrease"<<endl;
//            Graph g1=*this;
            runT1 = 0;
            for(int u=0;u<updateBatch;u++){
                ID1=testdata[u].first.first;
                ID2=testdata[u].first.second;
                oldW=testdata[u].second;
                newW=oldW*0.5;
//                ID1=testdata2[u].first.first;
//                ID2=testdata2[u].first.second;
//                oldW=get<0>(testdata2[u].second);
//                newW=oldW-(oldW-get<1>(testdata2[u].second))/5;
//                if(newW==oldW){
//                    newW = get<1>(testdata2[u].second);
//                }
//                cout<<oldW<<" "<<get<1>(testdata2[u].second)<< " "<<newW<<endl;
                if(newW == oldW){
                    newW = oldW-1;
                }
                if(newW < 1) {
                    cout<<"New edge weight is not positive! "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                    exit(1);
                }
                if(ifDebug){
                    cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW;//<<endl;
                }


                tt1.start();
//                g1.Decrease(ID1,ID2,oldW,newW);
                Decrease(ID1,ID2,oldW,newW);
                tt1.stop();

                runT1 += tt1.GetRuntime();



                if(ifDebug){
                    cout<<" Time: "<<tt1.GetRuntime()<< " s."<<endl;
                    CorrectnessCheck(100);
                }

//                testdata[u].second = newW;
            }
            cout<<"Average Decrease update Time: "<<runT1/updateBatch<<" s.\n"<<endl;
//            break;
            Neighbor = NeighborTemp; AdjaCore = AdjaCoreTemp; AdjaCoreMap = AdjaCoreMapTemp;
            Label = LabelTemp; Tree = TreeTemp; //LabelV.Labels = LabelVTemp;
//            break;
        }
        case 2:{
            cout<<"Update type: Increase"<<endl;
            runT1 = 0;
            for(int u=0;u<updateBatch;u++){
//            for(int u=406;u<updateBatch;u++){//core 30
                ID1=testdata[u].first.first;
                ID2=testdata[u].first.second;
                oldW=testdata[u].second;
                newW=oldW*1.5;
//                ID1=testdata2[u].first.first;
//                ID2=testdata2[u].first.second;
//                oldW=get<0>(testdata2[u].second);
//                newW=get<2>(testdata2[u].second);
                if(newW == oldW){
                    newW = oldW+1;
                }
                if(ifDebug){
                    cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW;//<<endl;
                }


                tt1.start();
                Increase(ID1,ID2,oldW,newW);
                tt1.stop();

                runT1 += tt1.GetRuntime();

                if(ifDebug){
                    cout<<" Time: "<<tt1.GetRuntime()<< " s."<<endl;
                    CorrectnessCheck(100);
                }

            }
            cout<<"Average Increase update Time: "<<runT1/updateBatch<<" s."<<endl;
            Neighbor = NeighborTemp; AdjaCore = AdjaCoreTemp; AdjaCoreMap = AdjaCoreMapTemp;
            Label = LabelTemp; Tree = TreeTemp; //LabelV.Labels = LabelVTemp;
            break;
        }
        default:
            break;
    }
}
//function for computing the index size
void Graph::IndexsizeCTH2H(){
    unsigned long long m=0,m1=0,m2=0,m3=0,m4=0,m5=0;

    //core index
    for(int k=0;k<Label.size();k++){
        m1+=Label[k].size()*2*sizeof(int);
    }

    for(int i=0;i<PruningPointSet.size();i++){//Order
        for(auto it=PruningPointSet[i].begin();it!=PruningPointSet[i].end();it++){//Order
            m2+=(1+(*it).second.size())*sizeof(int);
        }
    }

    //Periphery index
    if(algoTree==0){
        for(int i=0;i<Tree.size();i++){
            m3+=Tree[i].vert.size()*4*sizeof(int);//ID,dis,cnt,pos
        }
    }else if(algoTree==1){
        for(int i=0;i<Tree.size();i++){
            m3+=Tree[i].vert.size()*4*sizeof(int);//ID,dis,cnt,pos
            m3+=Tree[i].dis.size()*3*sizeof(int);//ID,dis,cnt
            m3+=Tree[i].disInf.size()*2*sizeof(int);//ID,dis
        }
    }


    for(int i=0;i< SCconNodesMT.size();i++){
        for(auto it=SCconNodesMT[i].begin(); it!=SCconNodesMT[i].end(); it++){
            m4+=sizeof(int)+(*it).second.size()*2*sizeof(int);
        }
    }

    //cout<<"Index size "<<(double)m1/1024/1024<<", Pruning point size "<<(double)m2/1024/1024<<endl;
    m=m1+m2+m3+m4+m5;
    cout<<"Core label size: "<<(double)m1/1024/1024<<" MB"<<endl;
    cout<<"Pruning point size: "<<(double)m2/1024/1024<<" MB"<<endl;
    cout<<"Core Index size: "<<(double)(m1+m2)/1024/1024<<" MB"<<endl;
    cout<<"Tree label size: "<<(double)m3/1024/1024<<" MB"<<endl;
    cout<<"Tree update info size: "<<(double)m4/1024/1024<<" MB"<<endl;
    cout<<"Tree Index size: "<<(double)(m3+m4)/1024/1024<<" MB"<<endl;
    cout<<"Minimum index size "<<(double)(m1+m3)/1024/1024<<" MB"<<endl;
    cout<<"Overall index size "<<(double)m/1024/1024<<" MB"<<endl;
}


/// Query Processing
//function for correctness check
void Graph::CorrectnessCheck(int runtimes){
    Timer tt;
    double runT=0;
    srand (time(NULL));
    int s, t, d1, d2, d3;
//    runtimes = 1;
    cout<<"Correctness check ("<<runtimes<<" rounds) ... ";
    for(int i=0;i<runtimes;i++){
//        if(i%100==0) cout<<i<<endl;
        s=rand()%node_num;
        t=rand()%node_num;
//        s=530835,t=802626;//GO
//        s=248272,t=242169;//NY

        if(runtimes == 1){
//            cout<<"s: "<<s<<" ; t: "<<t<<endl;
        }
        d1=Dijkstra(s,t,Neighbor);

        if(algoTree==0){
            tt.start();
            d2=Query_CH(s,t);
            tt.stop();
        }else if(algoTree==1){
            tt.start();
            d2=Query(s,t);
            tt.stop();
        }

        runT+=tt.GetRuntime();
//        cout<<s<<"("<<CoreTag[s]<<") "<<t<<"("<<CoreTag[t]<<") "<<d2<<" "<<d1<<endl;
        if(d1!=d2){
            cout<<"InCorrect! "<<i<<": "<<s<<"("<<NodeOrder[s]<<") "<<t<<"("<<NodeOrder[t]<<") "<<d2<<" "<<d1;
            cout<<" ; CoreTag: "<<CoreTag[s]<<" "<<CoreTag[t]<<endl;
            if(algoTree==1){
                QueryDebug(s,t);
            }

        }
    }
    cout<<"Average Query Time: "<<1000*runT/runtimes<<" ms."<<endl;
}
//function for core index correctness check
void Graph::CorrectnessCheckCore(){
    srand (time(NULL));
    int s, t, d1, d2, d3;
    vector<int> coreVertex;
    for(int i=0;i<node_num;++i){
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
                cout<<"InCorrect! "<<s<<"("<<NodeOrder[s]<<") "<<t<<"("<<NodeOrder[t]<<"): "<<d1<<" "<<d2<<endl;
//				DijkstraPath(s,t);
//				DijkstraCorePath(s,t);
            }
        }else
            i--;
    }
}
//function for core index correctness check, all pairs
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
//function for efficiency test
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
    cout<<"Efficiency test. Run times: "<<runtimes<<endl;
    Timer tt;
    double runT=0;
    int d1,d2;
    bool ifDebug=false;
//    ifDebug=true;

    if(ifDebug){
        cout<<"With correctness check."<<endl;
    }

    if(!LabelV.Labels.empty()){
        ifVectorQ=true;
    }else{
        ifVectorQ= false;
    }

    vector<pair<int,double>> queryTimes;
    queryTimes.assign(4, make_pair(0,0));
    vector<string> queryTypes;
    queryTypes.push_back("Both-core: ");
    queryTypes.push_back("Same-partition: ");
    queryTypes.push_back("Different-partition: ");
    queryTypes.push_back("Core-partition: ");

    if(algoTree==0){
        for(int i=0;i<runtimes;i++){
            ID1=ODpair[i].first,ID2=ODpair[i].second;
            tt.start();
            d1=Query_CH(ID1,ID2);
            tt.stop();
            runT+=tt.GetRuntime();
//            cout<<i<<", "<<ID1<<"("<<CoreTag[ID1]<<") "<<ID2<<"("<<CoreTag[ID2]<<") : "<<tt.GetRuntime()<<endl;
            if(CoreTag[ID1]==-1 && CoreTag[ID2]==-1){//both core
                queryTimes[0].first++;
                queryTimes[0].second+=tt.GetRuntime();
            }else if(CoreTag[ID1]!=-1 && CoreTag[ID1]==CoreTag[ID2]){//same-partition
                queryTimes[1].first++;
                queryTimes[1].second+=tt.GetRuntime();
            }else if(CoreTag[ID1]!=-1 && CoreTag[ID2]!=-1 && CoreTag[ID1]!=CoreTag[ID2]){//different-partition
                queryTimes[2].first++;
                queryTimes[2].second+=tt.GetRuntime();
            }else if((CoreTag[ID1]!=-1 && CoreTag[ID2]==-1) || (CoreTag[ID1]==-1 && CoreTag[ID2]!=-1)){//one core, one partition
                queryTimes[3].first++;
                queryTimes[3].second+=tt.GetRuntime();
            }else{
                cout<<"Wrong scenario. "<<i<<", "<<ID1<<"("<<CoreTag[ID1]<<") "<<ID2<<"("<<CoreTag[ID2]<<") : "<<tt.GetRuntime()<<endl; exit(1);
            }
            if(ifDebug){
                d2= Dijkstra(ID1,ID2,Neighbor);
                if(d1!=d2){
                    cout<<"Incorrect query result! "<<ID1<<" "<<ID2<<" "<<d1<<" "<<d2<<endl; exit(1);
                }
            }
        }
    }else if(algoTree==1){
        for(int i=0;i<runtimes;i++){
            ID1=ODpair[i].first,ID2=ODpair[i].second;
            tt.start();
            d1=Query(ID1,ID2);
            tt.stop();
            runT+=tt.GetRuntime();
//            cout<<i<<", "<<ID1<<"("<<CoreTag[ID1]<<") "<<ID2<<"("<<CoreTag[ID2]<<") : "<<tt.GetRuntime()<<endl;
            if(CoreTag[ID1]==-1 && CoreTag[ID2]==-1){//both core
                queryTimes[0].first++;
                queryTimes[0].second+=tt.GetRuntime();
            }else if(CoreTag[ID1]!=-1 && CoreTag[ID1]==CoreTag[ID2]){//same-partition
                queryTimes[1].first++;
                queryTimes[1].second+=tt.GetRuntime();
            }else if(CoreTag[ID1]!=-1 && CoreTag[ID2]!=-1 && CoreTag[ID1]!=CoreTag[ID2]){//different-partition
                queryTimes[2].first++;
                queryTimes[2].second+=tt.GetRuntime();
            }else if((CoreTag[ID1]!=-1 && CoreTag[ID2]==-1) || (CoreTag[ID1]==-1 && CoreTag[ID2]!=-1)){//one core, one partition
                queryTimes[3].first++;
                queryTimes[3].second+=tt.GetRuntime();
            }else{
                cout<<"Wrong scenario. "<<i<<", "<<ID1<<"("<<CoreTag[ID1]<<") "<<ID2<<"("<<CoreTag[ID2]<<") : "<<tt.GetRuntime()<<endl; exit(1);
            }
            if(ifDebug){
                d2= Dijkstra(ID1,ID2,Neighbor);
                if(d1!=d2){
                    cout<<"Incorrect query result! "<<ID1<<" "<<ID2<<" "<<d1<<" "<<d2<<endl; exit(1);
                }
            }
        }
    }
    cout<<"Query time: "<<endl;
    for(int i=0;i<queryTimes.size();++i){
        cout<<queryTypes[i]<<" ("<<queryTimes[i].first<<"): "<<1000*queryTimes[i].second/queryTimes[i].first<<endl;
    }


    ifVectorQ= false;

    cout<<"Average Query Time: "<<(double)runT*1000/runtimes<<" ms."<<endl;
}

void Graph::EffiCheckQueryTypeTest(int runtimes) {
    EffiCheckQueryType(sourcePath+dataset+"sameCore_CT"+ to_string(bandWidth)+".query", runtimes, 0);
    EffiCheckQueryType(sourcePath+dataset+"coreTree_CT"+ to_string(bandWidth)+".query", runtimes, 1);
    EffiCheckQueryType(sourcePath+dataset+"crossTree_CT"+ to_string(bandWidth)+".query", runtimes, 2);
    EffiCheckQueryType(sourcePath+dataset+"sameTree_CT"+ to_string(bandWidth)+".query", runtimes, 3);
}

//function for efficiency test
void Graph::EffiCheckQueryType(string filename, int runtimes, int type) {
    ifstream IF(filename);
    if(!IF){
        cout<<"Cannot open Map "<<filename<<endl;
        exit(1);
    }
    bool ifDebug=false;
//    ifDebug=true;

    int num, ID1, ID2;
    vector<pair<int,int>> ODpair;
    IF>>num;
    for(int k=0;k<num;k++){
        IF>>ID1>>ID2;
        if(type==0){// same-core
            if(CoreTag[ID1]==-1 && CoreTag[ID2]==-1){
                ODpair.push_back(make_pair(ID1, ID2));
            }else{
                cout<<"Wrong Same-Core OD pair! "<<ID1<<"("<<CoreTag[ID1]<<"), "<<ID2<<"("<<CoreTag[ID2]<<")."<<endl; exit(1);
            }
        }else if(type==1){// core-tree
            if(CoreTag[ID1]==-1 && CoreTag[ID2]!=-1){
                ODpair.push_back(make_pair(ID1, ID2));
            }else{
                cout<<"Wrong Core-Tree OD pair! "<<ID1<<"("<<CoreTag[ID1]<<"), "<<ID2<<"("<<CoreTag[ID2]<<")."<<endl; exit(1);
            }
        }else if(type==2){// tree-tree
            if(CoreTag[ID1] != CoreTag[ID2] && CoreTag[ID1]!=-1 && CoreTag[ID2]!=-1){
                ODpair.push_back(make_pair(ID1, ID2));
            }else{
                cout<<"Wrong Tree-Tree OD pair! "<<ID1<<"("<<CoreTag[ID1]<<"), "<<ID2<<"("<<CoreTag[ID2]<<")."<<endl; exit(1);
            }
        }else if(type==3){// same-tree
            if(CoreTag[ID1]==CoreTag[ID2] && CoreTag[ID1]!=-1){
                ODpair.push_back(make_pair(ID1, ID2));
            }else{
                cout<<"Wrong Same-Tree OD pair! "<<ID1<<"("<<CoreTag[ID1]<<"), "<<ID2<<"("<<CoreTag[ID2]<<")."<<endl; exit(1);
            }
        }
    }
    if(runtimes > num){
        runtimes = num;
    }

    cout<<"Efficiency test. "<<filename;
    cout<<"; Run times: "<<runtimes<<endl;
    int s, t;
//    std::chrono::high_resolution_clock::time_point t1, t2;
//    std::chrono::duration<double> time_span;
    Timer tt;
    double runT=0;

    if(ifDebug){
        cout<<"With correctness check."<<endl;
    }

    vector<int> LSize1, LSize2;
    LSize1.assign(runtimes,0);
    LSize2.assign(runtimes,0);
    double lsize1=0,lsize2=0;
    int d1,d2;

    ifVectorQ=true;
    for(int i=0;i<runtimes;i++){
        ID1=ODpair[i].first, ID2=ODpair[i].second;
        tt.start();
        d1=Query(ID1, ID2);
        tt.stop();
        if(ifDebug){
            d2= Dijkstra(ID1,ID2, Neighbor);
            if(d1!=d2){
                cout<<"Incorrect query result! "<<ID1<<" "<<ID2<<" "<<d1<<" "<<d2<<endl; exit(1);
            }
        }
        LSize1[i]=Label[ID1].size();
        LSize2[i]=Label[ID2].size();
        lsize1+=Label[ID1].size(); lsize2+=Label[ID2].size();
        runT += tt.GetRuntime();
    }
    ifVectorQ=false;
    cout<<"Average label size of ID1: "<<lsize1/runtimes<<" ; label size of ID2: "<<lsize2/runtimes<<endl;
    cout<<"Average Query Time: "<<runT*1000/runtimes<<" ms."<<endl;
}

void Graph::DFSTree(vector<int>& tNodes, int id){
    tNodes.push_back(Tree[id].uniqueVertex);
    for(int i=0;i<Tree[id].ch.size();++i){
        DFSTree(tNodes,Tree[id].ch[i]);
    }
}

void Graph::SameTreeUpdateGen(vector<int>& tRoots, int times){
    string STQuery = sourcePath + dataset + ".updateST";

    /*---OD pairs generation---*/
    int pairs = 0;
    int node_start, node_end, tid;
    int temp = 0;
    ofstream outFile(STQuery, ios::out);
    if (!outFile) {
        cout << "File opening failed." << endl;
        assert(outFile);
    }
    cout << "Same-tree update OD pairs file generating..." << endl;
    outFile << times << endl;
    //generate random OD pairs
    pairs = 0;
    int tSize=tRoots.size();
    map<int,vector<int>> tVertex; tVertex.clear();
    set<pair<int,int>> edgeIdSet;
    edgeIdSet.clear();

    while (pairs < times) {
        flag1:
        tid = tRoots[rand() % tSize];
        if(tVertex.find(tid) == tVertex.end()){//if not found
            vector<int> tNodes;
            DFSTree(tNodes,tid);
            tVertex.insert({tid,tNodes});
        }
        vector<int> tNodes = tVertex[tid];
        node_start = tNodes[rand() % tNodes.size()];
        int pos = rand()%Neighbor[node_start].size();
        node_end = Neighbor[node_start][pos].first;
        int ti=0;
        while(CoreTag[node_end]==-1){//if neighbor is core
            pos = rand()%Neighbor[node_start].size();
            node_end = Neighbor[node_start][pos].first;
            ti++;
            if(ti>=10){
                goto flag1;
            }
        }
        pair<int,int> OD;
        assert(node_start!=node_end);
        if(node_start<node_end){
            OD.first=node_start; OD.second=node_end;
        }else{
            OD.first=node_end; OD.second=node_start;
        }
        double prop = 0.01*(rand()%100+1);//[0.01,1]
        int oldW = Neighbor[node_start][pos].second;
        int newW1 = oldW * (1.0 - 0.5*prop);//[0.5,0.995]
        int newW2 = oldW * (1.0 + 0.5*prop);//[1.005,1.5]

        if(newW1 >= 2 && edgeIdSet.find(OD)==edgeIdSet.end()){//if edge weight is no smaller than 2, and it has not been added.
            edgeIdSet.insert(OD);
            outFile << node_start << ' ' << node_end << ' '<<oldW<<' '<<newW1<<' '<<newW2<< endl;
        }else{
            continue;
        }
        ++pairs;
    }

    outFile.close();
    cout << "Finished." << endl;
}

void Graph::SameTreeQueryGen(vector<int>& tRoots, int times){
    string STQuery = sourcePath + "tmp/"+ dataset + "sameTree_CT"+ to_string(bandWidth)+".query";

    /*---OD pairs generation---*/
    int pairs = 0;
    int node_start, node_end, tid;
    int temp = 0;
    ofstream outFile(STQuery, ios::out);
    if (!outFile) {
        cout << "File opening failed." << endl;
        assert(outFile);
    }
    cout << "Same-tree query OD pairs file generating..." << endl;
    outFile << times << endl;
    //generate random OD pairs
    pairs = 0;
    int tSize=tRoots.size();
    map<int,vector<int>> tVertex; tVertex.clear();

    while (pairs < times) {
        tid = tRoots[rand() % tSize];
        if(tVertex.find(tid) == tVertex.end()){//if not found
            vector<int> tNodes;
            DFSTree(tNodes,tid);
            tVertex.insert({tid,tNodes});
        }
        vector<int> tNodes = tVertex[tid];
        node_start = tNodes[rand() % tNodes.size()];
        node_end = tNodes[rand() % tNodes.size()];
        while(node_end == node_start){
            node_end = tNodes[rand() % tNodes.size()];
        }
//        outFile << node_start+1 << ' ' << node_end+1 << endl;
        outFile << node_start << ' ' << node_end << endl;
        ++pairs;
    }

    outFile.close();
    cout << "Finished." << endl;
}

void Graph::SameTreeQueryTest(string filename,int runtimes){
    filename = filename + "ST";
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
    cout<<"Same-tree query test. Run times: "<<runtimes<<endl;
    int s, t;
    std::chrono::high_resolution_clock::time_point t1, t2;
    std::chrono::duration<double> time_span;
    double runT;

    t1=std::chrono::high_resolution_clock::now();
    for(int i=0;i<runtimes;i++){
        if(CoreTag[ID1] == CoreTag[ID2] && CoreTag[ID1] != -1){
            Query(ODpair[i].first,ODpair[i].second);
        }else{
            cout<<"Wrong query type: "<<ODpair[i].first<<" "<<ODpair[i].second<<" "<<CoreTag[ODpair[i].first]<<" "<<CoreTag[ODpair[i].second]<<endl;

        }

    }
    t2=std::chrono::high_resolution_clock::now();

    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
    runT= time_span.count();
    cout<<"Average Query Time: "<<(double)runT*1000/runtimes<<" ms."<<endl;
}

//function for core graph debugging
void Graph::CoreGraphDebug(string graphfile) {
    cout<<"!!!!!!!!!!!!!!!!!!! Core graph test !!!!!!!!!!!!!!!!!!!"<<endl;
    std::chrono::high_resolution_clock::time_point t10;
    std::chrono::high_resolution_clock::time_point t11;
    std::chrono::duration<double> time_span1;
    double runT1=0;
    int s,t,d1,d2;
    int ID1,ID2,oldW,newW;

    ReadCoreGraph(graphfile);
//    set<int> setB;
//    DFS_CC(AdjaCoreMap,CoreVertex,setB,node_num);

    vector<int> coreVertex;
    for(auto it=CoreVertex.begin();it!=CoreVertex.end();++it){
        coreVertex.push_back(*it);
    }
    int num=coreVertex.size();

//    CompareLabel(graphfile+"BPLL.label", graphfile+"BPLLNew.label");
//    exit(0);

//    DijkstraPathCore(190161,181451);

    Timer tt1;
    tt1.start();

//    ReadLabels2(graphfile+"PLL",0);//read PLL
//    vSm.reserve(node_num);
//    for(int i = 0; i < node_num; i++)
//    {
//        Semaphore* s = new Semaphore(1);
//        vSm.push_back(s);
//    }
//    PruningPointBuild(ifParallel);

//    ReadLabels(graphfile+"PLL",0);//read PLL
//    ReadLabels(graphfile+"PSL",0);//read PSL
//    ReadLabels(graphfile,0);//read PLL
//    Construct_core(0);//PLL
//    Construct_core(1);//PSL
//    Construct_core(2);//LCC
//    Construct_core(3);//PCL
    Construct_core(4);//BPCL

//    Construct_core(2);//Batch PSL


    if(dataset == "Test"){
//        Construct_core(0);// PLL index construction
//    Construct_core(1);// PSL index construction
//    WriteLabels(graphfile+".labelPLL");
    }else{
//        ReadLabels(graphfile);
    }

    tt1.stop();
    cout<<"Time for index construction: "<<tt1.GetRuntime()<<" s."<<endl;

//    cout<<QueryCore(0,3)<<" "<<DijkstraCore(0,3)<<endl;
//    cout<<QueryCore(1,3)<<" "<<DijkstraCore(1,3)<<endl;
//    CleanLabel(Label, graphfile);

    CorrectnessCheckCore();

//    WriteLabels(graphfile+"PCL");
//    WriteLabels(graphfile+"LCC");
//    WriteLabels(graphfile+"BPLL-clean");
//    WriteLabels(graphfile+"PLL");
//    WriteLabels(graphfile+"PSL");
//    exit(0);


    //update
    string updateFile2 = graphfile+".update";
    vector<pair<pair<int,int>,pair<int,int>>> testdata2;
    ReadUpdate2(updateFile2, testdata2);
    cout<<"Batch number: "<<testdata2.size()<<endl;
    int BatchNum = 1;
//    for(int u=0;u<testdata2.size();++u){
    for(int u=0;u<BatchNum;++u){
        ID1=testdata2[u].first.first;
        ID2=testdata2[u].first.second;
        oldW=testdata2[u].second.first;
        newW=testdata2[u].second.second;
        cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;

//        cout<<"Before update: "<<QueryCore(164496,207450)<<" "<<DijkstraCore(164496,207450)<<endl;
//        cout<<QueryCore(207450,211814)<<" "<<DijkstraCore(207450,211814)<<endl;
//        if(Label[207450].find(211814) != Label[207450].end()){//if found
//            cout<<Label[207450][211814]<<" "<<QueryCore(207450,211814)<<" "<<DijkstraCore(207450,211814)<<endl;
//        }

        int tHub = 207559;//core 8
        tHub = 230142;

        int id1 = 147520, id2 = 207559;
        id1=54337,id2=226281;

//        if(u==0){
//            cout << "Before update" << endl;
//            for(auto it=Label[id1].begin();it!=Label[id1].end();++it){
//                if(it->first == tHub){
//                    cout<<tHub<<" is the hub of "<<id1<<", dis is "<<it->second<<"("<<DisQueryVallyDebug(id1,tHub,AdjaCore,Label)<<" "<<DijkstraCore(id1,tHub)<<")"<<endl;
//                }
//            }
//            for(auto it=Label[id2].begin();it!=Label[id2].end();++it){
//                if(it->first == tHub){
//                    cout<<tHub<<" is the hub of "<<id2<<", dis is "<<it->second<<"("<<DisQueryVallyDebug(id2,tHub,AdjaCore,Label)<<" "<<DijkstraCore(id2,tHub)<<")"<<endl;
//                }
//            }
////            cout<<DijkstraCore(id1,tHub)<<" "<<DisQueryPeak(id1,tHub,Label)<<" "<<DisQueryVallyDebug(id1,tHub,AdjaCore)<<endl;
////            cout<<DijkstraCore(id2,tHub)<<" "<<DisQueryPeak(id2,tHub,Label)<<" "<<DisQueryVallyDebug(id2,tHub,AdjaCore)<<endl;
////            cout<<DijkstraCore(id1,id2)<<" "<<DisQueryPeak(id1,id2,Label)<<" "<<DisQueryVallyDebug(id1,id2,AdjaCore)<<endl;
//            if(PruningPointSet2[id1].find(tHub) != PruningPointSet2[id1].end()){
//                cout<<"Pruning point "<<id1<<" "<<tHub<<" "<<PruningPointSet2[id1][tHub]<<endl;
//            }
//            if(PruningPointSet2[id2].find(tHub) != PruningPointSet2[id2].end()){
//                cout<<"Pruning point "<<id2<<" "<<tHub<<" "<<PruningPointSet2[id2][tHub]<<endl;
//            }
//            QueryCoreDebug(id1, id2);
//        }

        if(u==1){// || u==279
//            WriteCoreGraph(graphfile+"-1");
//            WriteLabels(graphfile+"-1");
//            cout<<"Done."<<endl;
//            QueryCoreDebug(142488,143850);
        }

//        DijkstraPathCore(212434,99110);
//        QueryCoreDebug(164496,207450);
        t10=std::chrono::high_resolution_clock::now();
//        IncreasePSL(ID1,ID2,oldW,newW,AdjaCore,Label,PruningPointNew,NoSupportedPair);//original vector version with NoSuportedPair
//        IncreasePSLNew(ID1,ID2,oldW,newW,AdjaCore,Label,PruningPointSetOrder);//set version with NoSuportedPair
        IncreasePSLNew(ID1,ID2,oldW,newW,AdjaCore,Label,PruningPointSet);//set version with NoSuportedPair
        t11=std::chrono::high_resolution_clock::now();

        time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
        runT1 += time_span1.count();

//        cout << "After update" << endl;
//        for(auto it=Label[id1].begin();it!=Label[id1].end();++it){
//            if(it->first == tHub){
//                cout<<tHub<<" is the hub of "<<id1<<", dis is "<<it->second<<"("<<DisQueryVallyDebug(id1,tHub,AdjaCore,Label)<<" "<<DijkstraCore(id1,tHub)<<")"<<endl;
//            }
//        }
//        for(auto it=Label[id2].begin();it!=Label[id2].end();++it){
//            if(it->first == tHub){
//                cout<<tHub<<" is the hub of "<<id2<<", dis is "<<it->second<<"("<<DisQueryVallyDebug(id2,tHub,AdjaCore,Label)<<" "<<DijkstraCore(id2,tHub)<<")"<<endl;
//            }
//        }
////        cout<<DijkstraCore(id1,tHub)<<" "<<DisQueryPeak(id1,tHub,Label)<<" "<<DisQueryVallyDebug(id1,tHub,AdjaCore)<<endl;
////        cout<<DijkstraCore(id2,tHub)<<" "<<DisQueryPeak(id2,tHub,Label)<<" "<<DisQueryVallyDebug(id2,tHub,AdjaCore)<<endl;
////        cout<<DijkstraCore(id1,id2)<<" "<<DisQueryPeak(id1,id2,Label)<<" "<<DisQueryVallyDebug(id1,id2,AdjaCore)<<endl;
//        if(PruningPointSet2[id1].find(tHub) != PruningPointSet2[id1].end()){
//            cout<<"Pruning point "<<id1<<" "<<tHub<<" "<<PruningPointSet2[id1][tHub]<<endl;
//        }
//        if(PruningPointSet2[id2].find(tHub) != PruningPointSet2[id2].end()){
//            cout<<"Pruning point "<<id2<<" "<<tHub<<" "<<PruningPointSet2[id2][tHub]<<endl;
//        }
//        QueryCoreDebug(id1, id2);



//        if(dataset == "Test"){
        if(false){
            CorrectnessCheckCoreAll();
        }
        else{
//            CorrectnessCheckCoreAll();
            cout<<"Correctness checking..."<<endl;
            for(int i=0;i<1;i++){
                s=coreVertex[rand()%num];
                t=coreVertex[rand()%num];
                s=151036,t=195578;//core 23
                s=142294,t=141603;//core 23-set int
                s=150990,t=196782;//core 24
                s=150990,t=144591;
//                s=150990,t=144575;
//                s=214227,t=140253;
                s=41076,t=72885;//core 25
                s=182561,t=79284;//core 26
                s=142868,t=79284;//core 26
//                s=142816,t=79284;
//                s=147051,t=196720;//core 27
//                s=147084,t=196720;//core 27
                s=211099,t=208714;//core 28
//                s=211099,t=211337;//core 28
//                s=203827,t=207273;//core 29-1
//                s=45828, t=63068;
                s=71976, t=230938;//NY, original
//                s=71972, t=230938;//NY, original
                s=223905, t=226281;//NY, original2
                s=54337, t=226281;
                s=1, t=2;

                assert(CoreVertex.find(s)!=CoreVertex.end());
                assert(CoreVertex.find(t)!=CoreVertex.end());
                d1=DijkstraCore(s,t);
                d2=QueryCore(s,t);

                cout<<s<<"("<<CoreTag[s]<<") "<<t<<"("<<CoreTag[t]<<") "<<d2<<" "<<d1<<endl;
                if(d1!=d2){
                    int d3=DijkstraCore(t,s);
                    cout<<"InCorrect! "<<i<<": "<<s<<"("<<NodeOrder[s]<<") "<<t<<"("<<NodeOrder[t]<<") : "<<d2<<" "<<d1<<" "<<d3<<endl;
                    QueryDebug(s,t);
                }
            }
        }

    }

//    WriteLabels(graphfile+"PLL2");
//    WriteLabels(graphfile+"PSL2");

    cout<<"Average time for update: "<<runT1/BatchNum<<" s."<<endl;
    exit(0);
}
//function for Query processing
int Graph::Query_CH(int ID1, int ID2){
    int dis=INF;

    if(CoreTag[ID1]==-1 && CoreTag[ID2]==-1){//Case 1: both in core
        //cout<<"Core-Core"<<endl;
        dis=QueryCore(ID1, ID2);
    }else if(CoreTag[ID1]==-1 && CoreTag[ID2]!=-1){//Case 2: ID2 in partition, ID1 in core
        //cout<<"Core-Parti"<<endl;
        dis=QueryPartiCore_CH(ID2,ID1);
    }else if(CoreTag[ID1]!=-1 && CoreTag[ID2]==-1){//Case 2: ID1 in partition, ID2 in core
        //cout<<"Parti-Core"<<endl;
        dis=QueryPartiCore_CH(ID1,ID2);
    }else if(CoreTag[ID1]!=-1 && CoreTag[ID2]!=-1){//both in partition
        if(CoreTag[ID1] != CoreTag[ID2]){//Case 3: in different peripheries
            dis= QueryPartiParti_CH(ID1,ID2);
        }else{//Case 4: in the same periphery
            dis= QuerySameParti_CH(ID1,ID2);
        }
    }
    return dis;
}

//function for Query processing
int Graph::Query(int ID1, int ID2){
    int dis=INF;

    if(CoreTag[ID1]==-1 && CoreTag[ID2]==-1){//Case 1: both in core
//        cout<<"Core-Core"<<endl;
        dis=QueryCore(ID1, ID2);
//        dis=LabelV.query(ID1, ID2);
    }else if(CoreTag[ID1]==-1 && CoreTag[ID2]!=-1){//Case 2: ID2 in partition, ID1 in core
//        cout<<"Core-Parti"<<endl;
        dis=QueryPartiCore(ID2, ID1);

    }else if(CoreTag[ID1]!=-1 && CoreTag[ID2]==-1){//Case 2: ID1 in partition, ID2 in core
//        cout<<"Parti-Core"<<endl;
        dis=QueryPartiCore(ID1, ID2);
    }else if(CoreTag[ID1]!=-1 && CoreTag[ID2]!=-1){//both in partition
        if(CoreTag[ID1] != CoreTag[ID2]){//Case 3: in different peripheries
            dis= QueryPartiParti(ID1,ID2);
        }else{//Case 4: in the same periphery
//                cout<<"Same partition!"<<endl;
            if(strategy==PostBoundary){
                dis= QuerySameParti2(ID1,ID2);
            }else{
                dis= QuerySameParti(ID1,ID2);
            }

        }
    }
    return dis;
}
//Case 1: both in core
int Graph::QueryCore(int ID1, int ID2){
    int d=INF;
    if(ifVectorQ){
        d=LabelV.query(ID1,ID2);
    }
    else{
        int hubfinal,dis1final,dis2final;
        int hub, dis1, dis2;
        if(CoreTag[ID1]!=-1 ||CoreTag[ID2]!=-1){
            cout<<"Core Tag wrong! "<<CoreTag[ID1]<< " "<<CoreTag[ID2]<<endl;
            exit(1);
        }
//    cout<<NodeOrder[ID1]<<" "<<NodeOrder[ID2]<<endl;
        for(auto it=Label[ID1].begin();it!=Label[ID1].end();it++){
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

//	cout<<"hub "<<hubfinal<<",dis "<<d<<endl;
//	cout<<"check labeling "<<dis1final<<" "<<DijkstraCore(ID1,hubfinal)<<" "<<dis2final<<" "<<DijkstraCore(ID2,hubfinal)<<endl;
    }


    return d;
}
//Case 2: one core, one tree
int Graph::QueryPartiCore(int ID1, int ID2){//ID1 partition, ID2 core
    int d=INF;

    int pid=CoreTag[ID1];
    int bid;
    int dis1,dis2;

    for(auto it=Tree[rank[ID1]].disInf.begin();it!=Tree[rank[ID1]].disInf.end();++it){
        bid=it->first;
        dis1=it->second;
        dis2= QueryCore(bid,ID2);
//        dis2= LabelV.query(bid,ID2);
        if(d>dis1+dis2)
            d=dis1+dis2;
    }
    return d;
}
//Case 3: Different trees
int Graph::QueryPartiParti(int ID1, int ID2){//both are within the partition
    int d=INF;

    int pid1=CoreTag[ID1];
    int pid2=CoreTag[ID2];
    if(pid1==pid2){//if in the same partition
        cout<<"Wrong for partition-partition query!"<<endl;
        exit(1);

    }else{//if in different partitions
//        cout<<"Parti-Parti"<<endl;
        vector<int> B1=BoundVertex[pid1];
        vector<int> B2=BoundVertex[pid2];

        map<int,int> m1,m2;
        m1.clear();
        m2.clear();
        int bID1, bID2, tempdis;
        int b1,b2,d1,d2;
        for(int i=0;i<B1.size()-1;i++){
            bID1=B1[i];
            assert(Tree[rank[ID1]].disInf.find(bID1)!=Tree[rank[ID1]].disInf.end());
//            m1.insert(make_pair(bID1,Tree[rank[ID1]].disInfV[i]));
            m1.insert(make_pair(bID1,Tree[rank[ID1]].disInf[bID1]));
        }
        for(int j=0;j<B2.size()-1;j++){
            bID2=B2[j];
            assert(Tree[rank[ID2]].disInf.find(bID2)!=Tree[rank[ID2]].disInf.end());
//            m2.insert(make_pair(bID2,Tree[rank[ID2]].disInfV[j]));
            m2.insert(make_pair(bID2,Tree[rank[ID2]].disInf[bID2]));
        }

        for(int k=0;k<B1.size()-1;k++){
            bID1=B1[k];

            if(m1[bID1]>d)
                continue;

            for(int z=0;z<B2.size()-1;z++){
                bID2=B2[z];

                if(m2[bID2]>d)
                    continue;

                tempdis=m1[bID1]+QueryCore(bID1,bID2)+m2[bID2];
//                tempdis=m1[bID1]+LabelV.query(bID1,bID2)+m2[bID2];
                if(tempdis<d){
                    d=tempdis;
                    d1=m1[bID1]; d2=m2[bID2];
                    b1=bID1; b2=bID2;
                }

            }
        }

//        cout<<"b1, b2, d1, d2: "<<b1<<" "<<b2<<" "<<d1<<" "<<d2<<endl;
    }

    return d;
}
//Case 4: Same tree, for original version
int Graph::QuerySameParti(int ID1, int ID2){//both are within the partition
    int d=INF;

    int pid1=CoreTag[ID1];
    int pid2=CoreTag[ID2];
    if(pid1==pid2){//if in the same partition
//        cout<<"Same-Parti"<<endl;
        int temp_dis = QueryPeripheryTree(ID1,ID2,pid1);/// d2 may be wrong sometimes
        if(temp_dis<d)//QueryH2HPartition(ID1,ID2,pid1)
            d=temp_dis;//QueryH2HPartition(ID1,ID2,pid1);
        vector<int> B=BoundVertex[pid1];
        map<int,int> m1,m2;
        m1.clear();
        m2.clear();
        vector<int> B1,B2;
        B1.clear();
        B2.clear();
        int bID,d1,d2;
        for(int i=0;i<B.size()-1;i++){
            bID=B[i];
            assert(Tree[rank[ID1]].disInf.find(bID)!=Tree[rank[ID1]].disInf.end());
            assert(Tree[rank[ID2]].disInf.find(bID)!=Tree[rank[ID2]].disInf.end());
//            d1=Tree[rank[ID1]].disInf[i];
//            d2=Tree[rank[ID2]].disInf[i];
            d1=Tree[rank[ID1]].disInf[bID];
            d2=Tree[rank[ID2]].disInf[bID];

            if(d1<d){
                B1.push_back(bID);
                m1.insert(make_pair(bID,d1));
            }
            if(d2<d){
                B2.push_back(bID);
                m2.insert(make_pair(bID,d2));
            }
        }

        int bID1, bID2, tempdis;
        if(!B1.empty() && !B2.empty()){
            for(int k=0;k<B1.size();k++){
                bID1=B1[k];
                if(m1[bID1]>d)
                    continue;
                for(int z=0;z<B2.size();z++){
                    bID2=B2[z];
                    if(m2[bID2]>d)
                        continue;
                    tempdis=m1[bID1]+QueryCore(bID1,bID2)+m2[bID2];
//                    tempdis=m1[bID1]+LabelV.query(bID1,bID2)+m2[bID2];
                    if(tempdis<d)
                        d=tempdis;
                }
            }
        }

    }else{//if in different partitions
        cout<<"Wrong for same partition query!"<<endl;
        exit(1);
    }

    return d;
}
//Case 4: Same tree, for query-orient version
int Graph::QuerySameParti2(int ID1, int ID2){//both are within the partition
    int d=INF;

    int pid1=CoreTag[ID1];
    int pid2=CoreTag[ID2];
    if(pid1==pid2){//if in the same partition
//        cout<<"Same-Parti"<<endl;
        d = QueryPeripheryTree2(ID1,ID2,pid1);

    }else{//if in different partitions
        cout<<"Wrong for same partition query!"<<endl;
        exit(1);
    }

    return d;
}

//CH query for OD within one periphery, single source
void Graph::QueryPeriphery_CH(int ID1, set<int> todo, map<int,int> & results){
    //priority queue
    benchmark::heap<2,int,int> pqueue(node_num);

    //closed or not
    vector<bool> vVisited(node_num, false);
    //the existing shortest distance
    vector<int>	vDistance(node_num, INF);

    vDistance[ID1] = 0;
    pqueue.update(ID1,0);
    int topID,topDis,tempID,tempDis,weight;

    while( ! todo.empty() && ! pqueue.empty() ){
        pqueue.extract_min(topID, topDis);

        vVisited[topID] = true;

        if(todo.find(topID) != todo.end()){//if found
            todo.erase(topID);
            results[topID] = topDis;
        }

//        for(auto out=NeighborCon[topID].begin();out!=NeighborCon[topID].end();++out){
        for(auto out=Tree[rank[topID]].vert.begin();out!=Tree[rank[topID]].vert.end();++out){
            tempID = (*out).first;
            weight = (*out).second.first;

            if(!vVisited[tempID]){
                tempDis = vDistance[topID] + weight;
                if(vDistance[tempID] > tempDis){
                    vDistance[tempID] = tempDis;
                    pqueue.update(tempID, tempDis);
                }
            }
        }
    }

}

//CH query for OD within one periphery
int Graph::QueryPeriphery_CH(int ID1, int ID2){
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

            if(vDistanceForward[topNodeIDForward] > d)
                bF = true;

            vVisitedF[topNodeIDForward] = true;

            if(vVisitedB[topNodeIDForward]){
                int distTmp=topDisForward+vDistanceBackward[topNodeIDForward];
                if(distTmp<d){
                    d=distTmp;
                }
            }

//            for(auto out=NeighborCon[topNodeIDForward].begin();out!=NeighborCon[topNodeIDForward].end();out++){
            for(auto out=Tree[rank[topNodeIDForward]].vert.begin();out!=Tree[rank[topNodeIDForward]].vert.end();out++){
                neighborNodeID = (*out).first;
                neighborLength = (*out).second.first;

                int df = vDistanceForward[topNodeIDForward] + neighborLength;
                if(!vVisitedF[neighborNodeID]){
                    if(vDistanceForward[neighborNodeID] > df){
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
                }
            }

//            for(auto in=NeighborCon[topNodeIDBackward].begin();in!=NeighborCon[topNodeIDBackward].end();in++){
            for(auto in=Tree[rank[topNodeIDBackward]].vert.begin();in!=Tree[rank[topNodeIDBackward]].vert.end();in++){
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

//Case 2: one core, one tree
int Graph::QueryPartiCore_CH(int ID1, int ID2){//ID1 partition, ID2 core
    int d=INF;

    int pid=CoreTag[ID1];
    vector<int> B=BoundVertex[pid];
    int bid;
    int dis1,dis2;

    set<int> todo;
    todo.clear();
    todo.insert(B.begin(), B.end());//get the partition vertex set
    todo.erase(B[B.size()-1]);

    map<int,int> results;
    QueryPeriphery_CH(ID1, todo, results);

    for(int i=0;i<B.size()-1;i++){
        bid=B[i];
        dis1=results[bid];
        dis2=QueryCore(bid,ID2);
        if(d>dis1+dis2)
            d=dis1+dis2;
    }

    return d;
}
//Case 3: Different trees
int Graph::QueryPartiParti_CH(int ID1, int ID2){//both are within the partition
    int d=INF;

    int pid1=CoreTag[ID1];
    int pid2=CoreTag[ID2];
    if(pid1==pid2){//if in the same partition
        cout<<"Wrong for partition-partition query!"<<endl;
        exit(1);

    }else{//if in different partitions
        //cout<<"Parti-Parti"<<endl;
        vector<int> B1=BoundVertex[pid1];
        vector<int> B2=BoundVertex[pid2];

        map<int,int> m1,m2;
        m1.clear();
        m2.clear();
        int bID1, bID2, tempdis;
        int b1,b2,d1,d2;
        set<int> todo1,todo2; todo1.clear(); todo2.clear();
        todo1.insert(B1.begin(),B1.end()); todo2.insert(B2.begin(),B2.end());
        todo1.erase(B1[B1.size()-1]); todo2.erase(B2[B2.size()-1]);
        map<int,int> result1,result2;

        QueryPeriphery_CH(ID1,todo1,result1);
        for(int i=0;i<B1.size()-1;i++){
            bID1=B1[i];
//            assert(Tree[rank[ID1]].disInf.find(bID1)!=Tree[rank[ID1]].disInf.end());
//            m1.insert(make_pair(bID1,Tree[rank[ID1]].disInf[i]));
//            m1.insert(make_pair(bID1,Tree[rank[ID1]].disInf[bID1]));
//            m1.insert({bID1,QueryPeriphery_CH(ID1, bID1)});
            m1.insert({bID1,result1[bID1]});
        }
        QueryPeriphery_CH(ID2,todo2,result2);
        for(int j=0;j<B2.size()-1;j++){
            bID2=B2[j];
//            assert(Tree[rank[ID2]].disInf.find(bID2)!=Tree[rank[ID2]].disInf.end());
//            m2.insert(make_pair(bID2,Tree[rank[ID2]].disInf[j]));
//            m2.insert(make_pair(bID2,Tree[rank[ID2]].disInf[bID2]));
//            m2.insert({bID2,QueryPeriphery_CH(ID2, bID2)});
            m2.insert({bID2,result2[bID2]});
        }

        for(int k=0;k<B1.size()-1;k++){
            bID1=B1[k];

            if(m1[bID1]>d)
                continue;

            for(int z=0;z<B2.size()-1;z++){
                bID2=B2[z];

                if(m2[bID2]>d)
                    continue;

                tempdis=m1[bID1]+QueryCore(bID1,bID2)+m2[bID2];
                if(tempdis<d){
                    d=tempdis;
                    d1=m1[bID1]; d2=m2[bID2];
                    b1=bID1; b2=bID2;
                }

            }
        }

    }

    return d;
}
//Case 4: Same tree, for original version
int Graph::QuerySameParti_CH(int ID1, int ID2){//both are within the partition
    int d=INF;

    int pid1=CoreTag[ID1];
    int pid2=CoreTag[ID2];
    if(pid1==pid2){//if in the same partition
        //cout<<"Same-Parti"<<endl;
        int temp_dis = QueryPeriphery_CH(ID1,ID2);/// d2 may be wrong sometimes
        if(temp_dis<d)//QueryH2HPartition(ID1,ID2,pid1)
            d=temp_dis;//QueryH2HPartition(ID1,ID2,pid1);
        vector<int> B=BoundVertex[pid1];
        map<int,int> m1,m2;
        m1.clear();
        m2.clear();
        vector<int> B1,B2;
        B1.clear();
        B2.clear();
        int bID,d1,d2;
        set<int> todo; todo.clear();
        todo.insert(B.begin(),B.end()); todo.erase(B[B.size()-1]);
        map<int,int> result1,result2;

        QueryPeriphery_CH(ID1,todo,result1);
        QueryPeriphery_CH(ID2,todo,result2);

        for(int i=0;i<B.size()-1;i++){
            bID=B[i];
//           assert(Tree[rank[ID1]].disInf.find(bID)!=Tree[rank[ID1]].disInf.end());
//           assert(Tree[rank[ID2]].disInf.find(bID)!=Tree[rank[ID2]].disInf.end());
//            d1=Tree[rank[ID1]].disInf[i];
//            d2=Tree[rank[ID2]].disInf[i];
//            d1=Tree[rank[ID1]].disInf[bID];
//            d2=Tree[rank[ID2]].disInf[bID];

//            d1=QueryPeriphery_CH(ID1, bID);
//            d2=QueryPeriphery_CH(ID2, bID);
            d1=result1[bID];
            d2=result2[bID];

            if(d1<d){
                B1.push_back(bID);
                m1.insert(make_pair(bID,d1));
            }
            if(d2<d){
                B2.push_back(bID);
                m2.insert(make_pair(bID,d2));
            }
        }

        int bID1, bID2, tempdis;
        if(!B1.empty() && !B2.empty()){
            for(int k=0;k<B1.size();k++){
                bID1=B1[k];
                if(m1[bID1]>d)
                    continue;
                for(int z=0;z<B2.size();z++){
                    bID2=B2[z];
                    if(m2[bID2]>d)
                        continue;
                    tempdis=m1[bID1]+QueryCore(bID1,bID2)+m2[bID2];
                    if(tempdis<d)
                        d=tempdis;
                }
            }
        }

    }else{//if in different partitions
        cout<<"Wrong for same partition query!"<<endl;
        exit(1);
    }

    return d;
}

//function for Query processing, debug version
int Graph::QueryDebug(int ID1, int ID2){
    int dis=INF;

    if(CoreTag[ID1]==-1 && CoreTag[ID2]==-1){//Case 1: both in core
        cout<<"Core-Core"<<endl;
        dis=QueryCoreDebug(ID1, ID2);

    }else if(CoreTag[ID1]==-1 && CoreTag[ID2]!=-1){//Case 2: ID2 in partition, ID1 in core
        cout<<"Core-Parti"<<endl;
        dis=QueryPartiCoreDebug(ID2, ID1);
    }else if(CoreTag[ID1]!=-1 && CoreTag[ID2]==-1){//Case 2: ID1 in partition, ID2 in core
        cout<<"Parti-Core"<<endl;
        dis=QueryPartiCoreDebug(ID1, ID2);

    }else if(CoreTag[ID1]!=-1 && CoreTag[ID2]!=-1){//both in partition

//            dis=QueryPartiParti(ID1,ID2);
        if(CoreTag[ID1] != CoreTag[ID2]){//Case 3: in different peripheries
            cout<<"Parti-Parti"<<endl;
            int d=INF;
            int b1,b2,d1,d2;//final results
            int pid1=CoreTag[ID1];
            int pid2=CoreTag[ID2];

            vector<int> B1=BoundVertex[pid1];
            vector<int> B2=BoundVertex[pid2];

            map<int,int> m1,m2;
            m1.clear();
            m2.clear();
            int bID1, bID2, tempdis;
            for(int i=0;i<B1.size()-1;i++){
                bID1=B1[i];
                assert(Tree[rank[ID1]].disInf.find(bID1)!=Tree[rank[ID1]].disInf.end());
//            m1.insert(make_pair(bID1,Tree[rank[ID1]].disInf[i]));
                m1.insert(make_pair(bID1,Tree[rank[ID1]].disInf[bID1]));
            }
            for(int j=0;j<B2.size()-1;j++){
                bID2=B2[j];
                assert(Tree[rank[ID2]].disInf.find(bID2)!=Tree[rank[ID2]].disInf.end());
//            m2.insert(make_pair(bID2,Tree[rank[ID2]].disInf[j]));
                m2.insert(make_pair(bID2,Tree[rank[ID2]].disInf[bID2]));
            }

            for(int k=0;k<B1.size()-1;k++){
                bID1=B1[k];

                if(m1[bID1]>d)
                    continue;

                for(int z=0;z<B2.size()-1;z++){
                    bID2=B2[z];

                    if(m2[bID2]>d)
                        continue;

                    tempdis=m1[bID1]+QueryCore(bID1,bID2)+m2[bID2];
                    if(tempdis<d){
                        d=tempdis;
                        b1=bID1; b2=bID2; d1=m1[bID1]; d2=m2[bID2];
                    }
                }
            }
            dis=d;
            int d_12=QueryCore(b1,b2), dDijk_s=Dijkstra(ID1,b1,Neighbor), dDijk_12=Dijkstra(b1,b2,Neighbor), dDijk_t=Dijkstra(b2,ID2,Neighbor);
            cout<<ID1<<" "<<b1<<"("<<NodeOrder[b1]<<") "<<b2<<"("<<NodeOrder[b2]<<") "<<ID2<<" : "<<d1<<" "<<d_12<<" "<<d2<<" ; "<<dDijk_s<<" "<<dDijk_12<<"("<<DijkstraCore(b1,b2)<<") "<<dDijk_t<<endl;

//                if(d1!=dDijk_s){
//                    DijkstraPath(ID1,b1);
//                }
//                if(d_12!=dDijk_12){
//                    DijkstraPath(b1,b2);
//                }
//                if(d2!=dDijk_t){
//                    DijkstraPath(b2,ID2);
//                }

        }else{//Case 4: in the same periphery
            cout<<"Same-Parti"<<endl;
//                dis= QuerySameParti(ID1,ID2);
            int d=INF;
            int b1,b2,df1,df2;
            int pid1=CoreTag[ID1];
            int pid2=CoreTag[ID2];

            int temp_dis = QueryPeripheryTree(ID1, ID2, pid1);/// d2 may be wrong sometimes
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
            for (int i = 0; i < B.size() - 1; i++) {
                bID = B[i];
                assert(Tree[rank[ID1]].disInf.find(bID) != Tree[rank[ID1]].disInf.end());
                assert(Tree[rank[ID2]].disInf.find(bID) != Tree[rank[ID2]].disInf.end());
                d1 = Tree[rank[ID1]].disInf[bID];
                d2 = Tree[rank[ID2]].disInf[bID];

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
//Case 1: both in core, debug version
int Graph::QueryCoreDebug(int ID1, int ID2){
    int hubfinal,dis1final,dis2final;
    int d=INF;

    int hub, dis1, dis2, temp;
    if(CoreTag[ID1]!=-1 ||CoreTag[ID2]!=-1){
        cout<<"Core Tag wrong! "<<CoreTag[ID1]<< " "<<CoreTag[ID2]<<endl;
        exit(1);
    }
//    cout<<NodeOrder[ID1]<<" "<<NodeOrder[ID2]<<endl;
    for(auto it=Label[ID1].begin();it!=Label[ID1].end();it++){
        hub=(*it).first;
        dis1=(*it).second;

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
//    DijkstraPathCore(ID1,ID2);
    return d;
}
//Case 2: one core, one tree, debug version
int Graph::QueryPartiCoreDebug(int ID1, int ID2){//ID1 partition, ID2 core
    int d=INF;

    assert(CoreTag[ID1]>=0 && CoreTag[ID2]==-1);

    int pid=CoreTag[ID1];
    int bid;
    int dis1,dis2;
    int bfinal=-1,d1=INF,d2=INF;

    for(auto it=Tree[rank[ID1]].disInf.begin();it!=Tree[rank[ID1]].disInf.end();++it){
        bid=it->first;
        dis1=it->second;
        dis2= QueryCore(bid,ID2);
        if(d>dis1+dis2){
            d=dis1+dis2;
            bfinal = bid;
            d1 = dis1; d2 = dis2;
        }

    }

    int dDijk_s=Dijkstra(ID1,bfinal,Neighbor);
    int dDijk_t=Dijkstra(bfinal,ID2,Neighbor);
    cout<<ID1<<" "<<bfinal<<"("<<NodeOrder[bfinal]<<") "<<ID2<<"("<<NodeOrder[ID2]<<") : "<<d1<<" "<<d2<<" ; "<<dDijk_s<<" "<<dDijk_t<<"("<<DijkstraCore(bfinal,ID2)<<") "<<endl;

//    if(d1!=dDijk_s){
//        DijkstraPath(ID1,bfinal);
//    }
//
//    if(d2!=dDijk_t){
//        DijkstraPath(bfinal,ID2);
//    }

    return d;
}

/// Index Maintenance

//function of maintaining edge weight decrease updates
void Graph::Decrease(int a, int b, int oldW, int newW){
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

    int pid = -1;

    AdjaCoreMapOld = AdjaCoreMap;
//    changedInterfaceEdges.clear();

    if((CoreTag[a]==-1 && BoundTag[a].first==0) || (CoreTag[b]==-1 && BoundTag[b].first==0)){// either endpoint is non-boundary core vertex
        //cout<<"******************change 1*******************"<<endl;
//        cout<<"Core-Core"<<endl;
        if(algoCoreU==0) {
//		DecreasePLL(a,b,oldW,newW,AdjaCore,Label);
            DecreasePSL(a, b, oldW, newW, AdjaCore, Label);
        }else if(algoCoreU==1){
            PLLdec(a, b, oldW, newW, AdjaCore);
        }
    }else if(CoreTag[a]!=-1 || CoreTag[b]!=-1){//either endpoint is in-partition vertex
        //cout<<"******************change 2*******************"<<endl;
        if(CoreTag[a]!=-1)//if a is periphery vertex
            pid=CoreTag[a];
        else
            pid=CoreTag[b];

        //cout<<"<<<<<<<<<<"<<endl;
        //cout<<CoreTag[a]<<" "<<BoundTag[a]<<" "<<CoreTag[b]<<" "<<BoundTag[b]<<" pid "<<pid<<endl;
        //cout<<"decrease "<<a<<" "<<b<<" "<<oldW<<" "<<newW<<endl;
        if(algoTree==0){
            DecreaseCHNew(a,b,newW,Neighbor,Tree,rank,heightMax);
        }else if(algoTree==1){
            DecreaseH2HNew(a,b,newW,Neighbor,Tree,rank,heightMax,true);
//        DecreaseH2HNew(a,b,newW,Neighbor,Tree,rank,heightMax,false);
        }


        //cout<<">>>>>>>>>>"<<endl;

        /// check whether the update of periphery has affected core
        int ID1,ID2,dis,olddis;
//        for(auto it=changedInterfaceEdges.begin();it!=changedInterfaceEdges.end();++it){
//            ID1=it->first;
//            for(auto it2=it->second.begin();it2!=it->second.end();++it2){
//                ID2=it2->first; olddis=it2->second;
//                dis = AdjaCoreMap[ID1][ID2];
//                if(olddis>dis){
//                    if(algoCoreU==0) {
////                   DecreasePLL(ID1,ID2,olddis,dis,AdjaCore,Label);
//                        DecreasePSL(ID1, ID2, olddis, dis, AdjaCore, Label);
//                    }else if(algoCoreU==1){
//                        PLLdec(ID1, ID2, olddis, dis, AdjaCore);
//                    }
//                    extUpdate = true;
//                }
//            }
//        }
        for(int i=0;i<BoundVertex[pid].size()-1;i++){
            ID1=BoundVertex[pid][i];

            //if(PartiVertexInverted[pid].find(ID1)!=PartiVertexInverted[pid].end()){//only if ID1 has neighbors in partition pid
            for(int j=i+1;j<BoundVertex[pid].size()-1;j++){
                ID2=BoundVertex[pid][j];

                //if(PartiVertexInverted[pid].find(ID2)!=PartiVertexInverted[pid].end()){
//                dis = QueryPeripheryTree(ID1,ID2,pid);///
                dis = AdjaCoreMap[ID1][ID2];
//                if(dis != temp_dis)
//                    cout<<dis<<" "<<temp_dis<<endl;
                olddis=AdjaCoreMapOld[ID1][ID2];
                if(olddis>dis){
                    if(algoCoreU==0) {
//                   DecreasePLL(ID1,ID2,olddis,dis,AdjaCore,Label);
                        DecreasePSL(ID1, ID2, olddis, dis, AdjaCore, Label);
                    }else if(algoCoreU==1){
                        PLLdec(ID1, ID2, olddis, dis, AdjaCore);
                    }

                }


            }
//
        }
    }else if(BoundTag[a].first==1 && BoundTag[b].first==1){//Both end points are interface vertex
        //cout<<"******************change 3*******************"<<endl;
        ///periphery update
//        cout<<"Both interface vertex!"<<endl;
        int olddis=AdjaCoreMap[a][b];
        if(olddis>newW){
            AdjaCoreMap[a][b]=newW;
            AdjaCoreMap[b][a]=newW;
            if(algoCoreU==0) {
//            DecreasePLL(a,b,olddis,newW,AdjaCore,Label);
                DecreasePSL(a, b, olddis, newW, AdjaCore, Label);
            }else if(algoCoreU==1){
                PLLdec(a, b, olddis, newW, AdjaCore);
            }
        }
    }

}
//New version, 2023-02-28
void Graph::Increase(int a, int b, int oldW, int newW){
    for(int i=0;i<Neighbor[a].size();i++){
        if(Neighbor[a][i].first==b){
            if(oldW != Neighbor[a][i].second){
                cout<<"Old edge weight is incorrect! "<<a<<" "<<b<<": "<<oldW<<" "<<Neighbor[a][i].second<<endl;
                oldW = Neighbor[a][i].second;
            }
//            cout<<a<<" "<<b<<" "<<Neighbor[a][i].second<<" "<<newW<<endl;
            Neighbor[a][i].second=newW;
            break;
        }
    }
    for(int i=0;i<Neighbor[b].size();i++){
        if(Neighbor[b][i].first==a){
            if(oldW != Neighbor[b][i].second){
                cout<<"Old edge weight is incorrect! "<<b<<" "<<a<<": "<<oldW<<" "<<Neighbor[b][i].second<<endl;
                oldW = Neighbor[b][i].second;
            }
//            cout<<b<<" "<<a<<" "<<Neighbor[b][i].second<<" "<<newW<<endl;
            Neighbor[b][i].second=newW;
            break;
        }
    }

//    AdjaCoreMapOld = AdjaCoreMap;
//    NodeOrder_ = NodeOrder;///

    int pid = -1;
    if((CoreTag[a]==-1 && !BoundTag[a].first) || (CoreTag[b]==-1 && !BoundTag[b].first)){//edges in the core
//        cout<<"edge in core"<<endl;
        if(algoCoreU==0){
    //        IncreasePSL(a,b,oldW,newW,AdjaCore,Label,PruningPointNew,NoSupportedPair);
//            IncreasePSLNew(a,b,oldW,newW,AdjaCore,Label,PruningPointSetOrder);//set version with NoSupportedPair, queue
            IncreasePSLNew(a,b,oldW,newW,AdjaCore,Label,PruningPointSet);//set version with NoSupportedPair, queue
        }else if(algoCoreU==1){
            PLLinc(a,b,oldW,newW,AdjaCore);
        }

    }else if(CoreTag[a]!=-1 || CoreTag[b]!=-1){//edges in one partition, either endpoint is in-partition vertex
//        cout<<"edge in partition"<<endl;

        if(CoreTag[a]!=-1)
            pid=CoreTag[a];
        else
            pid=CoreTag[b];

        int ID1,ID2,dis,olddis;
        AdjaCoreMapOld = AdjaCoreMap;


        //cout<<"<<<<<<<<<<"<<endl;
        //cout<<CoreTag[a]<<" "<<BoundTag[a]<<" "<<CoreTag[b]<<" "<<BoundTag[b]<<" pid "<<pid<<endl;
        if(algoTree==0){
            IncreaseCHNew(a,b,oldW,newW,Neighbor,Tree,rank,heightMax,SCconNodesMT,VidtoTNid);
        }else if(algoTree==1){
            IncreaseH2HNew(a,b,oldW,newW,Neighbor,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,true);
//        IncreaseH2HNew(a,b,oldW,newW,Neighbor,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,false);
        }

        //cout<<">>>>>>>>>>"<<endl;

//            cout<<"Flag 2"<<endl;

        for(int i=0;i<BoundVertex[pid].size()-1;i++){
            ID1=BoundVertex[pid][i];

            //if(PartiVertexInverted[pid].find(ID1)!=PartiVertexInverted[pid].end()){//only if ID1 has neighbors in partition pid
            for(int j=i+1;j<BoundVertex[pid].size()-1;j++){
                ID2=BoundVertex[pid][j];

                dis = AdjaCoreMap[ID1][ID2];
                olddis=AdjaCoreMapOld[ID1][ID2];
                if(olddis < dis){
//                    cout<<ID1<<" "<<ID2<<" "<<olddis<<" "<<dis<<endl;
                    if(algoCoreU==0) {
//                    IncreasePSL(ID1,ID2,olddis,dis,AdjaCore,Label,PruningPointNew,NoSupportedPair);//original version
//                        IncreasePSLNew(ID1, ID2, olddis, dis, AdjaCore, Label,PruningPointSetOrder);//set version with NoSupportedPair, queue
                        IncreasePSLNew(ID1, ID2, olddis, dis, AdjaCore, Label, PruningPointSet);//set version with NoSupportedPair, queue
                    }else if(algoCoreU==1){
                        PLLinc(ID1,ID2,olddis,dis,AdjaCore);
                    }
                }

            }
        }

    }else{//Both end points are boundary vertex
        //cout<<CoreTag[a]<<" "<<BoundTag[a]<<" "<<CoreTag[b]<<" "<<BoundTag[b]<<endl;
//        cout<<"edge between boundary vertex"<<endl;
        int olddis=AdjaCoreMap[a][b];
        assert(olddis <= oldW);
        if(olddis==oldW){//it indicates the update of e(a,b) may influence AdjaCoreMap[a][b]
            int newDis=INF;//the distance supported by partitions
            set<int> newSets;

            int lid,hid;
            if(NodeOrder[a]<NodeOrder[b]){
                lid=a;hid=b;
            }else{
                lid=b;hid=a;
            }
//            cout<<"Super edge: "<<lid<<" "<<hid<<endl;
//            assert(SuppPartiID[lid].find(hid)!=SuppPartiID[lid].end());
//            assert(!SuppPartiID[lid][hid].empty());

            int Cw=INF; int countwt=0;
            for(int i=0;i<Neighbor[lid].size();i++){
                if(Neighbor[lid][i].first==hid){
                    Cw=Neighbor[lid][i].second;//the weight value in the original graph
                    countwt=1;
                    break;
                }
            }
            int ssw=INF,wtt=INF,wid=-1;
            vector<pair<int,int>> Wnodes; //Wnodes.clear(); //Supportive vertices

            /// calculate the correct value by traversing supportive vertices
//            assert(SCconNodesMT[lid][hid].size() == SCconNodesMT[hid][lid].size());
            if(lid<hid){
                Wnodes=SCconNodesMT[lid][hid]; //cout<<"wid num "<<Wnodes.size()<<endl;
            }else{
                Wnodes=SCconNodesMT[hid][lid]; //cout<<"wid num "<<Wnodes.size()<<endl;
            }


            for(int i=0;i<Wnodes.size();i++){
                wid=Wnodes[i].first;//wid must be periphery vertex
                assert(CoreTag[wid] != -1);

                pid = CoreTag[wid];

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
                    newDis=Cw;
                    newSets.clear();
                    newSets.insert(pid);
                }else if(ssw+wtt==Cw){
                    countwt+=1;
                    newSets.insert(pid);
                }

                if(lid<hid){
                    if(SuppPartiID[lid][hid].find(pid) != SuppPartiID[lid][hid].end()){
                        if(SuppPartiID[lid][hid][pid] != ssw+wtt){//values of SuppPartiID may be wrong
                            SuppPartiID[lid][hid][pid] = ssw+wtt;
                        }
                    }
                }else{
                    if(SuppPartiID[hid][lid].find(pid) != SuppPartiID[hid][lid].end()){
                        if(SuppPartiID[hid][lid][pid] != ssw+wtt){//values of SuppPartiID may be wrong
                            SuppPartiID[hid][lid][pid] = ssw+wtt;
                        }
                    }
                }

            }

//            for(auto it=SuppPartiID[a][b].begin();it!=SuppPartiID[a][b].end();++it){//for the partitions that contains a and b
//                if(it->second < newDis){
//                    newDis = it->second;
//                    newSets.clear();
//                    newSets.insert(it->first);
//                }else if(it->second == newDis){
//                    newSets.insert(it->first);
//                }
//            }
            if(newDis > oldW){//trigger update
                if(newDis <= newW){//if the partition-supported edge weight is smaller than newW, we update edge weight from oldW to newDis, <
//                    cout<<"Supported by partition: "<<oldW<<" "<<newDis<<" "<<newW<<", "<<newSets.size()<<" "<<SuppPartiID[a][b].size()<<" "<<SuppPartiIDReal[a][b].first<<" "<<SuppPartiIDReal[a][b].second.size()<<endl;
                    AdjaCoreMap[a][b] = newDis;
                    AdjaCoreMap[b][a] = newDis;
                    if(a<b){
                        SuppPartiIDReal[a][b].first = newDis;
                        SuppPartiIDReal[a][b].second = newSets;
                        /// update interface entries
                        for(auto it=SuppPartiID[a][b].begin();it!=SuppPartiID[a][b].end();++it){//for the partitions that contains a and b
//                        cout<<"Partition "<<it->first<<" "<<it->second<<endl;
                            int Pid = it->first;
                            int rootID = BoundVertex[Pid][BoundVertex[Pid].size()-1];
                            vector<int> interfaceP;
                            for(auto it2=Tree[rank[rootID]].vert.begin();it2!=Tree[rank[rootID]].vert.end();++it2){
                                interfaceP.emplace_back(it2->first);
                            }
                            InterfacePropagate( rank[rootID], interfaceP,Tree,true);
                        }
                    }else{
                        SuppPartiIDReal[b][a].first = newDis;
                        SuppPartiIDReal[b][a].second = newSets;
                        /// update interface entries
                        for(auto it=SuppPartiID[b][a].begin();it!=SuppPartiID[b][a].end();++it){//for the partitions that contains a and b
//                        cout<<"Partition "<<it->first<<" "<<it->second<<endl;
                            int Pid = it->first;
                            int rootID = BoundVertex[Pid][BoundVertex[Pid].size()-1];
                            vector<int> interfaceP;

                            for(auto it2=Tree[rank[rootID]].vert.begin();it2!=Tree[rank[rootID]].vert.end();++it2){
                                interfaceP.emplace_back(it2->first);
                            }
                            InterfacePropagate( rank[rootID], interfaceP,Tree,true);
                        }
                    }


                    /// update the core index
                    if(algoCoreU==0) {
//                    IncreasePSL(a,b,oldW,newDis,AdjaCore,Label,PruningPointNew,NoSupportedPair);//original version
//                        IncreasePSLNew(a, b, oldW, newDis, AdjaCore, Label,PruningPointSetOrder);//set version with NoSupportedPair, queue
                        IncreasePSLNew(a, b, oldW, newDis, AdjaCore, Label, PruningPointSet);//set version with NoSupportedPair, queue
                    }else if (algoCoreU==1){
                        PLLinc(a, b, oldW, newDis, AdjaCore);
                    }

                }
                else{//if newW < newDis
//                    cout<<"Supported by newW: "<<oldW<<" "<<newDis<<" "<<newW<<", "<<newSets.size()<<" "<<SuppPartiID[a][b].size()<<" "<<SuppPartiIDReal[a][b].first<<" "<<SuppPartiIDReal[a][b].second.size()<<endl;
                    AdjaCoreMap[a][b] = newW;
                    AdjaCoreMap[b][a] = newW;
                    if(a<b){
                        SuppPartiIDReal[a][b].first = newW;
                        SuppPartiIDReal[a][b].second.clear();
                        /// update interface entries
//                        for(auto it=SuppPartiID[a][b].begin();it!=SuppPartiID[a][b].end();++it){//for the partitions that contains a and b
////                        cout<<"Partition "<<it->first<<" "<<it->second<<endl;
//                            int Pid = it->first;
//                            int rootID = BoundVertex[Pid][BoundVertex[Pid].size()-1];
//                            vector<int> interfaceP;
//
//                            for(auto it2=Tree[rank[rootID]].vert.begin();it2!=Tree[rank[rootID]].vert.end();++it2){
//                                interfaceP.emplace_back(it2->first);
//                            }
//                            InterfacePropagate( rank[rootID], interfaceP,Tree,true);
//                        }
                    }else{
                        SuppPartiIDReal[b][a].first = newW;
                        SuppPartiIDReal[b][a].second.clear();
                        /// update interface entries
//                        for(auto it=SuppPartiID[b][a].begin();it!=SuppPartiID[b][a].end();++it){//for the partitions that contains a and b
////                        cout<<"Partition "<<it->first<<" "<<it->second<<endl;
//                            int Pid = it->first;
//                            int rootID = BoundVertex[Pid][BoundVertex[Pid].size()-1];
//                            vector<int> interfaceP;
//
//                            for(auto it2=Tree[rank[rootID]].vert.begin();it2!=Tree[rank[rootID]].vert.end();++it2){
//                                interfaceP.emplace_back(it2->first);
//                            }
//                            InterfacePropagate( rank[rootID], interfaceP,Tree,true);
//                        }
                    }

                    /// update core index
                    if(algoCoreU==0) {
//                    IncreasePSL(a,b,oldW,newW,AdjaCore,Label,PruningPointNew,NoSupportedPair);//original version
//                        IncreasePSLNew(a, b, oldW, newW, AdjaCore, Label,PruningPointSetOrder);//set version with NoSupportedPair, queue
                        IncreasePSLNew(a, b, oldW, newW, AdjaCore, Label,PruningPointSet);//set version with NoSupportedPair, queue
                    }else if(algoCoreU==1){
                        PLLinc(a, b, oldW, newW, AdjaCore);
                    }
                }
            }

        }

    }


}









