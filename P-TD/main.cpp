/*
 * main.cpp
 *
 *  Created on: 19 Sep 2022
 *      Author: zhangmengxuan, Xinjie ZHOU
 */
#include "head.h"
int ifDebug = false;

int main(int argc, char** argv){
    if( argc != 5 && argc != 6){
        printf("usage:\n<arg1> source path, e.g. /export/project/xzhouby\n");
        printf("<arg2> dataset, e.g. NY\n");
        printf("<arg3> partition method, e.g. NC\n");
        printf("<arg4> partition number, e.g. 64\n");
        printf("<arg5> (optional) update type, e.g. 0: Decrease; 1: Increase.\n");
        exit(0);
    }
    cout<<"This is test for Edge-cut H2H!"<<endl;
    string partitionPath;
//    //read the partition files

//	string graphfile="/media/TraminerData/mengxuan/PartitionData/NY";
    string graphfile="/Users/zhouxj/Documents/1-Research/Datasets/NY";
    string sourcePath = "/Users/zhouxj/Documents/1-Research/Datasets/";

    int updateType = DECREASE;
    int updateBatch = 100;
    int updateVolume = 1;
    int runtimes = 1000;

	Graph g;
    g.algoName = "Bubble";
    //	g.pnum=75;
    g.pnum=6;
    g.threadnum=10;//thread number

    if(argc > 1){
        sourcePath = argv[1];
        cout<<"argv[1]: "<<argv[1]<<endl;
        g.dataset = argv[2];
        cout<<"argv[2]: "<<argv[2]<<endl;
        g.algoName = argv[3];
        cout<<"argv[3]: "<<argv[3]<<endl;
        g.pnum = atoi(argv[4]);
        cout<<"argv[4]: "<<argv[4]<<endl;
        if(argc > 5){
            updateType = atoi(argv[5]);
            cout<<"argv[5]: "<<argv[5]<<endl;
        }
    }
    cout<<"Dataset: "<<g.dataset<<endl;
    cout<<"Partition method: "<<g.algoName<<endl;
    cout<<"Thread number: "<<g.threadnum<<endl;
    if(ifDebug){
        cout<<"Debug mode!"<<endl;
    }else{
        cout<<"Test mode!"<<endl;
    }

    graphfile = sourcePath + "/" + g.dataset + "/" +g.dataset;
    partitionPath = graphfile + "_" + g.algoName+"_"+ to_string(g.pnum);

    g.GraphRead(graphfile);
    g.OrderRead(partitionPath+"/vertex_order");

	//Read partitioned subgraphs
    g.GraphPartitionRead(partitionPath);

    cout<<"Finished reading..."<<endl;

    std::chrono::high_resolution_clock::time_point t10;
    std::chrono::high_resolution_clock::time_point t11;
    std::chrono::duration<double> time_span1;
    double runT1=0;
    double runT2=0;
    double runT3=0;
    double maxrT=0;

    Timer tt;

    //for pre-boundary
//    g.PreBoundaryCompute(true);
//    exit(0);

    /// Task1: Index construction
    tt.start();
    cout<<"Index construction test..."<<endl;
    //Index construction of subgraphs
    g.Trees.clear();
    g.toRMQs.clear();
    g.RMQIndexs.clear();
    g.ranks.clear();
    g.heightMaxs.clear();
    //g.SCconNodess.clear();
    g.SCconNodesMTs.clear();
    g.VidtoTNids.clear();
    cout<<"Partition index construction..."<<endl;
    if(ifDebug)
        cout<<"partition ID ";
    for(int partiID=0;partiID<g.pnum;partiID++){
        if(ifDebug)
            cout<<partiID<<" ";
        Graph h;
        h.nodenum=g.nodenum;
        h.NodeOrder=g.NodeOrder;
        h.vNodeOrder=g.vNodeOrder;
        h.Neighbors=g.NeighborsParti[partiID];

        t10=std::chrono::high_resolution_clock::now();
        h.H2Hindex();
        t11=std::chrono::high_resolution_clock::now();
        time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
        runT1+= time_span1.count();
        if(maxrT < time_span1.count()){
            maxrT = time_span1.count();
        }
        g.Trees.push_back(h.Tree);
        g.toRMQs.push_back(h.toRMQ);
        g.RMQIndexs.push_back(h.RMQIndex);
        g.ranks.push_back(h.rank);
        g.heightMaxs.push_back(h.heightMax);
        //g.SCconNodess.push_back(h.SCconNodes);
        g.SCconNodesMTs.push_back(h.SCconNodesMT);
        g.VidtoTNids.push_back(h.VidtoTNid);
    }
    if(ifDebug)
        cout<<endl;
    cout<<"Partitions' index construction time: "<<runT1<< " s. maxrT: "<<maxrT<<" s."<<endl;

    Timer tt1;
    tt1.start();
    //Build overlay graph
    g.OverlayGraphConstructPost();
    tt1.stop();
    runT2 = tt1.GetRuntime();
    cout<<"Overlay graph construction time: "<< tt1.GetRuntime()<<" s."<<endl;
    if(ifDebug){
        g.CorrectnessCheckOverlay();
    }


    //Index construction of overlay graph
    Graph l;
    l.Neighbors=g.NeighborsOverlay;
    l.nodenum=g.nodenum;
    l.NodeOrder=g.NodeOrder;
    l.vNodeOrder=g.vNodeOrder;

    t10=std::chrono::high_resolution_clock::now();
    l.H2Hindex();
    t11=std::chrono::high_resolution_clock::now();
    time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
    runT3 = time_span1.count();
    cout<<"Overlay Index Construction Time: "<<runT3<<" s."<<endl;
    g.TreeOverlay=l.Tree;
    g.toRMQOverlay=l.toRMQ;
    g.RMQIndexOverlay=l.RMQIndex;
    g.rankOverlay=l.rank;
    g.heightMaxOverlay=l.heightMax;
    //g.SCconNodesOverlay=l.SCconNodes;
    g.SCconNodesOverlayMT=l.SCconNodesMT;
    g.VidtoTNidOverlay=l.VidtoTNid;

    if(ifDebug){
        //Correctness Check
        g.CorrectnessCheck(100);
    }
    tt.stop();

    cout<<"Overall Index Construction Time: "<<runT1+runT2+runT3<<" s."<<endl;
    cout<<"Minimum Overall Index Construction Time: "<<maxrT+runT3<<" s."<<endl;
    cout<<"All Time: "<<tt.GetRuntime()<<" s."<<endl;
    g.Indexsize();
    cout<<"-----------"<<endl;

    ///Task 2: query processing
    cout<<"Query processing test..."<<endl;
    //Query processing efficiency
    g.EffiCheck(graphfile+".query", runtimes);///Efficiency test

    g.EffiCheck();///Query distribution test
    cout<<"-----------"<<endl;
//    exit(0);
    ///Task3: Index update
    cout<<"Index update test..."<<endl;
    // read updates
    string file = graphfile + ".update";

    vector<pair<pair<int,int>,pair<int,int>>> wBatch;
    int ID1, ID2, oldW, newW;
    srand (time(NULL));
    cout<<"Update batch: "<<updateBatch<<endl;

    switch (updateType) {
        case 0:{
            //Decrease update
            cout<<"Update type: Decrease"<<endl;
            Graph g1=g;
            g1.ReadUpdates(file);
            t10=std::chrono::high_resolution_clock::now();
            for(int u=0;u<updateBatch;u++){
                ID1 = g1.updateEdges[u].first.first;
                ID2 = g1.updateEdges[u].first.second;
                oldW = g1.updateEdges[u].second;
                newW=oldW*0.5;
                if(newW < 1) {
                    cout<<"New edge weight is not positive! "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                    exit(1);
                }
                if(ifDebug){
                    cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                }

                g1.DecreaseSingle(ID1,ID2,newW);
                if(ifDebug){
                    g1.CorrectnessCheck(100);
                }

            }
            t11=std::chrono::high_resolution_clock::now();
            time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
            runT1= time_span1.count();
            cout<<"Average Decrease update Time: "<<runT1/updateBatch<<" s."<<endl;
            break;
        }
        case 1:{
            //Increase update
            cout<<"Update type: Increase"<<endl;
//            Graph g2=g;
            g.ReadUpdates(file);
            t10=std::chrono::high_resolution_clock::now();

            for(int u=0;u<updateBatch;u++){
                ID1 = g.updateEdges[u].first.first;
                ID2 = g.updateEdges[u].first.second;
                oldW = g.updateEdges[u].second;
                newW=oldW*2;
                if(ifDebug)
                    cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                g.IncreaseSingle(ID1,ID2,oldW,newW);
                if(ifDebug){
                    g.CorrectnessCheck(100);
                }

            }
            t11=std::chrono::high_resolution_clock::now();
            time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
            runT1= time_span1.count();
            cout<<"Average Increase update Time: "<<runT1/updateBatch<<" s."<<endl;
            break;
        }
        default:{
            cout<<"Wrong update type!"<<endl;
            break;
        }

    }

    cout<<"-----------"<<endl;

	return 0;
}




