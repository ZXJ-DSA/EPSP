/*
 * main.cpp
 *
 *  Created on: 19 Sep 2022
 *      Author: zhangmengxuan
 */
#include "head.h"
int ifDebug = false;

int main(int argc, char** argv){
    if( argc != 6 && argc != 7){
        printf("usage:\n<arg1> source path, e.g. /export/project/xzhouby\n");
        printf("<arg2> dataset, e.g. NY\n");
        printf("<arg3> partition method, e.g. Bubble\n");
        printf("<arg4> partition number, e.g. 4\n");
        printf("<arg5> task, e.g. 1: index construction; 2: query processing; 3: index update.\n");
        printf("<arg6> (only for update) update type, e.g. 0: Increase; 1: Decrease.\n");
        exit(0);
    }
    cout<<"This is test for Edge-cut H2H!"<<endl;
    string partitionPath;
//    //read the partition files

//	string graphfile="/media/TraminerData/mengxuan/PartitionData/NY";
    string graphfile="/Users/zhouxj/Documents/1-Research/Datasets/NY";
    string sourcePath = "/Users/zhouxj/Documents/1-Research/Datasets/";

    int updateType = INCREASE;
    int updateBatch = 100;
    int updateVolume = 1;
    int runtimes = 1000;

	Graph g;
    g.algoName = "Bubble";
    //	g.pnum=75;
    g.pnum=6;
    g.threadnum=150;//thread number
    int taskType = 1;
    if(argc > 1){
        sourcePath = argv[1];
        g.dataset = argv[2];
        g.algoName = argv[3];
        g.pnum = atoi(argv[4]);
        taskType = atoi(argv[5]);
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
//	g.GraphRead(graphfile+"/Graph");
//	g.OrderRead(graphfile+"/Order75");
//    g.GraphRead(graphfile+"/NY");
//    g.OrderRead(graphfile+"/NY_NC_6/vertex_order");
    g.GraphRead(graphfile);
    g.OrderRead(partitionPath+"/vertex_order");


	//Read partitioned subgraphs
//	g.GraphPartitionRead(graphfile+"/partition"+to_string(g.pnum));
//    g.GraphPartitionRead(graphfile+"/NY_NC_6");
    g.GraphPartitionRead(partitionPath);

    cout<<"Finished reading..."<<endl;

    std::chrono::high_resolution_clock::time_point t10;
    std::chrono::high_resolution_clock::time_point t11;
    std::chrono::duration<double> time_span1;
    double runT1=0;

    if(argc > 6){
        updateType = atoi(argv[6]);
        cout<<"updateType: "<<updateType<<endl;
    }
    Timer tt;
    switch (taskType) {
        case 1:/// Index construction
        {
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
            cout<<"Partitions' index construction time: "<<runT1<< " s."<<endl;

            Timer tt1;
            tt1.start();
            //Build overlay graph
            g.OverlayGraphConstructPost();
            tt1.stop();
            cout<<"Overlay graph construction time: "<< tt1.GetRuntime()<<" s."<<endl;
            if(ifDebug){
                g.CorrectnessCheckOverlay();
            }

//            cout<<"Overlay graph finish construction!"<<endl;

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
            runT1= time_span1.count();
            cout<<"Overlay Index Construction Time: "<<runT1<<" s."<<endl;
            g.TreeOverlay=l.Tree;
            g.toRMQOverlay=l.toRMQ;
            g.RMQIndexOverlay=l.RMQIndex;
            g.rankOverlay=l.rank;
            g.heightMaxOverlay=l.heightMax;
            //g.SCconNodesOverlay=l.SCconNodes;
            g.SCconNodesOverlayMT=l.SCconNodesMT;
            g.VidtoTNidOverlay=l.VidtoTNid;

//            if(ifDebug){
//                //Correctness Check
//                g.CorrectnessCheck();
//            }
            tt.stop();
            cout<<"Time used for index construction: "<<tt.GetRuntime()<<" s."<<endl;
            g.Indexsize();
            cout<<"-----------"<<endl;
//            break;
        }
        case 2:///query processing
        {
            cout<<"Query processing test..."<<endl;
            if(ifDebug){
                //Correctness Check
                g.CorrectnessCheck(100);
            }
            //Query processing efficiency
//        g.EffiCheck(graphfile+"/ODeffi");
            g.EffiCheck(graphfile+".query", runtimes);

            g.EffiCheck();
            cout<<"-----------"<<endl;
            //break;
        }
        case 3:///Index update
        {
            cout<<"Index update test..."<<endl;
//            if(argc != 7){
//                cout<<"Wrong input parameter!"<<endl;
//                exit(1);
//            }
            // read updates
            string file = graphfile + ".update";

            vector<pair<pair<int,int>,pair<int,int>>> wBatch;
            int ID1, ID2, oldW, newW;
            srand (time(NULL));
            cout<<"Update batch: "<<updateBatch<<endl;
            if(updateType == INCREASE){
                cout<<"Update type: Increase"<<endl;
                Graph g2=g;
                g2.ReadUpdates(file);
                t10=std::chrono::high_resolution_clock::now();

                /// for batch update
//                wBatch.clear();
//                for(int u=0;u<updateBatch;u++) {
//                    if (ifDebug){
//                        cout << "Batch " << u << ": " << ID1 << " " << ID2 << " " << oldW << " " << newW << endl;
//                    }
//                    for(int i=u*updateVolume;i<(u+1)*updateVolume;++i) {//for each edge update
//                        ID1 = g2.updateEdges[i].first.first;
//                        ID2 = g2.updateEdges[i].first.second;
//                        oldW = g2.updateEdges[i].second;
//                        newW = oldW * 1.5;
//                        wBatch.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
//                    }
//                    if (ifDebug)
//                        cout << "Batch " << u << ": " << ID1 << " " << ID2 << " " << oldW << " " << newW << endl;
//                    g2.IncreaseBatch(wBatch,g2.Neighbors,g2.Tree,g2.rank,g2.heightMax,g2.SCconNodesMT,g2.VidtoTNid);
//                    if (ifDebug) {
//                        g2.CorrectnessCheck(runtimes);
//                    }
//                }
                /// for single update
                for(int u=0;u<updateBatch;u++){
//                    ID1=rand()%g2.nodenum;
//                    ID2=g2.Neighbors[ID1][0].first;
//                    oldW=g2.Neighbors[ID1][0].second;
                    ID1 = g2.updateEdges[u].first.first;
                    ID2 = g2.updateEdges[u].first.second;
                    oldW = g2.updateEdges[u].second;
                    newW=oldW*2;
                    if(ifDebug)
                        cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                    g2.IncreaseSingle(ID1,ID2,oldW,newW);
                    if(ifDebug){
                        g2.CorrectnessCheck(100);
                    }

                }
                t11=std::chrono::high_resolution_clock::now();
                time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
                runT1= time_span1.count();
                cout<<"Average Increase update Time: "<<runT1/updateBatch<<" s."<<endl;
            }else if(updateType == DECREASE){
                cout<<"Update type: Decrease"<<endl;
                Graph g1=g;
                g1.ReadUpdates(file);
                t10=std::chrono::high_resolution_clock::now();
                for(int u=0;u<updateBatch;u++){
//                    ID1=rand()%g1.nodenum;
//                    ID2=g1.Neighbors[ID1][0].first;
//                    oldW=g1.Neighbors[ID1][0].second;
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
            }else{
                cout<<"Wrong update type!"<<endl;
            }
            cout<<"-----------"<<endl;
            break;
        }
        default:
            break;
    }

	return 0;
}




