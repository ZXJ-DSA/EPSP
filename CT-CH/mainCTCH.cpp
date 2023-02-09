/*
 * mainCTCH.cpp
 *
 *  Created on: 14 Jan 2023
 *      Author: Xinjie ZHOU
 */
#include "head3.h"

bool ifExtension = false;

int main(int argc, char** argv){

    if( argc != 4 && argc != 5 && argc != 6){//
        printf("usage:\n<arg1> source path, e.g. /export/project/xzhouby\n");
        printf("<arg2> name of dataset, e.g. NY\n");
        printf("<arg3> tree width, e.g. 20\n");
        printf("<arg4> (optional) update type, (0: Decrease; 1: Increase) e.g. 0\n");
        printf("<arg5> (optional) if use extended label, (0: false; 1: True), default: 0\n");
        exit(0);
    }
    cout<<"This is test for CT-CH!"<<endl;
    string DesFile="./data/";
    string dataset = "NY";
    int treeWidth = 20;
    int updateType = 0;
    int runtimes = 1000;
    int updateBatch = 100;//200;

    if(argc > 1) {
        cout << "argc: " << argc << endl;
        cout << "argv[1]: " << argv[1] << endl;//source path
        DesFile = argv[1];
        if(argc > 2){
            cout << "argv[2]: " << argv[2] << endl;//dataset
            dataset = argv[2];
        }
        if(argc > 3){
            cout << "argv[3]: " << argv[3] << endl;
            treeWidth = stoi(argv[3]);
        }
        if(argc > 4){
            cout << "argv[4]: " << argv[4] << endl;
            updateType = stoi(argv[4]);
        }
        if(argc > 5){
            cout << "argv[5]: " << argv[5] << endl;
            ifExtension = stoi(argv[5]);
        }
    }

	//used for running time calculation
	std::chrono::high_resolution_clock::time_point t10;
	std::chrono::high_resolution_clock::time_point t11;
	std::chrono::duration<double> time_span1;
	double runT1=0;
    double runT2=0;
    double runT3=0;

//    string graphfile="/media/TraminerData/mengxuan/MengxuanGraphWPSL/Cond/CondWeighted";
    string graphfile=DesFile+"/"+dataset+"/"+dataset;
//    string orderfile=graphfile+".order";
    string ODfile=graphfile+".query";
    string updateFile=graphfile+".update";

    Graph g;
    g.threadnum=10;//thread number of parallel computation (can be changed)
    g.ifParallel = true;
    g.Width=treeWidth;//15;
    g.dataset=dataset;
    cout<<"Dataset: "<<dataset<<endl;
    cout<<"Tree width: "<<g.Width<<endl;
    cout<<"Thread number: "<<g.threadnum<<endl;
    int ID1,ID2, oldW,newW;

//    g.CoreGraphDebug(graphfile+"Core2");

    g.ReadGraph(graphfile);
//    g.StainingMethod(0);


    ///Task 1: Index construction
    t10=std::chrono::high_resolution_clock::now();
    ///tree index construction
    //g.Construct_tree(true);
    //CH index construction
    g.MDEContract();
    t11=std::chrono::high_resolution_clock::now();
    time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
    runT1 = time_span1.count();

    g.Create_Partitions();
    cout<<"Partition's index construction time: "<<runT1<<" s"<<endl;

    t10=std::chrono::high_resolution_clock::now();
    ///core index construction
    g.Construct_core();

    t11=std::chrono::high_resolution_clock::now();
    time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
    runT2 = time_span1.count();
    cout<<"Core's index construction time: "<<runT2<<" s."<< endl;
	cout<<"-----------------"<<endl;


    cout<<"Overall Construction Time: "<<runT1+runT2+runT3<<" s."<<endl;
//    cout<<"Minimum Overall Construction Time: "<<runT3+runT4+maxrT<<" s."<<endl;

    g.indexsizeCHH2H();//index (labeling+pruning point) size computation
    g.CorrectnessCheck_CH(200);
//    exit(0);
    ///Task 2: Query processing
//    g.EffiCheck("/media/TraminerData/mengxuan/MengxuanGraphWPSL/Cond/ODpair");//query efficiency test
    g.EffiCheck_CH(ODfile,runtimes);//query efficiency test
//    exit(0);
    ///Task 3: Index update
    vector<pair<pair<int,int>,int>> testdata;
    g.ReadUpdate(updateFile, testdata);

    cout<<"Update Batch: "<<updateBatch<<endl;
//    int ID1, ID2, oldW, newW;
//    srand (time(NULL));
    switch (updateType) {
        case 0:{
            cout<<"Update type: Decrease"<<endl;
//            Graph g1=g;
            runT1 = 0;
            for(int u=0;u<updateBatch;u++){
                ID1=testdata[u].first.first;
                ID2=testdata[u].first.second;
                oldW=testdata[u].second;
                newW=oldW*0.5;
                if(newW < 1) {
                    cout<<"New edge weight is not positive! "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                    exit(1);
                }
                cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                t10=std::chrono::high_resolution_clock::now();
                g.Decrease_CH(ID1,ID2,oldW,newW);
                t11=std::chrono::high_resolution_clock::now();
                time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
                runT1 += time_span1.count();
//                g.CorrectnessCheck_CH(100);
//                testdata[u].second = newW;
            }
            cout<<"Average Decrease update Time: "<<runT1/updateBatch<<" s.\n"<<endl;
            break;
        }
        case 1:{
            cout<<"Update type: Increase"<<endl;
//            Graph g2=g;
            runT1 = 0;
            for(int u=0;u<updateBatch;u++){//46
                ID1=testdata[u].first.first;
                ID2=testdata[u].first.second;
                oldW=testdata[u].second;
                newW=oldW*2;
                cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
//                g.CorrectnessCheck(200);
//                g.DijkstraPath(142488,143850);
                //cout<<"Before update: "<<g.Query(261222,146486)<<" "<<g.Dijkstra(261222,146486,g.Neighbor)<<endl;

                //g.DijkstraPath(261222, 146486);
                //g.WriteCoreGraph(graphfile+"Core");
                t10=std::chrono::high_resolution_clock::now();
                g.Increase_CH(ID1,ID2,oldW,newW);
                t11=std::chrono::high_resolution_clock::now();
//                g.WriteCoreGraph(graphfile+"Core2");
                time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
                runT1 += time_span1.count();
                //cout<<"After update: "<<g.Query(261222,146486)<<" "<<g.Dijkstra(261222,146486,g.Neighbor)<<endl;
//                g.DijkstraPath(63055,233384);
//                g.QueryCoreDebug(142488,143850);
//                g.CorrectnessCheck_CH(100);
            }
            cout<<"Average Increase update Time: "<<runT1/updateBatch<<" s."<<endl;
            break;
        }
        default:
            break;
    }

    cout<<"------------------"<<endl;
	return 0;
}





