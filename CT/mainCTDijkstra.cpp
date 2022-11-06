/*
 * main.cpp
 *
 *  Created on: 13 Oct 2022
 *      Author: zhangmengxuan
 */
#include "head.h"

int main(int argc, char** argv){

    if( argc != 4 && argc != 5){// && argc != 6
        printf("usage:\n<arg1> source path, e.g. /export/project/xzhouby\n");
        printf("<arg2> name of dataset, e.g. NY\n");
        printf("<arg3> tree width, e.g. 20\n");
        printf("<arg4> (optional) update type, (0: Decrease; 1: Increase) e.g. 0\n");
        exit(0);
    }
    cout<<"This is test for CT-DS!"<<endl;
    string DesFile="./data/";
    string dataset = "NY";
    int treeWidth = 20;
    int updateType = 0;

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
    }

    int runtimes = 1000;
    int updateBatch = 10;

	//used for running time calculation
	std::chrono::high_resolution_clock::time_point t10;
	std::chrono::high_resolution_clock::time_point t11;
	std::chrono::duration<double> time_span1;
	double runT1=0;

//    string graphfile="/media/TraminerData/mengxuan/MengxuanGraphWPSL/Cond/CondWeighted";
    string graphfile=DesFile+"/"+dataset+"/"+dataset;
    string orderfile=graphfile+".order";
    string ODfile=graphfile+".query";
    string updateFile=graphfile+".update";

    Graph g;
    g.threadnum=10;//thread number of parallel computation (can be changed)
    g.Width=treeWidth;//20;
    cout<<"Dataset: "<<dataset<<endl;
    cout<<"Tree width: "<<g.Width<<endl;
    cout<<"Thread number: "<<g.threadnum<<endl;

    g.ReadGraph(graphfile);
    g.StainingMethod(0);

    g.PartitionPreProcess();
    g.PartitionPostProcess();

    ///Task 1: construct the index of Core
    cout<<"Begin Core's Index Construction.................."<<endl;
    Graph h;
    h.nodenum=g.nodenum;
    h.Neighbors=g.AdjaCore;
    h.NodeOrder=g.NodeOrder;
    h.vNodeOrder=g.vNodeOrder;
    h.threadnum=g.threadnum;
    t10=std::chrono::high_resolution_clock::now();
    //****************PSL construction***************************
    cout<<"PSL"<<endl;
    h.vSm.reserve(h.nodenum);
    for(int i = 0; i < h.nodenum; i++)
    {
        Semaphore* s = new Semaphore(1);
        h.vSm.push_back(s);
    }
    h.IndexConstructMThread2New();
    //***********************************************************

    //****************PLL construction***************************
//    cout<<"PLL"<<endl;
    //h.IndexConst1();
    //***********************************************************
    t11=std::chrono::high_resolution_clock::now();
    time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
    runT1 = time_span1.count();
    g.IndexCore=h.Label;
    g.PruningPointCore=h.PruningPointNew;
    g.NoSupportedPairCore=h.NoSupportedPair;
    cout<<"Core's Index Construction Finished"<<endl;
    cout<<"Construction Time: "<<runT1<<" s."<<endl;

    g.indexsizeCTDijk();//index (labeling+pruning point) size computation

    //g.CorrectnessCheck();
    ///Task 2: Query processing
    g.EffiCheck(ODfile,runtimes);//query efficiency test

    ///Task 3: Index update
    vector<pair<pair<int,int>,int>> testdata;
    g.ReadUpdate(updateFile, testdata);
    cout<<"Update Batch: "<<updateBatch<<endl;

    //index update
    int ID1, ID2, oldW, newW;
    srand (time(NULL));
    switch (updateType) {
        case 0:{
            //Decrease update
            cout<<"Update type: Decrease"<<endl;
            Graph g1=g;
            t10=std::chrono::high_resolution_clock::now();
            for(int u=0;u<updateBatch;u++){
                ID1=testdata[u].first.first;
                ID2=testdata[u].first.second;
                oldW=testdata[u].second;
                newW=oldW*0.5;
                if(newW < 1){
                    cout<<"Edge weight wrong! "<<newW<<endl;
                    exit(1);
                }
                g1.Decrease(ID1,ID2,oldW,newW);
            }
            t11=std::chrono::high_resolution_clock::now();
            time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
            runT1= time_span1.count();
            cout<<"Average Decrease update Time: "<<runT1/updateBatch<<" s."<<endl;
//            break;
        }
        case 1:{
            cout<<"Update type: Increase"<<endl;
//            Graph g2=g;
            t10=std::chrono::high_resolution_clock::now();
            for(int u=0;u<updateBatch;u++){
                ID1=testdata[u].first.first;
                ID2=testdata[u].first.second;
                oldW=testdata[u].second;
                newW=oldW*2;
                g.Increase(ID1,ID2,oldW,newW);
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

    cout<<"------------------"<<endl;
	return 0;
}



