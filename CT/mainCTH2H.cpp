/*
 * main.cpp
 *
 *  Created on: 14 Oct 2022
 *      Author: zhangmengxuan, Xinjie ZHOU
 */
#include "head.h"

bool ifExtensiion = true;

int main(int argc, char** argv){

    if( argc != 4 && argc != 5){// && argc != 6
        printf("usage:\n<arg1> source path, e.g. /export/project/xzhouby\n");
        printf("<arg2> name of dataset, e.g. NY\n");
        printf("<arg3> tree width, e.g. 20\n");
        printf("<arg4> (optional) update type, (0: Decrease; 1: Increase) e.g. 0\n");
        exit(0);
    }
    cout<<"This is test for CT-TD!"<<endl;
    string DesFile="./data/";
    string dataset = "NY";
    int treeWidth = 20;
    int updateType =0;

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
    double runT2=0;
    double runT3=0;
    double runT4=0;
    double maxrT=0;

//    string graphfile="/media/TraminerData/mengxuan/MengxuanGraphWPSL/Cond/CondWeighted";
    string graphfile=DesFile+"/"+dataset+"/"+dataset;
    string orderfile=graphfile+".order";
    string ODfile=graphfile+".query";
    string updateFile=graphfile+".update";

    Graph g;
    g.threadnum=10;//thread number of parallel computation (can be changed)
    g.ifParallel = true;
    g.Width=treeWidth;//15;
    cout<<"Dataset: "<<dataset<<endl;
    cout<<"Tree width: "<<g.Width<<endl;
    cout<<"Thread number: "<<g.threadnum<<endl;

    g.ReadGraph(graphfile);
    g.StainingMethod(0);

    g.PartitionPreProcess();
    g.PartitionPostProcess();

    ///Task 1: Index construction
    //construct index for peripheries
    g.Trees.clear();
    g.toRMQs.clear();
    g.RMQIndexs.clear();
    g.ranks.clear();
    g.heightMaxs.clear();
    g.SCconNodesMTs.clear();
    g.VidtoTNids.clear();
//    cout<<"partition ID ";
    int ID;
    cout<<"Periphery index construction..."<<endl;

    for(int partiID=0;partiID<g.partiNum;partiID++){
        if(partiID%1000 == 0)
            cout<<"partiID: "<<partiID<<endl;
        Graph h;
        h.nodenum=g.nodenum;
        h.NodeOrder=g.NodeOrder;
        h.vNodeOrder=g.vNodeOrder;
        h.Neighbors=g.AdjaParti[partiID];
        //h.Neighbors.assign(h.nodenum, vector<pair<int,int>>());

        t10=std::chrono::high_resolution_clock::now();
        h.H2Hindex();
        t11=std::chrono::high_resolution_clock::now();
        time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
        runT1+= time_span1.count();
        if(maxrT < time_span1.count())
            maxrT = time_span1.count();
        g.Trees.push_back(h.Tree);
        g.toRMQs.push_back(h.toRMQ);
        g.RMQIndexs.push_back(h.RMQIndex);
        g.ranks.push_back(h.rank);
        g.heightMaxs.push_back(h.heightMax);
        g.SCconNodesMTs.push_back(h.SCconNodesMT);
        g.VidtoTNids.push_back(h.VidtoTNid);
    }


    cout<<"Partition index construction time: "<<runT1<<" s. maxrT: "<<maxrT<<" s."<<endl;

    t10=std::chrono::high_resolution_clock::now();
    g.OverlayGraph();//finish the overlay graph construction
    t11=std::chrono::high_resolution_clock::now();
    time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
    runT2 = time_span1.count();
    cout<<"Overlay graph construction time: "<<runT2<<" s."<< endl;
    //construct the index of Core
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
    runT3 = time_span1.count();
    g.IndexCore=h.Label;
    g.PruningPointCore=h.PruningPointNew;
    g.NoSupportedPairCore=h.NoSupportedPair;
//    cout<<"Core's Index Construction Finished"<<endl;
    cout<<"Core's construction time: "<<runT3<<" s."<< endl;
    if(ifExtensiion){
        cout<<"Begin Extension Index Construction..."<<endl;
        t10=std::chrono::high_resolution_clock::now();
        g.ExtensionIndexConstruct(g.ifParallel);
        t11=std::chrono::high_resolution_clock::now();
        time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
        runT4 = time_span1.count();
        cout<<"Extension index construction time: "<<runT4<<" s."<<endl;
    }
    cout<<"Overall Construction Time: "<<runT1+runT2+runT3+runT4<<" s."<<endl;
    cout<<"Minimum Overall Construction Time: "<<runT3+runT4+maxrT<<" s."<<endl;

    g.indexsizeCTH2H();//index (labeling+pruning point) size computation
    //g.CorrectnessCheck();

    ///Task 2: Query processing
//    g.EffiCheck("/media/TraminerData/mengxuan/MengxuanGraphWPSL/Cond/ODpair");//query efficiency test
    g.EffiCheck(ODfile,runtimes);//query efficiency test

    ///Task 3: Index update
    vector<pair<pair<int,int>,int>> testdata;
    g.ReadUpdate(updateFile, testdata);

    cout<<"Update Batch: "<<updateBatch<<endl;
    int ID1, ID2, oldW, newW;
    srand (time(NULL));
    switch (updateType) {
        case 0:{
            cout<<"Update type: Decrease"<<endl;
            Graph g1=g;
            t10=std::chrono::high_resolution_clock::now();
            for(int u=0;u<updateBatch;u++){
                ID1=testdata[u].first.first;
                ID2=testdata[u].first.second;
                oldW=testdata[u].second;
                newW=oldW*0.5;
                if(newW < 1) {
                    cout<<"New edge weight is not positive! "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                    exit(1);
                }
//                cout<<u<<" "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                g1.Decrease(ID1,ID2,oldW,newW);
//                g1.CorrectnessCheck();
            }
            t11=std::chrono::high_resolution_clock::now();
            time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
            runT1= time_span1.count();
            cout<<"Average Decrease update Time: "<<runT1/updateBatch<<" s."<<endl;
//            break;
        }
        case 1:{
            cout<<"Update type: Increase"<<endl;
            //Graph g2=g;
            t10=std::chrono::high_resolution_clock::now();
            for(int u=0;u<updateBatch;u++){
                ID1=testdata[u].first.first;
                ID2=testdata[u].first.second;
                oldW=testdata[u].second;
                newW=oldW*2;
                g.Increase(ID1,ID2,oldW,newW);
        //        g.CorrectnessCheck();
            }
            t11=std::chrono::high_resolution_clock::now();
            time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
            runT1= time_span1.count();
            cout<<"Average Increase update Time: "<<runT1/updateBatch<<" s."<<endl;
        }
    }

    cout<<"------------------"<<endl;
	return 0;
}
