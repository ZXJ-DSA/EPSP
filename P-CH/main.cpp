/*
 * main.cpp
 *
 *  Created on: 19 Oct 2022
 *      Author: zhangmengxuan, Xinjie ZHOU
 */
#include "head.h"

int main(int argc, char** argv){
    if( argc != 5 && argc !=6){
        printf("usage:\n<arg1> source path, e.g. /export/project/xzhouby\n");
        printf("<arg2> dataset, e.g. NY\n");
        printf("<arg3> partition method, e.g. NC\n");
        printf("<arg4> partition number, e.g. 64\n");
        printf("<arg5> (optional) update type, e.g. 0: Decrease; 1: Increase.\n");
        exit(0);
    }
    cout<<"This is test for Edge-cut CH!"<<endl;

	string graphfile="/media/TraminerData/mengxuan/PartitionData/NY";
    string sourcefile;
    string dataset;
//    string partiName;
    int updateType = 0;
    int updateBatch = 100;
    int runtimes = 1000;

	Graph g;
    g.pnum=75;
    g.threadnum=10;

    if(argc > 1){
        sourcefile = argv[1];
        cout<<"argv[1]: "<<argv[1]<<endl;
        dataset = argv[2];
        cout<<"argv[2]: "<<argv[2]<<endl;
        g.partiName = argv[3];
        cout<<"argv[3]: "<<argv[3]<<endl;
        g.pnum = atoi(argv[4]);
        cout<<"argv[4]: "<<argv[4]<<endl;
        if(argc > 5){
            updateType = atoi(argv[5]);
            cout<<"argv[5]: "<<argv[5]<<endl;
        }
    }
    cout<<"Dataset: "<<dataset<<endl;
    cout<<"Partition method: "<<g.partiName<<endl;
    cout<<"Thread number: "<<g.threadnum<<endl;

    graphfile = sourcefile + "/"+dataset+"/"+dataset;

    string partitionPath = graphfile+"_"+g.partiName + "_"+to_string(g.pnum);

    g.GraphRead(graphfile);
    g.OrderRead(partitionPath+"/vertex_order");


	//Read partitioned subgraphs
    g.GraphPartitionRead(partitionPath);

    /// Task 1: Index construction
	std::chrono::high_resolution_clock::time_point t10;
	std::chrono::high_resolution_clock::time_point t11;
	std::chrono::duration<double> time_span1;
	double runT1=0;
    double runT2=0;
    double runT3=0;
    double maxrT=0;
	//Index construction of subgraphs
	g.NeighborCons.clear();
	g.SCconNodesMTs.clear();
//	cout<<"partition ID ";
	for(int partiID=0;partiID<g.pnum;partiID++){
//		cout<<partiID<<" ";
		Graph h;
		h.nodenum=g.nodenum;
		h.NodeOrder=g.NodeOrder;
		h.threadnum=g.threadnum;
		h.vNodeOrder=g.vNodeOrder;
		h.Neighbor=g.NeighborsParti[partiID];

		t10=std::chrono::high_resolution_clock::now();
		h.CHindexMT();
		t11=std::chrono::high_resolution_clock::now();
		time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
		runT1+= time_span1.count();
        if(maxrT < time_span1.count()){
            maxrT = time_span1.count();
        }
		g.NeighborCons.push_back(h.NeighborCon);
		g.SCconNodesMTs.push_back(h.SCconNodesMT);
	}
//	cout<<endl;
    cout<<"Partitions' index construction time: "<<runT1<< " s. maxrT: "<<maxrT<<" s."<<endl;

	//Build overlay graph
    t10=std::chrono::high_resolution_clock::now();
	g.OverlayGraphConstructPost();
    t11=std::chrono::high_resolution_clock::now();
    time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
    runT2 = time_span1.count();
	//g.CorrectnessCheckOverlay();
    cout<<"Overlay graph construction time: "<< runT2<<" s."<<endl;

	//Index construction of overlay graph
	Graph l;
	l.Neighbor=g.NeighborsOverlay;
	l.nodenum=g.nodenum;
	l.threadnum=g.threadnum;
	l.NodeOrder=g.NodeOrder;
	l.vNodeOrder=g.vNodeOrder;

	t10=std::chrono::high_resolution_clock::now();
	l.CHindexMT();
	t11=std::chrono::high_resolution_clock::now();
	time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
    runT3= time_span1.count();
    cout<<"Overall Index Construction Time: "<<runT1+runT2+runT3<<" s."<<endl;
    cout<<"Minimum Overall Index Construction Time: "<<maxrT+runT2+runT3<<" s."<<endl;
	g.NeighborConOverlay=l.NeighborCon;
	g.SCconNodesMTOverlay=l.SCconNodesMT;

//	g.CorrectnessCheck();
	//Index size
	g.Indexsize();

    ///Task2: Efficiency test
	//Query processing efficiency
//	g.EffiCheck(graphfile+"/ODeffi");
    g.EffiCheck(graphfile+".query", runtimes);

    g.EffiCheck();

//    exit(0);

    ///Task 3: Update test
	//Index update
    string updateFile = graphfile + ".update";
    vector<pair<pair<int,int>,int>> testdata;
    g.UpdateRead(updateFile, testdata);
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
                g1.DecreaseSingle(ID1,ID2,oldW,newW);
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
                g.IncreaseSingle(ID1,ID2,oldW,newW);
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



    cout<<"----------------\n"<<endl;
	return 0;
}



