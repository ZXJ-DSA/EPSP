/*
 * mainPost.cpp
 *
 *  Created on: 27 Oct 2022
 *      Author: zhangmengxuan, Xinjie ZHOU
 */
#include "head.h"

int main(int argc, char** argv){
    if( argc != 5){
        printf("usage:\n<arg1> source path, e.g. /export/project/xzhouby\n");
        printf("<arg2> dataset, e.g. NY\n");
        printf("<arg3> partition method, e.g. Bubble\n");
        printf("<arg4> partition number, e.g. 4\n");
        exit(0);
    }
    cout<<"This is test for Edge-cut CH (Post-boundary)!"<<endl;

    string graphfile="/media/TraminerData/mengxuan/PartitionData/NY";
    string sourcefile;
    string dataset;
//    string partiName;
    int updateBatch = 100;
    int runtimes = 1000;

	Graph g;
	g.pnum=75;
	g.threadnum=30;

    if(argc > 1){
        sourcefile = argv[1];
        cout<<"argv[1]: "<<argv[1]<<endl;
        dataset = argv[2];
        cout<<"argv[2]: "<<argv[2]<<endl;
        g.partiName = argv[3];
        cout<<"argv[3]: "<<argv[3]<<endl;
        g.pnum = atoi(argv[4]);
        cout<<"argv[4]: "<<argv[4]<<endl;
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

	std::chrono::high_resolution_clock::time_point t10;
	std::chrono::high_resolution_clock::time_point t11;
	std::chrono::duration<double> time_span1;
	double runT1=0,maxrunT=0;
    double runT2=0;
    double runT3=0;
    double runT4=0;
    double maxrunT2=0;
	//Index construction of subgraphs
	g.NeighborCons.clear();
	g.SCconNodesMTs.clear();
	//cout<<"partition ID ";
	for(int partiID=0;partiID<g.pnum;partiID++){
		//cout<<partiID<<" ";
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
		if(time_span1.count()>maxrunT)
			maxrunT=time_span1.count();
		g.NeighborCons.push_back(h.NeighborCon);
		g.SCconNodesMTs.push_back(h.SCconNodesMT);
	}
	//cout<<endl;
	cout<<"Index construction step1: Partitions' index construction time "<<runT1<<", maxrunT "<<maxrunT<<endl;

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
	cout<<"Index construction step2: Overlay Index Construction Time: "<<runT1<<" s."<<endl;
	g.NeighborConOverlay=l.NeighborCon;
	g.SCconNodesMTOverlay=l.SCconNodesMT;

	//index rebuild of subgraph
	g.NeighborsPartiPost=g.NeighborsParti;
	g.PartitionUpdate();
	for(int partiID=0;partiID<g.pnum;partiID++){
		//cout<<partiID<<" ";
		Graph h;
		h.nodenum=g.nodenum;
		h.NodeOrder=g.NodeOrder;
		h.threadnum=g.threadnum;
		h.vNodeOrder=g.vNodeOrder;
		h.Neighbor=g.NeighborsPartiPost[partiID];

		t10=std::chrono::high_resolution_clock::now();
		h.CHindexMT();
		t11=std::chrono::high_resolution_clock::now();
		time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
		runT4+= time_span1.count();
		if(time_span1.count()>maxrunT2)
			maxrunT2=time_span1.count();
		g.NeighborConsPost.push_back(h.NeighborCon);
	}
	//cout<<endl;
	cout<<"Index construction step3: Partitions Index Reconstruction Time: "<<runT4<<" s, maxrunT: "<<maxrunT2<<" s."<<endl;

    cout<<"Overall Index Construction Time: "<<runT1+runT2+runT3+runT4<<" s."<<endl;
    cout<<"Minimum Overall Index Construction Time: "<<maxrunT+runT3+maxrunT2<<" s."<<endl;

	//Index size
	g.IndexsizePost();

	//Query processing efficiency
	g.EffiCheckPost(graphfile+".query", runtimes);
	g.EffiCheckPost1();


	return 0;
}


