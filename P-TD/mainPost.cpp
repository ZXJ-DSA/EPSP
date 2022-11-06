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
        printf("<arg3> partition method, e.g. NC\n");
        printf("<arg4> partition number, e.g. 4\n");
        exit(0);
    }
    cout<<"This is test for Edge-cut PLL (Post-boundary)!"<<endl;

    string graphfile="/media/TraminerData/mengxuan/PartitionData/NY";
    string sourcefile;
    string dataset;
//    string partiName;
    int updateBatch = 100;
    int runtimes = 1000;

	Graph g;
	g.pnum=75;
	g.threadnum=150;

    if(argc > 1){
        sourcefile = argv[1];
        cout<<"argv[1]: "<<argv[1]<<endl;
        dataset = argv[2];
        cout<<"argv[2]: "<<argv[2]<<endl;
        g.algoName = argv[3];
        cout<<"argv[3]: "<<argv[3]<<endl;
        g.pnum = atoi(argv[4]);
        cout<<"argv[4]: "<<argv[4]<<endl;
    }
    cout<<"Dataset: "<<dataset<<endl;
    cout<<"Partition method: "<<g.algoName<<endl;
    cout<<"Thread number: "<<g.threadnum<<endl;
    graphfile = sourcefile + "/"+dataset+"/"+dataset;

    string partitionPath = graphfile+"_"+g.algoName + "_"+to_string(g.pnum);

//    g.GraphRead(graphfile+"/Graph");
//    g.OrderRead(graphfile+"/Order75");
    g.GraphRead(graphfile);
    g.OrderRead(partitionPath+"/vertex_order");

	//Read partitioned subgraphs
//	g.GraphPartitionRead(graphfile+"/partition"+to_string(g.pnum));
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
	g.Trees.clear();
	g.toRMQs.clear();
	g.RMQIndexs.clear();
	g.ranks.clear();
	g.heightMaxs.clear();
	//g.SCconNodess.clear();
	g.SCconNodesMTs.clear();
	g.VidtoTNids.clear();
	//cout<<"partition ID ";
	for(int partiID=0;partiID<g.pnum;partiID++){
		//cout<<partiID<<" ";
		Graph h;
		h.nodenum=g.nodenum;
		h.threadnum=g.threadnum;
		h.NodeOrder=g.NodeOrder;
		h.vNodeOrder=g.vNodeOrder;
		h.Neighbors=g.NeighborsParti[partiID];

		t10=std::chrono::high_resolution_clock::now();
		h.H2Hindex();
		t11=std::chrono::high_resolution_clock::now();
		time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
		runT1+= time_span1.count();
		if(time_span1.count()>maxrunT)
			maxrunT=time_span1.count();
		g.Trees.push_back(h.Tree);
		g.toRMQs.push_back(h.toRMQ);
		g.RMQIndexs.push_back(h.RMQIndex);
		g.ranks.push_back(h.rank);
		g.heightMaxs.push_back(h.heightMax);
		//g.SCconNodess.push_back(h.SCconNodes);
		g.SCconNodesMTs.push_back(h.SCconNodesMT);
		g.VidtoTNids.push_back(h.VidtoTNid);
	}
	//cout<<endl;
	cout<<"Index construction step1: Partitions' index construction time: "<<runT1<<" s, maxrunT: "<<maxrunT<<" s."<<endl;

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
	l.Neighbors=g.NeighborsOverlay;
	l.nodenum=g.nodenum;
	l.NodeOrder=g.NodeOrder;
	l.vNodeOrder=g.vNodeOrder;
	l.threadnum=g.threadnum;

	t10=std::chrono::high_resolution_clock::now();
	l.H2Hindex();
	t11=std::chrono::high_resolution_clock::now();
	time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
	runT3= time_span1.count();
	cout<<"Index construction step2: Overlay Index Construction Time: "<<runT3<<" s."<<endl;
	g.TreeOverlay=l.Tree;
	g.toRMQOverlay=l.toRMQ;
	g.RMQIndexOverlay=l.RMQIndex;
	g.rankOverlay=l.rank;
	g.heightMaxOverlay=l.heightMax;
	//g.SCconNodesOverlay=l.SCconNodes;
	g.SCconNodesOverlayMT=l.SCconNodesMT;
	g.VidtoTNidOverlay=l.VidtoTNid;

	//index rebuild of subgraphs
	g.NeighborsPartiPost=g.NeighborsParti;
	g.PartitionUpdate();
	g.Trees.clear();
	g.toRMQs.clear();
	g.RMQIndexs.clear();
	g.ranks.clear();
	g.heightMaxs.clear();
	//g.SCconNodess.clear();
	g.SCconNodesMTs.clear();
	g.VidtoTNids.clear();
	//cout<<"partition ID ";
	for(int partiID=0;partiID<g.pnum;partiID++){
		//cout<<partiID<<" ";
		Graph h;
		h.nodenum=g.nodenum;
		h.threadnum=g.threadnum;
		h.NodeOrder=g.NodeOrder;
		h.vNodeOrder=g.vNodeOrder;
		h.Neighbors=g.NeighborsPartiPost[partiID];

		t10=std::chrono::high_resolution_clock::now();
		h.H2Hindex();
		t11=std::chrono::high_resolution_clock::now();
		time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
		runT4+= time_span1.count();
		if(time_span1.count()>maxrunT2)
			maxrunT2=time_span1.count();
		g.Trees.push_back(h.Tree);
		g.toRMQs.push_back(h.toRMQ);
		g.RMQIndexs.push_back(h.RMQIndex);
		g.ranks.push_back(h.rank);
		g.heightMaxs.push_back(h.heightMax);
		//g.SCconNodess.push_back(h.SCconNodes);
		g.SCconNodesMTs.push_back(h.SCconNodesMT);
		g.VidtoTNids.push_back(h.VidtoTNid);
	}
	//cout<<endl;
	cout<<"Index construction step3: Partitions Index Reconstruction Time: "<<runT4<<" s, maxrunT: "<<maxrunT2<<" s."<<endl;

    cout<<"Overall Index Construction Time: "<<runT1+runT2+runT3+runT4<<" s."<<endl;
    cout<<"Minimum Overall Index Construction Time: "<<maxrunT+runT3+maxrunT2<<" s."<<endl;

    //Index size
    g.IndexsizePost();

	//Query processing efficiency
//	g.EffiCheckPost(graphfile+"/ODeffi");
    cout<<"Query processing test.."<<endl;
    g.EffiCheckPost(graphfile+".query", runtimes);

	g.EffiCheckPost1();

	cout<<"---------------"<<endl;

	return 0;
}



