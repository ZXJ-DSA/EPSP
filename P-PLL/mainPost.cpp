/*
 * mainPost.cpp
 *
 *  Created on: 26 Oct 2022
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

	//Index construction step1: subgraphs
	g.LabelParti.clear();
	g.PruningPointParti.clear();
	set<pair<int,int>> setpair;
	setpair.clear();
	g.NoSupportedParti.assign(g.pnum,setpair);
	g.NoSupportedPartiPost.assign(g.pnum,setpair);
	g.NoSupportedOverlay.clear();
	for(int partiID=0;partiID<g.pnum;partiID++){
        if(partiID%20==0)
            cout<<"partition ID "<<partiID<<endl;
		Graph h;
		h.nodenum=g.nodenum;
		h.NodeOrder=g.NodeOrder;
		h.vNodeOrder=g.vNodeOrder;
		h.Neighbors=g.NeighborsParti[partiID];
		h.threadnum=g.threadnum;

		t10=std::chrono::high_resolution_clock::now();
		//***************PSL construction**************
		h.vSm.reserve(h.nodenum);
		for(int i = 0; i < h.nodenum; i++)
		{
			Semaphore* s = new Semaphore(1);
			h.vSm.push_back(s);
		}
		h.IndexConstructMThread2New();
		//***************PLL construction**************
		//h.IndexConst();
		t11=std::chrono::high_resolution_clock::now();
		time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
		runT1+= time_span1.count();
		if(time_span1.count()>maxrunT)
			maxrunT=time_span1.count();
		g.LabelParti.push_back(h.Label);
		g.PruningPointParti.push_back(h.PruningPointNew);
	}
	cout<<"Index construction step1: Partitions' index construction time: "<<runT1<<" s, maxrunT: "<<maxrunT<<" s."<<endl;

	//Build overlay graph
    t10=std::chrono::high_resolution_clock::now();
    g.OverlayGraphConstructPost();
    t11=std::chrono::high_resolution_clock::now();
    time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
    runT2 = time_span1.count();
    //g.CorrectnessCheckOverlay();
    cout<<"Overlay graph construction time: "<< runT2<<" s."<<endl;

	//Index construction step2: overlay graph
	Graph l;
	l.Neighbors=g.NeighborsOverlay;
	l.nodenum=g.nodenum;
	l.NodeOrder=g.NodeOrder;
	l.vNodeOrder=g.vNodeOrder;
	l.threadnum=g.threadnum;
	t10=std::chrono::high_resolution_clock::now();
	//***************PSL construction**************
    cout<<"PSL"<<endl;
	l.vSm.reserve(l.nodenum);
	for(int i = 0; i < l.nodenum; i++)
	{
		Semaphore* s = new Semaphore(1);
		l.vSm.push_back(s);
	}
	l.IndexConstructMThread2New();
	//***************PLL construction**************
	//l.IndexConst();
	t11=std::chrono::high_resolution_clock::now();
	time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
	runT3 = time_span1.count();
	cout<<"Index construction step2: Overlay Index Construction Time: "<<runT3<<" s."<<endl;
	g.LabelOverlay=l.Label;
	g.PruningPointOverlay=l.PruningPointNew;

	//Index construction step3: update subgraph
	g.NeighborsPartiPost=g.NeighborsParti;
	g.LabelPartiPost=g.LabelParti;
	g.PruningPointPartiPost=g.PruningPointParti;
	t10=std::chrono::high_resolution_clock::now();
	g.PartitionUpdate();
	t11=std::chrono::high_resolution_clock::now();
	time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
	runT4 = time_span1.count();
	cout<<"Index construction step3: Partition's index update time: "<<runT4<<" s."<<endl;
    cout<<"Overall Index Construction Time: "<<runT1+runT2+runT3+runT4<<" s."<<endl;
    cout<<"Minimum Overall Index Construction Time: "<<maxrunT+runT3+runT4<<" s."<<endl;

    g.IndexsizePost();

	//Query processing efficiency
//	g.EffiCheckPost(graphfile+"/ODeffi");
    cout<<"Query processing test.."<<endl;
    g.EffiCheckPost(graphfile+".query", runtimes);

	//Query efficiency for different query distribution
	g.EffiCheckPost1();

	//cout<<"The update time of post is the sum of No-Boundary Update Time + Index construction step3"<<endl;

	//Index update
	//int ID1, ID2, oldW, newW;
	//srand (time(NULL));

	/*Graph g2=g;
	t10=std::chrono::high_resolution_clock::now();
	for(int u=0;u<100;u++){
		ID1=rand()%g2.nodenum;
		ID2=g2.Neighbors[ID1][0].first;
		oldW=g2.Neighbors[ID1][0].second;
		newW=oldW*0.5;
		g2.DecreasePost(ID1,ID2,oldW,newW);
	}
	t11=std::chrono::high_resolution_clock::now();
	time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
	runT1= time_span1.count();
	cout<<"Decrease update Time "<<runT1/100<<endl;

	Graph g1=g;
	t10=std::chrono::high_resolution_clock::now();
	for(int u=0;u<100;u++){
		ID1=rand()%g1.nodenum;
		ID2=g1.Neighbors[ID1][0].first;
		oldW=g1.Neighbors[ID1][0].second;
		newW=oldW*1.5;
		g1.IncreasePost(ID1,ID2,oldW,newW);
	}
	t11=std::chrono::high_resolution_clock::now();
	time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t11-t10);
	runT1= time_span1.count();
	cout<<"Increase update Time "<<runT1/100<<endl;*/

	return 0;
}



