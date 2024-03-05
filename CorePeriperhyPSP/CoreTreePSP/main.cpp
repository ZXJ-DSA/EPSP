/*
 * main.cpp
 *
 *  Created on: 14 Oct 2022
 *      Author: Xinjie ZHOU
 */
#include "head.h"

int main(int argc, char** argv){

    if( argc < 4 || argc > 10){//
        printf("usage:\n<arg1> source path, e.g. /export/project/xzhouby\n");
        printf("<arg2> name of dataset, e.g. NY\n");
        printf("<arg3> tree width, e.g. 20\n");
        printf("<arg4> (optional) Tree index type, (0: CH; 1: TD), default: 1\n");
        printf("<arg5> (optional) update type, (0: No Update Test; 1: Decrease; 2: Increase), default: 0\n");
        printf("<arg6> (optional) update number, eg. 1000\n");
        printf("<arg7> (optional) thread number, eg. 15\n");
        printf("<arg8> (optional) PSP strategy, (1: Pre-boundary; 2: No-boundary; 3: Post-boundary), default: 2\n");
        printf("<arg9> (optional) preprocessing task, (1: Same-tree query generation), default: 1\n");
        exit(0);
    }

    string DesFile="./data/";
    string dataset = "NY";
    int treeWidth = 10;
    int algoCoreC = 0;//0: BPCL; 1: PCL; 2: PLL; 3: WPSL; 4: GLL; 5: Read
    int algoCoreU = 0;//0: PDPLL; 1: SDPLL
    int algoTree = 1;//0: CH; 1: TD
    int updateType = 0;
    int runtimes = 10000;
//    runtimes = 100;
    int updateBatch = 10;
    updateBatch = 1000;
//    updateBatch = 50;
    int threadNum = 15;
    int batchSize = 15;
    int strategy = 2;//PSP strategy
    int preTask = 0;

    if(argc > 1) {
        cout << "argc: " << argc << endl;
        cout << "argv[1]: " << argv[1] << endl;//source path
        DesFile = argv[1];
        if(argc > 2){
            cout << "argv[2]: " << argv[2] << endl;//dataset
            dataset = argv[2];
        }
        if(argc > 3){
            cout << "argv[3]: " << argv[3] << endl;//treewidth
            treeWidth = stoi(argv[3]);
        }
        if(argc > 4){
            cout << "argv[4]: " << argv[4] << endl;//tree index type
            algoTree = stoi(argv[4]);
        }
        if(argc > 5){
            cout << "argv[5]: " << argv[5] << endl;//update type
            updateType = stoi(argv[5]);
        }
        if(argc > 6){
            cout << "argv[6]: " << argv[6] << endl;//update number
            updateBatch = stoi(argv[6]);
        }
        if(argc > 7){
            cout << "argv[7]: " << argv[7] << endl;//thread number
            threadNum = stoi(argv[7]);
        }
        if(argc > 8){
            cout << "argv[8]: " << argv[8] << endl;//PSP strategy
            strategy = stoi(argv[8]);
        }
        if(argc > 9){
            cout << "argv[9]: " << argv[9] << endl;//Preprocessing task
            preTask = stoi(argv[9]);
            if(preTask!=1){
                cout<<"Wrong preprocessing task!"<<preTask<<endl; exit(1);
            }
        }
    }


	//used for running time calculation
    Timer tt0;
    tt0.start();

    string sourcePath=DesFile+"/"+dataset+"/";
//    string orderfile=graphfile+".order";
    string ODfile=sourcePath+dataset+".query";


    string updateFile=sourcePath+dataset+".update";

    Graph g;
    g.threadnum=threadNum;//thread number of parallel computation (can be changed)
    if(batchSize>=threadNum){
        g.batchsize=batchSize;
    }else{
        g.batchsize=threadNum;
    }

    g.ifParallel = true;
    g.bandWidth=treeWidth;//15;
    g.dataset=dataset;
    g.strategy=strategy;
    g.algoCoreC=algoCoreC;
    g.algoCoreU=algoCoreU;
    g.algoTree=algoTree;
    cout<<"Dataset: "<<dataset<<endl;
    cout<<"Bandwidth: "<<g.bandWidth<<endl;
    cout<<"Thread number: "<<g.threadnum<<endl;
    cout<<"Batch size: "<<g.batchsize<<endl;
    cout<<"PSP strategy: "<<g.strategy<<endl;
    cout<<"Tree index: ";
    if(g.algoTree==0){
        cout<<"CH"<<endl;
    }else if(g.algoTree==1){
        cout<<"H2H"<<endl;
    }else{
        cout<<"wrong tree index type!"<<endl; exit(1);
    }



    if(g.algoCoreU==0){
        cout<<"Core Update algorithm: Propagation-based DPLL."<<endl;
    }else if(g.algoCoreU==1){
        cout<<"Core Update algorithm: Search-based DPLL."<<endl;
    }

//    g.CoreGraphDebug(graphfile+"Core3");//Test
//    g.CoreGraphDebug(graphfile+"CorePLL_148");//NY, original

    g.sourcePath=sourcePath;

    g.ReadGraph(sourcePath+dataset);//
//    g.StainingMethod(0);

    if(preTask==1){
        g.QueryGenerationSameParti();//same-partition query generation
    }

    ///Task 1: Index construction
    g.IndexConstruction();
//    g.WriteCTIndex(graphfile);

//    g.CTQueryUpdateGen(10000);

    ///Task 2: Query processing
    g.CorrectnessCheck(100);
    g.EffiCheck(ODfile,runtimes);//query efficiency test
    g.EffiCheck(sourcePath+"tmp/"+dataset+"_sameParti_CT"+ to_string(treeWidth)+".query",runtimes);//query efficiency test
//    exit(0);
//#ifdef __APPLE__
////    cout<<"The platform is macOS."<<endl;
//    g.EffiCheck(ODfile+"SamePartiPlanar",runtimes);//query efficiency test
//#else
//    g.EffiCheck("/home/data/xzhouby/datasets/"+g.dataset+"/"+g.dataset+".querySamePartiCoreTree",runtimes);//query efficiency test
//#endif

//    g.SameTreeQueryTest(ODfile,runtimes);
//    exit(0);
    ///Task 3: Index update
    g.IndexMaintenance(updateFile,updateType,updateBatch);//index maintenance
//    g.IndexMaintenanceTypeTest(updateType,updateBatch);
//    g.IndexMaintenance(updateFile+"2",updateType,updateBatch);//index maintenance
//    g.IndexMaintenance(updateFile+"ST",updateType,updateBatch);//same-tree index maintenance

    tt0.stop();
    cout<<"\nOverall runtime: "<<tt0.GetRuntime()<<" s."<<endl;
    cout<<"------------------\n"<<endl;
    exit(0);
    g.clear();
	return 0;
}
