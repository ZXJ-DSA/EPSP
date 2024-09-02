/*
 * main.cpp
 *
 *  Created on: 14 Oct 2022
 *      Author: Xinjie ZHOU
 */
#include "headPSP.h"

int main(int argc, char** argv){

    if( argc < 4 || argc > 15){//
        printf("Planar Dynamic PSP Test.\nusage:\n");
        printf("<arg1> source path, e.g. /export/project/xzhouby\n");
        printf("<arg2> name of dataset, e.g. NY\n");
        printf("<arg3> shortest path index, 1: CH; 2: TD; 3: PLL. default: 2\n");
        printf("<arg4> (optional) partitioned shortest path strategy, 1: Pre-boundary; 2: No-boundary; 3: Post-boundary. default: 2\n");
        printf("<arg5> (optional) partition method, (NC: PUNCH; MT: METIS; SC: SCOTCH; kahypar: KaHyPar; geometric: RCB; Bubble: Bubble; HEP: HEP; CLUGP: CLUGP), default: NC\n");
        printf("<arg6> (optional) partition number, e.g. 64\n");
        printf("<arg7> (optional) update type, (0: No Update Test; 1: Decrease; 2: Increase; 3: Edge insert; 4: Vertex insert; 5: Edge delete; 6: Vertex delete), default: 0\n");
        printf("<arg8> (optional) update batch number, default: 100\n");
        printf("<arg9> (optional) update batch volume, default: 1\n");
        printf("<arg10> (optional) edge weight change percentage, default: 0.5\n");
        printf("<arg11> (optional) overlay simplification optimization, 0 or 1, default: 1\n");
        printf("<arg12> (optional) thread number, default: 15\n");
        printf("<arg13> (optional) percent for scalability test, default: 0\n");
        printf("<arg14> (optional) preprocessing task (1: Partitioned MDE Ordering; 2: Partitioned Query Generation; 3: Structural Update Generation)\n");
        exit(0);
    }

    string DesFile="./data/";
    string dataset = "NY";
    int algoChoice = 2;//1: CH, 2: H2H, 3: PLL
    int partitionNum = 20;
    int PSPStrategy = 2;//1: Pre, 2: No, 3: Post
    string algoParti = "NC";
    int updateType = 0;
    int runtimes = 10000;
    runtimes = 1000;
    int batchNum = 10;
    batchNum = 100;
    int batchSize=1;
    int threadNum = 15;
    int preTask=0;
    int cutType=EdgeCut;
    int ifFullOpt=true;
    double changeRatio=0.5;
    int percent=0;

    if(argc > 1) {
        cout << "argc: " << argc << endl;
        cout << "argv[1] (Source Path): " << argv[1] << endl;//source path
        DesFile = argv[1];

        cout << "argv[2] (Dataset): " << argv[2] << endl;//dataset
        dataset = argv[2];

        cout << "argv[3] (SP Index): " << argv[3] << endl;//system index
        algoChoice = stoi(argv[3]);

        if(argc > 4){
            cout << "argv[4] (PSP Strategy): " << argv[4] << endl;//PSP Strategy
            PSPStrategy = stoi(argv[4]);
        }

        if(argc > 5){
            cout << "argv[5] (Partition Method): " << argv[5] << endl;//partition method
            algoParti = argv[5];
            if(algoParti == "HEP" || algoParti == "CLUGP"){
                cutType=VertexCut;
            }
        }

        if(argc > 6){
            cout << "argv[6] (Partition Number): " << argv[6] << endl;//partition number
            partitionNum = stoi(argv[6]);
        }

        if(argc > 7){
            cout << "argv[7] (Update Type): " << argv[7] << endl;//update type
            updateType = stoi(argv[7]);
        }

        if(argc > 8){
            cout << "argv[8] (Batch Number): " << argv[8] << endl;//update number
            batchNum = stoi(argv[8]);
        }

        if(argc > 9){
            cout << "argv[9] (Batch Size): " << argv[9] << endl;//batch size
            batchSize = stoi(argv[9]);
        }

        if(argc > 10){
            cout << "argv[10] (Weight Change Ratio): " << argv[10] << endl;//edge weight change ratio
            changeRatio = stod(argv[10]);
        }

        if(argc > 11){
            cout << "argv[11] (Overlay Optimization): " << argv[11] << endl;//overlay optimization
            ifFullOpt = stoi(argv[11]);
        }

        if(argc > 12){
            cout << "argv[12] (Thread Number): " << argv[12] << endl;//thread number
            threadNum = stoi(argv[12]);
        }

        if(argc > 13){
            cout << "argv[13] (Percent for Scalability): " << argv[13] << endl;//percent
            percent = stoi(argv[13]);
        }

        if(argc > 14){
            cout << "argv[14] (Preprocessing Task): " << argv[14] << endl;//preprocessing task
            preTask = stoi(argv[14]);
        }

    }

	//used for running time calculation
    Timer tt0;
    tt0.start();

//    string graphfile="/media/TraminerData/mengxuan/MengxuanGraphWPSL/Cond/CondWeighted";
    string sourcePath=DesFile+"/"+dataset+"/";
    string ODfile=sourcePath+dataset+".query";
    string updateFile=sourcePath+dataset+".update";


    Graph g;
    g.threadnum=threadNum;//thread number of parallel computation (can be changed)
    g.sourcePath=sourcePath;
    g.ifParallel = true;
    g.dataset=dataset;
    g.algoChoice=algoChoice;
    g.PSPStrategy=PSPStrategy;
    g.algoParti=algoParti;
    g.partiNum=partitionNum;
    g.pNum=partitionNum;
    g.cutType=cutType;
    g.updateType=updateType;
    g.ifFullOpt=ifFullOpt;
    cout<<"Planar Dynamic PSP Test."<<endl;
    cout<<"Dataset: "<<dataset<<endl;
    if(g.algoChoice==CH){
        cout<<"SP Index: CH"<<endl;
    }else if(g.algoChoice==H2H){
        cout<<"SP Index: H2H"<<endl;
    }else if(g.algoChoice==PLL){
        cout<<"SP Index: PLL"<<endl;
    }else{
        cout<<"Wrong SP index! "<<g.algoChoice<<endl; exit(1);
    }
    if(g.PSPStrategy==PreBoundary){
        cout<<"PSP Strategy: Pre-boundary strategy"<<endl;
    }else if(g.PSPStrategy==NoBoundary){
        cout<<"PSP Strategy: No-boundary strategy"<<endl;
    }else if(g.PSPStrategy==PostBoundary){
        cout<<"PSP Strategy: Post-boundary strategy"<<endl;
    }else{
        cout<<"Wrong PSP strategy! "<<g.PSPStrategy<<endl; exit(1);
    }
    if(g.cutType==EdgeCut){
        cout<<"Cut Type: Edge-cut"<<endl;
    }else if(g.cutType==VertexCut){
        cout<<"Cut Type: Vertex-cut"<<endl;
    }
    if(g.ifFullOpt && g.PSPStrategy!=PreBoundary){
        cout<<"With overlay simplification optimization."<<endl;
    }else{
        cout<<"Without overlay simplification optimization."<<endl;
    }
//    string pMethods[8]={"NC","MT", "SC", "kahypar", "geometric", "Bubble", "HEP", "CLUGP"};
//    bool ifFind=false;
//    for(int i=0;i<pMethods->size();++i){
//        if(g.algoParti==pMethods[i]){
//            ifFind=true;
//            break;
//        }
//    }
//    if(ifFind){
//        cout<<"Partition method: "<<g.algoParti<<endl;
//    }else{
//        cout<<"Incorrect partition method! "<<g.algoParti<<endl; exit(1);
//    }

    cout<<"Partition number: "<<g.partiNum<<endl;
    cout<<"Thread number: "<<g.threadnum<<endl;

    if(percent!=0){
        ODfile=sourcePath+dataset+"_"+ to_string(percent)+".query";
        updateFile=sourcePath+dataset+"_"+ to_string(percent)+".update";
        g.percent=percent;
    }

    if(preTask==1){
        if(cutType==EdgeCut){
            g.EdgeCutVertexOrdering(2);//Boundary-first MDE ordering
//            g.EdgeCutVertexOrdering(3);//Degree-based Boundary-first ordering
        }else if(cutType==VertexCut){
            g.VertexCutVertexOrdering();
        }

    }
    else if(preTask==2){
        g.QueryGenerationParti(true);//same-partition query generation

    }
    else if(preTask==3){
        g.StructuralUpdateGeneration(1000);
    }

    ///Task 1: Index construction
    g.IndexConstruction();



    ///Task 2: Query processing
//    g.CorrectnessCheck(100);
    g.EffiCheck(ODfile,runtimes);//query efficiency test
//    g.EffiCheck(sourcePath+"partitions/"+dataset+"_"+algoParti+"_"+to_string(partitionNum)+"/same_parti.query",runtimes);//same-partition query efficiency test
//    g.EffiCheck(sourcePath+"partitions/"+dataset+"_"+algoParti+"_"+to_string(partitionNum)+"/both_core.query",runtimes);//same-partition query efficiency test
//    g.EffiCheck(sourcePath+"partitions/"+dataset+"_"+algoParti+"_"+to_string(partitionNum)+"/core_parti.query",runtimes);//same-partition query efficiency test
//    g.EffiCheck(sourcePath+"partitions/"+dataset+"_"+algoParti+"_"+to_string(partitionNum)+"/parti_parti.query",runtimes);//same-partition query efficiency test
    for(int i=1;i<11;++i){
//        g.EffiCheck(sourcePath+dataset+"_spatial.Q"+ to_string(i),runtimes);//same-partition query efficiency test
//        g.EffiCheck(sourcePath+dataset+"_SP.Q"+ to_string(i),runtimes);//same-partition query efficiency test
    }

//    exit(0);
    ///Task 3: Index update
//    g.WriteCoreIndex(sourcePath+dataset+"Core");
//    g.WriteCoreGraph(sourcePath+dataset+"Core");
    g.IndexMaintenance(updateType,batchNum,batchSize,changeRatio);//index maintenance

    tt0.stop();
    cout<<"\nOverall runtime: "<<tt0.GetRuntime()<<" s."<<endl;
    cout<<"------------------\n"<<endl;
//    exit(0);
    g.clear();
	return 0;
}
