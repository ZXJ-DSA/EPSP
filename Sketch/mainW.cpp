#include <iostream>
#include "BatchHLW.h"
#include "BatchHLW.hpp"
#include <sys/stat.h>

using namespace std;

int main(int argc, char **argv) {
    if( argc != 5 && argc != 6){
        printf("usage:\n<arg1> source path, e.g. /export/project/xzhouby\n");
        printf("<arg2> dataset, e.g. NY\n");
        printf("<arg3> number of landmarks, e.g. 20\n");
        printf("<arg4> task, 1: construct index, 2: query processing, 3: index update. e.g. 1\n");
        printf("<arg5> (optional) update type, 0: Decrease, 1: Increase. e.g. 0\n");
        exit(0);
    }
    cout<<"This is test for Sketch!"<<endl;

    int runtimes = 100;
    int taskType = 1;
    int updateType = 0;//DECREASE INCREASE MIX
    int updateBatch = 10;
    int updateVolume = 1;
    bool ifParallel = true;

    string dataset;
    string sourcePath;
    string graphFile;
    int landmarkNum=10;

    if(argc > 1) {
        cout << "argc: " << argc << endl;
        cout << "argv[1]: " << argv[1] << endl;//source path
        sourcePath = argv[1];
        if(argc > 2){
            cout << "argv[2]: " << argv[2] << endl;//dataset
            dataset = argv[2];

        }
        if(argc > 3){
            cout << "argv[3]: " << argv[3] << endl;
            landmarkNum = stoi(argv[3]);
        }
        if(argc > 4){
            cout << "argv[4]: " << argv[4] << endl;
            taskType = stoi(argv[4]);
        }
        if(argc > 5){
            cout << "argv[5]: " << argv[5] << endl;
            updateType = stoi(argv[5]);
        }
    }

    graphFile = sourcePath + "/" + dataset + "/" + dataset;


    Timer tt0;
    tt0.start();

//    HighwayLabellingWW hl(graphFile, k);
    HighwayLabellingWW hl;
    hl.dataset = dataset;
    hl.threadnum = 10;/// thread number

    cout<<"Dataset: "<<hl.dataset<<endl;
    cout<<"Landmark number: "<<landmarkNum<<endl;
    cout<<"Thread number: "<<hl.threadnum<<endl;
    cout << "Loading Graph..." << std::endl;
    hl.ReadGraph(graphFile, landmarkNum);

    /// connectivity checking
    if(ifDebug){
        hl.CheckCC();
    }

    string dirname = sourcePath + "/" + hl.dataset + "/tmp";
    mkdir(dirname.c_str(), 0777);
    hl.landmarkTopK.resize(landmarkNum);

    switch (taskType) {
        case 1:{
            ///Task 1: Index construction
            hl.SelectLandmarks_HD(hl.landmarkTopK ,graphFile);
            hl.SaveLandmarks(dirname,hl.landmarkTopK);
            cout << "Constructing Highway Cover Labelling..." << std::endl;
            hl.BuildIndex(hl.landmarkTopK,ifParallel);
            hl.storeLabelling(dirname); //storing labelling to disk
//            break;
        }
        case 2:{
            ///Task 2: Query processing
            HighwayLabellingWW hl1=hl;
            if(taskType == 2){
                hl1.LoadLandmarks(dirname,hl1.landmarkTopK);
            }
            hl1.loadLabelling_Pruned(dirname); //loading labelling from disk
            hl1.RemoveLandmarks(hl1.landmarkTopK);
            cout << "Querying Highway Cover Labelling..." << std::endl;
            if(ifDebug){
                /// correctness checking
                hl1.CorrectnessCheck(runtimes, false);
            }
            hl1.EfficencyTest(graphFile,runtimes);
//            break;
        }
        case 3:{
            ///Task 3: Index update
            if(taskType == 3){
                hl.LoadLandmarks(dirname,hl.landmarkTopK);
            }
            hl.loadLabelling_Full(dirname, hl.landmarkTopK); //loading labelling from disk
            cout<<"Update batch: "<<updateBatch<<" ; Update volume: "<<updateVolume<<endl;
            switch (updateType) {
                case DECREASE:{
                    HighwayLabellingWW hl2=hl;
                    cout<<"Update type: Decrease"<<endl;
                    hl2.UpdateLabellingW(graphFile, 1,DECREASE,updateBatch,updateVolume,ifParallel);
//                    break;
                }
                case INCREASE:{
                    cout<<"Update type: Increase"<<endl;
                    hl.UpdateLabellingW(graphFile, 1,INCREASE,updateBatch,updateVolume,ifParallel);
                    break;
                }
                default:{
                    cout<<"Wrong update type!"<<endl;
                    break;
                }
            }

        }
    }



    tt0.stop();
    cout<<"Overall run time: "<<tt0.GetRuntime()<<endl;
    cout<<"--------------------\n"<<endl;
	return 0;
}
