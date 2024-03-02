#include <iostream>
#include "BatchHLWW.h"
#include "BatchHLWW.hpp"
#include <sys/stat.h>

using namespace std;

int main(int argc, char **argv) {
    int runtimes = 1000;
    int updateType = INCREASE;//DECREASE INCREASE MIX
    int updateBatch = 20;
	updateBatch=1;
    int updateVolume = 1;
    bool ifParallel = true;

    if( argc < 5 || argc > 7){
        printf("usage:\n<arg1> source path, e.g. /Users/zhouxj/Documents/1-Research/Datasets\n");
        printf("<arg2> number of landmarks, e.g. 20\n");
        printf("<arg3> dataset, e.g. NY\n");
        printf("<arg4> task, 1: construct index, 2: index update, 3: query processing. e.g. 1\n");
        printf("<arg5> (only for index update) method type, e.g. 0: BHL+, 1: BHL\n");
        printf("<arg6> (only for index update) update type, e.g. 0: Increase, 1: Decrease\n");
        exit(0);
    }

    Timer tt0;
    tt0.start();
  if (argc > 1) {
    int k = atoi(argv[2]);
    string sourcePath = argv[1];
    string graphFile = sourcePath + "/" + argv[3] + "/" + argv[3];
    cout << "Loading Graph..." << std::endl;
//    HighwayLabellingWW hl(graphFile, k);
      HighwayLabellingWW hl;
      hl.dataset = argv[3];
        hl.threadnum = 150;/// thread number

      cout<<"Dataset: "<<hl.dataset<<endl;
      cout<<"Landmark number: "<<k<<endl;
		cout<<"Thread number: "<<hl.threadnum<<endl;
	  hl.ReadGraph(graphFile, k);

    /// connectivity checking
    if(ifDebug){
//        hl.CheckCC();
    }

      string dirname = sourcePath + "/" + hl.dataset + "/tmp";
      mkdir(dirname.c_str(), 0777);
    hl.landmarkTopK.resize(k);

    if(string(argv[4]).compare("1") == 0) {///construct index
//        if(argc != 5){
//            cout<<"Wrong parameter number!"<<endl;
//            exit(1);
//        }
      hl.SelectLandmarks_HD(hl.landmarkTopK ,graphFile);
      hl.SaveLandmarks(dirname,hl.landmarkTopK);
      cout << "Constructing Highway Cover Labelling..." << std::endl;
      hl.BuildIndex(hl.landmarkTopK,ifParallel);
      hl.storeLabelling(dirname); //storing labelling to disk
    } else if(string(argv[4]).compare("2") == 0) {///index update
//        if(argc != 7){
//            cout<<"Wrong parameter number!"<<endl;
//            exit(1);
//        }
//        updateType = atoi(argv[6]);
      hl.LoadLandmarks(dirname,hl.landmarkTopK);
      hl.loadLabelling_Full(dirname, hl.landmarkTopK); //loading labelling from disk
      if(atoi(argv[5]) == 0){
          cout<<"Update method: BHL+"<<endl;
      }else if(atoi(argv[5]) == 1){
          cout<<"Update method: BHL"<<endl;
      }
	  updateType=stoi(argv[6]);
      if(updateType == DECREASE){
          cout<<"Update type: Decrease"<<endl;
      }else if(updateType == INCREASE){
          cout<<"Update type: Increase"<<endl;
      }else{
          cout<<"Wrong update type!"<<endl;
      }
      cout<<"Update batch: "<<updateBatch<<" ; Update volume: "<<updateVolume<<endl;
      cout << "Updating Highway Cover Labelling..." << std::endl;
      hl.UpdateLabellingW(graphFile, atoi(argv[5]),updateType,updateBatch,updateVolume,ifParallel);
//      if(updateType == INCREASE){
//          updateType = DECREASE;
//      }else if(updateType == DECREASE){
//          updateType = INCREASE;
//      }
//        if(updateType == DECREASE){
//            cout<<"Update type: Decrease"<<endl;
//        }else if(updateType == INCREASE){
//            cout<<"Update type: Increase"<<endl;
//        }else{
//            cout<<"Wrong update type!"<<endl;
//        }
//      hl.UpdateLabellingW(graphFile, atoi(argv[5]),updateType,updateBatch,updateVolume,ifParallel);
//      hl.UpdateLabelling(argv[4], atoi(argv[5]));

//      hl.storeLabelling(string(argv[4])+"new_"); //storing labelling to disk after update
//      hl.deallocate();
    } else if (string(argv[4]).compare("3") == 0) {///query processing
//        if(argc != 5){
//            cout<<"Wrong parameter number!"<<endl;
//            exit(1);
//        }
      hl.LoadLandmarks(dirname,hl.landmarkTopK);
      hl.loadLabelling_Pruned(dirname); //loading labelling from disk
      hl.RemoveLandmarks(hl.landmarkTopK);
      cout << "Querying Highway Cover Labelling..." << std::endl;
      if(ifDebug){
          /// correctness checking
          hl.CorrectnessCheck(runtimes, false);
      }
//      hl.QueryDistance(argv[2], argv[6], 100);
      hl.EfficencyTest(graphFile,runtimes);
    }
  }
  tt0.stop();
    cout<<"Overall run time: "<<tt0.GetRuntime()<<endl;
  cout<<"--------------------\n"<<endl;
	return 0;
}
