/*
 * Filename:    main.cpp
 * Description: execution entry for hierarchical PLL algorithm
 * Created:     25 August 2022
 * Authors:     Xinjie ZHOU
 */
#include <sys/stat.h>
#include "gtree.hpp"
#include "gstartree.hpp"

bool ifDebug = false;//

int main(int argc, char** argv){
    Timer tt;
    tt.start();
    Gstartree gt;
//    gt.dataset = "Rome"; // cal FLA
//    gt.dataset = "FLA"; // cal FLA
    gt.dataset = "NY"; // cal FLA
    gt.thread_num = 100;
    int runtimes = 100;
    task_type = 1; //1: G-Tree building; 2:G*-Tree (shortcuts) building; 3: G-Tree update; 4: G*-Tree query; 5: generate OD pairs.
    bool ifGstar = false;   //if use shortcuts of G*-tree
    int updateType = DECREASE; //update type DECREASE INCREASE
    int updateVolume = 1;//the update volume of each batch
    int updateBatch = 10;//the number of batch
    bool ifPrune = true;//whether to use pruning strategy for G-tree update
    bool ifParallel = true;//whether to use multi-thread to speed up

    if( argc != 5 && argc != 6){
        printf("usage:\n<arg1> source path, e.g. ~/datasets\n");
        printf("<arg2> name of dataset, e.g. NY\n");
        printf("<arg3> partition number, e.g. 256\n");
        printf("<arg4> task type, e.g. 1: G-Tree building; 2: Query processing; 3: G-Tree update; 4: generate query OD pairs; 5: generate update OD pairs. default: 1\n");
        printf("<arg5> (optional) update type, e.g. 0: Decrease; 1: Increase. default: 0\n");
        exit(0);
    }


    if(argc > 1) {
        cout << "argc: " << argc << endl;
        cout << "argv[1]: " << argv[1] << endl;//partition number
        DataPath = argv[1];
        if(argc > 2){
            cout << "argv[2]: " << argv[2] << endl;//dataset
            gt.dataset = argv[2];
        }

        if(argc > 3){
            cout << "argv[3]: " << argv[3] << endl;//data path

            gt.parti_num = atoi(argv[3]);
        }
        if(argc > 4){
            cout << "argv[4]: " << argv[4] << endl;//task type
            task_type = atoi(argv[4]);
        }
        if(argc > 5){
            cout << "argv[5]: " << argv[5] << endl;//task type
            updateType = atoi(argv[5]);
        }
    }
    if(ifDebug){
	    cout<<"Debug mode."<<endl;
    }else{
	    cout<<"Test mode."<<endl;
    }
    cout << "The dataset is " << gt.dataset << "." <<endl;
    cout<<"Dataset: "<<gt.dataset<<endl;
    cout<<"Leaf node size: "<<LEAF_CAP<<endl;
    cout<<"Thread number: "<<gt.thread_num<<endl;

    dirname = DataPath + "/" + argv[2] + "/tmp";
    mkdir(dirname.c_str(), 0777);

    gt.PathInit(false);
    if(ifParallel) {
	    cout<<"With multi-thread computation!"<<endl;
    }
    else cout<<"Without multi-thread computation!"<<endl;
    switch (task_type) {
        case 1: {
            /// G-tree building
            gt.gtree_build(ifParallel);
//            break;
        }
        case 2:{
            /// Query processing
            gt.dist_main(runtimes,ifGstar);
//            break;
        }
        case 3:{
            /// G-tree update
            cout<<"Update batch: "<<updateBatch<<endl;
            cout<<"Update volume of each batch: "<<updateVolume<<endl;
            switch (updateType) {
                case DECREASE:{
                    Gstartree gt1 = gt;
                    cout<<"Update type: Decrease!"<<endl;
                    gt1.GtreeUpdateParalel(updateType,updateVolume,updateBatch,ifParallel);
//                    break;
                }
                case INCREASE:{
                    cout<<"Update type: Increase!"<<endl;
                    gt.GtreeUpdateParalel(INCREASE,updateVolume,updateBatch,ifParallel);
                    break;
                }
                default:{
                    cout<<"Wrong update type!"<<endl;
                    break;
                }
            }
            break;
        }
        case 5:{
            /// generate query OD pairs
            gt.ODpairGenerate(10000);
            break;
        }
        case 6:{
            /// generate update OD pairs
            gt.UpdateGenerate(1000);
            break;
        }
        default:
            break;
    }

    tt.stop();
    cout << "\nOverall time cost:" << tt.GetRuntime() << " seconds." << endl;
    cout<<"--------------------\n"<<endl;
}
