/*
 * Filename:    main.cpp
 * Description: execution entry for G-tree algorithm
 * Created:     25 August 2022
 * Authors:     Xinjie ZHOU
 */
#include <sys/stat.h>
#include "gtree.hpp"
#include "gstartree.hpp"

int main(int argc, char** argv){
    if( argc < 6 || argc >8){
        printf("usage:\n<arg1> source path, e.g. ~/datasets\n");
        printf("<arg2> name of dataset, e.g. NY\n");
        printf("<arg3> leaf node size, e.g. 64\n");
        printf("<arg4> index type, e.g. 1: G-Tree. default: 1\n");
        printf("<arg5> task type, e.g. 1: Index construction and Querying; 2: Only Query processing; 3: Only Index update. default: 1\n");
        printf("<arg6> (optional) update type, e.g. 1: decrease update; 2: increase update. default: 1\n");
        printf("<arg7> (optional) update number. default: 1000\n");
        exit(0);
    }

    Timer tt;
    tt.start();
    Gstartree gt;
//    gt.dataset = "Rome"; // cal FLA
//    gt.dataset = "FLA"; // cal FLA
    gt.dataset = "NY"; // cal FLA
    int runtimes = 1000;
    int index_type = 1;
    int task_type = 1; //1: G-Tree building; 2:G*-Tree (shortcuts) building; 3: G-Tree update; 4: G*-Tree query; 5: generate OD pairs.
    bool ifGstar = false;   //if use shortcuts of G*-tree
//    ifGstar=true;
    int updateType = INCREASE; //update type DECREASE INCREASE
    int updateVolume = 1;//the update volume of each batch
    int updateBatch = 100;//the number of batch
    updateBatch=10;
    bool ifPrune = true;//whether to use pruning strategy for G-tree update
    gt.ifParallel = true;//whether to use multi-thread to speed up
//  gt.ifParallel = false;
//    ifMac = false;//false
    gt.ifDebug=true;
    gt.ifDebug= false;
    if(argc > 1) {
        cout << "argc: " << argc << endl;
        cout << "argv[1] (Source path): " << argv[1] << endl;//source path
        gt.DataPath = argv[1];
        if(argc > 2){
            cout << "argv[2] (Dataset): " << argv[2] << endl;//dataset
            gt.dataset = argv[2];
        }
        if(argc > 3){
            cout << "argv[3] (Leaf node size): " << argv[3] << endl;//leaf node size
            LEAF_CAP = atoi(argv[3]);
        }
        if(argc > 4){
            cout << "argv[4] (Index type): " << argv[4] << endl;//index type
            gt.indexType = atoi(argv[4]);
        }
        if(argc > 5){
            cout << "argv[5] (Task type): " << argv[5] << endl;//task type
            task_type = atoi(argv[5]);
        }
        if(argc > 6){
            cout << "argv[6] (Update type): " << argv[6] << endl;//update type
            updateType = atoi(argv[6]);
        }
        if(argc > 7){
            cout << "argv[7] (Update number): " << argv[7] << endl;//update number
            updateBatch = atoi(argv[7]);
        }
    }
    if(gt.ifDebug){
	    cout<<"Debug mode."<<endl;
    }else{
	    cout<<"Test mode."<<endl;
    }
    cout << "The dataset is " << gt.dataset << "." <<endl;

    gt.dirname = gt.DataPath + "/" + argv[2] + "/tmp";
    mkdir(gt.dirname.c_str(), 0777);

    gt.PathInit(gt.ifMac);

    switch (task_type) {
        case 1:
            // index construction
            gt.IndexConstruction();
            // query processing
            gt.dist_main(runtimes);
            break;
        case 2:
            /// querying
            gt.init_gstarQ();
            gt.dist_main(runtimes);
            break;
        case 3:
            /// G-tree update
            gt.IndexMaintenance(updateType,updateBatch,updateVolume);
            /*if(updateType == INCREASE){
		    updateType = DECREASE;
		    cout<<"Update type: Decrease!"<<endl;
	    }else if(updateType == DECREASE){
		    updateType = INCREASE;
		    cout<<"Update type: Increase!"<<endl;
	    }
            gt.GtreeUpdateParalel(updateType,updateVolume,updateBatch,ifParallel);*/
	        break;
        case 4: /// generate query OD pairs
            gt.ODpairGenerate(10000);
            break;
        case 5:/// generate update OD pairs
            gt.UpdateGenerate(1000);
            break;
        default:
            break;
    }

    tt.stop();
    cout << "\nOverall time cost:" << tt.GetRuntime() << " seconds." << endl;
}