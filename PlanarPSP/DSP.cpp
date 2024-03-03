/*
 * H2H.cpp
 *
 *  Created on: 22 Dec 2022
 *      Author: Xinjie Zhou
 */
#include "headSP.h"
#include "H2H.hpp"
#include "PLL.hpp"


///////////////////////////// Index Construction ////////////////////////////////

void Graph::IndexConstruction(int algo){
    algoIndex = algo;

    switch (algoIndex) {
        case 0:{
            cout<<"Dijkstra's algorithm."<<endl;

            break;
        }
        case 1:{
            cout<<"CH algorithm."<<endl;
            CHIndexConstruct();
            break;
        }
        case 2:{
            cout<<"H2H algorithm."<<endl;
            H2HIndexConstruct();
            break;
        }
        case 3:{
            cout<<"PLL algorithm."<<endl;
//            PLLIndexConstruction(1);//PLL
//            PLLIndexConstruction(2);//PSL
            PLLIndexConstruction(3);//BPCL
            break;
        }
        default:{
            cout<<"Wrong SP index! "<<algoIndex<<endl; exit(1);
        }
    }


}



///////////////////////////// Query Processing ////////////////////////////////

//function for correctness check
void Graph::CorrectnessCheck(int runtimes){
    if(algoIndex==0){
        return;
    }
    Timer tt;
    double runT=0;
    srand (time(NULL));
    int s, t, d1=0, d2=0, d3;
//    runtimes = 1;
    cout<<"Correctness check ("<<runtimes<<" rounds) ... ";
    for(int i=0;i<runtimes;i++){
//        if(i%100==0) cout<<i<<endl;
        s=rand()%node_num;
        t=rand()%node_num;
//        s=52413,t=21929;//NY
//        cout<<"Query "<<i<<": "<<s<<" "<<t<<endl;

//        if(runtimes == 1){
//            cout<<"s: "<<s<<" ; t: "<<t<<endl;
//        }
        tt.start();
        d1=Dijkstra(s,t,Neighbor);
        tt.stop();
        if(algoIndex==0){
            d2=d1;
        }
        if(algoIndex==1){
            tt.start();
            d2=QueryCHWP(s,t);
            tt.stop();
        }else if(algoIndex==2){
            tt.start();
            d2=QueryH2H(s,t);
            tt.stop();
        }else if(algoIndex==3){
            tt.start();
            d2=QueryPLL(s,t);
            tt.stop();
        }

        runT+=tt.GetRuntime();
//        cout<<s<<"("<<CoreTag[s]<<") "<<t<<"("<<CoreTag[t]<<") "<<d2<<" "<<d1<<endl;
        if(d1!=d2){
            cout<<"InCorrect! "<<i<<": "<<s<<"("<<NodeOrder[s]<<") "<<t<<"("<<NodeOrder[t]<<") "<<d2<<" "<<d1<<endl;
            exit(1);
        }
    }
    cout<<"Average Query Time: "<<1000*runT/runtimes<<" ms."<<endl;
}

void Graph::EffiCheck(int runtimes){
    string ODfile=sourcePath+dataset+".query";
    ifstream IF(ODfile);
    if(!IF){
        cout<<"Cannot open Map "<<ODfile<<endl;
        exit(1);
    }
    cout<<"Query file: "<<ODfile<<endl;
    int num, ID1, ID2;
    vector<pair<int,int>> ODpair;
    IF>>num;
    for(int k=0;k<num;k++){
        IF>>ID1>>ID2;
        ODpair.push_back(make_pair(ID1, ID2));
    }
    if(runtimes > num){
        runtimes = num;
    }
    cout<<"Efficiency test. Run times: "<<runtimes<<endl;
    int s, t;
    double runT=0;
    int d1, d2;
    Timer tt;
    clock_t start = clock();

    vector<int> results(runtimes,-1);
    for(int i=0;i<runtimes;i++){
        s=ODpair[i].first; t=ODpair[i].second;
//        d1=Dijkstra(s,t,Neighbor);
        if(algoIndex==1){
            tt.start();
            d2=QueryCHWP(s,t);
            tt.stop();
        }
        else if(algoIndex==2){
            tt.start();
            d2=QueryH2H(s,t);
            tt.stop();
        }
        else if(algoIndex==3){
            tt.start();
            d2=QueryPLL(s,t);
            tt.stop();
        }

        runT += tt.GetRuntime();
        results[i]=d2;
//        cout<<d2<<" ";//<<endl;
//        if(d1!=d2){
//            cout<<"Wrong! "<<s<<" "<<t<<" "<<d2<<" "<<d1<<endl; exit(1);
//        }
    }

//    cout<<endl;

    cout<<"Average Query Time: "<<1000*runT/runtimes<<" ms. "<<1000*(double)(clock() - start) / (CLOCKS_PER_SEC*runtimes)<<" ms."<<endl;
}



///////////////////////////// Index Maintenance ////////////////////////////////

void Graph::IndexMaintenance(int updateType, int batchNum, bool ifBatch, int batchSize) {
    cout<<"Index update test..."<<endl;
    // read updates
    string file = sourcePath+dataset + ".update";
    bool ifDebug=false;
//    ifDebug=true;
    vector<pair<pair<int,int>,pair<int,int>>> wBatch;
    int ID1, ID2, oldW, newW;
    srand (0);

    cout<<"Update batch: "<<batchNum<<" ; Batch size: "<<batchSize<<endl;
    vector<pair<pair<int,int>,int>> updateData;
    ReadUpdate(file, updateData);
    Timer tt;
    double runT1=0, runT2 = 0;
    switch (updateType) {
        case 0:{
            break;
        }
        case 1:{
            //Decrease update
            cout<<"Update type: Decrease"<<endl;
            Graph g2=*this;
            if(ifBatch && algoIndex!=3){//for batch update
                if(batchNum*batchSize>updateData.size()){
                    batchNum=floor(updateData.size()/batchSize);
                }
                vector<pair<pair<int,int>,pair<int,int>>> wBatch;
                for(int u=0;u<batchNum;u++){
                    wBatch.clear();
                    for(int i=0;i<batchSize;++i){
                        ID1 = updateData[u*batchSize+i].first.first;
                        ID2 = updateData[u*batchSize+i].first.second;
                        oldW = updateData[u*batchSize+i].second;
                        newW=oldW*0.5;
                        if(newW < 1) {
                            cout<<"New edge weight is not positive! "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                            exit(1);
                        }
                        if(ifDebug){
                            cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                        }
                        wBatch.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
                    }

                    if(algoIndex==1){
                        tt.start();
                        g2.CHdecBat(wBatch);
//                        CHdecBat(wBatch);
                        tt.stop();
                    }
                    else if(algoIndex==2){
                        tt.start();
                        g2.H2HdecBat(wBatch);
//                        H2HdecBat(wBatch);
                        tt.stop();
                    }

                    runT1 += tt.GetRuntime();
                    if(ifDebug){
                        g2.CorrectnessCheck(100);
//                        CorrectnessCheck(100);
                    }
                }
                cout<<"Average Decrease batch update Time: "<<runT1/(batchNum*batchSize)<<" s."<<endl;
            }
            else{//for single-edge update
                for(int u=0;u<batchNum;u++){
                    wBatch.clear();
                    ID1 = updateData[u].first.first;
                    ID2 = updateData[u].first.second;
                    oldW = updateData[u].second;
                    newW=oldW*0.5;
                    if(newW < 1) {
                        cout<<"New edge weight is not positive! "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                        exit(1);
                    }
                    if(ifDebug){
                        cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                    }
                    wBatch.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
                    if(algoIndex==1){
                        tt.start();
                        g2.CHdecBat(wBatch);
//                        CHdecBat(wBatch);
                        tt.stop();
                    }
                    else if(algoIndex==2){
                        tt.start();
                        g2.H2HdecBat(wBatch);
//                        H2HdecBat(wBatch);
                        tt.stop();
                    }
                    else if(algoIndex==3){
                        tt.start();
                        g2.DecreasePSL(ID1,ID2,oldW,newW,g2.Neighbor,g2.Label);
//                        DecreasePSL(ID1,ID2,oldW,newW,Neighbor,Label);
                        tt.stop();
                    }
                    runT1 += tt.GetRuntime();
                    if(ifDebug){
                        CorrectnessCheck(100);
                    }
                }
                cout<<"Average Decrease single-edge update Time: "<<runT1/batchNum<<" s."<<endl;
            }

//            break;
        }
        case 2:{
            //Increase update
            cout<<"Update type: Increase"<<endl;
            if(ifBatch && algoIndex!=3){//for batch update
                if(batchNum*batchSize>updateData.size()){
                    batchNum=floor(updateData.size()/batchSize);
                }
                vector<pair<pair<int,int>,pair<int,int>>> wBatch;
                for(int u=0;u<batchNum;u++){
                    wBatch.clear();
                    for(int i=0;i<batchSize;++i){
                        ID1 = updateData[u*batchSize+i].first.first;
                        ID2 = updateData[u*batchSize+i].first.second;
                        oldW = updateData[u*batchSize+i].second;
                        newW=oldW*1.5;
                        if(ifDebug){
                            cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                        }
                        wBatch.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
                    }

                    if(algoIndex==1){
                        tt.start();
                        CHincBatMT(wBatch);
                        tt.stop();
                    }
                    else if(algoIndex==2){
                        tt.start();
                        H2HincBatMT(wBatch);
                        tt.stop();
                    }

                    runT2 += tt.GetRuntime();
                    if(ifDebug){
//                        CorrectnessCheckH2H(100);
                    }

                }
                cout<<"Average Increase batch update Time: "<<runT2/(batchNum*batchSize)<<" s."<<endl;
            }
            else {//for single-edge update
                for(int u=0;u<batchNum;u++){
                    wBatch.clear();
                    ID1 = updateData[u].first.first;
                    ID2 = updateData[u].first.second;
                    oldW = updateData[u].second;
                    newW=oldW*2;
                    if(ifDebug){
                        cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<" "<<oldW<<" "<<newW<<endl;
                    }
                    wBatch.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
                    if(algoIndex==1){
                        tt.start();
                        CHincBatMT(wBatch);
                        tt.stop();
                    }
                    else if(algoIndex==2){
                        tt.start();
                        H2HincBatMT(wBatch);
                        tt.stop();
                    }
                    else if(algoIndex==3){
                        tt.start();
                        IncreasePSL(ID1,ID2,oldW,newW,Neighbor,Label,PruningPointSet,PruningPointSet2);
                        tt.stop();
                    }
                    runT2 += tt.GetRuntime();
                    if(ifDebug){
                        CorrectnessCheck(100);
                    }
                }
                cout<<"Average Increase single-edge update Time: "<<runT2/batchNum<<" s."<<endl;
            }
            break;
        }
        default:{
            cout<<"Wrong update type!"<<endl;
            break;
        }

    }
}

