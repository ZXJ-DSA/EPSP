//
// Created by Xinjie ZHOU on 25/8/2022.
//

#ifndef GSTARTREE_HPP
#define GSTARTREE_HPP

#include "gstartree.h"

void Gstartree::IndexConstruction(){
    switch (indexType) {
        case gtreeIndex:{//G-tree
            cout<<"Index type: G-Tree."<<endl;
            gtree_build(ifParallel);
            break;
        }
        case gstarIndex:{//G*-tree
            cout<<"Index type: G*-Tree."<<endl;
            gstartree_build();
            break;
        }
        case lgtreeIndex:{//LG-tree
            cout<<"Index type: LG-Tree."<<endl;
            LGTreeIndexBuild();
            break;
        }
        case tgtreeIndex:{
            cout<<"Index type: TG-Tree."<<endl;
            TDGTreeIndexBuild();
            break;
        }
        default:{
            cout<<"Wrong index type! "<<indexType<<endl; exit(1);
        }
    }

}

void Gstartree::TDGTreeIndexBuild(){
    // init
    init();
    double t1=0,t2=0,t3=0,t4=0,t5=0;

    /// Hierarchical tree building
    Timer tt;
    cout << "Start to build G-tree..."<<endl;
    tt.start();
    build();
//    cout <<"Done."<<endl;
    tt.stop();
    t1 = tt.GetRuntime();
    cout << "The time for G-tree building: " << t1 << " s." << endl;

    // dump gtree
//    cout << "Saving G-tree..."<<endl;
    if(percentScale==0){
        gtree_save(dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".TD.gtree.bin", dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".TD.paths.bin");
    }
    else{
        gtree_save(dirname+"/"+dataset+"_"+ to_string(percentScale) + "." + to_string(LEAF_CAP) + ".TD.gtree.bin", dirname+"/"+dataset+"_"+ to_string(percentScale) + "." + to_string(LEAF_CAP) + ".TD.paths.bin");
    }
//    gtree_save_txt();
//    cout <<"Done."<<endl;
//    exit(0);

    /// original overlay graph building
    tt.start();
    getOriginalOverlayGraph();
    tt.stop();
    t2 = tt.GetRuntime();
    cout << "The time for overlay graph computing: " << t2 << " s." << endl;
    // calculate distance matrix
    cout << "Start to constructing index..."<<endl;
    tt.start();
    hierarchy_shortest_path_calculation_NoBoundary(ifParallel);
//    hierarchy_shortest_path_calculation_NoBoundary(false);
//    cout << "Done."<<endl;
    tt.stop();
    t3 = tt.GetRuntime();
    cout << "The time for partition index building: " << t3 << " s." << endl;
    if(percentScale==0){
        hierarchy_shortest_path_save(dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".TD.minds.bin");
    }
    else{
        hierarchy_shortest_path_save(dirname+"/"+dataset+"_"+ to_string(percentScale) + "." + to_string(LEAF_CAP) + ".TD.minds.bin");
    }

    /// Construct L~
    if(ifHierarchicalOrdering){
        cout<<"Vertex ordering method: Hierarchical MDE Ordering."<<endl;
        tt.start();
        H2HconOrderMT(false);///
        tt.stop();
    }
    else{
        cout<<"Vertex ordering method: Naive MDE Ordering."<<endl;
        tt.start();
        H2HconOrderMT(true);///
        tt.stop();
    }
    t4=tt.GetRuntime();
    if(percentScale==0){
        graphInfoSave(dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".TD.graph.bin");
    }
    else{
        graphInfoSave(dirname+"/"+dataset+"_"+ to_string(percentScale) + "." + to_string(LEAF_CAP) + ".TD.graph.bin");
    }

    cout << "The time for overlay index building: " << t4 << " s." << endl;
    cout << "Overall time for index construction: " << t1+t2+t3+t4+t5 << " s." << endl;
    // dump distance matrix
//    cout << "Saving distance matrix..."<<endl;
//    H2HIndexSave();
//    cout << "Done."<<endl;
    IndexSizeH2H();
}

void Gstartree::LGTreeIndexBuild(){
    // init
    init();
    double t1=0,t2=0,t3=0,t4=0,t5=0;

    /// Hierarchical tree building
    Timer tt;
    cout << "Start to build G-tree..."<<endl;
    tt.start();
    build();
//    cout <<"Done."<<endl;
    tt.stop();
    t1 = tt.GetRuntime();
    cout << "The time for G-tree building: " << t1 << " s." << endl;

    // dump gtree
//    cout << "Saving G-tree..."<<endl;
    gtree_save(dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".LG.gtree.bin", dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".LG.paths.bin");
//    gtree_save_txt();
//    cout <<"Done."<<endl;
//    exit(0);

    /// original overlay graph building
    tt.start();
    getOriginalOverlayGraph();
    tt.stop();
    t2 = tt.GetRuntime();
    cout << "The time for overlay graph computing: " << t2 << " s." << endl;
    // calculate distance matrix
    cout << "Start to constructing index..."<<endl;
    tt.start();
    hierarchy_shortest_path_calculation_PreBoundary(ifParallel);
//    hierarchy_shortest_path_calculation_PreBoundary(false);
//    cout << "Done."<<endl;
    tt.stop();
    t3 = tt.GetRuntime();
    cout << "The time for partition index building: " << t3 << " s." << endl;
    hierarchy_shortest_path_save(dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".LG.minds.bin");
    /// Construct L~
    tt.start();
    VertexHierarchyBuild();//non-neighbor set based vertex hierarchy
    BorderLabelBuild();//border vertex label list computation
    tt.stop();

    t4=tt.GetRuntime();
    graphInfoSave(dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".LG.graph.bin");
    overlayIndexSave(dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".LG.PLL.bin");
    cout << "The time for overlay index building: " << t4 << " s." << endl;
    cout << "Overall time for index construction: " << t1+t2+t3+t4+t5 << " s." << endl;
    // dump distance matrix
//    cout << "Saving distance matrix..."<<endl;
//    H2HIndexSave();
//    cout << "Done."<<endl;
    IndexSizePLL();
}
///function of building the vertex hierarchy by non-neighbor set
void Gstartree::VertexHierarchyBuild(){
    Timer tt;
    tt.start();
    vertexLevels.clear();
    NeighborCon.assign(node_num,vector<pair<int,pair<int,int>>>());
    SCconNodesMT.assign(node_num, map<int, vector<int>>());
    levelTag.assign(node_num,-1);
    //get original overlay graph
    E.assign(node_num,unordered_map<int,pair<int,int>>());
    for(int i=0;i<OverlayGraph.size();i++){
        if(!OverlayGraph[i].empty()){
            for(auto it=OverlayGraph[i].begin();it!=OverlayGraph[i].end();++it){
                E[i].insert(make_pair(it->first,make_pair(it->second,1)));
            }
        }
    }
    _DD_.assign(node_num,0);

    set<DegComp1> Deg;//min first
    int degree;
    vector<bool> exist(node_num,false);//if in the core, all vertices is originally in core

    for(int ID=0;ID<OverlayGraph.size();++ID) {
        if (Nodes[ID].isborder) {
            degree = E[ID].size();
            exist[ID] = true;
            if (degree > 0) {//get degree
                _DD_[ID] = degree;
                Deg.insert(DegComp1(ID));
            } else {
                cout << "Wrong!! Degree of " << ID << " is " << degree << endl;
                exit(1);
            }
        }
    }

    vNodeOrder.clear();
    for(int i=0;i<node_num;++i){
        if(!Nodes[i].isborder){//if not border
            vNodeOrder.push_back(i);
        }
    }
    cout<<"overlay vertex number: "<<node_num-vNodeOrder.size()<<endl;
    vector<bool> change(node_num,false);//whether the neighbor (degree) has changed
//    unordered_set<int> neighborSet; neighborSet.clear();
    unordered_set<int> Temp; Temp.clear();
    vector<int> layerVertex; unordered_set<int> layerVertexSet; layerVertex.clear();

    bool CutLabel=false;
    int layer=0;
    int ID1,ID2;
    int count=0;

    int x;//minimum degree first

    while(!Deg.empty() || !Temp.empty()){
        if(!Temp.empty()){
            assert(Deg.empty());
            Deg.clear();
            for(auto it=Temp.begin();it!=Temp.end();++it) {
                int ID=*it;
                if (Nodes[ID].isborder) {
                    degree = E[ID].size();
//                    exist[ID] = true;
                    Deg.insert(DegComp1(ID));
                }
            }
            vertexLevels.push_back(layerVertex);
            layerVertex.clear(); layerVertexSet.clear();
            Temp.clear();
            layer++;
        }
        if(layer%500==0){
            cout<<"layer "<<layer<<": "<<count<<endl;
        }

        while(!Deg.empty()){
            x=(*Deg.begin()).x;//minimum degree first
            while(change[x]){//update the degree if it is changed
                Deg.erase(DegComp1(x));
                _DD_[x]=E[x].size();
                Deg.insert(DegComp1(x));
                change[x]=false;
                x=(*Deg.begin()).x;
            }
//        cout<<x<<" "<<E[x].size()<<endl;
            if(Temp.find(x)==Temp.end()){//if not found, add to layer vertex
                levelTag[x]=layer;
                layerVertex.push_back(x);
                layerVertexSet.insert(x);
                count++;
                vNodeOrder.push_back(x);
                Deg.erase(Deg.begin());

                vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();
                for(auto it=E[x].begin();it!=E[x].end();it++){
                    ID2=it->first;
                    if(exist[ID2]){
                        Neigh.push_back(*it);
                        if(Temp.find(ID2)==Temp.end()){
                            Temp.insert(ID2);
                        }
                    }else{
                        cout<<ID2<<" has been removed. "<<endl;
                    }
                }
                if(Neigh.empty() && !Deg.empty()){
                    cout<<"!!! Neigh is empty for "<<x<<" : "<<E[x].size()<<" "<<Nodes[x].adjnodes.size()<<" "<<count<<endl;
//            exit(1);
                }
                NeighborCon[x].assign(Neigh.begin(),Neigh.end());

                exist[x]=false;
                //delete the star
                for(int i=0;i<Neigh.size();i++){
                    int y=Neigh[i].first;
                    deleteE(x,y);//delete x from y's adjacency list
                    change[y]=true;
                }
                //add all-pair neighbors
                if(Neigh.size()<=100){
                    for(int i=0;i<Neigh.size();i++) {
                        ID1 = Neigh[i].first;
                        for (int j = i + 1; j < Neigh.size(); j++) {
                            ID2 = Neigh[j].first;
                            insertE(ID1, ID2, Neigh[i].second.first + Neigh[j].second.first);
                            /// For TD update

                            if(Neigh[i].first<Neigh[j].first)
                                SCconNodesMT[Neigh[i].first][Neigh[j].first].push_back(x);//no direction
                            else if(Neigh[j].first<Neigh[i].first)
                                SCconNodesMT[Neigh[j].first][Neigh[i].first].push_back(x);

                        }
                    }
                }
                else{
                    if(Neigh.size()>thread_num){
                        int step=Neigh.size()/thread_num;
                        boost::thread_group thread;
                        for(int i=0;i<thread_num;i++){
                            pair<int,int> p;
                            p.first=i*step;
                            if(i==thread_num-1)
                                p.second=Neigh.size();
                            else
                                p.second=(i+1)*step;
                            thread.add_thread(new boost::thread(&Gstartree::NeighborComorder, this, boost::ref(Neigh), p, x));
                        }
                        thread.join_all();
                    }else{
                        boost::thread_group thread;
                        for(int i=0;i<Neigh.size();i++){
                            pair<int,int> p;
                            p.first=i; p.second=(i+1);
                            thread.add_thread(new boost::thread(&Gstartree::NeighborComorder, this, boost::ref(Neigh), p, x));
                        }
                        thread.join_all();
                    }
                }
            }
            else{
                Deg.erase(Deg.begin());
            }
        }
    }
    if(vNodeOrder.size()!=node_num){
        cout<<"vNodeOrder size incorrect! "<<vNodeOrder.size()<<endl;
        exit(1);
    }
    NodeOrder.assign(node_num,-1);
    for(int k=0;k<vNodeOrder.size();k++){
        ID1=vNodeOrder[k];
        if(ID1>=0 && ID1<node_num){
            NodeOrder[ID1]=k;
        }else{
            cout<<"Wrong ID! "<<ID1<<" "<<k<<endl; exit(1);
        }
    }
    tt.stop();
    cout<<"Time for non-neighbor set based vertex ordering: "<<tt.GetRuntime()<<" s."<<endl;
    cout<<"Vertex hierarchy number: "<<vertexLevels.size()<<endl;
}

void Gstartree::BorderLabelBuild(){
    Timer tt;
    tt.start();
    cout<<"Begin border label construction..."<<endl;
    Label.assign(node_num,unordered_map<int,int>());
    int ID1,ID2,ID3,weight;
    for(int i=0;i<vertexLevels.size()-1;++i){
        for(int j=0;j<vertexLevels[i].size();++j){
            ID1=vertexLevels[i][j];
            Label[ID1].insert({ID1,0});
            for(auto it=NeighborCon[ID1].begin();it!=NeighborCon[ID1].end();++it){
                ID2=it->first; weight=it->second.first;
                Label[ID1].insert({ID2,weight});
            }
        }
    }
    int level=vertexLevels.size()-1;
//    cout<<"size of the highest level: "<<vertexLevels[level].size()<<endl;
    for(int j=0;j<vertexLevels[level].size();++j){
        ID1=vertexLevels[level][j];
        Label[ID1].insert({ID1,0});
    }
    for(int i=vertexLevels.size()-2;i>=0;--i){
//        if(i%500==0){
//            cout<<"level "<<i<<endl;
//        }
        for(int j=0;j<vertexLevels[i].size();++j){
            ID1=vertexLevels[i][j];
            for(auto it=NeighborCon[ID1].begin();it!=NeighborCon[ID1].end();++it){
                ID2=it->first; weight=it->second.first;
                for(auto it2=Label[ID2].begin();it2!=Label[ID2].end();++it2){
                    ID3=it2->first;
                    if(Label[ID1].find(ID3)==Label[ID1].end()){//if not found
                        Label[ID1].insert({ID3,Label[ID2][ID3]+weight});
                    }else{//if found
                        if(Label[ID1][ID3] > Label[ID2][ID3]+weight){
                            Label[ID1][ID3] = Label[ID2][ID3]+weight;
                        }
                    }
                }
            }
        }
    }
    tt.stop();
    cout<<"Time for label construction: "<<tt.GetRuntime()<<" s."<<endl;
}

//void Gstartree::leafNodeMatrix_save(){
//    char filename[300];
//    if(dataset == "cal"){
//        strcpy(filename, DataPath.c_str());
//        strcat(filename, FILE_ONTREE_MIND.c_str());
//    }else {
//        strcpy(filename, dirname.c_str());
//        strcat(filename, "/");
//        string temp=dataset + "." + to_string(LEAF_CAP) + ".TD.minds.bin";
//        strcat(filename, temp.c_str());
//    }
//    FILE* fout = fopen( filename, "wb" );
//    if(fout == NULL){
//        cout<<"Failed to open file "<<filename<<endl;
//        exit(1);
//    }
//    int* buf;
//    int count;
//    for ( int i = 0; i < GTree.size(); i++ ){
//        // union borders
//        count = GTree[i].union_borders.size();
//        fwrite( &count, sizeof(int), 1, fout );
//        buf = new int[count];
//        copy( GTree[i].union_borders.begin(), GTree[i].union_borders.end(), buf );
//        fwrite( buf, sizeof(int), count, fout );
//        delete[] buf;
//        // mind
//        count = GTree[i].mind.size();
//        fwrite( &count, sizeof(int), 1, fout );
//        buf = new int[count];
//        copy( GTree[i].mind.begin(), GTree[i].mind.end(), buf );
//        fwrite( buf, sizeof(int), count, fout );
//        delete[] buf;
//    }
//    fclose(fout);
//}

void Gstartree::H2HIndexSave(){
    char filename[300];

    strcpy(filename, dirname.c_str());
    strcat(filename, "/");
    string temp=dataset + "." + to_string(LEAF_CAP) + ".TD.H2H.bin";
    strcat(filename, temp.c_str());

    FILE* fout = fopen( filename, "wb" );
    if(fout == NULL){
        cout<<"Failed to open file "<<filename<<endl;
        exit(1);
    }
    int* buf;
    int count;
    /// Tree
    fwrite( &heightMax, sizeof(int), 1, fout );
    count = Tree.size();
    fwrite( &count, sizeof(int), 1, fout );
    for ( int i = 0; i < Tree.size(); i++ ){
        // uniqueVertex, pa, height, hdepth
        fwrite( &Tree[i].uniqueVertex, sizeof(int), 1, fout );
        fwrite( &Tree[i].pa, sizeof(int), 1, fout );
        fwrite( &Tree[i].height, sizeof(int), 1, fout );
        fwrite( &Tree[i].hdepth, sizeof(int), 1, fout );
        // vert
        count = Tree[i].vert.size();
        fwrite( &count, sizeof(int), 1, fout );
        for(int j=0;j<count;++j){
            fwrite( &Tree[i].vert[j].first, sizeof(int), 1, fout );
            fwrite( &Tree[i].vert[j].second.first, sizeof(int), 1, fout );
            fwrite( &Tree[i].vert[j].second.second, sizeof(int), 1, fout );
        }
        // dis
        count = Tree[i].dis.size();
        fwrite( &count, sizeof(int), 1, fout );
        buf = new int[count];
        copy( Tree[i].dis.begin(), Tree[i].dis.end(), buf );
        fwrite( buf, sizeof(int), count, fout );
        delete[] buf;
        // cnt
        count = Tree[i].cnt.size();
        fwrite( &count, sizeof(int), 1, fout );
        buf = new int[count];
        copy( Tree[i].cnt.begin(), Tree[i].cnt.end(), buf );
        fwrite( buf, sizeof(int), count, fout );
        delete[] buf;
        // pos
        count = Tree[i].pos.size();
        fwrite( &count, sizeof(int), 1, fout );
        buf = new int[count];
        copy( Tree[i].pos.begin(), Tree[i].pos.end(), buf );
        fwrite( buf, sizeof(int), count, fout );
        delete[] buf;
        // vAncestor
        count = Tree[i].vAncestor.size();
        fwrite( &count, sizeof(int), 1, fout );
        buf = new int[count];
        copy( Tree[i].vAncestor.begin(), Tree[i].vAncestor.end(), buf );
        fwrite( buf, sizeof(int), count, fout );
        delete[] buf;
        // ch
        count = Tree[i].ch.size();
        fwrite( &count, sizeof(int), 1, fout );
        buf = new int[count];
        copy( Tree[i].ch.begin(), Tree[i].ch.end(), buf );
        fwrite( buf, sizeof(int), count, fout );
        delete[] buf;
        // FN
        count = Tree[i].FN.size();
        fwrite( &count, sizeof(int), 1, fout );
        bool *buf2 = new bool[count];
        copy( Tree[i].FN.begin(), Tree[i].FN.end(), buf2 );
        fwrite( buf2, sizeof(int), count, fout );
        delete[] buf2;
    }
    /// rank
    count = rank.size();
    fwrite( &count, sizeof(int), 1, fout );
    buf = new int[count];
    copy( rank.begin(), rank.end(), buf );
    fwrite( buf, sizeof(int), count, fout );
    delete[] buf;
    /// SCconNodesMT
    count = SCconNodesMT.size();
    fwrite( &count, sizeof(int), 1, fout );
    for ( int i = 0; i < SCconNodesMT.size(); i++ ){
        // vert
        count = SCconNodesMT[i].size();
        fwrite( &count, sizeof(int), 1, fout );
        for(auto it=SCconNodesMT[i].begin();it!=SCconNodesMT[i].end();++it){
            fwrite( &it->first, sizeof(int), 1, fout );
            count = it->second.size();
            fwrite( &count, sizeof(int), 1, fout );
            buf = new int[count];
            copy( it->second.begin(), it->second.end(), buf );
            fwrite( buf, sizeof(int), count, fout );
            delete[] buf;
        }
    }
    /// VidtoTNid
    count = VidtoTNid.size();
    fwrite( &count, sizeof(int), 1, fout );
    for ( int i = 0; i < VidtoTNid.size(); i++ ){
        // vert
        count = VidtoTNid[i].size();
        fwrite( &count, sizeof(int), 1, fout );
        buf = new int[count];
        copy( VidtoTNid[i].begin(), VidtoTNid[i].end(), buf );
        fwrite( buf, sizeof(int), count, fout );
        delete[] buf;
    }

    fclose(fout);
}

void Gstartree::graphInfoSave(string filename){
//    char filename[300];
//    strcpy(filename, dirname.c_str());
//    strcat(filename, "/");
//    string temp=dataset + "." + to_string(LEAF_CAP) + ".TD.graph.bin";
//    strcat(filename, temp.c_str());

    FILE* fout = fopen( filename.c_str(), "wb" );
    if(fout == NULL){
        cout<<"Failed to open file "<<filename<<endl;
        exit(1);
    }
    int* buf;
    int count;
    /// Overlay graph
    for(int i=0;i<OverlayGraph.size();++i){
        count = OverlayGraph[i].size();
        fwrite( &count, sizeof(int), 1, fout );
        for(auto it=OverlayGraph[i].begin();it!=OverlayGraph[i].end();++it){
            fwrite( &it->first, sizeof(int), 1, fout );
            fwrite( &it->second, sizeof(int), 1, fout );
        }
    }
    /// Ordering
    count=node_num;
    buf = new int[count];
    copy( NodeOrder.begin(), NodeOrder.end(), buf );
    fwrite( buf, sizeof(int), count, fout );
    delete[] buf;
    fclose(fout);
}

void Gstartree::overlayIndexSave(string filename){
//    char filename[300];
//    strcpy(filename, dirname.c_str());
//    strcat(filename, "/");
//    string temp=dataset + "." + to_string(LEAF_CAP) + ".TD.graph.bin";
//    strcat(filename, temp.c_str());

    FILE* fout = fopen( filename.c_str(), "wb" );
    if(fout == NULL){
        cout<<"Failed to open file "<<filename<<endl;
        exit(1);
    }
    int* buf;
    int count;
    /// PLL
    for(int i=0;i<Label.size();++i){
        count = Label[i].size();
        fwrite( &count, sizeof(int), 1, fout );
        for(auto it=Label[i].begin();it!=Label[i].end();++it){
            fwrite( &it->first, sizeof(int), 1, fout );
            fwrite( &it->second, sizeof(int), 1, fout );
        }
    }
    /// NeighborCon
    for(int i=0;i<NeighborCon.size();++i){
        count = NeighborCon[i].size();
        fwrite( &count, sizeof(int), 1, fout );
        for(auto it=NeighborCon[i].begin();it!=NeighborCon[i].end();++it){
            fwrite( &it->first, sizeof(int), 1, fout );
            fwrite( &it->second.first, sizeof(int), 1, fout );
            fwrite( &it->second.second, sizeof(int), 1, fout );
        }
    }
    /// vertexLevels
    count = vertexLevels.size();
    fwrite( &count, sizeof(int), 1, fout );
    for(int i=0;i<vertexLevels.size();++i){
        count = vertexLevels[i].size();
        fwrite( &count, sizeof(int), 1, fout );
        buf = new int[count];
        copy( vertexLevels[i].begin(), vertexLevels[i].end(), buf );
        fwrite( buf, sizeof(int), count, fout );
        delete[] buf;
    }
    /// SCconNodesMT
    for(int i=0;i<SCconNodesMT.size();++i){
        count = SCconNodesMT[i].size();
        fwrite( &count, sizeof(int), 1, fout );
        for(auto it=SCconNodesMT[i].begin();it!=SCconNodesMT[i].end();++it){
            fwrite( &it->first, sizeof(int), 1, fout );
            count = it->second.size();
            fwrite( &count, sizeof(int), 1, fout );
            buf = new int[count];
            copy( it->second.begin(), it->second.end(), buf );
            fwrite( buf, sizeof(int), count, fout );
            delete[] buf;
        }
    }
    fclose(fout);
}

// load distance matrix from file
void Gstartree::load_overlayIndex(string filename) {
    FILE* fin = fopen( filename.c_str(), "rb" );
    if(fin == NULL){
        cout << "Failed to open file " << filename << endl;
        exit(1);
    }
    cout << "Loading overlay index... ";
    int* buf;
    int count, count2;
    int ID1, ID2, weight, num;
    /// PLL
    Label.assign(node_num,unordered_map<int,int>());
    for(ID1=0;ID1<Label.size();++ID1){
        fread( &count, sizeof(int), 1, fin );
        if(count>=1){
            for(int j=0;j<count;++j){
                fread( &ID2, sizeof(int), 1, fin );
                fread( &weight, sizeof(int), 1, fin );
                Label[ID1].insert({ID2,weight});
            }
        }
    }
    /// NeighborCon
    NeighborCon.assign(node_num,vector<pair<int,pair<int,int>>>());
    for(ID1=0;ID1<NeighborCon.size();++ID1) {
        fread(&count, sizeof(int), 1, fin);
        if (count >= 1) {
            for (int j = 0; j < count; ++j) {
                fread(&ID2, sizeof(int), 1, fin);
                fread(&weight, sizeof(int), 1, fin);
                fread(&num, sizeof(int), 1, fin);
                NeighborCon[ID1].emplace_back(ID2, make_pair(weight, num));
            }
        }
    }
    /// vertexLevels
    fread( &count, sizeof(int), 1, fin );
    vertexLevels.assign(count,vector<int>());
    for(int i=0;i<vertexLevels.size();++i){
        fread( &count, sizeof(int), 1, fin );
        buf = new int[count];
        fread( buf, sizeof(int), count, fin );
        for ( int j = 0; j < count; ++j ){
            vertexLevels[i].push_back(buf[j]);
        }
        delete[] buf;
    }
    /// SCconNodesMT
    SCconNodesMT.assign(node_num, map<int, vector<int>>());
    for(int i=0;i<SCconNodesMT.size();++i){
        fread( &count, sizeof(int), 1, fin );
        for(int j=0;j<count;++j){
            fread( &ID2, sizeof(int), 1, fin );
            SCconNodesMT[i].insert({ID2,vector<int>()});
            fread( &count2, sizeof(int), 1, fin );
            buf = new int[count2];
            fread( buf, sizeof(int), count2, fin );
            for ( int j = 0; j < count2; ++j ){
                SCconNodesMT[i][ID2].push_back(buf[j]);
            }
            delete[] buf;
        }
    }
    fclose(fin);
    cout<<"vertex levels: "<<vertexLevels.size()<<" Done."<<endl;
}

void Gstartree::IndexMaintenance( int updateType, int updateBatch,int updateVolume, bool ifRead){
    // init
    if(ifRead){
        if(indexType==gtreeIndex || indexType==gstarIndex){
            init_ReadIndex();
        }
        else if(indexType==lgtreeIndex){
            init_LGTreeIndex();
        }
        else if(indexType==tgtreeIndex){
            init_TGTreeIndex();
        }
        CorrectnessCheck(100);
    }

    // read updates
    string file = string(DataPath) + "/" + dataset + "/" + FILE_UPDATE;
    if(percentScale!=0){
        file = string(DataPath) + "/" + dataset + "/" + dataset+"_"+ to_string(percentScale)+".update";
    }
    vector<pair<pair<int,int>,int>> updateEdges;
    ReadUpdates(file,updateEdges);

//    if(updateType == INCREASE) cout<<"Update type: Increase!"<<endl;
//    if(updateType == DECREASE) cout<<"Update type: Decrease!"<<endl;
//    if(updateType == MIX) cout<<"Update type: Random!"<<endl;

    cout<<"Update batch: "<<updateBatch<<endl;
    cout<<"Update volume of each batch: "<<updateVolume<<endl;

    if(indexType==tgtreeIndex||indexType==lgtreeIndex){
        NodeOrder_=NodeOrder;
    }

//    set<OrderCompMax> OC;
//    OC.insert(OrderCompMax(30));
//    OC.insert(OrderCompMax(300));
//    OC.insert(OrderCompMax(30000));
//    while(!OC.empty()){
//        int x=OC.begin()->x;
//        cout<<x<<" "<<NodeOrder_[x]<<endl;
//        OC.erase(OC.begin());
//    }

    bool ifDebug= false;
//    ifDebug=true;

    double ave_time = 0;

    int ID1,ID2,oldW,newW;
    vector<pair<pair<int,int>,pair<int,int>>> updates;
    Timer tt;
    ///preprocess: get the update edges
    assert(updateVolume<=updateEdges.size());
    switch (updateType) {
        case DECREASE:{
            cout<<"\nUpdate type: Decrease"<<endl;
            ave_time = 0;
            auto GTreeTemp=GTree; auto NodesTemp=Nodes; auto OverlayGraphTemp=OverlayGraph; auto TreeTemp=Tree; auto LabelTemp=Label; auto NeighborConTemp=NeighborCon;
            for(int batch_i=0;batch_i<updateBatch;++batch_i){
                cout<<"Batch "<<batch_i<<" :";
                updates.clear();
                for(int i=batch_i;i<batch_i+updateVolume;++i) {//for each edge update
                    ID1 = updateEdges[i].first.first;
                    ID2 = updateEdges[i].first.second;
                    oldW = updateEdges[i].second;
                    newW = (int) (0.5 * oldW);
                    updates.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
                    if(ifDebug){
                        cout<<" "<<ID1<<"("<<Nodes[ID1].isborder<<","<<Nodes[ID1].inleaf<<") "<<ID2<<"("<<Nodes[ID2].isborder<<","<<Nodes[ID2].inleaf<<") "<<oldW<<" "<<newW;
                    }
                }
                if(ifDebug){
                    cout<<endl;
                }

//                cout<<"Before update: "<<endl;
//                int s=202763, t=73343;
//                s=196894, t=145058;
//                int dis1=Query_LGTree(s,t);
//                cout<<s<<"("<<Nodes[s].isborder<<","<<Nodes[s].inleaf<<") "<<t<<"("<<Nodes[t].isborder<<","<<Nodes[t].inleaf<<") "<<dis1<<" "<<Dijkstra(s,t,Nodes)<<endl;

                if(indexType==gtreeIndex){
                    GTreeO = GTree; //old GTree
                    NodesO = Nodes;//old Nodes
                    tt.start();
                    GtreeIndexUpdate(updates,ifParallel);
                    tt.stop();
                }else if(indexType==gstarIndex){
                    tt.start();
                    GstartreeIndexUpdate(updates,ifParallel);
                    tt.stop();
                }else if(indexType==lgtreeIndex){
                    tt.start();
                    LGTreeIndexUpdate(updates,ifParallel,DECREASE);
                    tt.stop();
                }
                else if(indexType==tgtreeIndex){
                    tt.start();
                    TGTreeIndexUpdate(updates,ifParallel,DECREASE);
                    tt.stop();
                }

                cout<<"update time: "<<tt.GetRuntime()<<" s."<<endl;
                ave_time += tt.GetRuntime();
                if(ifDebug){
                    CorrectnessCheck(100);
                }
            }
            cout<<"Average decrease update time: "<< ave_time/updateBatch <<" s."<<endl;
//            break;
            GTree=GTreeTemp; Nodes=NodesTemp; Tree=TreeTemp; OverlayGraph=OverlayGraphTemp; Label=LabelTemp; NeighborCon=NeighborConTemp;
        }
        case INCREASE:{
            cout<<"\nUpdate type: Increase"<<endl;
            ave_time = 0;
            for(int batch_i=0;batch_i<updateBatch;++batch_i){
                cout<<"Batch "<<batch_i<<" :";
                updates.clear();
                for(int i=batch_i;i<batch_i+updateVolume;++i) {//for each edge update
                    ID1 = updateEdges[i].first.first;
                    ID2 = updateEdges[i].first.second;
                    oldW = updateEdges[i].second;
                    newW = (int) (1.5 * oldW);//update
                    updates.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
                    if(ifDebug){
                        cout<<" "<<ID1<<"("<<Nodes[ID1].isborder<<","<<Nodes[ID1].inleaf<<") "<<ID2<<"("<<Nodes[ID2].isborder<<","<<Nodes[ID2].inleaf<<") "<<oldW<<" "<<newW;
                    }
                }
                if(ifDebug){
                    cout<<endl;
                }
//                cout<<"Before update: "<<endl;
//                int s=138682, t=146478;
//                int dis1=Query_TGTree(s,t);
//                cout<<s<<"("<<Nodes[s].isborder<<","<<Nodes[s].inleaf<<") "<<t<<"("<<Nodes[t].isborder<<","<<Nodes[t].inleaf<<") "<<dis1<<" "<<Dijkstra(s,t,Nodes)<<endl;

                if(indexType==gtreeIndex){
                    GTreeO = GTree; //old GTree
                    NodesO = Nodes;//old Nodes
                    tt.start();
                    GtreeIndexUpdate(updates,ifParallel);
                    tt.stop();
                }else if(indexType==gstarIndex){
                    tt.start();
                    GstartreeIndexUpdate(updates,ifParallel);
                    tt.stop();
                }else if(indexType==lgtreeIndex){
                    tt.start();
                    LGTreeIndexUpdate(updates,ifParallel,INCREASE);
                    tt.stop();
                }
                else if(indexType==tgtreeIndex){
                    tt.start();
                    TGTreeIndexUpdate(updates,ifParallel,INCREASE);
                    tt.stop();
                }

                cout<<"Update time: "<<tt.GetRuntime()<<" s."<<endl;
                ave_time += tt.GetRuntime();
                if(ifDebug){
                    CorrectnessCheck(100);
                }
            }
            cout<<"Average increase update time: "<< ave_time/updateBatch <<" s."<<endl;
            break;
        }
        case MIX:{
            cout<<"\nUpdate type: Random"<<endl;
            ave_time = 0;
            for(int batch_i=0;batch_i<updateBatch;++batch_i){
                cout<<"Batch "<<batch_i;

                updates.clear();
                for(int i=batch_i;i<batch_i+updateVolume;++i) {//for each edge update
                    ID1 = updateEdges[i].first.first;
                    ID2 = updateEdges[i].first.second;
                    oldW = updateEdges[i].second;
                    newW = (int) (2 * oldW * (rand() / double(RAND_MAX)));
                    updates.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
                }
                GTreeO = GTree; //old GTree
                NodesO = Nodes;//old Nodes
                tt.start();
                if(indexType==gtreeIndex){
                    GtreeIndexUpdate(updates,ifParallel);
                }else if(indexType==gstarIndex){
                    GstartreeIndexUpdate(updates,ifParallel);
                }
                tt.stop();
                cout<<" ; update time: "<<tt.GetRuntime()<<" s."<<endl;
                ave_time += tt.GetRuntime();
                if(ifDebug){
                    CorrectnessCheck(100);
                }
            }
            cout<<"The average time used for updating: "<< ave_time/updateBatch <<" s."<<endl;
            break;
        }
        default:
        {
            cout<<"Wrong update type! "<<updateType<<endl; exit(1);
        }
    }
}

//function of dealing with batch edge increase update of G-tree
void Gstartree::GstartreeIndexUpdate(vector<pair<pair<int,int>,pair<int,int> > > & updates, bool ifParallel){//update(<ID1,ID2>,oldW) vector<pair<pair<int,int>,int> > & updates,

//    double ave_time = 0;
    int ID1,ID2,oldW,newW;
    vector<pair<pair<int,int>,pair<int,int> > > uWithinLeafNode;// edge update within the same leaf node (<ID1,ID2>,<oldW,newW>)
    vector<pair<pair<int,int>,pair<int,int> > > uCrossLeafNode;// edge update among borders (<ID1,ID2>,<oldW,newW>)
    vector<pair<pair<int,int>,pair<int,int> > > uAmongBorders;// edge/shortcut update among borders (<ID1,ID2>,<oldW,newW>)
    bool flag_borderUpdate = false;//whether there is border update of leaf node
//    Timer tt;
//    tt.start();
    bool ifOnlyShortcut=false;
//    ifOnlyShortcut=true;

    if(!ifOnlyShortcut){
        for(int i=0;i<updates.size();++i){//for each edge update
            ID1 = updates[i].first.first;
            ID2 = updates[i].first.second;
            oldW = updates[i].second.first;
            newW = updates[i].second.second;

//        cout<<ID1 << " "<<ID2<<" ("<<oldW<<"->"<<newW<<") "<<endl;
            //update edge weight
            for(int j=0;j<Nodes[ID1].adjnodes.size();j++){
                if(Nodes[ID1].adjnodes[j]==ID2){
                    assert(Nodes[ID1].adjweight[j] == oldW);
                    Nodes[ID1].adjweight[j]=newW;
                    break;
                }
            }
            for(int j=0;j<Nodes[ID2].adjnodes.size();j++){
                if(Nodes[ID2].adjnodes[j]==ID1){
                    assert(Nodes[ID2].adjweight[j] == oldW);
                    Nodes[ID2].adjweight[j]=newW;
                    break;
                }
            }
            //identify the edge type
            if(Nodes[ID1].inleaf == Nodes[ID2].inleaf){//if it is an update within the leaf node
                uWithinLeafNode.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
            }else{//if it is an update among the leaf nodes
                uCrossLeafNode.emplace_back(make_pair(ID1,ID2),make_pair(oldW,newW));
                uAmongBorders.emplace_back(make_pair(ID1,ID2),make_pair(oldW,newW));
                flag_borderUpdate = true;
            }
        }
//        cout<<"After update: "<<Dijkstra(759,760,Nodes)<<endl;
        //// Method 1: directly reconstruct
//        for(int i=0;i<GTree.size();++i){
//            GTree[i].mind.clear();
//        }
//        hierarchy_shortest_path_calculation(ifParallel);

        //// Method 2: process updates within leaf node with pruning strategy, and reconstruct non-leaf nodes
        /// update processing
        int lnID;
//        int ID1_pos,ID2_pos;
        vector<Node> graph;//temp graph
        graph = Nodes;//original graph

        bool borderUpdate = false;
        int nleafnode;
        unordered_set<int> nodesToCheck;//the set of leaf nodes for update checking
        nodesToCheck.clear();
        unordered_set<int> leafsToCheck;//the set of leaf nodes for update checking
        leafsToCheck.clear();
        unordered_set<int> updatedNodes;//the set of all updated nodes
        updatedNodes.clear();
        unordered_set<int> borderAffectedNodes;
        borderAffectedNodes.clear();
        vector<int> cands;

        /// Update the leaf nodes
        Timer tt1;
        tt1.start();
        // for update within the same leaf node
        if(!uWithinLeafNode.empty()){
            for(auto it=uWithinLeafNode.begin();it!=uWithinLeafNode.end();++it){//deal with each update within leaf node
                ID1 = it->first.first; ID2 = it->first.second; oldW = it->second.first; newW = it->second.second;
                assert(Nodes[ID1].inleaf == Nodes[ID2].inleaf);
                lnID = Nodes[ID1].inleaf;//get the leaf node id
                /// insert the pending update leaf node
                if(leafsToCheck.find(lnID) == leafsToCheck.end()){//if not found
                    leafsToCheck.insert(lnID);
                }
            }
        }
        // for update among leaf nodes, a naive idea is that we can identify the affected area firstly, and propagate updates to other area (leaf nodes) if the border of this area is updated.
        int lca = -1;
        bool flag_cross = false;
        if(!uCrossLeafNode.empty()){
            flag_cross = true;
            vector<int> path_a, path_b;
            ///identify the affected area caused by cross-leaf edge updates
            for(int i=0;i<uCrossLeafNode.size();++i){
                ID1 = uCrossLeafNode[i].first.first; ID2 = uCrossLeafNode[i].first.second; oldW = uCrossLeafNode[i].second.first; newW = uCrossLeafNode[i].second.second;
                path_a = Nodes[ID1].gtreepath, path_b = Nodes[ID2].gtreepath;
                lca = path_a[find_LCA_pos(ID1,ID2)];
                //identify the affected area
                if(GTree[lca].isleaf){//if it is leaf node
                    if(leafsToCheck.find(lca) == leafsToCheck.end()){//if not found
                        leafsToCheck.insert(lca);
                    }
                }else{//if it is non-leaf node
                    nodesToCheck.insert(lca);
                    GetChildNodes(lca, nodesToCheck, leafsToCheck, updatedNodes);
                }

            }
        }
        // update leaf nodes
        borderAffectedNodes.clear();
        LeafLevelUpdate(graph,leafsToCheck,updatedNodes, borderAffectedNodes, ifParallel);//incorrect for parallel version
        tt1.stop();
//        cout<<"The time of leaf nodes updating: "<<tt1.GetRuntime()<<" s."<<endl;

        if(flag_cross){
//            cout<<nodesToCheck.size()<<endl;
            assert(!nodesToCheck.empty());
        }

        /// Update non-leaf nodes
        tt1.start();
        if(!borderAffectedNodes.empty() || !nodesToCheck.empty()){
            hierarchy_shortest_path_calculate(ifParallel);//correct. ifParallel
//            UpdateUpwardsPropagate(graph, updatedNodes,borderAffectedNodes, ifParallel, lca, nodesToCheck);//incorrect
        }
        tt1.stop();
//        cout<<"The time used for propagating updates to upper levels: "<<tt1.GetRuntime()<<" s."<<endl;

//        cout<<"Done."<<endl;
//    tt.stop();
//    cout<<"The time used for updating: "<<tt.GetRuntime()<<" s."<<endl;
//    ave_time += tt.GetRuntime();

    }


    if(indexType==gstarIndex){
        int idx_i, idx_j;
        multimap<int,pair<int,int> > shortcut_node_pairs;
        for (int i=0;i<leafNodePairs.size();++i) {
            idx_i = leafNodePairs[i].first;
            idx_j = leafNodePairs[i].second;
            shortcut_node_pairs.insert({compute_lca_level(idx_i,idx_j).first,make_pair(idx_i,idx_j)});
        }

        int top_level_;
        for(auto it=shortcut_node_pairs.begin();it!=shortcut_node_pairs.end();++it){
            top_level_ = it->first;
            idx_i = it->second.first; idx_j = it->second.second;
//        if(top_level_ == 1)
//            continue;
//        if(idx_i == 1759 && idx_j == 2176){
            build_shortcuts_for_nodes(idx_i, idx_j, top_level_, top_level_);
//        }
        }
    }
}


// load distance matrix from file
void Gstartree::load_graphOrdering(string filename){
//    char filename[300];
//    if(dataset == "cal"){
//        strcpy(filename, DataPath.c_str());
//        strcat(filename, FILE_ONTREE_MIND.c_str());
//    }else {
//        strcpy(filename, dirname.c_str());
//        strcat(filename, "/");
//        strcat(filename, FILE_ONTREE_MIND.c_str());
//    }
    FILE* fin = fopen( filename.c_str(), "rb" );
    if(fin == NULL){
        cout << "Failed to open file " << filename << endl;
        exit(1);
    }
    cout << "Loading overlay graph information... "<<filename<<endl;
    int* buf;
    int count, pos = 0;
    int ID1, ID2, weight;
    /// Overlay graph
    overlayNodeNum=0;
    OverlayGraph.assign(node_num,unordered_map<int,int>());
    for(ID1=0;ID1<OverlayGraph.size();++ID1){
        fread( &count, sizeof(int), 1, fin );
        if(count>=1){
            overlayNodeNum++;
            for(int j=0;j<count;++j){
                fread( &ID2, sizeof(int), 1, fin );
                fread( &weight, sizeof(int), 1, fin );
                OverlayGraph[ID1].insert({ID2,weight});
            }
        }

    }
    /// Ordering
    buf = new int[node_num];
    NodeOrder.clear();
    fread( buf, sizeof(int), node_num, fin );
    for ( int i = 0; i < node_num; i++ ){
        NodeOrder.push_back(buf[i]);
    }
    delete[] buf;
    fclose(fin);
    vNodeOrder.assign(node_num,-1);
    for(int i=0;i<node_num;++i){
        vNodeOrder[NodeOrder[i]]=i;
    }
//    cout<<"Done."<<endl;
    cout<<"Overlay graph size: "<<overlayNodeNum<<endl;
}

void Gstartree::load_shortcuts(string filename) {
    FILE *fin = fopen(filename.c_str(), "rb");
//    FILE *fin = fopen(FILE_SHORTCUT, "rb");
    if(fin == NULL){
        cout << "Failed to open file " << filename << endl;
        exit(1);
    }
    cout<<"Loading shortcuts... ";
    int *buf;
    int count;
    int shortcutPairNum=0;
    int leafNodePairNum;
    leafNodePairs.clear();

    fread(&leafNodePairNum, sizeof(int), 1, fin);//shortcuts.size()
    int kv1_size, kv2_size;
    fread(&kv1_size, sizeof(int), 1, fin);//shortcuts.size()

    int i, j;

    for (int x = 0; x < kv1_size; ++x) {
        fread(&i, sizeof(int), 1, fin);//it.first, n1
        fread(&kv2_size, sizeof(int), 1, fin);//it.second.size()

        for (int y = 0; y < kv2_size; ++y) {
            fread(&j, sizeof(int), 1, fin);//it.second.first, n2
            fread(&count, sizeof(int), 1, fin);//it.second.second.size()

            buf = new int[count];
            fread(buf, sizeof(int), count, fin);
            shortcuts[i][j].clear();
            for(int k=0;k<count;++k){
                shortcuts[i][j].push_back(buf[k]);
            }
//            shortcuts[i][j].insert(shortcuts[i][j].end(), buf, buf + count);
//            shortcuts[j][i].clear();
//            shortcuts[j][i].emplace_back(-1);
            ++shortcutPairNum;
            leafNodePairs.emplace_back(i,j);
            delete[] buf;
//            ///check correctness of shortcuts
//            int pb=GTree[i].borders.size();
//            int qb=GTree[j].borders.size();
//            int temp_dis;
//            for(int p=0;p<GTree[i].borders.size();++p){
//                for(int q=0;q<GTree[j].borders.size();++q){
//                    temp_dis = Dijkstra(GTree[i].borders[p],GTree[j].borders[q]);
////                    cout<<p*qb+q<<endl;
//                    if(shortcuts[i][j][p*qb+q] != temp_dis){
//                        cout << GTree[i].borders[p] << " "<< GTree[j].borders[q] << " " << shortcuts[i][j][p*qb+q] << " " <<  temp_dis << endl;
//                    }
//
//                }
//            }
        }
    }
    fclose(fin);
    cout <<"Shortcut pair number: "<<shortcutPairNum<<" "<<leafNodePairNum<<" "<<leafNodePairs.size()<<endl;

}

// up_pos & current_pos(used for quickly locating parent & child nodes)
/*void Gstartree::build_up_and_down_pos() {
    unordered_map<int, int> pos_map;
    for (int i = 1; i < GTree.size(); i++) {//for each tree node
        GTree[i].current_pos.clear();
        GTree[i].up_pos.clear();

        // current_pos
        pos_map.clear();
        for (int j = 0; j < GTree[i].union_borders.size(); j++) {
            pos_map[GTree[i].union_borders[j]] = j;//from vertex id to position id
        }
        for (int j = 0; j < GTree[i].borders.size(); j++) {
            GTree[i].current_pos.push_back(pos_map[GTree[i].borders[j]]);//get the borders' position id in this node
        }
        // up_pos
        pos_map.clear();
        for (int j = 0; j < GTree[GTree[i].father].union_borders.size(); j++) {
            pos_map[GTree[GTree[i].father].union_borders[j]] = j;//map the vertex id of parent's borders to position id
        }
        for (int j = 0; j < GTree[i].borders.size(); j++) {
            GTree[i].up_pos.push_back(pos_map[GTree[i].borders[j]]);// get the borders' position id in the parent node
        }
    }
}*/

//old function of computing the score of shortcut pair
double Gstartree::compute_value_of_leaf_node_pair(int i, int j) {//inline
    double w = double(GTree[i].leafnodes.size()) * GTree[j].leafnodes.size() / (GTree[i].borders.size() * GTree[j].borders.size());// larger value indicates denser area
    double v = GTree[i].borders.size() + GTree[j].borders.size();//v is function mu

    int father;

    while ((father = GTree[i].father) != 0) {//
        v += GTree[father].borders.size() * GTree[i].borders.size();
        i = father;
    }

    while ((father = GTree[j].father) != 0) {
        v += GTree[father].borders.size() * GTree[j].borders.size();
        j = father;
    }

    v += GTree[i].borders.size() * GTree[j].borders.size();

    return v * w;
}
//new function of computing the score of shortcut pair: correct one
/*double Gstartree::compute_value_of_leaf_node_pair(int i, int j) {//inline
    double w = double(GTree[i].leafnodes.size()) * GTree[j].leafnodes.size() / (GTree[i].borders.size() * GTree[j].borders.size());// larger value indicates denser area
    double v = GTree[i].borders.size() + GTree[j].borders.size();//v is function mu

    int father,lca=0;

    vector<int> treePath_i, treePath_j;
    int n1=i;
    int n2=j;
    while ((father = GTree[n1].father) != 0) {//
        treePath_i.push_back(father);
        n1=father;
    }
    while ((father = GTree[n2].father) != 0) {
        treePath_j.push_back(father);
        n2=father;
    }
    for(int p=0;p<treePath_i.size();++p){
        bool flag_find = false;
        for(int q=0;q<treePath_j.size();++q){
            if(treePath_i[p] == treePath_j[q]){
                lca = treePath_i[p];
                flag_find = true;
                break;
            }
        }
        if(flag_find)
            break;
    }
    n1=i;
    while ((father = GTree[n1].father) != lca) {//
        v += GTree[father].borders.size() * GTree[n1].borders.size();
        n1 = father;
    }
    n2=j;
    while ((father = GTree[n2].father) != lca) {
        v += GTree[father].borders.size() * GTree[n2].borders.size();
        n2 = father;
    }

    v += GTree[n1].borders.size() * GTree[n2].borders.size();

    return v * w;
}*/
//function to check whether two leaf nodes are adjacent
bool Gstartree::check_leaf_adjacent(TreeNode &li, TreeNode &lj) {

    for (const auto &bi : li.borders) {
        for (const auto &bj : lj.borders) {
            if (find(Nodes[bi].adjnodes.begin(), Nodes[bi].adjnodes.end(), bj) != Nodes[bi].adjnodes.end()) {//if li and lj have adjacent border vertex
                return true;
            }
        }
    }
    return false;
}
pair<int,int> Gstartree::compute_lca_level(int x, int y){
    const auto &path_x = Nodes[GTree[x].leafnodes[0]].gtreepath;
    const auto &path_y = Nodes[GTree[y].leafnodes[0]].gtreepath;

    /// Find LCA
    int top_level_x=1,top_level_y=1;
//    if(path_x.size() != path_y.size()){
//        cout<<"!"<<endl;
//    }
//    assert(path_x.size() == path_y.size());
    for(int i=0;i<path_x.size();++i){
        for(int j=0;j<path_y.size();++j){
            if(path_x[i]==path_y[j]){
                top_level_x=i+1;
                top_level_y=j+1;
            }
        }
    }
    if(top_level_x!=top_level_y){
        cout<<"Error in LCA level checking."<<endl;
        exit(1);
    }
    return make_pair(top_level_x,top_level_y);
}
//function of building shortcuts for given node pair, original version, incorrect
void Gstartree::build_shortcuts_for_nodes(int x, int y) {
    const auto &path_x = Nodes[GTree[x].leafnodes[0]].gtreepath;
    const auto &path_y = Nodes[GTree[y].leafnodes[0]].gtreepath;

    if (!GTree[x].is_cached) {
        const auto &path = path_x;
        int tn, cn, min, dist, posa, posb, posz, poszz;
        unsigned long union_border_size;

        for (int z = 0; z < GTree[path[1]].borders.size(); ++z) {

            if (!GTree[path[1]].is_cached) {
                posz = GTree[path[1]].up_pos[z];
                for (int zz = 0; zz < GTree[path[1]].borders.size(); ++zz) {
                    poszz = GTree[path[1]].up_pos[zz];
                    GTree[path[1]].cacheB[z][zz] = GTree[0].mind[posz * GTree[0].union_borders.size() + poszz];//the distance among the borders of GTree[path[1]]
                }
            }

            for (int i = 1; i < path.size() - 1; ++i) {
                tn = path[i];
                cn = path[i + 1];

                if (!GTree[cn].is_cached) {
                    for (int j = 0; j < GTree[cn].borders.size(); ++j) {
                        union_border_size = GTree[tn].union_borders.size();
                        min = INT_MAX;
                        posa = GTree[cn].up_pos[j];
                        for (int k = 0; k < GTree[tn].borders.size(); ++k) {
                            posb = GTree[tn].current_pos[k];
                            dist = GTree[tn].cacheB[z][k] + GTree[tn].mind[posb * union_border_size + posa];
                            if (dist < min) min = dist;
                        }
                        GTree[cn].cacheB[z][j] = min;
                    }
                }
            }
        }

        for (int i = 1; i < path.size(); ++i) {
            GTree[path[i]].is_cached = true;
        }
    }

    if (!GTree[y].is_cached) {
        const auto &path = path_y;
        int tn, cn, min, dist, posa, posb, posz, poszz;
        unsigned long union_border_size;

        for (int z = 0; z < GTree[path[1]].borders.size(); ++z) {

            if (!GTree[path[1]].is_cached) {
                posz = GTree[path[1]].up_pos[z];
                for (int zz = 0; zz < GTree[path[1]].borders.size(); ++zz) {
                    poszz = GTree[path[1]].up_pos[zz];
                    GTree[path[1]].cacheB[z][zz] = GTree[0].mind[posz * GTree[0].union_borders.size() + poszz];
                }
            }

            for (int i = 1; i < path.size() - 1; ++i) {
                tn = path[i];
                cn = path[i + 1];

                if (!GTree[cn].is_cached) {
                    for (int j = 0; j < GTree[cn].borders.size(); ++j) {
                        union_border_size = GTree[tn].union_borders.size();
                        min = INT_MAX;
                        posa = GTree[cn].up_pos[j];
                        for (int k = 0; k < GTree[tn].borders.size(); ++k) {
                            posb = GTree[tn].current_pos[k];
                            dist = GTree[tn].cacheB[z][k] + GTree[tn].mind[posb * union_border_size + posa];
                            if (dist < min) min = dist;
                        }
                        GTree[cn].cacheB[z][j] = min;
                    }
                }
            }
        }

        for (int i = 1; i < path.size(); ++i) {
            GTree[path[i]].is_cached = true;
        }
    }


    int ns_top = path_x[1];
    int nt_top = path_y[1];
    int ns_top_up_pos, nt_top_up_pos;
    auto union_border_size = GTree[0].union_borders.size();
    int min, dist;
    vector<vector<int>> mid(GTree[ns_top].borders.size(), vector<int>(GTree[nt_top].borders.size(), 0));

    for (int i = 0; i < GTree[ns_top].borders.size(); i++) {
        ns_top_up_pos = GTree[ns_top].up_pos[i];
        for (int j = 0; j < GTree[nt_top].borders.size(); j++) {
            nt_top_up_pos = GTree[nt_top].up_pos[j];
            mid[i][j] = GTree[0].mind[ns_top_up_pos * union_border_size + nt_top_up_pos];
        }
    }

    for (int i = 0; i < GTree[x].borders.size(); ++i) {
        for (int j = 0; j < GTree[y].borders.size(); ++j) {
            min = INT_MAX;
            for (int z = 0; z < GTree[ns_top].borders.size(); ++z) {
                for (int zz = 0; zz < GTree[nt_top].borders.size(); ++zz) {
                    dist = GTree[x].cacheB[z][i] + GTree[y].cacheB[zz][j] + mid[z][zz];
                    if (dist < min) {
                        min = dist;
                    }
                }
            }
            int temp_dist=Dijkstra(GTree[x].borders[i],GTree[y].borders[j],Nodes);
            if(min != temp_dist){
                cout << GTree[x].borders[i] << " "<< GTree[y].borders[j] << " " << min << " " << temp_dist << endl; exit(1);
            }
            shortcuts[x][y].emplace_back(min);
        }
    }

}

void Gstartree::build_shortcuts_for_nodes(int x, int y,int top_level_x,int top_level_y) {//inline
    const auto &path_x = Nodes[GTree[x].leafnodes[0]].gtreepath;
    const auto &path_y = Nodes[GTree[y].leafnodes[0]].gtreepath;

    /// Find LCA
//    int top_level_x=1,top_level_y=1;
//    top_level_x = compute_lca_level(x,y);
//    top_level_y = top_level_x;

    /// update cache flag
    if(top_level_x != top_level || top_level_y != top_level){
//        cout<<"!"<<top_level_x<<" " <<top_level_y << endl;
        for (int i = 0; i < GTree.size(); ++i) {
            GTree[i].is_cached = false;
            int j = i, father=j;
            vector<int> fathers;
            while (father > 0){
                fathers.emplace_back(father);
                father = GTree[father].father;
            }
            if(int(fathers.size())-top_level_x >= 0){
                j = fathers[fathers.size()-top_level_x];
                GTree[i].cacheB.resize(GTree[j].borders.size(), vector<int>(GTree[i].borders.size(), 0));
            }
        }
        top_level = top_level_x;
    }
//    if(top_level > 1){
//        cout<<x<<" "<<y<<endl;
//    }

    if (!GTree[x].is_cached) {//if not cached
        const auto &path = path_x;
        int tn, cn, min, dist, posa, posb, posz, poszz;
        unsigned long union_border_size;

        for (int z = 0; z < GTree[path[top_level_x]].borders.size(); ++z) {//for each border in level-1 tree node, compute the dis to all borders of children nodes

            if (!GTree[path[top_level_x]].is_cached) {
                posz = GTree[path[top_level_x]].up_pos[z];//position id in the root node
                for (int zz = 0; zz < GTree[path[top_level_x]].borders.size(); ++zz) {
                    poszz = GTree[path[top_level_x]].up_pos[zz];
                    GTree[path[top_level_x]].cacheB[z][zz] = GTree[path[top_level_x-1]].mind[posz * GTree[path[top_level_x-1]].union_borders.size() + poszz];//get dis from each border to other borders in level-1 node
//                    int temp_dis = Dijkstra(GTree[path[top_level_x]].borders[z],GTree[path[top_level_x]].borders[zz]);
//                    if(GTree[path[top_level_x]].cache[z][zz]!=temp_dis){
//                        cout << GTree[path[top_level_x]].borders[z] << " "<< GTree[path[top_level_x]].borders[zz] << " " << GTree[path[top_level_x]].cache[z][zz] << " " << temp_dis << endl;
//                    }
                }
            }

            for (int i = top_level_x; i < path.size() - 1; ++i) {
                tn = path[i];
                cn = path[i + 1];

                if (!GTree[cn].is_cached) {
                    for (int j = 0; j < GTree[cn].borders.size(); ++j) {//for each border in cn
                        union_border_size = GTree[tn].union_borders.size();
                        min = INT_MAX;
                        posa = GTree[cn].up_pos[j];//position id in parent node (tn)

                        for (int k = 0; k < GTree[tn].borders.size(); ++k) {//dis from borders of level cn to borders of level tn
                            posb = GTree[tn].current_pos[k];
                            dist = GTree[tn].cacheB[z][k] + GTree[tn].mind[posb * union_border_size + posa];
                            if (dist < min)
                                min = dist;
                        }

                        GTree[cn].cacheB[z][j] = min;
//                        int temp_dis = Dijkstra(GTree[path[top_level_x]].borders[z],GTree[cn].borders[j]);
//                        if(min!=temp_dis){
//                            cout << GTree[path[top_level_x]].borders[z] << " "<< GTree[cn].borders[j] << " " << min << " " << temp_dis << endl;
//                        }
                    }
                }
            }
        }

        for (int i = top_level_x; i < path.size(); ++i) {
            GTree[path[i]].is_cached = true;
        }
    }

    if (!GTree[y].is_cached) {
        const auto &path = path_y;
        int tn, cn, min, dist, posa, posb, posz, poszz;
        unsigned long union_border_size;

        for (int z = 0; z < GTree[path[top_level_y]].borders.size(); ++z) {//compute distance from each border of level-1 node to all borders of its children nodes

            if (!GTree[path[top_level_y]].is_cached) {
                posz = GTree[path[top_level_y]].up_pos[z];
                for (int zz = 0; zz < GTree[path[top_level_y]].borders.size(); ++zz) {
                    poszz = GTree[path[top_level_y]].up_pos[zz];
                    GTree[path[top_level_y]].cacheB[z][zz] = GTree[path[top_level_y-1]].mind[posz * GTree[path[top_level_y-1]].union_borders.size() + poszz];
//                    int temp_dis = Dijkstra(GTree[path[top_level_y]].borders[z],GTree[path[top_level_y]].borders[zz]);
//                    if(GTree[path[top_level_y]].cache[z][zz]!=temp_dis){
//                        cout << GTree[path[top_level_y]].borders[z] << " "<< GTree[path[top_level_y]].borders[zz] << " " << GTree[path[top_level_y]].cache[z][zz] << " " << temp_dis << endl;
//                    }
                }
            }

            for (int i = top_level_y; i < path.size() - 1; ++i) {
                tn = path[i];
                cn = path[i + 1];

                if (!GTree[cn].is_cached) {
                    for (int j = 0; j < GTree[cn].borders.size(); ++j) {
                        union_border_size = GTree[tn].union_borders.size();
                        min = INT_MAX;
                        posa = GTree[cn].up_pos[j];
                        for (int k = 0; k < GTree[tn].borders.size(); ++k) {
                            posb = GTree[tn].current_pos[k];
                            dist = GTree[tn].cacheB[z][k] + GTree[tn].mind[posb * union_border_size + posa];
                            if (dist < min)
                                min = dist;
                        }
                        GTree[cn].cacheB[z][j] = min;
//                        int temp_dis = Dijkstra(GTree[path[top_level_y]].borders[z],GTree[cn].borders[j]);
//                        if(min!=temp_dis){
//                            cout << GTree[path[top_level_y]].borders[z] << " "<< GTree[cn].borders[j] << " " << min << " " << temp_dis << endl;
//                        }
                    }
                }
            }
        }

        for (int i = top_level_y; i < path.size(); ++i) {
            GTree[path[i]].is_cached = true;
        }
    }


    int ns_top = path_x[top_level_x];
    int nt_top = path_y[top_level_y];
    int ns_top_up_pos, nt_top_up_pos;
    auto union_border_size = GTree[path_x[top_level_x-1]].union_borders.size();
    int min, dist;
    vector<vector<int>> mid(GTree[ns_top].borders.size(), vector<int>(GTree[nt_top].borders.size(), 0));

    for (int i = 0; i < GTree[ns_top].borders.size(); i++) {//compute the cross-branch distance
        ns_top_up_pos = GTree[ns_top].up_pos[i];
        for (int j = 0; j < GTree[nt_top].borders.size(); j++) {
            nt_top_up_pos = GTree[nt_top].up_pos[j];
            mid[i][j] = GTree[path_x[top_level_x-1]].mind[ns_top_up_pos * union_border_size + nt_top_up_pos];
//            int temp_dis = Dijkstra(GTree[ns_top].borders[i],GTree[nt_top].borders[j]);
//            if(mid[i][j]!=temp_dis){
//                cout << GTree[ns_top].borders[i] << " "<< GTree[nt_top].borders[j] << " " << mid[i][j] << " " << temp_dis << endl;
//            }
        }
    }
    /// compute final distance
    for (int i = 0; i < GTree[x].borders.size(); ++i) {//for each border in leaf node x
        for (int j = 0; j < GTree[y].borders.size(); ++j) {//for each border in leaf node y
            min = INT_MAX;
            for (int z = 0; z < GTree[ns_top].borders.size(); ++z) {//for each border in top-level
                for (int zz = 0; zz < GTree[nt_top].borders.size(); ++zz) {
                    dist = GTree[x].cacheB[z][i] + GTree[y].cacheB[zz][j] + mid[z][zz];
                    if (dist < min) {
                        min = dist;
                    }
                }
            }
//            int temp_dis = Dijkstra(GTree[x].borders[i],GTree[y].borders[j],Nodes);
//            if(min!=temp_dis){
//                cout << GTree[x].borders[i] << " "<< GTree[y].borders[j] << " " << min << " " << temp_dis << endl;
//            }
            shortcuts[x][y].emplace_back(min);
        }
    }

}
//function of building shortcuts
void Gstartree::build_shortcuts() {//inline

    priority_queue<leaf_node_pair> leaf_node_pairs;

    const auto n = double(SHORTCUT_THRESHOLD);//eta, the shortcut threshold
    double new_value;

    long count = 0;
    double value_sum = 0.0;

//    for (int i = 0; i < leaf_nodes.size() - 1; ++i) {
//        for (int j = i + 1; j < leaf_nodes.size(); ++j) {
//            if (GTree[leaf_nodes[i]].father != GTree[leaf_nodes[j]].father &&
//                check_leaf_adjacent(GTree[leaf_nodes[i]], GTree[leaf_nodes[j]])) {
//                new_value = compute_value_of_leaf_node_pair(leaf_nodes[i], leaf_nodes[j]);
//                if (leaf_node_pairs.size() < n) {
//                    leaf_node_pairs.emplace(leaf_nodes[i], leaf_nodes[j], new_value);
//                } else if (leaf_node_pairs.top().v < new_value) {
//                    leaf_node_pairs.pop();
//                    leaf_node_pairs.emplace(leaf_nodes[i], leaf_nodes[j], new_value);
//                }
//            }
//            printf("+%ld\n", count++);
//        }
//    }
    assert(!leaf_nodes.empty());
    for (int i = 0; i < leaf_nodes.size() - 1; ++i) {//for each leaf node
        for (int j = i + 1; j < leaf_nodes.size(); ++j) {
            if (GTree[leaf_nodes[i]].father != GTree[leaf_nodes[j]].father &&
                check_leaf_adjacent(GTree[leaf_nodes[i]], GTree[leaf_nodes[j]])) {//if they do not have the same father but adjacent in the graph
                new_value = compute_value_of_leaf_node_pair(leaf_nodes[i], leaf_nodes[j]);//compute mu, larger value mean higher cost
                if (value_sum < n) {
                    leaf_node_pairs.emplace(leaf_nodes[i], leaf_nodes[j], new_value);
                    value_sum += (GTree[leaf_nodes[i]].borders.size() * GTree[leaf_nodes[j]].borders.size());
                } else if (leaf_node_pairs.top().v < new_value) {//if
                    leaf_node_pairs.pop();
                    leaf_node_pairs.emplace(leaf_nodes[i], leaf_nodes[j], new_value);
                }
            }
        }
    }

    count = leaf_node_pairs.size();

    int idx_i, idx_j;
    multimap<int,pair<int,int> > shortcut_node_pairs;
    while (!leaf_node_pairs.empty()) {
        idx_i = leaf_node_pairs.top().x;
        idx_j = leaf_node_pairs.top().y;
        shortcut_node_pairs.insert({compute_lca_level(idx_i,idx_j).first,make_pair(idx_i,idx_j)});
//        cout << "x = " << idx_i << "; y = " << idx_j << endl;
//        build_shortcuts_for_nodes(idx_i, idx_j);//incorrect
        leafNodePairs.emplace_back(idx_i,idx_j);
        leaf_node_pairs.pop();
    }

    int top_level_;
    for(auto it=shortcut_node_pairs.begin();it!=shortcut_node_pairs.end();++it){
        top_level_ = it->first;
        idx_i = it->second.first; idx_j = it->second.second;
//        if(top_level_ == 1)
//            continue;
//        if(idx_i == 1759 && idx_j == 2176){
            build_shortcuts_for_nodes(idx_i, idx_j, top_level_, top_level_);
//        }
    }
    cout << "Shortcut number of leaf-node pair: " << count << " " <<leafNodePairs.size()<<endl;
}

void Gstartree::save_shortcuts(string FS) {
    cout << "Begin saving shortcuts to file: " << FS << endl;
    FILE *fout = fopen(FS.c_str(), "wb");
    if(fout == NULL){
        cout << "Failed to open file " << FS << endl;
        exit(1);
    }

    int *buf;
    unsigned long count;

    int temp=leafNodePairs.size();
    assert(temp>0);
    fwrite(&temp, sizeof(int), 1, fout);//shortcut pair number

    int kv1_size = shortcuts.size();
    fwrite(&kv1_size, sizeof(int), 1, fout);//shortcut size

    for (auto &kv1 : shortcuts) {
        int i = kv1.first;
        fwrite(&i, sizeof(int), 1, fout);//source node

        int kv2_size = kv1.second.size();
        fwrite(&kv2_size, sizeof(int), 1, fout);//size

        for (auto &kv2 : kv1.second) {
            auto j = kv2.first;
            fwrite(&j, sizeof(int), 1, fout);//target node

            count = kv2.second.size();
            fwrite(&count, sizeof(int), 1, fout);//size

            buf = new int[count];
            copy(kv2.second.begin(), kv2.second.end(), buf);//shortcuts
            fwrite(buf, sizeof(int), count, fout);

            delete[] buf;
        }
    }

    fclose(fout);
}

int Gstartree::gstartree_build() {
    Timer tt;
    double t1,t2,t3;

    bool ifRead=true;
//    ifRead=false;

    string filename=dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".Gstar.minds.bin";
    FILE *fin = fopen( filename.c_str(), "rb" );
    if(fin == NULL){//if not exist
        ifRead= false;
//        cout << "Failed to open file " << filename << endl;
//        exit(1);
    }else{//if exist
        ifRead=true;
        fclose(fin);
    }

    if(ifRead){
        ReadGraph_W();
        load_gtree(dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".Gstar.gtree.bin", dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".Gstar.paths.bin");
        load_minds(dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".Gstar.minds.bin");
    }else{
        // init
        init();

        // gtree_build
        cout << "Start to build G-tree..."<<endl;
        tt.start();
        build();
//    cout <<"Done."<<endl;
        tt.stop();
        t1 = tt.GetRuntime();
        cout << "The time for G-tree building: " << t1 << " s." << endl;


        // dump gtree
        gtree_save(dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".Gstar.gtree.bin", dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".Gstar.paths.bin");
        // calculate distance matrix
        cout << "Start to calculate distance matrix..."<<endl;
        tt.start();
        hierarchy_shortest_path_calculation(ifParallel);
        tt.stop();
        t2 = tt.GetRuntime();
        cout << "The time for distance matrix computing: " << t2 << " s." << endl;
        hierarchy_shortest_path_save(dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".Gstar.minds.bin");
    }

    build_up_and_down_pos();


    if(dataset=="NY"){
        SHORTCUT_THRESHOLD=100000;
    }else if(dataset=="FLA"){
        SHORTCUT_THRESHOLD=150000;
    }else if(dataset=="W"){
        SHORTCUT_THRESHOLD=1100000;
    }else if(dataset=="USA"){
        SHORTCUT_THRESHOLD==45000000;
    }else{
        cout<<"Wrong dataset! "<<dataset<<endl; exit(1);
    }

    cout << "Shortcut threshold: "<< SHORTCUT_THRESHOLD << endl;
    cout << "Building shortcuts..." << endl;

    tt.start();
//    for(int i=0;i<GTree.size();++i){
//        int j = i, father;
//        while ((father = GTree[j].father) > 0) j = father;
//        GTree[i].cacheB.resize(GTree[j].borders.size(), vector<int>(GTree[i].borders.size(), 0));
//    }
    build_shortcuts();
    tt.stop();
    t3=tt.GetRuntime();
    cout << "The time for shortcut building: " << t3 << " s." << endl;
    cout << "Overall time for index construction: " << t1+t2+t3 << " s." << endl;

//    char filename[300];
//    if(dataset == "cal"){
//        strcpy(filename, DataPath.c_str());
//        strcat(filename, FILE_SHORTCUT.c_str());
//        strcat(filename, "-");
//        strcat(filename, std::to_string(percent).c_str());
//        strcat(filename, ".bin");
//    }else {
//        strcpy(filename, dirname.c_str());
//        strcat(filename, "/");
//        strcat(filename, FILE_SHORTCUT.c_str());
//        strcat(filename, "-");
//        strcat(filename, std::to_string(percent).c_str());
//        strcat(filename, ".bin");
//    }
    cout << "Saving shortcuts..." << endl;
    filename=dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".Gstar.sc_"+to_string(SHORTCUT_THRESHOLD)+".bin";
    save_shortcuts(filename);
    cout << "Done." << endl;

    return 0;
}

void Gstartree::init_ReadIndex() {
    if(dataset == "cal"){
        load_graph();
    }else{
        ReadGraph_W();
    }

    if(indexType==gtreeIndex){
        load_gtree(dirname+"/"+FILE_GTREE, dirname+"/"+FILE_NODES_GTREE_PATH);
        load_minds(dirname+"/"+FILE_ONTREE_MIND);
    }else if(indexType==gstarIndex){
        load_gtree(dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".Gstar.gtree.bin", dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".Gstar.paths.bin");
        load_minds(dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".Gstar.minds.bin");
        shortcuts.clear();
        if(dataset=="NY"){
            SHORTCUT_THRESHOLD=100000;
        }else if(dataset=="FLA"){
            SHORTCUT_THRESHOLD=150000;
        }else if(dataset=="W"){
            SHORTCUT_THRESHOLD=1100000;
        }else if(dataset=="USA"){
            SHORTCUT_THRESHOLD==45000000;
        }else{
            cout<<"Wrong dataset! "<<dataset<<endl; exit(1);
        }
        string filename=dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".Gstar.sc_"+to_string(SHORTCUT_THRESHOLD)+".bin";
        load_shortcuts(filename);

    }
    build_up_and_down_pos();//
}

void Gstartree::init_LGTreeIndex() {
    if(dataset == "cal"){
        load_graph();
    }else{
        ReadGraph_W();
    }

    load_gtree(dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".LG.gtree.bin", dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".LG.paths.bin");
    load_minds(dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".LG.minds.bin");
    load_graphOrdering(dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".LG.graph.bin");
    load_overlayIndex(dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".LG.PLL.bin");

    /// Construct L~
//    Timer tt;
//    if(ifHierarchicalOrdering){
//        tt.start();
//        H2HconOrderMT(false);///
//        tt.stop();
//    }
//    else{
//        tt.start();
//        H2HconOrderMT(true);///
//        tt.stop();
//    }
//    cout << "The time for overlay index building: " << tt.GetRuntime() << " s." << endl;
    // dump distance matrix
//    cout << "Saving distance matrix..."<<endl;
//    H2HIndexSave();
//    cout << "Done."<<endl;
    IndexSizePLL();
}

void Gstartree::init_TGTreeIndex() {
    if(dataset == "cal"){
        load_graph();
    }else{
        ReadGraph_W();
    }
    if(percentScale!=0){
        load_gtree(dirname+"/"+dataset+"_"+ to_string(percentScale) + "." + to_string(LEAF_CAP) + ".TD.gtree.bin", dirname+"/"+dataset+"_"+ to_string(percentScale) + "." + to_string(LEAF_CAP) + ".TD.paths.bin");
        load_minds(dirname+"/"+dataset +"_"+ to_string(percentScale) + "." + to_string(LEAF_CAP) + ".TD.minds.bin");
        load_graphOrdering(dirname+"/"+dataset+"_"+ to_string(percentScale) + "." + to_string(LEAF_CAP) + ".TD.graph.bin");
    }
    else{
        load_gtree(dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".TD.gtree.bin", dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".TD.paths.bin");
        load_minds(dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".TD.minds.bin");
        load_graphOrdering(dirname+"/"+dataset + "." + to_string(LEAF_CAP) + ".TD.graph.bin");
    }


    /// Construct L~
    Timer tt;
    if(ifHierarchicalOrdering){
        tt.start();
        H2HconOrderMT(false);///
        tt.stop();
    }
    else{
        tt.start();
        H2HconOrderMT(true);///
        tt.stop();
    }
    cout << "The time for overlay index building: " << tt.GetRuntime() << " s." << endl;
    // dump distance matrix
//    cout << "Saving distance matrix..."<<endl;
//    H2HIndexSave();
//    cout << "Done."<<endl;
    IndexSizeH2H();
}

//void Gstartree::init_query() {
//    for (auto &tn: GTree) {
//        //tn.oclist.clear();
//        tn.cache = vector<int>(tn.borders.size(), 0);
//        tn.is_visited = false;
//    }
//}

vector<int> Gstartree::load_objects() {
    unordered_set<int> o;
    o.clear();
    char filename[300];
    if(dataset == "cal"){
        strcpy(filename, DataPath.c_str());
        strcat(filename, FILE_OBJECT.c_str());
    }else {
        strcpy(filename, dirname.c_str());
        strcat(filename, "/");
        strcat(filename, FILE_OBJECT.c_str());
    }
    FILE *fin = fopen(filename, "r");
    if(fin == NULL){
        cout << "Failed to open file " << filename << endl;
        exit(1);
    }
    int oid;
    while (fscanf(fin, "%d", &oid) == 1) {
        o.insert(oid);
    }
    fclose(fin);

    vector<int> res(o.begin(), o.end());

    return res;
}

/*int Gstartree::dijkstra_p2p(int s, int t) {//inline
    unordered_map<int, int> result;
    result[s] = 0;

    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    pq.emplace(0, s);

    int min, minpos, adjnode, weight;
    TIME_TICK_START

    while (!pq.empty()) {
        min = pq.top().first;
        minpos = pq.top().second;

        if (minpos == t) {
            TIME_TICK_END
            return min;
        }

        pq.pop();

        for (int i = 0; i < Nodes[minpos].adjnodes.size(); i++) {

            adjnode = Nodes[minpos].adjnodes[i];
            weight = Nodes[minpos].adjweight[i];

            if (result.find(adjnode) == result.end() || result[adjnode] > min + weight) {
                result[adjnode] = min + weight;
                pq.emplace(min + weight, adjnode);
            }
        }
    }

    return -1;
}*/

/*vector<int> Gstartree::dijkstra_candidate(int s, unordered_set<int> &cands) {//inline
    // init
    auto num_cands = cands.size();
    unordered_set<int> todo(cands.begin(), cands.end());
    unordered_map<int, int> result;
    result[s] = 0;

    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    pq.emplace(0, s);

    int min, minpos, adjnode, weight;
    while (!todo.empty() && !pq.empty()) {
        min = pq.top().first;
        minpos = pq.top().second;
        pq.pop();

        todo.erase(minpos);

        for (int i = 0; i < Nodes[minpos].adjnodes.size(); i++) {

            adjnode = Nodes[minpos].adjnodes[i];
            weight = Nodes[minpos].adjweight[i];

            if (result.find(adjnode) == result.end() || result[adjnode] > min + weight) {
                result[adjnode] = min + weight;
                pq.emplace(min + weight, adjnode);
            }
        }
    }

    // output
    vector<int> output;
    output.reserve(num_cands);
    for (const auto &iter : cands) {
        output.emplace_back(result[iter]);
    }

    // return
    return output;

}*/

/*vector<int> Gstartree::dijkstra_candidate(int s, vector<int> &cands) {//inline
    // init
    auto num_cands = cands.size();
    unordered_set<int> todo(cands.begin(), cands.end());
    unordered_map<int, int> result;
    result[s] = 0;

    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    pq.emplace(0, s);

    int min, minpos, adjnode, weight;
    while (!todo.empty() && !pq.empty()) {
        min = pq.top().first;
        minpos = pq.top().second;
        pq.pop();

        todo.erase(minpos);

        for (int i = 0; i < Nodes[minpos].adjnodes.size(); i++) {

            adjnode = Nodes[minpos].adjnodes[i];
            weight = Nodes[minpos].adjweight[i];

            if (result.find(adjnode) == result.end() || result[adjnode] > min + weight) {
                result[adjnode] = min + weight;
                pq.emplace(min + weight, adjnode);
            }
        }
    }

    // output
    vector<int> output;
    output.reserve(num_cands);
    for (const auto &iter : cands) {
        output.emplace_back(result[iter]);
    }

    // return
    return output;

}*/

inline bool Gstartree::check_shortcut(int ns, int nt) {
    return shortcuts.find(ns) != shortcuts.end() && shortcuts[ns].find(nt) != shortcuts[ns].end();
}

//inline int Gstartree::find_LCA_pos(int src, int dst) {
//    for (int i = 1; i < Nodes[src].deep && i < Nodes[dst].deep; ++i) {
//        if (Nodes[src].gtreepath[i] != Nodes[dst].gtreepath[i])
//            return i - 1;
//    }
//    return 0;
//}

void Gstartree::QueryProcessingTest(int runtimes) {
    switch (indexType) {
        case gtreeIndex:{
            init_ReadIndex();
            CorrectnessCheck(100);
            EfficiencyTest(runtimes);
            break;
        }
        case gstarIndex:{
            init_ReadIndex();
            CorrectnessCheck(100);
            EfficiencyTest(runtimes);
            break;
        }
        case lgtreeIndex:{
            init_LGTreeIndex();
            CorrectnessCheck(100);
            EfficiencyTest(runtimes);
            break;
        }
        case tgtreeIndex:{
            init_TGTreeIndex();
            CorrectnessCheck(100);
            EfficiencyTest(runtimes);
            break;
        }
        default:{
            cout<<"Wrong index type! "<<indexType<<endl; exit(1);
        }
    }
}

//stage-based query for LG-Tree
int Gstartree::Query_LGTree(int ID1, int ID2) {
    int d=INF;
    int Ns = Nodes[ID1].inleaf;
    int Nt = Nodes[ID2].inleaf;
    if(Nodes[ID1].isborder && Nodes[ID2].isborder){//if ID1 and ID2 are boundary vertices
        d = QueryPLL(ID1,ID2);
    }else if(Nodes[ID1].isborder && !Nodes[ID2].isborder ){//if only ID1 is boundary vertex
        d = QueryPartiBoundaryPLL(ID2,ID1);
    }else if(Nodes[ID2].isborder && !Nodes[ID1].isborder){
        d = QueryPartiBoundaryPLL(ID1,ID2);
    }else if(!Nodes[ID1].isborder && !Nodes[ID2].isborder){//if both are not boundary vertex
        if(Ns == Nt){//in the same partition
            cout<<"same-partition query! "<<ID1<<"("<<Ns<<") "<<ID2<<"("<<Nt<<")"<<endl;
//            d = Dijkstra(ID1,ID2,Nodes);
            d=dijkstra_p2p_leafNode(ID1, ID2);
            int dis1,dis2;
            int bid1,bid2;
            int leafSize1=GTree[Ns].leafnodes.size();
            int leafSize2=GTree[Nt].leafnodes.size();
            int posa = Nodes[ID1].inleafpos;//the position in leaf node
            int posb = Nodes[ID2].inleafpos;//the position in leaf node
            int disbb=INF;
            for(int i=0;i<GTree[Ns].borders.size();++i){
                bid1=GTree[Ns].borders[i];
                dis1=GTree[Ns].mind[i*leafSize1+posa];
                for(int j=0;j<GTree[Nt].borders.size();++j){
                    bid2=GTree[Nt].borders[j];
                    dis2=GTree[Nt].mind[j*leafSize2+posb];
                    int posb2=Nodes[bid2].inleafpos;
                    disbb=GTree[Nt].mind[i*leafSize1+posb2];
                    if(d>dis1+dis2+disbb){
                        d=dis1+dis2+disbb;
                    }
                }
            }
        }else{//if in different partitions
            d = QueryPartiPartiPLL(ID1,ID2);
        }
    }
    return d;
}

int Gstartree::Query_LGTreeDebug(int ID1, int ID2) {
    int d=INF;
    int Ns = Nodes[ID1].inleaf;
    int Nt = Nodes[ID2].inleaf;
    if(Nodes[ID1].isborder && Nodes[ID2].isborder){//if ID1 and ID2 are boundary vertices
        cout<<"boundary-boundary"<<endl;
        d = QueryPLLDebug(ID1,ID2);
    }else if(Nodes[ID1].isborder && !Nodes[ID2].isborder ){//if only ID1 is boundary vertex
        cout<<"boundary-parti"<<endl;
        d = QueryPartiBoundaryPLL(ID2,ID1);
    }else if(Nodes[ID2].isborder && !Nodes[ID1].isborder){
        cout<<"parti-boundary"<<endl;
        d = QueryPartiBoundaryPLL(ID1,ID2);
    }else if(!Nodes[ID1].isborder && !Nodes[ID2].isborder){//if both are not boundary vertex
        if(Ns == Nt){//in the same partition
            cout<<"same-partition query! "<<ID1<<"("<<Ns<<") "<<ID2<<"("<<Nt<<")"<<endl;
//            d = Dijkstra(ID1,ID2,Nodes);
            d=dijkstra_p2p_leafNode(ID1, ID2);
            int dis1,dis2;
            int bid1,bid2;
            int leafSize1=GTree[Ns].leafnodes.size();
            int leafSize2=GTree[Nt].leafnodes.size();
            int posa = Nodes[ID1].inleafpos;//the position in leaf node
            int posb = Nodes[ID2].inleafpos;//the position in leaf node
            int disbb=INF;
            for(int i=0;i<GTree[Ns].borders.size();++i){
                bid1=GTree[Ns].borders[i];
                dis1=GTree[Ns].mind[i*leafSize1+posa];
                for(int j=0;j<GTree[Nt].borders.size();++j){
                    bid2=GTree[Nt].borders[j];
                    dis2=GTree[Nt].mind[j*leafSize2+posb];
                    int posb2=Nodes[bid2].inleafpos;
                    disbb=GTree[Nt].mind[i*leafSize1+posb2];
                    if(d>dis1+dis2+disbb){
                        d=dis1+dis2+disbb;
                    }
                }
            }
        }
        else{//if in different partitions
            cout<<"parti-parti"<<endl;
            d = QueryPartiPartiPLLDebug(ID1,ID2);
        }
    }
    return d;
}

int Gstartree::QueryPLL(int ID1, int ID2){
    if(ID1==ID2) return 0;
    if(!Nodes[ID1].isborder || !Nodes[ID2].isborder) {
        cout<<"Wrong for this query pair. "<<ID1<<"("<<Nodes[ID1].isborder<<") "<<ID2<<"("<<Nodes[ID2].isborder<<")"<<endl; exit(1);
    }
    int d=INF;

    int hub, dis1, dis2;
    for(auto it=Label[ID1].begin();it!=Label[ID1].end();it++){
        hub=(*it).first;
        dis1=(*it).second;
        if(Label[ID2].find(hub)!=Label[ID2].end()){
            dis2=Label[ID2][hub];
            if(dis1+dis2<d){
                d=dis1+dis2;
                //cout<<"details "<<hub<<" "<<dis1<<" "<<dis2<<" "<<d<<endl;
            }
        }
    }
    return d;
}

int Gstartree::QueryPLLDebug(int ID1, int ID2){
    if(ID1==ID2) return 0;
    if(!Nodes[ID1].isborder || !Nodes[ID2].isborder) {
        cout<<"Wrong for this query pair. "<<ID1<<"("<<Nodes[ID1].isborder<<") "<<ID2<<"("<<Nodes[ID2].isborder<<")"<<endl; exit(1);
    }
    int d=INF;

    int hub, dis1, dis2;
    int hubf, dis1f, dis2f;
    for(auto it=Label[ID1].begin();it!=Label[ID1].end();it++){
        hub=(*it).first;
        dis1=(*it).second;
        if(Label[ID2].find(hub)!=Label[ID2].end()){
            dis2=Label[ID2][hub];
            if(dis1+dis2<d){
                hubf=hub, dis1f=dis1, dis2f=dis2;
                d=dis1+dis2;
                //cout<<"details "<<hub<<" "<<dis1<<" "<<dis2<<" "<<d<<endl;
            }
        }
    }
    int dis1_Dijk= Dijkstra(ID1,hubf,Nodes), dis2_Dijk= Dijkstra(hubf,ID2,Nodes);
    cout<<ID1<<"("<<Nodes[ID1].isborder<<","<<NodeOrder[ID1]<<") "<<hubf<<"("<<Nodes[hubf].isborder<<","<<NodeOrder[hubf]<<") "<<ID2<<"("<<Nodes[ID2].isborder<<","<<NodeOrder[ID2]<<"): "<<dis1f<<"("<<dis1_Dijk<<") "<<dis2f<<"("<<dis2_Dijk<<") "<<d<<"("<<Dijkstra(ID1,ID2,Nodes)<<")"<<endl;
    return d;
}

int Gstartree::QueryPartiBoundaryPLL(int ID1,int ID2){
    int d=INF;
    assert(Nodes[ID1].isborder==false && Nodes[ID2].isborder==true);
    int PID1=Nodes[ID1].inleaf;
    int PosID1=Nodes[ID1].inleafpos;
    vector<int> B1;
    int bID1,d1;
    int tempdis;
    B1 = GTree[PID1].borders;

    map<int,int> m1  ;
    m1.clear();

    int lnNum1 = GTree[PID1].leafnodes.size();

    //get distances to corresponding boundary vertex
    for(int i=0;i<B1.size();++i){
        bID1 = B1[i];
        d1 = GTree[PID1].mind[i*lnNum1 + PosID1];
        m1.insert({bID1,d1});
    }

    int b1f,b2f,d1f,d2f,dbb,dbbf;
    //get distances between boundary vertex
    for(int i=0;i<B1.size();++i){
        bID1 = B1[i];

        dbb = QueryPLL(bID1,ID2);
        tempdis = m1[bID1] + dbb;
        if(tempdis < d){
            d = tempdis;
            b1f = bID1;
            d1f = m1[bID1];
            dbbf = dbb;
        }

    }
    return d;
}

int Gstartree::QueryPartiPartiPLL(int ID1,int ID2){
    int d=INF;
    assert(Nodes[ID1].isborder==false && Nodes[ID2].isborder==false);
    int PID1=Nodes[ID1].inleaf, PID2=Nodes[ID2].inleaf;
    int PosID1=Nodes[ID1].inleafpos, PosID2=Nodes[ID2].inleafpos;
    vector<int> B1,B2;
    int bID1,bID2,d1,d2;
    int tempdis;
    B1 = GTree[PID1].borders; B2 = GTree[PID2].borders;

    map<int,int> m1,m2;
    m1.clear();
    m2.clear();
    int lnNum1 = GTree[PID1].leafnodes.size();
    int lnNum2 = GTree[PID2].leafnodes.size();
    //get distances to corresponding boundary vertex
    for(int i=0;i<B1.size();++i){
        bID1 = B1[i];
        d1 = GTree[PID1].mind[i*lnNum1 + PosID1];
        m1.insert({bID1,d1});
    }
    for(int i=0;i<B2.size();++i){
        bID2 = B2[i];
        d2 = GTree[PID2].mind[i*lnNum2 + PosID2];
        m2.insert({bID2,d2});
    }
    int b1f,b2f,d1f,d2f,dbb,dbbf;
    //get distances between boundary vertex
    for(int i=0;i<B1.size();++i){
        bID1 = B1[i];
        for(int j=0;j<B2.size();++j){
            bID2 = B2[j];
            dbb = QueryPLL(bID1,bID2);
            tempdis = m1[bID1] + dbb + m2[bID2];
            if(tempdis < d){
                d = tempdis;
                b1f = bID1; b2f = bID2;
                d1f = m1[bID1]; d2f =  m2[bID2];
                dbbf = dbb;
            }
        }
    }

    return d;
}

int Gstartree::QueryPartiPartiPLLDebug(int ID1,int ID2){
    int d=INF;
    assert(Nodes[ID1].isborder==false && Nodes[ID2].isborder==false);
    int PID1=Nodes[ID1].inleaf, PID2=Nodes[ID2].inleaf;
    int PosID1=Nodes[ID1].inleafpos, PosID2=Nodes[ID2].inleafpos;
    vector<int> B1,B2;
    int bID1,bID2,d1,d2;
    int tempdis;
    B1 = GTree[PID1].borders; B2 = GTree[PID2].borders;

    map<int,int> m1,m2;
    m1.clear();
    m2.clear();
    int lnNum1 = GTree[PID1].leafnodes.size();
    int lnNum2 = GTree[PID2].leafnodes.size();
    //get distances to corresponding boundary vertex
    for(int i=0;i<B1.size();++i){
        bID1 = B1[i];
        d1 = GTree[PID1].mind[i*lnNum1 + PosID1];
        m1.insert({bID1,d1});
    }
    for(int i=0;i<B2.size();++i){
        bID2 = B2[i];
        d2 = GTree[PID2].mind[i*lnNum2 + PosID2];
        m2.insert({bID2,d2});
    }
    int b1f,b2f,d1f,d2f,dbb,dbbf;
    //get distances between boundary vertex
    for(int i=0;i<B1.size();++i){
        bID1 = B1[i];
        for(int j=0;j<B2.size();++j){
            bID2 = B2[j];
            dbb = QueryPLL(bID1,bID2);
            tempdis = m1[bID1] + dbb + m2[bID2];
            if(tempdis < d){
                d = tempdis;
                b1f = bID1; b2f = bID2;
                d1f = m1[bID1]; d2f =  m2[bID2];
                dbbf = dbb;
            }
        }
    }

    int dis1_Dijk= Dijkstra(ID1,b1f,Nodes), dis2_Dijk= Dijkstra(b1f,b2f,Nodes),  dis3_Dijk= Dijkstra(b2f,ID2,Nodes);
    cout<<ID1<<"("<<Nodes[ID1].isborder<<","<<Nodes[ID1].inleaf<<") "<<b1f<<"("<<Nodes[b1f].isborder<<","<<Nodes[b1f].inleaf<<") "<<b2f<<"("<<Nodes[b2f].isborder<<","<<Nodes[b2f].inleaf<<") "<<ID2<<"("<<Nodes[ID2].isborder<<","<<Nodes[ID2].inleaf<<"): "<<d1f<<"("<<dis1_Dijk<<") "<<dbbf<<"("<<dis2_Dijk<<") "<<d2f<<"("<<dis3_Dijk<<") "<<d<<"("<<Dijkstra(ID1,ID2,Nodes)<<")"<<endl;
    return d;
}

//stage-based query for No-Boundary
int Gstartree::Query_TGTree(int ID1, int ID2) {
    int d=INF;
    int Ns = Nodes[ID1].inleaf;
    int Nt = Nodes[ID2].inleaf;
    if(Nodes[ID1].isborder && Nodes[ID2].isborder){//if ID1 and ID2 are boundary vertices
        d = QueryH2H(ID1,ID2);
    }else if(Nodes[ID1].isborder && !Nodes[ID2].isborder ){//if only ID1 is boundary vertex
        d = QueryPartiBoundary(ID2,ID1);
    }else if(Nodes[ID2].isborder && !Nodes[ID1].isborder){
        d = QueryPartiBoundary(ID1,ID2);
    }else if(!Nodes[ID1].isborder && !Nodes[ID2].isborder){//if both are not boundary vertex
        if(Ns == Nt){//in the same partition
            cout<<"same-partition query! "<<ID1<<"("<<Ns<<") "<<ID2<<"("<<Nt<<")"<<endl;
//            d = Dijkstra(ID1,ID2,Nodes);
            d=dijkstra_p2p_leafNode(ID1, ID2);
            int dis1,dis2;
            int bid1,bid2;
            int leafSize1=GTree[Ns].leafnodes.size();
            int leafSize2=GTree[Nt].leafnodes.size();
            int posa = Nodes[ID1].inleafpos;//the position in leaf node
            int posb = Nodes[ID2].inleafpos;//the position in leaf node
            int disbb=INF;
            for(int i=0;i<GTree[Ns].borders.size();++i){
                bid1=GTree[Ns].borders[i];
                dis1=GTree[Ns].mind[i*leafSize1+posa];
                for(int j=0;j<GTree[Nt].borders.size();++j){
                    bid2=GTree[Nt].borders[j];
                    dis2=GTree[Nt].mind[j*leafSize2+posb];
                    int posb2=Nodes[bid2].inleafpos;
                    disbb=GTree[Nt].mind[i*leafSize1+posb2];
                    if(d>dis1+dis2+disbb){
                        d=dis1+dis2+disbb;
                    }
                }
            }
        }else{//if in different partitions
            d = QueryPartiParti(ID1,ID2);
        }
    }
    return d;
}

int Gstartree::Query_TGTreeDebug(int ID1, int ID2) {
    int d=INF;
    int Ns = Nodes[ID1].inleaf;
    int Nt = Nodes[ID2].inleaf;
    if(Nodes[ID1].isborder && Nodes[ID2].isborder){//if ID1 and ID2 are boundary vertices
        cout<<"boundary-boundary."<<endl;
        d = QueryH2HDebug(ID1,ID2);
    }else if(Nodes[ID1].isborder && !Nodes[ID2].isborder ){//if only ID1 is boundary vertex
        cout<<"boundary-parti."<<endl;
        d = QueryPartiBoundaryDebug(ID2,ID1);
    }else if(Nodes[ID2].isborder && !Nodes[ID1].isborder){
        cout<<"parti-boundary."<<endl;
        d = QueryPartiBoundaryDebug(ID1,ID2);
    }else if(!Nodes[ID1].isborder && !Nodes[ID2].isborder){//if both are not boundary vertex
        if(Ns == Nt){//in the same partition
            cout<<"same-partition query! "<<ID1<<"("<<Ns<<") "<<ID2<<"("<<Nt<<")"<<endl;
//            d = Dijkstra(ID1,ID2,Nodes);
            d=dijkstra_p2p_leafNode(ID1, ID2);
            int dis1,dis2;
            int bid1,bid2;
            int borderSize=GTree[Ns].borders.size();
            int posa = Nodes[ID1].inleafpos;//the position in leaf node
            int posb = Nodes[ID2].inleafpos;//the position in leaf node
            for(int i=0;i<GTree[Ns].borders.size();++i){
                bid1=GTree[Ns].borders[i];
                dis1=GTree[Ns].mind[i*borderSize+posa];
                for(int j=0;j<GTree[Ns].borders.size();++j){
                    bid2=GTree[Ns].borders[j];
                    dis2=GTree[Ns].mind[j*borderSize+posb];
                    if(d>dis1+dis2){
                        d=dis1+dis2;
                    }
                }
            }
        }
        else{//if in different partitions
            cout<<"parti-parti."<<endl;
            d = QueryPartiParti(ID1,ID2);
        }
    }
    return d;
}

int Gstartree::QueryH2H(int ID1,int ID2){
    if(ID1==ID2) return 0;
    if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1)
        return INF;
    int r1=rank[ID1], r2=rank[ID2];
    int LCA=LCAQuery(r1,r2);

    if(LCA==r1)
        return Tree[r2].dis[Tree[r1].pos.back()];
    else if(LCA==r2)
        return Tree[r1].dis[Tree[r2].pos.back()];
    else{
        int tmp=INF;
        for(int i=0;i<Tree[LCA].pos.size();i++){
            if(tmp>Tree[r1].dis[Tree[LCA].pos[i]]+Tree[r2].dis[Tree[LCA].pos[i]])
                tmp=Tree[r1].dis[Tree[LCA].pos[i]]+Tree[r2].dis[Tree[LCA].pos[i]];
        }
        return tmp;
    }
}

int Gstartree::QueryH2HDebug(int ID1,int ID2){
    if(ID1==ID2) return 0;
    if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1)
        return INF;
    int r1=rank[ID1], r2=rank[ID2];
    int LCA=LCAQuery(r1,r2);
    int hub,dis1,dis2,dis1F,dis2F;
    int dis=INF;

    if(LCA==r1){
        cout<<"LCA=r1 "<<Tree[LCA].uniqueVertex<<endl;
        return Tree[r2].dis[Tree[r1].pos.back()];
    }
    else if(LCA==r2){
        cout<<"LCA=r2 "<<Tree[LCA].uniqueVertex<<endl;
        return Tree[r1].dis[Tree[r2].pos.back()];
    }
    else{
        cout<<"LCA= "<<Tree[LCA].uniqueVertex<<endl;
        for(int i=0;i<Tree[LCA].pos.size();i++){
            dis1=Tree[r1].dis[Tree[LCA].pos[i]]; dis2=Tree[r2].dis[Tree[LCA].pos[i]];
            if(dis>dis1+dis2){
                dis1F=dis1, dis2F=dis2; hub=Tree[LCA].vert[i].first;
                dis=dis1+dis2;
            }
        }

    }

    int dis1_Dijk= Dijkstra(ID1,hub,Nodes), dis2_Dijk= Dijkstra(hub,ID2,Nodes);
    cout<<ID1<<"("<<Nodes[ID1].isborder<<","<<Nodes[ID1].inleaf<<") "<<hub<<"("<<Nodes[hub].isborder<<","<<Nodes[hub].inleaf<<") "<<ID2<<"("<<Nodes[ID2].isborder<<","<<Nodes[ID2].inleaf<<"): "<<dis1F<<"("<<dis1_Dijk<<") "<<dis2F<<"("<<dis2_Dijk<<") "<<dis<<"("<<Dijkstra(ID1,ID2,Nodes)<<")"<<endl;
    return dis;
}

int Gstartree::QueryPartiBoundary(int ID1,int ID2){
    int d=INF;
    assert(Nodes[ID1].isborder==false && Nodes[ID2].isborder==true);
    int PID1=Nodes[ID1].inleaf;
    int PosID1=Nodes[ID1].inleafpos;
    vector<int> B1;
    int bID1,d1;
    int tempdis;
    B1 = GTree[PID1].borders;

    map<int,int> m1  ;
    m1.clear();

    int lnNum1 = GTree[PID1].leafnodes.size();

    //get distances to corresponding boundary vertex
    for(int i=0;i<B1.size();++i){
        bID1 = B1[i];
        d1 = GTree[PID1].mind[i*lnNum1 + PosID1];
        m1.insert({bID1,d1});
    }

    int b1f,b2f,d1f,d2f,dbb,dbbf;
    //get distances between boundary vertex
    for(int i=0;i<B1.size();++i){
        bID1 = B1[i];

        dbb = QueryH2H(bID1,ID2);
        tempdis = m1[bID1] + dbb;
        if(tempdis < d){
            d = tempdis;
            b1f = bID1;
            d1f = m1[bID1];
            dbbf = dbb;
        }

    }
    return d;
}

int Gstartree::QueryPartiBoundaryDebug(int ID1,int ID2){
    int d=INF;
    assert(Nodes[ID1].isborder==false && Nodes[ID2].isborder==true);
    int PID1=Nodes[ID1].inleaf;
    int PosID1=Nodes[ID1].inleafpos;
    vector<int> B1;
    int bID1,d1;
    int tempdis;
    B1 = GTree[PID1].borders;

    map<int,int> m1  ;
    m1.clear();

    int lnNum1 = GTree[PID1].leafnodes.size();

    //get distances to corresponding boundary vertex
    for(int i=0;i<B1.size();++i){
        bID1 = B1[i];
        d1 = GTree[PID1].mind[i*lnNum1 + PosID1];
        m1.insert({bID1,d1});
    }

    int b1f,b2f,d1f,d2f,dbb,dbbf;
    //get distances between boundary vertex
    for(int i=0;i<B1.size();++i){
        bID1 = B1[i];

        dbb = QueryH2H(bID1,ID2);
        tempdis = m1[bID1] + dbb;
        if(tempdis < d){
            d = tempdis;
            b1f = bID1;
            d1f = m1[bID1];
            dbbf = dbb;
        }

    }
    int dis1_Dijk= Dijkstra(ID1,b1f,Nodes), dis2_Dijk= Dijkstra(b1f,ID2,Nodes);
    cout<<ID1<<"("<<Nodes[ID1].isborder<<","<<Nodes[ID1].inleaf<<") "<<b1f<<"("<<Nodes[b1f].isborder<<","<<Nodes[b1f].inleaf<<") "<<ID2<<"("<<Nodes[ID2].isborder<<","<<Nodes[ID2].inleaf<<"): "<<d1f<<"("<<dis1_Dijk<<") "<<dbbf<<"("<<dis2_Dijk<<") "<<d<<"("<<Dijkstra(ID1,ID2,Nodes)<<")"<<endl;

    return d;
}

int Gstartree::QueryPartiParti(int ID1,int ID2){
    int d=INF;
    assert(Nodes[ID1].isborder==false && Nodes[ID2].isborder==false);
    int PID1=Nodes[ID1].inleaf, PID2=Nodes[ID2].inleaf;
    int PosID1=Nodes[ID1].inleafpos, PosID2=Nodes[ID2].inleafpos;
    vector<int> B1,B2;
    int bID1,bID2,d1,d2;
    int tempdis;
    B1 = GTree[PID1].borders; B2 = GTree[PID2].borders;

    map<int,int> m1,m2;
    m1.clear();
    m2.clear();
    int lnNum1 = GTree[PID1].leafnodes.size();
    int lnNum2 = GTree[PID2].leafnodes.size();
    //get distances to corresponding boundary vertex
    for(int i=0;i<B1.size();++i){
        bID1 = B1[i];
        d1 = GTree[PID1].mind[i*lnNum1 + PosID1];
        m1.insert({bID1,d1});
    }
    for(int i=0;i<B2.size();++i){
        bID2 = B2[i];
        d2 = GTree[PID2].mind[i*lnNum2 + PosID2];
        m2.insert({bID2,d2});
    }
    int b1f,b2f,d1f,d2f,dbb,dbbf;
    //get distances between boundary vertex
    for(int i=0;i<B1.size();++i){
        bID1 = B1[i];
        for(int j=0;j<B2.size();++j){
            bID2 = B2[j];
            dbb = QueryH2H(bID1,bID2);
            tempdis = m1[bID1] + dbb + m2[bID2];
            if(tempdis < d){
                d = tempdis;
                b1f = bID1; b2f = bID2;
                d1f = m1[bID1]; d2f =  m2[bID2];
                dbbf = dbb;
            }
        }
    }
    return d;
}

void Gstartree::CorrectnessCheck(int times){
    srand (time(NULL));
    int s, t;
    int dis1, dis2;
    Timer tt1,tt2;
    double time_Dijk=0,time_Q=0;
    double sameNodeT=0, diffNodeT=0;
    int times1=0, times2=0;

//    times=1;
    cout<<"Correctness check... Query number: "<<times<<endl;
    for(int i=0;i<times;i++){//times
        s=rand()%node_num;
        t=rand()%node_num;
//        s=235330, t=233696;
//        s=196894; t=145058;//b-b
//        s=196894; t=73422;//b-b
//        s=225621; t=145793;//b-b
//        int index_i=rand()%leafNodePairs.size();
//        s=GTree[leafNodePairs[index_i].first].leafnodes[rand()%GTree[leafNodePairs[index_i].first].leafnodes.size()], t=GTree[leafNodePairs[index_i].second].leafnodes[rand()%GTree[leafNodePairs[index_i].second].leafnodes.size()];

        if(indexType==gtreeIndex || indexType==gstarIndex){
            init_query(GTree);
            tt2.start();
//        dis_GTree=Distance_query(s,t);
            dis1=dist_query(s,t);
            tt2.stop();
        }
        else if(indexType==lgtreeIndex){
            tt2.start();
            dis1=Query_LGTree(s,t);
            tt2.stop();
        }
        else if(indexType==tgtreeIndex){
            tt2.start();
            dis1=Query_TGTree(s,t);
            tt2.stop();
        }

        time_Q += tt2.GetRuntime();
        if(Nodes[s].inleaf==Nodes[t].inleaf){
            sameNodeT+=tt2.GetRuntime();
            times1++;
        }else{
            diffNodeT+=tt2.GetRuntime();
            times2++;
        }

        tt1.start();
        dis2=Dijkstra(s,t,Nodes);
        tt1.stop();
        time_Dijk += tt1.GetRuntime();

        if(dis1!=dis2){
            cout<<"InCorrect!! "<<i<<": "<<s<<"("<<Nodes[s].isborder<<","<<Nodes[s].inleaf<<") "<<t<<"("<<Nodes[t].isborder<<","<<Nodes[t].inleaf<<") "<<dis1<<" "<<dis2<<endl;
            if(indexType==tgtreeIndex){
                Query_TGTreeDebug(s,t);
            }else if(indexType==lgtreeIndex){
                Query_LGTreeDebug(s,t);
            }

            exit(1);
        }

        //cout<<"PID "<<VtoParID[s]<<" "<<VtoParID[t]<<endl;
    }
//    cout << "Average run time of Dijkstra: "<<time_Dijk*1000/times<<" ms."<<endl;
    cout << "Average time: "<<time_Q*1000/times<<" ms."<<endl;
    cout<<"Same leaf node query time: "<<1000*sameNodeT/times1<<" ms ("<<times1<<") ; different leaf node query time: "<<1000*diffNodeT/times2<<" ms ("<<times2<<")"<<endl;
    //cout<<"Correctness finish!"<<endl;
}

inline int Gstartree::dist_query(int src, int dst) {
    if(src==dst){
        return 0;
    }
    int dis=INF;
    int Ns = Nodes[src].inleaf;
    int Nt = Nodes[dst].inleaf;

    // Case D1: vertices src and dst are in the same leaf node
    if (Ns == Nt) {
//        cout<<src << " " << dst << " : s and t are in the same leaf node."<<endl;
//        return dijkstra_p2p(src, dst);
        dis=dijkstra_p2p_leafNode(src, dst);
        int dis1,dis2;
        int bid1,bid2;
        int leafSize1=GTree[Ns].leafnodes.size();
        int leafSize2=GTree[Nt].leafnodes.size();
        int posa = Nodes[src].inleafpos;//the position in leaf node
        int posb = Nodes[dst].inleafpos;//the position in leaf node
        int disbb=INF;
        for(int i=0;i<GTree[Ns].borders.size();++i){
            bid1=GTree[Ns].borders[i];
            dis1=GTree[Ns].mind[i*leafSize1+posa];
            for(int j=0;j<GTree[Nt].borders.size();++j){
                bid2=GTree[Nt].borders[j];
                dis2=GTree[Nt].mind[j*leafSize2+posb];
                int posb2=Nodes[bid2].inleafpos;
                disbb=GTree[Nt].mind[i*leafSize1+posb2];
                if(dis>dis1+dis2+disbb){
                    dis=dis1+dis2+disbb;
                }
            }
        }
    }
    else{
        int posa = Nodes[src].inleafpos;//the position in leaf node
        int posb = Nodes[dst].inleafpos;//the position in leaf node
        auto num_border_ns = GTree[Ns].borders.size();
        auto num_leafnode_ns = GTree[Ns].leafnodes.size();
        auto num_border_nt = GTree[Nt].borders.size();
        auto num_leafnode_nt = GTree[Nt].leafnodes.size();

        // Case D2: there is a shortcut between Ns and Nt
        if(indexType==gstarIndex){
            if (check_shortcut(Ns, Nt)) {
//                cout<<src << " " << dst << " : there is a shortcut between Ns and Nt."<<endl;
                vector<int> sc = Ns < Nt ? shortcuts[Ns][Nt] : shortcuts[Nt][Ns];//chose the smaller one as the first dimension
                int temp_dist;
                for (int i = 0; i < num_border_ns; ++i) {
                    GTree[Ns].cache[i] = GTree[Ns].mind[i * num_leafnode_ns + posa];
//                    temp_dist=Dijkstra(src,GTree[Ns].borders[i],Nodes);
//                    if(GTree[Ns].cache[i] != temp_dist){
//                        cout << src << " "<< GTree[Ns].borders[i] << " " << GTree[Ns].cache[i] << " " << temp_dist << endl;
//                    }
                }

                for (int j = 0; j < num_border_nt; ++j) {
                    GTree[Nt].cache[j] = GTree[Nt].mind[j * num_leafnode_nt + posb];
//                    temp_dist = Dijkstra(dst,GTree[Nt].borders[j],Nodes);
//                    if(GTree[Nt].cache[j] != temp_dist){
//                        cout << dst << " "<< GTree[Nt].borders[j] << " " << GTree[Nt].cache[j] << " " << temp_dist << endl;
//                    }
                }

                int dist;

                if (Ns < Nt) {
                    for (int i = 0; i < num_border_ns; ++i) {
                        for (int j = 0; j < num_border_nt; ++j) {
                            assert(i * num_border_nt + j < sc.size());
                            dist = GTree[Ns].cache[i] + GTree[Nt].cache[j] + sc[i * num_border_nt + j];
//                            temp_dist = Dijkstra(GTree[Ns].borders[i],GTree[Nt].borders[j],Nodes);
//                            if(sc[i * num_border_nt + j] != temp_dist){
//                                cout << GTree[Ns].borders[i] << " "<< GTree[Nt].borders[j] << " " << sc[i * num_border_nt + j] << " " << temp_dist << endl;
//                                exit(1);
//                            }
                            if (dist < dis) {
                                dis = dist;
                            }
                        }
                    }
                } else {
                    for (int i = 0; i < num_border_ns; ++i) {
                        for (int j = 0; j < num_border_nt; ++j) {
                            assert(j * num_border_ns + i < sc.size());
                            dist = GTree[Ns].cache[i] + GTree[Nt].cache[j] +  sc[j * num_border_ns + i];
//                            temp_dist = Dijkstra(GTree[Ns].borders[i],GTree[Nt].borders[j],Nodes);
//                            if(sc[j * num_border_ns + i] != temp_dist){
//                                cout << GTree[Ns].borders[i] << " "<< GTree[Nt].borders[j] << " " << sc[j * num_border_ns + i] << " " << temp_dist << endl;
//                            }
                            if (dist < dis) {
                                dis = dist;
                            }
                        }
                    }
                }
                return dis;
            }
        }

        // Case D3: there is no shortcut between Ns and Nt

        // Find LCA index in gtreepath
        int LCA_pos = find_LCA_pos(src, dst);

        // Step out of leaf node Ns
        for (int i = 0; i < num_border_ns; ++i) {
//        cout << i*num_leafnode_ns+posa << " "<<GTree[Ns].mind.size() << endl;
            GTree[Ns].cache[i] = GTree[Ns].mind[i * num_leafnode_ns + posa];//store the distance from all borders to s
//        cout << src << " "<< GTree[Ns].borders[i] << " " << GTree[Ns].cache[i] << " " << Dijkstra(src,GTree[Ns].borders[i]) << endl;
        }

        // Init some variables
        const auto &up_path = Nodes[src].gtreepath;
        const auto &down_path = Nodes[dst].gtreepath;
        int cn, tn, min, dist, posx, posy;
        unsigned long union_border_size;

        // Step out of nodes until meeting LCA
        // The cache of each node 'tn' stores the distance from vertex src to node tn's child then to tn
        TIME_TICK_START
        for (auto i = up_path.size() - 2; i >= LCA_pos + 1; --i) {
            tn = up_path[i];
            cn = up_path[i + 1];  // child node
            union_border_size = GTree[tn].union_borders.size();
            for (int j = 0; j < GTree[tn].borders.size(); j++) {//
                min = INT_MAX;
                posx = GTree[tn].current_pos[j];//get the position id of tn's borders in this node
                for (int k = 0; k < GTree[cn].borders.size(); k++) {
                    posy = GTree[cn].up_pos[k];//get the position id of cn's borders in the parent node
                    dist = GTree[cn].cache[k] + GTree[tn].mind[posx * union_border_size + posy];
                    if (dist < min) {
                        min = dist;
                    }
                }
                GTree[tn].cache[j] = min;//update the distance from the borders to its children's borders
//            cout << src << " "<< GTree[tn].borders[j] << " " << GTree[tn].cache[j] << " " << Dijkstra(src,GTree[tn].borders[j]) << endl;
            }
        }


        // Step across LCA (from one branch to another)
        // The cache of Nt's top ancestor node 'nt_top' stores the distance
        // from vertex src to Ns's top ancestor 'ns_top' node then to 'nt_top'
//    cout<<"LCA"<<endl;
        int ns_top = up_path[LCA_pos + 1];
        int nt_top = down_path[LCA_pos + 1];
        int lca_node = up_path[LCA_pos];
        int ns_top_up_pos, nt_top_up_pos;
        union_border_size = GTree[lca_node].union_borders.size();
        for (int i = 0; i < GTree[nt_top].borders.size(); i++) {
            min = INT_MAX;
            nt_top_up_pos = GTree[nt_top].up_pos[i];//the position id in parent node (LCA)
            for (int j = 0; j < GTree[ns_top].borders.size(); j++) {
                ns_top_up_pos = GTree[ns_top].up_pos[j];//the position id in parent node (LCA)
//            dist = GTree[ns_top].cache[i] + GTree[lca_node].mind[nt_top_up_pos * union_border_size + ns_top_up_pos];
                dist = GTree[ns_top].cache[j] + GTree[lca_node].mind[nt_top_up_pos * union_border_size + ns_top_up_pos];
                if (dist < min) {
                    min = dist;
                }
            }
            GTree[nt_top].cache[i] = min;
//        cout << src << " "<< GTree[nt_top].borders[i] << " " << GTree[nt_top].cache[i] << " " << Dijkstra(src,GTree[nt_top].borders[i]) << endl;
        }


        // Step into nodes until meeting Nt
        // The cache of each node 'tn' stores the distance from vertex src to node tn's parent then to tn
        for (auto i = LCA_pos + 2; i < down_path.size(); ++i) {
            tn = down_path[i];
            cn = down_path[i - 1];   // parent node
            union_border_size = GTree[cn].union_borders.size();
            for (int j = 0; j < GTree[tn].borders.size(); j++) {
                min = INT_MAX;
                posx = GTree[tn].up_pos[j];
                for (int k = 0; k < GTree[cn].borders.size(); k++) {
                    posy = GTree[cn].current_pos[k];
                    dist = GTree[cn].cache[k] + GTree[cn].mind[posy * union_border_size + posx];
                    if (dist < min) {
                        min = dist;
                    }
                }
                // update
                GTree[tn].cache[j] = min;
//            cout << src << " "<< GTree[tn].borders[j] << " " << GTree[tn].cache[j] << " " << Dijkstra(src,GTree[tn].borders[j]) << endl;
            }
        }

        // Step into the leaf node Nt
        min = INT_MAX;
        for (int i = 0; i < num_border_nt; ++i) {
            dist = GTree[Nt].cache[i] + GTree[Nt].mind[i * num_leafnode_nt + posb];
            if (dist < min) {
                min = dist;
            }
        }
        TIME_TICK_END
        dis=min;
    }


    return dis;
}

void Gstartree::EfficiencyTest(int run_times) {
//    for(int i=0;i<GTree.size();++i){
//        if(GTree[i].isleaf){
//            cout<<i<<"("<<GTree[i].leafnodes.size()<<"):";
//            for(int j=0;j<GTree[i].leafnodes.size();++j){
//                cout<<" "<<GTree[i].leafnodes[j];
//            }
//            cout<<endl;
//            break;
//        }
//    }

    bool ifDebug=false;
//    ifDebug=true;

    int ID1, ID2, num;
    vector<pair<int, int>> ODpair;

    long double query_time;

    string filename;
    if(dataset == "cal"){
        filename = DataPath + FILE_QUERY;
    }else {
        filename = DataPath +  "/" + dataset + "/" + FILE_QUERY;
    }
    if(percentScale!=0){
        filename = DataPath +  "/" + dataset + "/" + dataset +"_"+ to_string(percentScale)+".query";
    }
    cout<<"Efficiency test. Query file: "<<filename<<endl;
    //Open file
    ifstream inFile(filename, ios::in);
    if (!inFile) {
        cout << "Failed to open file " << filename << endl;
        exit(1);
    }
    inFile >> num;
    for (int i = 0; i < run_times; ++i) {
        inFile >> ID1 >> ID2;
        ODpair.emplace_back(make_pair(ID1, ID2));
    }
    inFile.close();

    long double total_time = 0.0;
    int count = 0;
    int dis1 = 0;
    int dis2 = 0;

    cout << "Begin query processing... Query number: "<<run_times;
    if(ifDebug){
        cout<<" ; with correctness check.";
    }
    cout<<endl;

    Timer tt;
    double runT=0;
    double sameNodeT=0, diffNodeT=0;
    int times1=0, times2=0;
    for (int i = 0; i < run_times; ++i) {//0
        ID1 = ODpair[i].first;
        ID2 = ODpair[i].second;
        if(indexType==gtreeIndex || indexType==gstarIndex){
            init_query(GTree);
            tt.start();
            dis1 = dist_query(ID1, ID2);
            tt.stop();
        }else if(indexType==lgtreeIndex){
            tt.start();
            dis1 = Query_LGTree(ID1, ID2);
            tt.stop();
        }else if(indexType==tgtreeIndex){
            tt.start();
            dis1 = Query_TGTree(ID1, ID2);
            tt.stop();
        }

        if(Nodes[ID1].inleaf==Nodes[ID2].inleaf){
            sameNodeT+=tt.GetRuntime();
            times1++;
        }else{
            diffNodeT+=tt.GetRuntime();
            times2++;
        }
        runT+=tt.GetRuntime();
        if(ifDebug){
            dis2 = Dijkstra(ID1,ID2,Nodes);
            if(dis1 != dis2){
                cout << "Incorrect! " <<i<<": "<< ID1 << " " << ID2 << " : " << dis1 << " " << dis2 << endl;
            }
        }
        ++count;
    }

//    cout<<"Done."<<endl;

    cout << "The average time for querying: " << runT*1000/run_times << " ms." << endl;
    cout<<"Same leaf node query time: "<<1000*sameNodeT/times1<<" ms ("<<times1<<") ; different leaf node query time: "<<1000*diffNodeT/times2<<" ms ("<<times2<<")"<<endl;

//    cout << (total_time / count) << endl;
}


/// For N-TS-HP
//function of vertex ordering
/*void Gstartree::getOverlayOrder(){
    NodeOrder.assign(node_num,-1);
    vNodeOrder.assign(node_num,-1);
    int count=0;

    int curID=0;
    int curLevel=0;
    int tempLevel;
    //get boundary order by BFS
//    getBOrder(curID,HighToLow);
    queue<pair<int,int>> queue1;//vertex ID, tree level
    queue1.push(make_pair(curID,0));
    boundaryLevel.clear();
    pair<int,int> topPair;
    vector<int> currentB;
    //Get order in top-down manner
    while(!queue1.empty()){
        topPair = queue1.front();
        curID = topPair.first; tempLevel = topPair.second;
        queue1.pop();
        if(curLevel != tempLevel){
            boundaryLevel.push_back(currentB);
            currentB.clear();
            curLevel = tempLevel;
        }
        if(!GTree[curID].children.empty()){

            for(int i=0;i<GTree[curID].children.size();++i){
                int pID = GTree[curID].children[i];
                queue1.push(make_pair(pID,tempLevel+1));

                for(int j=0;j<GTree[pID].borders.size();++j){
                    int ID = GTree[pID].borders[j];
                    if(!Nodes[ID].isborder){
                        cout<<"Wrong border! "<<ID<<endl; exit(1);
                    }
                    if(Nodes[ID].isborder && NodeOrder[ID] == -1){//if ID is boundary vertex and it has never been assigned order
                        count++;
                        NodeOrder[ID] = node_num - count;
                        assert(NodeOrder[ID]>=0 && NodeOrder[ID]<node_num);
                        vNodeOrder[NodeOrder[ID]] = ID;
                        currentB.push_back(ID);
                    }
                }
            }
        }
    }

    cout<<"Number of layers: "<<boundaryLevel.size()<<endl;
//    cout<<"Node number of overlay graph: "<<OverlayGraph.size()<<endl;

    if(count != overlayNodeNum){
        cout<<"count inconsistent! "<<count<<" "<<overlayNodeNum<<endl; exit(1);
    }

//    for(int i=0;i<node_num;++i){
//        if(NodeOrder[i] == -1){
//            HighToLow.push_back(i);
//            NodeOrder[i] = node_num - HighToLow.size();
//            vNodeOrder[NodeOrder[i]] = i;
//        }
//    }
//    assert(HighToLow.size() == node_num);

}*/

void Gstartree::getOverlayLayers(){
    Timer tt;
    tt.start();
    vector<bool> ifDone(node_num,false);
    int count=0;
    int curID=0;
    int curLevel=0;
    int tempLevel;
    //get boundary order by BFS
//    getBOrder(curID,HighToLow);
    queue<pair<int,int>> queue1;//vertex ID, tree level
    queue1.push(make_pair(curID,0));
    boundaryLevel.clear();
    pair<int,int> topPair;
    vector<int> currentB;
    //Get order in top-down manner
    while(!queue1.empty()){
        topPair = queue1.front();
        curID = topPair.first; tempLevel = topPair.second;
        queue1.pop();
        if(curLevel != tempLevel){
            boundaryLevel.push_back(currentB);
            currentB.clear();
            curLevel = tempLevel;
        }
        if(!GTree[curID].children.empty()){

            for(int i=0;i<GTree[curID].children.size();++i){
                int pID = GTree[curID].children[i];
                queue1.push(make_pair(pID,tempLevel+1));

                for(int j=0;j<GTree[pID].borders.size();++j){
                    int ID = GTree[pID].borders[j];
                    if(!Nodes[ID].isborder){
                        cout<<"Wrong border! "<<ID<<endl; exit(1);
                    }
                    if(Nodes[ID].isborder && !ifDone[ID]){//if ID is boundary vertex and it has never been assigned order
                        count++;
                        currentB.push_back(ID);
                        ifDone[ID]=true;
                    }
                }
            }
        }
    }
    tt.stop();
    cout<<"Number of layers: "<<boundaryLevel.size()<<" ; time: "<<tt.GetRuntime()<<endl;
//    cout<<"Node number of overlay graph: "<<OverlayGraph.size()<<endl;

    if(count != overlayNodeNum){
        cout<<"count inconsistent! "<<count<<" "<<overlayNodeNum<<endl; exit(1);
    }


//    for(int i=0;i<node_num;++i){
//        if(NodeOrder[i] == -1){
//            HighToLow.push_back(i);
//            NodeOrder[i] = node_num - HighToLow.size();
//            vNodeOrder[NodeOrder[i]] = i;
//        }
//    }
//    assert(HighToLow.size() == node_num);
}




void Gstartree::getOriginalOverlayGraph(){
    int neiID,weight;
    OverlayGraph.assign(node_num,unordered_map<int,int>());
    overlayNodeNum=0;
    for(int i=0;i<Nodes.size();i++){
        if(Nodes[i].isborder){
            overlayNodeNum++;
            for(int j=0;j<Nodes[i].adjnodes.size();j++){
                neiID=Nodes[i].adjnodes[j], weight=Nodes[i].adjweight[j];
                if(Nodes[neiID].isborder){
                    OverlayGraph[i].insert({neiID,weight});
                }
            }
        }
    }
    cout<<"Overlay graph size: "<<overlayNodeNum<<endl;
}
//function of obtaining the partition graphs
/*void Gstartree::getPartiGraph(){
    map<int,unordered_set<int>> partiVSet;
    vector<int> leafNodes;
    for(int i=0;i<GTree.size();++i){
        if(GTree[i].isleaf){
            leafNodes.push_back(i);

            unordered_set<int> tempSet; tempSet.clear();
            tempSet.insert(GTree[i].leafnodes.begin(), GTree[i].leafnodes.end());//get the partition vertex set
            partiVSet.insert({i,tempSet});
        }
    }
    int pID,ID;
    int nid,wei;
    for(int i=0;i<leafNodes.size();++i){
        pID = leafNodes[i];
        PartiGraph.insert({pID,unordered_map<int,vector<pair<int,int>>>()});
        for(int j=0;j<GTree[pID].leafnodes.size();++j){
            ID = GTree[pID].leafnodes[j];
            vector<pair<int,int>> adjs;
            for(int k=0;k<Nodes[ID].adjnodes.size();++k){
                nid = Nodes[ID].adjnodes[k];
                if(partiVSet[pID].find(nid) != partiVSet[pID].end()){//if found
                    wei = Nodes[ID].adjweight[k];
                    adjs.emplace_back(nid,wei);
                }
            }
            PartiGraph[pID].insert({ID,adjs});
        }
    }
}*/
//calculate the distance matrix of leaf nodes by the no-boundary strategy
void Gstartree::hierarchy_shortest_path_calculation_NoBoundary(bool ifParallel){
    //get partition graph
    Timer tt;
//    tt.start();
//    getPartiGraph();
//    tt.stop();
//    cout<<"Time for partition graph generation: "<<tt.GetRuntime()<<" s."<<endl;

    vector<int> leafNodes;
    for(int i=0;i<GTree.size();++i){
        if(GTree[i].isleaf){
            leafNodes.push_back(i);
        }
    }
    cout<<"Leaf node number (partition number): "<<leafNodes.size()<<endl;

    tt.start();
    /// It seems that G-tree is an unbalanced tree
    // bottom up calculation
    // temp graph
//    vector<Node> graph;
//    graph = Nodes;//original graph
    vector<int> cands;
    vector<int> result;
    unordered_map<int, unordered_map<int,int> > vertex_pairs;//result of distance matrix

    // do dijkstra
    int s, t, tn, nid, cid, weight;
    vector<int> tnodes, tweight;
    set<int> nset;

    int partiNum = leafNodes.size();

    if(ifParallel){/// multi-thread
//        vSmNode.reserve(node_num);
//        for(int i = 0; i < node_num; i++)
//        {
//            Semaphore* s = new Semaphore(1);
//            vSmNode.push_back(s);
//        }
        cout<<"Multi-thread computation for partition index."<<endl;
        if(partiNum>thread_num){
            int step=partiNum/thread_num;
            boost::thread_group threadf;
            for(int i=0;i<thread_num;i++){
                pair<int,int> p;
                p.first=i*step;
                if(i==thread_num-1){
                    p.second=partiNum;
                }
                else{
                    p.second=(i+1)*step;
                }
                threadf.add_thread(new boost::thread(&Gstartree::leafNodeDisMatrixCompute, this, p, boost::ref(leafNodes), false));
            }
            threadf.join_all();
        }else{
            boost::thread_group threadf;
            for(int pid=0;pid<partiNum;++pid) {
                threadf.add_thread(new boost::thread(&Gstartree::leafNodeDisMatrixCompute, this, make_pair(pid,pid+1), boost::ref(leafNodes), false));
            }
            threadf.join_all();
        }

//        if(!vSmNode.empty()){
//            for(int i = 0; i < node_num; i++)
//                delete vSmNode[i];
//            vSmNode.clear();
//        }
    }
    else{/// single thread
        /// Construct Li
        cout<<"Single thread computation for partition index."<<endl;
        for(int i=0;i<leafNodes.size();++i){
            tn = leafNodes[i];
            // cands = leafnodes
            cands = GTree[tn].leafnodes;//for leaf node, the cands is the leafnodes
            for(int j=0;j<cands.size();++j){
                Nodes[cands[j]].inleaf = tn;
                Nodes[cands[j]].inleafpos = j;
            }
            // union borders = borders;
            GTree[tn].union_borders = GTree[tn].borders;

            // start to do min dis
            vertex_pairs.clear();
            // for each border, do min dis
            int cc = 0;
            for ( int k = 0; k < GTree[tn].union_borders.size(); k++ ){
                //printf("DIJKSTRA...LEAF=%d BORDER=%d\n", tn, GTree[tn].union_borders[k] );
//                result = dijkstra_candidate( GTree[tn].union_borders[k], cands, graph );//distance vector from s to all borders
//                result = dijkstra_candidate_No( GTree[tn].union_borders[k], cands, PartiGraph[tn] );//distance vector from s to all borders
                result = dijkstra_candidate_No( GTree[tn].union_borders[k], cands, Nodes);//distance vector from s to all borders
                //printf("DIJKSTRA...END\n");

                // save to map
                for ( int p = 0; p < result.size(); p ++ ){
                    GTree[tn].mind.push_back( result[p] );//store the distance vector of s
                    vertex_pairs[GTree[tn].union_borders[k]][cands[p]] = result[p];
                }
            }


            // second, add inter connected edges (shortcuts)
            for ( int k = 0; k < GTree[tn].borders.size(); k++ ){
                s = GTree[tn].borders[k];
                for ( int p = 0; p < GTree[tn].borders.size(); p++ ){
                    if ( k == p ) continue;

                    t = GTree[tn].borders[p];
//                    graph[s].adjnodes.push_back( t );
//                    graph[s].adjweight.push_back( vertex_pairs[s][t] );
                    /// For Overlay graph
                    if(OverlayGraph[s].find(t)==OverlayGraph[s].end()){//if not found
                        OverlayGraph[s].insert({t,vertex_pairs[s][t]});
                    }else{//if found
                        if(OverlayGraph[s][t] > vertex_pairs[s][t]){//if equal
                            OverlayGraph[s][t] = vertex_pairs[s][t];
                        }
                        else if(OverlayGraph[s][t] < vertex_pairs[s][t]){
                            cout<<"Wrong for this case! "<<OverlayGraph[s][t]<<" "<<vertex_pairs[s][t]<<endl; exit(1);
                        }
                    }

                }
            }
        }
    }

    tt.stop();
    cout<<"The time used for leaf node index generation: "<<tt.GetRuntime()<<" s."<<endl;

//    for(int i=0;i<vSmNode.size();++i){
//        delete vSmNode[i];
//    }
}
//calculate the distance matrix of leaf nodes by the pre-boundary strategy
void Gstartree::hierarchy_shortest_path_calculation_PreBoundary(bool ifParallel){
    //get partition graph
    Timer tt;
//    tt.start();
//    getPartiGraph();
//    tt.stop();
//    cout<<"Time for partition graph generation: "<<tt.GetRuntime()<<" s."<<endl;

    vector<int> leafNodes;
    for(int i=0;i<GTree.size();++i){
        if(GTree[i].isleaf){
            leafNodes.push_back(i);
        }
    }
    cout<<"Leaf node number (partition number): "<<leafNodes.size()<<endl;

    tt.start();
    /// It seems that G-tree is an unbalanced tree
    // bottom up calculation
    // temp graph
//    vector<Node> graph;
//    graph = Nodes;//original graph
    vector<int> cands;
    vector<int> result;
    unordered_map<int, unordered_map<int,int> > vertex_pairs;//result of distance matrix

    // do dijkstra
    int s, t, tn, nid, cid, weight;
    vector<int> tnodes, tweight;
    set<int> nset;

    int partiNum = leafNodes.size();

    if(ifParallel){/// multi-thread
//        vSmNode.reserve(node_num);
//        for(int i = 0; i < node_num; i++)
//        {
//            Semaphore* s = new Semaphore(1);
//            vSmNode.push_back(s);
//        }
        cout<<"Multi-thread computation for partition index."<<endl;
        if(partiNum>thread_num){
            int step=partiNum/thread_num;
            boost::thread_group threadf;
            for(int i=0;i<thread_num;i++){
                pair<int,int> p;
                p.first=i*step;
                if(i==thread_num-1)
                    p.second=partiNum;
                else
                    p.second=(i+1)*step;
//                threadf.add_thread(new boost::thread(&Gstartree::leafNodeAllPair_No, this, p, boost::ref(leafNodes),boost::ref(PartiGraph)));
                threadf.add_thread(new boost::thread(&Gstartree::leafNodeDisMatrixCompute, this, p, boost::ref(leafNodes), true));
            }
            threadf.join_all();
        }else{
            boost::thread_group threadf;
            for(int pid=0;pid<partiNum;++pid) {
//                threadf.add_thread(new boost::thread(&Gstartree::leafNodeAllPair_No, this, make_pair(pid,pid+1), boost::ref(leafNodes),boost::ref(PartiGraph)));
                threadf.add_thread(new boost::thread(&Gstartree::leafNodeDisMatrixCompute, this, make_pair(pid,pid+1), boost::ref(leafNodes), true));
            }
            threadf.join_all();
        }

//        if(!vSmNode.empty()){
//            for(int i = 0; i < node_num; i++)
//                delete vSmNode[i];
//            vSmNode.clear();
//        }
    }
    else{/// single thread
        /// Construct Li
        cout<<"Single thread computation for partition index."<<endl;
        for(int i=0;i<leafNodes.size();++i){
            tn = leafNodes[i];
            // cands = leafnodes
            cands = GTree[tn].leafnodes;//for leaf node, the cands is the leafnodes
            for(int j=0;j<cands.size();++j){
                Nodes[cands[j]].inleaf = tn;
                Nodes[cands[j]].inleafpos = j;
            }
            // union borders = borders;
            GTree[tn].union_borders = GTree[tn].borders;

            // start to do min dis
            vertex_pairs.clear();
            // for each border, do min dis
            int cc = 0;
            for ( int k = 0; k < GTree[tn].union_borders.size(); k++ ){
                //printf("DIJKSTRA...LEAF=%d BORDER=%d\n", tn, GTree[tn].union_borders[k] );
//                result = dijkstra_candidate( GTree[tn].union_borders[k], cands, graph );//distance vector from s to all borders
//                result = dijkstra_candidate_No( GTree[tn].union_borders[k], cands, PartiGraph[tn] );//distance vector from s to all borders
                result = dijkstra_candidate( GTree[tn].union_borders[k], cands, Nodes);//distance vector from s to all borders
                //printf("DIJKSTRA...END\n");

                // save to map
                for ( int p = 0; p < result.size(); p ++ ){
                    GTree[tn].mind.push_back( result[p] );//store the distance vector of s
                    vertex_pairs[GTree[tn].union_borders[k]][cands[p]] = result[p];
                }
            }


            // second, add inter connected edges (shortcuts)
            for ( int k = 0; k < GTree[tn].borders.size(); k++ ){
                s = GTree[tn].borders[k];
                for ( int p = 0; p < GTree[tn].borders.size(); p++ ){
                    if ( k == p ) continue;

                    t = GTree[tn].borders[p];
//                    graph[s].adjnodes.push_back( t );
//                    graph[s].adjweight.push_back( vertex_pairs[s][t] );
                    /// For Overlay graph
                    if(OverlayGraph[s].find(t)==OverlayGraph[s].end()){//if not found
                        OverlayGraph[s].insert({t,vertex_pairs[s][t]});
                    }else{//if found
                        if(OverlayGraph[s][t] > vertex_pairs[s][t]){//if equal
                            OverlayGraph[s][t] = vertex_pairs[s][t];
                        }
                        else if(OverlayGraph[s][t] < vertex_pairs[s][t]){
                            cout<<"Wrong for this case! "<<OverlayGraph[s][t]<<" "<<vertex_pairs[s][t]<<endl; exit(1);
                        }
                    }

                }
            }
        }
    }

    tt.stop();
    cout<<"The time used for leaf node index generation: "<<tt.GetRuntime()<<" s."<<endl;

//    for(int i=0;i<vSmNode.size();++i){
//        delete vSmNode[i];
//    }
}

/*void Gstartree::leafNodeAllPair_No(pair<int,int> p, vector<int> & leafNodes, map<int,unordered_map<int,vector<pair<int,int>>>> & graphs){
    vector<int> cands;
    vector<int> result;
    unordered_map<int, unordered_map<int,int> > vertex_pairs;//result of distance matrix
    int s, t, tn, nid, cid, weight;
    vector<int> tnodes, tweight;
    set<int> nset;

    for(int i=p.first;i<p.second;++i){
        int tn = leafNodes[i];
//        unordered_map<int,vector<pair<int,int>>> graph = graphs[tn];
//        unordered_map<int,vector<pair<int,int>>> graph = PartiGraph[tn];
        // cands = leafnodes
        cands = GTree[tn].leafnodes;//for leaf node, the cands is the leafnodes
        for(int j=0;j<cands.size();++j){
            Nodes[cands[j]].inleaf = tn;
            Nodes[cands[j]].inleafpos = j;
        }
        // union borders = borders;
        GTree[tn].union_borders = GTree[tn].borders;

        // start to do min dis
        vertex_pairs.clear();
        // for each border, do min dis
        int cc = 0;
        for ( int k = 0; k < GTree[tn].union_borders.size(); k++ ){
            //printf("DIJKSTRA...LEAF=%d BORDER=%d\n", tn, GTree[tn].union_borders[k] );
//            result = dijkstra_candidate( GTree[tn].union_borders[k], cands, graph );//distance vector from s to all borders
//            result = dijkstra_candidate_No( GTree[tn].union_borders[k], cands,  PartiGraph[tn] );//distance vector from s to all borders
            result = dijkstra_candidate_No( GTree[tn].union_borders[k], cands,  Nodes );//distance vector from s to all borders
            //printf("DIJKSTRA...END\n");

            // save to map
            for ( int p = 0; p < result.size(); p ++ ){
                GTree[tn].mind.push_back( result[p] );//store the distance vector of s
                vertex_pairs[GTree[tn].union_borders[k]][cands[p]] = result[p];
            }
        }

        // IMPORTANT! after all border finished, degenerate graph
        // first, remove inward edges
        for ( int k = 0; k < GTree[tn].borders.size(); k++ ){//for each border vertex
            s = GTree[tn].borders[k];
            tnodes.clear();
            tweight.clear();
            for ( int p = 0; p < Nodes[s].adjnodes.size(); p++ ){//for each adjacent vertex
                nid = Nodes[s].adjnodes[p];
                weight = Nodes[s].adjweight[p];
                // if adj node in same tree node

                if(BoundaryTag[nid]){
                    vSmNode[s]->wait();
                    E[s].insert({nid,make_pair(weight,1)});
                    vSmNode[s]->notify();
                }

//                if ( graph[nid].gtreepath.size() <= i || graph[nid].gtreepath[i] != tn ){// add the higher-level nodes or other node in the same level
//                    // only remain those useful vertices, i.e., borders
//                    tnodes.push_back(nid);
//                    tweight.push_back(weight);
//                }
            }
            // cut it
//                graph[s].adjnodes = tnodes;//update the adjacency lists of graph, only left the useful boundary vertices
//                graph[s].adjweight = tweight;
            /// For NeighborCon
//            for(int p=0;p<tnodes.size();++p){
//                E[s].insert({tnodes[p],make_pair(tweight[p],1)});
//            }

        }
        // second, add inter connected edges (shortcuts)
        for ( int k = 0; k < GTree[tn].borders.size(); k++ ){
            for ( int p = 0; p < GTree[tn].borders.size(); p++ ){
                if ( k == p ) continue;
                s = GTree[tn].borders[k];
                t = GTree[tn].borders[p];
//                    graph[s].adjnodes.push_back( t );
//                    graph[s].adjweight.push_back( vertex_pairs[s][t] );
                /// For NeighborCon
                vSmNode[s]->wait();
                E[s].insert({t,make_pair(INF,1)});
                vSmNode[s]->notify();

                if(vertex_pairs[s][t]<INF){
                    if(E[s].find(t) != E[s].end()){//if found
                        if(E[s][t].first == vertex_pairs[s][t]){//if equal
                            E[s][t].second += 1;
                        }else{//if not equal
                            assert(E[s][t].first > vertex_pairs[s][t]);
                            E[s][t].first = vertex_pairs[s][t];
                            E[s][t].second = 1;
                        }
                    }else{//if not found
                        vSmNode[s]->wait();
                        E[s].insert({t,make_pair(vertex_pairs[s][t],1)});
                        vSmNode[s]->notify();
                    }
                }
//                else{
//                    cout<<"INF! "<<s<<" "<<t<<" "<<vertex_pairs[s][t]<<endl;
//                }

            }
        }

    }

}*/

void Gstartree::leafNodeDisMatrixCompute(pair<int,int> p, vector<int> & leafNodes, bool ifGlobal){
    vector<int> cands;
    vector<int> result;
    unordered_map<int, unordered_map<int,int> > vertex_pairs;//result of distance matrix
    int s, t, tn, nid, cid, weight;
    vector<int> tnodes, tweight;
    set<int> nset;

    for(int i=p.first;i<p.second;++i){
        int tn = leafNodes[i];
//        unordered_map<int,vector<pair<int,int>>> graph = graphs[tn];
//        unordered_map<int,vector<pair<int,int>>> graph = PartiGraph[tn];
        // cands = leafnodes
        cands = GTree[tn].leafnodes;//for leaf node, the cands is the leafnodes
        for(int j=0;j<cands.size();++j){
            Nodes[cands[j]].inleaf = tn;
            Nodes[cands[j]].inleafpos = j;
        }
        // union borders = borders;
        GTree[tn].union_borders = GTree[tn].borders;

        // start to do min dis
        vertex_pairs.clear();
        // for each border, do min dis
        int cc = 0;
        if(ifGlobal){//pre-boundary strategy
            for ( int k = 0; k < GTree[tn].union_borders.size(); k++ ){
                //printf("DIJKSTRA...LEAF=%d BORDER=%d\n", tn, GTree[tn].union_borders[k] );
                result = dijkstra_candidate( GTree[tn].union_borders[k], cands,  Nodes );//distance vector from s to all borders
                //printf("DIJKSTRA...END\n");

                // save to map
                for ( int p = 0; p < result.size(); p ++ ){
                    GTree[tn].mind.push_back( result[p] );//store the distance vector of s
                    vertex_pairs[GTree[tn].union_borders[k]][cands[p]] = result[p];
                }
            }
        }
        else{//no-boundary strategy
            for ( int k = 0; k < GTree[tn].union_borders.size(); k++ ){
                //printf("DIJKSTRA...LEAF=%d BORDER=%d\n", tn, GTree[tn].union_borders[k] );
                result = dijkstra_candidate_No( GTree[tn].union_borders[k], cands,  Nodes );//distance vector from s to all borders
                //printf("DIJKSTRA...END\n");

                // save to map
                for ( int p = 0; p < result.size(); p ++ ){
                    GTree[tn].mind.push_back( result[p] );//store the distance vector of s
                    vertex_pairs[GTree[tn].union_borders[k]][cands[p]] = result[p];
                }
            }
        }

        // second, add inter connected edges (shortcuts)
        for ( int k = 0; k < GTree[tn].borders.size(); k++ ){
            for ( int p = 0; p < GTree[tn].borders.size(); p++ ){
                if ( k == p ) continue;
                s = GTree[tn].borders[k];
                t = GTree[tn].borders[p];
//                    graph[s].adjnodes.push_back( t );
//                    graph[s].adjweight.push_back( vertex_pairs[s][t] );
                if(OverlayGraph[s].find(t)==OverlayGraph[s].end()){//if not found
                    OverlayGraph[s].insert({t,vertex_pairs[s][t]});
                }else{//if found
                    if(OverlayGraph[s][t] > vertex_pairs[s][t]){//if equal
                        OverlayGraph[s][t] = vertex_pairs[s][t];
                    }
                    else if(OverlayGraph[s][t] < vertex_pairs[s][t]){
                        cout<<"Wrong for this case! "<<OverlayGraph[s][t]<<" "<<vertex_pairs[s][t]<<endl; exit(1);
                    }
                }

            }
        }

    }

}

void Gstartree::H2HconOrderMT(bool ifMDE){
    Timer tt;
    tt.start();
    if(ifMDE){
        MDEContract();//MDE-based contraction
    }
    else{
        getOverlayLayers();//get boundary layers
        StratifiedMDEOrdering(boundaryLevel);
//        VertexContraction();//Order-based contraction
    }
    tt.stop();
    cout<<"MDE contraction time: "<<tt.GetRuntime()<<" s."<<endl;
    makeTree();
    makeIndex();

}

//Function of contracting vertices by MDE
void Gstartree::MDEContract(){
    Timer tt;
    tt.start();
    cout<<"Start MDE-based contraction..."<<endl;

    NeighborCon.assign(node_num,vector<pair<int,pair<int,int>>>());
    //get original overlay graph
    E.assign(node_num,unordered_map<int,pair<int,int>>());
    for(int i=0;i<OverlayGraph.size();i++){
        if(!OverlayGraph[i].empty()){
            for(auto it=OverlayGraph[i].begin();it!=OverlayGraph[i].end();++it){
                E[i].insert(make_pair(it->first,make_pair(it->second,1)));
            }
        }
    }
    //for H2H update
    SCconNodesMT.assign(node_num, map<int, vector<int>>());

    _DD_.assign(node_num,0);
//    _DD2_.assign(node_num,0);
//    DD.assign(node_num,0); DD2.assign(node_num,0);

    set<DegComp1> Deg;//min first
//    vector<bool> active(node_num,false);//flag that indicate whether a vertex is active for contraction
    int degree;
    vector<bool> exist(node_num,false);//if in the core, all vertices is originally in core

    for(int ID=0;ID<OverlayGraph.size();++ID){
        if(Nodes[ID].isborder){
            degree=E[ID].size();
            exist[ID] = true;
            if(degree > 0){//get degree
                _DD_[ID]=degree;
//            _DD2_[ID]=degree;
                Deg.insert(DegComp1(ID));
//            active[i] = true;
            }else{
                cout<<"Wrong!! Degree of "<<ID<<" is "<<degree<<endl;
                exit(1);
            }
        }
    }

    vNodeOrder.clear();
    for(int i=0;i<node_num;++i){
        if(!Nodes[i].isborder){//if not border
            vNodeOrder.push_back(i);
        }
    }
//    while(vNodeOrder.size()<node_num-OverlayGraph.size()){
//        vNodeOrder.push_back(-1);
//    }
    vector<bool> change(node_num,false);//whether the neighbor (degree) has changed

    bool CutLabel=false;
    int count=0;
    int ID1,ID2;

//    cout<<"192402: "<<E[192402].size()<<" "<<Nodes[192402].inleaf<<" "<<Nodes[192402].isborder<<endl;
//    for(auto i1=E[192402].begin();i1!=E[192402].end();++i1){
//        cout<<"E "<<i1->first<<" "<<Nodes[i1->first].inleaf<<" "<<Nodes[i1->first].isborder<<endl;;
//    }
//    for(auto i1=Nodes[192402].adjnodes.begin();i1!=Nodes[192402].adjnodes.end();++i1){
//        cout<<"Adj "<<*i1<<" "<<Nodes[*i1].inleaf<<" "<<Nodes[*i1].isborder<<endl;
//    }

    //Get the order of all vertices by MDE
    while(!Deg.empty()){
        count+=1;
        int x=(*Deg.begin()).x;//minimum degree first
        while(change[x]){//update the degree if it is changed
            Deg.erase(DegComp1(x));
//            Deg.erase(Deg.begin());
            _DD_[x]=E[x].size();
//            _DD2_[x]=E[x].size();
            Deg.insert(DegComp1(x));
            change[x]=false;
            x=(*Deg.begin()).x;
        }
//        cout<<x<<" "<<E[x].size()<<endl;
        vNodeOrder.push_back(x);//least important vertex first
        Deg.erase(Deg.begin());

//        vector<pair<int,int>> Neigh; Neigh.clear();
//        for(auto it=E[x].begin();it!=E[x].end();it++){//for the neighbors of x
//            if(existCore[(*it).first]){//if in the core
//                Neigh.push_back(*it);
//            }
//        }
//        NeighborConCore[x].assign(Neigh.begin(),Neigh.end());//for core vertex, the neighbors are all truly-core vertices
        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();
        for(auto it=E[x].begin();it!=E[x].end();it++){
            if(exist[(*it).first]){
                Neigh.push_back(*it);
            }
        }
        if(Neigh.empty() && !Deg.empty()){
            cout<<"!!! Neigh is empty for "<<x<<" : "<<E[x].size()<<" "<<Nodes[x].adjnodes.size()<<" "<<count<<endl;
//            exit(1);
        }
        NeighborCon[x].assign(Neigh.begin(),Neigh.end());

        exist[x]=false;
        //delete the star
        for(int i=0;i<Neigh.size();i++){
            int y=Neigh[i].first;
            deleteE(x,y);//delete x from y's adjacency list
            change[y]=true;
        }
        //add all-pair neighbors
        if(Neigh.size()<=100){
            for(int i=0;i<Neigh.size();i++) {
                ID1 = Neigh[i].first;
                for (int j = i + 1; j < Neigh.size(); j++) {
                    ID2 = Neigh[j].first;
                    insertE(ID1, ID2, Neigh[i].second.first + Neigh[j].second.first);
                    /// For TD update

                    if(Neigh[i].first<Neigh[j].first)
                        SCconNodesMT[Neigh[i].first][Neigh[j].first].push_back(x);//no direction
                    else if(Neigh[j].first<Neigh[i].first)
                        SCconNodesMT[Neigh[j].first][Neigh[i].first].push_back(x);

                }
            }
        }
        else{
            if(Neigh.size()>thread_num){
                int step=Neigh.size()/thread_num;
                boost::thread_group thread;
                for(int i=0;i<thread_num;i++){
                    pair<int,int> p;
                    p.first=i*step;
                    if(i==thread_num-1)
                        p.second=Neigh.size();
                    else
                        p.second=(i+1)*step;
                    thread.add_thread(new boost::thread(&Gstartree::NeighborComorder, this, boost::ref(Neigh), p, x));
                }
                thread.join_all();
            }else{
                boost::thread_group thread;
                for(int i=0;i<Neigh.size();i++){
                    pair<int,int> p;
                    p.first=i; p.second=(i+1);
                    thread.add_thread(new boost::thread(&Gstartree::NeighborComorder, this, boost::ref(Neigh), p, x));
                }
                thread.join_all();
            }
        }

    }
    tt.stop();
    cout<<"The time used for contraction: "<<tt.GetRuntime()<<" s."<<endl;

    if(vNodeOrder.size()!=node_num){
        cout<<"vNodeOrder size incorrect! "<<vNodeOrder.size()<<endl;
        exit(1);
    }
    NodeOrder.assign(node_num,-1);
    for(int k=0;k<vNodeOrder.size();k++){
        ID1=vNodeOrder[k];
        if(ID1>=0 && ID1<node_num){
            NodeOrder[ID1]=k;
        }else{
            cout<<"Wrong ID! "<<ID1<<" "<<k<<endl; exit(1);
        }

    }
//    cout << "???" << endl;
//    while(vNodeOrder.size()<node_num){
//        vNodeOrder.push_back(-1);
//    }

}

void Gstartree::StratifiedMDEOrdering(vector<vector<int>> & boundaryLevel){
    NeighborCon.assign(node_num,vector<pair<int,pair<int,int>>>());
    SCconNodesMT.assign(node_num, map<int, vector<int>>());
    _DD_.assign(node_num,0);

    set<DegComp1> Deg;//min first
//    vector<bool> active(node_num,false);//flag that indicate whether a vertex is active for contraction
    int degree;
    vector<bool> exist(node_num,false);//if in the core, all vertices is originally in core

    E.assign(node_num,unordered_map<int,pair<int,int>>());
    for(int i=0;i<OverlayGraph.size();i++){
        if(!OverlayGraph[i].empty()){
            for(auto it=OverlayGraph[i].begin();it!=OverlayGraph[i].end();++it){
                E[i].insert(make_pair(it->first,make_pair(it->second,1)));
            }
            if(Nodes[i].isborder){
                exist[i] = true;
//                degree=E[i].size();
//                if(degree > 0){//get degree
//                    _DD_[i]=degree;
//                }else{
//                    cout<<"Wrong!! Degree of "<<i<<" is "<<degree<<endl;
//                    exit(1);
//                }
            }
            else{
                cout<<"Wrong for this vertex! "<<i<<endl; exit(1);
            }
        }
    }

    vNodeOrder.clear();
    for(int i=0;i<node_num;++i){
        if(!Nodes[i].isborder){//if not border
            vNodeOrder.push_back(i);
        }
    }
    cout<<"overlay vertex number: "<<node_num-vNodeOrder.size()<<endl;

    vector<bool> change(node_num,false);//whether the neighbor (degree) has changed

    bool CutLabel=false;
    int count=0;
    int ID1,ID2;
    //Get the order of all vertices by MDE
    Timer tt;
    cout<<"Start Hierarchical MDE-based contraction..."<<endl;
    for(int layer_i=boundaryLevel.size()-1;layer_i>=0;--layer_i){
        Deg.clear();
        change.assign(node_num,false);
        for(int j=0;j<boundaryLevel[layer_i].size();++j){
            ID1=boundaryLevel[layer_i][j];
            if(Nodes[ID1].isborder){
                degree=E[ID1].size();
                if(degree > 0){//get degree
                    _DD_[ID1]=degree;
                    Deg.insert(DegComp1(ID1));
                }else{
                    cout<<"Wrong!! Degree of "<<ID1<<" is "<<degree<<endl;
                    exit(1);
                }
            }
            else{
                cout<<"Wrong for this vertex! "<<ID1<<endl; exit(1);
            }
        }

        tt.start();
        while(!Deg.empty()){
            count+=1;
            int x=(*Deg.begin()).x;//minimum degree first
            while(change[x]){//update the degree if it is changed
                Deg.erase(DegComp1(x));
//            Deg.erase(Deg.begin());
                _DD_[x]=E[x].size();
//            _DD2_[x]=E[x].size();
                Deg.insert(DegComp1(x));
                change[x]=false;
                x=(*Deg.begin()).x;
            }
//        cout<<x<<" "<<E[x].size()<<endl;
            vNodeOrder.push_back(x);//least important vertex first
            Deg.erase(Deg.begin());

            vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();
            for(auto it=E[x].begin();it!=E[x].end();it++){
                if(exist[(*it).first]){
                    Neigh.push_back(*it);
                }
            }
            if(Neigh.empty() && !Deg.empty()){
                cout<<"!!! Neigh is empty for "<<x<<" : "<<E[x].size()<<" "<<Nodes[x].adjnodes.size()<<" "<<count<<endl;
//            exit(1);
            }
            NeighborCon[x].assign(Neigh.begin(),Neigh.end());

            exist[x]=false;
            //delete the star
            for(int i=0;i<Neigh.size();i++){
                int y=Neigh[i].first;
                deleteE(x,y);//delete x from y's adjacency list
                change[y]=true;
            }
            //add all-pair neighbors
            if(Neigh.size()<=100){
                for(int i=0;i<Neigh.size();i++) {
                    ID1 = Neigh[i].first;
                    for (int j = i + 1; j < Neigh.size(); j++) {
                        ID2 = Neigh[j].first;
                        insertE(ID1, ID2, Neigh[i].second.first + Neigh[j].second.first);
                        /// For TD update
                        if(Neigh[i].first<Neigh[j].first)
                            SCconNodesMT[Neigh[i].first][Neigh[j].first].push_back(x);//no direction
                        else if(Neigh[j].first<Neigh[i].first)
                            SCconNodesMT[Neigh[j].first][Neigh[i].first].push_back(x);

                    }
                }
            }
            else{
                if(Neigh.size()>thread_num){
                    int step=Neigh.size()/thread_num;
                    boost::thread_group thread;
                    for(int i=0;i<thread_num;i++){
                        pair<int,int> p;
                        p.first=i*step;
                        if(i==thread_num-1)
                            p.second=Neigh.size();
                        else
                            p.second=(i+1)*step;
                        thread.add_thread(new boost::thread(&Gstartree::NeighborComorder, this, boost::ref(Neigh), p, x));
                    }
                    thread.join_all();
                }else{
                    boost::thread_group thread;
                    for(int i=0;i<Neigh.size();i++){
                        pair<int,int> p;
                        p.first=i; p.second=(i+1);
                        thread.add_thread(new boost::thread(&Gstartree::NeighborComorder, this, boost::ref(Neigh), p, x));
                    }
                    thread.join_all();
                }
            }


        }
        tt.stop();
        cout<<"Layer "<<layer_i<<" : "<<boundaryLevel[layer_i].size()<<" ; time: "<<tt.GetRuntime()<<" s."<<endl;
    }

    if(vNodeOrder.size()!=node_num){
        cout<<"vNodeOrder size incorrect! "<<vNodeOrder.size()<<endl;
        exit(1);
    }
    NodeOrder.assign(node_num,-1);
    for(int k=0;k<vNodeOrder.size();k++){
        ID1=vNodeOrder[k];
        if(ID1>=0 && ID1<node_num){
            NodeOrder[ID1]=k;
        }else{
            cout<<"Wrong ID! "<<ID1<<" "<<k<<endl; exit(1);
        }
    }

}

//Function of contracting vertices by given order
void Gstartree::VertexContraction(){
    Timer tt;
    tt.start();
    cout<<"Start contraction..."<<endl;

    //for H2H update
    vector<bool> exist(node_num,false);//if in the core, all vertices is originally in core
    SCconNodesMT.assign(node_num, map<int, vector<int>>());
    NeighborCon.assign(node_num,vector<pair<int,pair<int,int>>>());//temporal graph to store Neighbors in the core, for graph contraction
    //get original overlay graph
    E.assign(node_num,unordered_map<int,pair<int,int>>());
    for(int i=0;i<OverlayGraph.size();i++){
        if(!OverlayGraph[i].empty()){
            exist[i]=true;
            for(auto it=OverlayGraph[i].begin();it!=OverlayGraph[i].end();++it){
                E[i].insert(make_pair(it->first,make_pair(it->second,1)));
            }
        }
    }

    _DD_.assign(node_num,0);
//    _DD2_.assign(node_num,0);
//    DD.assign(node_num,0); DD2.assign(node_num,0);

//    vector<bool> active(node_num,false);//flag that indicate whether a vertex is active for contraction
    int degree;



//    while(vNodeOrder.size()<node_num-OverlayGraph.size()){
//        vNodeOrder.push_back(-1);
//    }
    vector<bool> change(node_num,false);//whether the neighbor (degree) has changed



    bool CutLabel=false;
    int count=0;
    int ID1,ID2;

//    cout<<"192402: "<<E[192402].size()<<" "<<Nodes[192402].inleaf<<" "<<Nodes[192402].isborder<<endl;
//    for(auto i1=E[192402].begin();i1!=E[192402].end();++i1){
//        cout<<"E "<<i1->first<<" "<<Nodes[i1->first].inleaf<<" "<<Nodes[i1->first].isborder<<endl;;
//    }
//    for(auto i1=Nodes[192402].adjnodes.begin();i1!=Nodes[192402].adjnodes.end();++i1){
//        cout<<"Adj "<<*i1<<" "<<Nodes[*i1].inleaf<<" "<<Nodes[*i1].isborder<<endl;
//    }

    //Get the order of all vertices by MDE
    for(int i=node_num-overlayNodeNum;i<node_num;++i){
        count+=1;
        int x=vNodeOrder[i];
        assert(x>=0 && x<node_num);
        if(!Nodes[x].isborder){
            cout<<"wrong for border! "<<x<<endl; exit(1);
        }
        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();
        for(auto it=E[x].begin();it!=E[x].end();it++){
            if(exist[(*it).first]){
                Neigh.push_back(*it);
            }
        }

        NeighborCon[x].assign(Neigh.begin(),Neigh.end());

        exist[x]=false;
        //delete the star
        for(int i=0;i<Neigh.size();i++){
            int y=Neigh[i].first;
            deleteE(x,y);//delete x from y's adjacency list
            change[y]=true;
        }
        //add all-pair neighbors
        for(int i=0;i<Neigh.size();i++) {
            ID1 = Neigh[i].first;
            for (int j = i + 1; j < Neigh.size(); j++) {
                ID2 = Neigh[j].first;
                insertE(ID1, ID2, Neigh[i].second.first + Neigh[j].second.first);
                /// For TD update

                if(Neigh[i].first<Neigh[j].first)
                    SCconNodesMT[Neigh[i].first][Neigh[j].first].push_back(x);//no direction
                else if(Neigh[j].first<Neigh[i].first)
                    SCconNodesMT[Neigh[j].first][Neigh[i].first].push_back(x);

                change[ID1] = true;
                change[ID2] = true;
            }
        }

    }
    tt.stop();
    cout<<"The time used for contraction: "<<tt.GetRuntime()<<" s."<<endl;


//    cout << "???" << endl;
//    while(vNodeOrder.size()<node_num){
//        vNodeOrder.push_back(-1);
//    }

}

void Gstartree::makeTree(){
    vector<int> vecemp; //vecemp.clear();
    VidtoTNid.assign(node_num,vecemp);

    rank.assign(node_num,0);
    //Tree.clear();
    int len=vNodeOrder.size()-1;
    heightMax=0;
    if(vNodeOrder.size()!=node_num){
        cout<<"vNodeOrder size wrong! "<<vNodeOrder.size()<<" "<<node_num<<endl; exit(1);
    }

    TDNode rootn;
    int x=vNodeOrder[len];
    //cout<<"len "<<len<<" , ID "<<x<<endl;
    assert(x>=0 && x<node_num);
    rootn.vert=NeighborCon[x];
    rootn.uniqueVertex=x;
    rootn.pa=-1;
    rootn.height=1;
    rank[x]=0;
    Tree.push_back(rootn);
    len--;


    int nn;
    for(;len>=0;len--){
//        if(len%100==0){
//            cout<<"len: "<<len<<endl;
//        }
        x=vNodeOrder[len];
        assert(x>=0 && x<node_num);
        if(!Nodes[x].isborder){//if not boundary vertex
            break;
        }
        if(treeWidth<NeighborCon[x].size()){
            treeWidth=NeighborCon[x].size();
        }
        TDNode nod;
        nod.vert=NeighborCon[x];
        nod.uniqueVertex=x;
        int pa=match(x,NeighborCon[x]);
        Tree[pa].ch.push_back(Tree.size());
        nod.pa=pa;
        nod.height=Tree[pa].height+1;

        nod.hdepth=Tree[pa].height+1;
        for(int i=0;i<NeighborCon[x].size();i++){
            nn=NeighborCon[x][i].first;
            VidtoTNid[nn].push_back(Tree.size());
            if(Tree[rank[nn]].hdepth<Tree[pa].height+1)
                Tree[rank[nn]].hdepth=Tree[pa].height+1;
        }
        if(nod.height>heightMax)
            heightMax=nod.height;
        rank[x]=Tree.size();
        Tree.push_back(nod);
//        cout<<"len "<<len<<" , ID "<<x<<" ; "<<NeighborCon[x].size()<<endl;
//        if(len == 263982){
//            cout<<len<<" "<<NeighborCon[x].size()<<endl;
//        }
    }
    cout<<"Tree size: "<<Tree.size()<<" ; Treewidth: "<<treeWidth<<" ; Tree height: "<<heightMax<<endl;
    if(Tree.size()!=overlayNodeNum){
        cout<<"Tree number inconsistent! "<<Tree.size()<<" "<<overlayNodeNum<<endl; exit(1);
    }
}

void Gstartree::makeIndex(){
    makeRMQ();

    //initialize
    vector<int> list; //list.clear();
    list.push_back(Tree[0].uniqueVertex);
    Tree[0].pos.clear();
    Tree[0].pos.push_back(0);

    for(int i=0;i<Tree[0].ch.size();i++){
        makeIndexDFS(Tree[0].ch[i],list);
    }

}

int Gstartree::match(int x,vector<pair<int,pair<int,int>>> &vert) {
    if(vert.empty()){
        cout<<"Empty vert! "<<x<<endl;
        exit(1);
    }
    int nearest = vert[0].first;
    for (int i = 1; i < vert.size(); i++) {
        if (rank[vert[i].first] > rank[nearest])
            nearest = vert[i].first;
    }
    int p = rank[nearest];
    return p;
}

void Gstartree::makeRMQDFS(int p, int height){
    toRMQ[p] = EulerSeq.size();
    EulerSeq.push_back(p);
    for (int i = 0; i < Tree[p].ch.size(); i++){
        makeRMQDFS(Tree[p].ch[i], height + 1);
        EulerSeq.push_back(p);
    }
}

void Gstartree::makeRMQ(){
    //EulerSeq.clear();
    toRMQ.assign(node_num,0);
    //RMQIndex.clear();
    makeRMQDFS(0, 1);
    RMQIndex.push_back(EulerSeq);

    int m = EulerSeq.size();
    for (int i = 2, k = 1; i < m; i = i * 2, k++){
        vector<int> tmp;
        //tmp.clear();
        tmp.assign(m,0);
        for (int j = 0; j < m - i; j++){
            int x = RMQIndex[k - 1][j], y = RMQIndex[k - 1][j + i / 2];
            if (Tree[x].height < Tree[y].height)
                tmp[j] = x;
            else tmp[j] = y;
        }
        RMQIndex.push_back(tmp);
    }
}

int Gstartree::LCAQuery(int _p, int _q){
    int p = toRMQ[_p], q = toRMQ[_q];
    if (p > q){
        int x = p;
        p = q;
        q = x;
    }
    int len = q - p + 1;
    int i = 1, k = 0;
    while (i * 2 < len){
        i *= 2;
        k++;
    }
    q = q - i + 1;
    if (Tree[RMQIndex[k][p]].height < Tree[RMQIndex[k][q]].height)
        return RMQIndex[k][p];
    else return RMQIndex[k][q];
}

void Gstartree::makeIndexDFS(int p, vector<int>& list){
    //initialize
    int NeiNum=Tree[p].vert.size();
    Tree[p].pos.assign(NeiNum+1,0);
    Tree[p].dis.assign(list.size(),INF);
    Tree[p].cnt.assign(list.size(),0);
    Tree[p].FN.assign(list.size(),true);
    Tree[p].vAncestor=list;

    //pos
    //map<int,Nei> Nmap; Nmap.clear();//shortcut infor ordered by the pos ID
    for(int i=0;i<NeiNum;i++){
        for(int j=0;j<list.size();j++){
            if(Tree[p].vert[i].first==list[j]){
                Tree[p].pos[i]=j;
                Tree[p].dis[j]=Tree[p].vert[i].second.first;
                Tree[p].cnt[j]=1;
                break;
            }
        }
    }
    Tree[p].pos[NeiNum]=list.size();

    //dis
    for(int i=0;i<NeiNum;i++){
        int x=Tree[p].vert[i].first;
        int disvb=Tree[p].vert[i].second.first;
        int k=Tree[p].pos[i];//the kth ancestor is x

        for(int j=0;j<list.size();j++){
            int y=list[j];//the jth ancestor is y

            int z;//the distance from x to y
            if(k!=j){
                if(k<j)
                    z=Tree[rank[y]].dis[k];
                else if(k>j)
                    z=Tree[rank[x]].dis[j];

                if(Tree[p].dis[j]>z+disvb){
                    Tree[p].dis[j]=z+disvb;
                    Tree[p].FN[j]=false;
                    Tree[p].cnt[j]=1;
                }else if(Tree[p].dis[j]==z+disvb){
                    Tree[p].cnt[j]+=1;
                }
            }
        }
    }

    //nested loop
    list.push_back(Tree[p].uniqueVertex);
    for(int i=0;i<Tree[p].ch.size();i++){
        makeIndexDFS(Tree[p].ch[i],list);
    }
    list.pop_back();
}


//function of deleting u from v's neighbor
void Gstartree::deleteE(int u,int v){
//    if(E[u].find(v)!=E[u].end()){
//        E[u].erase(E[u].find(v));
//        //DD[u]--;
//    }

    if(E[v].find(u)!=E[v].end()){
        E[v].erase(E[v].find(u));
        //DD[v]--;
    }
}
//function of inserting u to v's neighbor and verse vice.
void Gstartree::insertE(int u,int v,int w){
    if(E[u].find(v)==E[u].end()){
        E[u].insert(make_pair(v,make_pair(w,1)));
        //DD[u]++;
        //DD2[u]++;
    }
    else{
        if(E[u][v].first>w)
            E[u][v]=make_pair(w,1);
        else if(E[u][v].first==w)
            E[u][v].second+=1;
    }

    if(E[v].find(u)==E[v].end()){
        E[v].insert(make_pair(u,make_pair(w,1)));
        //DD[v]++;
        //DD2[v]++;
    }
    else{
        if(E[v][u].first>w)
            E[v][u]=make_pair(w,1);
        else if(E[v][u].first==w)
            E[v][u].second+=1;
    }
}

void Gstartree::insertEMTorder(int u,int v,int w){
    if(E[u].find(v)==E[u].end()){
        E[u].insert(make_pair(v,make_pair(w,1)));
    }
    else{
        if(E[u][v].first>w)
            E[u][v]=make_pair(w,1);
        else if(E[u][v].first==w)
            E[u][v].second+=1;
    }
}

void Gstartree::NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x){
//    sm->wait();
    int ID1, w1;
    int ID2, w2;
    for(int k=p.first;k<p.second;k++){
        ID1=Neighvec[k].first;
        w1=Neighvec[k].second.first;
        for(int h=0;h<Neighvec.size();h++){
            ID2=Neighvec[h].first;
            w2=Neighvec[h].second.first;
            if(ID1==ID2){
                continue;
            }
            insertEMTorder(ID1, ID2, w1+w2);
            if(ID1<ID2){
                if(SCconNodesMT[ID1].find(ID2)==SCconNodesMT[ID1].end()){//if not found
                    SCconNodesMT[ID1].insert({ID2,vector<int>()});
                }
                SCconNodesMT[ID1][ID2].push_back(x);
            }

        }
    }
//    sm->notify();
}

//compute the distance by no-boundary strategy
vector<int> Gstartree::dijkstra_candidate_No( int s, vector<int> &cands, vector<Node> &graph){
    // init
    int PID=Nodes[s].inleaf;

    set<int> todo;
    todo.clear();
    todo.insert(cands.begin(), cands.end());//get the partition vertex set

    map<int,int> result;
    result.clear();
    set<int> visited;
    visited.clear();
//    unordered_map<int,int> q;//map used for mapping vertex id to its distance from source vertex
//    benchmark::heap<2, int, int> q(node_num);
    priority_queue<pair<int,int>,vector<pair<int,int>>,greater<pair<int,int>>> q;
//    vector<bool> closed(node_num, false); //flag vector of whether closed
    vector<Distance> cost(node_num, INF);   //vector of cost
    int temp_dis;
    q.push(make_pair(0,s));
    cost[s] = 0;
//    q[s] = 0;//map of source vertex

    // start
    int min, minpos, adjnode, weight;
    pair<int,int> topE;
    while( ! todo.empty() && ! q.empty() ){
        min = -1;
        topE=q.top();
        min=topE.first; minpos=topE.second;
        q.pop();
        // put min to result, add to visited
        while(visited.find(minpos)!=visited.end()){//if found
            if(!q.empty()){
                topE=q.top();
                min=topE.first; minpos=topE.second;
                q.pop();
            }else{
                cout<<"priority queue empty."<<endl;
                goto Results;
            }
        }

        if(visited.find(minpos)==visited.end()){//if not found
            result[minpos] = min;
            visited.insert( minpos );
        }else{
            cout<<"Wrong!"<<endl; exit(1);
        }

        if ( todo.find( minpos ) != todo.end() ){//if found, erase visited vertex
            todo.erase( minpos );
        }

        // expand on graph (the original graph)
        for ( int i = 0; i < graph[minpos].adjnodes.size(); i++ ){
            adjnode = graph[minpos].adjnodes[i];
//        for ( int i = 0; i < graph[minpos].size(); i++ ){
//            adjnode = graph[minpos][i].first;
            if(Nodes[adjnode].inleaf != PID){
                continue;
            }
            if ( visited.find( adjnode ) != visited.end() ){//if found, ie, it is visited
                continue;
            }
//            weight = graph[minpos][i].second;//edge weight
            weight = graph[minpos].adjweight[i];//edge weight
            temp_dis = min + weight;
            if ( temp_dis < cost[adjnode] ){
                q.push(make_pair(temp_dis,adjnode));
                cost[adjnode] = temp_dis;
            }
        }
    }

    Results:
    // output
    vector<int> output;
    for ( int i = 0; i < cands.size(); i++ ){
        if(result.find(cands[i]) != result.end()){//if found
            output.push_back( result[cands[i]] );//only push the distance result of vertex in cands
        }else{
            output.push_back( INF );//only push the distance result of vertex in cands
        }

    }

    // return
    return output;//vector of distance matrix values
}

void Gstartree::IndexSizeH2H() {
    long long m=0,m1=0,m2=0;

    for(int i=0;i<GTree.size();++i){
        m1+=GTree[i].mind.size()*sizeof(int);
    }

    for(int i=0;i<Tree.size();i++){
        m2+=Tree[i].dis.size()*sizeof(int);//dis
        m2+=Tree[i].pos.size()*sizeof(int);//pos
        m2+=Tree[i].dis.size()*sizeof(int);//cnt
        m2+=Tree[i].vert.size()*3*sizeof(int);//neighID/weight/count
    }

    m=m1+m2;
    cout<<"Partition Index size: "<<(double)m1/1024/1024<<" MB."<<endl;
    cout<<"Overlay Index size: "<<(double)m2/1024/1024<<" MB."<<endl;
    cout<<"Overall Index size: "<<(double)m/1024/1024<<" MB."<<endl;
}


void Gstartree::IndexSizePLL() {
    long long m=0,m1=0,m2=0;

    for(int i=0;i<GTree.size();++i){
        m1+=GTree[i].mind.size()*sizeof(int);
    }

    for(int i=0;i<Label.size();i++){
        m2+=Label[i].size()*2*sizeof(int);
    }

    m=m1+m2;
    cout<<"Partition Index size: "<<(double)m1/1024/1024<<" MB."<<endl;
    cout<<"Overlay Index size: "<<(double)m2/1024/1024<<" MB."<<endl;
    cout<<"Overall Index size: "<<(double)m/1024/1024<<" MB."<<endl;
}



//// Index maintenance

void Gstartree::LGTreeIndexUpdate(vector<pair<pair<int,int>,pair<int,int> > > & updates, bool ifParallel, int updateType){//update(<ID1,ID2>,oldW) vector<pair<pair<int,int>,int> > & updates,
//    double ave_time = 0;
    int ID1,ID2,oldW,newW;
    vector<pair<pair<int,int>,pair<int,int> > > uWithinLeafNode;// edge update within the same leaf node (<ID1,ID2>,<oldW,newW>)
    vector<pair<pair<int,int>,pair<int,int> > > uCrossLeafNode;// edge/shortcut update among borders (<ID1,ID2>,<oldW,newW>)
    bool flag_borderUpdate = false;//whether there is border update of leaf node
//    Timer tt;
//    tt.start();

    for(int i=0;i<updates.size();++i){//for each edge update
        ID1 = updates[i].first.first;
        ID2 = updates[i].first.second;
        oldW = updates[i].second.first;
        newW = updates[i].second.second;

//        cout<<ID1 << " "<<ID2<<" ("<<oldW<<"->"<<newW<<") "<<endl;
        //update edge weight
        for(int j=0;j<Nodes[ID1].adjnodes.size();j++){
            if(Nodes[ID1].adjnodes[j]==ID2){
                assert(Nodes[ID1].adjweight[j] == oldW);
                Nodes[ID1].adjweight[j]=newW;
                break;
            }
        }
        for(int j=0;j<Nodes[ID2].adjnodes.size();j++){
            if(Nodes[ID2].adjnodes[j]==ID1){
                assert(Nodes[ID2].adjweight[j] == oldW);
                Nodes[ID2].adjweight[j]=newW;
                break;
            }
        }
        //identify the edge type
        if(Nodes[ID1].inleaf == Nodes[ID2].inleaf){//if it is an update within the leaf node
            uWithinLeafNode.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
        }else{//if it is an update among the leaf nodes
            uCrossLeafNode.emplace_back(make_pair(ID1,ID2),make_pair(oldW,newW));
            flag_borderUpdate = true;
        }
    }

    /// update processing
    int lnID;
//        int ID1_pos,ID2_pos;
//    vector<Node> graph;//temp graph
//    graph = Nodes;//original graph

    /// Update the leaf nodes
    Timer tt1;
    tt1.start();
    // for update within the same leaf node
    vector<pair<pair<int,int>,pair<int,int> > > testData;
    map<pair<int,int>,pair<int,int>> affectedBPairs;
    if(!uWithinLeafNode.empty()){//if same-leaf node update
        for(auto it=uWithinLeafNode.begin();it!=uWithinLeafNode.end();++it){//deal with each update within leaf node
            ID1 = it->first.first; ID2 = it->first.second; oldW = it->second.first; newW = it->second.second;
            assert(Nodes[ID1].inleaf == Nodes[ID2].inleaf);
            int PID = Nodes[ID1].inleaf;

//            PartitionUpdate_Pre(PID,affectedBPairs,false);
            PartitionUpdate_Pre(PID,affectedBPairs,true);
            //update overlay graph index
            if(!affectedBPairs.empty()){
//                cout<<" affectedBPairs size: "<<affectedBPairs.size()<<endl;
                if(!uCrossLeafNode.empty()){//if there is cross-leaf updates
                    cout<<"Multiple update: in-leaf node and cross-leaf node."<<endl;
                    for(auto it=uCrossLeafNode.begin();it!=uCrossLeafNode.end();++it){
                        oldW=it->second.first; newW=it->second.second;
                        if(it->first.first<it->first.second){
                            ID1=it->first.first, ID2=it->first.second;
                        }else{
                            ID1=it->first.second, ID2=it->first.first;
                        }
                        if(affectedBPairs.find(make_pair(ID1,ID2))!=affectedBPairs.end()){//if found
                            if(newW<affectedBPairs[make_pair(ID1,ID2)].second){
                                affectedBPairs[make_pair(ID1,ID2)].second=newW;
                            }
                        }else{//if not found
                            affectedBPairs.insert({make_pair(ID1,ID2), make_pair(oldW,newW)});
                        }
                    }
                    for(auto it=affectedBPairs.begin();it!=affectedBPairs.end();++it){
                        testData.emplace_back(it->first,it->second);
                    }
                }
                else{//if there is no cross-leaf updates
                    int olddis, newdis;
                    for(auto it=affectedBPairs.begin();it!=affectedBPairs.end();++it){
                        ID1=it->first.first; ID2=it->first.second; newdis=it->second.second;
                        if(OverlayGraph[ID1].find(ID2) != OverlayGraph[ID1].end()){//if found
                            olddis=OverlayGraph[ID1][ID2];
                        }else{//if not found
                            cout<<"Not found edge e("<<ID1<<","<<ID2<<") in overlay graph!"<<endl; exit(1);
                        }
                        if(updateType==DECREASE){
                            if(newdis<olddis){
                                for(int i=0;i<Nodes[ID1].adjnodes.size();++i){
                                    if(Nodes[ID1].adjnodes[i]==ID2){
                                        cout<<"Exist original edge. "<<ID1<<" "<<ID2<<" "<<Nodes[ID1].adjweight[i]<<" "<<olddis<<" "<<newdis<<endl;
                                    }
                                }
                                testData.emplace_back(make_pair(ID1,ID2), make_pair(olddis,newdis));
                            }else if(newdis>olddis){
                                cout<<"Something wrong happens. "<<ID1<<"("<<Nodes[ID1].isborder<<") "<<ID2<<"("<<Nodes[ID2].isborder<<") : "<<newdis<<" "<<olddis<< endl;
                                exit(1);
                            }
                        }else if(updateType==INCREASE){
                            if(newdis>olddis){
                                for(int i=0;i<Nodes[ID1].adjnodes.size();++i){
                                    if(Nodes[ID1].adjnodes[i]==ID2){
                                        cout<<"Exist original edge. "<<ID1<<" "<<ID2<<" "<<Nodes[ID1].adjweight[i]<<" "<<olddis<<" "<<newdis<<endl;
                                    }
                                }
                                testData.emplace_back(make_pair(ID1,ID2), make_pair(olddis,newdis));
                            }else if(newdis<olddis){
                                cout<<"Something wrong happens. "<<ID1<<"("<<Nodes[ID1].isborder<<") "<<ID2<<"("<<Nodes[ID2].isborder<<") : "<<newdis<<" "<<olddis<< endl;
                                exit(1);
                            }
                        }

                    }
                }
//                for(auto it=affectedBPairs.begin();it!=affectedBPairs.end();++it){
//                    testData.emplace_back(it->first,it->second);
//                }
            }
        }
    }
    else{
        testData=uCrossLeafNode;
    }



    cout<<"affectedBPairs size: "<<affectedBPairs.size()<<" ; testData size: "<<testData.size()<< endl;
//    for(int i=0;i<testData.size();++i){
//        ID1=testData[i].first.first, ID2=testData[i].first.second, oldW=testData[i].second.first, newW=testData[i].second.second;
//        cout<<i<<": "<<ID1<<"("<<Nodes[ID1].isborder<<","<<Nodes[ID1].inleaf<<") "<<ID2<<"("<<Nodes[ID2].isborder<<","<<Nodes[ID2].inleaf<<") "<<oldW<<" "<<newW<<endl;
//    }
    /// Update overlay index
    if(updateType == DECREASE){
        BorderLabelUpdateDec(testData);
    }else if(updateType == INCREASE){
        BorderLabelUpdateInc(testData);
    }

}

void Gstartree::PartitionUpdate_Pre(int tn, map<pair<int,int>,pair<int,int>> & affectedBPairs, bool ifParallel){
    /// Construct Li
    vector<int> cands;
    vector<int> result;
    int s, t, nid, cid, weight, ID;


    // cands = leafnodes
    cands = GTree[tn].leafnodes;//for leaf node, the cands is the leafnodes
    int vNum = cands.size();
    // union borders = borders;
//    GTree[tn].union_borders = GTree[tn].borders;
    assert(GTree[tn].union_borders.size() == GTree[tn].borders.size());
    int bNum = GTree[tn].borders.size();

    if(ifParallel){/// multi-thread
        if(bNum>thread_num){
            int step=bNum/thread_num;
            boost::thread_group threadf;
            for(int i=0;i<thread_num;i++){
                pair<int,int> p;
                p.first=i*step;
                if(i==thread_num-1)
                    p.second=bNum;
                else
                    p.second=(i+1)*step;
                threadf.add_thread(new boost::thread(&Gstartree::boundaryVUpdate_Pre, this, p, tn, boost::ref(cands),boost::ref(affectedBPairs)));
            }
            threadf.join_all();
        }else{
            boost::thread_group threadf;
            for(int pid=0;pid<bNum;++pid) {
                threadf.add_thread(new boost::thread(&Gstartree::boundaryVUpdate_Pre, this, make_pair(pid,pid+1), tn, boost::ref(cands),boost::ref(affectedBPairs)));
            }
            threadf.join_all();
        }
    }
    else{///single thread
        // for each border, do min dis
        int cc = 0;
        for ( int k = 0; k < GTree[tn].union_borders.size(); k++ ){
            ID = GTree[tn].union_borders[k];
            //printf("DIJKSTRA...LEAF=%d BORDER=%d\n", tn, GTree[tn].union_borders[k] );
            result = dijkstra_candidate( ID, cands, Nodes );//distance vector from s to all borders
            //printf("DIJKSTRA...END\n");

            // save to map
            for ( int p = 0; p < result.size(); ++p ){
                if(result[p] != GTree[tn].mind[k*vNum + p]){
                    if(Nodes[cands[p]].isborder){//if it is boundary vertex
                        if(ID<cands[p]){
                            affectedBPairs.insert({make_pair(ID,cands[p]), make_pair(GTree[tn].mind[k*vNum + p],result[p])});
                        }else{
                            affectedBPairs.insert({make_pair(cands[p],ID), make_pair(GTree[tn].mind[k*vNum + p],result[p])});
                        }

                    }
                    GTree[tn].mind[k*vNum + p] = result[p];//update the distance matrix
                }

            }
        }
    }

//    int bid1,bid2;
//    for(int i=0;i<GTree[tn].borders.size();++i){
//        bid1=GTree[tn].borders[i];
//        for(int j=i+1;j<GTree[tn].borders.size();++j){
//            bid2=GTree[tn].borders[j];
//        }
//    }
}

void Gstartree::boundaryVUpdate_Pre(pair<int,int> p, int tn, vector<int> & cands, map<pair<int,int>,pair<int,int>> & affectedBPairs){
    int vNum = cands.size();
    int ID;

    for(int k=p.first;k<p.second;++k){
        ID = GTree[tn].union_borders[k];
        vector<int> result = dijkstra_candidate( ID, cands, Nodes );//distance vector from s to all borders

        // save to map
        for ( int p = 0; p < result.size(); ++p ){
            if(result[p] != GTree[tn].mind[k*vNum + p]){
                if(Nodes[cands[p]].isborder){//if it is boundary vertex
                    if(ID<cands[p]){
                        affectedBPairs.insert({make_pair(ID,cands[p]), make_pair(GTree[tn].mind[k*vNum + p],result[p])});
                    }
//                    else{
//                        affectedBPairs.insert({make_pair(cands[p],ID), make_pair(GTree[tn].mind[k*vNum + p],result[p])});
//                    }
                }
                GTree[tn].mind[k*vNum + p] = result[p];//update the distance matrix
            }
        }
    }

}

void Gstartree::BorderLabelUpdateDec(vector<pair<pair<int,int>,pair<int,int>>>& wBatch){
    /// Stage 1: update the shortcut graph
    vector<pair<pair<int,int>,pair<int,int>>> updateEdges;
    CHdecBat(wBatch,updateEdges);
//    cout<<"update shortcut pair number: "<<updateEdges.size()<<endl;

    /// Stage 2: update the border labels
    bool ifDebug=false;
    if(!updateEdges.empty()){
        int ID1,ID2,ID3,oldW,newW;
        unordered_set<int> LabelC;//label changed
        unordered_set<int> WeightC;//weight changed
        vector<unordered_map<int,int>> LabelCdis;//{s,(t,d)} maintain the fresh distance and avoid search in the adjacent list
        vector<unordered_map<int,int>> WeightCdis;//{s,(t,d)} maintain the fresh distance and avoid search in the adjacent list
        LabelCdis.assign(node_num,unordered_map<int,int>()); WeightCdis.assign(node_num,unordered_map<int,int>());
        for(auto it=updateEdges.begin();it!=updateEdges.end();++it){
            ID1=it->first.first, ID2=it->first.second, oldW=it->second.first, newW=it->second.second;
            if(Label[ID1].find(ID2)!=Label[ID1].end()){//if found
                WeightC.insert(ID1);
                WeightCdis[ID1].insert({ID2,oldW});
                if(Label[ID1][ID2]>newW){
                    LabelC.insert(ID1);
                    LabelCdis[ID1].insert({ID2,Label[ID1][ID2]});
                    Label[ID1][ID2]=newW;
                }
//                else{
//                    cout<<"Seems wrong for this label update. "<<ID1<<" "<<ID2<<" "<<Label[ID1][ID2]<<" "<<weight<<endl; exit(1);
//                }
            }else{
                cout<<"Wrong for this shortcut. "<<ID1<<"("<<NodeOrder[ID1]<<") "<<ID2<<"("<<NodeOrder[ID2]<<")"<<endl; exit(1);
            }
        }

        for(int i=vertexLevels.size()-2;i>=0;--i){
//        if(i%500==0){
//            cout<<"level "<<i<<endl;
//        }
            for(int j=0;j<vertexLevels[i].size();++j){
                ID1=vertexLevels[i][j];
//                if(ID1==196894){
//                    cout<<ID1<<endl;
//                    ifDebug=true;
//                }
                if(WeightC.find(ID1) == WeightC.end()){//if not found, which means that no edge weight changed for ID1
                    bool ifUpdate=false;
                    for(auto it=NeighborCon[ID1].begin();it!=NeighborCon[ID1].end();++it){
                        ID2=it->first; newW=it->second.first;
                        //only check the affected labels
                        if(LabelC.find(ID2)!=LabelC.end()){//if found
                            for(auto it2=LabelCdis[ID2].begin();it2!=LabelCdis[ID2].end();++it2){
                                ID3=it2->first;
                                if(Label[ID1][ID3]>Label[ID2][ID3]+newW){
                                    LabelCdis[ID1].insert({ID3,Label[ID1][ID3]});
                                    Label[ID1][ID3]=Label[ID2][ID3]+newW;
                                    ifUpdate=true;
                                }
                            }
                        }
                    }
                    if(ifUpdate){
                        LabelC.insert(ID1);
                    }
                }
                else{//if found, which means that there is edge weight changed for ID1
                    bool ifUpdate=false;
                    for(auto it=NeighborCon[ID1].begin();it!=NeighborCon[ID1].end();++it){
                        ID2=it->first; newW=it->second.first;
                        if(WeightCdis[ID1].find(ID2)!=WeightCdis[ID1].end()){//if found
                            for(auto it2=Label[ID2].begin();it2!=Label[ID2].end();++it2){
                                ID3=it2->first;
                                if(Label[ID1].find(ID3)==Label[ID1].end()){//if not found
                                    cout<<"Wrong! Label not exists. "<<ID1<<" "<<ID3<<endl; exit(1);
                                }else{//if found
                                    if(Label[ID1][ID3] > Label[ID2][ID3]+newW){
                                        LabelCdis[ID1].insert({ID3,Label[ID1][ID3]});
                                        Label[ID1][ID3] = Label[ID2][ID3]+newW;
                                        ifUpdate=true;
                                    }
                                }
                            }
                        }
                        else{
                            //only check the affected labels
                            if(LabelC.find(ID2)!=LabelC.end()){//if found
                                for(auto it2=LabelCdis[ID2].begin();it2!=LabelCdis[ID2].end();++it2){
                                    ID3=it2->first;
                                    if(Label[ID1][ID3]>Label[ID2][ID3]+newW){
                                        LabelCdis[ID1].insert({ID3,Label[ID1][ID3]});
                                        Label[ID1][ID3]=Label[ID2][ID3]+newW;
                                        ifUpdate=true;
                                    }
                                }
                            }
                        }
                    }
                    if(ifUpdate){
                        LabelC.insert(ID1);
                    }
                }

            }
        }
    }

}

int Gstartree::DisQueryVally(int ID1, int ID2){
    int neiID;
    int neiDis,d=INF;
    for(int i=0;i<NeighborCon[ID1].size();i++){
        neiID=NeighborCon[ID1][i].first;
        neiDis=NeighborCon[ID1][i].second.first;
        if(NodeOrder[neiID]<=NodeOrder[ID2] && Label[neiID].find(ID2)!=Label[neiID].end()){//if found
            if(neiDis+Label[neiID][ID2]<d){
                d=neiDis+Label[neiID][ID2];
            }
        }
    }
    return d;
}

void Gstartree::BorderLabelUpdateInc(vector<pair<pair<int,int>,pair<int,int>>>& wBatch){
    /// Stage 1: update the shortcut graph
    vector<pair<pair<int,int>,pair<int,int>>> updateEdges;
    CHincBatMT(wBatch,updateEdges);
    cout<<"update shortcut pair number: "<<updateEdges.size()<<endl;

    /// Stage 2: update the border labels
    bool ifDebug=false;
    if(!updateEdges.empty()){
        int ID1,ID2,ID3,oldW,newW,olddis,weight,newdis;
        unordered_set<int> LabelC;//label changed
        unordered_set<int> WeightC;//weight changed
        vector<unordered_map<int,int>> LabelCdis;//{s,(t,d)} maintain the fresh distance and avoid search in the adjacent list
        vector<unordered_map<int,int>> WeightCdis;//{s,(t,d)} maintain the fresh distance and avoid search in the adjacent list
        LabelCdis.assign(node_num,unordered_map<int,int>()); WeightCdis.assign(node_num,unordered_map<int,int>());
        for(auto it=updateEdges.begin();it!=updateEdges.end();++it){
            ID1=it->first.first, ID2=it->first.second, oldW=it->second.first, newW=it->second.second;
            if(Label[ID1].find(ID2)!=Label[ID1].end()){//if found
                WeightC.insert(ID1);
                WeightCdis[ID1].insert({ID2,oldW});
                if(Label[ID1][ID2]==oldW){
                    newdis=DisQueryVally(ID1,ID2);//may be incorrect
                    if(Label[ID1][ID2]<newdis){
                        LabelC.insert(ID1);
                        LabelCdis[ID1].insert({ID2,Label[ID1][ID2]});
//                    Label[ID1][ID2]=newW;
                        Label[ID1][ID2]=newdis;
                    }
                    else{
                        assert(Label[ID1][ID2]==newdis);
                    }

                }
//                else{
//                    cout<<"Seems wrong for this label update. "<<ID1<<" "<<ID2<<" "<<Label[ID1][ID2]<<" "<<weight<<endl; exit(1);
//                }
            }else{
                cout<<"Wrong for this shortcut. "<<ID1<<"("<<NodeOrder[ID1]<<") "<<ID2<<"("<<NodeOrder[ID2]<<")"<<endl; exit(1);
            }
        }

        for(int i=vertexLevels.size()-2;i>=0;--i){
//        if(i%500==0){
//            cout<<"level "<<i<<endl;
//        }
            for(int j=0;j<vertexLevels[i].size();++j){
                ID1=vertexLevels[i][j];
//                if(ID1==196894){
//                    cout<<ID1<<endl;
//                    ifDebug=true;
//                }
                if(WeightC.find(ID1) == WeightC.end()){//if not found, which means that no edge weight changed for ID1
                    bool ifUpdate=false;
                    for(auto it=NeighborCon[ID1].begin();it!=NeighborCon[ID1].end();++it){
                        ID2=it->first; weight=it->second.first;
                        //only check the affected labels
                        if(LabelC.find(ID2)!=LabelC.end()){//if found
                            for(auto it2=LabelCdis[ID2].begin();it2!=LabelCdis[ID2].end();++it2){
                                ID3=it2->first; olddis=it2->second;
                                if(Label[ID1][ID3]==olddis+weight){
                                    newdis=DisQueryVally(ID1,ID3);//may be incorrect
                                    if(Label[ID1][ID3]<newdis){
                                        LabelCdis[ID1].insert({ID3,Label[ID1][ID3]});
                                        Label[ID1][ID3]=newdis;
                                        ifUpdate=true;
                                    }
                                }
                            }
                        }
                    }
                    if(ifUpdate){
                        LabelC.insert(ID1);
                    }
                }
                else{//if found, which means that there is edge weight changed for ID1
                    bool ifUpdate=false;
                    for(auto it=NeighborCon[ID1].begin();it!=NeighborCon[ID1].end();++it){
                        ID2=it->first; newW=it->second.first;
                        if(WeightCdis[ID1].find(ID2)!=WeightCdis[ID1].end()){//if found, edge weight from ID1 to ID2 has changed
                            oldW=WeightCdis[ID1][ID2];
                            for(auto it2=Label[ID2].begin();it2!=Label[ID2].end();++it2){
                                ID3=it2->first;
                                if(Label[ID1].find(ID3)==Label[ID1].end()){//if not found
                                    cout<<"Wrong! Label not exists. "<<ID1<<" "<<ID3<<endl; exit(1);
                                }else{//if found
                                    if(LabelCdis[ID2].find(ID3)!=LabelCdis[ID2].end()){//if found
                                        olddis=LabelCdis[ID2][ID3];
                                    }else{
                                        olddis=Label[ID2][ID3];
                                    }
                                    if(Label[ID1][ID3] == olddis+oldW){
                                        newdis=DisQueryVally(ID1,ID3);//may be incorrect
                                        if(Label[ID1][ID3]<newdis){
                                            LabelCdis[ID1].insert({ID3,Label[ID1][ID3]});
                                            Label[ID1][ID3]=newdis;
                                            ifUpdate=true;
                                        }
                                    }
                                }
                            }


                        }
                        else{//if edge (ID1,ID2) does not change
                            //only check the affected labels
                            if(LabelC.find(ID2)!=LabelC.end()){//if found
                                for(auto it2=LabelCdis[ID2].begin();it2!=LabelCdis[ID2].end();++it2){
                                    ID3=it2->first; olddis=LabelCdis[ID2][ID3];
                                    if(Label[ID1][ID3]==olddis+newW){
                                        newdis=DisQueryVally(ID1,ID3);//may be incorrect
                                        if(Label[ID1][ID3]<newdis){
                                            LabelCdis[ID1].insert({ID3,Label[ID1][ID3]});
                                            Label[ID1][ID3]=newdis;
                                            ifUpdate=true;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if(ifUpdate){
                        LabelC.insert(ID1);
                    }
                }

            }
        }
    }
}

void Gstartree::CHdecBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& updateEdges){
    //maintain the index caused by the weight change
    //NodeOrders.clear();
    set<OrderComp> OC;
    map<pair<int,int>,int> OCdis;//{(s,t),d} maintain the fresh distance and avoid search in the adjacent list
    //OC.clear(); OCdis.clear();

    int a,b,newW;//the weight of (a,b) decrease to newW
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first;
        b=wBatch[k].first.second;
        newW=wBatch[k].second.second;

        //modify the information in original graph
        if(OverlayGraph[a].find(b)!=OverlayGraph[a].end()){
            OverlayGraph[a][b]=newW;
        }else{
            cout<<"Wrong for Neighbors!"<<endl; exit(1);
        }
        if(OverlayGraph[b].find(a)!=OverlayGraph[b].end()){
            OverlayGraph[b][a]=newW;
        }else{
            cout<<"Wrong for Neighbors!"<<endl; exit(1);
        }

        if(NodeOrder[a]<NodeOrder[b]){
            for(int i=0;i<NeighborCon[a].size();i++){
                if(NeighborCon[a][i].first==b){
                    if(NeighborCon[a][i].second.first>newW){
                        //cout<<OutNeighborCon[a][i].second.first<<"..........."<<newW<<endl;
                        updateEdges.emplace_back(make_pair(a,b),make_pair(NeighborCon[a][i].second.first,newW));
                        NeighborCon[a][i].second.first=newW;
                        NeighborCon[a][i].second.second=1;

                        OCdis[make_pair(a,b)]=newW;
//                        cout<<a<<"("<<NodeOrder_[a]<<") "<<b<<"("<<NodeOrder_[b]<<")"<<endl;
                        OC.insert(OrderComp(a,b));
                    }else if(NeighborCon[a][i].second.first==newW)
                        NeighborCon[a][i].second.second+=1;
                    break;
                }
            }
        }else{
            for(int i=0;i<NeighborCon[b].size();i++){
                if(NeighborCon[b][i].first==a){
                    if(NeighborCon[b][i].second.first>newW){
                        updateEdges.emplace_back(make_pair(b,a),make_pair(NeighborCon[b][i].second.first,newW));
                        NeighborCon[b][i].second.first=newW;
                        NeighborCon[b][i].second.second=1;

                        OCdis[make_pair(b,a)]=newW;
                        OC.insert(OrderComp(b,a));
                    }else if(NeighborCon[b][i].second.first==newW)
                        NeighborCon[b][i].second.second+=1;
                    break;
                }
            }
        }
    }


    while(!OC.empty()){
        int s=(*OC.begin()).x; int t=(*OC.begin()).y;
        int wt;
        OC.erase(OC.begin());
        wt=OCdis[make_pair(s,t)];
        map<int,int> InM2t; //InM2t.clear();
        vector<pair<int,int>> InMLower; //InMLower.clear();
        for(int i=0;i<NeighborCon[s].size();i++){
            if(NodeOrder[NeighborCon[s][i].first]>NodeOrder[t])
                InM2t.insert(make_pair(NeighborCon[s][i].first,NeighborCon[s][i].second.first));
            else if(NodeOrder[NeighborCon[s][i].first]<NodeOrder[t])
                InMLower.push_back(make_pair(NeighborCon[s][i].first,NeighborCon[s][i].second.first));
        }
        int inID,inW,inWt;
        for(int i=0;i<NeighborCon[t].size();i++){
            inID=NeighborCon[t][i].first;
            if(InM2t.find(inID)!=InM2t.end()){
                inW=InM2t[inID];
                inWt=NeighborCon[t][i].second.first;
                if(inWt>inW+wt){
                    updateEdges.emplace_back(make_pair(t,inID),make_pair(NeighborCon[t][i].second.first,inW+wt));
                    NeighborCon[t][i].second.first=inW+wt;
                    NeighborCon[t][i].second.second=1;

                    OCdis[make_pair(t,inID)]=inW+wt;
                    OrderComp oc={t,inID};
                    OC.insert(oc);
                }else if(inWt==inW+wt){
                    NeighborCon[t][i].second.second+=1;
                }
            }
        }

        for(int i=0;i<InMLower.size();i++){
            inID=InMLower[i].first; inW=InMLower[i].second;
            for(int j=0;j<NeighborCon[inID].size();j++){
                if(NeighborCon[inID][j].first==t){
                    inWt=NeighborCon[inID][j].second.first;
                    if(inWt>inW+wt){
                        updateEdges.emplace_back(make_pair(inID,t),make_pair(NeighborCon[inID][j].second.first,inW+wt));
                        NeighborCon[inID][j].second.first=inW+wt;
                        NeighborCon[inID][j].second.second=1;

                        OCdis[make_pair(inID,t)]=inW+wt;
                        OrderComp oc={inID,t};
                        OC.insert(oc);
                    }else if(inWt==inW+wt)
                        NeighborCon[inID][j].second.second+=1;
                    break;
                }
            }
        }
    }//finish change index
}

void Gstartree::CHincBatMT(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& updateEdges){
    //NodeOrders.clear();
    set<OrderComp> OC; //OC.clear();
    map<pair<int,int>,int> OCdis;//{(s,t),d} maintain the old distance before refreshed and avoid search in the adjacent list
    //OCdis.clear();

    for(int wb=0;wb<wBatch.size();wb++){
        int a=wBatch[wb].first.first;
        int b=wBatch[wb].first.second;
        int oldW=wBatch[wb].second.first;
        int newW=wBatch[wb].second.second;

        //modify the original graph information
        if(OverlayGraph[a].find(b)!=OverlayGraph[a].end()){
            OverlayGraph[a][b]=newW;
        }else{
            cout<<"Wrong for Neighbors!"<<endl; exit(1);
        }
        if(OverlayGraph[b].find(a)!=OverlayGraph[b].end()){
            OverlayGraph[b][a]=newW;
        }else{
            cout<<"Wrong for Neighbors!"<<endl; exit(1);
        }


        if(NodeOrder[a]<NodeOrder[b]){
            for(int i=0;i<NeighborCon[a].size();i++){
                if(NeighborCon[a][i].first==b){
                    if(NeighborCon[a][i].second.first==oldW){
                        NeighborCon[a][i].second.second-=1;
                        if(NeighborCon[a][i].second.second<1){
                            OrderComp oc={a,b};
                            OC.insert(oc);
                            OCdis[make_pair(a,b)]=oldW;
                        }
                    }
                    break;
                }
            }
        }else{
            for(int i=0;i<NeighborCon[b].size();i++){
                if(NeighborCon[b][i].first==a){
                    if(NeighborCon[b][i].second.first==oldW){
                        NeighborCon[b][i].second.second-=1;
                        if(NeighborCon[b][i].second.second<1){
                            OrderComp oc={b,a};
                            OC.insert(oc);
                            OCdis[make_pair(b,a)]=oldW;
                        }
                    }
                    break;
                }
            }
        }
    }

    while(!OC.empty()){
        int s=(*OC.begin()).x; int t=(*OC.begin()).y;
        int wt;
        OC.erase(OC.begin());
        wt=OCdis[make_pair(s,t)];//distance of s--->t before change
        int inID,inW;
        map<int,int> HigherIn; vector<pair<int,int>> LowerIn;
        //HigherIn.clear(); LowerIn.clear();
        //the shortcuts infected by s-->t
        for(int i=0;i<NeighborCon[s].size();i++){
            inID=NeighborCon[s][i].first;
            inW=NeighborCon[s][i].second.first;
            if(NodeOrder[inID]<NodeOrder[t]){
                LowerIn.emplace_back(inID,inW);
            }else if(NodeOrder[inID]>NodeOrder[t]){
                HigherIn.insert(make_pair(inID,inW));
            }
        }
        for(int i=0;i<NeighborCon[t].size();i++){
            inID=NeighborCon[t][i].first;
            if(HigherIn.find(inID)!=HigherIn.end()){
                inW=HigherIn[inID];
                if(NeighborCon[t][i].second.first==wt+inW){
                    NeighborCon[t][i].second.second-=1;
                    if(NeighborCon[t][i].second.second<1){
                        OrderComp oc={t,inID};
                        OC.insert(oc);
                        OCdis[make_pair(t,inID)]=wt+inW;
                    }
                }
            }
        }
        for(int i=0;i<LowerIn.size();i++){
            inID=LowerIn[i].first; inW=LowerIn[i].second;
            for(int j=0;j<NeighborCon[inID].size();j++){
                if(NeighborCon[inID][j].first==t){
                    if(NeighborCon[inID][j].second.first==inW+wt){
                        NeighborCon[inID][j].second.second-=1;
                        if(NeighborCon[inID][j].second.second<1){
                            OrderComp oc={inID,t};
                            OC.insert(oc);
                            OCdis[make_pair(inID,t)]=wt+inW;
                        }
                    }
                    break;
                }
            }
        }

        //get the new weight value of s-->t
        wt=INF; int countwt=0;
//        for(int i=0;i<Neighbor[s].size();i++){
//            if(Neighbor[s][i].first==t){
//                wt=Neighbor[s][i].second;//the weight value in the original graph
//                countwt=1;
//                break;
//            }
//        }
        for (auto it = OverlayGraph[s].begin(); it != OverlayGraph[s].end(); ++it) {///Neighbor of overlay graph
            if (it->first == t) {
                wt = it->second;//the weight value in the original graph
                countwt = 1;
                break;
            }
        }
        int ssw,wtt,wid;
        vector<int> Wnodes; //Wnodes.clear();
        if(s<t){
            //Wnodes=SCconNodes[make_pair(s,t)]; //cout<<"wid num "<<Wnodes.size()<<endl;
            Wnodes=SCconNodesMT[s][t];
        }else{
            //Wnodes=SCconNodes[make_pair(t,s)];
            Wnodes=SCconNodesMT[t][s];
        }

        for(int i=0;i<Wnodes.size();i++){
            wid=Wnodes[i];
            for(int j=0;j<NeighborCon[wid].size();j++){
                if(NeighborCon[wid][j].first==s){
                    ssw=NeighborCon[wid][j].second.first;
                }
                if(NeighborCon[wid][j].first==t){
                    wtt=NeighborCon[wid][j].second.first;
                }
            }

            if(ssw+wtt<wt){
                wt=ssw+wtt;
                countwt=1;
            }else if(ssw+wtt==wt){
                countwt+=1;
            }
        }

        //refresh the weight value of s--t in the index
        for(int i=0;i<NeighborCon[s].size();i++){
            if(NeighborCon[s][i].first==t){
//                cout<<"Refresh shortcut: "<<s<<" "<<t<<" "<<NeighborCon[s][i].second.first<<" "<<wt<<endl;
                updateEdges.emplace_back(make_pair(s,t), make_pair(NeighborCon[s][i].second.first,wt));
                NeighborCon[s][i].second.first=wt;
                NeighborCon[s][i].second.second=countwt;
                break;
            }
        }
    }
}

void Gstartree::TGTreeIndexUpdate(vector<pair<pair<int,int>,pair<int,int> > > & updates, bool ifParallel, int updateType){//update(<ID1,ID2>,oldW) vector<pair<pair<int,int>,int> > & updates,
//    double ave_time = 0;
    int ID1,ID2,oldW,newW;
    vector<pair<pair<int,int>,pair<int,int> > > uWithinLeafNode;// edge update within the same leaf node (<ID1,ID2>,<oldW,newW>)
    vector<pair<pair<int,int>,pair<int,int> > > uCrossLeafNode;// edge/shortcut update among borders (<ID1,ID2>,<oldW,newW>)
    bool flag_borderUpdate = false;//whether there is border update of leaf node
//    Timer tt;
//    tt.start();

    for(int i=0;i<updates.size();++i){//for each edge update
        ID1 = updates[i].first.first;
        ID2 = updates[i].first.second;
        oldW = updates[i].second.first;
        newW = updates[i].second.second;

//        cout<<ID1 << " "<<ID2<<" ("<<oldW<<"->"<<newW<<") "<<endl;
        //update edge weight
        for(int j=0;j<Nodes[ID1].adjnodes.size();j++){
            if(Nodes[ID1].adjnodes[j]==ID2){
                assert(Nodes[ID1].adjweight[j] == oldW);
                Nodes[ID1].adjweight[j]=newW;
                break;
            }
        }
        for(int j=0;j<Nodes[ID2].adjnodes.size();j++){
            if(Nodes[ID2].adjnodes[j]==ID1){
                assert(Nodes[ID2].adjweight[j] == oldW);
                Nodes[ID2].adjweight[j]=newW;
                break;
            }
        }
        //identify the edge type
        if(Nodes[ID1].inleaf == Nodes[ID2].inleaf){//if it is an update within the leaf node
            uWithinLeafNode.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
        }else{//if it is an update among the leaf nodes
            uCrossLeafNode.emplace_back(make_pair(ID1,ID2),make_pair(oldW,newW));
            flag_borderUpdate = true;
        }
    }

    /// update processing
    int lnID;
//        int ID1_pos,ID2_pos;
//    vector<Node> graph;//temp graph
//    graph = Nodes;//original graph

    /// Update the leaf nodes
    Timer tt1;
    tt1.start();
    // for update within the same leaf node
    vector<pair<pair<int,int>,pair<int,int> > > testData;
    map<pair<int,int>,pair<int,int>> affectedBPairs;
    if(!uWithinLeafNode.empty()){//if same-leaf node update
        for(auto it=uWithinLeafNode.begin();it!=uWithinLeafNode.end();++it){//deal with each update within leaf node
            ID1 = it->first.first; ID2 = it->first.second; oldW = it->second.first; newW = it->second.second;
            assert(Nodes[ID1].inleaf == Nodes[ID2].inleaf);
            int PID = Nodes[ID1].inleaf;

//            PartitionUpdate_No(PID,affectedBPairs,false);
            PartitionUpdate_No(PID,affectedBPairs,true);
            //update overlay graph index
            if(!affectedBPairs.empty()){
//                cout<<" affectedBPairs size: "<<affectedBPairs.size()<<endl;
                if(!uCrossLeafNode.empty()){//if there is cross-leaf updates
                    cout<<"Multiple update: in-leaf node and cross-leaf node."<<endl;
                    for(auto it=uCrossLeafNode.begin();it!=uCrossLeafNode.end();++it){
                        oldW=it->second.first; newW=it->second.second;
                        if(it->first.first<it->first.second){
                            ID1=it->first.first, ID2=it->first.second;
                        }else{
                            ID1=it->first.second, ID2=it->first.first;
                        }
                        if(affectedBPairs.find(make_pair(ID1,ID2))!=affectedBPairs.end()){//if found
                            if(newW<affectedBPairs[make_pair(ID1,ID2)].second){
                                affectedBPairs[make_pair(ID1,ID2)].second=newW;
                            }
                        }else{//if not found
                            affectedBPairs.insert({make_pair(ID1,ID2), make_pair(oldW,newW)});
                        }
                    }
                    for(auto it=affectedBPairs.begin();it!=affectedBPairs.end();++it){
                        testData.emplace_back(it->first,it->second);
                    }
                }
                else{//if there is no cross-leaf updates
                    int olddis, newdis;
                    for(auto it=affectedBPairs.begin();it!=affectedBPairs.end();++it){
                        ID1=it->first.first; ID2=it->first.second; newdis=it->second.second;
                        if(OverlayGraph[ID1].find(ID2) != OverlayGraph[ID1].end()){//if found
                            olddis=OverlayGraph[ID1][ID2];
                        }else{//if not found
                            cout<<"Not found edge e("<<ID1<<","<<ID2<<") in overlay graph!"<<endl; exit(1);
                        }
                        if(updateType==DECREASE){
                            if(newdis<olddis){
                                for(int i=0;i<Nodes[ID1].adjnodes.size();++i){
                                    if(Nodes[ID1].adjnodes[i]==ID2){
                                        cout<<"Exist original edge. "<<ID1<<" "<<ID2<<" "<<Nodes[ID1].adjweight[i]<<" "<<olddis<<" "<<newdis<<endl;
                                    }
                                }
                                testData.emplace_back(make_pair(ID1,ID2), make_pair(olddis,newdis));
                            }else if(newdis>olddis){
                                cout<<"Something wrong happens. "<<ID1<<"("<<Nodes[ID1].isborder<<") "<<ID2<<"("<<Nodes[ID2].isborder<<") : "<<newdis<<" "<<olddis<< endl;
                                exit(1);
                            }
                        }else if(updateType==INCREASE){
                            if(newdis>olddis){
                                for(int i=0;i<Nodes[ID1].adjnodes.size();++i){
                                    if(Nodes[ID1].adjnodes[i]==ID2){
                                        cout<<"Exist original edge. "<<ID1<<" "<<ID2<<" "<<Nodes[ID1].adjweight[i]<<" "<<olddis<<" "<<newdis<<endl;
                                    }
                                }
                                testData.emplace_back(make_pair(ID1,ID2), make_pair(olddis,newdis));
                            }else if(newdis<olddis){
                                cout<<"Something wrong happens. "<<ID1<<"("<<Nodes[ID1].isborder<<") "<<ID2<<"("<<Nodes[ID2].isborder<<") : "<<newdis<<" "<<olddis<< endl;
                                exit(1);
                            }
                        }

                    }
                }
//                for(auto it=affectedBPairs.begin();it!=affectedBPairs.end();++it){
//                    testData.emplace_back(it->first,it->second);
//                }
            }
        }
    }
    else{
        testData=uCrossLeafNode;
    }



//    cout<<"affectedBPairs size: "<<affectedBPairs.size()<<" ; testData size: "<<testData.size()<< endl;
//    for(int i=0;i<testData.size();++i){
//        ID1=testData[i].first.first, ID2=testData[i].first.second, oldW=testData[i].second.first, newW=testData[i].second.second;
//        cout<<i<<": "<<ID1<<"("<<Nodes[ID1].isborder<<","<<Nodes[ID1].inleaf<<") "<<ID2<<"("<<Nodes[ID2].isborder<<","<<Nodes[ID2].inleaf<<") "<<oldW<<" "<<newW<<endl;
//    }
    /// Update overlay index
    if(updateType == DECREASE){
        H2HdecBat(testData);
    }else if(updateType == INCREASE){
        H2HincBatMT(testData);
    }

}


void Gstartree::PartitionUpdate_No(int tn, map<pair<int,int>,pair<int,int>> & affectedBPairs, bool ifParallel){
    /// Construct Li
    vector<int> cands;
    vector<int> result;
    int s, t, nid, cid, weight, ID;


    // cands = leafnodes
    cands = GTree[tn].leafnodes;//for leaf node, the cands is the leafnodes
    int vNum = cands.size();
    // union borders = borders;
//    GTree[tn].union_borders = GTree[tn].borders;
    assert(GTree[tn].union_borders.size() == GTree[tn].borders.size());
    int bNum = GTree[tn].borders.size();

    if(ifParallel){/// multi-thread
        if(bNum>thread_num){
            int step=bNum/thread_num;
            boost::thread_group threadf;
            for(int i=0;i<thread_num;i++){
                pair<int,int> p;
                p.first=i*step;
                if(i==thread_num-1)
                    p.second=bNum;
                else
                    p.second=(i+1)*step;
                threadf.add_thread(new boost::thread(&Gstartree::boundaryVUpdate_No, this, p, tn, boost::ref(cands),boost::ref(affectedBPairs)));
            }
            threadf.join_all();
        }else{
            boost::thread_group threadf;
            for(int pid=0;pid<bNum;++pid) {
                threadf.add_thread(new boost::thread(&Gstartree::boundaryVUpdate_No, this, make_pair(pid,pid+1), tn, boost::ref(cands),boost::ref(affectedBPairs)));
            }
            threadf.join_all();
        }
    }
    else{///single thread
        // for each border, do min dis
        int cc = 0;
        for ( int k = 0; k < GTree[tn].union_borders.size(); k++ ){
            ID = GTree[tn].union_borders[k];
            //printf("DIJKSTRA...LEAF=%d BORDER=%d\n", tn, GTree[tn].union_borders[k] );
//                result = dijkstra_candidate( GTree[tn].union_borders[k], cands, graph );//distance vector from s to all borders
//            result = dijkstra_candidate_No( ID, cands, PartiGraph[tn] );//distance vector from s to all borders
            result = dijkstra_candidate_No( ID, cands, Nodes );//distance vector from s to all borders
            //printf("DIJKSTRA...END\n");

            // save to map
            for ( int p = 0; p < result.size(); ++p ){
                if(result[p] != GTree[tn].mind[k*vNum + p]){
                    if(Nodes[cands[p]].isborder){//if it is boundary vertex
                        if(ID<cands[p]){
                            affectedBPairs.insert({make_pair(ID,cands[p]), make_pair(GTree[tn].mind[k*vNum + p],result[p])});
                        }else{
                            affectedBPairs.insert({make_pair(cands[p],ID), make_pair(GTree[tn].mind[k*vNum + p],result[p])});
                        }

                    }
                    GTree[tn].mind[k*vNum + p] = result[p];//update the distance matrix
                }

            }
        }
    }

//    int bid1,bid2;
//    for(int i=0;i<GTree[tn].borders.size();++i){
//        bid1=GTree[tn].borders[i];
//        for(int j=i+1;j<GTree[tn].borders.size();++j){
//            bid2=GTree[tn].borders[j];
//        }
//    }
}

void Gstartree::boundaryVUpdate_No(pair<int,int> p, int tn, vector<int> & cands, map<pair<int,int>,pair<int,int>> & affectedBPairs){
    int vNum = cands.size();
    int ID;

    for(int k=p.first;k<p.second;++k){
        ID = GTree[tn].union_borders[k];
//        vector<int> result = dijkstra_candidate_No( ID, cands, PartiGraph[tn] );//distance vector from s to all borders
        vector<int> result = dijkstra_candidate_No( ID, cands, Nodes );//distance vector from s to all borders

        // save to map
        for ( int p = 0; p < result.size(); ++p ){
            if(result[p] != GTree[tn].mind[k*vNum + p]){
                if(Nodes[cands[p]].isborder){//if it is boundary vertex
                    if(ID<cands[p]){
                        affectedBPairs.insert({make_pair(ID,cands[p]), make_pair(GTree[tn].mind[k*vNum + p],result[p])});
                    }
//                    else{
//                        affectedBPairs.insert({make_pair(cands[p],ID), make_pair(GTree[tn].mind[k*vNum + p],result[p])});
//                    }
                }
                GTree[tn].mind[k*vNum + p] = result[p];//update the distance matrix
            }
        }
    }

}


void Gstartree::H2HdecBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch) {
    map<int, int> checkedDis;

    for (int i = 0; i < Tree.size(); i++) {
        Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
    }

//NodeOrderss.clear();
    vector<set<int>> SCre; //SCre.clear();
    set<int> ss; //ss.clear();
    SCre.assign(node_num, ss);//{vertexID, set<int>}
    set<OrderCompMin> OC; //OC.clear();//vertexID in decreasing node order

    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed

    int a, b, oldW, newW, lid, hid;
    for (int k = 0; k < wBatch.size(); k++) {
        a = wBatch[k].first.first;
        b = wBatch[k].first.second;
        oldW = wBatch[k].second.first;
        newW = wBatch[k].second.second;
        if (NodeOrder[a] < NodeOrder[b]) {
            lid = a;
            hid = b;
        } else {
            lid = b;
            hid = a;
        }

        if(OverlayGraph[a].find(b)!=OverlayGraph[a].end()){
            OverlayGraph[a][b]=newW;
        }else{
            cout<<"Wrong for Neighbors!"<<endl; exit(1);
        }
        if(OverlayGraph[b].find(a)!=OverlayGraph[b].end()){
            OverlayGraph[b][a]=newW;
        }else{
            cout<<"Wrong for Neighbors!"<<endl; exit(1);
        }


        for (int i = 0; i < Tree[rank[lid]].vert.size(); i++) {
            if (Tree[rank[lid]].vert[i].first == hid) {
                if (Tree[rank[lid]].vert[i].second.first > newW) {
                    Tree[rank[lid]].vert[i].second.first = newW;
                    Tree[rank[lid]].vert[i].second.second = 1;
                    SCre[lid].insert(hid);
                    OC.insert(OrderCompMin(lid));
                } else if (Tree[rank[lid]].vert[i].second.first == newW) {
                    Tree[rank[lid]].vert[i].second.second += 1;
                }
                break;
            }
        }

    }

    vector<int> ProBeginVertexSet; //ProBeginVertexSet.clear();
    vector<int> ProBeginVertexSetNew;
    int ProBeginVertexID;
    int ProID;
//processing the stars
    while (!OC.empty()) {
        ProID = (*OC.begin()).x;
        OC.erase(OC.begin());
        vector<pair<int, pair<int, int>>> Vert = Tree[rank[ProID]].vert;
        bool ProIDdisCha = false;//to see if the distance labeling of proID change or not
        for (auto it = SCre[ProID].begin(); it != SCre[ProID].end(); it++) {
            int Cid = *it;
            int Cw;
            int cidH = Tree[rank[Cid]].height - 1;

            map<int, int> Hnei; //Hnei.clear();
            vector<pair<int, int>> Lnei; //Lnei.clear();
            for (int j = 0; j < Vert.size(); j++) {
                if (NodeOrder[Vert[j].first] > NodeOrder[Cid]) {
                    Hnei[Vert[j].first] = Vert[j].second.first;
                } else if (NodeOrder[Vert[j].first] < NodeOrder[Cid]) {
                    Lnei.push_back(make_pair(Vert[j].first, Vert[j].second.first));
                } else {
                    Cw = Vert[j].second.first;
                }
            }

            if (Tree[rank[ProID]].dis[cidH] > Cw) {
                Tree[rank[ProID]].dis[cidH] = Cw;
                Tree[rank[ProID]].FN[cidH] = true;
                ProIDdisCha = true;
                Tree[rank[ProID]].DisRe.insert(Cid);
            } else if (Tree[rank[ProID]].dis[cidH] == Cw) {
                Tree[rank[ProID]].FN[cidH] = true;
            }

            int hid, hidHeight, lid, lidHeight, wsum;
            for (int j = 0; j < Tree[rank[Cid]].vert.size(); j++) {
                hid = Tree[rank[Cid]].vert[j].first;
                hidHeight = Tree[rank[hid]].height - 1;
                if (Hnei.find(hid) != Hnei.end()) {
                    wsum = Cw + Hnei[hid];
                    if (wsum < Tree[rank[Cid]].vert[j].second.first) {
                        Tree[rank[Cid]].vert[j].second.first = wsum;
                        Tree[rank[Cid]].vert[j].second.second = 1;
                        SCre[Cid].insert(hid);
                        OC.insert(OrderCompMin(Cid));
                    } else if (wsum == Tree[rank[Cid]].vert[j].second.first) {
                        Tree[rank[Cid]].vert[j].second.second += 1;
                    }

                }
            }
            for (int j = 0; j < Lnei.size(); j++) {
                lid = Lnei[j].first;
                lidHeight = Tree[rank[lid]].height - 1;
                for (int k = 0; k < Tree[rank[lid]].vert.size(); k++) {
                    if (Tree[rank[lid]].vert[k].first == Cid) {
                        wsum = Cw + Lnei[j].second;
                        if (Tree[rank[lid]].vert[k].second.first > wsum) {
                            Tree[rank[lid]].vert[k].second.first = wsum;
                            Tree[rank[lid]].vert[k].second.second = 1;
                            SCre[lid].insert(Cid);
                            OC.insert(OrderCompMin(lid));
                        } else if (Tree[rank[lid]].vert[k].second.first == wsum) {
                            Tree[rank[lid]].vert[k].second.second += 1;
                        }

                        break;
                    }
                }
            }
        }

        if (ProIDdisCha) {//if the distance labeling is dectected changed
            vertexIDChL.insert(ProID);
            ProBeginVertexSetNew.clear();
            ProBeginVertexSetNew.reserve(ProBeginVertexSet.size() + 1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew = rank[ProID], r;
            for (int i = 0; i < ProBeginVertexSet.size(); i++) {
                r = rank[ProBeginVertexSet[i]];
                if (LCAQuery(rnew, r) != rnew) {
                    ProBeginVertexSetNew.push_back(ProBeginVertexSet[i]);
                }
            }
            ProBeginVertexSet = ProBeginVertexSetNew;
        }
    }

//cout<<"Finish bottom-up refresh"<<endl;
    for (int i = 0; i < ProBeginVertexSet.size(); i++) {
        ProBeginVertexID = ProBeginVertexSet[i];
        vector<int> linee; //linee.clear();
        linee.reserve(heightMax);
        int pachidd = Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
        while (Tree[rank[pachidd]].height > 1) {
            linee.insert(linee.begin(), pachidd);
            pachidd = Tree[Tree[rank[pachidd]].pa].uniqueVertex;
        }
        linee.insert(linee.begin(), pachidd);
        EachNodeProBDis5(rank[ProBeginVertexID], linee, vertexIDChL, checkedDis);
    }
//return checkedDis.size();
}

void Gstartree::EachNodeProBDis5(int child,vector<int>& line,set<int>& vertexIDChL, map<int,int>& checkedDis){
    bool ProIDdisCha=false;

    if(Tree[child].DisRe.size()!=0){
        for(int k=0;k<Tree[child].vert.size();k++){
            int b=Tree[child].vert[k].first, bH=Tree[rank[b]].height-1,vbW=Tree[child].vert[k].second.first;
            if(Tree[child].FN[bH]){
                if(Tree[child].DisRe.find(b)!=Tree[child].DisRe.end()){//all ancestor check
                    for(int i=0;i<bH;i++){
                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
                            Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
                            Tree[child].FN[i]=false;
                            ProIDdisCha=true;
                        }
                    }
                    for(int i=bH+1;i<line.size();i++){
                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
                            Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
                            Tree[child].FN[i]=false;
                            ProIDdisCha=true;
                        }
                    }

                }else{//partial ancestor check

                    if(vertexIDChL.find(b)!=vertexIDChL.end()){
                        for(int i=0;i<bH;i++){
                            checkedDis.insert(make_pair(child,i));
                            if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
                                Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
                                Tree[child].FN[i]=false;
                                ProIDdisCha=true;
                            }
                        }
                    }
                    for(int i=bH+1;i<line.size();i++){
                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
                            Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
                            Tree[child].FN[i]=false;
                            ProIDdisCha=true;
                        }
                    }

                }
            }
        }
    }else{
        for(int k=0;k<Tree[child].vert.size();k++){
            int b=Tree[child].vert[k].first, bH=Tree[rank[b]].height-1,vbW=Tree[child].vert[k].second.first;
            if(Tree[child].FN[bH]){
                if(vertexIDChL.find(b)!=vertexIDChL.end()){
                    for(int i=0;i<bH;i++){
                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
                            Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
                            Tree[child].FN[i]=false;
                            ProIDdisCha=true;
                        }
                    }
                }
                for(int i=bH+1;i<line.size();i++){
                    checkedDis.insert(make_pair(child,i));
                    if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
                        Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
                        Tree[child].FN[i]=false;
                        ProIDdisCha=true;
                    }
                }
            }
        }
    }

    if(ProIDdisCha){
        vertexIDChL.insert(Tree[child].uniqueVertex);
    }

    line.push_back(Tree[child].uniqueVertex);
    for(int i=0;i<Tree[child].ch.size();i++){
        EachNodeProBDis5(Tree[child].ch[i], line, vertexIDChL,checkedDis);
    }
    line.pop_back();

}

void Gstartree::H2HincBatMT(vector<pair<pair<int,int>,pair<int,int>>>& wBatch) {
    int checknum = 0;
    map<pair<int, int>, int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
//OCdis.clear();

//NodeOrderss.clear();
    vector<set<int>> SCre; //SCre.clear();
    set<int> ss; //ss.clear();
    SCre.assign(node_num, ss);//{vertexID, set<int>}
    set<OrderCompMin> OC;
    OC.clear();//vertexID in decreasing node order

    for (int k = 0; k < wBatch.size(); k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int oldW = wBatch[k].second.first;
        int newW = wBatch[k].second.second;

        if (oldW < newW) {
            if(OverlayGraph[a].find(b)!=OverlayGraph[a].end()){
                if(OverlayGraph[a][b]!=oldW){//only works for no-boundary
                    cout<<"Inconsistent! "<<OverlayGraph[a][b]<<" "<<oldW<<endl; exit(1);
                }
                OverlayGraph[a][b]=newW;
            }else{
                cout<<"Wrong for Neighbors!"<<endl; exit(1);
            }
            if(OverlayGraph[b].find(a)!=OverlayGraph[b].end()){
                if(OverlayGraph[b][a]!=oldW){
                    cout<<"Inconsistent! "<<OverlayGraph[b][a]<<" "<<oldW<<endl; exit(1);
                }
                OverlayGraph[b][a]=newW;
            }else{
                cout<<"Wrong for Neighbors!"<<endl; exit(1);
            }

            int lid, hid;
            if (NodeOrder[a] < NodeOrder[b]) {
                lid = a; hid = b;
            } else {
                lid = b; hid = a;
            }

            for (int i = 0; i < Tree[rank[lid]].vert.size(); i++) {
                if (Tree[rank[lid]].vert[i].first == hid) {
                    if (Tree[rank[lid]].vert[i].second.first == oldW) {
                        Tree[rank[lid]].vert[i].second.second -= 1;
//                        cout<<lid<<" "<<hid<<" "<<Tree[rank[lid]].vert[i].second.first<<" "<<Tree[rank[lid]].vert[i].second.second<<endl;
                        if (Tree[rank[lid]].vert[i].second.second < 1) {
                            OCdis[make_pair(lid, hid)] = oldW;
                            SCre[lid].insert(hid);
                            OC.insert(OrderCompMin(lid));
                        }
                    }
                    break;
                }
            }
        }else{
            cout<<"Wrong for this edge update. "<<a<<" "<<b<<" "<<oldW<<" "<<newW<<endl; exit(1);
        }
    }

    vector<int> ProBeginVertexSet;
    ProBeginVertexSet.clear();
    vector<int> ProBeginVertexSetNew;
    bool influence;
    int ProID;
    vector<int> line;
    while (!OC.empty()) {
        ProID = (*OC.begin()).x;
        OC.erase(OC.begin());
        vector<pair<int, pair<int, int>>> Vert = Tree[rank[ProID]].vert;
        influence = false;

//each ProID corresponds to a line
        line.clear();
        line.reserve(heightMax);
        int pachid = ProID;
        while (Tree[rank[pachid]].height > 1) {
            line.insert(line.begin(), pachid);
            pachid = Tree[Tree[rank[pachid]].pa].uniqueVertex;
        }
        line.insert(line.begin(), pachid);

        for (auto it = SCre[ProID].begin(); it != SCre[ProID].end(); it++) {
            int Cid = *it;
            int Cw = OCdis[make_pair(ProID, Cid)];
            int cidH = Tree[rank[Cid]].height - 1;

            map<int, int> Hnei; //Hnei.clear();
            vector<pair<int, int>> Lnei; //Lnei.clear();
            for (int j = 0; j < Vert.size(); j++) {
                if (NodeOrder[Vert[j].first] > NodeOrder[Cid]) {
                    Hnei[Vert[j].first] = Vert[j].second.first;
                } else if (NodeOrder[Vert[j].first] < NodeOrder[Cid]) {
                    Lnei.push_back(make_pair(Vert[j].first, Vert[j].second.first));
                }
            }
//check the affected shortcuts
            int hid, lid;
            for (int j = 0; j < Tree[rank[Cid]].vert.size(); j++) {
                hid = Tree[rank[Cid]].vert[j].first;
                if (Hnei.find(hid) != Hnei.end()) {
                    if (Cw + Hnei[hid] == Tree[rank[Cid]].vert[j].second.first) {
                        Tree[rank[Cid]].vert[j].second.second -= 1;
                        if (Tree[rank[Cid]].vert[j].second.second < 1) {
                            SCre[Cid].insert(hid);
                            OC.insert(OrderCompMin(Cid));
                            OCdis[make_pair(Cid, hid)] = Cw + Hnei[hid];
                        }
                    }
                }
            }
            for (int j = 0; j < Lnei.size(); j++) {
                lid = Lnei[j].first;
                for (int k = 0; k < Tree[rank[lid]].vert.size(); k++) {
                    if (Tree[rank[lid]].vert[k].first == Cid) {
                        if (Tree[rank[lid]].vert[k].second.first == Cw + Lnei[j].second) {
                            Tree[rank[lid]].vert[k].second.second -= 1;
                            if (Tree[rank[lid]].vert[k].second.second < 1) {
                                SCre[lid].insert(Cid);
                                OC.insert(OrderCompMin(lid));
                                OCdis[make_pair(lid, Cid)] = Cw + Lnei[j].second;
                            }
                        }
                        break;
                    }
                }
            }


//before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
            if (Tree[rank[ProID]].FN[cidH]) {
                influence = true;
//higher than Cid
                for (int i = 0; i < cidH; i++) {
                    if (Tree[rank[ProID]].dis[i] == Cw + Tree[rank[Cid]].dis[i]) {
                        Tree[rank[ProID]].cnt[i] -= 1;
                    }
                }

//equal to Cid
                Tree[rank[ProID]].FN[cidH] = false;
                Tree[rank[ProID]].cnt[cidH] -= 1;

//lower than Cid
                for (int i = cidH + 1; i < Tree[rank[ProID]].dis.size(); i++) {
                    if (Tree[rank[ProID]].dis[i] == Cw + Tree[rank[line[i]]].dis[cidH]) {
                        Tree[rank[ProID]].cnt[i] -= 1;
                    }
                }
            }

//get the new value of shortcut
//	cout<<Cw<<" increase to ";
            Cw = INF;
            int countwt = 0;

//            for (int i = 0; i < Nodes[ProID].adjnodes.size(); i++) {///Neighbor
//                if (Nodes[ProID].adjnodes[i] == Cid) {
//                    Cw = Nodes[ProID].adjweight[i];//the weight value in the original graph
//                    countwt = 1;
//                    break;
//                }
//            }

            for (auto it = OverlayGraph[ProID].begin(); it != OverlayGraph[ProID].end(); ++it) {///Neighbor of overlay graph
                if (it->first == Cid) {
                    Cw = it->second;//the weight value in the original graph
                    countwt = 1;
                    break;
                }
            }

            int ssw, wtt, wid;
            vector<int> Wnodes;
            Wnodes.clear();

            if (ProID < Cid)
                Wnodes = SCconNodesMT[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
            else
                Wnodes = SCconNodesMT[Cid][ProID];
            if (Wnodes.size() > 0) {
                for (int i = 0; i < Wnodes.size(); i++) {
                    wid = Wnodes[i];
                    for (int j = 0; j < Tree[rank[wid]].vert.size(); j++) {
                        if (Tree[rank[wid]].vert[j].first == ProID) {
                            ssw = Tree[rank[wid]].vert[j].second.first;
                        }
                        if (Tree[rank[wid]].vert[j].first == Cid) {
                            wtt = Tree[rank[wid]].vert[j].second.first;
                        }
                    }

                    if (ssw + wtt < Cw) {
                        Cw = ssw + wtt;
                        countwt = 1;
                    } else if (ssw + wtt == Cw) {
                        countwt += 1;
                    }
                }
            }

//cout<<Cw<<endl;
//refresh the shortcut to the new value
            for (int i = 0; i < Tree[rank[ProID]].vert.size(); i++) {
                if (Tree[rank[ProID]].vert[i].first == Cid) {
                    Tree[rank[ProID]].vert[i].second.first = Cw;
                    Tree[rank[ProID]].vert[i].second.second = countwt;
                    break;
                }
            }
        }

        if (influence) {
            ProBeginVertexSetNew.clear();
            ProBeginVertexSetNew.reserve(ProBeginVertexSet.size() + 1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew = rank[ProID], r;
            for (int i = 0; i < ProBeginVertexSet.size(); i++) {
                r = rank[ProBeginVertexSet[i]];
                if (LCAQuery(rnew, r) != rnew) {
                    ProBeginVertexSetNew.push_back(ProBeginVertexSet[i]);
                }
            }
            ProBeginVertexSet = ProBeginVertexSetNew;
        }

    }

    int ProBeginVertexID;
    for (int i = 0; i < ProBeginVertexSet.size(); i++) {
        ProBeginVertexID = ProBeginVertexSet[i];
        vector<int> linee; //linee.clear();
        linee.reserve(heightMax);
        int pachidd = Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
        while (Tree[rank[pachidd]].height > 1) {
            linee.insert(linee.begin(), pachidd);
            pachidd = Tree[Tree[rank[pachidd]].pa].uniqueVertex;
        }
        linee.insert(linee.begin(), pachidd);

        eachNodeProcessIncrease1(rank[ProBeginVertexID], linee, checknum);
    }
//return checknum;
}

void Gstartree::eachNodeProcessIncrease1(int children, vector<int>& line, int& changelabel){
    int childID=Tree[children].uniqueVertex;
    int childH=Tree[children].height-1;
    for(int i=0;i<Tree[children].dis.size();i++){
        if(Tree[children].cnt[i]==0){
            changelabel+=1;
            //firstly, check which dis can be infected
            int disBF=Tree[children].dis[i];
            int PID;
            //chidlID
            for(int k=0;k<VidtoTNid[childID].size();k++){
                PID=VidtoTNid[childID][k];
                if(Tree[PID].FN[childH] && Tree[PID].dis[i]==disBF+Tree[PID].dis[childH]){
                    Tree[PID].cnt[i]-=1;
                }
            }

            //line[i]
            for(int k=0;k<VidtoTNid[line[i]].size();k++){
                PID=VidtoTNid[line[i]][k];
                if(Tree[PID].height>Tree[children].height && Tree[PID].vAncestor[childH] == childID){
                    if(Tree[PID].FN[i] && Tree[PID].dis[childH]==disBF+Tree[PID].dis[i]){
                        Tree[PID].cnt[childH]-=1;
                    }
                }
            }

            //secondly, calculate the actual distance
            int dis=INF; int count=0;
            int Dvb; int b,bH; int DDvb=INF;
            for(int j=0;j<Tree[children].vert.size();j++){
                Dvb=Tree[children].vert[j].second.first;
                b=Tree[children].vert[j].first;
                bH=Tree[rank[b]].height-1;
                if(bH<i){
                    if(Dvb+Tree[rank[line[i]]].dis[bH]<dis){
                        dis=Dvb+Tree[rank[line[i]]].dis[bH];
                        count=1;
                    }else if(Dvb+Tree[rank[line[i]]].dis[bH]==dis){
                        count+=1;
                    }
                }else if(bH==i){
                    DDvb=Dvb;
                    if(Dvb<dis){
                        dis=Dvb;
                        count=1;
                    }else if(Dvb==dis){
                        count+=1;
                    }
                }else{
                    if(Dvb+Tree[rank[b]].dis[i]<dis){
                        dis=Dvb+Tree[rank[b]].dis[i];
                        count=1;
                    }else if(Dvb+Tree[rank[b]].dis[i]==dis){
                        count+=1;
                    }
                }
            }
            if(DDvb==dis) Tree[children].FN[i]=true;
            Tree[children].dis[i]=dis;
            Tree[children].cnt[i]=count;
        }
    }

    line.push_back(childID);
    for(int i=0;i<Tree[children].ch.size();i++){
        eachNodeProcessIncrease1(Tree[children].ch[i],line,changelabel);
    }
    line.pop_back();
}

#endif //GSTARTREE_HPP
