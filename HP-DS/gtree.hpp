//
// Created by Xinjie ZHOU on 25/8/2022.
//

#ifndef GTREE_HPP
#define GTREE_HPP
// macro for 64 bits file, larger than 2G
#define _FILE_OFFSET_BITS 64

#include <metis.h>
#include "gtree.h"

extern bool ifDebug;

// use for metis
// idx_t = int64_t / real_t = double
idx_t nvtxs; // |vertices|
idx_t ncon; // number of weight per vertex
idx_t* xadj; // array of adjacency of indices
idx_t* adjncy; // array of adjacency nodes
idx_t* vwgt; // array of weight of nodes
idx_t* adjwgt; // array of weight of edges in adjncy
idx_t nparts; // number of parts to partition
idx_t objval; // edge cut for partitioning solution
idx_t* part; // array of partition vector
idx_t options[METIS_NOPTIONS]; // option array

void Gtree::PathInit(bool ifMac){
    if(dataset == "cal"){
        if(ifMac){
            DataPath = "/Users/zhouxj/Documents/1-Research/1-Papers/0-My_Papers/2-EPSP/GTree-master/src/gtree/";
        }
//        else{
//            DataPath = "/home/s4451682/Xinjie/Hier_PLL/";
//        }
    }else{
        if(ifMac){
            DataPath = "/Users/zhouxj/Documents/1-Research/Datasets/";
        }
//        else{
//            DataPath = "/media/TraminerData/s4451682/GraphDataforPartition/Processed/";
//        }
    }
    FILE_NODE = dataset + ".cnode";//const
    FILE_EDGE = dataset + ".cedge";
// gtree index disk storage
    FILE_NODES_GTREE_PATH = dataset + "." + to_string(parti_num) + ".paths.bin";
    FILE_GTREE = dataset + "." + to_string(parti_num) +".gtree.bin";
    FILE_ONTREE_MIND = dataset + "." + to_string(parti_num) + ".minds.bin";
// input
    FILE_OBJECT = dataset + "." + to_string(parti_num) + ".object";
    FILE_SHORTCUT = dataset + "." + to_string(parti_num)+ ".shortcuts";
    FILE_QUERY = dataset + ".query";
    FILE_UPDATE = dataset + ".update";
}

// METIS setting options
void Gtree::options_setting(){
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY; // Multilevel k-way partitioning
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT; // Edge-cut minimization
    options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM; // Sorted heavy-edge matching
    options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_GROW; // Computes a bisection at random followed by a refinement METIS_IPTYPE_RANDOM
    options[METIS_OPTION_RTYPE] = METIS_RTYPE_FM; // FM-based cut refinement
    // options[METIS_OPTION_NCUTS] = 1; // the number of different partitionings that it will compute.
    // options[METIS_OPTION_NITER] = 10; // the number of iterations for the refinement algorithms at each stage of the uncoarsening process.
    /* balance factor, used to be 500 */
//	options[METIS_OPTION_UFACTOR] = 5;//500; // Specifies the maximum allowed load imbalance among the partitions. default value is 30
    // options[METIS_OPTION_MINCONN];
    options[METIS_OPTION_CONTIG] = 1; // the partitioning routines should try to produce partitions that are contiguous
    // options[METIS_OPTION_SEED];
    options[METIS_OPTION_NUMBERING] = 0; // indicate which numbering scheme is used for the adjacency structure of a graph, start from 0
    // options[METIS_OPTION_DBGLVL] = 0; // Specifies the amount of progress/debugging information will be printed during the execution of the algorithms.
}

// input init
void Gtree::load_graph(){
    FILE *fin;
    // load node
    printf("LOADING NODE...\t");
    char filename[200];
    if(dataset == "cal"){
        strcpy(filename, DataPath.c_str());
        strcat(filename, FILE_NODE.c_str());
    }else {
        strcpy(filename, dirname.c_str());
        strcat(filename, "/");
        strcat(filename, FILE_NODE.c_str());
    }
    fin = fopen(filename,"r");
    if(fin == NULL){
        cout<<"Failed to open file "<<filename<<endl;
        assert(fin != NULL);
        exit(1);
    }
    int nid;
    double x,y;
    while( fscanf(fin, "%d %lf %lf", &nid, &x, &y ) == 3 ){
        Node node = { x, y };
        Nodes.push_back(node);
    }
    fclose(fin);
    printf("COMPLETE. Node number = %d\n", (int)Nodes.size());
    node_num = Nodes.size();

    // load edge
    printf("LOADING EDGE...\t");
    char filename2[200];
    if(dataset == "cal"){
        strcpy(filename2, DataPath.c_str());
        strcat(filename2, FILE_EDGE.c_str());
    }else {
        strcpy(filename2, dirname.c_str());
        strcat(filename2, "/");
        strcat(filename2, FILE_EDGE.c_str());
    }
//	fin = fopen(FILE_EDGE, "r");
    fin = fopen(filename2,"r");
    if(fin == NULL){
        cout<<"Failed to open file "<<filename2<<endl;
        assert(fin != NULL);
        exit(1);
    }
    int eid;
    int snid, enid;
    double weight;
    int iweight;
    noe = 0;
    while( fscanf(fin,"%d %d %d %lf", &eid, &snid, &enid, &weight ) == 4 ){
        ++noe;
        assert(enid<node_num &&snid <node_num);
        iweight = (int) (weight * WEIGHT_INFLATE_FACTOR );
        Nodes[snid].adjnodes.push_back( enid );
        Nodes[snid].adjweight.push_back( iweight );
        Nodes[enid].adjnodes.push_back( snid );
        Nodes[enid].adjweight.push_back( iweight );
    }
    fclose(fin);
    printf("COMPLETE. Edge number = %d\n", noe);
    edge_num = noe;
}
//read weighted graph
void Gtree::ReadGraph_W() {

    //normal distribution generator
    Timer tt;
    tt.start();
    string lineStr;
    string line_symbol;
    string temp_str;
    int ID1, ID2, weight;

    string r_file = string(DataPath)  + "/" + dataset + "/" + dataset;
//    string r_file = string(DataPath)  + dataset + "/" + dataset + "_Processed";

    ifstream inFile(r_file, ios::in);
    if (!inFile) { // if not exist
        cout << "Fail to open file" << r_file << endl;
        exit(1);
    }

    /// read graph and recording the degree of vertices
    cout << dataset<<" graph edges Data loading..." << endl;
    inFile >> node_num >> edge_num;

    LEAF_CAP = node_num/parti_num;
    cout<<"Leaf node maximum size: "<<LEAF_CAP<<endl;
    Nodes.resize(node_num);
    while (inFile) {//read each line
        inFile >> ID1 >> ID2 >> weight;
        assert(weight > 0);
        if(dataset != "FRIEND"){
            Nodes[ID1].adjnodes.push_back( ID2 );
            Nodes[ID1].adjweight.push_back( weight );
        }else{
            Nodes[ID1].adjnodes.push_back( ID2 );
            Nodes[ID1].adjweight.push_back( weight );
            Nodes[ID2].adjnodes.push_back( ID1 );
            Nodes[ID2].adjweight.push_back( weight );
        }
    }
    inFile.close();

    tt.stop();
    cout << "The time used for data reading: " << tt.GetRuntime() << " s." << endl;
    cout << "Node number: " << node_num << ", Edge number: " << edge_num << endl;
    cout << "--------------------" << endl;

}

// transform original data format to that suitable for METIS
void Gtree::data_transform_init( set<int> &nset ){
    // nvtxs, ncon
    nvtxs = nset.size();
    ncon = 1; // The number of balancing constraints. It should be at least 1.
//    cout << "test 1.1" << endl;
    xadj = new idx_t[nset.size() + 1]; // index array
//    cout << "test 1.2" << endl;
    if(noe != -1){
        adjncy = new idx_t[noe * 2]; // adjacency array
        adjwgt = new idx_t[noe * 2]; // adjacency weight
    }else{
        adjncy = new idx_t[edge_num+1]; // adjacency array
        adjwgt = new idx_t[edge_num+1]; // adjacency weight
    }
//    cout << "test 1.3" << endl;
    int xadj_pos = 1;
    int xadj_accum = 0;
    int adjncy_pos = 0;

    // xadj, adjncy, adjwgt
    unordered_map<int,int> nodemap;
    nodemap.clear();

    xadj[0] = 0;
    int i = 0;
    int nid, fanout, enid;
    for ( set<int>::iterator it = nset.begin(); it != nset.end(); ++it, ++i ){
        // init node map
        nodemap[*it] = i; //new vertex map

        nid = *it;
        fanout = Nodes[nid].adjnodes.size();
        for ( int j = 0; j < fanout; ++j ){
            enid = Nodes[nid].adjnodes[j];
            // ensure edges within
            if ( nset.find( enid ) != nset.end() ){ // if not found
                ++xadj_accum;

                adjncy[adjncy_pos] = enid; // record original adjacency id
                adjwgt[adjncy_pos] = Nodes[nid].adjweight[j];
                ++adjncy_pos;
            }
        }
        xadj[xadj_pos++] = xadj_accum;
    }
//    cout << "test 1.4" << endl;
    // adjust nodes number started by 0
    for ( int j = 0; j < adjncy_pos; ++j ){
        if(nodemap[adjncy[j]] != adjncy[j]){
            adjncy[j] = nodemap[adjncy[j]]; // update the adjacency node id to the mapped one
        }
    }

    // adjwgt -> 1
    if (ADJWEIGHT_SET_TO_ALL_ONE){ // if it is unweighted graph
        for ( int j = 0; j < adjncy_pos; ++j ){
            adjwgt[j] = 1;
        }
    }

    // nparts
    nparts = PARTITION_PART;
//    cout << "test 1.5" << endl;
    // part
    part = new idx_t[nset.size()];
//    cout << "test 1.6" << endl;
}

void Gtree::init(){
    if(dataset == "cal"){
        load_graph();
    }else{
        ReadGraph_W();
    }
    options_setting();
}

void Gtree::finalize(){
    delete xadj;
    delete adjncy;
    delete adjwgt;
    delete part;
}

// graph partition
// input: nset = a set of node id
// output: <node, node belong to partition id>
unordered_map<int,int> Gtree::graph_partition( set<int> &nset ){
    unordered_map<int,int> result;

    // transform data to metis
    data_transform_init( nset );
//    cout << "test 2." << endl;
    // partition, result -> part
    // k way partition, k=nparts=4
    METIS_PartGraphKway(
            &nvtxs,
            &ncon,
            xadj,
            adjncy,
            NULL,
            NULL,
            adjwgt,
            &nparts,
            NULL,
            NULL,
            options,
            &objval,
            part
    );
//    cout << "test 3." << endl;
    // push to result
    result.clear();
    int i = 0;
    for ( set<int>::iterator it = nset.begin(); it != nset.end(); it++, i++ ){
        result[*it] = part[i];
    }

    // finalize
    finalize();

    return result;
}

// gtree construction
void Gtree::build(){
    // init root
    TreeNode root;
    root.isleaf = false;
    root.father = -1;
    GTree.push_back(root);

    // init stack
    stack<Status> buildstack;//last in, first out; first in, last out
    Status rootstatus;
    rootstatus.tnid = 0;
    rootstatus.nset.clear();
    for ( int i = 0; i < Nodes.size(); i++ ){//root node
        rootstatus.nset.insert(i);
    }
    buildstack.push( rootstatus );

    // start to build
    unordered_map<int,int> presult;
    set<int> childset[PARTITION_PART];


    while( buildstack.size() > 0 ){
        // pop top
        Status current = buildstack.top();//current status
        buildstack.pop();

        // update gtreepath
        for ( set<int>::iterator it = current.nset.begin(); it != current.nset.end(); it++ ){
            Nodes[*it].gtreepath.push_back( current.tnid );
        }

        // check cardinality
        if ( current.nset.size() <= LEAF_CAP ){//if the size is smaller than leaf node limitation
            // build leaf node
            GTree[current.tnid].isleaf = true;
            GTree[current.tnid].leafnodes.clear();
            for ( set<int>::iterator it = current.nset.begin(); it != current.nset.end(); it++ ){
                GTree[current.tnid].leafnodes.push_back( *it );
            }
            continue;
        }
//        cout << "test 1." << endl;
        /// graph partitioning by METIS
//		printf("PARTITIONING...NID=%d...SIZE=%d...", current.tnid, (int)current.nset.size() );
        presult = graph_partition( current.nset );
//		printf("COMPLETE.\n");

        // construct child node set
        for ( int i = 0; i < PARTITION_PART; i++ ){
            childset[i].clear();
        }
        int slot;
        for ( set<int>::iterator it = current.nset.begin(); it != current.nset.end(); it++ ){
            slot = presult[*it];//partition id
            childset[slot].insert(*it);
        }

        // generate child tree nodes
        int childpos;
        for ( int i = 0; i < PARTITION_PART; i++ ){
            TreeNode tnode;
            tnode.isleaf = false;
            tnode.father = current.tnid;

            // insert to GTree first
            GTree.push_back(tnode);
            childpos = GTree.size() - 1; // tree node position
            GTree[current.tnid].children.push_back( childpos );

            // calculate border nodes
            GTree[childpos].borders.clear();
            for ( set<int>::iterator it = childset[i].begin(); it != childset[i].end(); it++ ){//for each child partition

                bool isborder = false;
                for ( int j = 0; j < Nodes[*it].adjnodes.size(); j++ ){
                    if ( childset[i].find( Nodes[*it].adjnodes[j] ) == childset[i].end() ){//if exist a neighbor vertex that is not in the same tree node
                        isborder = true;
                        break;
                    }
                }
                if ( isborder ){//if it is border vertex
                    GTree[childpos].borders.push_back(*it);
                    // update globally
                    Nodes[*it].isborder = true;
                }
            }

            // add to stack
            Status ongoingstatus;
            ongoingstatus.tnid = childpos;//new tree id
            ongoingstatus.nset = childset[i];
            buildstack.push(ongoingstatus);

        }

    }


}

// dump gtree index to file
void Gtree::gtree_save(){
    // FILE_GTREE
    char filename[200];
    if(dataset == "cal"){
        strcpy(filename, DataPath.c_str());
        strcat(filename, FILE_GTREE.c_str());
    }else {
        strcpy(filename, dirname.c_str());
        strcat(filename, "/");
        strcat(filename, FILE_GTREE.c_str());
    }
    FILE *fout = fopen( filename, "wb" );
//	FILE *fout = fopen( FILE_GTREE, "wb" );
    if(fout == NULL){
        cout<<"Failed to open file "<<filename<<endl;
        exit(1);
    }
    int *buf = new int[ Nodes.size() ];
    for ( int i = 0; i < GTree.size(); i++ ){
        // borders
        int count_borders = GTree[i].borders.size();
        fwrite( &count_borders, sizeof(int), 1, fout );
        copy( GTree[i].borders.begin(), GTree[i].borders.end(), buf );
        fwrite( buf, sizeof(int), count_borders, fout );
        // children
        int count_children = GTree[i].children.size();
        fwrite( &count_children, sizeof(int), 1, fout );
        copy( GTree[i].children.begin(), GTree[i].children.end(), buf );
        fwrite( buf, sizeof(int), count_children, fout );
        // isleaf
        fwrite( &GTree[i].isleaf, sizeof(bool), 1, fout );
        // leafnodes
        int count_leafnodes = GTree[i].leafnodes.size();
        fwrite( &count_leafnodes, sizeof(int), 1, fout );
        copy( GTree[i].leafnodes.begin(), GTree[i].leafnodes.end(), buf );
        fwrite( buf, sizeof(int), count_leafnodes, fout );
        // father
        fwrite( &GTree[i].father, sizeof(int), 1, fout );
    }
    fclose(fout);

    // FILE_NODES_GTREE_PATH
    char filename2[200];
    if(dataset == "cal"){
        strcpy(filename2, DataPath.c_str());
        strcat(filename2, FILE_NODES_GTREE_PATH.c_str());
    }else {
        strcpy(filename2, dirname.c_str());
        strcat(filename2, "/");
        strcat(filename2, FILE_NODES_GTREE_PATH.c_str());
    }
    fout = fopen( filename2, "wb" );
//    fout = fopen( (DataPath+FILE_NODES_GTREE_PATH).c_str(), "wb" );
    if(fout == NULL){
        cout<<"Failed to open file "<<filename2<<endl;
        exit(1);
    }
    for ( int i = 0; i < Nodes.size(); i++ ){
        int count = Nodes[i].gtreepath.size();
        fwrite( &count, sizeof(int), 1, fout );
        copy( Nodes[i].gtreepath.begin(), Nodes[i].gtreepath.end(), buf );
        fwrite( buf, sizeof(int), count, fout );
    }
    fclose(fout);
    delete[] buf;
}

// load gtree index from file
/*void Gtree::gtree_load(){
    // FILE_GTREE
    char filename[200];
    if(dataset == "cal"){
        strcpy(filename, DataPath.c_str());
        strcat(filename, FILE_GTREE.c_str());
    }else {
        strcpy(filename, DataPath.c_str());
        strcat(filename, dataset.c_str());
        strcat(filename, "/");
        strcat(filename, FILE_GTREE.c_str());
    }
    FILE *fin = fopen( filename, "rb" );
//	FILE *fin = fopen( FILE_GTREE, "rb" );
    if(fin == NULL){
        cout<<"Failed to open file "<<filename<<endl;
        exit(1);
    }
    cout << "Loading G-tree... ";
    int *buf = new int[ Nodes.size() ];
    int count_borders, count_children, count_leafnodes;
    bool isleaf;
    int father;

    // clear gtree
    GTree.clear();

    while( fread( &count_borders, sizeof(int), 1, fin ) ){
        TreeNode tn;
        // borders
        tn.borders.clear();
        fread( buf, sizeof(int), count_borders, fin );
        for ( int i = 0; i < count_borders; i++ ){
            tn.borders.push_back(buf[i]);
        }
        // children
        fread( &count_children, sizeof(int), 1, fin );
        fread( buf, sizeof(int), count_children, fin );
        for ( int i = 0; i < count_children; i++ ){
            tn.children.push_back(buf[i]);
        }
        // isleaf
        fread( &isleaf, sizeof(bool), 1, fin );
        tn.isleaf = isleaf;
        // leafnodes
        fread( &count_leafnodes, sizeof(int), 1, fin );
        fread( buf, sizeof(int), count_leafnodes, fin );
        for ( int i = 0; i < count_leafnodes; i++ ){
            tn.leafnodes.push_back(buf[i]);
        }
        // father
        fread( &father, sizeof(int), 1, fin );
        tn.father = father;

        GTree.push_back(tn);
    }
    fclose(fin);

    // FILE_NODES_GTREE_PATH
    int count;
    char filename2[200];
    if(dataset == "cal"){
        strcpy(filename2, DataPath.c_str());
        strcat(filename2, FILE_NODES_GTREE_PATH.c_str());
    }else {
        strcpy(filename2, DataPath.c_str());
        strcat(filename2, dataset.c_str());
        strcat(filename2, "/");
        strcat(filename2, FILE_NODES_GTREE_PATH.c_str());
    }
    fin = fopen( filename2, "wb" );
//    fin = fopen( (DataPath+FILE_NODES_GTREE_PATH).c_str(), "rb" );
    if(fin == NULL){
        cout<<"Failed to open file "<<filename2<<endl;
        exit(1);
    }
    int pos = 0;
    while( fread( &count, sizeof(int), 1, fin ) ){
        fread( buf, sizeof(int), count, fin );
        // clear gtreepath
        Nodes[pos].gtreepath.clear();
        for ( int i = 0; i < count; i++ ){
            Nodes[pos].gtreepath.push_back( buf[i] );
        }
        // pos increase
        pos ++;
    }
    fclose(fin);
    delete[] buf;
    cout << "Done." << endl;
}*/

// load gtree index from file
void Gtree::load_gtreeQ() {
    // FILE_GTREE
    char filename[200];
    if(dataset == "cal"){
        strcpy(filename, DataPath.c_str());
        strcat(filename, FILE_GTREE.c_str());
    }else {
        strcpy(filename, dirname.c_str());
        strcat(filename, "/");
        strcat(filename, FILE_GTREE.c_str());
    }
    FILE *fin = fopen( filename, "rb" );
//    FILE *fin = fopen(FILE_GTREE, "rb");
    if(fin == NULL){
        cout << "Failed to open file " << filename << endl;
        exit(1);
    }
    cout << "Loading G-tree... ";
    int *buf = new int[Nodes.size()];
    int count_borders, count_children, count_leafnodes;
    bool isleaf;
    int father;

    // clear gtree
    GTree.clear();

    int node_count = 0;

    while (fread(&count_borders, sizeof(int), 1, fin)) {
        TreeNode tn;
        // borders
        tn.borders.clear();
        fread(buf, sizeof(int), count_borders, fin);
        for (int i = 0; i < count_borders; i++) {
            tn.borders.push_back(buf[i]);
        }
        // children
        fread(&count_children, sizeof(int), 1, fin);
        fread(buf, sizeof(int), count_children, fin);
        tn.children.clear();///
        for (int i = 0; i < count_children; i++) {
            tn.children.push_back(buf[i]);
        }
        // isleaf
        fread(&isleaf, sizeof(bool), 1, fin);
        tn.isleaf = isleaf;
        // leafnodes
        fread(&count_leafnodes, sizeof(int), 1, fin);
        fread(buf, sizeof(int), count_leafnodes, fin);
        tn.leafnodes.clear();///
        for (int i = 0; i < count_leafnodes; i++) {
            tn.leafnodes.push_back(buf[i]);
            Nodes[buf[i]].inleafpos = i;
        }
        // father
        fread(&father, sizeof(int), 1, fin);
        tn.father = father;

        GTree.push_back(tn);
    }
    fclose(fin);

    // FILE_NODES_GTREE_PATH
    int count;
    char filename2[200];
    if(dataset == "cal"){
        strcpy(filename2, DataPath.c_str());
        strcat(filename2, FILE_NODES_GTREE_PATH.c_str());
    }else {
        strcpy(filename2, dirname.c_str());
        strcat(filename2, "/");
        strcat(filename2, FILE_NODES_GTREE_PATH.c_str());
    }
    fin = fopen( filename2, "rb" );
    if(fin == NULL){
        cout << "Failed to open file " << filename2 << endl;
        exit(1);
    }
    int pos = 0;
    while (fread(&count, sizeof(int), 1, fin)) {
        fread(buf, sizeof(int), count, fin);
        // clear gtreepath
        Nodes[pos].gtreepath.clear();
        for (int i = 0; i < count; i++) {
            Nodes[pos].gtreepath.push_back(buf[i]);
        }
        Nodes[pos].deep = count;
        Nodes[pos].inleaf = buf[count - 1];
        // pos increase
        pos++;
    }
    fclose(fin);
    delete[] buf;
    cout << "Done."<<endl;
}

// up_pos & current_pos(used for quickly locating parent & child nodes)
void Gtree::build_up_and_down_pos() {
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
}

// dijkstra search, used for single-source shortest path search WITHIN one gtree leaf node!//old version
// input: s = source node
//        cands = candidate node list
//        graph = search graph(this can be set to subgraph)
/*vector<int> Gtree::dijkstra_candidate( int s, vector<int> &cands, vector<Node> &graph ){
    // init
    set<int> todo;
    todo.clear();
    todo.insert(cands.begin(), cands.end());//get the partition vertex set

    unordered_map<int,int> result;
    result.clear();
    set<int> visited;
    visited.clear();
    unordered_map<int,int> q;//map used for mapping vertex id to its distance from source vertex
    q.clear();
    q[s] = 0;//map of source vertex

    // start
    int min, minpos, adjnode, weight;
    while( ! todo.empty() && ! q.empty() ){
        min = -1;
        //linear scan for finding vertex with minimal distance
        for ( unordered_map<int,int>::iterator it = q.begin(); it != q.end(); it ++ ){
            if ( min == -1 ){
                minpos = it -> first; //vertex id
                min = it -> second; //distance from source vertex
            }
            else{
                if ( it -> second < min ){
                    min = it -> second;
                    minpos = it -> first;
                }
            }
        }

        // put min to result, add to visited
        result[minpos] = min;
        visited.insert( minpos );
        q.erase(minpos);

        if ( todo.find( minpos ) != todo.end() ){//if found, erase visited vertex
            todo.erase( minpos );
        }

        // expand on graph (the original graph)
        for ( int i = 0; i < graph[minpos].adjnodes.size(); i++ ){
            adjnode = graph[minpos].adjnodes[i];
            if ( visited.find( adjnode ) != visited.end() ){//if found, ie, it is visited
                continue;
            }
            weight = graph[minpos].adjweight[i];//edge weight

            if ( q.find(adjnode) != q.end() ){//if found in q
                if ( min + weight < q[adjnode] ){
                    q[adjnode] = min + weight;
                }
            }
            else{
                q[adjnode] = min + weight;
            }

        }
    }

    // output
    vector<int> output;
    for ( int i = 0; i < cands.size(); i++ ){
        output.push_back( result[cands[i]] );//only push the distance result of vertex in cands (i.e., in this partition)
    }

    // return
    return output;//vector of distance matrix values
}*/
// heap-based dijkstra search, used for single-source shortest path search WITHIN one gtree leaf node!
// input: s = source node
//        cands = candidate node list
//        graph = search graph(this can be set to subgraph)
// output: the distance vector from s to all borders
vector<int> Gtree::dijkstra_candidate( int s, vector<int> &cands, vector<Node> &graph ){
    // init
    set<int> todo;
    todo.clear();
    todo.insert(cands.begin(), cands.end());//get the partition vertex set

    unordered_map<int,int> result;
    result.clear();
    set<int> visited;
    visited.clear();
//    unordered_map<int,int> q;//map used for mapping vertex id to its distance from source vertex
    benchmark::heap<2, int, int> q(node_num);
//    vector<bool> closed(node_num, false); //flag vector of whether closed
    vector<Distance> cost(node_num, INF);   //vector of cost
    int temp_dis;
    q.clear();
    q.update(s,0);
    cost[s] = 0;
//    q[s] = 0;//map of source vertex

    // start
    int min, minpos, adjnode, weight;
    while( ! todo.empty() && ! q.empty() ){
        min = -1;
        q.extract_min(minpos,min);

        // put min to result, add to visited
        result[minpos] = min;
        visited.insert( minpos );
//        q.erase(minpos);

        if ( todo.find( minpos ) != todo.end() ){//if found, erase visited vertex
            todo.erase( minpos );
        }

        // expand on graph (the original graph)
        for ( int i = 0; i < graph[minpos].adjnodes.size(); i++ ){
            adjnode = graph[minpos].adjnodes[i];
            if ( visited.find( adjnode ) != visited.end() ){//if found, ie, it is visited
                continue;
            }
            weight = graph[minpos].adjweight[i];//edge weight
            temp_dis = min + weight;
            if ( temp_dis < cost[adjnode] ){
                q.update(adjnode,temp_dis);
                cost[adjnode] = temp_dis;
            }

        }
    }

    // output
    vector<int> output;
    for ( int i = 0; i < cands.size(); i++ ){
        output.push_back( result[cands[i]] );//only push the distance result of vertex in cands
    }

    // return
    return output;//vector of distance matrix values
}
// calculate the distance matrix, algorithm shown in section 5.2 of paper
void Gtree::hierarchy_shortest_path_calculation(bool ifParallel){
    // level traversal
    vector< vector<int> > treenodelevel;//hierarchical tree

    vector<int> current;
    current.clear();
    current.push_back(0);//root node
    treenodelevel.push_back(current);

    vector<int> mid;
    while( current.size() != 0 ){
        mid = current;// intermediate tree node
        current.clear();
        for ( int i = 0; i < mid.size(); i++ ){
            for ( int j = 0; j < GTree[mid[i]].children.size(); j++ ){
                current.push_back( GTree[mid[i]].children[j] );
            }
        }
        if ( current.size() == 0 ) break;
        treenodelevel.push_back( current );
    }
    /// It seems that G-tree is an unbalanced tree
    // bottom up calculation
    // temp graph
    vector<Node> graph;
    graph = Nodes;//original graph
    vector<int> cands;
    vector<int> result;
    unordered_map<int, unordered_map<int,int> > vertex_pairs;//result of distance matrix

    // do dijkstra
    int s, t, tn, nid, cid, weight;
    vector<int> tnodes, tweight;
    set<int> nset;

    if(ifParallel){/// multi-thread
        int threadnum = 0;
        vector<vector<int>> tnVector(threadnum);
        vector<vector<unordered_map<int, unordered_map<int,int> >>> vertex_pairsVV;//result of distance matrix
        for ( int i = treenodelevel.size() - 1; i >= 0; i-- ){//start from the lowest level
            boost::thread_group threadf;
            if(treenodelevel[i].size() > thread_num){
                threadnum = thread_num;
            }else{
                threadnum = treenodelevel[i].size();
            }
            tnVector.assign(threadnum,vector<int>());
            vertex_pairsVV.assign(threadnum,vector<unordered_map<int, unordered_map<int,int> >>());
            int thread_i = 0;
            for ( int j = 0; j < treenodelevel[i].size(); j++ ) {//for each partition in this level
                tn = treenodelevel[i][j];
                tnVector[thread_i].push_back(tn);
                if(thread_i >= threadnum-1){
                    thread_i = 0;
                }else{
                    thread_i++;
                }
            }
            for(int ti=0;ti<threadnum;++ti){
                vertex_pairsVV[ti].assign(tnVector[ti].size(),unordered_map<int, unordered_map<int,int> >());
                /// multi-thread
                threadf.add_thread(new boost::thread(&Gtree::hierarchy_shortest_path_compute, this, i, boost::ref(tnVector[ti]), boost::ref(graph), boost::ref(vertex_pairsVV), ti));
            }
            threadf.join_all();
            for(int ti=0;ti<threadnum;++ti) {
                int ii = 0;
                for (auto it = tnVector[ti].begin(); it != tnVector[ti].end(); ++it, ++ii) {
                    vertex_pairs = vertex_pairsVV[ti][ii];
                    int tn = *it;
                    // IMPORTANT! after all border finished, degenerate graph
                    // first, remove inward edges
                    for (int k = 0; k < GTree[tn].borders.size(); k++) {//for each border vertex
                        s = GTree[tn].borders[k];
                        tnodes.clear();
                        tweight.clear();
                        for (int p = 0; p < graph[s].adjnodes.size(); p++) {//for each adjacent vertex
                            nid = graph[s].adjnodes[p];
                            weight = graph[s].adjweight[p];
                            // if adj node in same tree node

                            if (graph[nid].gtreepath.size() <= i || graph[nid].gtreepath[i] !=
                                                                    tn) {// add the higher-level nodes or other node in the same level
                                // only remain those useful vertices, i.e., borders
                                tnodes.push_back(nid);
                                tweight.push_back(weight);
                            }
                        }
                        // cut it
                        graph[s].adjnodes = tnodes;//update the adjacency lists of graph, only left the useful boundary vertices
                        graph[s].adjweight = tweight;
                    }
                    // second, add inter connected edges (shortcuts)
                    for (int k = 0; k < GTree[tn].borders.size(); k++) {
                        for (int p = 0; p < GTree[tn].borders.size(); p++) {
                            if (k == p) continue;
                            s = GTree[tn].borders[k];
                            t = GTree[tn].borders[p];
                            graph[s].adjnodes.push_back(t);
                            graph[s].adjweight.push_back(vertex_pairs[s][t]);
                        }
                    }
                }
            }
        }
    }
    else{/// single thread
        for ( int i = treenodelevel.size() - 1; i >= 0; i-- ){//start from the lowest level
            for ( int j = 0; j < treenodelevel[i].size(); j++ ){//for each partition in this level
                tn = treenodelevel[i][j];

                cands.clear();
                if ( GTree[tn].isleaf ){//for leaf node
                    // cands = leafnodes
                    cands = GTree[tn].leafnodes;//for leaf node, the cands is the leafnodes
                    // union borders = borders;
                    GTree[tn].union_borders = GTree[tn].borders;
                }
                else{//for non-leaf node
                    nset.clear();
                    for ( int k = 0; k < GTree[tn].children.size(); k++ ){
                        cid = GTree[tn].children[k];
                        nset.insert( GTree[cid].borders.begin(), GTree[cid].borders.end() );//get the borders of all children
                    }
                    // union borders = cands;

                    cands.clear();
                    for ( set<int>::iterator it = nset.begin(); it != nset.end(); it ++ ){//push each vertex of this partition into cands
                        cands.push_back( *it );//for non-leaf node, the cands is the borders of all child nodes
                    }
                    GTree[tn].union_borders = cands;//for non-leaf node, the union_borders contains all the borders of children
                }

                // start to do min dis
                vertex_pairs.clear();

                // for each border, do min dis
                int cc = 0;

                for ( int k = 0; k < GTree[tn].union_borders.size(); k++ ){
                    //printf("DIJKSTRA...LEAF=%d BORDER=%d\n", tn, GTree[tn].union_borders[k] );
                    result = dijkstra_candidate( GTree[tn].union_borders[k], cands, graph );//distance vector from s to all borders
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
                    for ( int p = 0; p < graph[s].adjnodes.size(); p++ ){//for each adjacent vertex
                        nid = graph[s].adjnodes[p];
                        weight = graph[s].adjweight[p];
                        // if adj node in same tree node

                        if ( graph[nid].gtreepath.size() <= i || graph[nid].gtreepath[i] != tn ){// add the higher-level nodes or other node in the same level
                            // only remain those useful vertices, i.e., borders
                            tnodes.push_back(nid);
                            tweight.push_back(weight);
                        }
                    }
                    // cut it
                    graph[s].adjnodes = tnodes;//update the adjacency lists of graph, only left the useful boundary vertices
                    graph[s].adjweight = tweight;
                }
                // second, add inter connected edges (shortcuts)
                for ( int k = 0; k < GTree[tn].borders.size(); k++ ){
                    for ( int p = 0; p < GTree[tn].borders.size(); p++ ){
                        if ( k == p ) continue;
                        s = GTree[tn].borders[k];
                        t = GTree[tn].borders[p];
                        graph[s].adjnodes.push_back( t );
                        graph[s].adjweight.push_back( vertex_pairs[s][t] );
                    }
                }
            }
        }
    }


}

void Gtree::hierarchy_shortest_path_calculate( bool ifParallel){
    // level traversal
    vector< vector<int> > treenodelevel;//hierarchical tree

    vector<int> current;
    current.clear();
    current.push_back(0);//root node
    treenodelevel.push_back(current);

    vector<int> mid;
    while( current.size() != 0 ){
        mid = current;// intermediate tree node
        current.clear();
        for ( int i = 0; i < mid.size(); i++ ){
            for ( int j = 0; j < GTree[mid[i]].children.size(); j++ ){
                current.push_back( GTree[mid[i]].children[j] );
            }
        }
        if ( current.size() == 0 ) break;
        treenodelevel.push_back( current );
    }
    /// It seems that G-tree is an unbalanced tree
    // bottom up calculation
    // temp graph
    vector<Node> graph;
    graph = Nodes;//original graph
    vector<int> cands;
    vector<int> result;
    unordered_map<int, unordered_map<int,int> > vertex_pairs;//result of distance matrix

    // do dijkstra
    int s, t, tn, nid, cid, weight;
    vector<int> tnodes, tweight;
    set<int> nset;

    if(ifParallel){
        int threadnum = 0;
        vector<vector<int>> tnVector(threadnum);
        vector<vector<unordered_map<int, unordered_map<int,int> >>> vertex_pairsVV;//result of distance matrix
        for ( int i = treenodelevel.size() - 2; i >= 0; i-- ){//start from the lowest level
            boost::thread_group threadf;
            if(treenodelevel[i].size() > thread_num){
                threadnum = thread_num;
            }else{
                threadnum = treenodelevel[i].size();
            }
            tnVector.assign(threadnum,vector<int>());
            vertex_pairsVV.assign(threadnum,vector<unordered_map<int, unordered_map<int,int> >>());
            int thread_i = 0;
            for ( int j = 0; j < treenodelevel[i].size(); j++ ) {//for each partition in this level
                tn = treenodelevel[i][j];
                tnVector[thread_i].push_back(tn);
                if(thread_i >= threadnum-1){
                    thread_i = 0;
                }else{
                    thread_i++;
                }
            }
            for(int ti=0;ti<threadnum;++ti){
                vertex_pairsVV[ti].assign(tnVector[ti].size(),unordered_map<int, unordered_map<int,int> >());
                /// multi-thread
                threadf.add_thread(new boost::thread(&Gtree::hierarchy_shortest_path_compute_update, this, i, tnVector[ti], boost::ref(graph), boost::ref(vertex_pairsVV), ti));
            }
            threadf.join_all();
            for(int ti=0;ti<threadnum;++ti) {
                int ii=0;
                for (auto it = tnVector[ti].begin(); it != tnVector[ti].end(); ++it, ++ii) {
                    vertex_pairs = vertex_pairsVV[ti][ii];
                    int tn = *it;
                    // IMPORTANT! after all border finished, degenerate graph
                    // first, remove inward edges
                    for (int k = 0; k < GTree[tn].borders.size(); k++) {//for each border vertex
                        s = GTree[tn].borders[k];
                        tnodes.clear();
                        tweight.clear();
                        for (int p = 0; p < graph[s].adjnodes.size(); p++) {//for each adjacent vertex
                            nid = graph[s].adjnodes[p];
                            weight = graph[s].adjweight[p];
                            // if adj node in same tree node

                            if (graph[nid].gtreepath.size() <= i || graph[nid].gtreepath[i] !=
                                                                    tn) {// add the higher-level nodes or other node in the same level
                                // only remain those useful vertices, i.e., borders
                                tnodes.push_back(nid);
                                tweight.push_back(weight);
                            }
                        }
                        // cut it
                        graph[s].adjnodes = tnodes;//update the adjacency lists of graph, only left the useful boundary vertices
                        graph[s].adjweight = tweight;
                    }
                    // second, add inter connected edges (shortcuts)
                    for (int k = 0; k < GTree[tn].borders.size(); k++) {
                        for (int p = 0; p < GTree[tn].borders.size(); p++) {
                            if (k == p) continue;
                            s = GTree[tn].borders[k];
                            t = GTree[tn].borders[p];
                            graph[s].adjnodes.push_back(t);
                            graph[s].adjweight.push_back(vertex_pairs[s][t]);
                        }
                    }
                }
            }

        }
    }
    else{
        for ( int i = treenodelevel.size() - 2; i >= 0; i-- ){//start from the lowest level
            for ( int j = 0; j < treenodelevel[i].size(); j++ ){//for each partition in this level
                tn = treenodelevel[i][j];

                cands.clear();
                if ( GTree[tn].isleaf ){//for leaf node
                    // cands = leafnodes
                    cands = GTree[tn].leafnodes;//for leaf node, the cands is the leafnodes
                    // union borders = borders;
//                    GTree[tn].union_borders = GTree[tn].borders;
                }
                else{//for non-leaf node
//                    nset.clear();
//                    for ( int k = 0; k < GTree[tn].children.size(); k++ ){
//                        cid = GTree[tn].children[k];
//                        nset.insert( GTree[cid].borders.begin(), GTree[cid].borders.end() );//get the borders of all children
//                    }
                    // union borders = cands;

                    cands.clear();
//                    for ( set<int>::iterator it = nset.begin(); it != nset.end(); it ++ ){//push each vertex of this partition into cands
//                        cands.push_back( *it );//for non-leaf node, the cands is the borders of all child nodes
//                    }
//                    GTree[tn].union_borders = cands;//for non-leaf node, the union_borders contains all the borders of children
                    cands = GTree[tn].union_borders;
                }

                // start to do min dis
                vertex_pairs.clear();

                // for each border, do min dis
                int cc = 0;

                for ( int k = 0; k < GTree[tn].union_borders.size(); k++ ){
                    //printf("DIJKSTRA...LEAF=%d BORDER=%d\n", tn, GTree[tn].union_borders[k] );
                    result = dijkstra_candidate( GTree[tn].union_borders[k], cands, graph );//distance vector from s to all borders
                    //printf("DIJKSTRA...END\n");

                    // save to map
                    for ( int p = 0; p < result.size(); p ++ ){
                        GTree[tn].mind[k*cands.size()+p] =  result[p];//store the distance vector of s
                        vertex_pairs[GTree[tn].union_borders[k]][cands[p]] = result[p];
                    }
                }

                // IMPORTANT! after all border finished, degenerate graph
                // first, remove inward edges
                for ( int k = 0; k < GTree[tn].borders.size(); k++ ){//for each border vertex
                    s = GTree[tn].borders[k];
                    tnodes.clear();
                    tweight.clear();
                    for ( int p = 0; p < graph[s].adjnodes.size(); p++ ){//for each adjacent vertex
                        nid = graph[s].adjnodes[p];
                        weight = graph[s].adjweight[p];
                        // if adj node in same tree node

                        if ( graph[nid].gtreepath.size() <= i || graph[nid].gtreepath[i] != tn ){// add the higher-level nodes or other node in the same level
                            // only remain those useful vertices, i.e., borders
                            tnodes.push_back(nid);
                            tweight.push_back(weight);
                        }
                    }
                    // cut it
                    graph[s].adjnodes = tnodes;//update the adjacency lists of graph, only left the useful boundary vertices
                    graph[s].adjweight = tweight;
                }
                // second, add inter connected edges (shortcuts)
                for ( int k = 0; k < GTree[tn].borders.size(); k++ ){
                    for ( int p = 0; p < GTree[tn].borders.size(); p++ ){
                        if ( k == p ) continue;
                        s = GTree[tn].borders[k];
                        t = GTree[tn].borders[p];
                        graph[s].adjnodes.push_back( t );
                        graph[s].adjweight.push_back( vertex_pairs[s][t] );
                    }
                }
            }
        }
    }


}

void Gtree::hierarchy_shortest_path_compute_update(int i, vector<int> tn_v, vector<Node> & graph, vector<vector<unordered_map<int, unordered_map<int,int> >>> & vertex_pairsVV, int thread_i){
    vector<int> cands;
    vector<int> result;
    unordered_map<int, unordered_map<int,int> > vertex_pairs;//result of distance matrix

    // do dijkstra
    int s, t, nid, cid, weight;
    vector<int> tnodes, tweight;
    set<int> nset;
    int ii=0;
    for(auto it=tn_v.begin();it!=tn_v.end();++it,++ii){
        int tn = *it;
        cands.clear();
        if ( GTree[tn].isleaf ){//for leaf node
            // cands = leafnodes
            cands = GTree[tn].leafnodes;//for leaf node, the cands is the leafnodes
            // union borders = borders;
            GTree[tn].union_borders = GTree[tn].borders;
        }
        else{//for non-leaf node
            cands = GTree[tn].union_borders;//for non-leaf node, the union_borders contains all the borders of children
        }

        // start to do min dis
        vertex_pairs.clear();

        // for each border, do min dis
        int cc = 0;

        for ( int k = 0; k < GTree[tn].union_borders.size(); k++ ){
            //printf("DIJKSTRA...LEAF=%d BORDER=%d\n", tn, GTree[tn].union_borders[k] );
            result = dijkstra_candidate( GTree[tn].union_borders[k], cands, graph );//distance vector from s to all borders
            //printf("DIJKSTRA...END\n");

            // save to map
            for ( int p = 0; p < result.size(); p++ ){
                GTree[tn].mind[k*cands.size() + p] = result[p];//store the distance vector of s
                vertex_pairs[GTree[tn].union_borders[k]][cands[p]] = result[p];
            }
        }

        vertex_pairsVV[thread_i][ii] = vertex_pairs;

//        // IMPORTANT! after all border finished, degenerate graph
//        // first, remove inward edges
//        for ( int k = 0; k < GTree[tn].borders.size(); k++ ){//for each border vertex
//            s = GTree[tn].borders[k];
//            tnodes.clear();
//            tweight.clear();
//            for ( int p = 0; p < graph[s].adjnodes.size(); p++ ){//for each adjacent vertex
//                nid = graph[s].adjnodes[p];
//                weight = graph[s].adjweight[p];
//                // if adj node in same tree node
//
//                if ( graph[nid].gtreepath.size() <= i || graph[nid].gtreepath[i] != tn ){// add the higher-level nodes or other node in the same level
//                    // only remain those useful vertices, i.e., borders
//                    tnodes.push_back(nid);
//                    tweight.push_back(weight);
//                }
//            }
//            // cut it
//            graph[s].adjnodes = tnodes;//update the adjacency lists of graph, only left the useful boundary vertices
//            graph[s].adjweight = tweight;
//        }
//        // second, add inter connected edges (shortcuts)
//        for ( int k = 0; k < GTree[tn].borders.size(); k++ ){
//            for ( int p = 0; p < GTree[tn].borders.size(); p++ ){
//                if ( k == p ) continue;
//                s = GTree[tn].borders[k];
//                t = GTree[tn].borders[p];
//                graph[s].adjnodes.push_back( t );
//                graph[s].adjweight.push_back( vertex_pairs[s][t] );
//            }
//        }
    }

}

void Gtree::hierarchy_shortest_path_compute(int i, vector<int> & tn_v, vector<Node> & graph, vector<vector<unordered_map<int, unordered_map<int,int> >>> & vertex_pairsVV, int thread_i){
    vector<int> cands;
    vector<int> result;
    unordered_map<int, unordered_map<int,int> > vertex_pairs;//result of distance matrix

    // do dijkstra
    int s, t, nid, cid, weight;
    vector<int> tnodes, tweight;
    set<int> nset;
    int ii=0;
    for(auto it=tn_v.begin();it!=tn_v.end();++it,++ii){
        int tn = *it;
        cands.clear();
        if ( GTree[tn].isleaf ){//for leaf node
            // cands = leafnodes
            cands = GTree[tn].leafnodes;//for leaf node, the cands is the leafnodes
            // union borders = borders;
            GTree[tn].union_borders = GTree[tn].borders;
        }
        else{//for non-leaf node
            nset.clear();
            for ( int k = 0; k < GTree[tn].children.size(); k++ ){
                cid = GTree[tn].children[k];
                nset.insert( GTree[cid].borders.begin(), GTree[cid].borders.end() );//get the borders of all children
            }
            // union borders = cands;

            cands.clear();
            for ( set<int>::iterator it = nset.begin(); it != nset.end(); it ++ ){//push each vertex of this partition into cands
                cands.push_back( *it );//for non-leaf node, the cands is the borders of all child nodes
            }
            GTree[tn].union_borders = cands;//for non-leaf node, the union_borders contains all the borders of children
        }

        // start to do min dis
        vertex_pairs.clear();

        // for each border, do min dis
        int cc = 0;

        for ( int k = 0; k < GTree[tn].union_borders.size(); k++ ){
            //printf("DIJKSTRA...LEAF=%d BORDER=%d\n", tn, GTree[tn].union_borders[k] );
            result = dijkstra_candidate( GTree[tn].union_borders[k], cands, graph );//distance vector from s to all borders
            //printf("DIJKSTRA...END\n");

            // save to map
            for ( int p = 0; p < result.size(); p ++ ){
                GTree[tn].mind.push_back( result[p] );//store the distance vector of s
                vertex_pairs[GTree[tn].union_borders[k]][cands[p]] = result[p];
            }
        }

        vertex_pairsVV[thread_i][ii] = vertex_pairs;
//        // IMPORTANT! after all border finished, degenerate graph
//        // first, remove inward edges
//        for ( int k = 0; k < GTree[tn].borders.size(); k++ ){//for each border vertex
//            s = GTree[tn].borders[k];
//            tnodes.clear();
//            tweight.clear();
//            for ( int p = 0; p < graph[s].adjnodes.size(); p++ ){//for each adjacent vertex
//                nid = graph[s].adjnodes[p];
//                weight = graph[s].adjweight[p];
//                // if adj node in same tree node
//
//                if ( graph[nid].gtreepath.size() <= i || graph[nid].gtreepath[i] != tn ){// add the higher-level nodes or other node in the same level
//                    // only remain those useful vertices, i.e., borders
//                    tnodes.push_back(nid);
//                    tweight.push_back(weight);
//                }
//            }
//            // cut it
//            graph[s].adjnodes = tnodes;//update the adjacency lists of graph, only left the useful boundary vertices
//            graph[s].adjweight = tweight;
//        }
//        // second, add inter connected edges (shortcuts)
//        for ( int k = 0; k < GTree[tn].borders.size(); k++ ){
//            for ( int p = 0; p < GTree[tn].borders.size(); p++ ){
//                if ( k == p ) continue;
//                s = GTree[tn].borders[k];
//                t = GTree[tn].borders[p];
//                graph[s].adjnodes.push_back( t );
//                graph[s].adjweight.push_back( vertex_pairs[s][t] );
//            }
//        }
    }

}

// dump distance matrix into file
void Gtree::hierarchy_shortest_path_save(){
    char filename[200];
    if(dataset == "cal"){
        strcpy(filename, DataPath.c_str());
        strcat(filename, FILE_ONTREE_MIND.c_str());
    }else {
        strcpy(filename, dirname.c_str());
        strcat(filename, "/");
        strcat(filename, FILE_ONTREE_MIND.c_str());
    }
    FILE* fout = fopen( filename, "wb" );
    if(fout == NULL){
        cout<<"Failed to open file "<<filename<<endl;
        exit(1);
    }
    int* buf;
    int count;
    for ( int i = 0; i < GTree.size(); i++ ){
        // union borders
        count = GTree[i].union_borders.size();
        fwrite( &count, sizeof(int), 1, fout );
        buf = new int[count];
        copy( GTree[i].union_borders.begin(), GTree[i].union_borders.end(), buf );
        fwrite( buf, sizeof(int), count, fout );
        delete[] buf;
        // mind
        count = GTree[i].mind.size();
        fwrite( &count, sizeof(int), 1, fout );
        buf = new int[count];
        copy( GTree[i].mind.begin(), GTree[i].mind.end(), buf );
        fwrite( buf, sizeof(int), count, fout );
        delete[] buf;
    }
    fclose(fout);
}

// load distance matrix from file
void Gtree::load_minds(){
    char filename[200];
    if(dataset == "cal"){
        strcpy(filename, DataPath.c_str());
        strcat(filename, FILE_ONTREE_MIND.c_str());
    }else {
        strcpy(filename, dirname.c_str());
        strcat(filename, "/");
        strcat(filename, FILE_ONTREE_MIND.c_str());
    }
    FILE* fin = fopen( filename, "rb" );
//	FILE* fin = fopen( FILE_ONTREE_MIND, "rb" );
    if(fin == NULL){
        cout << "Failed to open file " << filename << endl;
        exit(1);
    }
    cout << "Loading distance matrix... ";
    int* buf;
    int count, pos = 0;
    while( fread( &count, sizeof(int), 1, fin ) ){
        // union borders
        buf = new int[count];
        fread( buf, sizeof(int), count, fin );
        GTree[pos].union_borders.clear();
        for ( int i = 0; i < count; i++ ){
            GTree[pos].union_borders.push_back(buf[i]);
        }
        delete[] buf;
        // mind
        fread( &count, sizeof(int), 1, fin );
        buf = new int[count];
        fread( buf, sizeof(int), count, fin );
        GTree[pos].mind.clear();
        for ( int i = 0; i < count; i++ ){
            GTree[pos].mind.push_back(buf[i]);
        }
        if(pos==5061){
//            int i=0;
//            int span = GTree[pos].leafnodes.size();
//            for(auto it=GTree[pos].mind.begin();it!=GTree[pos].mind.end();++it){
//                if(i%span==0)
//                    cout<<endl;
//                cout<<*it<<" ";
//                ++i;
//            }
//            cout<<endl;
//            cout<<GTree[pos].isleaf<<endl;
//            cout<<GTree[pos].mind.size()<<endl;
//            int leafposID1=4,leafposID2=0;
//            cout<<GTree[pos].mind[leafposID2*span + leafposID1]<<" "<<GTree[pos].mind[leafposID1*span + leafposID2]<<" (old); "<<endl;
        }
        pos++;
        delete[] buf;
    }
    fclose(fin);
    cout << "Done." << endl;
}

//build G-Tree
int Gtree::gtree_build(bool ifParallel){
    // init
//    TIME_TICK_START
    init();
//    TIME_TICK_END
//    TIME_TICK_PRINT("INIT")
    double t1,t2;

    // gtree_build
    Timer tt;
    tt.start();
    cout << "Start to build G-tree..."<<endl;
//    TIME_TICK_START
    build();
//    TIME_TICK_END
//    TIME_TICK_PRINT("BUILD")
//    cout <<"Done."<<endl;
    tt.stop();
    t1 = tt.GetRuntime();
    cout << "The time for G-tree building: " << t1 << " s." << endl;

    // dump gtree
//    cout << "Saving G-tree..."<<endl;
    gtree_save();
//    cout <<"Done."<<endl;

    // calculate distance matrix
    cout << "Start to calculate distance matrix..."<<endl;
//    tt.start();
//    TIME_TICK_START
    hierarchy_shortest_path_calculation(ifParallel);
//    TIME_TICK_END
//    TIME_TICK_PRINT("MIND")
//    cout << "Done."<<endl;
    tt.stop();
    t2 = tt.GetRuntime();
    cout << "The time for distance matrix computing: " << t2-t1 << " s." << endl;
    cout << "Overall time for index construction: " << t2 << " s." << endl;
    // dump distance matrix
//    cout << "Saving distance matrix..."<<endl;
    hierarchy_shortest_path_save();
//    cout << "Done."<<endl;

    return 0;
}

//Function for generating random queries
void Gtree::ODpairGenerate(int times){
    string RandomDis = string(DataPath) + "/" + dataset + "/" + FILE_QUERY;

    if(node_num == 0){
        if(dataset == "cal"){
            RandomDis = string(DataPath) + FILE_QUERY;
            load_graph();
        }else{
            string r_file = string(DataPath)  + "/" +dataset + "/" + dataset;
//    string r_file = string(DataPath)  + dataset + "/" + dataset + "_Processed";

            ifstream inFile(r_file, ios::in);
            if (!inFile) { // if not exist
                cout << "Fail to open file" << r_file << endl;
                exit(1);
            }

            /// read graph and recording the degree of vertices
            inFile >> node_num >> edge_num;
            inFile.close();
        }
    }

    /*---OD pairs generation---*/
    int pairs = 0;
    int node_start, node_end;
    int temp = 0;
    ofstream outFile(RandomDis, ios::out);
    if (!outFile) {
        cout << "File opening failed." << endl;
        assert(outFile);
    }
    cout << "Query OD pairs file generating..." << endl;
    outFile << times << endl;
    //generate random OD pairs
    pairs = 0;
    while (pairs < times) {
        node_start = rand() % node_num;
        node_end = rand() % node_num;
        while(node_end == node_start){
            node_end = rand() % node_num;
        }
//        outFile << node_start+1 << ' ' << node_end+1 << endl;
        outFile << node_start << ' ' << node_end << endl;
        ++pairs;
    }

    outFile.close();
    cout << "Finished." << endl;
}

//Function for generating update edges
void Gtree::UpdateGenerate(int times){
    string filename = string(DataPath) + "/" + dataset + "/" + FILE_UPDATE;
    /// read graph
    string r_file = string(DataPath)  + "/" + dataset + "/" + dataset;
//    string r_file = string(DataPath)  + dataset + "/" + dataset + "_Processed";
    ifstream inFile(r_file, ios::in);
    if (!inFile) { // if not exist
        cout << "Fail to open file" << r_file << endl;
        exit(1);
    }
    inFile >> node_num >> edge_num;
    int ID1,ID2,weight;
    vector<pair<pair<int,int>,int>> edges;
    while(inFile){
        inFile >> ID1 >> ID2 >> weight;
        if(ID1 < ID2){
            edges.emplace_back(make_pair(ID1,ID2),weight);
        }
    }
    inFile.close();

    int size = edges.size();
    assert(2*size == edge_num);
    /*---OD pairs generation---*/
    int pairs = 0;
    int id = 0;

    ofstream outFile(filename, ios::out);
    if (!outFile) {
        cout << "File opening failed." << endl;
        assert(outFile);
    }
    cout << "Update OD pairs file generating..." << endl;
    outFile << times << endl;
    //generate random OD pairs
    pairs = 0;
    unordered_set<int> edgeIdSet;
    edgeIdSet.clear();
    while (pairs < times) {
        id = rand() % size;
        if(edges[id].second >= 2 && edgeIdSet.find(id)==edgeIdSet.end()){//if edge weight is no smaller than 2, and it has not been added.
            edgeIdSet.insert(id);
            outFile << edges[id].first.first << ' ' << edges[id].first.second << ' '<<edges[id].second<< endl;
        }else{
            continue;
        }

        ++pairs;
    }
    outFile.close();
    cout << "Finished." << endl;
}

//Function of in-memory Dijkstra's algorithm
Distance Gtree::Dijkstra(NodeId s, NodeId t, vector<Node> & Nodes) { // second version, powered by benchmark::heap
    if (s == t) return 0;
    benchmark::heap<2, NodeId, Distance> pqueue(node_num);
    int item_id, temp_id;
    Distance item_dis, temp_dis;
    vector<bool> closed(node_num, false); //flag vector of whether closed
    vector<Distance> cost(node_num, INF);   //vector of cost
//    vector<NodeId> pre(node_num, -1);       //vector of predecessor id
    Distance min_cost = INF;
    //Initiation of start node
    cost[s] = 0;//cost of start node
    pqueue.update(s, 0);

    //Iteration
    while (!pqueue.empty()) {//for every node in pqueue
        pqueue.extract_min(item_id, item_dis);// top and delete min item
//        cout<<item_id<<" ";
        if (item_id == t) {//if reach target node
            min_cost = cost[item_id];
//            cout<<"Minimal cost: " << min_cost<<endl;
            break;
        }
        //relaxation
//        cout <<Nodes[item_id].adjnodes.size()<<endl;
        for (int i = 0; i < Nodes[item_id].adjnodes.size(); ++i) {
            temp_id = Nodes[item_id].adjnodes[i];
            if (closed[temp_id])//if closed
                continue;
//            cout << item_id << " " << temp_id << endl;
            assert(Nodes[item_id].adjnodes.size() == Nodes[item_id].adjweight.size());
            temp_dis = item_dis + Nodes[item_id].adjweight[i];
            if (cost[temp_id] > temp_dis) {//slack operation
                cost[temp_id] = temp_dis;
//                pre[temp_id] = item_id;
                pqueue.update(temp_id, temp_dis);
            }
        }
        closed[item_id] = true;
    }
//    cout<<endl;
//    Dij_getPath(pre,s,t);//traverse the pre vector to get the shortest path
    return min_cost;
}
//function of reading update edges
void Gtree::ReadUpdates(string filename){
    int ID1, ID2, weight;
    ifstream inFile(filename, ios::in);
    if (!inFile) { // if not exist
        cout << "Fail to open file" << filename << endl;
        exit(1);
    }
    int num;
    inFile >> num;
    for(int i=0;i<num;i++){
        inFile>>ID1>>ID2>>weight;
        updateEdges.emplace_back(make_pair(ID1, ID2), weight);
    }
    inFile.close();
}
//function of checking if nid is an ancestor of ID2
bool Gtree::CheckIfAncestor(int ID2, int nid){
    vector<int> fathers = Nodes[ID2].gtreepath;
    for(int i=fathers.size()-1;i>=0;--i){
        if(fathers[i] == nid){
            return true;
        }
    }
    return false;
}
//function of initiation of update
void Gtree::init_update() {
    if(task_type == 3){
        ReadGraph_W();
    }
    load_gtreeQ();
    load_minds();
    build_up_and_down_pos();//
    init_borders();//set border flag
}
//function of getting isborder flag
void Gtree::init_borders(){
    int leaf_num=0;
    for(auto it=GTree.begin();it!=GTree.end();++it){
        if(it->isleaf){
            for(int i=0;i<it->borders.size();++i){
                Nodes[it->borders[i]].isborder = true;
            }
            ++leaf_num;
        }
    }
    cout<<"The leaf node number is: "<<leaf_num<<" ; The non-leaf node number is: "<<GTree.size()-leaf_num<<endl;
}
//function of initiating query processing
void Gtree::init_query(vector<TreeNode> & GTree) {
    for (auto &tn: GTree) {
        //tn.oclist.clear();
        tn.cache_q = vector<int>(tn.borders.size(), 0);
        tn.is_visited = false;
    }
}
//function of finding the LCA position in gtreepath
inline int Gtree::find_LCA_pos(int src, int dst) {
    for (int i = 1; i < Nodes[src].deep && i < Nodes[dst].deep; ++i) {
        if (Nodes[src].gtreepath[i] != Nodes[dst].gtreepath[i])
            return i - 1;
    }
    return 0;
}
// heap-based dijkstra search for updating the leaf nodes
bool Gtree::dijkstra_candidate_update(int ID1, int ID1_pos, vector<int> &cands, vector<Node> &graph){//vector<pair<pair<int,int>,pair<int,int>>> & borderShortcuts
    // init
    bool borderUpdate = false;
    set<int> todo;
    todo.clear();
    todo.insert(cands.begin(), cands.end());//get the partition vertex set

    map<int,int> result;
    result.clear();
    unordered_set<int> visited;
    visited.clear();
//    unordered_map<int,int> q;//map used for mapping vertex id to its distance from source vertex
    benchmark::heap<2, int, int> q(node_num);
//    vector<bool> closed(node_num, false); //flag vector of whether closed
    vector<Distance> cost(node_num, INF);   //vector of cost
    int temp_dis;
    q.clear();
    q.update(ID1,0);
    cost[ID1] = 0;
    // start
    int min, minpos, ID2, weight;
    while( ! todo.empty() && ! q.empty() ){
        min = -1;
        q.extract_min(minpos,min);
        // put min to result, add to visited
        result[minpos] = min;
        visited.insert( minpos );
//        q.erase(minpos);
        if ( todo.find( minpos ) != todo.end() ){//if found, erase visited vertex
            todo.erase( minpos );
        }
        // expand on graph (the original graph)
        for ( int i = 0; i < graph[minpos].adjnodes.size(); i++ ){
            ID2 = graph[minpos].adjnodes[i];
            if ( visited.find( ID2 ) != visited.end() ){//if found, ie, it is visited
                continue;
            }
            weight = graph[minpos].adjweight[i];//edge weight
            temp_dis = min + weight;
            if ( temp_dis < cost[ID2] ){
                q.update(ID2,temp_dis);
                cost[ID2] = temp_dis;
            }
        }
    }
    /// update
    int nid, ID2_pos;
    int oldW;
    nid = graph[ID1].inleaf;
    assert(GTree[nid].isleaf);
    int span = GTree[nid].leafnodes.size();
    for ( int i = 0; i < cands.size(); i++ ){
        ID2 = cands[i];
        ID2_pos = graph[ID2].inleafpos;
        assert(ID2_pos == i);
        if(result[ID2] != GTree[nid].mind[ID1_pos*span + ID2_pos]){//if there is update
            if(ifDebug){
                oldW = GTree[nid].mind[ID1_pos*span + ID2_pos];
                int temp_dis = Dijkstra(ID1,ID2,NodesO);
                if(temp_dis != oldW)
                    cout<<"Leaf node old incorrect: "<<ID1<<" "<<ID2<<": "<<oldW<<" "<<temp_dis<<endl;
                temp_dis = Dijkstra(ID1,ID2,graph);
                if(temp_dis != result[ID2])
                    cout<<"Leaf node new incorrect: "<<ID1<<" "<<ID2<<": "<<result[ID2]<<" "<<temp_dis<<endl;
//            cout<<ID1<<" "<<ID2<<": "<<oldW<<" (old); "<<result[ID2]<<" (new)."<<endl;
            }
            GTree[nid].mind[ID1_pos*span + ID2_pos] = result[ID2];
            if(graph[ID2].isborder){//if it is border
//                borderShortcuts.emplace_back(make_pair(ID1,ID2), make_pair(oldW,result[ID2]));//add
                borderUpdate = true;
//                    cout<<"!!! "<<ID1<<" "<<ID2<<endl;
            }
        }
    }
    // return
    return borderUpdate;//vector of distance matrix values
}
// heap-based dijkstra search for updating the non-leaf nodes
bool Gtree::dijkstra_candidate_update(int ID1, int ID1_pos, int nid, vector<int> &cands, vector<Node> &graph, unordered_set<int> & borderSet){// vector<pair<pair<int,int>,pair<int,int>>> & borderShortcuts,
    // init
    bool borderUpdate = false;
    set<int> todo;
    todo.clear();
    todo.insert(cands.begin(), cands.end());//get the partition vertex set

    map<int,int> result;
    result.clear();
    unordered_set<int> visited;
    visited.clear();
//    unordered_map<int,int> q;//map used for mapping vertex id to its distance from source vertex
    benchmark::heap<2, int, int> q(node_num);
//    vector<bool> closed(node_num, false); //flag vector of whether closed
    vector<Distance> cost(node_num, INF);   //vector of cost
    int temp_dis;
    q.clear();
    q.update(ID1,0);
    cost[ID1] = 0;
    // start
    int min, minpos, ID2, weight;
    while( ! todo.empty() && ! q.empty() ){
        min = -1;
        q.extract_min(minpos,min);
        // put min to result, add to visited
        result[minpos] = min;
        visited.insert( minpos );
//        q.erase(minpos);
        if ( todo.find( minpos ) != todo.end() ){//if found, erase visited vertex
            todo.erase( minpos );
        }
        // expand on graph (the original graph)
        for ( int i = 0; i < graph[minpos].adjnodes.size(); i++ ){
            ID2 = graph[minpos].adjnodes[i];
            if ( visited.find( ID2 ) != visited.end() ){//if found, ie, it is visited
                continue;
            }
            weight = graph[minpos].adjweight[i];//edge weight
            temp_dis = min + weight;
            if ( temp_dis < cost[ID2] ){
                q.update(ID2,temp_dis);
                cost[ID2] = temp_dis;
            }
        }
    }
    /// update
    int ID2_pos;
    int oldW;
    bool flag_border = false;
    if(borderSet.find(ID1)!=borderSet.end()){
        flag_border = true;
    }
    int span = GTree[nid].union_borders.size();
//    if(span!=cands.size())
//        cout<<"!!!"<<endl;
    for ( int i = 0; i < cands.size(); i++ ){
        ID2 = cands[i];
        ID2_pos = ID1_pos+1+i;
//        if(ID1==ID2)
//            continue;
        if(result[ID2] != GTree[nid].mind[ID1_pos*span + ID2_pos]){//if there is update
            if(ifDebug){
                oldW = GTree[nid].mind[ID1_pos*span + ID2_pos];
                int temp_dis = Dijkstra(ID1,ID2,NodesO);
                if(temp_dis != oldW)
                    cout<<"Non leaf node old incorrect: "<<ID1<<" "<<ID2<<": "<<oldW<<" "<<temp_dis<<endl;
                temp_dis = Dijkstra(ID1,ID2,graph);
                if(temp_dis != result[ID2])
                    cout<<"Non leaf node new incorrect: "<<ID1<<" "<<ID2<<": "<<result[ID2]<<" "<<temp_dis<<endl;
//                cout<<ID1<<" "<<ID2<<": "<<oldW<<" (old); "<<result[ID2]<<" (new)."<<endl;
                if(GTree[nid].mind[ID1_pos*span + ID2_pos] != GTree[nid].mind[ID2_pos*span + ID1_pos])
                    cout<<"!!! Original Distance matrix incorrect"<<endl;
            }

            GTree[nid].mind[ID1_pos*span + ID2_pos] = result[ID2];
            GTree[nid].mind[ID2_pos*span + ID1_pos] = result[ID2];
            if(flag_border){// if ID1 is border
                if(borderSet.find(ID2)!=borderSet.end()){//if ID2 is border too
//                    borderShortcuts.emplace_back(make_pair(ID1,ID2), make_pair(oldW,result[ID2]));//add
                    borderUpdate = true;
//                    cout<<"!!! "<<ID1<<" "<<ID2<<endl;
                }
            }
        }
    }
    // return
    return borderUpdate;//vector of distance matrix values
}
//function of contracting leaf node by adding shortcuts and removing inner nodes
void Gtree::LeafNodeContract(int lnID, vector<Node> & graph, unordered_set<int> & leafsToCheck, unordered_set<int> & updatedLeafs){//  unordered_set<int> & updatedLeafs,
    /// deal with current leaf node
    // IMPORTANT! after all border finished, degenerate graph
    // first, remove inward edges
    int s, t, nid, weight;
    vector<int> tnodes, tweight;
    set<int> nset;
    for ( int k = 0; k < GTree[lnID].borders.size(); k++ ){//for each border vertex
        s = GTree[lnID].borders[k];
        tnodes.clear();
        tweight.clear();
        for ( int p = 0; p < graph[s].adjnodes.size(); p++ ){//for each adjacent vertex
            nid = graph[s].adjnodes[p];//adjacent vertex id
            weight = graph[s].adjweight[p];
            // if adj node is border
            if ( graph[nid].inleaf != lnID && graph[nid].isborder){//
                // only remain those useful vertices, i.e., borders
                tnodes.push_back(nid);
                tweight.push_back(weight);
                int temp_id = graph[nid].inleaf;
                if(temp_id != lnID){
                    if(leafsToCheck.find(temp_id)==leafsToCheck.end() && updatedLeafs.find(temp_id) == updatedLeafs.end()){// if not found in updatedLeafs and leafsToCheck
                        leafsToCheck.insert(temp_id);
                    }
                }
            }
        }
        // cut it
        graph[s].adjnodes = tnodes;//update the adjacency lists of graph, only left the useful boundary vertices
        graph[s].adjweight = tweight;
    }
    // second, add inter connected edges (shortcuts)
    int nleafnode = GTree[lnID].leafnodes.size();
    int ID2_pos;
    for ( int k = 0; k < GTree[lnID].borders.size(); k++ ){
        for ( int p = 0; p < GTree[lnID].borders.size(); p++ ){
            if ( k == p ) continue;
            s = GTree[lnID].borders[k];
            t = GTree[lnID].borders[p];
            ID2_pos = graph[t].inleafpos;
            graph[s].adjnodes.push_back( t );
            graph[s].adjweight.push_back( GTree[lnID].mind[k*nleafnode + ID2_pos]);
            if(ifDebug){
                int temp_dis = Dijkstra(s,t,Nodes);
                if(temp_dis != GTree[lnID].mind[k*nleafnode + ID2_pos]){
                    cout<<"Leaf node contraction error: "<<s<<" "<<t<<GTree[lnID].mind[k*nleafnode + ID2_pos]<<" "<<temp_dis<<endl;
                }
            }
        }
    }
}
//function of contracting node by adding shortcuts and removing inner nodes
void Gtree::NonLeafNodeContract(int tn, int level_i, vector<Node> & graph, unordered_set<int> & nodesToCheck, unordered_set<int> & updatedNodes, vector< vector<int> > & treenodelevel, unordered_set<int> & borderSet){//  unordered_set<int> & updatedLeafs,
    // IMPORTANT! after all border finished, degenerate graph
    // first, remove inward edges
    vector<int> cands;
    vector<int> result;
//    unordered_map<int, unordered_map<int,int> > vertex_pairs;//result of distance matrix
    int s, t, nid, cid, weight;
    vector<int> tnodes, tweight;
    set<int> nset;
//    unordered_set<int> unionBorderSet; unionBorderSet.clear();
//    for(auto it=GTree[tn].union_borders.begin();it!=GTree[tn].union_borders.end();++it){
//        unionBorderSet.insert(*it);
//    }

    // IMPORTANT! after all border finished, degenerate graph
    // first, remove inward edges
    for ( int k = 0; k < GTree[tn].borders.size(); k++ ){//for each border vertex
        s = GTree[tn].borders[k];
        tnodes.clear();
        tweight.clear();
        for ( int p = 0; p < graph[s].adjnodes.size(); p++ ){//for each adjacent vertex
            t = graph[s].adjnodes[p];
            weight = graph[s].adjweight[p];

            tnodes.push_back(t);
            tweight.push_back(weight);
            // prune
            cid = graph[t].gtreepath[level_i];
            //if t is in the neighbor tree node of this level
            if(cid != tn){
                if(nodesToCheck.find(cid)==nodesToCheck.end() && updatedNodes.find(cid) == updatedNodes.end()){// if not found in updatedLeafs and leafsToCheck
                    nodesToCheck.insert(cid);
                }
            }
            /*if ( unionBorderSet.find(t)!=unionBorderSet.end() || cid != tn ){// add the higher-level nodes or other node in the same level
                // only remain those useful vertices, i.e., borders
//                tnodes.push_back(t);
//                tweight.push_back(weight);
                if(cid != tn){
                    if(nodesToCheck.find(cid)==nodesToCheck.end() && updatedNodes.find(cid) == updatedNodes.end()){// if not found in updatedLeafs and leafsToCheck
                        nodesToCheck.insert(cid);
                    }
                }
            }*/
        }
        // cut it
        graph[s].adjnodes = tnodes;//update the adjacency lists of graph, only left the useful boundary vertices
        graph[s].adjweight = tweight;
    }
    // second, add inter connected edges (shortcuts)
    int span = GTree[tn].union_borders.size();
    for ( int k = 0; k < GTree[tn].union_borders.size(); k++ ){
        s = GTree[tn].union_borders[k];
        if(borderSet.find(s) == borderSet.end())//if not found, i.e., it is not border
            continue;
        for ( int p = k+1; p < GTree[tn].union_borders.size(); p++ ){
            if ( k == p ) continue;
            t = GTree[tn].union_borders[p];
            if(borderSet.find(t) == borderSet.end())//if not found, i.e., it is not border
                continue;
            graph[s].adjnodes.push_back( t );
            graph[s].adjweight.push_back( GTree[tn].mind[k*span + p] );
            if(ifDebug){
                int temp_dis = Dijkstra(s,t,graph);
                if(temp_dis != GTree[tn].mind[k*span + p])
                    cout<<"Non-Leaf node contraction error:  "<<s<<" "<<t<<": "<<temp_dis<<" "<<GTree[tn].mind[k*span + p]<<endl;
            }
        }
    }
}
//function of dealing with batch edge increase update of G-tree
/*void Gtree::GtreeUpdate(int updateType, int updateVolume,bool ifPrue){//update(<ID1,ID2>,oldW) vector<pair<pair<int,int>,int> > & updates,
    // init
    init_update();
    // read updates
    string file = string(DataPath) + dataset + "/" + FILE_UPDATE;
    ReadUpdates(file);

//    cout<<"Processing updates... Volume of update: "<<updateVolume<<endl;
    //duplicate G-tree
    GTreeO = GTree; //old GTree
    NodesO = Nodes;//old Nodes

    int ID1,ID2,oldW,newW;
    vector<pair<pair<int,int>,pair<int,int> > > uWithinLeafNode;// edge update within the same leaf node (<ID1,ID2>,<oldW,newW>)
    vector<pair<pair<int,int>,pair<int,int> > > uCrossLeafNode;// edge update among borders (<ID1,ID2>,<oldW,newW>)
    vector<pair<pair<int,int>,pair<int,int> > > uAmongBorders;// edge/shortcut update among borders (<ID1,ID2>,<oldW,newW>)
    bool flag_borderUpdate = false;//whether there is border update of leaf node
    Timer tt;
    tt.start();
    ///preprocess: update edges
    assert(updateVolume<=updateEdges.size());
    for(int i=0;i<updateVolume;++i){//for each edge update
        ID1 = updateEdges[i].first.first;
        ID2 = updateEdges[i].first.second;
        oldW = updateEdges[i].second;
        switch (updateType) {
            case INCREASE:
                newW = 1.5*oldW;//update
                break;
            case DECREASE:
                newW = 0.5*oldW;
                break;
            case MIX:
                newW = 2*oldW * (rand()/double(RAND_MAX));
                break;
            default:
                cout<<"Wrong update type!"<<endl; break;
        }
        //update edge weight
        for(int i=0;i<Nodes[ID1].adjnodes.size();i++){
            if(Nodes[ID1].adjnodes[i]==ID2){
                assert(Nodes[ID1].adjweight[i] == oldW);
                Nodes[ID1].adjweight[i]=newW;
                break;
            }
        }
        for(int i=0;i<Nodes[ID2].adjnodes.size();i++){
            if(Nodes[ID2].adjnodes[i]==ID1){
                assert(Nodes[ID2].adjweight[i] == oldW);
                Nodes[ID2].adjweight[i]=newW;
                break;
            }
        }
        //identify the edge type
        if(Nodes[ID1].inleaf == Nodes[ID2].inleaf){
            uWithinLeafNode.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
        }else{
            uCrossLeafNode.emplace_back(make_pair(ID1,ID2),make_pair(oldW,newW));
            uAmongBorders.emplace_back(make_pair(ID1,ID2),make_pair(oldW,newW));
            flag_borderUpdate = true;
        }
    }
    /// update processing
    int lnID;
    int ID1_pos,ID2_pos;
    vector<Node> graph;//temp graph
    graph = Nodes;//original graph

    bool borderUpdate = false;
    int nleafnode;
    unordered_set<int> leafsToCheck;//the set of leaf nodes for update checking
    leafsToCheck.clear();
    unordered_set<int> updatedLeafs;
    updatedLeafs.clear();
    unordered_set<int> borderAffectedNodes;
    borderAffectedNodes.clear();
//    list<int> neighborNodes;
    vector<int> cands;
//    vector<pair<pair<int,int>,pair<int,int>>> aBorderShortcuts;//affected shortcuts among borders of leaf node
//    set<int> aLeafNodes;// set of affected leaf nodes

    int temp_pos;
    /// deal with edge update within the same leaf node
    if(!uWithinLeafNode.empty()){
        for(auto it=uWithinLeafNode.begin();it!=uWithinLeafNode.end();++it){//deal with each update within leaf node
            ID1 = it->first.first; ID2 = it->first.second; oldW = it->second.first; newW = it->second.second;
//            ID1_pos = Nodes[ID1].inleafpos; ID2_pos = Nodes[ID2].inleafpos;
            lnID = Nodes[ID1].inleaf;//get the leaf node id
            if(updatedLeafs.find(lnID)!=updatedLeafs.end())//if found
                continue;
            /// update leaf node
            updatedLeafs.insert(lnID);
            if(leafsToCheck.find(lnID) != leafsToCheck.end()){//if found in leafsToCheck
                leafsToCheck.erase(lnID);
            }
            //update of distance matrix of lnID: re-evaluate from each border
            cands = GTree[lnID].leafnodes;
            bool borderUpdate = false;
//            temp_pos = aBorderShortcuts.size();
            for(int i=0;i<GTree[lnID].borders.size();++i){// re-evaluation from each border to other vertices in the leaf node
                if(dijkstra_candidate_update(GTree[lnID].borders[i], i, cands,graph)){//if there is update among borders
                    borderUpdate = true;
                }
            }
            if(borderUpdate) {//if there is update among borders of the leaf node
//                aLeafNodes.insert({lnID, make_pair(temp_pos,aBorderShortcuts.size())});//insert to affected leaf nodes map
                LeafNodeContract(lnID, graph, leafsToCheck, updatedLeafs);//contract the leaf node
                borderAffectedNodes.insert(lnID);
                flag_borderUpdate = true;
            }
        }
    }
    /// deal with edge/shortcut update among borders
    // If there is update borders, we need to check all non-updated leaf nodes at first.
    // Since we do not know whether the index within leaf node will be affected, we have to check all leaf nodes.
    // So far, the distance matrix of some leaf nodes have been updated
    if(flag_borderUpdate){//!uAmongBorders.empty()
        Timer tt1;
        tt1.start();
        if(!ifPrue){//if not use prune strategy
            cout<<"Without pruning strategy!"<<endl;
            /// Method 1: The naive solution is to check all non-updated leaf node from very beginning
            for(int nid=0;nid<GTree.size();++nid){
                if(!GTree[nid].isleaf)//if it is not leaf node
                    continue;
                if( updatedLeafs.find(nid) == updatedLeafs.end()){ // if the leaf node has not been updated yet
                    borderUpdate = false;
                    /// Way 1: directly reconstruct
                    updatedLeafs.insert(nid);
                    //update of distance matrix of lnID: re-evaluate from each border
                    cands = GTree[nid].leafnodes;
                    bool borderUpdate = false;
    //            temp_pos = aBorderShortcuts.size();
                    for(int i=0;i<GTree[nid].borders.size();++i){// re-evaluation from each border to other vertices in the leaf node
                        if(dijkstra_candidate_update(GTree[nid].borders[i], i, cands,graph)){//if there is update among borders
                            borderUpdate = true;
                        }
                    }
                    if(borderUpdate) {//if there is update in the leaf node
    //                aLeafNodes.insert({lnID, make_pair(temp_pos,aBorderShortcuts.size())});//insert to affected leaf nodes map
                        LeafNodeContract(lnID, graph, updatedLeafs, leafsToCheck);//contract the leaf node
                        borderAffectedNodes.insert(lnID);
                    }
                    /// Way 2: check the border shortcuts and update
                    *//*nleafnode = GTree[nid].leafnodes.size();
                    int b1,b2;
                    int disbb,dis1,dis2;
                    for(int i=0;i<GTree[nid].borders.size()-1;++i){//for each border pair
                        b1 = GTree[nid].borders[i];
                        for(int j=i+1;j<GTree[nid].borders.size();++j){
                            b2 = GTree[nid].borders[j];
                            ID2_pos = Nodes[b2].inleafpos;
                            //check each updated border edge/shortcut
                            for(auto it=uAmongBorders.begin();it!=uAmongBorders.end();++it){
                                ID1 = it->first.first; ID2 = it->first.second; oldW = it->second.first; newW = it->second.second;
                                disbb = GTreeO[nid].mind[i*nleafnode + ID2_pos];
    //                            int dis = Dijkstra(b1,b2,NodesO);
    //                            if(disbb != dis)
    //                                cout<<b1<<" "<<b2<<": "<<disbb<<" "<<dis<<endl;
                                dis1 = ComputeDisByTree(b1,ID1,GTreeO);//search on old G-tree
    //                            dis = Dijkstra(b1,ID1,NodesO);
    //                            if(dis1 != dis)
    //                                cout<<b1<<" "<<ID1<<": "<<dis1<<" "<<dis<<endl;
                                dis2 = ComputeDisByTree(ID2,b2,GTreeO);
    //                            dis = Dijkstra(b2,ID2,NodesO);
    //                            if(dis2 != dis)
    //                                cout<<b2<<" "<<ID2<<": "<<dis2<<" "<<dis<<endl;
                                if(abs(disbb - dis1 - dis2) == oldW){//if the shortest path pass e(ID1,ID2), we need to re-compute the index of this leaf node
                                    borderUpdate = true;
                                    goto flag1;
                                }
                            }
                        }
                    }
                    flag1:
                    if(borderUpdate){
                        cout<<"!!!There is update."<<endl;
                        /// update leaf node
                        updatedLeafs.insert(lnID);
    //                    temp_pos = aBorderShortcuts.size();
                        for(int i=0;i<GTree[lnID].borders.size();++i){// re-evaluation from each border to other vertices in the leaf node
                            dijkstra_candidate_update(GTree[lnID].borders[i],i,cands,graph, uAmongBorders);
                        }
    //                    aLeafNodes.insert({lnID, make_pair(temp_pos,aBorderShortcuts.size())});//insert to affected leaf nodes map
                        LeafNodeContract(nid, graph, updatedLeafs, neighborNodesSet);
                    }*//*
                }
            }
            tt1.stop();
            cout<<"The number of updated leaf nodes: "<<updatedLeafs.size()<<endl;
            cout<<"The time used for checking leaf nodes: "<<tt1.GetRuntime()<<" s."<<endl;
            /// propagate to upper levels
            tt1.start();
            /// Way 1: Naive solution, simply update all upper-level nodes
            UpdatePropagateToUpperLevels(graph, updatedLeafs, true);
        }else{
            cout<<"With pruning strategy!"<<endl;
            /// Method 2: Propagate with pruning strategy
            // we can identify the affected area firstly, and propagate updates to other area (leaf nodes) if the border of this area is updated.
            vector<int> path_a, path_b;
            int lca;
            ///identify the affected area caused by cross-leaf edge updates
            for(int i=0;i<uCrossLeafNode.size();++i){
                ID1 = uCrossLeafNode[i].first.first; ID2 = uCrossLeafNode[i].first.second; oldW = uCrossLeafNode[i].second.first; newW = uCrossLeafNode[i].second.second;
                path_a = Nodes[ID1].gtreepath, path_b = Nodes[ID2].gtreepath;
                lca = path_a[find_LCA_pos(ID1,ID2)];
                //identify the affected area
                if(GTree[lca].isleaf){//if it is leaf node
                    if(updatedLeafs.find(lca) == updatedLeafs.end()){//if it has not been updated
                        leafsToCheck.insert(lca);
                    }
                }else{//if it is non-leaf node
                    GetChildNodes(lca, leafsToCheck, updatedLeafs);
                }
            }
            // then we update the affected area, i.e., we update the leaf nodes and propagate the updates if any
            while(!leafsToCheck.empty()){//if not empty
                lnID = *leafsToCheck.begin();
                leafsToCheck.erase(lnID);
                /// update leaf node
                updatedLeafs.insert(lnID);
                //update of distance matrix of lnID: re-evaluate from each border
                cands = GTree[lnID].leafnodes;
                bool borderUpdate = false;
                for(int i=0;i<GTree[lnID].borders.size();++i){// re-evaluation from each border to other vertices in the leaf node
                    if(dijkstra_candidate_update(GTree[lnID].borders[i], i, cands, graph)){//if there is update among borders
                        borderUpdate = true;
                    }
                }
                if(borderUpdate) {//if there is update in the leaf node
                    LeafNodeContract(lnID, graph, leafsToCheck, updatedLeafs);//contract the leaf node
                    borderAffectedNodes.insert(lnID);
                }
            }
            tt1.stop();
            cout<<"The number of updated leaf nodes: "<<updatedLeafs.size()<<endl;
            cout<<"The time used for checking leaf nodes: "<<tt1.GetRuntime()<<" s."<<endl;
            /// propagate to upper levels
            tt1.start();
            if(!borderAffectedNodes.empty()){
                cout<<"There is propagation to upper levels!"<<endl;
                /// Way 2: Propagate with pruning strategy
                UpdatePropagateToUpperLevelsPrune(graph, updatedLeafs,borderAffectedNodes);
            }
        }

        tt1.stop();
        cout<<"The time used for propagating updates to upper levels: "<<tt1.GetRuntime()<<" s."<<endl;
    }
    cout<<"Done."<<endl;
    tt.stop();
    cout<<"The time used for updating: "<<tt.GetRuntime()<<" s."<<endl;
}*/
void Gtree::LeafLevelUpdate(vector<Node> & graph, unordered_set<int> & leafsToCheck, unordered_set<int> & updatedNodes, unordered_set<int> & borderAffectedNodes, bool ifParallel){
    int ID1,ID2,oldW,newW;
    int lnID;
    vector<int> cands;

    int threadnum=1;
    unordered_map<int, bool> flag_update; flag_update.clear();
    // If there is update borders, we need to check all non-updated neighbor leaf nodes first.
    while(!leafsToCheck.empty()){//if not empty
        if(ifParallel){//multi-thread
            /// parallel computing by boost::thread_group
            if(leafsToCheck.size()>thread_num){
                threadnum = thread_num;
            }else{
                threadnum = leafsToCheck.size();
            }

            flag_update.clear();
            /// multi-thread
//              cout<<"Size of leafsToCheck: "<<leafsToCheck.size()<<endl;
            boost::thread_group threadf;
            auto it = leafsToCheck.begin();
            for(int i=0;i<threadnum;i++){
                lnID = *leafsToCheck.begin();
                leafsToCheck.erase(lnID);
                updatedNodes.insert(lnID);// update leaf node
                flag_update.insert({lnID,false});
                threadf.add_thread(new boost::thread(&Gtree::LeafUpdate, this, lnID, boost::ref(graph), boost::ref(leafsToCheck), boost::ref(updatedNodes), boost::ref(borderAffectedNodes), boost::ref(flag_update)));
            }
            threadf.join_all();
            // contraction cannot be paralleled
            for(auto it=flag_update.begin();it!=flag_update.end();++it){
                lnID = it->first;
                if(it->second){
                    LeafNodeContract(lnID, graph, leafsToCheck, updatedNodes);//contract the leaf node
                    borderAffectedNodes.insert(lnID);
                }
//                borderAffectedNodes.insert(lnID);//we have to check all higher-level tree nodes within the affected area
            }
        }
        else{//single thread
            /// single thread for updating leaf nodes
            lnID = *leafsToCheck.begin();
            leafsToCheck.erase(lnID);
            updatedNodes.insert(lnID);
            //update of distance matrix of lnID: re-evaluate from each border
            cands = GTree[lnID].leafnodes;
            bool borderUpdate = false;
            for(int i=0;i<GTree[lnID].borders.size();++i){// re-evaluation from each border to other vertices in the leaf node
                if(dijkstra_candidate_update(GTree[lnID].borders[i], i, cands,graph)){//if there is update among borders
                    borderUpdate = true;
                }
            }
            if(borderUpdate) {//if there is update among borders of the leaf node
                LeafNodeContract(lnID, graph, leafsToCheck, updatedNodes);//contract the leaf node
                borderAffectedNodes.insert(lnID);
            }
        }
    }

}

void Gtree::NonLeafLevelUpdate(vector<Node> & graph, int level_i,  vector< vector<int> > & treenodelevel, unordered_set<int> & nodesToCheck, unordered_set<int> & updatedNodes, unordered_set<int> & borderAffectedNodes, bool ifParallel){
    int ID1,ID2,oldW,newW;
    int lnID;
    vector<int> cands;

    int threadnum=1;
    unordered_map<int, bool> flag_update; flag_update.clear();
    // If there is update borders, we need to check all non-updated neighbor leaf nodes first.
    while(!nodesToCheck.empty()){//if not empty
        if(ifParallel){//multi-thread
            /// parallel computing by boost::thread_group
            if(nodesToCheck.size()>thread_num){
                threadnum = thread_num;
            }else{
                threadnum = nodesToCheck.size();
            }

            flag_update.clear();
            /// multi-thread
//              cout<<"Size of leafsToCheck: "<<leafsToCheck.size()<<endl;
            boost::thread_group threadf;
            auto it = nodesToCheck.begin();
            for(int i=0;i<threadnum;i++){
                lnID = *nodesToCheck.begin();
                nodesToCheck.erase(lnID);
                updatedNodes.insert(lnID);// update leaf node
                flag_update.insert({lnID,false});
                threadf.add_thread(new boost::thread(&Gtree::NonLeafUpdate, this, lnID, boost::ref(graph), boost::ref(flag_update)));
            }
            threadf.join_all();
            for(auto it=flag_update.begin();it!=flag_update.end();++it){
                lnID = it->first;
                unordered_set<int> borderSet;
                borderSet.clear();
                for(auto it=GTree[lnID].borders.begin();it!=GTree[lnID].borders.end();++it)
                    borderSet.insert(*it);
                if(it->second){
                    NonLeafNodeContract(lnID, level_i, graph, nodesToCheck, updatedNodes, treenodelevel, borderSet);//contract the leaf node
                    borderAffectedNodes.insert(lnID);
                }
//                borderAffectedNodes.insert(lnID);//we have to check all higher-level tree nodes within the affected area
            }
        }
        else{//single thread
            /// single thread for updating leaf nodes
            lnID = *nodesToCheck.begin();
            nodesToCheck.erase(lnID);
            updatedNodes.insert(lnID);
            //update of distance matrix of lnID: re-evaluate from each border

            unordered_set<int> borderSet;
            borderSet.clear();
            for(auto it=GTree[lnID].borders.begin();it!=GTree[lnID].borders.end();++it)
                borderSet.insert(*it);
//            cands = GTree[lnID].union_borders;
            bool borderUpdate = false;
            for(int i=0;i<GTree[lnID].union_borders.size();++i){// re-evaluation from each border to other vertices in the leaf node
                cands.clear();
                for(int j=i+1;j<GTree[lnID].union_borders.size();++j){
                    cands.emplace_back(GTree[lnID].union_borders[j]);
                }
                if(dijkstra_candidate_update(GTree[lnID].union_borders[i], i, lnID, cands,graph, borderSet)){//if there is update among borders
                    borderUpdate = true;
                }
            }
            if(borderUpdate) {//if there is update among borders of the leaf node
                NonLeafNodeContract(lnID, level_i, graph, nodesToCheck, updatedNodes, treenodelevel, borderSet);//contract the leaf node
                borderAffectedNodes.insert(lnID);
            }
        }
    }
}
//function of propagate update upwards
void Gtree::UpdateUpwardsPropagate(vector<Node> & graph, unordered_set<int> & updatedNodes, unordered_set<int> & borderAffectedNodes, bool ifParallel, int lca, unordered_set<int> & nodesToCheck){
    /// level traversal to get hierarchical tree
    vector< vector<int> > treenodelevel;//hierarchical tree, note that level 0 is root node
    vector<int> nodeToLevel;//map from node id to level id
    nodeToLevel.resize(GTree.size());
    vector<int> current;
    current.clear();
    current.push_back(0);//root node
    treenodelevel.push_back(current);
    nodeToLevel[0] = 0;
    vector<int> mid;
    while( current.size() != 0 ){
        mid = current;// intermediate tree node
        current.clear();
        for ( int i = 0; i < mid.size(); i++ ){
            for ( int j = 0; j < GTree[mid[i]].children.size(); j++ ){
                current.push_back( GTree[mid[i]].children[j] );
                nodeToLevel[GTree[mid[i]].children[j]] = treenodelevel.size();
            }
        }
        if ( current.size() == 0 ) break;
        treenodelevel.push_back( current );
    }

    int nid, level_i;
    vector<int> cands;
//    unordered_set<int> nodesToCheck;//the set of tree nodes for update checking
//    nodesToCheck.clear();
    unordered_set<int> bAffectedNodes;//store all tree nodes that have updated border shortcut
    /// propagate level by level until there is no affected tree node
    vector<unordered_set<int>> levelsToCheck;
    levelsToCheck.assign(treenodelevel.size(),unordered_set<int>());
    int level_h=treenodelevel.size();
    int level_l=0;
    // for the border affected nodes, we first need to make sure they are in the same level
    for(auto it: borderAffectedNodes){
        level_h = min(level_h, nodeToLevel[it]);
        level_l = max(level_l, nodeToLevel[it]-1);
        int father = GTree[it].father;
        levelsToCheck[nodeToLevel[father]].insert(father);
    }
    for(auto it: nodesToCheck){
//        level_h = min(level_h, nodeToLevel[it]);
        level_l = max(level_l, nodeToLevel[it]);
        levelsToCheck[nodeToLevel[it]].insert(it);
    }
    int pa = lca;
    while(pa>=0){
        levelsToCheck[nodeToLevel[pa]].insert(pa);
        if(pa != 0){
            pa = GTree[pa].father;
        }else{
            break;
        }
    }
    if(ifDebug){
        if(level_h != level_l){
            cout<<"The leaf nodes are in different tree levels! "<<level_h <<" "<<level_l<<endl;
        }
    }
    for (auto it = borderAffectedNodes.begin(); it != borderAffectedNodes.end(); ++it) {//get the nodesToCheck of father level
        int father = GTree[*it].father;
        level_i = nodeToLevel[*it];
        if(levelsToCheck[level_i-1].find(father) == levelsToCheck[level_i-1].end()){//if not found
            levelsToCheck[level_i-1].insert(father);
        }
    }

    for(level_i = level_l;level_i >= 0;--level_i){
        nodesToCheck = levelsToCheck[level_i];
        if(nodesToCheck.empty()){
            continue;
//            if(nodeToLevel[lca]>=level_i){
//                cout<<level_i<<" "<<nodeToLevel[lca]<<endl;
//            }else{
//                break;
//            }
        }// &&
        if(level_i > 0){
            bAffectedNodes.clear();
            /// update level_i
            NonLeafLevelUpdate(graph, level_i, treenodelevel, nodesToCheck, updatedNodes, bAffectedNodes, ifParallel);
            for (auto it = bAffectedNodes.begin(); it != bAffectedNodes.end(); ++it) {//get the nodesToCheck of father level
                int father = GTree[*it].father;
                if(levelsToCheck[level_i-1].find(father) == levelsToCheck[level_i-1].end()){//if not found
                    levelsToCheck[level_i-1].insert(father);
                }
            }
        }
        else{//for root node
            assert(nodesToCheck.size() == 1);
            /// single thread for updating leaf nodes
            int lnID = *nodesToCheck.begin();
            nodesToCheck.erase(lnID);
            updatedNodes.insert(lnID);
            //update of distance matrix of lnID: re-evaluate from each border

            unordered_set<int> borderSet;
            borderSet.clear();
            for(auto it=GTree[lnID].borders.begin();it!=GTree[lnID].borders.end();++it)
                borderSet.insert(*it);
//            cands = GTree[lnID].union_borders;
            bool borderUpdate = false;
            for(int i=0;i<GTree[lnID].union_borders.size();++i){// re-evaluation from each border to other vertices in the leaf node
                cands.clear();
                for(int j=i+1;j<GTree[lnID].union_borders.size();++j){
                    cands.emplace_back(GTree[lnID].union_borders[j]);
                }
                if(dijkstra_candidate_update(GTree[lnID].union_borders[i], i, lnID, cands,graph, borderSet)){//if there is update among borders
                    borderUpdate = true;
                }
            }
        }


    }

   /* /// process the first non-leaf node level
    bAffectedNodes.clear();
    for(auto it=levelsToCheck.begin();it!=levelsToCheck.end();++it){
        level_i = it->first; nodesToCheck = it->second;
        NonLeafLevelUpdate(graph, level_i, nodesToCheck, bAffectedNodes);
    }
    levelsToCheck.clear();

    /// process other non-leaf node level
    while(!bAffectedNodes.empty()) {// if bAffectedNodes is not empty, go to the next level
//        cout<<"Level "<<nodeToLevel[*bAffectedNodes.begin()]<<". Number of affected nodes: "<<bAffectedNodes.size()<<endl;
        for (auto it = bAffectedNodes.begin(); it != bAffectedNodes.end(); ++it) {//get the nodesToCheck of father level
            if(*it != 0)
                nodesToCheck.insert(GTree[*it].father);
        }
        bAffectedNodes.clear();

        NonLeafLevelUpdate(graph, level_i, nodesToCheck, bAffectedNodes);


        if(ifParallel){
            //deal with father level
            int threadnum=1;
            while (!nodesToCheck.empty()) {
                /// parallel computing by boost::thread_group
                if(nodesToCheck.size()>thread_num){
                    threadnum = thread_num;
                }else{
                    threadnum = nodesToCheck.size();
                }
//            threadnum = nodesToCheck.size();
                unordered_map<int, bool> flag_update; flag_update.clear();
                /// multi-thread
//            cout<<"Size of nodesToCheck: "<<nodesToCheck.size()<<endl;
                boost::thread_group threadf;
                for(int i=0;i<threadnum;i++){
                    //deal with the affected tree node: for each border-affected node, we need to check and update it
                    nid = *nodesToCheck.begin();
                    level_i = nodeToLevel[nid];
                    nodesToCheck.erase(nid);
                    updatedNodes.insert(nid);
//                if(nid == 0)
//                    cout<<nid<<endl;
                    flag_update.insert({nid,false});
                    threadf.add_thread(new boost::thread(&Gtree::NodeUpdate, this, nid, boost::ref(graph), boost::ref(flag_update)));
                }
                threadf.join_all();
                for(auto it=flag_update.begin();it!=flag_update.end();++it){
                    nid = it->first;
                    if(it->second){
                        unordered_set<int> borderSet;
                        borderSet.clear();
                        for(auto it=GTree[nid].borders.begin();it!=GTree[nid].borders.end();++it)
                            borderSet.insert(*it);
                        NonLeafNodeContract(nid, level_i, graph, nodesToCheck, updatedNodes, treenodelevel, borderSet);//contract the leaf node
//                    bAffectedNodes.insert(nid);
                    }
                    if(nid!=0)
                        bAffectedNodes.insert(nid);
                }
            }
        }else{

        }

    }*/
}
//function of dealing with batch edge increase update of G-tree
void Gtree::GtreeUpdateParalel(int updateType, int updateVolume, int updateBatch, bool ifParallel){//update(<ID1,ID2>,oldW) vector<pair<pair<int,int>,int> > & updates,
    // init
    init_update();
    // read updates
    string file = string(DataPath) + "/" + dataset + "/" + FILE_UPDATE;
    ReadUpdates(file);

    cout<<"Processing updates... "<<endl;
    //duplicate G-tree

    if(ifDebug){
        CorrectnessCheck(100);
    }
    double ave_time = 0;

    int ID1,ID2,oldW,newW;
    vector<pair<pair<int,int>,pair<int,int> > > uWithinLeafNode;// edge update within the same leaf node (<ID1,ID2>,<oldW,newW>)
    vector<pair<pair<int,int>,pair<int,int> > > uCrossLeafNode;// edge update among borders (<ID1,ID2>,<oldW,newW>)
    vector<pair<pair<int,int>,pair<int,int> > > uAmongBorders;// edge/shortcut update among borders (<ID1,ID2>,<oldW,newW>)
    bool flag_borderUpdate = false;//whether there is border update of leaf node
    Timer tt;
    ///preprocess: get the update edges
    assert(updateVolume<=updateEdges.size());
    for(int batch_i=0;batch_i<updateBatch;++batch_i){
        cout<<"Batch "<<batch_i<<":";
        GTreeO = GTree; //old GTree
        NodesO = Nodes;//old Nodes
        uWithinLeafNode.clear();
        uCrossLeafNode.clear();
        uAmongBorders.clear();
        flag_borderUpdate = false;
        tt.start();
//        cout<<"Before update: "<<Dijkstra(759,760,Nodes)<<endl;
        for(int i=batch_i;i<batch_i+updateVolume;++i){//for each edge update
            ID1 = updateEdges[i].first.first;
            ID2 = updateEdges[i].first.second;
//            if(ID1 == 759 && ID2 == 760){
//                cout<<"!!"<<endl;
//                cout << ID1<<" "<<Nodes[ID1].adjnodes.size()<< " "<<Nodes[ID1].adjweight.size()<<endl;
//                cout << ID2<<" "<<Nodes[ID2].adjnodes.size()<< " "<<Nodes[ID2].adjweight.size()<<endl;
//            }
            cout<<" " <<ID1 << " "<<ID2;
            oldW = updateEdges[i].second;
            switch (updateType) {
                case INCREASE:
                    newW = (int)(2*oldW);//update
                    break;
                case DECREASE:
                    newW = (int)(0.5*oldW);
                    break;
                case MIX:
                    newW = (int)(2*oldW * (rand()/double(RAND_MAX)));
                    break;
                default:
                    cout<<"Wrong update type!"<<endl; break;
            }
            cout<<" ("<<oldW<<"->"<<newW<<") ";
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
        cout<<endl;
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
        tt.stop();
        cout<<"The time used for updating: "<<tt.GetRuntime()<<" s."<<endl;
        ave_time += tt.GetRuntime();
        if(ifDebug){
            CorrectnessCheck(100);
        }

    }
    cout<<"The average time used for updating: "<< ave_time/updateBatch <<" s."<<endl;
}
 /// old version-2022.10.11
/* void Gtree::GtreeUpdateParalel(int updateType, int updateVolume, int updateBatch, bool ifParallel){//update(<ID1,ID2>,oldW) vector<pair<pair<int,int>,int> > & updates,
    // init
    init_update();
    // read updates
    string file = string(DataPath) + dataset + "/" + FILE_UPDATE;
    ReadUpdates(file);

    cout<<"Processing updates... "<<endl;
    //duplicate G-tree

    CorrectnessCheck(100);
    int ID1,ID2,oldW,newW;
    vector<pair<pair<int,int>,pair<int,int> > > uWithinLeafNode;// edge update within the same leaf node (<ID1,ID2>,<oldW,newW>)
    vector<pair<pair<int,int>,pair<int,int> > > uCrossLeafNode;// edge update among borders (<ID1,ID2>,<oldW,newW>)
    vector<pair<pair<int,int>,pair<int,int> > > uAmongBorders;// edge/shortcut update among borders (<ID1,ID2>,<oldW,newW>)
    bool flag_borderUpdate = false;//whether there is border update of leaf node
    Timer tt;
    ///preprocess: update edges
    assert(updateVolume<=updateEdges.size());
    for(int batch_i=0;batch_i<updateBatch;++batch_i){
        GTreeO = GTree; //old GTree
        NodesO = Nodes;//old Nodes
        tt.start();
        for(int i=batch_i;i<batch_i+updateVolume;++i){//for each edge update
            ID1 = updateEdges[i].first.first;
            ID2 = updateEdges[i].first.second;
            oldW = updateEdges[i].second;
            switch (updateType) {
                case INCREASE:
                    newW = (int)(2*oldW);//update
                    break;
                case DECREASE:
                    newW = (int)(0.5*oldW);
                    break;
                case MIX:
                    newW = (int)(2*oldW * (rand()/double(RAND_MAX)));
                    break;
                default:
                    cout<<"Wrong update type!"<<endl; break;
            }
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
        /// update processing
        int lnID;
        int ID1_pos,ID2_pos;
        vector<Node> graph;//temp graph
        graph = Nodes;//original graph

        bool borderUpdate = false;
        int nleafnode;
        unordered_set<int> leafsToCheck;//the set of leaf nodes for update checking
        leafsToCheck.clear();
        unordered_set<int> updatedLeafs;
        updatedLeafs.clear();
        unordered_set<int> borderAffectedNodes;
        borderAffectedNodes.clear();
//    list<int> neighborNodes;
        vector<int> cands;
//    vector<pair<pair<int,int>,pair<int,int>>> aBorderShortcuts;//affected shortcuts among borders of leaf node
//    set<int> aLeafNodes;// set of affected leaf nodes

        int temp_pos;
        /// deal with edge update within the same leaf node
        if(!uWithinLeafNode.empty()){
            for(auto it=uWithinLeafNode.begin();it!=uWithinLeafNode.end();++it){//deal with each update within leaf node
                ID1 = it->first.first; ID2 = it->first.second; oldW = it->second.first; newW = it->second.second;
//            ID1_pos = Nodes[ID1].inleafpos; ID2_pos = Nodes[ID2].inleafpos;
                lnID = Nodes[ID1].inleaf;//get the leaf node id
                if(updatedLeafs.find(lnID)!=updatedLeafs.end())//if found
                    continue;
                /// update leaf node
                updatedLeafs.insert(lnID);
                if(leafsToCheck.find(lnID) != leafsToCheck.end()){//if found in leafsToCheck
                    leafsToCheck.erase(lnID);
                }
                //update of distance matrix of lnID: re-evaluate from each border
                cands = GTree[lnID].leafnodes;
                bool borderUpdate = false;
//            temp_pos = aBorderShortcuts.size();
                for(int i=0;i<GTree[lnID].borders.size();++i){// re-evaluation from each border to other vertices in the leaf node
                    if(dijkstra_candidate_update(GTree[lnID].borders[i], i, cands,graph)){//if there is update among borders
                        borderUpdate = true;
                    }
                }
                if(borderUpdate) {//if there is update among borders of the leaf node
//                aLeafNodes.insert({lnID, make_pair(temp_pos,aBorderShortcuts.size())});//insert to affected leaf nodes map
                    LeafNodeContract(lnID, graph, leafsToCheck, updatedLeafs);//contract the leaf node
                    borderAffectedNodes.insert(lnID);
                    flag_borderUpdate = true;
                }
            }
        }
        /// deal with edge/shortcut update among borders
        // If there is update borders, we need to check all non-updated leaf nodes at first.
        // Since we do not know whether the index within leaf node will be affected, we have to check all leaf nodes.
        // So far, the distance matrix of some leaf nodes have been updated
        if(flag_borderUpdate){//!uAmongBorders.empty()
            Timer tt1;
            tt1.start();
//            cout<<"With pruning strategy!"<<endl;
//            if(ifParallel)
//                cout<<"With multi-thread computation!"<<endl;
//            else{
//                cout<<"Without multi-thread computation!"<<endl;
//            }
            /// Method 2: Propagate with pruning strategy
            // we can identify the affected area firstly, and propagate updates to other area (leaf nodes) if the border of this area is updated.
            vector<int> path_a, path_b;
            int lca;
            ///identify the affected area caused by cross-leaf edge updates
            for(int i=0;i<uCrossLeafNode.size();++i){
                ID1 = uCrossLeafNode[i].first.first; ID2 = uCrossLeafNode[i].first.second; oldW = uCrossLeafNode[i].second.first; newW = uCrossLeafNode[i].second.second;
                path_a = Nodes[ID1].gtreepath, path_b = Nodes[ID2].gtreepath;
                lca = path_a[find_LCA_pos(ID1,ID2)];
                //identify the affected area
                if(GTree[lca].isleaf){//if it is leaf node
                    if(updatedLeafs.find(lca) == updatedLeafs.end()){//if it has not been updated
                        leafsToCheck.insert(lca);
                    }
                }else{//if it is non-leaf node
                    GetChildLeafs(lca, leafsToCheck, updatedLeafs);
                }
            }
            // then we update the affected area, i.e., we update the leaf nodes and propagate the updates if any
            int threadnum=1;
            unordered_map<int, bool> flag_update; flag_update.clear();
            while(!leafsToCheck.empty()){//if not empty
                if(ifParallel){//multi-thread
                    /// parallel computing by boost::thread_group
                    if(leafsToCheck.size()>thread_num){
                        threadnum = thread_num;
                    }else{
                        threadnum = leafsToCheck.size();
                    }

                    flag_update.clear();
                    /// multi-thread
//              cout<<"Size of leafsToCheck: "<<leafsToCheck.size()<<endl;
                    boost::thread_group threadf;
                    auto it = leafsToCheck.begin();
                    for(int i=0;i<threadnum;i++){
                        lnID = *leafsToCheck.begin();
                        leafsToCheck.erase(lnID);
                        updatedLeafs.insert(lnID);// update leaf node
                        flag_update.insert({lnID,false});
                        threadf.add_thread(new boost::thread(&Gtree::LeafUpdate, this, lnID, boost::ref(graph), boost::ref(leafsToCheck), boost::ref(updatedLeafs), boost::ref(borderAffectedNodes), boost::ref(flag_update)));
                    }
                    threadf.join_all();
                    for(auto it=flag_update.begin();it!=flag_update.end();++it){
                        lnID = it->first;
                        if(it->second){
                            LeafNodeContract(lnID, graph, leafsToCheck, updatedLeafs);//contract the leaf node
//                            borderAffectedNodes.insert(lnID);
                        }
                        borderAffectedNodes.insert(lnID);//we have to check all higher-level tree nodes within the affected area
                    }
                }else{//single thread
                    /// single thread
                    lnID = *leafsToCheck.begin();
                    leafsToCheck.erase(lnID);
                    updatedLeafs.insert(lnID);// update leaf node
                    //update of distance matrix of lnID: re-evaluate from each border
                    cands = GTree[lnID].leafnodes;
                    bool borderUpdate = false;
//            temp_pos = aBorderShortcuts.size();
                    for(int i=0;i<GTree[lnID].borders.size();++i){// re-evaluation from each border to other vertices in the leaf node
                        if(dijkstra_candidate_update(GTree[lnID].borders[i], i, cands,graph)){//if there is update among borders
                            borderUpdate = true;
                        }
                    }
                    if(borderUpdate) {//if there is update among borders of the leaf node
//                aLeafNodes.insert({lnID, make_pair(temp_pos,aBorderShortcuts.size())});//insert to affected leaf nodes map
                        LeafNodeContract(lnID, graph, leafsToCheck, updatedLeafs);//contract the leaf node
                        borderAffectedNodes.insert(lnID);
                    }
                }
            }
            tt1.stop();
//            cout<<"The number of updated leaf nodes: "<<updatedLeafs.size()<<endl;
//            cout<<"The time used for checking leaf nodes: "<<tt1.GetRuntime()<<" s."<<endl;

            /// propagate to upper levels
            tt1.start();
            if(!borderAffectedNodes.empty()){
//                cout<<"There is propagation to upper levels!"<<endl;
                /// Way 2: Propagate with pruning strategy
                if(ifParallel){
                    UpdatePropagateToUpperLevelsPruneP(graph, updatedLeafs,borderAffectedNodes);//multi-thread
                }else{
                    UpdatePropagateToUpperLevelsPrune(graph, updatedLeafs,borderAffectedNodes);//single thread
                }
            }
            tt1.stop();
//            cout<<"The time used for propagating updates to upper levels: "<<tt1.GetRuntime()<<" s."<<endl;
        }
//        cout<<"Done."<<endl;
        tt.stop();
        cout<<"Batch "<<batch_i <<". The time used for updating: "<<tt.GetRuntime()<<" s."<<endl;
        CorrectnessCheck(100);
    }
}*/
/*void Gtree::GtreeUpdateParalel(int updateType, int updateVolume, int updateBatch, bool ifParallel){//update(<ID1,ID2>,oldW) vector<pair<pair<int,int>,int> > & updates,
    // init
    init_update();
    // read updates
    string file = string(DataPath) + dataset + "/" + FILE_UPDATE;
    ReadUpdates(file);

    cout<<"Processing updates... Volume of update: "<<updateVolume<<endl;
    //duplicate G-tree
    GTreeO = GTree; //old GTree
    NodesO = Nodes;//old Nodes

    int ID1,ID2,oldW,newW;
    vector<pair<pair<int,int>,pair<int,int> > > uWithinLeafNode;// edge update within the same leaf node (<ID1,ID2>,<oldW,newW>)
    vector<pair<pair<int,int>,pair<int,int> > > uCrossLeafNode;// edge update among borders (<ID1,ID2>,<oldW,newW>)
    vector<pair<pair<int,int>,pair<int,int> > > uAmongBorders;// edge/shortcut update among borders (<ID1,ID2>,<oldW,newW>)
    bool flag_borderUpdate = false;//whether there is border update of leaf node
    Timer tt;
    tt.start();
    ///preprocess: update edges
    assert(updateVolume<=updateEdges.size());
    for(int i=0;i<updateVolume;++i){//for each edge update
        ID1 = updateEdges[i].first.first;
        ID2 = updateEdges[i].first.second;
        oldW = updateEdges[i].second;
        switch (updateType) {
            case INCREASE:
                newW = 1.5*oldW;//update
                break;
            case DECREASE:
                newW = 0.5*oldW;
                break;
            case MIX:
                newW = 2*oldW * (rand()/double(RAND_MAX));
                break;
            default:
                cout<<"Wrong update type!"<<endl; break;
        }
        //update edge weight
        for(int i=0;i<Nodes[ID1].adjnodes.size();i++){
            if(Nodes[ID1].adjnodes[i]==ID2){
                assert(Nodes[ID1].adjweight[i] == oldW);
                Nodes[ID1].adjweight[i]=newW;
                break;
            }
        }
        for(int i=0;i<Nodes[ID2].adjnodes.size();i++){
            if(Nodes[ID2].adjnodes[i]==ID1){
                assert(Nodes[ID2].adjweight[i] == oldW);
                Nodes[ID2].adjweight[i]=newW;
                break;
            }
        }
        //identify the edge type
        if(Nodes[ID1].inleaf == Nodes[ID2].inleaf){
            uWithinLeafNode.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
        }else{
            uCrossLeafNode.emplace_back(make_pair(ID1,ID2),make_pair(oldW,newW));
            uAmongBorders.emplace_back(make_pair(ID1,ID2),make_pair(oldW,newW));
            flag_borderUpdate = true;
        }
    }
    /// update processing
    int lnID;
    int ID1_pos,ID2_pos;
    vector<Node> graph;//temp graph
    graph = Nodes;//original graph

    bool borderUpdate = false;
    int nleafnode;
    unordered_set<int> leafsToCheck;//the set of leaf nodes for update checking
    leafsToCheck.clear();
    unordered_set<int> updatedLeafs;
    updatedLeafs.clear();
    unordered_set<int> borderAffectedNodes;
    borderAffectedNodes.clear();
//    list<int> neighborNodes;
    vector<int> cands;
//    vector<pair<pair<int,int>,pair<int,int>>> aBorderShortcuts;//affected shortcuts among borders of leaf node
//    set<int> aLeafNodes;// set of affected leaf nodes

    int temp_pos;
    /// deal with edge update within the same leaf node
    if(!uWithinLeafNode.empty()){
        for(auto it=uWithinLeafNode.begin();it!=uWithinLeafNode.end();++it){//deal with each update within leaf node
            ID1 = it->first.first; ID2 = it->first.second; oldW = it->second.first; newW = it->second.second;
//            ID1_pos = Nodes[ID1].inleafpos; ID2_pos = Nodes[ID2].inleafpos;
            lnID = Nodes[ID1].inleaf;//get the leaf node id
            if(updatedLeafs.find(lnID)!=updatedLeafs.end())//if found
                continue;
            /// update leaf node
            updatedLeafs.insert(lnID);
            if(leafsToCheck.find(lnID) != leafsToCheck.end()){//if found in leafsToCheck
                leafsToCheck.erase(lnID);
            }
            //update of distance matrix of lnID: re-evaluate from each border
            cands = GTree[lnID].leafnodes;
            bool borderUpdate = false;
//            temp_pos = aBorderShortcuts.size();
            for(int i=0;i<GTree[lnID].borders.size();++i){// re-evaluation from each border to other vertices in the leaf node
                if(dijkstra_candidate_update(GTree[lnID].borders[i], i, cands,graph)){//if there is update among borders
                    borderUpdate = true;
                }
            }
            if(borderUpdate) {//if there is update among borders of the leaf node
//                aLeafNodes.insert({lnID, make_pair(temp_pos,aBorderShortcuts.size())});//insert to affected leaf nodes map
                LeafNodeContract(lnID, graph, leafsToCheck, updatedLeafs);//contract the leaf node
                borderAffectedNodes.insert(lnID);
                flag_borderUpdate = true;
            }
        }
    }
    /// deal with edge/shortcut update among borders
    // If there is update borders, we need to check all non-updated leaf nodes at first.
    // Since we do not know whether the index within leaf node will be affected, we have to check all leaf nodes.
    // So far, the distance matrix of some leaf nodes have been updated
    if(flag_borderUpdate){//!uAmongBorders.empty()
        Timer tt1;
        tt1.start();
        cout<<"With pruning strategy!"<<endl;
        if(ifParallel)
            cout<<"With multi-thread computation!"<<endl;
        else{
            cout<<"Without multi-thread computation!"<<endl;
        }
        /// Method 2: Propagate with pruning strategy
        // we can identify the affected area firstly, and propagate updates to other area (leaf nodes) if the border of this area is updated.
        vector<int> path_a, path_b;
        int lca;
        ///identify the affected area caused by cross-leaf edge updates
        for(int i=0;i<uCrossLeafNode.size();++i){
            ID1 = uCrossLeafNode[i].first.first; ID2 = uCrossLeafNode[i].first.second; oldW = uCrossLeafNode[i].second.first; newW = uCrossLeafNode[i].second.second;
            path_a = Nodes[ID1].gtreepath, path_b = Nodes[ID2].gtreepath;
            lca = path_a[find_LCA_pos(ID1,ID2)];
            //identify the affected area
            if(GTree[lca].isleaf){//if it is leaf node
                if(updatedLeafs.find(lca) == updatedLeafs.end()){//if it has not been updated
                    leafsToCheck.insert(lca);
                }
            }else{//if it is non-leaf node
                GetChildLeafs(lca, leafsToCheck, updatedLeafs);
            }
        }
        // then we update the affected area, i.e., we update the leaf nodes and propagate the updates if any
        int threadnum=1;
        unordered_map<int, bool> flag_update;
        while(!leafsToCheck.empty()){//if not empty
            if(ifParallel){//multi-thread
                /// parallel computing by boost::thread_group
                if(leafsToCheck.size()>thread_num){
                    threadnum = thread_num;
                }else{
                    threadnum = leafsToCheck.size();
                }
//                threadnum = leafsToCheck.size();
                flag_update.clear();
                /// multi-thread
//              cout<<"Size of leafsToCheck: "<<leafsToCheck.size()<<endl;
                boost::thread_group threadf;
                auto it = leafsToCheck.begin();
                for(int i=0;i<threadnum;i++){
                    lnID = *leafsToCheck.begin();
                    leafsToCheck.erase(lnID);
                    updatedLeafs.insert(lnID);// update leaf node
                    flag_update.insert({lnID,false});
                    threadf.add_thread(new boost::thread(&Gtree::LeafUpdate, this, lnID, boost::ref(graph), boost::ref(leafsToCheck), boost::ref(updatedLeafs), boost::ref(borderAffectedNodes), boost::ref(flag_update)));
                }
                threadf.join_all();
                for(auto it=flag_update.begin();it!=flag_update.end();++it){
                    lnID = it->first;
                    if(it->second){
                        LeafNodeContract(lnID, graph, leafsToCheck, updatedLeafs);//contract the leaf node
                        borderAffectedNodes.insert(lnID);
                    }
                }
            }else{//single thread
                /// single thread
                lnID = *leafsToCheck.begin();
                leafsToCheck.erase(lnID);
                updatedLeafs.insert(lnID);// update leaf node
                LeafUpdate(lnID, graph, leafsToCheck, updatedLeafs, borderAffectedNodes,flag_update );
            }
        }
        tt1.stop();
        cout<<"The number of updated leaf nodes: "<<updatedLeafs.size()<<endl;
        cout<<"The time used for checking leaf nodes: "<<tt1.GetRuntime()<<" s."<<endl;
        /// propagate to upper levels
        tt1.start();
        if(!borderAffectedNodes.empty()){
            cout<<"There is propagation to upper levels!"<<endl;
            /// Way 2: Propagate with pruning strategy
            if(ifParallel){
                UpdatePropagateToUpperLevelsPruneP(graph, updatedLeafs,borderAffectedNodes);//multi-thread
            }else{
                UpdatePropagateToUpperLevelsPrune(graph, updatedLeafs,borderAffectedNodes);//single thread
            }
        }
        tt1.stop();
        cout<<"The time used for propagating updates to upper levels: "<<tt1.GetRuntime()<<" s."<<endl;
    }
    cout<<"Done."<<endl;
    tt.stop();
    cout<<"The time used for updating: "<<tt.GetRuntime()<<" s."<<endl;
}*/
//function for leaf node updating
void Gtree::LeafUpdate(int lnID, vector<Node> & graph, unordered_set<int> & leafsToCheck, unordered_set<int> & updatedLeafs, unordered_set<int> & borderAffectedNodes, unordered_map<int, bool> & flag_update){
//    lnID = *leafsToCheck.begin();
//    leafsToCheck.erase(lnID);
//    /// update leaf node
//    updatedLeafs.insert(lnID);
    //update of distance matrix of lnID: re-evaluate from each border
    vector<int> cands = GTree[lnID].leafnodes;
    bool borderUpdate = false;
    for(int i=0;i<GTree[lnID].borders.size();++i){// re-evaluation from each border to other vertices in the leaf node
        if(dijkstra_candidate_update(GTree[lnID].borders[i], i, cands, graph)){//if there is update among borders
            borderUpdate = true;
        }
    }
    if(borderUpdate) {//if there is update in the leaf node
//        if(flag_update.empty()){
//            LeafNodeContract(lnID, graph, leafsToCheck, updatedLeafs);//contract the leaf node
//            borderAffectedNodes.insert(lnID);
//        }else{
//            flag_update[lnID] = true;
//        }
//        LeafNodeContract(lnID, graph, leafsToCheck, updatedLeafs);//contract the leaf node
//        borderAffectedNodes.insert(lnID);
        if(flag_update.find(lnID)!=flag_update.end()){//if found
            if(!flag_update[lnID])
                flag_update[lnID] = true;
        }else{
            cout<<"Flag update error!"<<endl;
        }
    }
}
//function of non-leaf node updating
void Gtree::NonLeafUpdate(int nid, vector<Node> & graph, unordered_map<int, bool> & flag_update){
    unordered_set<int> borderSet; borderSet.clear();
    for(auto it=GTree[nid].borders.begin();it!=GTree[nid].borders.end();++it)
        borderSet.insert(*it);
    //update of distance matrix of lnID: re-evaluate from each border
    vector<int> cands = GTree[nid].union_borders;
    bool borderUpdate = false;
    for (int i = 0; i < GTree[nid].union_borders.size()-1; ++i) {// re-evaluation from each border to other vertices in the leaf node
        cands.clear();
        for(int j=i+1;j<GTree[nid].union_borders.size();++j){
            cands.emplace_back(GTree[nid].union_borders[j]);
        }
        if (dijkstra_candidate_update(GTree[nid].union_borders[i], i, nid, cands, graph,borderSet)) {//if there is update among borders
            borderUpdate = true;
        }
    }
    if(borderUpdate) {//if there is update in the leaf node
        flag_update[nid] = true;
    }
}
//function of propagating updates to upper levels
void Gtree::UpdatePropagateToUpperLevels(vector<Node> & graph, unordered_set<int> & updatedLeafs, bool ifLeaf){
    // level traversal
    vector< vector<int> > treenodelevel;//hierarchical tree

    vector<int> current;
    current.clear();
    current.push_back(0);//root node
    treenodelevel.push_back(current);

    vector<int> mid;
    while( current.size() != 0 ){
        mid = current;// intermediate tree node
        current.clear();
        for ( int i = 0; i < mid.size(); i++ ){
            for ( int j = 0; j < GTree[mid[i]].children.size(); j++ ){
                current.push_back( GTree[mid[i]].children[j] );
            }
        }
        if ( current.size() == 0 ) break;
        treenodelevel.push_back( current );
    }
    /// It seems that G-tree is an unbalanced tree
    // bottom up calculation
    // temp graph
    vector<int> cands;
    vector<int> result;
    unordered_map<int, unordered_map<int,int> > vertex_pairs;//result of distance matrix

    // do dijkstra
    int s, t, tn, nid, cid, weight;
    vector<int> tnodes, tweight;
    set<int> nset;

    for ( int i = treenodelevel.size() - 1; i >= 0; i-- ){//start from the lowest level
        for ( int j = 0; j < treenodelevel[i].size(); j++ ){//for each partition in this level
            tn = treenodelevel[i][j];
            if(!ifLeaf){//if no need to check leaf node
                if(GTree[tn].isleaf)//if this is leaf node, skip it
                    continue;
            }else{//if we need to check leaf node
                if(updatedLeafs.find(tn) != updatedLeafs.end()){//if found, ie, the leaf node has been updated
                    continue;
                }
            }

            cands.clear();
            if ( GTree[tn].isleaf ){//for leaf node
                // cands = leafnodes
                cands = GTree[tn].leafnodes;
                // union borders = borders;
                GTree[tn].union_borders = GTree[tn].borders;
            }
            else{//for non-leaf node
//                nset.clear();
//                for ( int k = 0; k < GTree[tn].children.size(); k++ ){
//                    cid = GTree[tn].children[k];
//                    nset.insert( GTree[cid].borders.begin(), GTree[cid].borders.end() );//get the borders of all children
//                }
//                // union borders = cands;
//
//                cands.clear();
//                for ( set<int>::iterator it = nset.begin(); it != nset.end(); it ++ ){//push each vertex of this partition into cands
//                    cands.push_back( *it );
//                }
//                GTree[tn].union_borders = cands;//for non-leaf node, the union_borders contains all the borders of children
                cands = GTree[tn].union_borders;
            }

            // start to do min dis
            vertex_pairs.clear();

            // for each border, do min dis
            int cc = 0;
            int span = GTree[tn].union_borders.size();
            for ( int k = 0; k < GTree[tn].union_borders.size(); k++ ){
                result = dijkstra_candidate( GTree[tn].union_borders[k], cands, graph );//distance vector from s to all borders
                // save to map
                for ( int p = 0; p < result.size(); p ++ ){
                    if(GTree[tn].mind[k*span + p] != result[p]){
//                        int temp_dis = Dijkstra(GTree[tn].union_borders[k],cands[p],NodesO);
//                        if(temp_dis != GTree[tn].mind[k*span + p])
//                            cout<<temp_dis<<" "<<GTree[tn].mind[k*span + p]<<endl;
//                        temp_dis = Dijkstra(GTree[tn].union_borders[k],cands[p],graph);
//                        if(temp_dis != result[p])
//                            cout<<temp_dis<<" "<<result[p]<<endl;
//                        cout<<GTree[tn].mind[k*span + p]<<" "<<result[p]<<endl;
                        GTree[tn].mind[k*span + p] = result[p];
                    }
                    vertex_pairs[GTree[tn].union_borders[k]][cands[p]] = result[p];
                }
            }

            // IMPORTANT! after all border finished, degenerate graph
            // first, remove inward edges
            for ( int k = 0; k < GTree[tn].borders.size(); k++ ){//for each border vertex
                s = GTree[tn].borders[k];
                tnodes.clear();
                tweight.clear();
                for ( int p = 0; p < graph[s].adjnodes.size(); p++ ){//for each adjacent vertex
                    nid = graph[s].adjnodes[p];
                    weight = graph[s].adjweight[p];
                    // if adj node in same tree node

                    if ( graph[nid].gtreepath.size() <= i || graph[nid].gtreepath[i] != tn ){///
                        // only remain those useful vertices, i.e., borders
                        tnodes.push_back(nid);
                        tweight.push_back(weight);
                    }
                }
                // cut it
                graph[s].adjnodes = tnodes;//update the adjacency lists of graph, only left the useful boundary vertices
                graph[s].adjweight = tweight;
            }
            // second, add inter connected edges (shortcuts)
            for ( int k = 0; k < GTree[tn].borders.size(); k++ ){
                for ( int p = 0; p < GTree[tn].borders.size(); p++ ){
                    if ( k == p ) continue;
                    s = GTree[tn].borders[k];
                    t = GTree[tn].borders[p];
                    graph[s].adjnodes.push_back( t );
                    graph[s].adjweight.push_back( vertex_pairs[s][t] );
                }
            }
        }
    }
}
//function of propagating updates to upper levels with pruning strategy
void Gtree::UpdatePropagateToUpperLevelsPrune(vector<Node> & graph, unordered_set<int> & updatedLeafs, unordered_set<int> & borderAffectedNodes){
    /// level traversal to get hierarchical tree
    vector< vector<int> > treenodelevel;//hierarchical tree, note that level 0 is root node
    vector<int> nodeToLevel;//map from node id to level id
    nodeToLevel.resize(GTree.size());
    vector<int> current;
    current.clear();
    current.push_back(0);//root node
    treenodelevel.push_back(current);
    nodeToLevel[0] = 0;
    vector<int> mid;


    while( current.size() != 0 ){
        mid = current;// intermediate tree node
        current.clear();
        for ( int i = 0; i < mid.size(); i++ ){
            for ( int j = 0; j < GTree[mid[i]].children.size(); j++ ){
                current.push_back( GTree[mid[i]].children[j] );
                nodeToLevel[GTree[mid[i]].children[j]] = treenodelevel.size();
            }
        }
        if ( current.size() == 0 ) break;
        treenodelevel.push_back( current );
    }

    int nid, level_i;
    vector<int> cands;
    unordered_set<int> nodesToCheck;//the set of tree nodes for update checking
    nodesToCheck.clear();
    unordered_set<int> bAffectedNodes = borderAffectedNodes;//store all tree nodes that have updated border shortcut
    unordered_set<int> updatedNodes;//set of updated non-leaf nodes
    updatedNodes.clear();
    /// propagate level by level until there is no affected tree node
    while(!bAffectedNodes.empty()) {// if bAffectedNodes is not empty, go to the next level
//        cout<<"Level "<<nodeToLevel[*bAffectedNodes.begin()]<<". Number of affected nodes: "<<bAffectedNodes.size()<<endl;
        for (auto it = bAffectedNodes.begin(); it != bAffectedNodes.end(); ++it) {//get the nodesToCheck of father level
            if(*it!=0)
                nodesToCheck.insert(GTree[*it].father);
        }
        bAffectedNodes.clear();
        //deal with father level
        while (!nodesToCheck.empty()) {
            //deal with the affected tree node: for each border-affected node, we need to check and update it
            nid = *nodesToCheck.begin();
            level_i = nodeToLevel[nid];
            nodesToCheck.erase(nid);
            updatedNodes.insert(nid);
//            if(nid == 0)
//                cout<<nid<<endl;
            unordered_set<int> borderSet; borderSet.clear();
            for(auto it=GTree[nid].borders.begin();it!=GTree[nid].borders.end();++it)
                borderSet.insert(*it);
            //update of distance matrix of lnID: re-evaluate from each border
            cands = GTree[nid].union_borders;
            bool borderUpdate = false;
            for (int i = 0; i < GTree[nid].union_borders.size()-1; ++i) {// re-evaluation from each border to other vertices in the leaf node
                cands.clear();
                for(int j=i+1;j<GTree[nid].union_borders.size();++j){
                    cands.emplace_back(GTree[nid].union_borders[j]);
                }
                if (dijkstra_candidate_update(GTree[nid].union_borders[i], i, nid, cands, graph,borderSet)) {//if there is update among borders
                    borderUpdate = true;
                }
            }
            if (borderUpdate) {//if there is update in the tree node
//                NonLeafNodeContract(nid, level_i, graph, nodesToCheck, updatedNodes, treenodelevel, borderSet);//contract the leaf node
//                bAffectedNodes.insert(nid);
            }
            if(nid!=0)
                bAffectedNodes.insert(nid);
        }
    }
}
//function of propagating updates to upper levels with pruning strategy parallelly
void Gtree::UpdatePropagateToUpperLevelsPruneP(vector<Node> & graph, unordered_set<int> & updatedLeafs, unordered_set<int> & borderAffectedNodes){
    /// level traversal to get hierarchical tree
    vector< vector<int> > treenodelevel;//hierarchical tree, note that level 0 is root node
    vector<int> nodeToLevel;//map from node id to level id
    nodeToLevel.resize(GTree.size());
    vector<int> current;
    current.clear();
    current.push_back(0);//root node
    treenodelevel.push_back(current);
    nodeToLevel[0] = 0;
    vector<int> mid;
    while( current.size() != 0 ){
        mid = current;// intermediate tree node
        current.clear();
        for ( int i = 0; i < mid.size(); i++ ){
            for ( int j = 0; j < GTree[mid[i]].children.size(); j++ ){
                current.push_back( GTree[mid[i]].children[j] );
                nodeToLevel[GTree[mid[i]].children[j]] = treenodelevel.size();
            }
        }
        if ( current.size() == 0 ) break;
        treenodelevel.push_back( current );
    }

    int nid, level_i;
    vector<int> cands;
    unordered_set<int> nodesToCheck;//the set of tree nodes for update checking
    nodesToCheck.clear();
    unordered_set<int> bAffectedNodes = borderAffectedNodes;//store all tree nodes that have updated border shortcut
    unordered_set<int> updatedNodes;//set of updated non-leaf nodes
    updatedNodes.clear();
    /// propagate level by level until there is no affected tree node
    while(!bAffectedNodes.empty()) {// if bAffectedNodes is not empty, go to the next level
//        cout<<"Level "<<nodeToLevel[*bAffectedNodes.begin()]<<". Number of affected nodes: "<<bAffectedNodes.size()<<endl;
        for (auto it = bAffectedNodes.begin(); it != bAffectedNodes.end(); ++it) {//get the nodesToCheck of father level
            if(*it != 0)
                nodesToCheck.insert(GTree[*it].father);
        }
        bAffectedNodes.clear();
        //deal with father level
        int threadnum=1;
        while (!nodesToCheck.empty()) {
            /// parallel computing by boost::thread_group
            if(nodesToCheck.size()>thread_num){
                threadnum = thread_num;
            }else{
                threadnum = nodesToCheck.size();
            }
//            threadnum = nodesToCheck.size();
            unordered_map<int, bool> flag_update; flag_update.clear();
            /// multi-thread
//            cout<<"Size of nodesToCheck: "<<nodesToCheck.size()<<endl;
            boost::thread_group threadf;
            for(int i=0;i<threadnum;i++){
                //deal with the affected tree node: for each border-affected node, we need to check and update it
                nid = *nodesToCheck.begin();
                level_i = nodeToLevel[nid];
                nodesToCheck.erase(nid);
                updatedNodes.insert(nid);
//                if(nid == 0)
//                    cout<<nid<<endl;
                flag_update.insert({nid,false});
                threadf.add_thread(new boost::thread(&Gtree::NonLeafUpdate, this, nid, boost::ref(graph), boost::ref(flag_update)));
            }
            threadf.join_all();
            for(auto it=flag_update.begin();it!=flag_update.end();++it){
                nid = it->first;
                if(it->second){
                    unordered_set<int> borderSet;
                    borderSet.clear();
                    for(auto it=GTree[nid].borders.begin();it!=GTree[nid].borders.end();++it)
                        borderSet.insert(*it);
                    NonLeafNodeContract(nid, level_i, graph, nodesToCheck, updatedNodes, treenodelevel, borderSet);//contract the leaf node
//                    bAffectedNodes.insert(nid);
                }
                if(nid!=0)
                    bAffectedNodes.insert(nid);
            }
        }
    }
}
//function of getting the children nodes of offspring
void Gtree::GetChildNodes(int lca, unordered_set<int> & nodesToCheck, unordered_set<int> & leafsToCheck, unordered_set<int> & updatedLeafs){
    for(auto it=GTree[lca].children.begin();it!=GTree[lca].children.end();++it){
        if(GTree[*it].isleaf){//if it is leaf node
            if(updatedLeafs.find(*it) == updatedLeafs.end())//if not found
                leafsToCheck.insert(*it);
        }else{
            if(updatedLeafs.find(*it) == updatedLeafs.end())//if not found
                nodesToCheck.insert(*it);
            GetChildNodes(*it,nodesToCheck,leafsToCheck,updatedLeafs);
        }
    }
}
void Gtree::GetChildLeafs(int lca, unordered_set<int> & leafsToCheck, unordered_set<int> & updatedLeafs){
    for(auto it=GTree[lca].children.begin();it!=GTree[lca].children.end();++it){
        if(GTree[*it].isleaf){//if it is leaf node
            if(updatedLeafs.find(*it) == updatedLeafs.end())//if not found
                leafsToCheck.insert(*it);
        }else{
            GetChildNodes(*it,leafsToCheck,leafsToCheck,updatedLeafs);
        }
    }
}
//function of computing the distance between two vertices according to G-tree
int Gtree::ComputeDisByTree(int src, int dst, vector<TreeNode> & GTree) {
    int Ns = Nodes[src].inleaf;
    int Nt = Nodes[dst].inleaf;
    int posa = Nodes[src].inleafpos;//the position in leaf node
    int posb = Nodes[dst].inleafpos;//the position in leaf node
    auto num_border_ns = GTree[Ns].borders.size();
    auto num_leafnode_ns = GTree[Ns].leafnodes.size();
    auto num_border_nt = GTree[Nt].borders.size();
    auto num_leafnode_nt = GTree[Nt].leafnodes.size();

    //initiate cache_q and is_visited of each tree node
    init_query(GTree);

    // Find LCA index in gtreepath
    int LCA_pos = find_LCA_pos(src, dst);

    // Step out of leaf node Ns
    for (int i = 0; i < num_border_ns; ++i) {
//        cout << i*num_leafnode_ns+posa << " "<<GTree[Ns].mind.size() << endl;
        GTree[Ns].cache_q[i] = GTree[Ns].mind[i * num_leafnode_ns + posa];//store the distance from all borders to s
    }

    // Init some variables
    const auto &up_path = Nodes[src].gtreepath;
    const auto &down_path = Nodes[dst].gtreepath;
    int cn, tn, min, dist, posx, posy;
    unsigned long union_border_size;

    // Step out of nodes until meeting LCA
    // The cache_q of each node 'tn' stores the distance from vertex src to node tn's child then to tn
    for (auto i = up_path.size() - 2; i >= LCA_pos + 1; --i) {
        tn = up_path[i];
        cn = up_path[i + 1];  // child node
        union_border_size = GTree[tn].union_borders.size();
        for (int j = 0; j < GTree[tn].borders.size(); j++) {
            min = INT_MAX;
            posx = GTree[tn].current_pos[j];//get the position id of tn's borders in this node
            for (int k = 0; k < GTree[cn].borders.size(); k++) {
                posy = GTree[cn].up_pos[k];//get the position id of cn's borders in the parent node
                dist = GTree[cn].cache_q[k] + GTree[tn].mind[posx * union_border_size + posy];
                if (dist < min) {
                    min = dist;
                }
            }
            GTree[tn].cache_q[j] = min;//update the distance from the borders to its children's borders
        }
    }

    // Step across LCA (from one branch to another)
    // The cache_q of Nt's top ancestor node 'nt_top' stores the distance
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
//            dist = GTree[ns_top].cache_q[i] + GTree[lca_node].mind[nt_top_up_pos * union_border_size + ns_top_up_pos];
            dist = GTree[ns_top].cache_q[j] + GTree[lca_node].mind[nt_top_up_pos * union_border_size + ns_top_up_pos];
            if (dist < min) {
                min = dist;
            }
        }
        GTree[nt_top].cache_q[i] = min;
    }


    // Step into nodes until meeting Nt
    // The cache_q of each node 'tn' stores the distance from vertex src to node tn's parent then to tn
    for (auto i = LCA_pos + 2; i < down_path.size(); ++i) {
        tn = down_path[i];
        cn = down_path[i - 1];   // parent node
        union_border_size = GTree[cn].union_borders.size();
        for (int j = 0; j < GTree[tn].borders.size(); j++) {
            min = INT_MAX;
            posx = GTree[tn].up_pos[j];
            for (int k = 0; k < GTree[cn].borders.size(); k++) {
                posy = GTree[cn].current_pos[k];
                dist = GTree[cn].cache_q[k] + GTree[cn].mind[posy * union_border_size + posx];
                if (dist < min) {
                    min = dist;
                }
            }
            // update
            GTree[tn].cache_q[j] = min;
        }
    }
    // Step into the leaf node Nt
    min = INT_MAX;
    for (int i = 0; i < num_border_nt; ++i) {
        dist = GTree[Nt].cache_q[i] + GTree[Nt].mind[i * num_leafnode_nt + posb];
        if (dist < min) {
            min = dist;
        }
    }
    return min;
}
//Dijkstra's algorithm from point to point
int Gtree::dijkstra_p2p(int s, int t) {//inline
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
}

int Gtree::Distance_query(int src, int dst) {
    int Ns = Nodes[src].inleaf;
    int Nt = Nodes[dst].inleaf;

    // Case D1: vertices src and dst are in the same leaf node
    if (Ns == Nt) {
//        cout<<src << " " << dst << " : s and t are in the same leaf node."<<endl;
        return dijkstra_p2p(src, dst);
    }

    int posa = Nodes[src].inleafpos;//the position in leaf node
    int posb = Nodes[dst].inleafpos;//the position in leaf node
    auto num_border_ns = GTree[Ns].borders.size();
    auto num_leafnode_ns = GTree[Ns].leafnodes.size();
    auto num_border_nt = GTree[Nt].borders.size();
    auto num_leafnode_nt = GTree[Nt].leafnodes.size();

    // Case D3: there is no shortcut between Ns and Nt

    // Find LCA index in gtreepath
    int LCA_pos = find_LCA_pos(src, dst);

    // Step out of leaf node Ns
    for (int i = 0; i < num_border_ns; ++i) {
//        cout << i*num_leafnode_ns+posa << " "<<GTree[Ns].mind.size() << endl;
        GTree[Ns].cache_q[i] = GTree[Ns].mind[i * num_leafnode_ns + posa];//store the distance from all borders to s
//        cout << src << " "<< GTree[Ns].borders[i] << " " << GTree[Ns].cache_q[i] << " " << Dijkstra(src,GTree[Ns].borders[i]) << endl;
    }

    // Init some variables
    const auto &up_path = Nodes[src].gtreepath;
    const auto &down_path = Nodes[dst].gtreepath;
    int cn, tn, min, dist, posx, posy;
    unsigned long union_border_size;

    // Step out of nodes until meeting LCA
    // The cache_q of each node 'tn' stores the distance from vertex src to node tn's child then to tn
    TIME_TICK_START
    for (auto i = up_path.size() - 2; i >= LCA_pos + 1; --i) {
        tn = up_path[i];
        cn = up_path[i + 1];  // child node
        union_border_size = GTree[tn].union_borders.size();
        for (int j = 0; j < GTree[tn].borders.size(); j++) {
            min = INT_MAX;
            posx = GTree[tn].current_pos[j];//get the position id of tn's borders in this node
            for (int k = 0; k < GTree[cn].borders.size(); k++) {
                posy = GTree[cn].up_pos[k];//get the position id of cn's borders in the parent node
                dist = GTree[cn].cache_q[k] + GTree[tn].mind[posx * union_border_size + posy];
                if (dist < min) {
                    min = dist;
                }
            }
            GTree[tn].cache_q[j] = min;//update the distance from the borders to its children's borders
//            if(ifDebug){
//                int temp_dis = Dijkstra(src,GTree[tn].borders[j],Nodes);
//                if(temp_dis != GTree[tn].cache_q[j]){
//                    cout << "Branch of s incorrect: " <<src << " "<< GTree[tn].borders[j] << " " << GTree[tn].cache_q[j] << " " << temp_dis << endl;
//                }
//            }

        }
    }


    // Step across LCA (from one branch to another)
    // The cache_q of Nt's top ancestor node 'nt_top' stores the distance
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
//            dist = GTree[ns_top].cache_q[i] + GTree[lca_node].mind[nt_top_up_pos * union_border_size + ns_top_up_pos];
            dist = GTree[ns_top].cache_q[j] + GTree[lca_node].mind[nt_top_up_pos * union_border_size + ns_top_up_pos];
            if (dist < min) {
                min = dist;
            }
        }
        GTree[nt_top].cache_q[i] = min;
//        if(ifDebug){
//            int temp_dis = Dijkstra(src,GTree[nt_top].borders[i],Nodes);
//            if(temp_dis != GTree[nt_top].cache_q[i]){
//                cout << "Step across LCA incorrect: " <<src << " "<< GTree[nt_top].borders[i] << " " << GTree[nt_top].cache_q[i] << " " << temp_dis << endl;
//            }
//        }

    }


    // Step into nodes until meeting Nt
    // The cache_q of each node 'tn' stores the distance from vertex src to node tn's parent then to tn
    for (auto i = LCA_pos + 2; i < down_path.size(); ++i) {
        tn = down_path[i];
        cn = down_path[i - 1];   // parent node
        union_border_size = GTree[cn].union_borders.size();
        for (int j = 0; j < GTree[tn].borders.size(); j++) {
            min = INT_MAX;
            posx = GTree[tn].up_pos[j];
            for (int k = 0; k < GTree[cn].borders.size(); k++) {
                posy = GTree[cn].current_pos[k];
                dist = GTree[cn].cache_q[k] + GTree[cn].mind[posy * union_border_size + posx];
                if (dist < min) {
                    min = dist;
                }
            }
            // update
            GTree[tn].cache_q[j] = min;
//            if(ifDebug){
//                int temp_dis = Dijkstra(src,GTree[tn].borders[j],Nodes);
//                if(temp_dis != GTree[tn].cache_q[j]){
//                    cout << "Branch of t incorrect: " << src << " "<< GTree[tn].borders[j] << " " << GTree[tn].cache_q[j] << " " << temp_dis << endl;
//                }
//            }

        }
    }

    // Step into the leaf node Nt
    min = INT_MAX;
    for (int i = 0; i < num_border_nt; ++i) {
        dist = GTree[Nt].cache_q[i] + GTree[Nt].mind[i * num_leafnode_nt + posb];
        if (dist < min) {
            min = dist;
        }
    }
    TIME_TICK_END
    return min;
}

void Gtree::CorrectnessCheck(int times){
    srand (time(NULL));
    int s, t;
    int dis_Dijk, dis_GTree;
    Timer tt1,tt2;
    double time_Dijk=0,time_PLL=0;

    for(int i=0;i<times;i++){//times
        s=rand()%node_num;
        t=rand()%node_num;
//        s=759; t=760;
        tt1.start();
        dis_Dijk=Dijkstra(s,t,Nodes);
        tt1.stop();
        time_Dijk += tt1.GetRuntime();
        init_query(GTree);
        tt2.start();
        dis_GTree=Distance_query(s,t);
        tt2.stop();
        time_PLL += tt2.GetRuntime();
        if(dis_Dijk!=dis_GTree){
            cout<<"InCorrect!! "<<s<<" "<<t<<" "<<dis_GTree<<" "<<dis_Dijk<<endl;
        }

        //cout<<"PID "<<VtoParID[s]<<" "<<VtoParID[t]<<endl;
    }
//    cout << "Average run time of Dijkstra: "<<time_Dijk*1000/times<<" ms."<<endl;
//    cout << "Average run time of Plannar_PLL: "<<time_PLL*1000/times<<" ms."<<endl;
    //cout<<"Correctness finish!"<<endl;
}


#endif //GTREE_HPP

