//
// Created by Xinjie ZHOU on 25/8/2022.
//

#ifndef GSTARTREE_HPP
#define GSTARTREE_HPP

#include "gstartree.h"

extern bool ifDebug;
// load gtree index from file
void Gstartree::load_gtree() {
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
    if(fin == NULL){
        cout << "Failed to open file " << filename << endl;
        exit(1);
    }
//	FILE *fin = fopen( FILE_GTREE, "rb" );
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
        }
        // father
        fread(&father, sizeof(int), 1, fin);
        tn.father = father;

        GTree.push_back(tn);

        if (tn.isleaf) {
            leaf_nodes.push_back(node_count);
        }

        int i = node_count;

        GTree[i].is_cached = false;
        int j = i, father;
        while ((father = GTree[j].father) > 0) j = father;

        GTree[i].cache.resize(GTree[j].borders.size(), vector<int>(GTree[i].borders.size(), 0));
//        if(i== 1362){
//            cout<<j<<endl;
//            cout<<GTree[j].borders.size()<<" "<<GTree[i].borders.size()<<endl;
//            cout<<GTree[i].cache.size()<<" "<<GTree[i].cache[0].size()<<endl;
//        }

        ++node_count;

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
//    fin = fopen(FILE_NODES_GTREE_PATH, "rb");
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
        // pos increase
        pos++;
    }
    fclose(fin);
    delete[] buf;
}

// load gtree index from file
/*void Gstartree::load_gtreeQ() {
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
        strcpy(filename2, DataPath.c_str());
        strcat(filename2, dataset.c_str());
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
}*/

void Gstartree::load_shortcuts() {
    unsigned long eta;
    if(dataset == "cal"){
        eta = (unsigned long)2*noe*percent/100;
    }else{
        eta = (unsigned long)edge_num*percent/100;
    }

    char filename[200];
    if(dataset == "cal"){
        strcpy(filename, DataPath.c_str());
        strcat(filename, FILE_SHORTCUT.c_str());
        strcat(filename, "-");
        strcat(filename, std::to_string(percent).c_str());
        strcat(filename, ".bin");
    }else {
        strcpy(filename, dirname.c_str());
        strcat(filename, "/");
        strcat(filename, FILE_SHORTCUT.c_str());
        strcat(filename, "-");
        strcat(filename, std::to_string(percent).c_str());
        strcat(filename, ".bin");
    }
    FILE *fin = fopen(filename, "rb");
//    FILE *fin = fopen(FILE_SHORTCUT, "rb");
    if(fin == NULL){
        cout << "Failed to open file " << filename << endl;
        exit(1);
    }
    cout<<"Loading shortcuts... ";
    int *buf;
    int count;

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
            shortcuts[i][j].insert(shortcuts[i][j].end(), buf, buf + count);
            shortcuts[j][i].clear();
            shortcuts[j][i].emplace_back(-1);
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
    cout <<"Done."<<endl;

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
/*double Gstartree::compute_value_of_leaf_node_pair(int i, int j) {//inline
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
}*/
//new function of computing the score of shortcut pair: correct one
double Gstartree::compute_value_of_leaf_node_pair(int i, int j) {//inline
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
}
//function to check whether two leaf nodes are adjacent
bool Gstartree::check_leaf_adjacent(TreeNode &li, TreeNode &lj) {

    for (const auto &bi : li.borders) {
        for (const auto &bj : lj.borders) {
            if (find(Nodes[bi].adjnodes.begin(), Nodes[bi].adjnodes.end(), bj) != Nodes[bi].adjnodes.end()) {//if there is an adjacent in lj
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
//function of building shortcuts for given node pair
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
                GTree[i].cache.resize(GTree[j].borders.size(), vector<int>(GTree[i].borders.size(), 0));
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
                    GTree[path[top_level_x]].cache[z][zz] = GTree[path[top_level_x-1]].mind[posz * GTree[path[top_level_x-1]].union_borders.size() + poszz];//get dis from each border to other borders in level-1 node
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
                            dist = GTree[tn].cache[z][k] + GTree[tn].mind[posb * union_border_size + posa];
                            if (dist < min)
                                min = dist;
                        }

                        GTree[cn].cache[z][j] = min;
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
//                    if(z==203){
//                        cout<<GTree[path[top_level_y]].borders.size()<<endl;
//                        cout<<path[top_level_y]<<endl;
//                        cout<<GTree[path[top_level_y]].cache.size()<<" "<<GTree[path[top_level_y]].cache[z].size() <<endl;
//                        cout<<GTree[path[top_level_y-1]].mind.size()<<" "<<posz * GTree[path[top_level_y-1]].union_borders.size() + poszz<<endl;
//                    }
                    GTree[path[top_level_y]].cache[z][zz] = GTree[path[top_level_y-1]].mind[posz * GTree[path[top_level_y-1]].union_borders.size() + poszz];
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
                            dist = GTree[tn].cache[z][k] + GTree[tn].mind[posb * union_border_size + posa];
                            if (dist < min)
                                min = dist;
                        }
                        GTree[cn].cache[z][j] = min;
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
                    dist = GTree[x].cache[z][i] + GTree[y].cache[zz][j] + mid[z][zz];
                    if (dist < min) {
                        min = dist;
                    }
                }
            }
//            int temp_dis = Dijkstra(GTree[x].borders[i],GTree[y].borders[j]);
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

    for (int i = 0; i < leaf_nodes.size() - 1; ++i) {//for each leaf node
        for (int j = i + 1; j < leaf_nodes.size(); ++j) {
            if (GTree[leaf_nodes[i]].father != GTree[leaf_nodes[j]].father &&
                check_leaf_adjacent(GTree[leaf_nodes[i]], GTree[leaf_nodes[j]])) {//if they are not sibling but adjacent
                new_value = compute_value_of_leaf_node_pair(leaf_nodes[i], leaf_nodes[j]);//larger value mean higher cost
                if (value_sum < n) {
                    leaf_node_pairs.emplace(leaf_nodes[i], leaf_nodes[j], new_value);
//                    value_sum += (GTree[i].borders.size() * GTree[j].borders.size());
                    value_sum += (GTree[leaf_nodes[i]].borders.size() * GTree[leaf_nodes[j]].borders.size());
                } else if (leaf_node_pairs.top().v < new_value) {//if
                    leaf_node_pairs.pop();
                    leaf_node_pairs.emplace(leaf_nodes[i], leaf_nodes[j], new_value);
                }
            }
        }
    }

    count = leaf_node_pairs.size();

    cout << "Shortcut number of leaf-node pair: " << count << endl;

    int idx_i, idx_j;
    multimap<int,pair<int,int> > shortcut_node_pairs;
    while (!leaf_node_pairs.empty()) {
        idx_i = leaf_node_pairs.top().x;
        idx_j = leaf_node_pairs.top().y;
        shortcut_node_pairs.insert({compute_lca_level(idx_i,idx_j).first,make_pair(idx_i,idx_j)});
//        cout << "x = " << idx_i << "; y = " << idx_j << endl;
//        build_shortcuts_for_nodes(idx_i, idx_j);
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
    if(task_type == 2){
        if(dataset == "cal"){
            load_graph();
        }else{
            ReadGraph_W();
        }
    }
    /// data loading
    load_gtree();
    load_minds();
    build_up_and_down_pos();

//    cout << Dijkstra(20979,22302) << endl;

//    vector<unsigned long> eta_list = {100000, 500000, 1000000, 1500000, 3000000};
    vector<unsigned long> eta_list;//the shortcut thresholds
    if(dataset == "cal"){
        eta_list = {(unsigned long)2*noe*percent/100, (unsigned long)2*noe/10, (unsigned long)2*noe/5};
    }else{
        eta_list = {(unsigned long)edge_num*percent/100, (unsigned long)edge_num/10, (unsigned long)edge_num/5};
    }


    for (auto eta : eta_list) {
        SHORTCUT_THRESHOLD = eta;
        Timer tt;
        tt.start();
        cout << "Building shortcuts..." << endl;
        cout << "Percentage of shortcuts: "<<percent<<" %; eta: " << eta << endl;
        build_shortcuts();
        cout << "Done." << endl;
        tt.stop();
        cout << "The time for G*-tree building: " << tt.GetRuntime() << " seconds." << endl;

        char filename[200];
        if(dataset == "cal"){
            strcpy(filename, DataPath.c_str());
            strcat(filename, FILE_SHORTCUT.c_str());
            strcat(filename, "-");
            strcat(filename, std::to_string(percent).c_str());
            strcat(filename, ".bin");
        }else {
            strcpy(filename, dirname.c_str());
            strcat(filename, "/");
            strcat(filename, FILE_SHORTCUT.c_str());
            strcat(filename, "-");
            strcat(filename, std::to_string(percent).c_str());
            strcat(filename, ".bin");
        }
        cout << "Saving shortcuts..." << endl;
        save_shortcuts(filename);
//        save_shortcuts(string(argv[6]) + "-" + std::to_string(eta) + ".sc");
        cout << "Done." << endl;
        break;
    }

    return 0;
}

void Gstartree::init_gstarQ(bool ifGstar) {
    if(task_type != 1){
        if(dataset == "cal"){
            load_graph();
        }else{
            ReadGraph_W();
        }
    }
    load_gtreeQ();
    load_minds();
    if(ifGstar){
        shortcuts.clear();
        load_shortcuts();
    }
    build_up_and_down_pos();//
}

//void Gstartree::init_query() {
//    for (auto &tn: GTree) {
//        //tn.oclist.clear();
//        tn.cache_q = vector<int>(tn.borders.size(), 0);
//        tn.is_visited = false;
//    }
//}

vector<int> Gstartree::load_objects() {
    unordered_set<int> o;
    o.clear();
    char filename[200];
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

vector<int> Gstartree::dijkstra_candidate(int s, unordered_set<int> &cands) {//inline
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

}

vector<int> Gstartree::dijkstra_candidate(int s, vector<int> &cands) {//inline
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

}

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

inline int Gstartree::dist_query(int src, int dst, bool ifGstar) {
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

    // Case D2: there is a shortcut between Ns and Nt
    if(ifGstar){
        if (check_shortcut(Ns, Nt)) {
            cout<<src << " " << dst << " : there is a shortcut between Ns and Nt."<<endl;
            vector<int> sc = Ns < Nt ? shortcuts[Ns][Nt] : shortcuts[Nt][Ns];//chose the smaller one as the first dimension
            int temp_dist;
            for (int i = 0; i < num_border_ns; ++i) {
                GTree[Ns].cache_q[i] = GTree[Ns].mind[i * num_leafnode_ns + posa];
//            temp_dist=Dijkstra(src,GTree[Ns].borders[i]);
//            if(GTree[Ns].cache_q[i] != temp_dist){
//                cout << src << " "<< GTree[Ns].borders[i] << " " << GTree[Ns].cache_q[i] << " " << temp_dist << endl;
//            }
            }

            for (int i = 0; i < num_border_nt; ++i) {
                GTree[Nt].cache_q[i] = GTree[Nt].mind[i * num_leafnode_nt + posb];
//            temp_dist = Dijkstra(dst,GTree[Nt].borders[i]);
//            if(GTree[Nt].cache_q[i] != temp_dist){
//                cout << dst << " "<< GTree[Nt].borders[i] << " " << GTree[Nt].cache_q[i] << " " << temp_dist << endl;
//            }
            }

            int min = INT_MAX, dist;

            if (Ns < Nt) {
                TIME_TICK_START
                for (int i = 0; i < num_border_ns; ++i) {
                    for (int j = 0; j < num_border_nt; ++j) {
//                    dist = GTree[Ns].cache_q[i] + GTree[Nt].cache_q[j] + sc[i * num_border_ns + j];
                        dist = GTree[Ns].cache_q[i] + GTree[Nt].cache_q[j] + sc[i * num_border_nt + j];
                        if (dist < min) {
                            min = dist;
                        }
                    }
                }
                TIME_TICK_END
            } else {
                TIME_TICK_START
                for (int i = 0; i < num_border_ns; ++i) {
                    for (int j = 0; j < num_border_nt; ++j) {
                        dist = GTree[Ns].cache_q[i] + GTree[Nt].cache_q[j] +  sc[j * num_border_ns + i];
//                    temp_dist = Dijkstra(GTree[Ns].borders[i],GTree[Nt].borders[j]);
//                    if(sc[j * num_border_ns + i] != temp_dist){
//                        cout << GTree[Ns].borders[i] << " "<< GTree[Nt].borders[j] << " " << sc[j * num_border_ns + i] << " " << temp_dist << endl;
//                    }
                        if (dist < min) {
                            min = dist;
                        }
                    }
                }
                TIME_TICK_END
            }
            return min;
        }
    }

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
//            cout << src << " "<< GTree[tn].borders[j] << " " << GTree[tn].cache_q[j] << " " << Dijkstra(src,GTree[tn].borders[j]) << endl;
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
//        cout << src << " "<< GTree[nt_top].borders[i] << " " << GTree[nt_top].cache_q[i] << " " << Dijkstra(src,GTree[nt_top].borders[i]) << endl;
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
//            cout << src << " "<< GTree[tn].borders[j] << " " << GTree[tn].cache_q[j] << " " << Dijkstra(src,GTree[tn].borders[j]) << endl;
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

void Gstartree::dist_main(int run_times, bool ifGstar) {
//    FILE_NODE = argv[2];
//    FILE_EDGE = argv[3];
//
//    FILE_NODES_GTREE_PATH = argv[4];
//    FILE_GTREE = argv[5];
//    FILE_ONTREE_MIND = argv[6];
//
//    FILE_QUERY = argv[7];
//    FILE_SHORTCUT = argv[8];

    // init
    if(task_type == 2){
        init_gstarQ(ifGstar);
    }

    if(ifDebug){
        CorrectnessCheck(100);
    }

    int ID1, ID2, num;
    vector<pair<int, int>> ODpair;

    long double query_time;

    string filename;
    if(dataset == "cal"){
        filename = DataPath + FILE_QUERY;
    }else {
        filename = DataPath +  "/" + dataset + "/" + FILE_QUERY;
    }
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
    int dis_Dijk = 0;
    int dis_Gtree = 0;

    cout << "Run times: "<<run_times<<endl;

    Timer tt;
    tt.start();
    for (int i = 0; i < run_times; ++i) {//0
        ID1 = ODpair[i].first;
        ID2 = ODpair[i].second;
        init_query(GTree);
        dis_Gtree = dist_query(ID1, ID2,ifGstar);
        if(ifDebug){
            dis_Dijk = Dijkstra(ID1,ID2,Nodes);
            if(dis_Dijk != dis_Gtree){
                cout << "Distance Error!!! " << ID1 << " " << ID2 << " : " << dis_Gtree << " " << dis_Dijk << endl;
            }
        }
        ++count;
    }

//    cout<<"Done."<<endl;

    tt.stop();
    cout << "The average time for querying: " << tt.GetRuntime()*1000/run_times << " ms." << endl;


//    cout << (total_time / count) << endl;


}
#endif //GSTARTREE_HPP
