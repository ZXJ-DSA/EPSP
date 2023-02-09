//
// Created by Xinjie ZHOU on 25/8/2022.
//

#ifndef GSTARTREE_H
#define GSTARTREE_H

#include <queue>
#include "gtree.h"

struct leaf_node_pair {
    int x, y;
    double v;

    leaf_node_pair(int _x, int _y, double _v) {
        x = _x;
        y = _y;
        v = _v;
    }

    bool operator<(const leaf_node_pair &rhs) const {
        return v < rhs.v;   // top item is largest
    }
};

class Gstartree:public Gtree{
public:
    unsigned long SHORTCUT_THRESHOLD = 0;
    vector<int> leaf_nodes;
    unordered_map<int, unordered_map<int, vector<int> > > shortcuts;
    vector<int> query_objects;

    void load_gtree();
//    void load_gtreeQ();
    void load_shortcuts();
//    void build_up_and_down_pos();
    double compute_value_of_leaf_node_pair(int i, int j);//inline
    bool check_leaf_adjacent(TreeNode &li, TreeNode &lj);
    pair<int,int> compute_lca_level(int x, int y);
    void build_shortcuts_for_nodes(int x, int y, int top_level_x,int top_level_y);//inline
//    void build_shortcuts_for_nodes(int x, int y);//inline
    void build_shortcuts();
    void save_shortcuts(string FS);

    int gstartree_build(); // build G*-Tree
    void init_gstarQ(bool ifGstar);
    vector<int> load_objects();
//    int dijkstra_p2p(int s, int t);
    vector<int> dijkstra_candidate(int s, unordered_set<int> &cands);
    vector<int> dijkstra_candidate(int s, vector<int> &cands);
    inline bool check_shortcut(int ns, int nt);
//    inline int find_LCA_pos(int src, int dst);
    inline int dist_query(int src, int dst, bool ifGstar);
    void dist_main(int run_times, bool ifGstar); // query based on G*-Tree
};

#endif //GSTARTREE_H
