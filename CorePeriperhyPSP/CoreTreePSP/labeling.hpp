// Hub Labels are lists of hubs and distances to them attached to every vertex in a graph.
// This file contains the class to store labels.
// Available methods include making a query, write labels to file, read labels from file,
// merge label tables, clean label tables.
//
//  Author: Xinjie ZHOU

#pragma once

#include <vector>
#include <limits>
#include <cassert>
#include <algorithm>
#include <utility>
#include <fstream>
#include <istream>
#include <omp.h>
#include "Heap.h"

namespace hl {

#define INF 99999999
//typedef unsigned int vertex;
typedef int vertex;
//typedef unsigned long long ull;

// Class to store labels
class Labeling {

public:
    std::vector< std::vector< std::pair<vertex,int> > > Labels;     // Lists of forward/reverse hubs

    vertex n;
    bool *lck;

    Labeling(size_t n = 0) :
        Labels(n, std::vector< std::pair<vertex,int> >()),
        n(n),
        lck(new bool[n]()) {}

    void resize(size_t n_=0){
        Labels.assign(n_, std::vector< std::pair<vertex,int> >()),
        n=n_,
        lck=new bool[n]();
    }

    void clear(){
        clearLabel();
//        delete[] lck;
    }

    // Find u-v distance
    int query(vertex u, vertex v) {
        int r = INF;
        for (size_t i=0, j=0; i < Labels[u].size() && j < Labels[v].size();) {
            if (Labels[u][i].first == Labels[v][j].first) {
                r = std::min(r, Labels[u][i++].second + Labels[v][j++].second);
            } else if (Labels[u][i].first < Labels[v][j].first) ++i;
            else ++j;
        }
        return r;
    }
    
    bool cover(vertex u, vertex v, uint d=INF) {
        for (size_t i=0, j=0; i < Labels[u].size() && j < Labels[v].size();) {
            if (Labels[u][i].first == Labels[v][j].first) {
                if (d >= Labels[u][i++].second + Labels[v][j++].second) return true;
            } else if (Labels[u][i].first < Labels[v][j].first) ++i;
            else ++j;
        }
        return false;
    }


//    inline bool clean_cover(Vertex u, Vertex v, unsigned int f, Distance d=infty, size_t hub_order=0) {
//        for (size_t i=0, j=0; i < label_v[u][f].size() && j < label_v[v][!f].size();) {
//            if (label_v[u][f][i] >= hub_order || label_v[v][!f][j] >= hub_order) return false;
//            if (label_v[u][f][i] == label_v[v][!f][j]) {
//                if (d >= label_d[u][f][i++] + label_d[v][!f][j++]) return true;
//            }
//            else if (label_v[u][f][i] < label_v[v][!f][j]) ++i;
//            else ++j;
//        }
//        return false;
//    }
//
//    inline void clean_roots(Vertex v, std::vector<Vertex> &order, unsigned int side)
//    {
//        std::vector<Vertex> temp_v;
//        std::vector<Distance> temp_d;
//        for (size_t i=0; i<label_v[v][side].size(); i++)
//        {
//            size_t hub_order = label_v[v][side][i];
//            Vertex hub =  order[hub_order];
//            Distance hub_dist = label_d[v][side][i];
//            if (!clean_cover(hub, v, side, hub_dist, hub_order))//;
//            {
//                temp_v.push_back(hub_order);
//                temp_d.push_back(hub_dist);
//            }
//        }
//        temp_v.swap(label_v[v][side]);
//        temp_d.swap(label_d[v][side]);
//    }

    // Add hub (v,d) to forward or reverse label of u
    inline void add(vertex u, vertex v, uint d) {
        
        while (!__sync_bool_compare_and_swap(&lck[u], false, true)) {} 
        Labels[u].emplace_back(v,d);
        lck[u]=false;
    }
    void add_lockfree(vertex u, vertex v, uint d) {
        Labels[u].emplace_back(v,d);
    }

    // Get labels
//    std::vector< std::vector<Vertex> > &get_label_hubs(Vertex u) { return label_v[u]; }
//    std::vector< std::vector<Distance> > &get_label_distances(Vertex u) { return label_d[u]; }

    // Get maximum label size
    size_t get_max() const {
        size_t max = 0;
        for (vertex v = 0; v < n; ++v){
            max = std::max(max, Labels[v].size());
        }

        return max;

        //size_t maxVal = 0;
        ////#pragma omp parallel for num_threads (NUM_THREAD) reduction (max: maxVal)
        ////{
        //    for (Vertex v=0; v<n; v++)
        //    {
        //        for (unsigned side = 0; side < 2; side++)
        //        {
        //            if (label_v[v][side].size() < maxVal)
        //                maxVal = label_v[v][side].size;
        //        }
        //    }
        //}
        //return maxVal;



    }

    // Get average label size
    double get_avg() const {
        long long total = 0;
        for (vertex v = 0; v < n; ++v)
            total += Labels[v].size() + Labels[v].size();
        return static_cast<double>(total)/n/2;
    }
    long long get_total() const {
        long long total = 0;
        for (vertex v = 0; v < n; ++v)
            total += Labels[v].size() + Labels[v].size();
        return total;
    }

    // Write labels to file
    bool write(char *filename) {
        std::ofstream file;
        file.open(filename);
        file << n << std::endl;
        for (vertex v = 0; v < n; ++v) {
            file << Labels[v].size();
            for (size_t i = 0; i < Labels[v].size(); ++i) {
                file << " " << Labels[v][i].first;
                file << " " << Labels[v][i].second;
            }
            file << std::endl;
        }
        file.close();
        return file.good();
    }

    // Read labels from file
    bool read(char *filename, vertex check_n = 0) {
        std::ifstream file;
        file.open(filename);
        file >> n;
        if (check_n && n != check_n) return false;
        Labels.resize(n, std::vector< std::pair<vertex,int> >());

        for (vertex v = 0; v < n; ++v) {
            size_t s;
            file >> s;
            Labels[v].resize(s);
            Labels[v].resize(s);
            for (size_t i = 0; i < s; ++i) {
                file >> Labels[v][i].first;
                file >> Labels[v][i].second;
            }
        }
        file >> std::ws;
        file.close();
        return file.eof() && !file.fail();
    }

    // Clear labels
    void clearLabel(unsigned NUM_THREAD=72) {
        #pragma omp parallel for num_threads(NUM_THREAD)
        for (vertex v = 0; v < n; ++v) {
            Labels[v].clear();
        }
    }

    // Sort labels before making queries
    void sort(unsigned NUM_THREAD) {

        // std::vector<std::vector<std::pair<Vertex, Distance>>> label(NUM_THREAD)
        // maxSize = get_max();
        // for (int i=0; i<NUM_THREAD; i++)
        //     label[i] = std::vector<std::pair<Vertex, Distance>>.reserve(maxSize)
        //#pragma omp parallel for num_threads(NUM_THREAD) schedule (dynamic)
        //for (Vertex v = 0; v < n; ++v) {
        //    for (int side = 0; side < 2; ++side) {
        //        for (size_t i = 0; i < label_v[v][side].size(); ++i)
        //            label[i].push_back(std::make_pair(label_v[v][side][i], label_d[v][side][i]));
        //        std::sort(label.begin(),label.end());
        //        for (size_t i = 0; i < label_v[v][side].size(); ++i) {
        //            label_v[v][side][i] = label[i].first;
        //            label_d[v][side][i] = label[i].second;
        //        }
        //        label[i].clear();
        //    }
        //}
        Timer tt;
        tt.start();
        #pragma omp parallel for num_threads(NUM_THREAD) schedule (dynamic, NUM_THREAD) 
        for (vertex v = 0; v < n; ++v) {
            std::sort(Labels[v].begin(),Labels[v].end());
        }
        tt.stop();
        std::cout<<"Time for label sorting: "<<tt.GetRuntime()<<" s."<<std::endl;
    }

    void postProcess(std::vector<std::unordered_map<vertex,int>>& Label, int ThreadNum){
        Timer tt;
        tt.start();
        Label.clear();
        Label.assign(n,std::unordered_map<vertex,int>());

        #pragma omp parallel for num_threads(ThreadNum) schedule (dynamic, ThreadNum)
        for(int id=0;id<n;++id){
            vertex id2;
            uint dis;
            for(auto it=Labels[id].begin();it!=Labels[id].end();++it){
                id2=it->first; dis=it->second;
                Label[id].insert({id2,dis});
            }
//            Labels[id].clear();
        }
//        Labels.clear();
        tt.stop();
        std::cout<<"Time for label post-processing: "<<tt.GetRuntime()<<" s."<<std::endl;
    }



};

struct Node{
    int ID;
    omp_lock_t lock;
};

// Class to store pruning point records
class PPR{

public:
    std::vector< std::vector< std::pair<vertex,vertex> > > PPRs;     // Lists of pruning point records
//    std::vector< std::unordered_map< vertex, std::unordered_set<vertex> > > PPRSet; // Hash table-based
    std::vector<Semaphore*> vSm;//locks of boost

    vertex n;
    bool *lck;

    PPR(size_t n = 0) : n(n), lck(new bool[n]()) {
        PPRs.assign(n, std::vector< std::pair<vertex,vertex> >());
    }

    void resize(size_t n_=0){
        PPRs.assign(n_, std::vector< std::pair<vertex,vertex> >()),
                n=n_,
                lck=new bool[n]();
    }

    void clear(){
        PPRs.clear();
//        delete[] lck;
        vSm.clear();
    }

    // Add hub (v,d) to forward or reverse label of u
    inline void add(vertex u, vertex v, uint d) {

        while (!__sync_bool_compare_and_swap(&lck[u], false, true)) {}
        PPRs[u].emplace_back(v,d);
        lck[u]=false;
    }
    void add_lockfree(vertex u, vertex v, uint d) {
        PPRs[u].emplace_back(v,d);
    }

    void write(std::string filename){
        std::ofstream OF2(filename);
        if(!OF2){
            std::cout<<"Cannot open "<<filename<<std::endl;
            exit(1);
        }
        std::cout<<"Write PPR..."<<std::endl;
//        Timer tt;
//        tt.start();
        int hub,ID2;
        for(int ID1=0;ID1<n;ID1++){
            if(!PPRs[ID1].empty()){//Order
                OF2<<ID1;
                for(auto it=PPRs[ID1].begin();it!=PPRs[ID1].end();++it){//Order
                    hub=it->first; ID2=it->second;
                    OF2<<" "<<hub<<" "<<ID2;
                }
                OF2<<std::endl;
            }

        }
        OF2.close();
//        tt.stop();
//        std::cout<<"Done."<<std::endl;
    }

    void pprInsert(std::vector<vertex>& p,std::vector<std::unordered_map<vertex,std::unordered_set<vertex>>>& PruningPointSet, std::vector<std::unordered_map<vertex,vertex>>& PruningPointSet2){
        for(auto it=p.begin();it!=p.end();++it){
            vertex id=*it;
            vertex hub,id2;
            for(auto it=PPRs[id].begin();it!=PPRs[id].end();++it){
                hub=it->first; id2=it->second;

                vSm[id]->wait();
                PruningPointSet[id][hub].insert(id2);
                PruningPointSet2[id][id2]=hub;
                vSm[id]->notify();

                vSm[id2]->wait();
                PruningPointSet[id2][hub].insert(id);
                PruningPointSet2[id2][id]=hub;
                vSm[id2]->notify();
            }
            PPRs[id].clear();
        }

    }

    // post-processing
    void postProcess(std::vector<std::unordered_map<vertex,std::unordered_set<vertex>>>& PruningPointSet, std::vector<std::unordered_map<vertex,vertex>>& PruningPointSet2, std::vector<vertex>& vertices){
        Timer tt;
        tt.start();
        PruningPointSet.assign(n,std::unordered_map<vertex,std::unordered_set<vertex>>());
        PruningPointSet2.assign(n,std::unordered_map<vertex,vertex>());

        vertex id,hub,id2;
        for(int i=0;i<vertices.size();++i){
            id=vertices[i];

            for(auto it=PPRs[id].begin();it!=PPRs[id].end();++it){
                hub=it->first; id2=it->second;

                PruningPointSet[id][hub].insert(id2);
                PruningPointSet2[id][id2]=hub;

                PruningPointSet[id2][hub].insert(id);
                PruningPointSet2[id2][id]=hub;
            }
            PPRs[id].clear();
        }

        PPRs.clear();

        tt.stop();
        std::cout<<"Time for PPR post-processing: "<<tt.GetRuntime()<<" s."<<std::endl;
    }

    void postProcess(std::vector<std::unordered_map<vertex,std::unordered_set<vertex>>>& PruningPointSet, std::vector<std::unordered_map<vertex,vertex>>& PruningPointSet2, std::vector<vertex>& vertices, int ThreadNum){
        Timer tt;
        tt.start();
        PruningPointSet.assign(n,std::unordered_map<vertex,std::unordered_set<vertex>>());
        PruningPointSet2.assign(n,std::unordered_map<vertex,vertex>());

        /// boost-based implementation
        vSm.reserve(n);
        for(int i = 0; i < n; i++)
        {
            Semaphore* s = new Semaphore(1);
            vSm.push_back(s);
        }
        std::vector<std::vector<vertex>> processID(ThreadNum);
        int pid;
        for(int i=0;i<vertices.size();++i){
            pid=i%ThreadNum;
            processID[pid].push_back(vertices[i]);
        }
        boost::thread_group thread;
        for(int i=0;i<processID.size();i++){
            thread.add_thread(new boost::thread(&PPR::pprInsert, this, boost::ref(processID[i]), boost::ref(PruningPointSet), boost::ref(PruningPointSet2)));
        }
        thread.join_all();

//        /// multiple threads: openmp-based implementation
//        std::vector<Node> vs(n);
//        for(int i=0;i<n;++i){
//            vs[i].ID=i;
//            omp_init_lock(&vs[i].lock);
//        }
//        #pragma omp parallel for num_threads(ThreadNum) schedule (dynamic, ThreadNum)
//        for(int id=0;id<n;++id){
//            vertex hub,id2;
//            for(auto it=PPRs[id].begin();it!=PPRs[id].end();++it){
//                hub=it->first; id2=it->second;
//
//                omp_set_lock(&vs[id].lock);
//                PruningPointSet[id][hub].insert(id2);
//                PruningPointSet2[id][id2]=hub;
//                omp_unset_lock(&vs[id].lock);
//
//                omp_set_lock(&vs[id2].lock);
//                PruningPointSet[id2][hub].insert(id);
//                PruningPointSet2[id2][id]=hub;
//                omp_unset_lock(&vs[id2].lock);
//            }
//            PPRs[id].clear();
//        }
//        for (int i = 0; i < n; ++i) {
//            omp_destroy_lock(&vs[i].lock);// Destroy locks
//        }



        PPRs.clear();

        tt.stop();
        std::cout<<"Time for PPR post-processing: "<<tt.GetRuntime()<<" s."<<std::endl;
    }


};

}
