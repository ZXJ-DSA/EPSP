#ifndef BATCHHLWW_HPP_
#define BATCHHLWW_HPP_

#include <sys/time.h>
#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <sstream>
#include <fstream>

#include "BatchHLWW.h"
#include "two_layer_queue.h"

//
// Implementation for weighted and undirected graphs.
//

using namespace std;

void HighwayLabellingWW::allocate() {
    //initiate original distance value from any vertex to landmarks
    distances.assign(V, vector<int>(K,INF));
    distances_1.assign(V, vector<int>(K,INF));

    highway.assign(K,vector<int>(K,INF));
    highway_1.assign(K,vector<int>(K,INF));
}

void HighwayLabellingWW::ReadGraph(string filename, int k) {
    V = 0; E = 0; K = k;

    std::ifstream ifs(filename);
    if (ifs.is_open()) {

        int ID1, ID2, weight; std::string query;
        std::unordered_map<int, int> vertex2id;//map vertex id

        if(getline(ifs, query)){
            std::istringstream iss(query);
            iss >> V >> E;//input edge
            cout<<"Node number: "<<V<<" ; Edge number: "<<E<<endl;
        }
        adj.assign(V,unordered_map<int,int>());

        if(dataset == "FRIEND"){
            cout<<"Dataset: "<<dataset<<endl;
        }else {
            cout<<"Dataset: "<<dataset<<endl;
        }

        bool flag= true;
        while (getline(ifs, query)){
            std::istringstream iss(query);
            iss >> ID1 >> ID2 >> weight;//input edge

            if (ID1 != ID2) {
                if(dataset == "FRIEND"){
                    adj[ID1][ID2] = weight;
                    adj[ID2][ID1] = weight;
                    if(flag){
                        cout<<"Flag 1"<<endl;
                        flag = false;
                    }
                }else{
//        adj[ID1][ID2] = 1;
                    adj[ID1][ID2] = weight;
//        adj[ID2].push_back(ID1);
                }

            }
        }
        ifs.close();

//    for (int v = 0; v < V; v++) {
//      std::sort(adj[v].begin(), adj[v].end());
//      adj[v].erase(std::unique(adj[v].begin(), adj[v].end()), adj[v].end());
//    }
        long long int edgenum=0;
        for(int i = 0; i < V; i++)
            edgenum += adj[i].size();
        if(edgenum != E){
            cout<<"Edge number is wrong!"<<endl;
            E = edgenum;
            std::cout << "V : " << V << " E : " << E << std::endl << std::endl;
            exit(1);
        }
    } else
        std::cout << "Unable to open file" << std::endl;

//  adjOld = adj;
}
//Correct version
/*HighwayLabellingWW::HighwayLabellingWW(std::string filename, int k) {
  V = 0; E = 0; K = k;

  std::ifstream ifs(filename);
  if (ifs.is_open()) {

    int ID1, ID2, weight; std::string query;
    std::unordered_map<int, int> vertex2id;//map vertex id

    if(getline(ifs, query)){
        std::istringstream iss(query);
        iss >> V >> E;//input edge
        cout<<"Node number: "<<V<<" ; Edge number: "<<E<<endl;
    }
    adj.assign(V,unordered_map<int,int>());

      if(dataset == "FRIEND"){
          cout<<"Correct! Dataset: "<<dataset<<endl;
      }else {
          cout<<"Incorrect! Dataset: "<<dataset<<endl;
      }

    bool flag= true;
    while (getline(ifs, query)){
      std::istringstream iss(query);
      iss >> ID1 >> ID2 >> weight;//input edge

      if (ID1 != ID2) {
          if(dataset == "FRIEND"){
              adj[ID1][ID2] = weight;
              adj[ID2][ID1] = weight;
              if(flag){
                  cout<<"Flag 1"<<endl;
                  flag = false;
              }
          }else{
//        adj[ID1][ID2] = 1;
              adj[ID1][ID2] = weight;
//        adj[ID2].push_back(ID1);
          }

      }
    }
    ifs.close();

//    for (int v = 0; v < V; v++) {
//      std::sort(adj[v].begin(), adj[v].end());
//      adj[v].erase(std::unique(adj[v].begin(), adj[v].end()), adj[v].end());
//    }
    long long int edgenum=0;
    for(int i = 0; i < V; i++)
        edgenum += adj[i].size();
    if(edgenum != E){
        cout<<"Edge number is wrong!"<<endl;
        E = edgenum;
        std::cout << "V : " << V << " E : " << E << std::endl << std::endl;
        exit(1);
    }
  } else
      std::cout << "Unable to open file" << std::endl;

//  adjOld = adj;
}*/

int HighwayLabellingWW::LabellingSize() {
  long size = 0;
  for (int i = 0; i < V; i++) {
    for (int j = 0; j < K; j++) {
      if(distances_1[i][j] != INF)
        size++;
    }
  }

  return (V + 2 *size) / (1024 * 1024);
}

void HighwayLabellingWW::ConstructHighwayLabellingNewP(pair<int,int> vi, vector<int> & topk) {
    for(int i=vi.first;i<vi.second;++i){
        ConstructHighwayLabellingNew(i,topk);
    }
}

void HighwayLabellingWW::ConstructHighwayLabellingNew(int i, vector<int> & topk) {
    int *P = new int[V];
    for(int j = 0; j < V; j++)
        P[j] = INF;

    vector<int> distance(V, INF);
    vector<bool> closed(V, false);//closed means the shortest distance has been calculated
    vector<bool> colored(V, false);//colored means the shortest path from curID to v passes through one landmark
    int topNodeID, topNodeDis, NNodeID, NNodeWei;

    benchmark::heap<2, int, int> Q(V);
    Q.update(topk[i],0);
    distance[topk[i]]=0;

    while(!Q.empty()){
        Q.extract_min(topNodeID, topNodeDis);
        closed[topNodeID]=true;
        P[topNodeID] = topNodeDis;

        if(landmarks.find(topNodeID)!=landmarks.end()){//if it is a landmark
            if(topNodeID!=topk[i]){
                colored[topNodeID]=true;
            }
        }

        for(auto it=adj[topNodeID].begin();it!=adj[topNodeID].end();it++){
            NNodeID = it->first; NNodeWei = it->second;
//            NNodeID=(*it).first; NNodeWei=(*it).second;
            if(!closed[NNodeID] && distance[NNodeID]>topNodeDis+NNodeWei){
                distance[NNodeID]=topNodeDis+NNodeWei;
                Q.update(NNodeID, topNodeDis+NNodeWei);
                if(colored[topNodeID])
                    colored[NNodeID]=true;
                else{
                    colored[NNodeID]=false;
                }
            }
        }
    }

    for(int j=0;j<V;++j){
        if(!colored[j]){
            distances[j][i] = P[j];
            distances_1[j][i] = P[j];
        }
        if(P[j] == INF){
            cout<<"INF: "<<i<<" "<<j<<endl;
        }
    }

    for(int j = 0; j < K; j++) {
        if(P[topk[j]] == INF){
            cout<<"Landmark INF: "<<i<<" "<<j<<endl;
        }
        highway[i][j] = P[topk[j]];
        highway_1[i][j] = P[topk[j]];
    }

    delete [] P;
}

void HighwayLabellingWW::ConstructHighwayLabelling(int i, vector<int> & topk) {

  int *P = new int[V];
  for(int j = 0; j < V; j++)
    P[j] = INF;

  std::queue<int> que[2];//first in first out, two queues: queue 1 stores the unaffected vertex while queue 2 stores the affected vertex.
  que[0].push(topk[i]); que[0].push(-1);
  distances[topk[i]][i] = 0; distances_1[topk[i]][i] = 0; P[topk[i]] = 0; int use = 0;
  while (!que[0].empty()) {
    int u = que[use].front();
    que[use].pop();

    if(u == -1) {//used to switch between queue 1 and queue 2
      use = 1 - use;
      que[use].push(-1);
      continue;
    }

    for (auto w : adj[u]) {
      if (P[w.first] == INF) {
        P[w.first] = P[u] + 1;
        if(use == 1 || landmarks.count(w.first) > 0)//for landmarks
          que[1].push(w.first);
        else {
          que[0].push(w.first);
          distances[w.first][i] = P[w.first];
          distances_1[w.first][i] = P[w.first];
        }
      }
    }
  }
  for(int j=0;j<V;++j){
      if(P[j] == INF){
          cout<<"INF: "<<i<<" "<<j<<endl;
      }
  }

  for(int j = 0; j < K; j++) {
      if(P[topk[j]] == INF){
          cout<<"Landmark INF: "<<i<<" "<<j<<endl;
      }
    highway[i][j] = P[topk[j]];
    highway_1[i][j] = P[topk[j]];
  }

  delete [] P;
}

void HighwayLabellingWW::BuildIndex(vector<int> & topk, bool ifParallel) {

  allocate();//allocate memory
  for(int i = 0; i < K; i++)
    landmarks[topk[i]] = i;//map from landmark vertex id to simple id

  // Start computing Highway Labelling (HL)
  time_ = -GetCurrentTimeSec();
  Timer tt;
  tt.start();
  if(ifParallel){
      /// multi thread
      //multiple thread
      if(K>threadnum){
          int step=K/threadnum;
          boost::thread_group thread;
          for(int i=0;i<threadnum;i++){
              pair<int,int> p;
              p.first=i*step;
              if(i==threadnum-1)
                  p.second=K;
              else
                  p.second=(i+1)*step;
              thread.add_thread(new boost::thread(&HighwayLabellingWW::ConstructHighwayLabellingNewP, this, p, boost::ref(topk)));
          }
          thread.join_all();
      }else{
          boost::thread_group thread;
          for(int i=0;i<K;i++){
              thread.add_thread(new boost::thread(&HighwayLabellingWW::ConstructHighwayLabellingNew, this, i, boost::ref(topk)));
          }
          thread.join_all();
      }

  }else{
      /// single thread
      for (int i = 0; i < K; i++){//for each landmark
          ConstructHighwayLabellingNew(i, topk);
//      ConstructHighwayLabelling(i, topk);
      }
  }

    tt.stop();
  time_ += GetCurrentTimeSec();
    std::cout << "Construction Time : " << tt.GetRuntime() << " s. Labelling Size: " << LabellingSize() << " MB" << std::endl;
//  std::cout << "Construction Time (sec.): " << time_ << " Labelling Size: " << LabellingSize() << " MB" << std::endl;
}

void HighwayLabellingWW::UpdateLabelling(std::string filename, int m) {
    adjOldO = adj;
  std::ifstream ifs(filename);
  int a, b; std::string op;
  std::vector<std::pair<std::string, std::pair<int, int> > > updates;

  while(ifs >> op >> a >> b) {
    if(op == "EI") {
      adj[a][b]=1;
      adj[b][a]=1;
    } else if(op == "ED") {
        if(adj[a].find(b)!=adj[a].end())
            adj[a].erase(b);
        if(adj[b].find(a)!=adj[b].end())
            adj[b].erase(a);
    }
    updates.emplace_back(op, std::make_pair(a, b));
  }
  ifs.close();

  adjOld = adj;

  Timer tt;
  tt.start();

  if(m == 0){
//      cout<<"Update method: BHL+"<<endl;
      BHL_Plus(updates);
  }else if (m == 1){
//      cout<<"Update method: BHL"<<endl;
//      BHL(updates);
      BHLW(updates);
  }

  tt.stop();

  for(a = 0; a < V; a++) {
    for(b = 0; b < K; b++)
      distances_1[a][b] = distances[a][b];
  }

  for(a = 0; a < K; a++) {
    for(b = 0; b < K; b++) {
      highway_1[a][b] = highway[a][b];
    }
  }
  /// update of vertices
  verticesUpdate();

  std::cout << "Batch Update Time (sec.): " << tt.GetRuntime() << " s. Updated Labelling Size: " << LabellingSize() << " MB" << std::endl;
}
//function of index updating for weighted graph
void HighwayLabellingWW::UpdateLabellingW(std::string filename, int m, int updateType, int updateBatch, int updateVolume, bool ifParallel) {
    if(ifDebug){
        adjOldO = adj;
    }
    /// read updates
    ReadUpdates(filename+".update");
    int ID1,ID2,oldW,newW;
    double runtime = 0;
//    cout<<"Update batch: "<<updateBatch<<" ; Update volume: "<<updateVolume<<endl;
    Timer tt;
    for(int batch_i=0;batch_i<updateBatch*updateVolume;batch_i+=updateVolume) {
        if(ifDebug){
            adj = adjOldO;
            cout<<"Batch "<<batch_i<<endl;
        }
        for (int i = batch_i; i < batch_i + updateVolume; ++i) {//for each edge update
            ID1 = updateEdges[i].first.first;
            ID2 = updateEdges[i].first.second;
            oldW = updateEdges[i].second;
            switch (updateType) {
                case INCREASE:
                    newW = (int) (2 * oldW);//update
                    break;
                case DECREASE:
                    newW = (int) (0.5 * oldW);
                    if(newW < 1){
                        cout<<"Decrease wrong! "<<newW<<endl;
                        exit(1);
                    }
                    break;
                case MIX:
                    newW = (int) (2 * oldW * (rand() / double(RAND_MAX)));
                    break;
                default:
                    cout << "Wrong update type!" << endl;
                    break;
            }
            if(newW <1){
                cout<<"New edge weight is too small! "<<newW<<endl;
                exit(1);
            }
            //update edge weight
            assert(adj[ID1][ID2] == oldW);
            assert(adj[ID2][ID1] == oldW);
            adj[ID1][ID2] = newW;
            adj[ID2][ID1] = newW;
        }
        /// update
        adjOld = adj;
        tt.start();
        if(updateType == INCREASE){//for increase update
            if(ifParallel){
                /// multi thread
                if(K>threadnum){
                    //multiple thread
                    int step=K/threadnum;
                    boost::thread_group thread;
                    for(int i=0;i<threadnum;i++){
                        pair<int,int> p;
                        p.first=i*step;
                        if(i==threadnum-1)
                            p.second=K;
                        else
                            p.second=(i+1)*step;
                        thread.add_thread(new boost::thread(&HighwayLabellingWW::ConstructHighwayLabellingNewP, this, p, boost::ref(landmarkTopK)));
                    }
                    thread.join_all();
                }else{
                    boost::thread_group thread;
                    for(int i=0;i<K;i++){
                        thread.add_thread(new boost::thread(&HighwayLabellingWW::ConstructHighwayLabellingNew, this, i, boost::ref(landmarkTopK)));
                    }
                    thread.join_all();
                }

            }else{
                for (int i = 0; i < K; i++){//for each landmark
                    ConstructHighwayLabellingNew(i, landmarkTopK);
                }
            }

        }else if(updateType == DECREASE){//for decrease update
            if(m == 0){
//                cout<<"Update method: BHL+"<<endl;
//            BHL_Plus(updates);
            }else if (m == 1){
//                cout<<"Update method: BHL"<<endl;
                BHLWP(updateEdges,ifParallel);
//                BHLW(updateEdges);
            }
        }

        tt.stop();
        runtime += tt.GetRuntime();
        if(ifDebug){
            adjOldO = adj;
            /// update of vertices
            verticesUpdate();
            RemoveLandmarks(landmarkTopK);
            CorrectnessCheck(20,true);
        }
    }
    std::ifstream ifs(filename);
    int a, b; std::string op;

    for(a = 0; a < V; a++) {
        for(b = 0; b < K; b++)
            distances_1[a][b] = distances[a][b];
    }

    for(a = 0; a < K; a++) {
        for(b = 0; b < K; b++) {
            highway_1[a][b] = highway[a][b];
        }
    }
    /// update of vertices
    verticesUpdate();

    cout << "Average Batch Update Time: " << runtime/updateBatch << " s. Updated Labelling Size: " << LabellingSize() << " MB" << endl;
}

int HighwayLabellingWW::ldPair(int d, bool landmark) { return landmark ? d << 1 : (d << 1) | 1; }
int HighwayLabellingWW::ldDist(int ld) { return ld >> 1; }
int HighwayLabellingWW::ldAdd(int ld, bool landmark) { return landmark ? (ld + 2) & ~1 : ld + 2; }

/** BHL^+ **/
//original correct version
void HighwayLabellingWW::BHL_Plus(std::vector<std::pair<std::string, std::pair<int, int> > > & updates) {
  int *A = new int[V];
  int *temp = new int[V];
  int *AFF_VERTS = new int[V]; int c; int br, beta;

  for(int i = 0; i < V; i++)
    A[i] = INF;

  for(int i = 0; i < K; i++) {

    // computing affected vertices
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int> >, std::greater<std::pair<int, int> > > que;
    std::queue<std::pair<int, int> > que1; c = 0;
    for(std::pair<std::string, std::pair<int, int> > iter : updates) {

      bool e = false;
      if(iter.first == "ED") e = true;

      temp[iter.second.first] = query(i, iter.second.first);
      temp[iter.second.second] = query(i, iter.second.second);

      if(temp[iter.second.first] > temp[iter.second.second]) {
        br = ldPair(ldAdd(ldPair(temp[iter.second.second], (distances_1[iter.second.second][i] == INF)), (landmarks.count(iter.second.first) > 0)), e); beta = ldPair(ldPair(temp[iter.second.first], (distances_1[iter.second.first][i] == INF)), true);
        if(br <= beta)
          que.push(std::make_pair(br, iter.second.first));
      } else if(temp[iter.second.first] < temp[iter.second.second]) {
        br = ldPair(ldAdd(ldPair(temp[iter.second.first], (distances_1[iter.second.first][i] == INF)), (landmarks.count(iter.second.second) > 0)), e); beta = ldPair(ldPair(temp[iter.second.second], (distances_1[iter.second.second][i] == INF)), true);
        if(br <= beta)
          que.push(std::make_pair(br, iter.second.second));
      }
    }

    while (!que1.empty() || !que.empty()) {

      std::pair<int, int> p;
      if(!que1.empty() && !que.empty()) {
        if(que.top().first < que1.front().first) {
          p = que.top(); que.pop();
        } else {
          p = que1.front(); que1.pop();
        }
      } else if(!que1.empty()) {
        p = que1.front(); que1.pop();
      } else {
        p = que.top(); que.pop();
      }

      if(A[p.second] == INF) {

        for (auto temp_ : adj[p.second]) {
            int w = temp_.first;
          temp[w] = query(i, w);

      br = ldPair(ldAdd(ldDist(p.first), (landmarks.count(w) > 0)), (p.first & 1)==0?true:false); beta = ldPair(ldPair(temp[w], (distances_1[w][i] == INF)), true);
          if(br <= beta)
            que1.push(std::make_pair(br, w));
        }

        A[p.second] = p.first;
        AFF_VERTS[c] = p.second; c++;
      }
    }

    // computing boundary vertices
    std::vector<std::pair<int, int> > V_aff;
    for(int j = 0; j < c; j++) {
      temp[AFF_VERTS[j]] = INF;
      for (auto temp_ : adj[AFF_VERTS[j]]) {
          int w = temp_.first;
        if(A[w] == INF)
          temp[AFF_VERTS[j]] = min(temp[AFF_VERTS[j]], temp[w] + 1);
      }

      if(temp[AFF_VERTS[j]] != INF)
        V_aff.push_back(std::make_pair(temp[AFF_VERTS[j]], AFF_VERTS[j]));
      landmarks.count(AFF_VERTS[j])>0?highway[i][landmarks[AFF_VERTS[j]]] = INF:distances[AFF_VERTS[j]][i] = INF;
    }

    // updating the labelling
    if(V_aff.size() > 0) {
      std::sort(V_aff.begin(), V_aff.end());

      std::queue<int> quee[2]; int use = 0;
      int d = V_aff[0].first; int x = 0;
      while(x < V_aff.size() && V_aff[x].first == d) {
        if(prunable(i, V_aff[x].second, temp, A)){
          quee[1].push(V_aff[x].second);
        } else {
          distances[V_aff[x].second][i] = temp[V_aff[x].second];
          quee[0].push(V_aff[x].second);
        }
        A[V_aff[x].second] = INF;
        x++;
      }

      quee[use].push(-1);
      while (!quee[0].empty() || x < V_aff.size()) {
        int u = quee[use].front();
        quee[use].pop();

        if(u == -1) {
          if(use == 0) { d++;
            while(x < V_aff.size() && V_aff[x].first == d) {
              if(A[V_aff[x].second] != INF) {
                if(prunable(i, V_aff[x].second, temp, A)) {
                  quee[1].push(V_aff[x].second);
                } else {
                  distances[V_aff[x].second][i] = temp[V_aff[x].second];
                  quee[0].push(V_aff[x].second);
                }
                A[V_aff[x].second] = INF;
              }
              x++;
            }
          }
          use = 1 - use;
          quee[use].push(-1);
          continue;
        }

        for (auto temp_ : adj[u]) {
            int w = temp_.first;
          if(A[w] != INF) {
            temp[w] = temp[u] + 1;
            if(use == 1 || prunable(i, w, temp, A)) {
              quee[1].push(w);
            } else {
              distances[w][i] = temp[w];
              quee[0].push(w);
            }
            A[w] = INF;
          }
        }
      }
    }

    if(i == K - 1)
      continue;

    for(int j = 0; j < c; j++) {
      if(A[AFF_VERTS[j]] != INF)
        A[AFF_VERTS[j]] = INF;
    }
  }

  delete [] A;
  delete [] temp;
  delete [] AFF_VERTS;
}
//modified version
void HighwayLabellingWW::BHL_PlusW(std::vector<std::pair<std::string, std::pair<int, int> > > & updates) {
    int *A = new int[V];
    int *temp = new int[V];
    int *AFF_VERTS = new int[V]; int c; int br, beta;

    for(int i = 0; i < V; i++)
        A[i] = INF;

    for(int i = 0; i < K; i++) {

        // computing affected vertices
        std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int> >, std::greater<std::pair<int, int> > > que;
        std::queue<std::pair<int, int> > que1; c = 0;
        for(std::pair<std::string, std::pair<int, int> > iter : updates) {

            bool e = false;
            if(iter.first == "ED") e = true;

            temp[iter.second.first] = query(i, iter.second.first);
            temp[iter.second.second] = query(i, iter.second.second);

            if(temp[iter.second.first] > temp[iter.second.second]) {
                br = ldPair(ldAdd(ldPair(temp[iter.second.second], (distances_1[iter.second.second][i] == INF)), (landmarks.count(iter.second.first) > 0)), e); beta = ldPair(ldPair(temp[iter.second.first], (distances_1[iter.second.first][i] == INF)), true);
                if(br <= beta)
                    que.push(std::make_pair(br, iter.second.first));
            } else if(temp[iter.second.first] < temp[iter.second.second]) {
                br = ldPair(ldAdd(ldPair(temp[iter.second.first], (distances_1[iter.second.first][i] == INF)), (landmarks.count(iter.second.second) > 0)), e); beta = ldPair(ldPair(temp[iter.second.second], (distances_1[iter.second.second][i] == INF)), true);
                if(br <= beta)
                    que.push(std::make_pair(br, iter.second.second));
            }
        }

        while (!que1.empty() || !que.empty()) {

            std::pair<int, int> p;
            if(!que1.empty() && !que.empty()) {
                if(que.top().first < que1.front().first) {
                    p = que.top(); que.pop();
                } else {
                    p = que1.front(); que1.pop();
                }
            } else if(!que1.empty()) {
                p = que1.front(); que1.pop();
            } else {
                p = que.top(); que.pop();
            }

            if(A[p.second] == INF) {

                for (auto temp_ : adj[p.second]) {
                    int w = temp_.first;
                    temp[w] = query(i, w);

                    br = ldPair(ldAdd(ldDist(p.first), (landmarks.count(w) > 0)), (p.first & 1)==0?true:false); beta = ldPair(ldPair(temp[w], (distances_1[w][i] == INF)), true);
                    if(br <= beta)
                        que1.push(std::make_pair(br, w));
                }

                A[p.second] = p.first;
                AFF_VERTS[c] = p.second; c++;
            }
        }

        // computing boundary vertices
        std::vector<std::pair<int, int> > V_aff;
        for(int j = 0; j < c; j++) {
            temp[AFF_VERTS[j]] = INF;
            for (auto temp_ : adj[AFF_VERTS[j]]) {
                int w = temp_.first;
                if(A[w] == INF)
                    temp[AFF_VERTS[j]] = min(temp[AFF_VERTS[j]], temp[w] + 1);
            }

            if(temp[AFF_VERTS[j]] != INF)
                V_aff.push_back(std::make_pair(temp[AFF_VERTS[j]], AFF_VERTS[j]));
            landmarks.count(AFF_VERTS[j])>0?highway[i][landmarks[AFF_VERTS[j]]] = INF:distances[AFF_VERTS[j]][i] = INF;
        }

        // updating the labelling
        if(V_aff.size() > 0) {
            std::sort(V_aff.begin(), V_aff.end());

            std::queue<int> quee[2]; int use = 0;
            int d = V_aff[0].first; int x = 0;
            while(x < V_aff.size() && V_aff[x].first == d) {
                if(prunable(i, V_aff[x].second, temp, A)){
                    quee[1].push(V_aff[x].second);
                } else {
                    distances[V_aff[x].second][i] = temp[V_aff[x].second];
                    quee[0].push(V_aff[x].second);
                }
                A[V_aff[x].second] = INF;
                x++;
            }

            quee[use].push(-1);
            while (!quee[0].empty() || x < V_aff.size()) {
                int u = quee[use].front();
                quee[use].pop();

                if(u == -1) {
                    if(use == 0) { d++;
                        while(x < V_aff.size() && V_aff[x].first == d) {
                            if(A[V_aff[x].second] != INF) {
                                if(prunable(i, V_aff[x].second, temp, A)) {
                                    quee[1].push(V_aff[x].second);
                                } else {
                                    distances[V_aff[x].second][i] = temp[V_aff[x].second];
                                    quee[0].push(V_aff[x].second);
                                }
                                A[V_aff[x].second] = INF;
                            }
                            x++;
                        }
                    }
                    use = 1 - use;
                    quee[use].push(-1);
                    continue;
                }

                for (auto temp_ : adj[u]) {
                    int w = temp_.first;
                    if(A[w] != INF) {
                        temp[w] = temp[u] + 1;
                        if(use == 1 || prunable(i, w, temp, A)) {
                            quee[1].push(w);
                        } else {
                            distances[w][i] = temp[w];
                            quee[0].push(w);
                        }
                        A[w] = INF;
                    }
                }
            }
        }

        if(i == K - 1)
            continue;

        for(int j = 0; j < c; j++) {
            if(A[AFF_VERTS[j]] != INF)
                A[AFF_VERTS[j]] = INF;
        }
    }

    delete [] A;
    delete [] temp;
    delete [] AFF_VERTS;
}
//original correct version
void HighwayLabellingWW::BHL(std::vector<std::pair<std::string, std::pair<int, int> > > & updates) {
    vector<int> A(V,INF);//used to judge whether the vertex belongs to AFF_VERTS
    vector<int> temp(V,INF);// distance from all vertex to current landmark
    int *AFF_VERTS = new int[V]; //affected vertex
    int c;//index of affected vertex
    int a, b;

    for(int i = 0; i < K; i++) {//for each landmark

        /// Algorithm 2: computing affected vertices
        std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int> >, std::greater<std::pair<int, int> > > que;//priority queue
        std::queue<std::pair<int, int> > que1; //queue
        c = 0;
        for(std::pair<std::string, std::pair<int, int> > iter : updates) {//for each update
            a = iter.second.first; b = iter.second.second;
            temp[a] = query(i, a);
            temp[b] = query(i, b);

            if(temp[a] > temp[b]) {
                que.push(std::make_pair(temp[b] + 1, a));
            } else if(temp[a] < temp[b]) {
                que.push(std::make_pair(temp[a] + 1, b));
            }
        }

        while (!que1.empty() || !que.empty()) {

            std::pair<int, int> p;//<distance, id>
            if(!que.empty() && !que1.empty()) {
                if(que.top().first < que1.front().first) {
                    p = que.top(); que.pop();
                } else {
                    p = que1.front(); que1.pop();
                }
            } else if(!que1.empty()) {
                p = que1.front(); que1.pop();
            } else {//if que1 is empty
                p = que.top(); que.pop();
            }

            if(A[p.second] == INF) {

                for (auto temp_ : adj[p.second]) {
                    int w = temp_.first;
                    temp[w] = query(i, w);
                    if(p.first + 1 <= temp[w]) {
                        que1.push(std::make_pair(p.first + 1, w));
                    }
                }

                A[p.second] = p.first;
                AFF_VERTS[c] = p.second; c++;
            }
        }
//        set<int> temp_set; temp_set.clear();
//        for(int j=0;j<c;j++){
//            temp_set.insert(AFF_VERTS[j]);
//        }
//        cout<<"Size: "<<temp_set.size()<<endl;
        /// computing boundary vertices
        std::vector<std::pair<int, int> > V_aff;//
        for(int j = 0; j < c; j++) {//for each vertex in AFF_VERTS
            temp[AFF_VERTS[j]] = INF;
            for (auto temp_ : adj[AFF_VERTS[j]]) {
                int w = temp_.first;
                if(A[w] == INF)//if w is not affected
                    temp[AFF_VERTS[j]] = min(temp[AFF_VERTS[j]], temp[w]+1);
            }

            if(temp[AFF_VERTS[j]] != INF)
                V_aff.emplace_back(temp[AFF_VERTS[j]], AFF_VERTS[j]);
            distances[AFF_VERTS[j]][i] = INF;
//            int temp_dis = Dijkstra(AFF_VERTS[j],landmarkTopK[i],adj);
//            if(temp_dis != distances[AFF_VERTS[j]][i])
//                cout<<"Wrong!! "<<AFF_VERTS[j]<<" "<<i<<": "<<distances[AFF_VERTS[j]][i]<<" "<<temp_dis<<endl;
        }

        /// Algorithm 4: updating the labelling
        if(V_aff.size() > 0) {
            std::sort(V_aff.begin(), V_aff.end());

            std::queue<int> quee[2]; int use = 0;
            int d = V_aff[0].first; int x = 0;//get vertex with minimal distance bounds
            while(x < V_aff.size() && V_aff[x].first == d) {//to get the set of Vmin
                if(prunable(i, V_aff[x].second, temp, A)) {//if it is prunable
                    quee[1].push(V_aff[x].second);
                } else {
                    distances[V_aff[x].second][i] = temp[V_aff[x].second];
                    int temp_dis = Dijkstra(V_aff[x].second,landmarkTopK[i],adj);
                    if(temp_dis != distances[V_aff[x].second][i])
                        cout<<"Wrong! "<<V_aff[x].second<<" "<<i<<": "<<distances[V_aff[x].second][i]<<" "<<temp_dis<<endl;
                    quee[0].push(V_aff[x].second);
                }
                A[V_aff[x].second] = INF;
                x++;
            }

            quee[use].push(-1);
            while (!quee[0].empty() || x < V_aff.size()) {
                int u = quee[use].front();
                quee[use].pop();

                if(u == -1) {
                    if(use == 0) { d++;
                        while(x < V_aff.size() && V_aff[x].first == d) {
                            if(A[V_aff[x].second] != INF) {
                                if(prunable(i, V_aff[x].second, temp, A)) {
                                    quee[1].push(V_aff[x].second);
                                } else {
                                    distances[V_aff[x].second][i] = temp[V_aff[x].second];
                                    int temp_dis = Dijkstra(V_aff[x].second,landmarkTopK[i],adj);
                                    if(temp_dis != distances[V_aff[x].second][i])
                                        cout<<"Wrong! "<<V_aff[x].second<<" "<<i<<": "<<distances[V_aff[x].second][i]<<" "<<temp_dis<<endl;
                                    quee[0].push(V_aff[x].second);
                                }
                                A[V_aff[x].second] = INF;
                            }
                            x++;
                        }
                    }
                    use = 1 - use;
                    quee[use].push(-1);
                    continue;
                }

                for (auto ww : adj[u]) {
                    int w = ww.first;
                    if(A[w] != INF) {
                        temp[w] = temp[u] + 1;
                        if(use == 1 || prunable(i, w, temp, A)) {
                            quee[1].push(w);
                        } else {
                            distances[w][i] = temp[w];
                              int temp_dis = Dijkstra(w,landmarkTopK[i],adj);
                              if(temp_dis != distances[w][i]){
                                  cout<<"Wrong! "<<w<<" "<<i<<": "<<distances[w][i]<<" "<<temp_dis<<endl;
                              }
                            quee[0].push(w);
                        }
                        A[w] = INF;
                    }
                }
            }
        }

        if(i == K - 1)
            continue;

        for(int j = 0; j < c; j++) {
            if(A[AFF_VERTS[j]] != INF)
                A[AFF_VERTS[j]] = INF;
        }
    }
    delete [] AFF_VERTS;
}
void HighwayLabellingWW::BatchSearch(vector<pair<string, pair<int, int> > > & updates, int rid , vector<int> & d_G, vector<int> & A, set<int> & V_AFF){//
    /// Algorithm 2: Batch Search, computing affected vertices
//    benchmark::heap<2, int, int> Q(V);
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int> >, std::greater<std::pair<int, int> > > Q;//priority queue
    int v,d;
    int a,b;
    int w,weight;
//    set<int> V_AFF;
    V_AFF.clear();

    for(std::pair<std::string, std::pair<int, int> > iter : updates) {//for each update
        a = iter.second.first; b = iter.second.second;
        d_G[a] = query(rid, a);//d_G(r,a)
        d_G[b] = query(rid, b);
        if(d_G[a] > d_G[b]) {
            Q.push(std::make_pair(d_G[b] + 1, a));
//            Q.update(a,d_G[b] + 1);//a is anchor while b is pre-anchor adj[b][a]
        } else if(d_G[a] < d_G[b]) {
            Q.push(std::make_pair(d_G[a] + 1, b));
//            Q.update(b,d_G[a] + 1);//b is anchor while a is pre-anchor adj[a][b]
        }
    }

    while (!Q.empty()) {//use two queues to get the minimal item
        v = Q.top().second; d = Q.top().first;
        Q.pop();

        if(A[v] == INF){//if not found, i.e., v does not belong to V_AFF+
            A[v] = d;
            V_AFF.insert(v);
            for (auto temp_ : adj[v]) {//for each neighbor of N_G'(v)
                w = temp_.first; weight = temp_.second;
                assert(weight == 1);
                d_G[w] = query(rid, w);//from landmark i to w
                if(d + weight <= d_G[w]) {//identify the composite-path-affected vertices
                    Q.push(std::make_pair(d + weight, w));
                    d_G[w] = d + weight;
                }
            }
        }

    }
//    return V_AFF;
}
//modified version
void HighwayLabellingWW::BHLW(std::vector<std::pair<std::string, std::pair<int, int> > > & updates) {
    vector<int> A;//used to judge whether the vertex belongs to AFF_VERTS
    vector<int> temp;// distance from all vertex to current landmark
    int *AFF_VERTS = new int[V]; //affected vertex
    int c=0;//index of affected vertex
    int a, b;
    int v;
    int w, weight;
    vector<pair<int,bool>> D_BOU(V,make_pair(INF,false));//(d,l)
    benchmark::heap<2, int, int> V_AFF(V);//affected vertices
    set<int> V_AFFPlus;

    for(int i = 0; i < K; i++) {//for each landmark
        A.assign(V,INF); temp.assign(V,INF); D_BOU.assign(V,make_pair(INF,false));
        /// Algorithm 2: Batch Search, computing affected vertices
        BatchSearch(updates,i,temp,A,V_AFFPlus);

        /// Algorithm 4: Batch Repair
        // Step 1: computing landmark distance bound
        for(auto it=V_AFFPlus.begin();it!=V_AFFPlus.end();++it){
            v = *it;
            temp[v] = INF;
            for (auto temp_ : adj[v]) {
                w = temp_.first; weight = temp_.second;
                if(A[w] == INF){//if w is not affected
                    temp[v] = min(temp[v], temp[w]+weight);
                    D_BOU[v].first = temp[v];
                    if(landmarkSet.find(v) != landmarkSet.end()){//if v is landmark
                        D_BOU[v].second = true;
                    }
                }
            }
            V_AFF.update(v,temp[v]);
//            distances[v][i] = INF;
        }

        /*for(int j = 0; j < c; j++) {//for each vertex in AFF_VERTS
            int v = AFF_VERTS[j];
            temp[v] = INF;
            for (auto temp_ : adj[v]) {
                int w = temp_.first; int weight = temp_.second;
                if(A[w] == INF){//if w is not affected
                    temp[v] = min(temp[v], temp[w]+weight);
                    D_BOU[v].first = temp[v];
                    if(landmarkSet.find(v) != landmarkSet.end()){//if v is landmark
                        D_BOU[v].second = true;
                    }
                }
            }

            if(temp[v] != INF){
//                V_Aff.insert(temp[AFF_VERTS[j]], AFF_VERTS[j]);
                V_aff.emplace_back(temp[v], v);
                V_AFF.update(v,D_BOU[v].first);
            }
//            distances[AFF_VERTS[j]][i] = INF;
        }*/

        // Step 2: updating the labels
        int topId,topDis; set<int> V_min;
        /*while(!V_Aff.empty()){
            // get V_min
            topDis = V_Aff.begin()->first;
            auto Vmin = V_Aff.equal_range(topDis);
            topIDs.clear();
            for(auto it=Vmin.first; it!=Vmin.second; ++it){
                topIDs.push_back(it->second);
            }
            for(int id = 0;id<topIDs.size();++id){//for each element in Vmin
                int v=topIDs[id];
                if(temp[v] == INF || landmarks.find(v)!=landmarks.end()){
                    distances[v][i] = INF;
                }else{
                    distances[v][i] = temp[v];
                }
                if(landmarks.find(v)!=landmarks.end()){
                    highway[i][landmarks[v]] = temp[v];
                    highway[landmarks[v]][i] = temp[v];
                    int temp_dis = Dijkstra(landmarkTopK[i],v,adj);
                    if(temp_dis != highway[i][landmarks[v]])
                        cout<<"Prunable Wrong! "<<i<<" "<<landmarks[v]<<": "<<highway[i][landmarks[v]]<<" "<<temp_dis<<" "<<Dijkstra(landmarkTopK[i],v,adjOldO)<<endl;
                }

                for (auto ww : adj[v]) {
                    int w = ww.first; int weight = ww.second;
                    if(A[w] != INF) {
                        temp[w] = temp[v] + weight;
                        if(use == 1 || prunableW(i, w, temp, A)) {
                            quee[1].push(w);
                        } else {
                            distances[w][i] = temp[w];
//                            int temp_dis = Dijkstra(w,landmarkTopK[i],adj);
//                            if(temp_dis != distances[w][i]){
//                                cout<<"Wrong! "<<w<<" "<<i<<": "<<distances[w][i]<<" "<<temp_dis<<endl;
//                            }
                        }
                        A[w] = INF;
                    }
                }


                if(prunableW(i, topIDs[id], temp, A)){
                    que2.update(V_aff[x].second,V_aff[x].first);
                } else {
                    distances[V_aff[x].second][i] = temp[V_aff[x].second];
//                    int temp_dis = Dijkstra(V_aff[x].second,landmarkTopK[i],adj);
//                    if(temp_dis != distances[V_aff[x].second][i])
//                        cout<<"Wrong! "<<V_aff[x].second<<" "<<i<<": "<<distances[V_aff[x].second][i]<<" "<<temp_dis<<endl;
//                    quee[0].push(V_aff[x].second);
                    que1.update(V_aff[x].second,V_aff[x].first);
                }
                A[V_aff[x].second] = INF;
            }
        }*/
        while(!V_AFF.empty()){
            //get V_min
            V_min.clear();
            V_AFF.extract_min(topId,topDis);
            V_min.insert(topId);
            while(!V_AFF.empty()){
                int next_dis = V_AFF.top_key();
                if(next_dis == topDis){
                    V_AFF.extract_min(topId,topDis);
                    V_min.insert(topId);
                }else{
                    break;
                }
            }
            //deal with every element in V_min
            for(auto it=V_min.begin();it!=V_min.end();++it){
                v = *it;
                A[v] = INF;
                bool ifLandmark = false;
                //set distance label
                if(D_BOU[v].first == INF || D_BOU[v].second){
                    distances[v][i] = INF;
                }else{
                    distances[v][i] = D_BOU[v].first;
//                    int temp_di = Dijkstra(v,landmarkTopK[i],adj);
//                    if(distances[v][i] != temp_di){
//                        cout<<"Distance incorrect! "<<landmarkTopK[i]<<" "<<v<<" "<<distances[v][i]<<" "<< temp_di<<endl;
//                    }
                }
                //set highway label
                if(landmarkSet.find(v) != landmarkSet.end()){//if it is a landmark
                    highway[landmarks[v]][i] = D_BOU[v].first;
                    highway[i][landmarks[v]] = D_BOU[v].first;
//                    int temp_di = Dijkstra(v,landmarkTopK[i],adj);
//                    if(highway[landmarks[v]][i] != temp_di){
//                        cout<<"Highway incorrect! "<<landmarkTopK[i]<<" "<<v<<" "<<highway[landmarks[v]][i]<<" "<< temp_di<<endl;
//                    }
                    ifLandmark = true;
                }
                //expansion
                for (auto ww : adj[v]) {
                    w = ww.first; weight = ww.second;
                    if(A[w] != INF && V_min.find(w) == V_min.end()) {//if w is in V_AFF and w is not in V_min
                        D_BOU[w].first = D_BOU[v].first + weight;
                        D_BOU[w].second = D_BOU[v].second;
                        if(ifLandmark){
                            D_BOU[w].second = true;
                        }
                        V_AFF.update(w,D_BOU[w].first);
//                        distances[w][i] = temp[w];
                        A[w] = INF;
                    }
                }
            }

        }

    }
    delete [] AFF_VERTS;
}
// algorithm 1
void HighwayLabellingWW::BatchSearch(vector<pair<pair<int,int>,int>> & updates, int rid , vector<int> & d_G, vector<int> & A, set<int> & V_AFF){//
    /// Algorithm 2: Batch Search, computing affected vertices
//    benchmark::heap<2, int, int> Q(V);
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int> >, std::greater<std::pair<int, int> > > Q;//priority queue
    int v,d;
    int a,b;
    int w,weight;
//    set<int> V_AFF;
    V_AFF.clear();

    for(auto iter : updates) {//for each update
        a = iter.first.first; b = iter.first.second;
        d_G[a] = query(rid, a);//d_G(r,a)
        d_G[b] = query(rid, b);
        if(d_G[a] > d_G[b]) {
            Q.push(std::make_pair(d_G[b] + adj[b][a], a));
//            Q.update(a,d_G[b] + 1);//a is anchor while b is pre-anchor adj[b][a]
        } else if(d_G[a] < d_G[b]) {
            Q.push(std::make_pair(d_G[a] + adj[a][b], b));
//            Q.update(b,d_G[a] + 1);//b is anchor while a is pre-anchor adj[a][b]
        }
    }

    while (!Q.empty()) {//use two queues to get the minimal item
        v = Q.top().second; d = Q.top().first;
        Q.pop();

        if(A[v] == INF){//if not found, i.e., v does not belong to V_AFF+
            A[v] = d;
            V_AFF.insert(v);
            for (auto temp_ : adj[v]) {//for each neighbor of N_G'(v)
                w = temp_.first; weight = temp_.second;
//                assert(weight == 1);
                d_G[w] = query(rid, w);//from landmark i to w
                if(d + weight <= d_G[w]) {//identify the composite-path-affected vertices
                    Q.push(std::make_pair(d + weight, w));
                    d_G[w] = d + weight;
                }
            }
        }

    }
//    return V_AFF;
}
//function of BHL for weighted graph
void HighwayLabellingWW::BHLW(vector<pair<pair<int,int>,int>> & updates) {
    vector<int> A;//used to judge whether the vertex belongs to AFF_VERTS
    vector<int> temp;// distance from all vertex to current landmark
    int *AFF_VERTS = new int[V]; //affected vertex
    int c=0;//index of affected vertex
    int a, b;
    int v;
    int w, weight;
    vector<pair<int,bool>> D_BOU(V,make_pair(INF,false));//(d,l)
    benchmark::heap<2, int, int> V_AFF(V);//affected vertices
    set<int> V_AFFPlus;

    for(int i = 0; i < K; i++) {//for each landmark
        A.assign(V,INF); temp.assign(V,INF); D_BOU.assign(V,make_pair(INF,false)); V_AFFPlus.clear();
        /// Algorithm 2: Batch Search, computing affected vertices
        BatchSearch(updates,i,temp,A,V_AFFPlus);
//        int cc = 0;
//        for(auto it=A.begin();it!=A.end();++it){
//            if(*it != INF){
//                cc++;
//            }
//        }
//        cout<<"cc: "<<cc<<endl;
        /// Algorithm 4: Batch Repair
//        std::vector<std::pair<int, int> > V_aff;
        // Step 1: computing landmark distance bound
        for(auto it=V_AFFPlus.begin();it!=V_AFFPlus.end();++it){
            v = *it;
//            if(v == 1602)
//                cout<<"affected "<<v<<endl;
            temp[v] = INF;
            for (auto temp_ : adj[v]) {
                w = temp_.first; weight = temp_.second;
                if(A[w] == INF){//if w is not affected
                    temp[v] = min(temp[v], temp[w]+weight);
                    D_BOU[v].first = temp[v];
                    if(landmarkSet.find(v) != landmarkSet.end()){//if v is landmark
                        D_BOU[v].second = true;
                    }
                }
            }
            V_AFF.update(v,temp[v]);
//            distances[v][i] = INF;
        }

        // Step 2: updating the labels
        int topId,topDis; set<int> V_min;
        /*while(!V_Aff.empty()){
            // get V_min
            topDis = V_Aff.begin()->first;
            auto Vmin = V_Aff.equal_range(topDis);
            topIDs.clear();
            for(auto it=Vmin.first; it!=Vmin.second; ++it){
                topIDs.push_back(it->second);
            }
            for(int id = 0;id<topIDs.size();++id){//for each element in Vmin
                int v=topIDs[id];
                if(temp[v] == INF || landmarks.find(v)!=landmarks.end()){
                    distances[v][i] = INF;
                }else{
                    distances[v][i] = temp[v];
                }
                if(landmarks.find(v)!=landmarks.end()){
                    highway[i][landmarks[v]] = temp[v];
                    highway[landmarks[v]][i] = temp[v];
                    int temp_dis = Dijkstra(landmarkTopK[i],v,adj);
                    if(temp_dis != highway[i][landmarks[v]])
                        cout<<"Prunable Wrong! "<<i<<" "<<landmarks[v]<<": "<<highway[i][landmarks[v]]<<" "<<temp_dis<<" "<<Dijkstra(landmarkTopK[i],v,adjOldO)<<endl;
                }

                for (auto ww : adj[v]) {
                    int w = ww.first; int weight = ww.second;
                    if(A[w] != INF) {
                        temp[w] = temp[v] + weight;
                        if(use == 1 || prunableW(i, w, temp, A)) {
                            quee[1].push(w);
                        } else {
                            distances[w][i] = temp[w];
//                            int temp_dis = Dijkstra(w,landmarkTopK[i],adj);
//                            if(temp_dis != distances[w][i]){
//                                cout<<"Wrong! "<<w<<" "<<i<<": "<<distances[w][i]<<" "<<temp_dis<<endl;
//                            }
                        }
                        A[w] = INF;
                    }
                }


                if(prunableW(i, topIDs[id], temp, A)){
                    que2.update(V_aff[x].second,V_aff[x].first);
                } else {
                    distances[V_aff[x].second][i] = temp[V_aff[x].second];
//                    int temp_dis = Dijkstra(V_aff[x].second,landmarkTopK[i],adj);
//                    if(temp_dis != distances[V_aff[x].second][i])
//                        cout<<"Wrong! "<<V_aff[x].second<<" "<<i<<": "<<distances[V_aff[x].second][i]<<" "<<temp_dis<<endl;
//                    quee[0].push(V_aff[x].second);
                    que1.update(V_aff[x].second,V_aff[x].first);
                }
                A[V_aff[x].second] = INF;
            }
        }*/
        while(!V_AFF.empty()){
            //get V_min
            V_min.clear();
            V_AFF.extract_min(topId,topDis);
            V_min.insert(topId);
            while(!V_AFF.empty()){
                int next_dis = V_AFF.top_key();
                if(next_dis == topDis){
                    V_AFF.extract_min(topId,topDis);
                    V_min.insert(topId);
                }else{
                    break;
                }
            }
            //deal with every element in V_min
            for(auto it=V_min.begin();it!=V_min.end();++it){
                v = *it;
                A[v] = INF;
                bool ifLandmark = false;
                //set distance label
                if(D_BOU[v].first == INF || D_BOU[v].second){//if v pass another landmark or distance bound of v is INF
                    if(distances[v][i]<INF){
                        distances[v][i] = INF;
                    }
//                    distances[v][i] = D_BOU[v].first;
                }else{
//                    cout<<distances[v][i]<<" "<<D_BOU[v].first<<endl;
                    distances[v][i] = D_BOU[v].first;
//                    if(ifDebug){
//                        int temp_di = Dijkstra(v,landmarkTopK[i],adj);
//                        if(distances[v][i] != temp_di){
//                            cout<<"Distance incorrect! "<<landmarkTopK[i]<<" "<<v<<" "<<distances[v][i]<<" "<< temp_di<<" "<<Dijkstra(v,landmarkTopK[i],adjOldO)<<endl;
//                        }
//                    }

                }
                //set highway label
                if(landmarkSet.find(v) != landmarkSet.end()){//if it is a landmark
                    highway[landmarks[v]][i] = D_BOU[v].first;
                    highway[i][landmarks[v]] = D_BOU[v].first;
//                    if(ifDebug){
//                        int temp_di = Dijkstra(v,landmarkTopK[i],adj);
//                        if(highway[landmarks[v]][i] != temp_di){
//                            cout<<"Highway incorrect! "<<landmarkTopK[i]<<" "<<v<<" "<<highway[landmarks[v]][i]<<" "<< temp_di<<endl;
//                        }
//                    }

                    ifLandmark = true;
                }
                //expansion
                for (auto ww : adj[v]) {
                    w = ww.first; weight = ww.second;
//                    if(w == 1602)
//                        cout<<w<<endl;
                    if(A[w] != INF && V_min.find(w) == V_min.end()) {//if w is in V_AFF
                        if(D_BOU[w].first >= D_BOU[v].first + weight){
                            D_BOU[w].first = D_BOU[v].first + weight;
                            D_BOU[w].second = D_BOU[v].second;
                            V_AFF.update(w,D_BOU[w].first);
//                        distances[w][i] = D_BOU[w].first;
//                        A[w] = INF;
                            if(ifLandmark){
                                D_BOU[w].second = true;
                            }
                        }
                    }
                }
            }

        }

    }
    delete [] AFF_VERTS;
}

//function of BHL for weighted graph
void HighwayLabellingWW::BHLWP(vector<pair<pair<int,int>,int>> & updates, bool ifParallel) {

    if(ifParallel){
        //multiple thread
        if(K>threadnum){
            int step=K/threadnum;
            boost::thread_group thread;
            for(int i=0;i<threadnum;i++){
                pair<int,int> p;
                p.first=i*step;
                if(i==threadnum-1)
                    p.second=K;
                else
                    p.second=(i+1)*step;
                thread.add_thread(new boost::thread(&HighwayLabellingWW::UpdateLandmark, this, p, boost::ref(updates)));
            }
            thread.join_all();
        }else{
            boost::thread_group thread;
            for(int i=0;i<K;i++){
                thread.add_thread(new boost::thread(&HighwayLabellingWW::UpdateLandmark, this, make_pair(i,i+1), boost::ref(updates)));
            }
            thread.join_all();
        }
    }
    else{
        for(int i = 0; i < K; i++) {//for each landmark
            UpdateLandmark(make_pair(i,i+1),updates);
        }
    }


}

void HighwayLabellingWW::UpdateLandmark(pair<int,int> landmarkRange, vector<pair<pair<int,int>,int>> & updates){
    vector<int> A;//used to judge whether the vertex belongs to AFF_VERTS
    vector<int> temp;// distance from all vertex to current landmark
    int *AFF_VERTS = new int[V]; //affected vertex
    vector<pair<int,bool>> D_BOU;//(d,l)
    set<int> V_AFFPlus;
    int v;
    int w, weight;
    benchmark::heap<2, int, int> V_AFF(V);//affected vertices

    for(int landmark_i = landmarkRange.first; landmark_i < landmarkRange.second; landmark_i++) {//for each landmark
        A.assign(V,INF); temp.assign(V,INF); D_BOU.assign(V,make_pair(INF,false)); V_AFFPlus.clear();
        /// Algorithm 2: Batch Search, computing affected vertices
        BatchSearch(updates,landmark_i,temp,A,V_AFFPlus);

        /// Algorithm 4: Batch Repair
//        std::vector<std::pair<int, int> > V_aff;
        // Step 1: computing landmark distance bound
        for(auto it=V_AFFPlus.begin();it!=V_AFFPlus.end();++it){
            v = *it;
//            if(v == 1602)
//                cout<<"affected "<<v<<endl;
            temp[v] = INF;
            for (auto temp_ : adj[v]) {
                w = temp_.first; weight = temp_.second;
                if(A[w] == INF){//if w is not affected
                    temp[v] = min(temp[v], temp[w]+weight);
                    D_BOU[v].first = temp[v];
                    if(landmarkSet.find(v) != landmarkSet.end()){//if v is landmark
                        D_BOU[v].second = true;
                    }
                }
            }
            V_AFF.update(v,temp[v]);
//            distances[v][i] = INF;
        }

        // Step 2: updating the labels
        int topId,topDis; set<int> V_min;

        while(!V_AFF.empty()){
            //get V_min
            V_min.clear();
            V_AFF.extract_min(topId,topDis);
            V_min.insert(topId);
            while(!V_AFF.empty()){
                int next_dis = V_AFF.top_key();
                if(next_dis == topDis){
                    V_AFF.extract_min(topId,topDis);
                    V_min.insert(topId);
                }else{
                    break;
                }
            }
            //deal with every element in V_min
            for(auto it=V_min.begin();it!=V_min.end();++it){
                v = *it;
                A[v] = INF;
                bool ifLandmark = false;
                //set distance label
                if(D_BOU[v].first == INF || D_BOU[v].second){//if v pass another landmark or distance bound of v is INF
                    if(distances[v][landmark_i]<INF){
                        distances[v][landmark_i] = INF;
                    }
//                    distances[v][i] = D_BOU[v].first;
                }else{
//                    cout<<distances[v][i]<<" "<<D_BOU[v].first<<endl;
                    distances[v][landmark_i] = D_BOU[v].first;
//                    if(ifDebug){
//                        int temp_di = Dijkstra(v,landmarkTopK[i],adj);
//                        if(distances[v][i] != temp_di){
//                            cout<<"Distance incorrect! "<<landmarkTopK[i]<<" "<<v<<" "<<distances[v][i]<<" "<< temp_di<<" "<<Dijkstra(v,landmarkTopK[i],adjOldO)<<endl;
//                        }
//                    }

                }
                //set highway label
                if(landmarkSet.find(v) != landmarkSet.end()){//if it is a landmark
                    highway[landmarks[v]][landmark_i] = D_BOU[v].first;
                    highway[landmark_i][landmarks[v]] = D_BOU[v].first;
//                    if(ifDebug){
//                        int temp_di = Dijkstra(v,landmarkTopK[i],adj);
//                        if(highway[landmarks[v]][i] != temp_di){
//                            cout<<"Highway incorrect! "<<landmarkTopK[i]<<" "<<v<<" "<<highway[landmarks[v]][i]<<" "<< temp_di<<endl;
//                        }
//                    }

                    ifLandmark = true;
                }
                //expansion
                for (auto ww : adj[v]) {
                    w = ww.first; weight = ww.second;
//                    if(w == 1602)
//                        cout<<w<<endl;
                    if(A[w] != INF && V_min.find(w) == V_min.end()) {//if w is in V_AFF
                        if(D_BOU[w].first >= D_BOU[v].first + weight){
                            D_BOU[w].first = D_BOU[v].first + weight;
                            D_BOU[w].second = D_BOU[v].second;
                            V_AFF.update(w,D_BOU[w].first);
//                        distances[w][i] = D_BOU[w].first;
//                        A[w] = INF;
                            if(ifLandmark){
                                D_BOU[w].second = true;
                            }
                        }
                    }
                }
            }

        }
    }

    delete [] AFF_VERTS;
}

bool HighwayLabellingWW::prunable(int i, int u, int *temp, int *A) {

  if(landmarks.count(u) > 0) {//if this is landmark
    highway[i][landmarks[u]] = temp[u];
    highway[landmarks[u]][i] = temp[u];
      int temp_dis = Dijkstra(landmarkTopK[i],u,adjOld);
      if(temp_dis != highway[i][landmarks[u]])
          cout<<"Wrong! "<<i<<" "<<landmarks[u]<<": "<<highway[i][landmarks[u]]<<" "<<temp_dis<<endl;
    return true;
  } else {//if this is not landmark
    for (auto ww : adj[u]) {
        int w = ww.first;
      if(A[w] == INF) {//if the neighbor has not been explored
        if(temp[w] == temp[u] - 1 && distances[w][i] == INF)
          return true;
      }
    }
  }
  return false;
}
bool HighwayLabellingWW::prunable(int i, int u, vector<int> & temp, vector<int> & A) {

    if(landmarks.count(u) > 0) {//if this is landmark
        highway[i][landmarks[u]] = temp[u];
        highway[landmarks[u]][i] = temp[u];
        int temp_dis = Dijkstra(landmarkTopK[i],u,adjOld);
        if(temp_dis != highway[i][landmarks[u]])
            cout<<"Wrong! "<<i<<" "<<landmarks[u]<<": "<<highway[i][landmarks[u]]<<" "<<temp_dis<<endl;
        return true;
    } else {//if this is not landmark
        for (auto ww : adj[u]) {
            int w = ww.first;
            if(A[w] == INF) {//if the neighbor has not been explored
                if(temp[w] == temp[u] - 1 && distances[w][i] == INF)
                    return true;
            }
        }
    }
    return false;
}
bool HighwayLabellingWW::prunableW(int i, int u, vector<int> & temp, vector<int> & A) {//i: landmark, u: affected vertex

    if(landmarks.count(u) > 0) {//if this is landmark
        highway[i][landmarks[u]] = temp[u];
        highway[landmarks[u]][i] = temp[u];
        int temp_dis = Dijkstra(landmarkTopK[i],u,adj);
        if(temp_dis != highway[i][landmarks[u]])
            cout<<"Prunable Wrong! "<<i<<" "<<landmarks[u]<<": "<<highway[i][landmarks[u]]<<" "<<temp_dis<<" "<<Dijkstra(landmarkTopK[i],u,adjOldO)<<endl;
        return true;
    } else {//if this is not landmark
        for (auto ww : adj[u]) {
            int w = ww.first;
            if(A[w] == INF) {//if the neighbor has not been explored
                if(temp[w] == temp[u] - ww.second && distances[w][i] == INF)//
                    return true;
            }
        }
    }
    return false;
}
//function of querying the distance from landmark r to vertex v
int HighwayLabellingWW::query(int r, int v) {

  int m = INF;
  for(int i = 0; i < K; i++) {
    m = min(m, distances_1[v][i] + highway_1[r][i]);
  }
  return m;
}

int HighwayLabellingWW::min(int a, int b) {
  return (a < b) ? a : b;
}

vector<int> HighwayLabellingWW::OrderRead(string filename){
    vector<int> vNodeOrder;
    vNodeOrder.assign(V,-1);

    ifstream IF(filename);
    if(!IF){
        cout<<"Cannot open file "<<filename<<endl;
        exit(1);
    }

    int num,ID,order;
    IF>>num;
    for(int i=0;i<num;i++){
        IF>>ID>>order;
        vNodeOrder[order]=ID;
    }
    return vNodeOrder;
}

//Function of selecting landmarks (high-degree vertices)
void HighwayLabellingWW::SelectLandmarks_HD(vector<int> & topk,string graphFile) {
    if(dataset == "NY" || dataset == "W" || dataset == "FLA" || dataset =="USA"){
        vector<int> vNodeOrder = OrderRead(graphFile+".order");
        for (int v = 0; v < K; v++){
            topk[v] = vNodeOrder[V - v - 1];
            //cout<<vNodeOrder[V - v - 1] <<" "<<adj[vNodeOrder[V - v - 1]].size()<<endl;
            landmarkSet.insert(topk[v]);
        }
    }else{
        std::vector<std::pair<int, int> > deg(V); long sum = 0;
        for (int v = 0; v < V; v++) {
            deg[v] = std::make_pair(adj[v].size(), v);
            sum = sum + adj[v].size();
        }
        std::sort(deg.rbegin(), deg.rend());

        for (int v = 0; v < K; v++){
            topk[v] = deg[v].second;
            landmarkSet.insert(topk[v]);
//            cout<<topk[v]<<" "<<adj[topk[v]].size()<<endl;
        }
    }

  assert(landmarkSet.size() == K);
}

void HighwayLabellingWW::RemoveLandmarks(vector<int> & topk) {

    if(ifDebug){
        adjOld = adj;
    }
  for(int i = 0; i < K; i++) {
    for (auto vv : adj[topk[i]]) {
        int v = vv.first;
      adj[v].erase(topk[i]);
    }
    adj[topk[i]].clear();
  }
}

void HighwayLabellingWW::HC_UB_naiveP(int s, int t, bool ifFull, int & dis) {
    if(!ifFull) {
        dis = HC_UB_naive(s, t);//distance bound
    }else{
        dis = HC_UB_naive_Full(s,t);
    }
}

int HighwayLabellingWW::HC_UB_naive(int s, int t) {

  int m = INF; int i, j;
  int r1,r2;
  int temp_d;
  for(i = 0; i < C[s]; i++) {
      r1 = vertices[s][i];
//      temp_d = distances[s][i];
//      int temp_dis = Dijkstra(s,landmarkTopK[r1],adjOld);
//      if(temp_dis != temp_d)
//          cout<<"Wrong! "<<s<<" "<<r1<<": "<<temp_d<<" "<<temp_dis<<endl;
      for (j = 0; j < C[t]; j++){
          r2 = vertices[t][j];
//          temp_d = distances[t][j];
//          cout<<"Dis Bound: "<<distances[s][r1]<<" "<<highway[r1][r2]<<" "<<distances[t][r2]<<endl;
//          cout<<"("<<s<<","<<vertices[s][r1]<<") ("<<vertices[t][r2]<<","<<t<<")"<<endl;
//          temp_dis = Dijkstra(t,landmarkTopK[r2],adjOld);
//          if(temp_dis != temp_d)
//              cout<<"Wrong!! "<<t<<" "<<r2<<": "<<temp_d<<" "<<temp_dis<<endl;
//          temp_dis = Dijkstra(landmarkTopK[r1],landmarkTopK[r2],adjOld);
//          temp_d = highway[r1][r2];
//          if(temp_dis != temp_d){
//              cout<<"Wrong!!! "<<r1<<" "<<r2<<": "<<temp_d<<" "<<temp_dis<<endl;
//          }
          m = min(m, distances[s][i] + highway[r1][r2] + distances[t][j]);
      }
  }
  assert(m<=INF);
  return m;
}

int HighwayLabellingWW::HC_UB_naive_Full(int s, int t) {

    int m = INF; int i, j;
    int r1,r2;
    for(i = 0; i < C[s]; i++) {
        r1 = vertices[s][i];
//        int temp_dis = Dijkstra(s,landmarkTopK[r1],adjOld);
//        if(temp_dis != distances[s][r1])
//            cout<<"Wrong! "<<s<<" "<<r1<<": "<<distances[s][r1]<<" "<<temp_dis<<endl;
        for (j = 0; j < C[t]; j++){
            r2 = vertices[t][j];
//            cout<<"Dis Bound: "<<distances[s][r1]<<" "<<highway[r1][r2]<<" "<<distances[t][r2]<<endl;
//            cout<<"("<<s<<","<<vertices[s][r1]<<") ("<<vertices[t][r2]<<","<<t<<")"<<endl;
//            temp_dis = Dijkstra(t,landmarkTopK[r2],adjOld);
//            if(temp_dis != distances[t][r2])
//                cout<<"Wrong! "<<t<<" "<<r2<<": "<<distances[t][r2]<<" "<<temp_dis<<endl;
//            temp_dis = Dijkstra(landmarkTopK[r1],landmarkTopK[r2],adjOld);
//            if(temp_dis != highway[r1][r2]){
//                cout<<"Wrong! "<<r1<<" "<<r2<<": "<<highway[r1][r2]<<" "<<temp_dis<<endl;
//            }
            m = min(m, distances[s][r1] + highway[r1][r2] + distances[t][r2]);
        }
    }
    assert(m<=INF);
    return m;
}

/*void HighwayLabellingWW::QueryDistance(std::string pairs, std::string output, int runtimes) {
  std::vector<TwoLayerQueue> qque; std::vector<int> qdist[2];

  qdist[0].resize(V, INF); qdist[1].resize(V, INF);
  qque.push_back(TwoLayerQueue(V)); qque.push_back(TwoLayerQueue(V));

  time_querying_millisec_ = 0; int s = 0, t = 0; int total = 0;
  std::ifstream ifs(pairs); std::ofstream ofs(output);
  if (!ifs) {
      cout<<"Cannot open file "<<pairs<<endl;
      exit(1);
  }
  if (!ofs) {
      cout<<"Cannot open file "<<output<<endl;
      exit(1);
  }
  vector<pair<int,int>> ODs;
  int times;
  ifs >> times;
  if(runtimes > times){
      runtimes = times;
  }
  cout<<"Run times: "<<runtimes<<endl;
  while(ifs >> s >> t){
      ODs.emplace_back(s,t);
  }
  ifs.close();
  for(int i=0;i<runtimes;++i){
      s = ODs[i].first; t = ODs[i].second;
      double a = -GetCurrentTimeMilliSec();

      int dist_upper = HC_UB_naive(s, t);//distance bound
      int res = dist_upper, dis[2] = {0, 0};
      for (int dir = 0; dir < 2; dir++){
          int v = dir == 0 ? s : t;
          qque[dir].clear();
          qque[dir].push(v);
          qque[dir].next();
          qdist[dir][v] = 0;
      }

      while (!qque[0].empty() && !qque[1].empty()) {
          int use = (qque[0].size() <= qque[1].size()) ? 0 : 1;
          dis[use]++;

          if (dis[0] + dis[1] == dist_upper) {
              res = dis[0] + dis[1];
              goto LOOP_END;
          }

          while (!qque[use].empty()) {

              int v = qque[use].front();
              qque[use].pop();

              for (auto ww : adj[v]) {
                    int w = ww.first;
                  int &src_d = qdist[    use][w];
                  int &dst_d = qdist[1 - use][w];
                  if (src_d != INF) continue;
                  if (dst_d != INF) {
                      res = qdist[use][v] + 1 + dst_d;
                      goto LOOP_END;
                  } else {
                      qque[use].push(w);
                      qdist[use][w] = qdist[use][v] + 1;
                  }
              }
          }
          qque[use].next();
      }
      LOOP_END:

      a += GetCurrentTimeMilliSec();
      time_querying_millisec_ += a;

      for (int dir = 0; dir < 2; dir++) {
          for (int v : qque[dir]) {
              qdist[dir][v] = INF;
          }
          qque[dir].clear();
      }
      ofs << s << " " << t << " " << (int) min(res, dist_upper) << "\n";
  }
  *//*while(ifs >> s >> t) { total++;

    double a = -GetCurrentTimeMilliSec();

    int dist_upper = HC_UB_naive(s, t);//distance bound
    int res = dist_upper, dis[2] = {0, 0};
    for (int dir = 0; dir < 2; dir++){
      int v = dir == 0 ? s : t;
      qque[dir].clear();
      qque[dir].push(v);
      qque[dir].next();
      qdist[dir][v] = 0;
    }

    while (!qque[0].empty() && !qque[1].empty()) {
      int use = (qque[0].size() <= qque[1].size()) ? 0 : 1;
      dis[use]++;

      if (dis[0] + dis[1] == dist_upper) {
        res = dis[0] + dis[1];
        goto LOOP_END;
      }

      while (!qque[use].empty()) {

        int v = qque[use].front();
        qque[use].pop();

        for (int w : adj[v]) {

          int &src_d = qdist[    use][w];
          int &dst_d = qdist[1 - use][w];
          if (src_d != INF) continue;
          if (dst_d != INF) {
            res = qdist[use][v] + 1 + dst_d;
            goto LOOP_END;
          } else {
            qque[use].push(w);
            qdist[use][w] = qdist[use][v] + 1;
          }
    }
      }
      qque[use].next();
    }
    LOOP_END:

    a += GetCurrentTimeMilliSec();
    time_querying_millisec_ += a;

    for (int dir = 0; dir < 2; dir++) {
      for (int v : qque[dir]) {
        qdist[dir][v] = INF;
      }
      qque[dir].clear();
    }

    ofs << s << " " << t << " " << (int) min(res, dist_upper) << "\n";
  }*//*

  ofs.close();

  std::cout << "Average Query Time (ms) : " << (double) time_querying_millisec_ / runtimes << std::endl;

}*/

void HighwayLabellingWW::EfficencyTest(string pairs,int runtimes){
    string filename = pairs + ".query";
    //Open file
    ifstream inFile(filename, ios::in);
    if (!inFile) {
        cout << "Failed to open file " << filename << endl;
        exit(1);
    }
    int num,ID1,ID2;
    vector<pair<int,int>> ODpair;
    inFile >> num;
    if(runtimes > num)
        runtimes = num;
    cout<<"Run times: "<<runtimes<<endl;
    for (int i = 0; i < runtimes; ++i) {
        inFile >> ID1 >> ID2;
        ODpair.emplace_back(make_pair(ID1, ID2));
    }
    inFile.close();
    // test
    double time=0;
    int disBatchHL;
    Timer tt;
    tt.start();
    for(int i=0;i<runtimes;++i){//runtimes
        ID1 = ODpair[i].first; ID2 = ODpair[i].second;
//        s = 3291; t = 1;
        disBatchHL = BatchHL_Query(ID1,ID2,adj, false);
    }
    tt.stop();
    time = tt.GetRuntime();
    std::cout << "Average Query Time: " << time*1000 / runtimes << " ms."<<std::endl;
}

void HighwayLabellingWW::storeLabelling(std::string filename) {
    string file = std::string(filename) + "/" + dataset + "." +std::to_string(K) + std::string(".index");
    FILE *fout = fopen( file.c_str(), "wb" );
    if(fout == NULL){
        cout<<"Failed to open file "<<file<<endl;
        exit(1);
    }

    for(int i = 0; i < V; i++) {
        int C = 0;
        for(int j = 0; j < K; j++) {
            if(distances[i][j] != INF){//if there is distance label from vertex i to landmark j
                C++;
            }
        }
        if(C==0){
            cout<<"Error! C == 0 : "<<i<<endl;
        }
        assert(C>0);
        fwrite( &C, sizeof(int), 1, fout );//count number
        for(int j = 0; j < K; j++) {
            if(distances[i][j] != INF) {
                fwrite( &j, sizeof(int), 1, fout );//mapped landmark id
                fwrite( &distances[i][j], sizeof(int), 1, fout );//distance
//                int temp_dis = Dijkstra(i,landmarkTopK[j],adj);
//                if(temp_dis!=distances[i][j]){
//                    cout<<"Distance wrong! "<<i<<" "<<j<<": "<<distances[i][j]<<" "<<temp_dis<<endl;
//                }
            }
        }
    }
    fclose(fout);
    file = std::string(filename) + "/" + dataset + "." + std::to_string(K) + std::string(".highway");
    fout = fopen( file.c_str(), "wb" );
    if(fout == NULL){
        cout<<"Failed to open file "<<file<<endl;
        exit(1);
    }
    for(int i = 0; i < K; i++) {
        for(int j = 0; j < K; j++) {
            if(highway[i][j] == INF){
                cout<<"!!!INF dis between landmark: "<<i<<" "<<j<<endl;
            }
            fwrite( &highway[i][j], sizeof(int), 1, fout );//distance from landmark i to landmark j
        }
    }
    fclose(fout);
}

void HighwayLabellingWW::loadLabelling_Full(std::string filename, vector<int> & topk) {
    string file = std::string(filename) + "/" + dataset + "." + std::to_string(K) + std::string(".index");
    FILE *fin = fopen( file.c_str(), "rb" );
    if(fin == NULL){
        cout << "Failed to open file " << file << endl;
        exit(1);
    }
  distances.assign(V,vector<int>(K,INF)); distances_1.assign(V,vector<int>(K,INF));

  int Count, idx;
  for(int i = 0; i < V; i++) {
      fread( &Count, sizeof(int), 1, fin );

    for(int j = 0; j < Count; j++) {
        fread( &idx, sizeof(int), 1, fin );//mapped landmark id
        fread( &distances[i][idx], sizeof(int), 1, fin );
        distances_1[i][idx] = distances[i][idx];
    }
  }
    fclose(fin);

  file = std::string(filename) + "/" + dataset + "." + std::to_string(K) + std::string(".highway");
    fin = fopen( file.c_str(), "rb" );
    if(fin == NULL){
        cout << "Failed to open file " << file << endl;
        exit(1);
    }
  highway.assign(K,vector<int>(K,INF)); highway_1.assign(K,vector<int>(K,INF));
  for(int i = 0; i < K; i++) {
    for(int j = 0; j < K; j++) {
        fread( &highway[i][j], sizeof(int), 1, fin );
        highway_1[i][j] = highway[i][j];
    }
  }
    fclose(fin);
}

void HighwayLabellingWW::loadLabelling_Pruned(std::string filename) {
    string file = std::string(filename) + "/" + dataset + "." + std::to_string(K) + std::string(".index");
    FILE *fin = fopen( file.c_str(), "rb" );
    if(fin == NULL){
        cout << "Failed to open file " << file << endl;
        exit(1);
    }
  C.resize(V);
  vertices.resize(V);
  distances.resize(V);

  for(int i = 0; i < V; i++) {
      fread( &C[i], sizeof(int), 1, fin );//landmark label number
    vertices[i].resize(C[i]);
    distances[i].resize(C[i]);
//    if(i == 2485)
//        cout<<i<<endl;
    for(int j = 0; j < C[i]; j++) {
        fread( &vertices[i][j], sizeof(int), 1, fin );//mapped landmark id
        fread( &distances[i][j], sizeof(int), 1, fin );//distance
//        int temp_dis = Dijkstra(i,landmarkTopK[vertices[i][j]],adj);
//        if(temp_dis!=distances[i][j]){
//            cout<<"Distance wrong! "<<i<<" "<<j<<": "<<distances[i][j]<<" "<<temp_dis<<endl;
//        }
    }
  }
    fclose(fin);

  file = std::string(filename) + "/" + dataset + "."+ std::to_string(K) + std::string(".highway");
    fin = fopen( file.c_str(), "rb" );
    if(fin == NULL){
        cout << "Failed to open file " << file << endl;
        exit(1);
    }
  highway.assign(K,vector<int>(K,INF));
  for(int i = 0; i < K; i++) {
    for(int j = 0; j < K; j++)
        fread( &highway[i][j], sizeof(int), 1, fin );
  }
    fclose(fin);

}
//function of saving landmarks
void HighwayLabellingWW::SaveLandmarks(string filename, vector<int> & topk){
    string file = filename + "/" + dataset + "." + to_string(K) + ".landmark";
    ofstream outFile(file, ios::out);
    if (!outFile) {
        cout << "File opening failed." << endl;
        exit(0);
    }
    cout << "Saving landmarks..." << endl;

    int p_i = 0;
    outFile << K << endl;
    for(int i=0;i<K;++i){
        outFile << topk[i] << endl;
        ++p_i;
    }
    outFile.close();
}
//function of saving landmarks
void HighwayLabellingWW::LoadLandmarks(string filename, vector<int> & topk){
    string file = filename + "/" + dataset + "." + to_string(K) + ".landmark";
    ifstream inFile(file, ios::in);
    if (!inFile) {
        cout << "File opening failed." << endl;
        exit(0);
    }
    cout << "Loading landmarks..." << endl;

    int p_i;
    inFile >> p_i;
    for(int i=0;i<p_i;++i){
        inFile >> topk[i];
        landmarkTopK[i]=topk[i];
        landmarks[topk[i]] = i;
        landmarkSet.insert(topk[i]);
    }
    inFile.close();
}
//function of BFS
int HighwayLabellingWW::BFS(int ID1, int ID2, vector<vector<int> > & graph){
    if(ID1 == ID2){
        return 0;
    }
    queue<int> queue;
    int item_id=0, temp_id=0;
    int item_dis,temp_dis;
    vector<bool> visited(V, false); //flag vector of whether closed
    vector<int> distance(V,INF);
    //Initiation of start node
    queue.push(ID1);
    distance[ID1]=0;
    visited[ID1]=true;
    assert(ID1 < V);
    int result_dis = INF;

    //Iteration
    while (!queue.empty()) {//for every node in queue !queue.empty()
        item_id = queue.front();// top and delete min item
        item_dis = distance[item_id];
        queue.pop();
        if(item_id == ID2){
            result_dis = item_dis;
            break;
        }
        //relaxation
        for (auto it = graph[item_id].begin(); it != graph[item_id].end(); ++it) {
            temp_id = *it; temp_dis = item_dis+1;
            if (!visited[temp_id]) {//if closed
                queue.push(temp_id);
                visited[temp_id] = true;
//                if(distance[temp_id] > temp_dis)
                    distance[temp_id] = temp_dis;
            }
        }
    }
    return result_dis;
}
//function of BFS
int HighwayLabellingWW::BFS(int ID1, int ID2, vector<unordered_map<int,int> > & graph){
    if(ID1 == ID2){
        return 0;
    }
    queue<int> queue;
    int item_id=0, temp_id=0;
    int item_dis,temp_dis;
    vector<bool> visited(V, false); //flag vector of whether closed
    vector<int> distance(V,INF);
    //Initiation of start node
    queue.push(ID1);
    distance[ID1]=0;
    visited[ID1]=true;
    assert(ID1 < V);
    int result_dis = INF;

    //Iteration
    while (!queue.empty()) {//for every node in queue !queue.empty()
        item_id = queue.front();// top and delete min item
        item_dis = distance[item_id];
        queue.pop();
        if(item_id == ID2){
            result_dis = item_dis;
            break;
        }
        //relaxation
        for (auto it = graph[item_id].begin(); it != graph[item_id].end(); ++it) {
            temp_id = it->first; temp_dis = item_dis+it->second;
            if (!visited[temp_id]) {//if closed
                queue.push(temp_id);
                visited[temp_id] = true;
//                if(distance[temp_id] > temp_dis)
                distance[temp_id] = temp_dis;
            }
        }
    }
    return result_dis;
}
//function of Bi-directional BFS
int HighwayLabellingWW::Bi_BFS(int ID1, int ID2, vector<vector<int> > & graph){
    if(ID1 == ID2){
        return 0;
    }
    queue<int> queue;
    int item_id=0, temp_id=0;
    int item_dis,temp_dis;
    vector<bool> visited(V, false); //flag vector of whether closed
    vector<int> distance(V,INF);
    //Initiation of start node
    queue.push(ID1);
    distance[ID1]=0;
    visited[ID1]=true;
    assert(ID1 < V);
    int result_dis = INF;

    //Iteration
    while (!queue.empty()) {//for every node in queue !queue.empty()
        item_id = queue.front();// top and delete min item
        item_dis = distance[item_id];
        queue.pop();
        if(item_id == ID2){
            result_dis = item_dis;
            break;
        }
        //relaxation
        for (auto it = graph[item_id].begin(); it != graph[item_id].end(); ++it) {
            temp_id = *it; temp_dis = item_dis+1;
            if (!visited[temp_id]) {//if closed
                queue.push(temp_id);
                visited[temp_id] = true;
//                if(distance[temp_id] > temp_dis)
                distance[temp_id] = temp_dis;
            }
        }
    }
    return result_dis;
}
//function of Dijkstra's algorithm
int HighwayLabellingWW::Dijkstra(int ID1, int ID2, vector<unordered_map<int,int>> & Nodes){
    benchmark::heap<2, int, int> pqueue(V);
    pqueue.update(ID1,0);

    vector<bool> closed(V, false);
    vector<int> distance(V, INF);
//	vector<int> prece(nodenum, 0);
    distance[ID1]=0;
    int topNodeID, topNodeDis;
    int NNodeID,NWeigh;

    int d=INF;//initialize d to infinite for the unreachable case

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        if(topNodeID==ID2){
            d=distance[ID2];
            break;
        }
        closed[topNodeID]=true;

        for(auto it=Nodes[topNodeID].begin();it!=Nodes[topNodeID].end();it++){
            NNodeID=(*it).first;
            NWeigh=(*it).second+topNodeDis;
            if(!closed[NNodeID]){
                if(distance[NNodeID]>NWeigh){
                    distance[NNodeID]=NWeigh;
                    pqueue.update(NNodeID, NWeigh);
                    //	prece[NNodeID]=topNodeID;
                }
            }
        }
    }

    return d;
}
void HighwayLabellingWW::BiDijkstraP(int ID1, int ID2, vector<unordered_map<int,int>> & Nodes, int & dis){// second version, powered by benchmark::heap
    if (ID1 == ID2){
        dis = 0;
        return;
    }
    benchmark::heap<2, int, int> pqueue(V), pqueue_b(V);
    int item_id, item_id_b, temp_id, temp_id_b;
    int item_dis, item_dis_b, temp_dis, temp_dis_b;
    int terminate_id;//termination id of bi-dijkstra
    vector<int> cost(V, INF);         //cost for current node to target node
    vector<bool> closed(V, false); //flag of whether having been closed
    vector<int> cost_b(V, INF);         //cost for current node to target node
    vector<bool> closed_b(V, false); //flag of whether having been closed
    vector<int> pre(V, -1); //forward search vector of predecessor vertex id
    vector<int> pre_b(V, -1); //backward search vector of predecessor vertex id
    int min_cost = INF;

    //Initiation of start node
    cost[ID1] = 0;//cost of start node
    closed[ID1] = true;
    cost_b[ID2] = 0;//cost of start node
    closed_b[ID2] = true;
    pre[ID1] = ID1;
    pre_b[ID2] = ID2;
    pqueue.update(ID1, 0);
    pqueue_b.update(ID2, 0);

    //Iteration
    while (!pqueue.empty() && !pqueue_b.empty()) {//if either priority queue is not empty
        if (pqueue.top_key() + pqueue_b.top_key() >= min_cost) {//condition of termination
            break;
        }
        /*-------Forward Search-------*/
        pqueue.extract_min(item_id, item_dis);// top min item
        closed[item_id] = true;
        //relax the adjacent nodes of node item
        for (auto it = Nodes[item_id].begin(); it != Nodes[item_id].end(); ++it) {
            temp_id = it->first;
            temp_dis = item_dis + it->second;
            if (!closed[temp_id]) {//if not visited
                if (cost[temp_id] > temp_dis) {//relaxation
                    cost[temp_id] = temp_dis;
                    pre[temp_id] = item_id;
                    pqueue.update(temp_id, cost[temp_id]);
                }
            }
            if (closed_b[temp_id] && temp_dis + cost_b[temp_id] < min_cost) {
                min_cost = temp_dis + cost_b[temp_id];
                terminate_id = temp_id;
//                cout<<"Forward! min_cost: "<<min_cost<<", terminate_id: "<<terminate_id<<endl;
            }
        }
        /*-------Reverse Search-------*/
        pqueue_b.extract_min(item_id_b, item_dis_b);// top min item
        closed_b[item_id_b] = true;

        //relax the adjacent nodes of node item
        for (auto it = Nodes[item_id_b].begin(); it != Nodes[item_id_b].end(); ++it) {
            temp_id_b = it->first;
            temp_dis_b = item_dis_b + it->second;
            if (!closed_b[temp_id_b]) {//if not closed
                if (cost_b[temp_id_b] > temp_dis_b) {//slack operation
                    cost_b[temp_id_b] = temp_dis_b;
                    pre_b[temp_id_b] = item_id_b;
                    pqueue_b.update(temp_id_b, cost_b[temp_id_b]);
                }
            }
            if (closed[temp_id_b] && temp_dis_b + cost[temp_id_b] < min_cost) {
                min_cost = temp_dis_b + cost[temp_id_b];
                terminate_id = temp_id_b;
//                cout<<"Reverse! min_cost: "<<min_cost<<", terminate_id: "<<terminate_id<<endl;
            }
        }
    }
    dis = min_cost;
}
//Function of in-memory Bidirectional Dijkstra's algorithm
int HighwayLabellingWW::BiDijkstra(int ID1, int ID2, vector<unordered_map<int,int>> & Nodes) {// second version, powered by benchmark::heap
    if (ID1 == ID2) return 0;
    benchmark::heap<2, int, int> pqueue(V), pqueue_b(V);
    int item_id, item_id_b, temp_id, temp_id_b;
    int item_dis, item_dis_b, temp_dis, temp_dis_b;
    int terminate_id;//termination id of bi-dijkstra
    vector<int> cost(V, INF);         //cost for current node to target node
    vector<bool> closed(V, false); //flag of whether having been closed
    vector<int> cost_b(V, INF);         //cost for current node to target node
    vector<bool> closed_b(V, false); //flag of whether having been closed
    vector<int> pre(V, -1); //forward search vector of predecessor vertex id
    vector<int> pre_b(V, -1); //backward search vector of predecessor vertex id
    int min_cost = INF;

    //Initiation of start node
    cost[ID1] = 0;//cost of start node
    closed[ID1] = true;
    cost_b[ID2] = 0;//cost of start node
    closed_b[ID2] = true;
    pre[ID1] = ID1;
    pre_b[ID2] = ID2;
    pqueue.update(ID1, 0);
    pqueue_b.update(ID2, 0);

    //Iteration
    while (!pqueue.empty() && !pqueue_b.empty()) {//if either priority queue is not empty
        if (pqueue.top_key() + pqueue_b.top_key() >= min_cost) {//condition of termination
            break;
        }
        /*-------Forward Search-------*/
        pqueue.extract_min(item_id, item_dis);// top min item
        closed[item_id] = true;
        //relax the adjacent nodes of node item
        for (auto it = Nodes[item_id].begin(); it != Nodes[item_id].end(); ++it) {
            temp_id = it->first;
            temp_dis = item_dis + it->second;
            if (!closed[temp_id]) {//if not visited
                if (cost[temp_id] > temp_dis) {//relaxation
                    cost[temp_id] = temp_dis;
                    pre[temp_id] = item_id;
                    pqueue.update(temp_id, cost[temp_id]);
                }
            }
            if (closed_b[temp_id] && temp_dis + cost_b[temp_id] < min_cost) {
                min_cost = temp_dis + cost_b[temp_id];
                terminate_id = temp_id;
//                cout<<"Forward! min_cost: "<<min_cost<<", terminate_id: "<<terminate_id<<endl;
            }
        }
        /*-------Reverse Search-------*/
        pqueue_b.extract_min(item_id_b, item_dis_b);// top min item
        closed_b[item_id_b] = true;

        //relax the adjacent nodes of node item
        for (auto it = Nodes[item_id_b].begin(); it != Nodes[item_id_b].end(); ++it) {
            temp_id_b = it->first;
            temp_dis_b = item_dis_b + it->second;
            if (!closed_b[temp_id_b]) {//if not closed
                if (cost_b[temp_id_b] > temp_dis_b) {//slack operation
                    cost_b[temp_id_b] = temp_dis_b;
                    pre_b[temp_id_b] = item_id_b;
                    pqueue_b.update(temp_id_b, cost_b[temp_id_b]);
                }
            }
            if (closed[temp_id_b] && temp_dis_b + cost[temp_id_b] < min_cost) {
                min_cost = temp_dis_b + cost[temp_id_b];
                terminate_id = temp_id_b;
//                cout<<"Reverse! min_cost: "<<min_cost<<", terminate_id: "<<terminate_id<<endl;
            }
        }
    }
    return min_cost;
}
//Function of in-memory Bidirectional Dijkstra's algorithm
int HighwayLabellingWW::BiDijkstraBounded(int ID1, int ID2, vector<unordered_map<int,int>> & Nodes, int dis_upper) {// second version, powered by benchmark::heap
    if (ID1 == ID2) return 0;
    benchmark::heap<2, int, int> pqueue(V), pqueue_b(V);
    int item_id, item_id_b, temp_id, temp_id_b;
    int item_dis, item_dis_b, temp_dis, temp_dis_b;
    int terminate_id;//termination id of bi-dijkstra
    vector<int> cost(V, INF);         //cost for current node to target node
    vector<bool> closed(V, false); //flag of whether having been closed
    vector<int> cost_b(V, INF);         //cost for current node to target node
    vector<bool> closed_b(V, false); //flag of whether having been closed
    vector<int> pre(V, -1); //forward search vector of predecessor vertex id
    vector<int> pre_b(V, -1); //backward search vector of predecessor vertex id
    int min_cost = INF;

    //Initiation of start node
    cost[ID1] = 0;//cost of start node
    closed[ID1] = true;
    cost_b[ID2] = 0;//cost of start node
    closed_b[ID2] = true;
    pre[ID1] = ID1;
    pre_b[ID2] = ID2;
    pqueue.update(ID1, 0);
    pqueue_b.update(ID2, 0);

    //Iteration
    while (!pqueue.empty() && !pqueue_b.empty()) {//if either priority queue is not empty
        if (pqueue.top_key() + pqueue_b.top_key() >= min_cost) {//condition of termination
            break;
        }
        if (pqueue.top_key() + pqueue_b.top_key() >= dis_upper) {//condition of termination
            min_cost = pqueue.top_key() + pqueue_b.top_key();
            break;
        }
        /*-------Forward Search-------*/
        pqueue.extract_min(item_id, item_dis);// top min item
        closed[item_id] = true;
        //relax the adjacent nodes of node item
        for (auto it = Nodes[item_id].begin(); it != Nodes[item_id].end(); ++it) {
            temp_id = it->first;
            temp_dis = item_dis + it->second;
            if (!closed[temp_id]) {//if not visited
                if (cost[temp_id] > temp_dis) {//relaxation
                    cost[temp_id] = temp_dis;
                    pre[temp_id] = item_id;
                    pqueue.update(temp_id, cost[temp_id]);
                }
            }
            if (closed_b[temp_id] && temp_dis + cost_b[temp_id] < min_cost) {
                min_cost = temp_dis + cost_b[temp_id];
                terminate_id = temp_id;
//                cout<<"Forward! min_cost: "<<min_cost<<", terminate_id: "<<terminate_id<<endl;
            }
        }
        /*-------Reverse Search-------*/
        pqueue_b.extract_min(item_id_b, item_dis_b);// top min item
        closed_b[item_id_b] = true;

        //relax the adjacent nodes of node item
        for (auto it = Nodes[item_id_b].begin(); it != Nodes[item_id_b].end(); ++it) {
            temp_id_b = it->first;
            temp_dis_b = item_dis_b + it->second;
            if (!closed_b[temp_id_b]) {//if not closed
                if (cost_b[temp_id_b] > temp_dis_b) {//slack operation
                    cost_b[temp_id_b] = temp_dis_b;
                    pre_b[temp_id_b] = item_id_b;
                    pqueue_b.update(temp_id_b, cost_b[temp_id_b]);
                }
            }
            if (closed[temp_id_b] && temp_dis_b + cost[temp_id_b] < min_cost) {
                min_cost = temp_dis_b + cost[temp_id_b];
                terminate_id = temp_id_b;
//                cout<<"Reverse! min_cost: "<<min_cost<<", terminate_id: "<<terminate_id<<endl;
            }
        }
    }
    return min_cost;
}
//function to check the correctness of BatchHL
void HighwayLabellingWW::CorrectnessCheck(int runtimes,bool ifFull){
    int s = 0, t = 0; int total = 0;
    cout<<"Run times: "<<runtimes<<endl;
    int disDijk, disBatchHL;
    double time=0;
    Timer tt;
    for(int i=0;i<runtimes;++i){//runtimes
        s = rand()%V; t = rand()%V;
//        s = 3291; t = 1;
        tt.start();
        disBatchHL = BatchHL_Query(s,t,adj,ifFull);
        tt.stop();
        time += tt.GetRuntime();
        disDijk = Dijkstra(s,t,adjOld);
//        disDijk = Dijkstra(s,t,adj);
        if(disBatchHL != disDijk){
            cout << "Incorrect: " << s << " " << t << " " << disBatchHL << " " << disDijk <<" "<<Dijkstra(s,t,adjOldO)<<endl;
        }
    }

    std::cout << "Average Query Time (ms) : " << time*1000 / runtimes << " ms."<<std::endl;

}
//function to check the correctness of BatchHL
void HighwayLabellingWW::CorrectnessCheckW(int runtimes,bool ifFull){
    int s = 0, t = 0; int total = 0;
    cout<<"Run times: "<<runtimes<<endl;
    int disDijk, disBatchHL;
    double time=0;
    Timer tt;
    for(int i=0;i<runtimes;++i){//runtimes
        s = rand()%V; t = rand()%V;
//        s = 2255; t = 3267;
        tt.start();
        disBatchHL = BatchHL_QueryW(s,t,adj,ifFull);
        tt.stop();
        time += tt.GetRuntime();
        disDijk = Dijkstra(s,t,adjOld);
        if(disBatchHL != disDijk){
            cout << "Incorrect: " << s << " " << t << " " << disBatchHL << " " << disDijk <<" "<<Dijkstra(s,t,adjOldO)<<endl;
        }
    }

    std::cout << "Average Query Time (ms) : " << time*1000 / runtimes << " ms."<<std::endl;

}

//function of checking the connectivity of graph: map version
void HighwayLabellingWW::DFS_CC(vector<unordered_map<int,int>> & Edges, set<int> & set_A) {
    /// DFS for connected component
    stack<int> stack_A;
//    set<int> set_A;//nodes waiting for visited
    set<int> set_B;//nodes visited for current component
    int item_id,temp_id;
    vector<bool> flag_visited(V,false);
    bool flag_finish = false;
    int temp_num = 0;
    int component_i = 0;
    pair<set<int>,long long int> mcc;

//    for(int i=0;i<nodenum;++i){
//        set_A.insert(i);
//    }
    int seed = *set_A.begin();
    stack_A.push(seed);
    set_A.erase(seed);
    set_B.insert(seed);
    flag_visited[seed] = true;
    //Get the connected components by DFS
    while(!set_A.empty()) {//if not finish
        temp_num = 0;
        while (!stack_A.empty()) {
            item_id = stack_A.top();
            stack_A.pop();
            for (auto it = Edges[item_id].begin(); it != Edges[item_id].end(); ++it) {
                temp_id = it->first;
                temp_num += 1;
                if (!flag_visited[temp_id]) {//if not visited
                    stack_A.push(temp_id);
                    set_A.erase(temp_id);
                    set_B.insert(temp_id);
                    flag_visited[temp_id] = true;
                }
            }
        }
        if (set_B.size() > mcc.first.size()) {
            mcc.first.clear();
            mcc.first = set_B;
            mcc.second = temp_num;// /2
        }
        if(!set_B.empty() && set_B.size() < mcc.first.size()){
//            cout<<"Smaller connected component with vertex size "<<set_B.size()<<": ";
//            for(auto it=set_B.begin();it!=set_B.end();++it){
//                cout<<*it<<" ";
//            }
//            cout<<"; degree: ";
//            for(auto it=set_B.begin();it!=set_B.end();++it){
//                cout<<Edges[*it].size()<<" ";
//            }
//            cout<<endl;
        }
        ++component_i;
        set_B.clear();
        if (!set_A.empty()) {
            stack_A.push(*set_A.begin());
            set_B.insert(*set_A.begin());
            flag_visited[*set_A.begin()] = true;
            set_A.erase(*set_A.begin());
        } else {
            break;
        }
    }
    if(component_i==1){
        cout<<"This graph has only one connected component."<<endl;
        cout<<"Nodes size of graph: "<< mcc.first.size() << endl;
        cout<<"Edges size of graph: "<< mcc.second << endl;
    }else{
        cout<<"This graph has "<< component_i <<" connected component!"<<endl;
        cout<<"Nodes size of the largest connected component is: "<<mcc.first.size()<<endl;
        cout<<"Edges size of the largest connected component is: "<<mcc.second<<endl;
    }
}

//function of checking the connectivity of graph
void HighwayLabellingWW::DFS_CC(vector<vector<int>> & Edges, set<int> & set_A) {
    /// DFS for connected component
    stack<int> stack_A;
//    set<int> set_A;//nodes waiting for visited
    set<int> set_B;//nodes visited for current component
    int item_id,temp_id;
    vector<bool> flag_visited(V,false);
    bool flag_finish = false;
    int temp_num = 0;
    int component_i = 0;
    pair<set<int>,int> mcc;

//    for(int i=0;i<nodenum;++i){
//        set_A.insert(i);
//    }
    int seed = *set_A.begin();
    stack_A.push(seed);
    set_A.erase(seed);
    set_B.insert(seed);
    flag_visited[seed] = true;
    //Get the connected components by DFS
    while(!set_A.empty()) {//if not finish
        temp_num = 0;
        while (!stack_A.empty()) {
            item_id = stack_A.top();
            stack_A.pop();
            for (auto it = Edges[item_id].begin(); it != Edges[item_id].end(); ++it) {
                temp_id = *it;
                temp_num += 1;
                if (!flag_visited[temp_id]) {//if not visited
                    stack_A.push(temp_id);
                    set_A.erase(temp_id);
                    set_B.insert(temp_id);
                    flag_visited[temp_id] = true;
                }
            }
        }
        if (set_B.size() > mcc.first.size()) {
            mcc.first.clear();
            mcc.first = set_B;
            mcc.second = temp_num;// /2
        }
        if(!set_B.empty() && set_B.size() < mcc.first.size()){
            cout<<"Smaller connected component with vertex size "<<set_B.size()<<": ";
            for(auto it=set_B.begin();it!=set_B.end();++it){
                cout<<*it<<" ";
            }
            cout<<"; degree: ";
            for(auto it=set_B.begin();it!=set_B.end();++it){
                cout<<Edges[*it].size()<<" ";
            }
            cout<<endl;
        }
        ++component_i;
        set_B.clear();
        if (!set_A.empty()) {
            stack_A.push(*set_A.begin());
            set_B.insert(*set_A.begin());
            flag_visited[*set_A.begin()] = true;
            set_A.erase(*set_A.begin());
        } else {
            break;
        }
    }
    if(component_i==1){
        cout<<"This graph has only one connected component."<<endl;
        cout<<"Nodes size of graph: "<< mcc.first.size() << endl;
        cout<<"Edges size of graph: "<< mcc.second << endl;
    }else{
        cout<<"This graph has "<< component_i <<" connected component!"<<endl;
        cout<<"Nodes size of the largest connected component is: "<<mcc.first.size()<<endl;
        cout<<"Edges size of the largest connected component is: "<<mcc.second<<endl;
    }
}

void HighwayLabellingWW::CheckCC() {
    /// check graph connectivity
    set<int> vSet;
    for(int i=0;i<V;++i){
        vSet.insert(i);
    }
    DFS_CC(adj,vSet);
}

/*int HighwayLabellingWW::BatchHL_Query(int s, int t,vector<vector<int> > & graph,bool ifFull) {
    if(s == t)
        return 0;
    std::vector<TwoLayerQueue> qque; std::vector<int> qdist[2];

    qdist[0].resize(V, INF); qdist[1].resize(V, INF);
    qque.push_back(TwoLayerQueue(V)); qque.push_back(TwoLayerQueue(V));

    int dist_upper;
    if(!ifFull) {
        dist_upper = HC_UB_naive(s, t);//distance bound
    }else{
        dist_upper = HC_UB_naive_Full(s,t);
    }

//    cout<<"dist_upper: "<<dist_upper<<endl;
//    cout<<"BFS on sparsified graph: "<<BFS(s,t,adj)<<endl;

    int res = dist_upper, dis[2] = {0, 0};
    for (int dir = 0; dir < 2; dir++){
        int v = dir == 0 ? s : t;
        qque[dir].clear();
        qque[dir].push(v);
        qque[dir].next();
        qdist[dir][v] = 0;
    }

    while (!qque[0].empty() && !qque[1].empty()) {
        int use = (qque[0].size() <= qque[1].size()) ? 0 : 1;
        dis[use]++;

        if (dis[0] + dis[1] == dist_upper) {
            res = dis[0] + dis[1];
            goto LOOP_END;
        }

        while (!qque[use].empty()) {

            int v = qque[use].front();
            qque[use].pop();

            for (int w : graph[v]) {

                int &src_d = qdist[    use][w];
                int &dst_d = qdist[1 - use][w];
                if (src_d != INF) continue;
                if (dst_d != INF) {
                    res = qdist[use][v] + 1 + dst_d;
                    goto LOOP_END;
                } else {
                    qque[use].push(w);
                    qdist[use][w] = qdist[use][v] + 1;
                }
            }
        }
        qque[use].next();
    }
    LOOP_END:
    for (int dir = 0; dir < 2; dir++) {
        for (int v : qque[dir]) {
            qdist[dir][v] = INF;
        }
        qque[dir].clear();
    }
    return (int) min(res, dist_upper);
}*/

int HighwayLabellingWW::BatchHL_Query(int s, int t,vector<unordered_map<int,int> > & graph,bool ifFull) {
    if(s == t)
        return 0;
    int dist_upper, res = INF;
    bool flag_ll = false;
    if(landmarkSet.find(s)!=landmarkSet.end()){
//        cout<<s<<" is a landmark; ";
        flag_ll = true;
    }else{
//        cout<<s<<" is not a landmark; ";
        flag_ll = false;
    }
    if(landmarkSet.find(t)!=landmarkSet.end()){
//        cout<<t<<" is a landmark."<<endl;
        if(!flag_ll)
            flag_ll = false;
    }else{
//        cout<<t<<" is not a landmark."<<endl;
        flag_ll = false;
    }
    if(flag_ll){
        assert(landmarks.find(s) != landmarks.end());
        assert(landmarks.find(t) != landmarks.end());
        res = highway[landmarks[s]][landmarks[t]];
    }else{
        int bi_dis;
//        boost::thread_group threadf;
//        threadf.add_thread(new boost::thread(&HighwayLabellingWW::HC_UB_naiveP, this, s,  t, ifFull, dist_upper));
//        threadf.add_thread(new boost::thread(&HighwayLabellingWW::BiDijkstraP, this, s, t, boost::ref(adj),boost::ref(bi_dis)));
//        threadf.join_all();
        if(!ifFull) {
            dist_upper = HC_UB_naive(s, t);//distance bound
        }else{
            dist_upper = HC_UB_naive_Full(s,t);
        }
//        bi_dis = BiDijkstra(s,t,adj);
        bi_dis = BiDijkstraBounded(s,t,adj,dist_upper);
//        cout<<"dist_upper: "<<dist_upper<<endl;
//        cout<<"BiDijkstra on sparsified graph: "<<bi_dis<<endl;
        res =  min(bi_dis, dist_upper);
    }

//    cout<<"Dijkstra on sparsified graph: "<<Dijkstra(s,t,adj)<<endl;
//    cout<<"BFS on sparsified graph: "<<BFS(s,t,adj)<<endl;

    return res;
}

int HighwayLabellingWW::BatchHL_QueryW(int s, int t,vector<unordered_map<int,int> > & graph,bool ifFull) {
    if(s == t)
        return 0;
    int dist_upper, res = INF;
    bool flag_ll = false;
    if(landmarkSet.find(s)!=landmarkSet.end()){
//        cout<<s<<" is a landmark; ";
        flag_ll = true;
    }else{
//        cout<<s<<" is not a landmark; ";
        flag_ll = false;
    }
    if(landmarkSet.find(t)!=landmarkSet.end()){
//        cout<<t<<" is a landmark."<<endl;
        if(!flag_ll)
            flag_ll = false;
    }else{
//        cout<<t<<" is not a landmark."<<endl;
        flag_ll = false;
    }
    if(flag_ll){
        assert(landmarks.find(s) != landmarks.end());
        assert(landmarks.find(t) != landmarks.end());
        res = highway[landmarks[s]][landmarks[t]];
    }else{
        if(!ifFull) {
            dist_upper = HC_UB_naive(s, t);//distance bound
        }else{
            dist_upper = HC_UB_naive_Full(s,t);
        }
        int bi_dis = BiDijkstra(s,t,adj);
//        cout<<"dist_upper: "<<dist_upper<<endl;
//        cout<<"BiDijkstra on sparsified graph: "<<bi_dis<<endl;
        res =  min(bi_dis, dist_upper);
    }

//    cout<<"Dijkstra on sparsified graph: "<<Dijkstra(s,t,adj)<<endl;
//    cout<<"BFS on sparsified graph: "<<BFS(s,t,adj)<<endl;

    return res;
}

//function of reading update edges
void HighwayLabellingWW::ReadUpdates(string filename){
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

void HighwayLabellingWW::verticesUpdate(){
    vertices.resize(V);
    C.resize(V);
    for(int i = 0; i < V; i++) {
        int Coun = 0;
        for (int j = 0; j < K; j++) {
            if (distances[i][j] != INF)//if there is distance label from vertex i to landmark j
                Coun++;
        }
        if (Coun == 0) {
            cout << "Error! C == 0 : " << i << endl;
        }
        C[i] = Coun;
        vertices[i].resize(Coun);
        int count = 0;
        for (int j = 0; j < K; j++) {
            if (distances[i][j] != INF) {//if there is distance label from vertex i to landmark j
                vertices[i][count] = j;
                ++count;
            }
        }
        assert(count == Coun);
    }

}


#endif  // BATCHHLWW_HPP_
