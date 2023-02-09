/*
 * CH.cpp
 *
 *  Created on: 14 Jan 2023
 *      Author: zhangmengxuan, Xinjie ZHOU
 */
#include "head3.h"
//CH query for OD within one periphery, single source
void Graph::QueryPeriphery_CH(int ID1, set<int> todo, map<int,int> & results){
    //priority queue
    benchmark::heap<2,int,int> pqueue(nodenum);

    //closed or not
    vector<bool> vVisited(nodenum, false);
    //the existing shortest distance
    vector<int>	vDistance(nodenum, INF);

    vDistance[ID1] = 0;
    pqueue.update(ID1,0);
    int topID,topDis,tempID,tempDis,weight;

    while( ! todo.empty() && ! pqueue.empty() ){
        pqueue.extract_min(topID, topDis);

        vVisited[topID] = true;

        if(todo.find(topID) != todo.end()){//if found
            todo.erase(topID);
            results[topID] = topDis;
        }

        for(auto out=NeighborConCH[topID].begin();out!=NeighborConCH[topID].end();++out){
            tempID = (*out).first;
            weight = (*out).second.first;


            if(!vVisited[tempID]){
                tempDis = vDistance[topID] + weight;
                if(vDistance[tempID] > tempDis){
                    vDistance[tempID] = tempDis;
                    pqueue.update(tempID, tempDis);
                }
            }
        }
    }

}

//CH query for OD within one periphery
int Graph::QueryPeriphery_CH(int ID1, int ID2){
	if(ID1==ID2) return 0;
	if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
	int d=INF;
	benchmark::heap<2,int,int> fHeapForward(nodenum);
	benchmark::heap<2, int, int> fHeapBackward(nodenum);

	//closed or not
	vector<bool> vVisitedF(nodenum, false);
	vector<bool> vVisitedB(nodenum, false);
	//the existing shortest distance
	vector<int>	vDistanceForward(nodenum, INF);
	vector<int>	vDistanceBackward(nodenum, INF);
	//stop search or not
	bool bF = false;
	bool bB = false;
	vDistanceForward[ID1] = 0;
	vDistanceBackward[ID2] = 0;
	fHeapForward.update(ID1,0);
	fHeapBackward.update(ID2,0);

	int topNodeIDForward, topNodeIDBackward,topDisForward,topDisBackward, neighborNodeID, neighborLength;

	while(!fHeapForward.empty() || !fHeapBackward.empty() )
	{
		if(bF && bB)
			break;
		if(bF && fHeapBackward.empty())
			break;
		if(bB && fHeapForward.empty())
			break;
		//Forward Search
		if(!fHeapForward.empty() && !bF)
		{
			fHeapForward.extract_min(topNodeIDForward, topDisForward);

			if(vDistanceForward[topNodeIDForward] > d)
				bF = true;

			vVisitedF[topNodeIDForward] = true;

			if(vVisitedB[topNodeIDForward]){
				int distTmp=topDisForward+vDistanceBackward[topNodeIDForward];
				if(distTmp<d){
					d=distTmp;
				}
			}

			for(auto out=NeighborConCH[topNodeIDForward].begin();out!=NeighborConCH[topNodeIDForward].end();out++){
				neighborNodeID = (*out).first;
				neighborLength = (*out).second.first;

				int df = vDistanceForward[topNodeIDForward] + neighborLength;
				if(!vVisitedF[neighborNodeID]){
					if(vDistanceForward[neighborNodeID] > df){
						vDistanceForward[neighborNodeID] = df;
						fHeapForward.update(neighborNodeID, df);
					}
				}
			}
		}

		//Backward Search
		if(!fHeapBackward.empty() && !bB)
		{
			fHeapBackward.extract_min(topNodeIDBackward, topDisBackward);

			if(vDistanceBackward[topNodeIDBackward] > d)
				bB = true;

			vVisitedB[topNodeIDBackward] = true;

			if(vVisitedF[topNodeIDBackward]){
				int distTmp=topDisBackward+vDistanceForward[topNodeIDBackward];
				if(distTmp<d){
					d=distTmp;
				}
			}

			for(auto in=NeighborConCH[topNodeIDBackward].begin();in!=NeighborConCH[topNodeIDBackward].end();in++){
				neighborNodeID = (*in).first;
				neighborLength = (*in).second.first;

				int db = vDistanceBackward[topNodeIDBackward] + neighborLength;
				if(!vVisitedB[neighborNodeID]){
					if(vDistanceBackward[neighborNodeID]>db){
						vDistanceBackward[neighborNodeID] = db;
						fHeapBackward.update(neighborNodeID, db);
					}
				}
			}
		}
	}
	return d;
}
//new function
int Graph::QuerySameParti_CH(int ID1, int ID2){//both are within the partition
    int d=INF;

    int pid1=CoreTag[ID1];
    int pid2=CoreTag[ID2];
    if(pid1==pid2){//if in the same partition
        //cout<<"Same-Parti"<<endl;
        int temp_dis = QueryPeriphery_CH(ID1,ID2);/// d2 may be wrong sometimes
        if(temp_dis<d)//QueryH2HPartition(ID1,ID2,pid1)
            d=temp_dis;//QueryH2HPartition(ID1,ID2,pid1);
        vector<int> B=BoundVertex[pid1];
        map<int,int> m1,m2;
        m1.clear();
        m2.clear();
        vector<int> B1,B2;
        B1.clear();
        B2.clear();
        int bID,d1,d2;
        set<int> todo; todo.clear();
        todo.insert(B.begin(),B.end()); todo.erase(B[B.size()-1]);
        map<int,int> result1,result2;

        QueryPeriphery_CH(ID1,todo,result1);
        QueryPeriphery_CH(ID2,todo,result2);

        for(int i=0;i<B.size()-1;i++){
            bID=B[i];
//           assert(Tree[rank[ID1]].disInf.find(bID)!=Tree[rank[ID1]].disInf.end());
//           assert(Tree[rank[ID2]].disInf.find(bID)!=Tree[rank[ID2]].disInf.end());
//            d1=Tree[rank[ID1]].disInf[i];
//            d2=Tree[rank[ID2]].disInf[i];
//            d1=Tree[rank[ID1]].disInf[bID];
//            d2=Tree[rank[ID2]].disInf[bID];

//            d1=QueryPeriphery_CH(ID1, bID);
//            d2=QueryPeriphery_CH(ID2, bID);
            d1=result1[bID];
            d2=result2[bID];

            if(d1<d){
                B1.push_back(bID);
                m1.insert(make_pair(bID,d1));
            }
            if(d2<d){
                B2.push_back(bID);
                m2.insert(make_pair(bID,d2));
            }
        }

        int bID1, bID2, tempdis;
        if(!B1.empty() && !B2.empty()){
            for(int k=0;k<B1.size();k++){
                bID1=B1[k];
                if(m1[bID1]>d)
                    continue;
                for(int z=0;z<B2.size();z++){
                    bID2=B2[z];
                    if(m2[bID2]>d)
                        continue;
                    tempdis=m1[bID1]+QueryCore(bID1,bID2)+m2[bID2];
                    if(tempdis<d)
                        d=tempdis;
                }
            }
        }

    }else{//if in different partitions
        cout<<"Wrong for same partition query!"<<endl;
        exit(1);
    }

    return d;
}
//old function
/*int Graph::QuerySameParti_CH(int ID1, int ID2){//both are within the partition
    int d=INF;

    int pid1=CoreTag[ID1];
    int pid2=CoreTag[ID2];
    if(pid1==pid2){//if in the same partition
        //cout<<"Same-Parti"<<endl;
        int temp_dis = QueryPeriphery_CH(ID1,ID2);/// d2 may be wrong sometimes
        if(temp_dis<d)//QueryH2HPartition(ID1,ID2,pid1)
            d=temp_dis;//QueryH2HPartition(ID1,ID2,pid1);
        vector<int> B=BoundVertex[pid1];
        map<int,int> m1,m2;
        m1.clear();
        m2.clear();
        vector<int> B1,B2;
        B1.clear();
        B2.clear();
        int bID,d1,d2;
        for(int i=0;i<B.size()-1;i++){
            bID=B[i];
//           assert(Tree[rank[ID1]].disInf.find(bID)!=Tree[rank[ID1]].disInf.end());
//           assert(Tree[rank[ID2]].disInf.find(bID)!=Tree[rank[ID2]].disInf.end());
//            d1=Tree[rank[ID1]].disInf[i];
//            d2=Tree[rank[ID2]].disInf[i];
//            d1=Tree[rank[ID1]].disInf[bID];
//            d2=Tree[rank[ID2]].disInf[bID];

            d1=QueryPeriphery_CH(ID1, bID);
            d2=QueryPeriphery_CH(ID2, bID);

            if(d1<d){
                B1.push_back(bID);
                m1.insert(make_pair(bID,d1));
            }
            if(d2<d){
                B2.push_back(bID);
                m2.insert(make_pair(bID,d2));
            }
        }

        int bID1, bID2, tempdis;
        if(!B1.empty() && !B2.empty()){
            for(int k=0;k<B1.size();k++){
                bID1=B1[k];
                if(m1[bID1]>d)
                    continue;
                for(int z=0;z<B2.size();z++){
                    bID2=B2[z];
                    if(m2[bID2]>d)
                        continue;
                    tempdis=m1[bID1]+QueryCore(bID1,bID2)+m2[bID2];
                    if(tempdis<d)
                        d=tempdis;
                }
            }
        }

    }else{//if in different partitions
        cout<<"Wrong for same partition query!"<<endl;
        exit(1);
    }

    return d;
}*/
//new function
int Graph::QueryPartiParti_CH(int ID1, int ID2){//both are within the partition
    int d=INF;

    int pid1=CoreTag[ID1];
    int pid2=CoreTag[ID2];
    if(pid1==pid2){//if in the same partition
        cout<<"Wrong for partition-partition query!"<<endl;
        exit(1);

    }else{//if in different partitions
        //cout<<"Parti-Parti"<<endl;
        vector<int> B1=BoundVertex[pid1];
        vector<int> B2=BoundVertex[pid2];

        map<int,int> m1,m2;
        m1.clear();
        m2.clear();
        int bID1, bID2, tempdis;
        int b1,b2,d1,d2;
        set<int> todo1,todo2; todo1.clear(); todo2.clear();
        todo1.insert(B1.begin(),B1.end()); todo2.insert(B2.begin(),B2.end());
        todo1.erase(B1[B1.size()-1]); todo2.erase(B2[B2.size()-1]);
        map<int,int> result1,result2;

        QueryPeriphery_CH(ID1,todo1,result1);
        for(int i=0;i<B1.size()-1;i++){
            bID1=B1[i];
//            assert(Tree[rank[ID1]].disInf.find(bID1)!=Tree[rank[ID1]].disInf.end());
//            m1.insert(make_pair(bID1,Tree[rank[ID1]].disInf[i]));
//            m1.insert(make_pair(bID1,Tree[rank[ID1]].disInf[bID1]));
//            m1.insert({bID1,QueryPeriphery_CH(ID1, bID1)});
            m1.insert({bID1,result1[bID1]});
        }
        QueryPeriphery_CH(ID2,todo2,result2);
        for(int j=0;j<B2.size()-1;j++){
            bID2=B2[j];
//            assert(Tree[rank[ID2]].disInf.find(bID2)!=Tree[rank[ID2]].disInf.end());
//            m2.insert(make_pair(bID2,Tree[rank[ID2]].disInf[j]));
//            m2.insert(make_pair(bID2,Tree[rank[ID2]].disInf[bID2]));
//            m2.insert({bID2,QueryPeriphery_CH(ID2, bID2)});
            m2.insert({bID2,result2[bID2]});
        }

        for(int k=0;k<B1.size()-1;k++){
            bID1=B1[k];

            if(m1[bID1]>d)
                continue;

            for(int z=0;z<B2.size()-1;z++){
                bID2=B2[z];

                if(m2[bID2]>d)
                    continue;

                tempdis=m1[bID1]+QueryCore(bID1,bID2)+m2[bID2];
                if(tempdis<d){
                    d=tempdis;
                    d1=m1[bID1]; d2=m2[bID2];
                    b1=bID1; b2=bID2;
                }

            }
        }

    }

    return d;
}
//old function
/*int Graph::QueryPartiParti_CH(int ID1, int ID2){//both are within the partition
    int d=INF;

    int pid1=CoreTag[ID1];
    int pid2=CoreTag[ID2];
    if(pid1==pid2){//if in the same partition
        cout<<"Wrong for partition-partition query!"<<endl;
        exit(1);

    }else{//if in different partitions
        //cout<<"Parti-Parti"<<endl;
        vector<int> B1=BoundVertex[pid1];
        vector<int> B2=BoundVertex[pid2];

        map<int,int> m1,m2;
        m1.clear();
        m2.clear();
        int bID1, bID2, tempdis;
        int b1,b2,d1,d2;
        for(int i=0;i<B1.size()-1;i++){
            bID1=B1[i];
//            assert(Tree[rank[ID1]].disInf.find(bID1)!=Tree[rank[ID1]].disInf.end());
//            m1.insert(make_pair(bID1,Tree[rank[ID1]].disInf[i]));
//            m1.insert(make_pair(bID1,Tree[rank[ID1]].disInf[bID1]));
            m1.insert({bID1,QueryPeriphery_CH(ID1, bID1)});
        }
        for(int j=0;j<B2.size()-1;j++){
            bID2=B2[j];
//            assert(Tree[rank[ID2]].disInf.find(bID2)!=Tree[rank[ID2]].disInf.end());
//            m2.insert(make_pair(bID2,Tree[rank[ID2]].disInf[j]));
//            m2.insert(make_pair(bID2,Tree[rank[ID2]].disInf[bID2]));
            m2.insert({bID2,QueryPeriphery_CH(ID2, bID2)});
        }

        for(int k=0;k<B1.size()-1;k++){
            bID1=B1[k];

            if(m1[bID1]>d)
                continue;

            for(int z=0;z<B2.size()-1;z++){
                bID2=B2[z];

                if(m2[bID2]>d)
                    continue;

                tempdis=m1[bID1]+QueryCore(bID1,bID2)+m2[bID2];
                if(tempdis<d){
                    d=tempdis;
                    d1=m1[bID1]; d2=m2[bID2];
                    b1=bID1; b2=bID2;
                }

            }
        }

    }

    return d;
}*/
//new function
int Graph::QueryPartiCore_CH(int ID1, int ID2){//ID1 partition, ID2 core
    int d=INF;

    int pid=CoreTag[ID1];
    vector<int> B=BoundVertex[pid];
    int bid;
    int dis1,dis2;

    set<int> todo;
    todo.clear();
    todo.insert(B.begin(), B.end());//get the partition vertex set
    todo.erase(B[B.size()-1]);

    map<int,int> results;
    QueryPeriphery_CH(ID1, todo, results);

    for(int i=0;i<B.size()-1;i++){
        bid=B[i];
        dis1=results[bid];
        dis2=QueryCore(bid,ID2);
        if(d>dis1+dis2)
            d=dis1+dis2;
    }

    return d;
}
//old function
/*int Graph::QueryPartiCore_CH(int ID1, int ID2){//ID1 partition, ID2 core
	int d=INF;

	int pid=CoreTag[ID1];
	vector<int> B=BoundVertex[pid];
	int bid;
	int dis1,dis2;

	for(int i=0;i<B.size()-1;i++){
		bid=B[i];
		dis1=QueryPeriphery_CH(ID1, bid);
		dis2=QueryCore(bid,ID2);
		if(d>dis1+dis2)
			d=dis1+dis2;
	}

	return d;
}*/

int Graph::Query_CH(int ID1, int ID2){
    int dis=INF;

    if(CoreTag[ID1]==-1 && CoreTag[ID2]==-1){//Case 1: both in core
        //cout<<"Core-Core"<<endl;
        dis=QueryCore(ID1, ID2);
    }else if(CoreTag[ID1]==-1 && CoreTag[ID2]!=-1){//Case 2: ID2 in partition, ID1 in core
        //cout<<"Core-Parti"<<endl;
        dis=QueryPartiCore_CH(ID2,ID1);
    }else if(CoreTag[ID1]!=-1 && CoreTag[ID2]==-1){//Case 2: ID1 in partition, ID2 in core
        //cout<<"Parti-Core"<<endl;
    	dis=QueryPartiCore_CH(ID1,ID2);
    }else if(CoreTag[ID1]!=-1 && CoreTag[ID2]!=-1){//both in partition
		if(CoreTag[ID1] != CoreTag[ID2]){//Case 3: in different peripheries
			dis= QueryPartiParti_CH(ID1,ID2);
		}else{//Case 4: in the same periphery
			dis= QuerySameParti_CH(ID1,ID2);
		}
    }
    return dis;
}

void Graph::EffiCheck_CH(string filename,int runtimes){
	ifstream IF(filename);
	if(!IF){
		cout<<"Cannot open Map "<<filename<<endl;
        exit(1);
	}
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
    cout<<"Run times: "<<runtimes<<endl;
	int s, t;
	std::chrono::high_resolution_clock::time_point t1, t2;
	std::chrono::duration<double> time_span;
	double runT;

	t1=std::chrono::high_resolution_clock::now();
	for(int i=0;i<runtimes;i++){
		Query_CH(ODpair[i].first,ODpair[i].second);
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"Average Query Time: "<<(double)runT*1000/runtimes<<" ms."<<endl;
}

void Graph::indexsizeCHH2H(){
	long long m=0,m1=0,m2=0,m3=0;

	//core index
	for(int k=0;k<Label.size();k++){
		m1+=Label[k].size()*2*sizeof(int);
	}


    for(int i=0;i<PruningPoint.size();i++){
        for(auto it=PruningPoint[i].begin();it!=PruningPoint[i].end();it++){
            m2+=(1+(*it).second.size())*sizeof(int);
        }
    }

	//Periphery index
   /* for(int i=0;i<Tree.size();i++){
        m3+=Tree[i].dis.size()*sizeof(int);
        m3+=Tree[i].disInf.size()*sizeof(int);
    }*/

    for(int i=0;i<NeighborCon.size();i++){
    	m3+=NeighborCon[i].size()*3*sizeof(int);
    }

    for(int i=0;i< SCconNodesMT.size();i++){
        for(auto it=SCconNodesMT[i].begin(); it!=SCconNodesMT[i].end(); it++){
            m3+=sizeof(int)+(*it).second.size()*sizeof(int);
        }
    }

    long long m4=0;
    //extended index
   /* if(ifExtension){
        assert(!IndexExt.empty());
        for(int i=0;i<nodenum;++i){
            if(!IndexExt[i].empty()){
                m4+=sizeof(int)+IndexExt[i].size()*2*sizeof(int);
            }
        }
    }*/


	//cout<<"Index size "<<(double)m1/1024/1024<<", Pruning point size "<<(double)m2/1024/1024<<endl;
	m=m1+m2+m3+m4;
	cout<<"Index size "<<(double)m/1024/1024<<" MB"<<endl;
}

void Graph::Create_Partitions() {
    cout<<"Creating partitions..."<<endl;
    //// Get tree
    vector<int> vecemp; //vecemp.clear();
    VidtoTNid.assign(nodenum,vecemp);

    //rank.assign(HighestOrder+2,0);
    rank.clear();
    rank.assign(nodenum,0);//the vector index of tree nodes, map from vertex to tree node
    int len=HighestOrder-1;
    heightMax=0;

    Node root;//virtual root node
    int x = vNodeOrder[len];
    if(NeighborCon[x].empty()){
        cout<<"There exist non-virtual root!"<<endl;
        root.uniqueVertex=x;
        len--;
    }else{
        root.uniqueVertex=-1;
    }
    root.height=1;
    Tree.push_back(root);

    int nn;
    for(;len>=0;len--){//check the vertices with order lower than HighestOrder
        x=vNodeOrder[len];
        if(existCore[x]){
            cout<<"Wrong: should be out of core"<<endl; exit(1);
        }
        Node nod;
        nod.vert=NeighborCon[x];
        nod.uniqueVertex=x;
        int pa=matchCore(x,NeighborCon[x]);

        //cout<<"pa "<<pa<<endl;

        Tree[pa].ch.push_back(Tree.size());
        nod.pa=pa;
        nod.height=Tree[pa].height+1;
        /// for update
        nod.hdepth=Tree[pa].height+1;
        for(int i=0;i<NeighborCon[x].size();i++){//for the neighbors which have higher order
            nn=NeighborCon[x][i].first;
            if(existCore[nn])//skip core vertex
                continue;
            VidtoTNid[nn].push_back(Tree.size());//record the child tree node rank who has direct super edge to nn
            if(Tree[rank[nn]].hdepth<Tree[pa].height+1)
                Tree[rank[nn]].hdepth=Tree[pa].height+1;
        }
        if(nod.height>heightMax)
            heightMax=nod.height;

        rank[x]=Tree.size();//the position of tree, higher-order vertex has lower rank
        if(pa==0){//root node of this tree
            nod.treeroot=rank[x];//tree root is itself
//            Tree[0].ch.push_back(Tree.size());//get the children of virtual root node
        }
        else{
            nod.treeroot=Tree[pa].treeroot;
        }
        Tree.push_back(nod);
    }

    //// Get partitions
    vector<int> cr;//the vertex rank of root's children
    cr.clear();
    int temp_ch=0;
    int temp_vert=0;
    int sum_singular=0;
    for(int k=0;k<Tree[0].ch.size();k++){//Tree[0] is the core
        int childrank=Tree[0].ch[k];
        if(Tree[childrank].ch.size()>0){///only the root which has children is regarded as a tree
            cr.push_back(childrank);
            temp_ch += Tree[childrank].ch.size();
            temp_vert += Tree[childrank].vert.size();
        }
        else{
            ++sum_singular;
//            cout<<"Single vertex in tree "<< childrank<<"!!!"<<endl;
        }
    }
//    cout<<"Tree[0].ch.size(): "<<Tree[0].ch.size()<<endl;
//    cout<<"Accumulated Tree[childrank].ch.size(): "<<temp_ch<<endl;
//    cout<<"Accumulated Tree[childrank].vert.size(): "<<temp_vert<<endl;
    partiNum=cr.size();

    cout<<"Periphery Number: "<<partiNum<<endl;
    cout<<"Nominal periphery root number: "<<sum_singular+partiNum<<endl;
    cout<<"Vertex number in core: "<<nodenum - HighestOrder + sum_singular <<endl;

    ///Get boundary vertex
    PartiUpdateExt.assign(partiNum, false);//initiation of partition update vector
    vector<int> vec;
    vec.clear();
    BoundVertex.assign(partiNum,vec);
    set<int> sset;
    sset.clear();
    BoundVertexSet.assign(partiNum,sset);
    BoundTag.assign(nodenum, make_pair(false,set<int>()));//the boundary vertex is consisted by interface vertex and root vertex
    map<int,int> PartiRoot;//partition root & partition ID
    PartiRoot.clear();
    unordered_set<int> BoundarySet;//the set of all boundary vertices
    BoundarySet.clear();

    int ID1,ID2;

    for(int PID=0;PID<partiNum;PID++){
        int childrank=cr[PID];
        for(int i=0;i<Tree[childrank].vert.size();i++){
            ID1 = Tree[childrank].vert[i].first;
            BoundVertex[PID].push_back(ID1);
            BoundVertexSet[PID].insert(ID1);
            BoundTag[ID1].first=true;//interface vertex
            BoundTag[ID1].second.insert(PID);

            BoundarySet.insert(ID1);
        }
        BoundVertex[PID].push_back(Tree[childrank].uniqueVertex);
//        BoundVertexSet[PID].insert(Tree[childrank].uniqueVertex);
//        BoundTag[Tree[childrank].uniqueVertex]=true;//root vertex

        PartiRoot.insert(make_pair(childrank,PID));//map from tree id to partition id
    }

    //get the adjacency list of interface vertex
    for(auto it=BoundarySet.begin();it!=BoundarySet.end();++it){
        ID1 = *it;
        assert(NeighborConCH[ID1].empty());
        for(auto it2=NeighborCon[ID1].begin();it2!=NeighborCon[ID1].end();++it2){
            ID2 = it2->first;
            if(BoundTag[ID2].first){//if the neighbor is boundary too
                NeighborConCH[ID1].emplace_back(ID2,it2->second);
            }
//            else{
//                cout<<"Not boundary: "<<ID1<<" "<<ID2<<endl;
//            }
        }
    }

    //if(PartiRoot.size()!=cr.size())
    //cout<<"something wrong with the boundary node"<<endl;

    CoreTag.assign(nodenum,-1);//-1 indicates core vertex or root vertex, i>=0 indicates non-core vertex (i.e., the inner-partition vertex) and which partition it belongs to
    int NodeID,RootNode,parentNode;
    int count=0;
    for(int len=HighestOrder-1;len>=0;len--){
        NodeID=vNodeOrder[len];
        RootNode=Tree[rank[NodeID]].treeroot;
        parentNode=Tree[rank[NodeID]].pa;
        if(parentNode!=0){//if it is not root vertex
            CoreTag[NodeID]=PartiRoot[RootNode];
        }
        else if(!Tree[rank[NodeID]].ch.empty()){//if the root vertex has children
            CoreTag[NodeID]=PartiRoot[RootNode];
        }else{
            count++;
        }
    }
    cout<<"Single vertex periphery: "<<sum_singular<<" "<<count<<endl;

    /// Get partition info: AdjaGraph (only for CT-DS) and AdjaCore
    AdjaCoreMap.clear();
    map<int,int> mii;
    mii.clear();
    AdjaCoreMap.assign(nodenum,mii);

    /// for AdjaCore
    for(int NodeID=0;NodeID<nodenum;NodeID++){
        if(CoreTag[NodeID]==-1){//core vertex
            for(int nei=0;nei<NeighborCon[NodeID].size();nei++) {//for each neighbor
                assert(CoreTag[NodeID]==-1);
                int neiID=NeighborCon[NodeID][nei].first;
                if(CoreTag[neiID] != -1){
                    cout<<"Wrong! The contracted neighbor "<<neiID<<" of "<<NodeID<<" is not core vertex!!!"<<endl;
                    exit(1);
                }
                AdjaCoreMap[NodeID][neiID]=NeighborCon[NodeID][nei].second.first;//the root vertex is regarded as core vertex
                AdjaCoreMap[neiID][NodeID]=NeighborCon[NodeID][nei].second.first;//the root vertex is regarded as core vertex
            }
        }
    }
    vector<pair<int,int>> vecpair;
    vecpair.clear();
    AdjaCore.assign(nodenum,vecpair);
    for(int i=0;i<AdjaCoreMap.size();i++){
        if(AdjaCoreMap[i].empty())
            continue;

        for(map<int,int>::iterator it=AdjaCoreMap[i].begin();it!=AdjaCoreMap[i].end();it++){
            AdjaCore[i].push_back(make_pair((*it).first, (*it).second));
        }
//        if(AdjaGraph[i].empty()){
//            AdjaGraph[i] = AdjaCore[i];
//        }else{
//            cout<<"Wrong!! AdjaGraph["<<i<<"] is not empty!"<<endl;
//            exit(1);
//        }
    }

    /// for AdjaGraph: current version of AdjaGraph is not truly correct
//    AdjaGraph.clear();
//    AdjaGraph.resize(nodenum);
    map<int,map<int,int>> mapvec;
    mapvec.clear();
    SuppPartiID.assign(nodenum, mapvec);
    map<int,pair<int,set<int>>> mapset;
    mapset.clear();
    SuppPartiIDReal.assign(nodenum, mapset);
    PartiVertex.assign(partiNum,unordered_set<int>());//only insert in-partition vertex in this case

    for(int id=0;id<nodenum;++id){
        if(CoreTag[id]!=-1){//for periphery vertex
            PartiVertex[CoreTag[id]].insert(id);
        }
    }

//    cout<<"Vertex number of each periphery:";
//    vector<int> PartiVertexNum(partiNum,0);
//    double aveVnum = 0.0;
//    for(int i=0;i<partiNum;++i){
//        PartiVertexNum[i] = PartiVertex[i].size();
//        cout<<" "<<PartiVertexNum[i];
//        aveVnum += PartiVertexNum[i];
//    }
//    cout<<" ; Average vertex number: "<<aveVnum / partiNum<<endl;

    for(int PID=0;PID<partiNum;PID++){
        int childrank=cr[PID];
        for(int i=0;i<Tree[childrank].vert.size();i++){
            ID1 = Tree[childrank].vert[i].first;
//            PartiVertex[PID].insert(ID1);//insert interface vertex

            for(int j=i+1;j<Tree[childrank].vert.size();++j){
                ID2 = Tree[childrank].vert[j].first;
                if(SuppPartiID[ID1].find(ID2)==SuppPartiID[ID1].end()){//if we cannot find ID2
                    SuppPartiID[ID1][ID2]=map<int,int>();
                    SuppPartiID[ID2][ID1]=map<int,int>();
                }
                SuppPartiID[ID1][ID2].insert({PID,INF});
                SuppPartiID[ID2][ID1].insert({PID,INF});
                for(auto it=SCconNodesMT[ID1][ID2].begin();it!=SCconNodesMT[ID1][ID2].end();++it){
                    assert(AdjaCoreMap[ID1].find(ID2)!=AdjaCoreMap[ID1].end());//must exist
                    assert(AdjaCoreMap[ID2].find(ID1)!=AdjaCoreMap[ID2].end());//must exist
                    if(CoreTag[it->first] == PID){
                        if(SuppPartiID[ID1][ID2][PID] > it->second){
                            SuppPartiID[ID1][ID2][PID] = it->second;
                            SuppPartiID[ID2][ID1][PID] = it->second;
                        }

                        if(SuppPartiIDReal[ID1].find(ID2) == SuppPartiIDReal[ID1].end()){
                            SuppPartiIDReal[ID1][ID2]=make_pair(AdjaCoreMap[ID1][ID2],set<int>());
                            SuppPartiIDReal[ID2][ID1]=make_pair(AdjaCoreMap[ID1][ID2],set<int>());
                        }
                        if(AdjaCoreMap[ID1][ID2] == it->second){
                            SuppPartiIDReal[ID1][ID2].second.insert(CoreTag[it->first]);
                            SuppPartiIDReal[ID2][ID1].second.insert(CoreTag[it->first]);
                        }
                    }else{
                        CoreTag[it->first];
                    }

                }
            }
        }
    }
    //clear useless variables
    existCore.clear();
}

void Graph::CorrectnessCheck_CH(int runtimes){
	srand (time(NULL));
	int s, t, d1, d2, d3;
//    runtimes = 1;
    cout<<"Correctness check ("<<runtimes<<" rounds) ..."<<endl;
	for(int i=0;i<runtimes;i++){
//        if(i%100==0)
//          cout<<i<<endl;
		s=rand()%nodenum;
		t=rand()%nodenum;
//        s=243736;t=146541;//parti-parti
//        s=142488;t=143850;//core-core
//        s=208312;t=210695;//same-parti
//        s=327243,t=752625;//parti-core
//        s=92495,t=401841;//core-parti
		d1=Dijkstra(s,t,Neighbor);
		d2=Query_CH(s,t);
//        d1=DijkstraCore(s,t);
//        d2=QueryCore(s,t);
//        cout<<s<<"("<<CoreTag[s]<<") "<<t<<"("<<CoreTag[t]<<") "<<d2<<" "<<d1<<endl;
		if(d1!=d2){
			cout<<"InCorrect! "<<i<<": "<<s<<" "<<t<<" "<<d2<<" "<<d1;
            cout<<" ; CoreTag: "<<CoreTag[s]<<" "<<CoreTag[t]<<endl;
            //QueryDebug(s,t);
		}
	}
}

vector<int> NodeOrders;

struct OrderComp{
    int x;
    int y;//order(x)<order(y)
    OrderComp(int _x, int _y){
        x=_x; y=_y;
    }
    bool operator< (const OrderComp& d) const{
        if(x==d.x && y==d.y){//avoid the redundant
            return false;
        }else{
            if(x!=d.x)
                return NodeOrders[x]<NodeOrders[d.x];
            if(y!=d.y)
                return NodeOrders[y]<NodeOrders[d.y];
        }
    }
};

void Graph::CHdecreaseBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch){
    //maintain the index caused by the weight change
    //NodeOrders.clear();
    NodeOrders.assign(NodeOrder.begin(),NodeOrder.end());
    set<OrderComp> OC;
    map<pair<int,int>,int> OCdis;//{(s,t),d} maintain the fresh distance and avoid search in the adjacent list
    //OC.clear(); OCdis.clear();

    int a,b,newW;//the weight of (a,b) decrease to newW
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first;
        b=wBatch[k].first.second;
        newW=wBatch[k].second.second;

        //modify the information in original graph
        /*for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
                Neighbor[b][i].second=newW;
                break;
            }
        }*/

        if(NodeOrder[a]<NodeOrder[b]){
            for(int i=0;i<NeighborConCH[a].size();i++){
                if(NeighborConCH[a][i].first==b){
                    if(NeighborConCH[a][i].second.first>newW){//if new edge weight is smaller than old edge weight
                        //cout<<OutNeighborCon[a][i].second.first<<"..........."<<newW<<endl;
                        NeighborConCH[a][i].second.first=newW;
                        NeighborConCH[a][i].second.second=1;

                        OCdis[make_pair(a,b)]=newW;
                        OC.insert(OrderComp(a,b));
                    }else if(NeighborConCH[a][i].second.first==newW)
                        NeighborConCH[a][i].second.second+=1;
                    break;
                }
            }
        }else{
            for(int i=0;i<NeighborConCH[b].size();i++){
                if(NeighborConCH[b][i].first==a){
                    if(NeighborConCH[b][i].second.first>newW){//if new edge weight is smaller than old edge weight
                        NeighborConCH[b][i].second.first=newW;
                        NeighborConCH[b][i].second.second=1;

                        OCdis[make_pair(b,a)]=newW;
                        OC.insert(OrderComp(b,a));
                    }else if(NeighborConCH[b][i].second.first==newW)
                        NeighborConCH[b][i].second.second+=1;
                    break;
                }
            }
        }
    }

//    cout<<"Flag 1"<<endl;
    while(!OC.empty()){
        int s=(*OC.begin()).x; int t=(*OC.begin()).y;
        int wt;
        OC.erase(OC.begin());
        wt=OCdis[make_pair(s,t)];
        /// deal with the shortcut update between interface vertices
        if(CoreTag[s]==-1 && CoreTag[t]==-1){//if both are core vertices
//            cout<<"Flag 2"<<endl;
            assert(BoundTag[s].first && BoundTag[t].first);
//            cout<<"Flag 3"<<endl;
            if(wt < AdjaCoreMap[s][t]){
//                cout<<"Flag 4"<<endl;
//                cout<<"core edge decrease update: "<<s<<" "<<t<<" "<<AdjaCoreMap[s][t]<<" "<<wt<<endl;
                AdjaCoreMap[s][t] = wt;
                AdjaCoreMap[t][s] = wt;
            }
        }


        map<int,int> InM2t; //InM2t.clear();
        vector<pair<int,int>> InMLower; //InMLower.clear();
        for(int i=0;i<NeighborConCH[s].size();i++){
            if(NodeOrder[NeighborConCH[s][i].first]>NodeOrder[t])//if the neighbor of has higher order than t
                InM2t.insert(make_pair(NeighborConCH[s][i].first,NeighborConCH[s][i].second.first));
            else if(NodeOrder[NeighborConCH[s][i].first]<NodeOrder[t])
                InMLower.push_back(make_pair(NeighborConCH[s][i].first,NeighborConCH[s][i].second.first));
        }
        int inID,inW,inWt;
        for(int i=0;i<NeighborConCH[t].size();i++){
            inID=NeighborConCH[t][i].first;
            if(InM2t.find(inID)!=InM2t.end()){//if found, i.e., inID is also t's neighbor
                inW=InM2t[inID];
                inWt=NeighborConCH[t][i].second.first;
                if(inWt>inW+wt){//if edge e(t, inID) needs to be updated
                    NeighborConCH[t][i].second.first=inW+wt;
                    NeighborConCH[t][i].second.second=1;
                    OCdis[make_pair(t,inID)]=inW+wt;
                    OrderComp oc={t,inID};
                    OC.insert(oc);
                }else if(inWt==inW+wt){
                    NeighborConCH[t][i].second.second+=1;
                }
            }
        }

        for(int i=0;i<InMLower.size();i++){
            inID=InMLower[i].first; inW=InMLower[i].second;
            for(int j=0;j<NeighborConCH[inID].size();j++){
                if(NeighborConCH[inID][j].first==t){//if t is also inID's neighbor
                    inWt=NeighborConCH[inID][j].second.first;
                    if(inWt>inW+wt){//if edge e(inID, t) needs to be updated
                        NeighborConCH[inID][j].second.first=inW+wt;
                        NeighborConCH[inID][j].second.second=1;

                        OCdis[make_pair(inID,t)]=inW+wt;
                        OrderComp oc={inID,t};
                        OC.insert(oc);
                    }else if(inWt==inW+wt)
                        NeighborConCH[inID][j].second.second+=1;
                    break;
                }
            }
        }
    }//finish change index
}

void Graph::CHincreaseBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch){
    //NodeOrders.clear();
    NodeOrders.assign(NodeOrder.begin(),NodeOrder.end());
    set<OrderComp> OC; //OC.clear();
    map<pair<int,int>,int> OCdis;//{(s,t),d} maintain the old distance before refreshed and avoid search in the adjacent list
    //OCdis.clear();

    for(int wb=0;wb<wBatch.size();wb++){
        int a=wBatch[wb].first.first;
        int b=wBatch[wb].first.second;
        int oldW=wBatch[wb].second.first;
        int newW=wBatch[wb].second.second;

        //modify the original graph information
        /*for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
                Neighbor[b][i].second=newW;
                break;
            }
        }*/


        if(NodeOrder[a]<NodeOrder[b]){
            for(int i=0;i<NeighborConCH[a].size();i++){
                if(NeighborConCH[a][i].first==b){
                    if(NeighborConCH[a][i].second.first==oldW){
                        NeighborConCH[a][i].second.second-=1;
                        if(NeighborConCH[a][i].second.second<1){//if the weight increase affect the shortcut weight
                            OrderComp oc={a,b};
                            OC.insert(oc);
                            OCdis[make_pair(a,b)]=oldW;
                        }
                    }
                    break;
                }
            }
        }else{
            for(int i=0;i<NeighborConCH[b].size();i++){
                if(NeighborConCH[b][i].first==a){
                    if(NeighborConCH[b][i].second.first==oldW){
                        NeighborConCH[b][i].second.second-=1;
                        if(NeighborConCH[b][i].second.second<1){//if the weight increase affect the shortcut weight
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
        /// deal with the shortcut update between interface vertices
        if(CoreTag[s]==-1 && CoreTag[t]==-1){//if both are core vertices
            assert(BoundTag[s].first && BoundTag[t].first);
            if(wt > AdjaCoreMap[s][t]){
                cout<<"core edge increase update: "<<s<<" "<<t<<" "<<AdjaCoreMap[s][t]<<" "<<wt<<endl;
                AdjaCoreMap[s][t] = wt;
                AdjaCoreMap[t][s] = wt;
            }
        }
        int inID,inW;
        map<int,int> HigherIn; vector<pair<int,int>> LowerIn;
        //HigherIn.clear(); LowerIn.clear();
        //the shortcuts infected by s-->t
        for(int i=0;i<NeighborConCH[s].size();i++){
            inID=NeighborConCH[s][i].first;
            inW=NeighborConCH[s][i].second.first;
            if(NodeOrder[inID]<NodeOrder[t]){
                LowerIn.push_back(make_pair(inID,inW));
            }else if(NodeOrder[inID]>NodeOrder[t]){
                HigherIn.insert(make_pair(inID,inW));
            }
        }
        for(int i=0;i<NeighborConCH[t].size();i++){
            inID=NeighborConCH[t][i].first;
            if(HigherIn.find(inID)!=HigherIn.end()){//if inID is the neighbor of both s and t
                inW=HigherIn[inID];
                if(NeighborConCH[t][i].second.first==wt+inW){
                    NeighborConCH[t][i].second.second-=1;
                    if(NeighborConCH[t][i].second.second<1){
                        OrderComp oc={t,inID};
                        OC.insert(oc);
                        OCdis[make_pair(t,inID)]=wt+inW;
                    }
                }
            }
        }
        for(int i=0;i<LowerIn.size();i++){
            inID=LowerIn[i].first; inW=LowerIn[i].second;
            for(int j=0;j<NeighborConCH[inID].size();j++){
                if(NeighborConCH[inID][j].first==t){//if inID is the neighbor of both s and t
                    if(NeighborConCH[inID][j].second.first==inW+wt){
                        NeighborConCH[inID][j].second.second-=1;
                        if(NeighborConCH[inID][j].second.second<1){
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
        for(int i=0;i<Neighbor[s].size();i++){
            if(Neighbor[s][i].first==t){
                wt=Neighbor[s][i].second;//the weight value in the original graph
                countwt=1;
                break;
            }
        }
        int ssw,wtt,wid;
        vector<pair<int,int>> Wnodes;
//        vector<int> Wnodes; //Wnodes.clear();
        Wnodes=SCconNodesMT[s][t];
        /*if(s<t){
            //Wnodes=SCconNodes[make_pair(s,t)]; //cout<<"wid num "<<Wnodes.size()<<endl;
            Wnodes=SCconNodesMT[s][t];
        }else{
            //Wnodes=SCconNodes[make_pair(t,s)];
            Wnodes=SCconNodesMT[s][t];
        }*/

        for(int i=0;i<Wnodes.size();i++){
            wid=Wnodes[i].first;
            for(int j=0;j<NeighborConCH[wid].size();j++){
                if(NeighborConCH[wid][j].first==s){
                    ssw=NeighborConCH[wid][j].second.first;
                }
                if(NeighborConCH[wid][j].first==t){
                    wtt=NeighborConCH[wid][j].second.first;
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
        for(int i=0;i<NeighborConCH[s].size();i++){
            if(NeighborConCH[s][i].first==t){
                NeighborConCH[s][i].second.first=wt;
                NeighborConCH[s][i].second.second=countwt;
                break;
            }
        }
    }
}

void Graph::Decrease_CH(int a, int b, int oldW, int newW){
    for(int i=0;i<Neighbor[a].size();i++){
        if(Neighbor[a][i].first==b){
            Neighbor[a][i].second=newW;
            break;
        }
    }
    for(int i=0;i<Neighbor[b].size();i++){
        if(Neighbor[b][i].first==a){
            Neighbor[b][i].second=newW;
            break;
        }
    }

    int pid = -1;
    //extUpdate = false;

    AdjaCoreMapOld = AdjaCoreMap;

    if((CoreTag[a]==-1 && BoundTag[a].first==0) || (CoreTag[b]==-1 && BoundTag[b].first==0)){// either endpoint is non-boundary core vertex
        //cout<<"******************change 1*******************"<<endl;
//        cout<<"Core-Core"<<endl;
//		DecreasePLL(a,b,oldW,newW,AdjaCore,Label);
        DecreasePSL(a,b,oldW,newW,AdjaCore,Label);
        //extUpdate = true;
    }else if(CoreTag[a]!=-1 || CoreTag[b]!=-1){//either endpoint is in-partition vertex
        //cout<<"******************change 2*******************"<<endl;
        if(CoreTag[a]!=-1)//if a is periphery vertex
            pid=CoreTag[a];
        else
            pid=CoreTag[b];

        //cout<<"<<<<<<<<<<"<<endl;
        //cout<<CoreTag[a]<<" "<<BoundTag[a]<<" "<<CoreTag[b]<<" "<<BoundTag[b]<<" pid "<<pid<<endl;
        //cout<<"decrease "<<a<<" "<<b<<" "<<oldW<<" "<<newW<<endl;
//        DecreaseH2HNew(a,b,newW,Neighbor,Tree,rank,heightMax,true);
//        DecreaseH2HNew(a,b,newW,Neighbor,Tree,rank,heightMax,false);
        vector<pair<pair<int,int>,pair<int,int>>> wBatch;
        wBatch.emplace_back(make_pair(make_pair(a,b), make_pair(oldW,newW)));
        CHdecreaseBat(wBatch);

        //cout<<">>>>>>>>>>"<<endl;

        /// check whether the update of periphery has affected core
        int ID1,ID2,dis,olddis;
        for(int i=0;i<BoundVertex[pid].size()-1;i++){
            ID1=BoundVertex[pid][i];

            //if(PartiVertexInverted[pid].find(ID1)!=PartiVertexInverted[pid].end()){//only if ID1 has neighbors in partition pid
            for(int j=i+1;j<BoundVertex[pid].size()-1;j++){
                ID2=BoundVertex[pid][j];

                //if(PartiVertexInverted[pid].find(ID2)!=PartiVertexInverted[pid].end()){
//                dis = QueryPeripheryTree(ID1,ID2,pid);///
                dis = AdjaCoreMap[ID1][ID2];
//                if(dis != temp_dis)
//                    cout<<dis<<" "<<temp_dis<<endl;
                olddis=AdjaCoreMapOld[ID1][ID2];
                if(olddis>dis){
//                   DecreasePLL(ID1,ID2,olddis,dis,AdjaCore,Label);
                    DecreasePSL(ID1,ID2,olddis,dis,AdjaCore,Label);
                    //extUpdate = true;
                }


            }

        }
    }else if(BoundTag[a].first==1 && BoundTag[b].first==1){//Both end points are interface vertex
        //cout<<"******************change 3*******************"<<endl;
        ///periphery update
//        cout<<"Both interface vertex!"<<endl;
        int olddis=AdjaCoreMap[a][b];
        if(olddis>newW){
            AdjaCoreMap[a][b]=newW;
            AdjaCoreMap[b][a]=newW;
//            DecreasePLL(a,b,olddis,newW,AdjaCore,Label);
            DecreasePSL(a,b,olddis,newW,AdjaCore,Label);
            //extUpdate = true;
        }

    }

}

//New version
void Graph::Increase_CH(int a, int b, int oldW, int newW){
    for(int i=0;i<Neighbor[a].size();i++){
        if(Neighbor[a][i].first==b){
//            if(oldW != Neighbor[a][i].second){
//                cout<<"Old edge weight is incorrect! "<<a<<" "<<b<<": "<<oldW<<" "<<Neighbor[a][i].second<<endl;
//                oldW = Neighbor[a][i].second;
//            }
            Neighbor[a][i].second=newW;
            break;
        }
    }
    for(int i=0;i<Neighbor[b].size();i++){
        if(Neighbor[b][i].first==a){
//            if(oldW != Neighbor[b][i].second){
//                cout<<"Old edge weight is incorrect! "<<b<<" "<<a<<": "<<oldW<<" "<<Neighbor[b][i].second<<endl;
//                oldW = Neighbor[b][i].second;
//            }
            Neighbor[b][i].second=newW;
            break;
        }
    }

    AdjaCoreMapOld = AdjaCoreMap;

    extUpdate = false;
    int pid = -1;
    if((CoreTag[a]==-1 && !BoundTag[a].first) || (CoreTag[b]==-1 && !BoundTag[b].first)){//edges in the core
        //cout<<"edge in core"<<endl;
//		IncreasePLL(a,b,oldW,newW,AdjaCore,Label,PruningPointNew,NoSupportedPair);
//        IncreasePSL(a,b,oldW,newW,AdjaCore,Label,PruningPointNew,NoSupportedPair);
        IncreasePSL(a,b,oldW,newW,AdjaCore,Label,PruningPoint);
        extUpdate = true;
    }else if(CoreTag[a]!=-1 || CoreTag[b]!=-1){//edges in one partition
        //cout<<"edge in partition"<<endl;

        if(CoreTag[a]!=-1)
            pid=CoreTag[a];
        else
            pid=CoreTag[b];

        int ID1,ID2,dis,olddis;

        //cout<<"<<<<<<<<<<"<<endl;
        //cout<<CoreTag[a]<<" "<<BoundTag[a]<<" "<<CoreTag[b]<<" "<<BoundTag[b]<<" pid "<<pid<<endl;
//        IncreaseH2HNew(a,b,oldW,newW,Neighbor,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,true);
//        IncreaseH2HNew(a,b,oldW,newW,Neighbor,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,false);
        vector<pair<pair<int,int>,pair<int,int>>> wBatch;
        wBatch.emplace_back(make_pair(make_pair(a,b), make_pair(oldW,newW)));
        CHincreaseBat(wBatch);
        //cout<<">>>>>>>>>>"<<endl;

//        ofstream OF("/Users/zhouxj/Documents/1-Research/Datasets/NY/NYCore.update");
        for(int i=0;i<BoundVertex[pid].size()-1;i++){
            ID1=BoundVertex[pid][i];

            //if(PartiVertexInverted[pid].find(ID1)!=PartiVertexInverted[pid].end()){//only if ID1 has neighbors in partition pid
            for(int j=i+1;j<BoundVertex[pid].size()-1;j++){
                ID2=BoundVertex[pid][j];

                dis = AdjaCoreMap[ID1][ID2];
                olddis=AdjaCoreMapOld[ID1][ID2];
                if(olddis < dis){
//                    cout<<ID1<<" "<<ID2<<" "<<olddis<<" "<<dis<<endl;

//                    OF<<ID1<<" "<<ID2<<" "<<olddis<<" "<<dis<<endl;
//                    IncreasePSL(ID1,ID2,olddis,dis,AdjaCore,Label,PruningPointNew,NoSupportedPair);
                    IncreasePSL(ID1,ID2,olddis,dis,AdjaCore,Label,PruningPoint);
                    extUpdate = true;

                }
            }
        }
//        OF.close();

    }else{//Both end points are boundary vertex
        //cout<<CoreTag[a]<<" "<<BoundTag[a]<<" "<<CoreTag[b]<<" "<<BoundTag[b]<<endl;
        //cout<<"edge between boundary vertex"<<endl;
//        cout<<"Both interface vertex!"<<endl;
        int olddis=AdjaCoreMap[a][b];
        assert(olddis <= oldW);
        if(olddis==oldW){//it indicates the update of e(a,b) may influence AdjaCoreMap[a][b]
            int newDis=INF;
            set<int> newSets;
            for(auto it=SuppPartiID[a][b].begin();it!=SuppPartiID[a][b].end();++it){
                if(it->second <= newW){
                    if(it->second < newDis){
                        newDis = it->second;
                        newSets.clear();
                        newSets.insert(it->first);
                    }else if(it->second == newDis){
                        newSets.insert(it->first);
                    }
                }
            }
            if(newDis == INF){//no supportive vertex can obtain super edge weight lower than newW
                AdjaCoreMap[a][b] = newW;
                AdjaCoreMap[b][a] = newW;
                SuppPartiIDReal[a][b].first = newW;
                SuppPartiIDReal[b][a].first = newW;
                SuppPartiIDReal[a][b].second.clear();
                SuppPartiIDReal[b][a].second.clear();
//                IncreasePSL(a,b,oldW,newW,AdjaCore,Label,PruningPointNew,NoSupportedPair);
                IncreasePSL(a,b,oldW,newW,AdjaCore,Label,PruningPoint);
                extUpdate = true;
            }else{
                AdjaCoreMap[a][b] = newDis;
                AdjaCoreMap[b][a] = newDis;
                SuppPartiIDReal[a][b].first = newDis;
                SuppPartiIDReal[b][a].first = newDis;
                SuppPartiIDReal[a][b].second = newSets;
                SuppPartiIDReal[b][a].second = newSets;
//                IncreasePSL(a,b,oldW,newDis,AdjaCore,Label,PruningPointNew,NoSupportedPair);
                IncreasePSL(a,b,oldW,newDis,AdjaCore,Label,PruningPoint);
                extUpdate = true;
            }

        }

    }


}

