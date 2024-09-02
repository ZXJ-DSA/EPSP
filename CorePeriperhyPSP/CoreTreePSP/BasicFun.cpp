/*
 * BasicFun.cpp
 *
 *  Created on: 13 Oct 2022
 *      Author: zhangmengxuan
 */
#include "head.h"

extern vector<int> NodeOrder_;//nodeID order

//// Graph RW
//function for reading graph
void Graph::ReadGraph(string filename){
	ifstream inGraph(filename);
	if(!inGraph){
		cout<<"Cannot open Map "<<filename<<endl;
        exit(1);
	}
    cout<<"Reading graph... "<<filename<<endl;
    Timer tt;
    tt.start();
	string line;
    getline(inGraph,line);
    vector<string> vs;
    boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
    node_num=stoi(vs[0]); edge_num=0;
    int tempENum=stoi(vs[1]);
    getline(inGraph,line);
	//graph g initialize
	Neighbor.assign(node_num, vector<pair<vertex,int>>());
	set<int> m; m.clear();
	vector<set<int>> E;
	E.assign(node_num,m);

	int ID1,ID2, weight;
	while(!line.empty()){
		vector<string> vs;
		boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
		ID1=stoi(vs[0]); ID2=stoi(vs[1]); weight=stoi(vs[2]);
//        weight=1;

		if(E[ID1].find(ID2)==E[ID1].end()){//if not found
			edge_num+=2;
			Neighbor[ID1].push_back(make_pair(ID2,weight));
			Neighbor[ID2].push_back(make_pair(ID1,weight));
			E[ID1].insert(ID2);
			E[ID2].insert(ID1);
		}
		if(inGraph.eof()) break;
		getline(inGraph,line);
	}
    inGraph.close();
    tt.stop();
	cout<<"Finish Reading! Vertex number: "<<node_num<<"; Edge number: "<<edge_num<<". Time: "<<tt.GetRuntime()<<" s."<< endl;
    assert(edge_num == tempENum);
}
//function of writing the whole graph into disk
void Graph::WriteGraph(string graphfile){
    cout<<"Writing graph..."<<endl;
    ofstream OF(graphfile);
    if(!OF){
        cout<<"Cannot open "<<graphfile<<endl;
        exit(1);
    }
    OF<<node_num<<" "<<edge_num<<endl;
    int tempE=0;
    for(int ID1=0;ID1<node_num;ID1++){
        for(int j=0;j<Neighbor[ID1].size();j++){
            int ID2=Neighbor[ID1][j].first;
            int wei=Neighbor[ID1][j].second;
            OF<<ID1<<" "<<ID2<<" "<<wei<<endl;
            tempE++;
        }
    }
    OF.close();
    cout<<"Done."<<endl;
}
//function of writing core graph into disk
void Graph::WriteCoreGraph(string graphfile){
    cout<<"Writing core graph..."<<endl;
    ofstream OF(graphfile);
    if(!OF){
        cout<<"Cannot open "<<graphfile<<endl;
        exit(1);
    }
//    OF<<node_num<<endl;
    int tempE=0;
    for(int ID1=0;ID1<node_num;ID1++){
        for(int j=0;j<AdjaCore[ID1].size();j++){
            int ID2=AdjaCore[ID1][j].first;
            int wei=AdjaCore[ID1][j].second;
            OF<<ID1<<" "<<ID2<<" "<<wei<<endl;
            tempE++;
        }
    }
    OF.close();
    ofstream OF2(graphfile+".order");
    if(!OF2){
        cout<<"Cannot open "<<graphfile+".order"<<endl;
        exit(1);
    }
//    OF2<<node_num<<" "<<tempE<<endl;
//    for(int ID1=0;ID1<node_num;ID1++){
//        OF2<<ID1<<" "<<NodeOrder[ID1]<<endl;
//    }
    int ID;
    OF2<<node_num<<endl;
    for(int i=node_num-1;i>=0;--i){
        ID=vNodeOrder[i];
        if(CoreTag[ID]==-1){
            OF2<<vNodeOrder[i]<<endl;
        }else{
            break;
        }
    }
    OF2.close();
    cout<<"Done."<<endl;
}
//function of reading core graph into disk
void Graph::ReadCoreGraph(string filename){
    ifstream inGraph(filename);
    if(!inGraph){
        cout<<"Cannot open Map "<<filename<<endl;
        exit(1);
    }

    string line;

    //graph g initialize


    getline(inGraph,line);
    vector<string> vs;
    boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
    node_num=stoi(vs[0]); edge_num=0;
    cout<<"Node number: "<<node_num<<endl;
    AdjaCore.assign(node_num, vector<pair<vertex,int>>());
    AdjaCoreMap.assign(node_num,map<int,int>());
    getline(inGraph,line);

    int ID1,ID2, weight;
    while(!line.empty()){
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
        ID1=stoi(vs[0]); ID2=stoi(vs[1]); weight=stoi(vs[2]);

        CoreVertex.insert(ID1);
        AdjaCore[ID1].push_back(make_pair(ID2,weight));
        AdjaCoreMap[ID1].insert({ID2,weight});
//        AdjaCore[ID1].push_back(make_pair(ID2,1));
//        AdjaCoreMap[ID1].insert({ID2,1});
        edge_num++;
        if(inGraph.eof())
            break;
        getline(inGraph,line);
    }

    cout<<"Finish Graph Reading! Nnum "<<CoreVertex.size()<<" Enum "<<edge_num<<endl;

    for(ID1=0;ID1<node_num;++ID1){
        for(auto it=AdjaCore[ID1].begin();it!=AdjaCore[ID1].end();++it){
            ID2=it->first; weight=it->second;
            if(AdjaCoreMap[ID2].find(ID1) == AdjaCoreMap[ID2].end()){//if not found
                cout<<"Wrong! "<<ID1<<" "<<ID2<<endl;
            }else if(weight != AdjaCoreMap[ID2][ID1])
            {
                cout<<"Wrong! "<<ID1<<" "<<ID2<<weight<<AdjaCoreMap[ID2][ID1]<<endl;
            }

        }
    }

    CoreTag.assign(node_num,0);
    for(auto it=CoreVertex.begin();it!=CoreVertex.end();++it){
        CoreTag[*it] = -1;
    }


    ifstream inGraph2(filename+".order");
    if(!inGraph2){
        cout<<"Cannot open Map "<<filename+".order"<<endl;
    }
    getline(inGraph2,line);
    vs.clear();
    boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
    int node_num=stoi(vs[0]);
//    int edge_num=stoi(vs[1]);
    assert(node_num == node_num);
//    assert(edge_num == edge_num);
    getline(inGraph2,line);
    NodeOrder.assign(node_num,-1);
    vNodeOrder.assign(node_num,-1);

    int order;
    while(!line.empty()){
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
        ID1=stoi(vs[0]); order=stoi(vs[1]);
        NodeOrder[ID1] = order;
        vNodeOrder[order] = ID1;
        if(inGraph2.eof())
            break;
        getline(inGraph2,line);
    }
    cout<<"Node finish ordering!"<<endl;
    NodeOrder_ = NodeOrder;

}
//function of vertex order reading
void Graph::WriteOrder(string filename){
    //Write order file to disk
    ofstream OF(filename);
    if(!OF){
        cout<<"Cannot open Map "<<filename<<endl;
        exit(1);
    }
    cout<<"Writing vertex order..."<<endl;
    OF<<node_num<<endl;
    for(int i = 0; i < NodeOrder.size(); i++){
        OF << i << "\t" << NodeOrder[i] << endl;//ID, order
    }

    OF.close();
    cout<<"Write done."<<endl;
}

void Graph::ReadOrder(string filename){
    ifstream inFile(filename, ios::in);
    if (!inFile) { // if not exist
        cout << "Fail to open file" << filename << endl;
        exit(1);
    }
    int nodeNum;
    string line;
    getline(inFile,line);
    vector<string> vs;
    boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
    nodeNum=stoi(vs[0]);

    if(nodeNum != node_num){
        cout<<"Wrong vertex number: "<<nodeNum<<" "<<node_num<<endl;
        exit(1);
    }
    NodeOrder.assign(nodeNum,-1);
    vNodeOrder.assign(nodeNum,-1);

    getline(inFile,line);

    int ID, order, num=0;
    while(!line.empty()){
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
        ID=stoi(vs[0]); order=stoi(vs[1]);

        NodeOrder[ID] = order;
        vNodeOrder[order] = ID;
        num++;
        if(inFile.eof())
            break;
        getline(inFile,line);
    }
    if(num!=nodeNum){
        cout<<"Inconsistent! "<<num<< " "<<nodeNum<<endl;
    }
}

//function for comparing the orders
void Graph::CompareOrder(string filename1, string filename2){
    cout<<"Comparing orders..."<<endl;
    ifstream IF1(filename1);
    if(!IF1){
        cout<<"Cannot open file "<<filename1<<endl;
        exit(1);
    }
    ifstream IF2(filename2);
    if(!IF2){
        cout<<"Cannot open file "<<filename2<<endl;
        exit(1);
    }
    vector<int> order1, order2;//Label1 is the ground truth
    order1.assign(node_num,-1);
    order2.assign(node_num,-1);
    string line;
    int ID,ord;

    //read label 1
    getline(IF1,line);
    vector<string> vs;
    boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
    int node_num=stoi(vs[0]);
    assert(node_num == node_num);
    getline(IF1,line);
    while(!line.empty()){
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
        if(vs.size()==2){
            ID=stoi(vs[0]), ord=stoi(vs[1]);
            order1[ID] = ord;
        }else{
            cout<<"Wrong syntax! vs.size(): "<<vs.size() <<" "<< line<<endl;
            exit(1);
        }


        if(IF1.eof())
            break;
        getline(IF1,line);
    }
    IF1.close();

    //read label 2
    getline(IF2,line);
    vs.clear();
    boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
    node_num=stoi(vs[0]);
    assert(node_num == node_num);
    getline(IF2,line);
    while(!line.empty()){
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
        if(vs.size()==2){
            ID=stoi(vs[0]), ord=stoi(vs[1]);
            order2[ID] = ord;
        }else{
            cout<<"Wrong syntax! vs.size(): "<<vs.size() <<" "<< line<<endl;
            exit(1);
        }
        if(IF2.eof())
            break;
        getline(IF2,line);
    }
    IF2.close();

    int ord1,ord2;

    for(int i=0;i<node_num;++i){
        if(order1[i]!=order2[i]){
            cout<<"Inconsistent! "<<i<<" "<<order1[i]<<" "<<order2[i]<<endl;
        }
    }
    cout<<"Done."<<endl;
}



//// Index RW
//function of writing core label to disk
void Graph::WriteCoreIndex(string file){
    cout<<"Writing core index into disk..."<<endl;
    Timer tt;
    tt.start();
    int tempE=0;
    int ID2,wei;
    ofstream OF(file+".core");
    if(!OF){
        cout<<"Cannot open "<<file+".core"<<endl;
        exit(1);
    }
    OF<<node_num<<endl;
    for(int ID1=0;ID1<node_num;ID1++){
        OF<<ID1<<" "<<Label[ID1].size();
        for(auto it=Label[ID1].begin();it!=Label[ID1].end();++it){
            ID2=it->first;
            wei=it->second;
            OF<<" "<<ID2<<" "<<wei;
            tempE++;
        }
        OF<<endl;
    }
    OF.close();
    /// Write pruning point
//    PPRV.write(file+".PPR");

    cout<<"Writing PPR..."<<endl;
    ofstream OF2(file+".prune");
    if(!OF2){
        cout<<"Cannot open "<<file+".prune"<<endl;
        exit(1);
    }
    int hub;
    for(int ID1=0;ID1<node_num;ID1++){
        if(!PruningPointSet[ID1].empty()){//Order
            OF2<<ID1<<" "<<PruningPointSet[ID1].size();
            for(auto it=PruningPointSet[ID1].begin();it!=PruningPointSet[ID1].end();++it){//Order
                hub=it->first;
                OF2<<" "<<hub<<" "<<it->second.size();

                for(auto it2=it->second.begin();it2!=it->second.end();++it2){
                    ID2 = *it2;//it2->ID;
                    OF2<<" "<<ID2;
                }
            }
            OF2<<endl;
        }

    }
    OF2.close();
//    cout<<"Done."<<endl;
    tt.stop();
    cout<<"Time for core index writing: "<<tt.GetRuntime()<<" s."<<endl;
}
//function of reading core label to disk
void Graph::ReadCoreIndex(string file) {
    cout<<"Reading core index..."<<endl;
    /// read label
    ifstream IF(file+".core");
    if(!IF){//if the label file does not exist, construct it
        cout<<"Cannot open file "<<file+".core"<<endl;
        Construct_core(algoCoreC);
        WriteCoreIndex(file);
//        exit(1);
    }
    else{
        Timer tt;
        tt.start();
        Label.assign(node_num, unordered_map<vertex,int>());

        PruningPointSet.clear();
        PruningPointSet.assign(node_num,unordered_map<vertex,unordered_set<vertex>>());
        PruningPointSet2.clear();
        PruningPointSet2.assign(node_num,unordered_map<vertex,vertex>());


        int nodeNum;
        int tempE=0;
        int ID1,ID2,weight,sz;
        string line;

        getline(IF,line);
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
        nodeNum = stoi(vs[0]);
        assert(nodeNum == node_num);
        getline(IF,line);

        while(!line.empty()){
            vector<string> vs;
            boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
            ID1=stoi(vs[0]); sz=stoi(vs[1]);
            int index=2;

            for(int i=0;i<sz;++i){
                ID2= stoi(vs[index++]);
                weight=stoi(vs[index++]);
                Label[ID1].insert({ID2,weight});
            }


            if(IF.eof())
                break;
            getline(IF,line);
        }

        IF.close();
        /// read pruning point
        ifstream IF2(file+".prune");
        if(!IF2){
            cout<<"Cannot open file "<<file+".prune"<<endl;
            exit(1);
        }

        getline(IF2,line);
        int hub,sz1,sz2;
        while(!line.empty()){
            vector<string> vs;
            boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
            ID1=stoi(vs[0]); sz1=stoi(vs[1]);
            int index=2;
            for(int i=0;i<sz1;++i){
                hub= stoi(vs[index++]); sz2=stoi(vs[index++]);
                for(int j=0;j<sz2;++j){
                    ID2=stoi(vs[index++]);
                    PruningPointSet[ID1][hub].insert(ID2);
                    PruningPointSet2[ID1][ID2] = hub;
                }

            }

            if(IF2.eof())
                break;
            getline(IF2,line);
        }

        IF2.close();
        tt.stop();
        cout<<"Time for index reading: "<<tt.GetRuntime()<<" s."<<endl;
    }

}
//function of writing core label to disk
void Graph::WriteTreeIndex(string file){
    cout<<"Writing tree index into disk..."<<endl;
    ofstream OF(file+".tree");
    if(!OF){
        cout<<"Cannot open "<<file+".tree"<<endl;
        exit(1);
    }

    OF<<Tree.size()<<endl;
    int tempE=0;
    int ID2,wei;
    for(int ID1=0;ID1<Tree.size();ID1++){
        OF<<Tree[ID1].uniqueVertex;
        OF<<" "<<Tree[ID1].treeroot;
        OF<<" "<<Tree[ID1].pa;
        OF<<" "<<Tree[ID1].height;
        OF<<" "<<Tree[ID1].hdepth;
        OF<<" "<<Tree[ID1].ch.size();
        for(int i=0;i<Tree[ID1].ch.size();++i){
            OF<<" "<<Tree[ID1].ch[i];
        }
        OF<<" "<<Tree[ID1].vert.size();//vert
        for(int i=0;i<Tree[ID1].vert.size();++i){
            OF<<" "<<Tree[ID1].vert[i].first<<" "<<Tree[ID1].vert[i].second.first<<" "<<Tree[ID1].vert[i].second.second<<" "<<Tree[ID1].pos[i];
        }
        OF<<" "<<Tree[ID1].vAncestor.size()-1;//dis
        for(int i=0;i<Tree[ID1].vAncestor.size()-1;++i){
            OF<<" "<<Tree[ID1].vAncestor[i]<<" "<<Tree[ID1].dis[i]<<" "<<Tree[ID1].cnt[i]<<" "<<Tree[ID1].FN[i];
        }
        OF<<" "<<Tree[ID1].disInf.size();//disInf
        for(auto it=Tree[ID1].disInf.begin();it!=Tree[ID1].disInf.end();++it){
            OF<<" "<<it->first<<" "<<it->second<<" "<<Tree[ID1].FNInf[it->first];
        }
        OF<<endl;
    }
    OF.close();
    /// Write SCconNodesMT
    ofstream OF2(file+".sc");
    if(!OF2){
        cout<<"Cannot open "<<file+".sc"<<endl;
        exit(1);
    }
    int cid;
//    OF2<<node_num<<endl;
    for(int ID1=0;ID1<node_num;ID1++){
        if(!SCconNodesMT[ID1].empty()){
            OF2<<ID1<<" "<<SCconNodesMT[ID1].size();
            for(auto it=SCconNodesMT[ID1].begin();it!=SCconNodesMT[ID1].end();++it){
                ID2=it->first;
                OF2<<" "<<ID2<<" "<<it->second.size();
                for(auto it2=it->second.begin();it2!=it->second.end();++it2){
                    cid = it2->first; wei = it2->second;
                    OF2<<" "<<cid<<" "<<wei;
                }
            }
        }

    }
    OF2.close();
    cout<<"Done."<<endl;
}
//function of reading core label to disk
void Graph::ReadTreeIndex(string file) {
    cout<<"Reading tree index..."<<endl;
    /// read label
    ifstream IF(file+".tree");
    if(!IF){//if the label file does not exist, construct it
        cout<<"Cannot open file "<<file+".tree"<<endl;
        exit(1);
    }
    Timer tt;
    tt.start();


    int tree_num;
    int tempE=0;
    int ID1,ID2,weight,sz;
    string line;

    getline(IF,line);
    vector<string> vs;
    boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
    tree_num = stoi(vs[0]);

    Tree.reserve(tree_num);

    for(int id=0;id<tree_num;++id){
        getline(IF,line);
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
        Tree[id].uniqueVertex=stoi(vs[0]);
        Tree[id].treeroot=stoi(vs[1]);
        Tree[id].pa=stoi(vs[2]);
        Tree[id].height=stoi(vs[3]);
        Tree[id].hdepth=stoi(vs[4]);
        int index=5;
        //ch
        int sz=stoi(vs[index++]);
        for(int i=0;i<sz;++i){
            Tree[id].ch.emplace_back(stoi(vs[index++]));
        }
        //vert
        int vertS=stoi(vs[index++]);
        for(int i=0;i<vertS;++i){
            pair<int,pair<int,int>> temp;
            temp.first=stoi(vs[index++]);
            temp.second.first=stoi(vs[index++]);
            temp.second.second=stoi(vs[index++]);
            Tree[id].vert.emplace_back(temp);
            Tree[id].pos.emplace_back(stoi(vs[index++]));
        }
        //dis
        int ancS=stoi(vs[index++]);
        Tree[id].pos.emplace_back(ancS);
        for(int i=0;i<ancS;++i){
            Tree[id].vAncestor.push_back(stoi(vs[index++]));
            Tree[id].dis.push_back(stoi(vs[index++]));
            Tree[id].cnt.push_back(stoi(vs[index++]));
            Tree[id].FN.push_back(stoi(vs[index++]));
        }
        Tree[id].vAncestor.emplace_back(Tree[id].uniqueVertex);
        Tree[id].dis.emplace_back(0);
        //disInf
        int intS=stoi(vs[index++]);

        for(int i=0;i<ancS;++i){
            ID2=stoi(vs[index++]);
            Tree[id].disInf.insert({ID2,stoi(vs[index++])});
            Tree[id].FNInf.insert({ID2,stoi(vs[index++])});
        }
        Tree[id].vAncestor.emplace_back(Tree[id].uniqueVertex);
        Tree[id].dis.emplace_back(0);
    }

    IF.close();
    /// read SuppPartiID
    ifstream IF2(file+".sc");
    if(!IF2){
        cout<<"Cannot open file "<<file+".sc"<<endl;
        exit(1);
    }

    SCconNodesMT.assign(node_num, map<int, vector<pair<int,int>>>());
    getline(IF2,line);
    int sz1,sz2,cid,wei;
    while(!line.empty()){
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
        ID1=stoi(vs[0]); sz1=stoi(vs[1]);
        int index=2;
        for(int i=0;i<sz1;++i){
            ID2= stoi(vs[index++]); sz2=stoi(vs[index++]);
            vector<pair<int,int>> tempV;
            for(int j=0;j<sz2;++j){
                cid=stoi(vs[index++]); wei= stoi(vs[index++]);
                tempV.emplace_back(cid,wei);
            }
            SCconNodesMT[ID1].insert({ID2,tempV});
        }

        if(IF2.eof())
            break;
        getline(IF2,line);
    }

    IF2.close();


    tt.stop();
    cout<<"Time for index reading: "<<tt.GetRuntime()<<" s."<<endl;
}
//function of writing core label to disk
void Graph::WriteTreeIndex2(string file){
    cout<<"Writing tree index into disk..."<<endl;
    ofstream OF(file+".tree");
    if(!OF){
        cout<<"Cannot open "<<file+".tree"<<endl;
        exit(1);
    }

    OF<<Tree.size()<<endl;
    int tempE=0;
    int ID2,wei;
    for(int ID1=0;ID1<Tree.size();ID1++){
        OF<<Tree[ID1].uniqueVertex;
        OF<<" "<<Tree[ID1].treeroot;
        OF<<" "<<Tree[ID1].pa;
        OF<<" "<<Tree[ID1].height;
        OF<<" "<<Tree[ID1].hdepth;
        OF<<" "<<Tree[ID1].ch.size();
        for(int i=0;i<Tree[ID1].ch.size();++i){
            OF<<" "<<Tree[ID1].ch[i];
        }
        OF<<" "<<Tree[ID1].vert.size();//vert
        for(int i=0;i<Tree[ID1].vert.size();++i){
            OF<<" "<<Tree[ID1].vert[i].first<<" "<<Tree[ID1].vert[i].second.first<<" "<<Tree[ID1].vert[i].second.second<<" "<<Tree[ID1].pos[i];
        }
        OF<<" "<<Tree[ID1].vAncestor.size()-1;//dis
        for(int i=0;i<Tree[ID1].vAncestor.size()-1;++i){
            OF<<" "<<Tree[ID1].vAncestor[i]<<" "<<Tree[ID1].dis[i]<<" "<<Tree[ID1].cnt[i]<<" "<<Tree[ID1].FN[i];
        }
        OF<<" "<<Tree[ID1].disInf.size();//disInf
        for(auto it=Tree[ID1].disInf.begin();it!=Tree[ID1].disInf.end();++it){
            OF<<" "<<it->first<<" "<<it->second<<" "<<Tree[ID1].FNInf[it->first];
        }
        OF<<endl;
    }
    OF.close();
    /// Write SCconNodesMT
    ofstream OF2(file+".sc");
    if(!OF2){
        cout<<"Cannot open "<<file+".sc"<<endl;
        exit(1);
    }
    int cid;
//    OF2<<node_num<<endl;
    for(int ID1=0;ID1<node_num;ID1++){
        if(!SCconNodesMT[ID1].empty()){
            OF2<<ID1<<" "<<SCconNodesMT[ID1].size();
            for(auto it=SCconNodesMT[ID1].begin();it!=SCconNodesMT[ID1].end();++it){
                ID2=it->first;
                OF2<<" "<<ID2<<" "<<it->second.size();
                for(auto it2=it->second.begin();it2!=it->second.end();++it2){
                    cid = it2->first; wei = it2->second;
                    OF2<<" "<<cid<<" "<<wei;
                }
            }
        }

    }
    OF2.close();
    cout<<"Done."<<endl;
}
//// Update RW
void Graph::ReadUpdate(string filename,vector<pair<pair<int,int>,int>>& TestData){
    TestData.clear();

    int num, ID1, ID2, oldw;
    ifstream IF(filename);
    if(!IF){
        cout<<"Cannot open file "<<filename<<endl;
        exit(1);
    }
    IF>>num;
    for(int i=0;i<num;i++){
        IF>>ID1>>ID2>>oldw;
        TestData.push_back(make_pair(make_pair(ID1, ID2), oldw));
    }
    IF.close();
}
void Graph::ReadUpdate2(string filename,vector<pair<pair<int,int>,pair<int,int>>>& TestData){
    TestData.clear();

    int num, ID1, ID2, oldW,newW;
    ifstream IF(filename);
    if(!IF){
        cout<<"Cannot open file "<<filename<<endl;
        exit(1);
    }
    string line;
    getline(IF,line);

    while(!line.empty()){
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
        ID1=stoi(vs[0]); ID2=stoi(vs[1]); oldW=stoi(vs[2]); newW=stoi(vs[3]);
        TestData.push_back(make_pair(make_pair(ID1, ID2), make_pair(oldW,newW)));
        if(IF.eof())
            break;
        getline(IF,line);
    }

    IF.close();
}

void Graph::ReadUpdate3(string filename,vector<pair<pair<int,int>,tuple<int,int,int>>>& TestData){
    TestData.clear();

    int num, ID1, ID2, oldW,newW1,newW2;
    ifstream IF(filename);
    if(!IF){
        cout<<"Cannot open file "<<filename<<endl;
        exit(1);
    }
    string line;
    getline(IF,line);
    vector<string> vs;
    boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
    num=stoi(vs[0]);
    getline(IF,line);

    while(!line.empty()){
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
        if(vs.size()==5){
            ID1=stoi(vs[0]); ID2=stoi(vs[1]); oldW=stoi(vs[2]); newW1=stoi(vs[3]); newW2=stoi(vs[4]);
            TestData.push_back(make_pair(make_pair(ID1, ID2), make_tuple(oldW,newW1, newW2)));
        }else{
            cout<<"Wrong input! vs.size: "<<vs.size()<<" "<<line<<endl;
        }

        if(IF.eof())
            break;
        getline(IF,line);
    }

    IF.close();
}

//// Dijkstra
//Dijkstra's algorithm
int Graph::Dijkstra(int ID1, int ID2,vector<vector<pair<vertex,int>>> &Neighbor){
	benchmark::heap<2, int, int> pqueue(node_num);
	pqueue.update(ID1,0);

	vector<bool> closed(node_num, false);
	vector<int> distance(node_num, INF);
	vector<int> prece(node_num, 0);
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

		for(auto it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
			NNodeID=(*it).first;
			NWeigh=(*it).second+topNodeDis;
			if(!closed[NNodeID]){
				if(distance[NNodeID]>NWeigh){
					distance[NNodeID]=NWeigh;
					pqueue.update(NNodeID, NWeigh);
					prece[NNodeID]=topNodeID;
				}
			}
		}
	}
    //retrieve path
//    RetrievePath(ID1, ID2, prece);

	return d;
}
//function of retrieving the shortest path
void Graph::RetrievePath(int ID1, int ID2, vector<int> & prece){

	//path retrieval
	vector<int> path;
	path.clear();
	path.push_back(ID2);
	int preID=prece[ID2];
	while(preID!=ID1){
		path.push_back(preID);
		preID=prece[preID];
	}
	path.push_back(ID1);

    pair<int,int> highestVertex(-1,0);//ID, order
	cout<<"path from "<<ID1<<" to "<<ID2<<": "<<endl;
	for(int i=path.size()-1;i>-1;i--){
		cout<<" "<<path[i]<<"("<<CoreTag[path[i]]<<","<<BoundTag[path[i]].first<<","<<NodeOrder[path[i]]<<") ";//<<endl;
        if(NodeOrder[path[i]] > highestVertex.second){
            highestVertex.second = NodeOrder[path[i]];
            highestVertex.first = path[i];
        }
        if(i>0){
            for(int j=0;j<Neighbor[path[i]].size();++j){
                if(Neighbor[path[i]][j].first == path[i-1]){
                    cout<<Neighbor[path[i]][j].second<<endl;
                    break;
                }
            }
        }
	}
	cout<<endl;
    cout<<"Highest-order vertex: "<<highestVertex.first<<" ("<<highestVertex.second<<")"<<endl;
}
//Dijkstra's search in core
int Graph::DijkstraCore(int ID1, int ID2){
    if(CoreTag[ID1] != -1 || CoreTag[ID2] != -1){
        cout<<"CoreTag wrong! "<<ID1<<"("<<CoreTag[ID1]<<") "<<ID2<<"("<<CoreTag[ID2]<<")"<<endl;
        return -1;
    }
    benchmark::heap<2, int, int> pqueue(node_num);
    pqueue.update(ID1,0);

    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);
    vector<int> prece(node_num, 0);
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

        for(auto it=AdjaCore[topNodeID].begin();it!=AdjaCore[topNodeID].end();it++){
            NNodeID=(*it).first;
            NWeigh=(*it).second+topNodeDis;
            if(!closed[NNodeID]){
                if(distance[NNodeID]>NWeigh){
                    distance[NNodeID]=NWeigh;
                    pqueue.update(NNodeID, NWeigh);
                    prece[NNodeID]=topNodeID;
                }
            }
        }
    }
    //retrieve path
//    RetrievePath(ID1, ID2, prece);

    return d;
}
//Dijkstra's search in core, debug version
pair<int,int> Graph::DijkstraCoreDebug(int ID1, int ID2){
    if(CoreTag[ID1] != -1 || CoreTag[ID2] != -1){
        cout<<"CoreTag wrong! "<<ID1<<"("<<CoreTag[ID1]<<") "<<ID2<<"("<<CoreTag[ID2]<<")"<<endl;
        return make_pair(-1,-1);
    }
    benchmark::heap<2, int, int> pqueue(node_num);
    pqueue.update(ID1,0);

    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);
    vector<int> prece(node_num, 0);
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

        for(auto it=AdjaCore[topNodeID].begin();it!=AdjaCore[topNodeID].end();it++){
            NNodeID=(*it).first;
            NWeigh=(*it).second+topNodeDis;
            if(!closed[NNodeID]){
                if(distance[NNodeID]>NWeigh){
                    distance[NNodeID]=NWeigh;
                    pqueue.update(NNodeID, NWeigh);
                    prece[NNodeID]=topNodeID;
                }
//                else if(distance[NNodeID]==NWeigh){
//                    if(NodeOrder[topNodeID] > NodeOrder[prece[NNodeID]]){
//                        prece[NNodeID]=topNodeID;
//                    }
//                }
            }
        }
    }

    //path retrieval
    vector<int> path;
    path.clear();
    path.push_back(ID2);
    int preID=prece[ID2];
    while(preID!=ID1){
        path.push_back(preID);
        preID=prece[preID];
    }
    path.push_back(ID1);

    int peakhub=-1;
    for(int i=path.size()-1;i>-1;i--){
        if(peakhub == -1){
            peakhub = path[i];
        }else{
            if(NodeOrder[path[i]] > NodeOrder[peakhub]){
                peakhub = path[i];
            }
        }

    }

    return make_pair(d,peakhub);
}

//// Other preprocessing
//function of generating query and update core label to disk
void Graph::CTQueryUpdateGen(int num) {
    cout<<"Writing CT queries into disk..."<<endl;
    string filename=sourcePath+"tmp/"+dataset;
    string QSameCore=filename+"_sameCore_CT" + to_string(bandWidth) + ".query";
    string QCoreTree=filename+"_coreTree_CT" + to_string(bandWidth) + ".query";
    string QTreeTree=filename+"_crossTree_CT" + to_string(bandWidth) + ".query";
    string QSameTree=filename+"_sameTree_CT" + to_string(bandWidth) + ".query";
    string USameCore=filename+"_sameCore_CT" + to_string(bandWidth) + ".update";
    string USameTree=filename+"_sameTree_CT" + to_string(bandWidth) + ".update";

    Timer tt;
    tt.start();
    int tempE=0;
    int ID1,ID2,wei;

    vector<pair<int,int>> vQSameCore;
    vector<pair<int,int>> vQCoreCore;
    vector<pair<int,int>> vQTreeTree;
    vector<pair<int,int>> vQSameTree;

    int baseCore = node_num-HighestOrder;

//    cout<<HighestOrder<<" "<<baseCore<<endl;
//    ID1 = vNodeOrder[HighestOrder-1];
//    ID2 = vNodeOrder[HighestOrder];
//    cout<<ID1<<" "<<CoreTag[ID1]<<" "<<NodeOrder[ID1]<<"; "<<ID2<<" "<<CoreTag[ID2]<<" "<<NodeOrder[ID2]<<endl;

    /// Same-core query
    ofstream outFile;
    outFile.open(QSameCore);
    if(!outFile){
        cout<<"Cannot open "<<QSameCore<<endl; exit(1);
    }
    tt.start();
    outFile<<num<<endl;
    for(int i=0;i<num;++i){
        ID1 = vNodeOrder[rand()%baseCore+HighestOrder];
        ID2 = vNodeOrder[rand()%baseCore+HighestOrder];
        while(ID1 == ID2){
            ID2 = vNodeOrder[rand()%baseCore+HighestOrder];
        }
        if(CoreTag[ID1]==-1 && CoreTag[ID2]==-1){
            outFile<<ID1<<" "<<ID2<<endl;
        }else{
            cout<<"Wrong Same-Core OD pair! "<<ID1<<"("<<CoreTag[ID1]<<"), "<<ID2<<"("<<CoreTag[ID2]<<")."<<endl; exit(1);
        }
    }
    outFile.close();
    tt.stop();
    cout<<"Same-core Query Generation Time: "<<tt.GetRuntime()<<" s."<<endl;

    map<int,vector<int>> tVertex;
    int pairs=0,tSize=tRoots.size();
    int tid;
    /// Same-tree query
    outFile.open(QSameTree);
    if(!outFile){
        cout<<"Cannot open "<<QSameTree<<endl; exit(1);
    }
    tt.start();
    outFile<<num<<endl;
    while (pairs < num) {
        tid = tRoots[rand() % tSize];
        if(tVertex.find(tid) == tVertex.end()){//if not found
            vector<int> tNodes;
            DFSTree(tNodes,tid);
            tVertex.insert({tid,tNodes});
        }
        vector<int> tNodes = tVertex[tid];
        ID1 = tNodes[rand() % tNodes.size()];
        ID2 = tNodes[rand() % tNodes.size()];
        while(ID1 == ID2){
            ID2 = tNodes[rand() % tNodes.size()];
        }
        if(CoreTag[ID1]==CoreTag[ID2] && CoreTag[ID1]!=-1){
            outFile << ID1 << " " << ID2 << endl;
        }else{
            cout<<"Wrong Same-Tree OD pair! "<<ID1<<"("<<CoreTag[ID1]<<"), "<<ID2<<"("<<CoreTag[ID2]<<")."<<endl; exit(1);
        }
        ++pairs;
    }
    outFile.close();
    tt.stop();
    cout<<"Same-tree Query Generation Time: "<<tt.GetRuntime()<<" s."<<endl;

    /// Core-Tree query
    outFile.open(QCoreTree);
    if(!outFile){
        cout<<"Cannot open "<<QCoreTree<<endl; exit(1);
    }
    tt.start();
    outFile<<num<<endl;
    for(int i=0;i<num;++i){
        ID1 = vNodeOrder[rand()%baseCore+HighestOrder];//core
        ID2 = vNodeOrder[rand()%HighestOrder];//tree
        if(CoreTag[ID1] ==-1 && CoreTag[ID2]!=-1){
            outFile << ID1 << " " << ID2 << endl;
        }else{
            cout<<"Wrong Same-Tree OD pair! "<<ID1<<"("<<CoreTag[ID1]<<"), "<<ID2<<"("<<CoreTag[ID2]<<")."<<endl; exit(1);
        }
    }
    outFile.close();
    tt.stop();
    cout<<"Core-tree Query Generation Time: "<<tt.GetRuntime()<<" s."<<endl;

    /// Tree-Tree query
    outFile.open(QTreeTree);
    if(!outFile){
        cout<<"Cannot open "<<QTreeTree<<endl; exit(1);
    }
    tt.start();
    outFile<<num<<endl;
    for(int i=0;i<num;++i){
        ID1 = vNodeOrder[rand()%HighestOrder];
        ID2 = vNodeOrder[rand()%HighestOrder];
        while(ID1 == ID2 || CoreTag[ID1]==CoreTag[ID2]){
            ID2 = vNodeOrder[rand()%HighestOrder];
        }
        if(CoreTag[ID1] != CoreTag[ID2] && CoreTag[ID1]!=-1 && CoreTag[ID2]!=-1){
            outFile << ID1 << " " << ID2 << endl;
        }else{
            cout<<"Wrong Same-Tree OD pair! "<<ID1<<"("<<CoreTag[ID1]<<"), "<<ID2<<"("<<CoreTag[ID2]<<")."<<endl; exit(1);
        }
    }
    outFile.close();
    tt.stop();
    cout<<"Tree-tree Query Generation Time: "<<tt.GetRuntime()<<" s."<<endl;

    int index_i=0;
    /// Same-Core Update
    outFile.open(USameCore);
    if(!outFile){
        cout<<"Cannot open "<<USameCore<<endl; exit(1);
    }
    tt.start();
    outFile<<num<<endl;
    int loop_i;
    pairs=0;
    while(pairs<num){
        ID1 = vNodeOrder[rand()%baseCore+HighestOrder];
        index_i = rand()%Neighbor[ID1].size();
        ID2 = Neighbor[ID1][index_i].first;
//        loop_i=0;
//        while(CoreTag[ID2]!=-1){
//            index_i = rand()%Neighbor[ID1].size();
//            ID2 = Neighbor[ID1][index_i].first;
//            ++loop_i;
//            if(loop_i>=5){
//                break;
//            }
//        }
        if(CoreTag[ID1]==-1 && CoreTag[ID2]==-1 && Neighbor[ID1][index_i].second>2){
            outFile<<ID1<<" "<<ID2<<" "<<Neighbor[ID1][index_i].second<< endl;
            ++pairs;
        }
//        else{
//            cout<<"Wrong Same-Core Update pair! "<<ID1<<"("<<CoreTag[ID1]<<"), "<<ID2<<"("<<CoreTag[ID2]<<")."<<endl; exit(1);
//        }
    }
    outFile.close();
    tt.stop();
    cout<<"Same-core Update Generation Time: "<<tt.GetRuntime()<<" s."<<endl;

    /// Same-Tree update
    outFile.open(USameTree);
    if(!outFile){
        cout<<"Cannot open "<<USameTree<<endl; exit(1);
    }
    tt.start();
    outFile<<num<<endl;
    pairs=0;
    while(pairs<num){
        ID1 = vNodeOrder[rand()%HighestOrder];
        index_i = rand()%Neighbor[ID1].size();
        ID2 = Neighbor[ID1][index_i].first;
//        while(CoreTag[ID2]==-1){
//            ID2 = Neighbor[ID1][rand()%Neighbor[ID1].size()].first;
//            index_i = rand()%Neighbor[ID1].size();
//        }
        if(CoreTag[ID1]!=-1 && Neighbor[ID1][index_i].second>2){
            outFile<<ID1<<" "<<ID2<<" "<<Neighbor[ID1][index_i].second<< endl;
            ++pairs;
        }
//        else{
//            cout<<"Wrong Same-Tree Update pair! "<<ID1<<"("<<CoreTag[ID1]<<"), "<<ID2<<"("<<CoreTag[ID2]<<")."<<endl; exit(1);
//        }
    }
    outFile.close();
    tt.stop();
    cout<<"Same-tree Update Generation Time: "<<tt.GetRuntime()<<" s."<<endl;


//    cout<<"Done."<<endl;
    tt.stop();
    cout<<"Time for core index writing: "<<tt.GetRuntime()<<" s."<<endl;
}
//function of generating update edges
void Graph::UpdateGene(int num, string filename){
    vector<pair<pair<int,int>, pair<int,int>>> UpdateData;

    set<pair<int,int>> Edges;
    vector<pair<pair<int,int>,int>> ENodeID;
    int ID1,ID2,wei;
    for(int i=0;i<Neighbor.size();i++){
        ID1=i;
        for(int j=0;j<Neighbor[i].size();j++){
            ID2=Neighbor[i][j].first;
            wei=Neighbor[i][j].second;
            if(ID1<ID2 && Edges.find(make_pair(ID1,ID2))==Edges.end()){
                Edges.insert(make_pair(ID1,ID2));
                ENodeID.push_back(make_pair(make_pair(ID1,ID2),wei));
            }
            else if(ID2<ID1 && Edges.find(make_pair(ID2,ID1))==Edges.end()){
                Edges.insert(make_pair(ID2,ID1));
                ENodeID.push_back(make_pair(make_pair(ID2,ID1),wei));
            }
        }
    }

    ofstream OF(filename);
    OF<<num<<endl;
    set<int> eid;
    for(int k=0;k<num;k++){
        int edgeid=rand()%ENodeID.size();
        if(eid.find(edgeid)==eid.end()){
            OF<<ENodeID[edgeid].first.first<<" "<<ENodeID[edgeid].first.second<<" "<<ENodeID[edgeid].second<<endl;
            eid.insert(edgeid);
        }else{
            k--;
        }
    }
    OF.close();
}
//function of generating OD pairs
void Graph::ODGene(int num, string filename){
    set<pair<int,int>> ODpair;
    vector<pair<int,int>> ODpairVec;

    srand (time(NULL));
    int s, t;
    for(int i=0;i<num;i++){
        s=rand()%node_num;
        t=rand()%node_num;
        if(ODpair.find(make_pair(s,t))==ODpair.end()){
            ODpairVec.push_back(make_pair(s,t));
            ODpair.insert(make_pair(s,t));
            ODpair.insert(make_pair(t,s));
        }else{
            i--;
        }
    }
    cout<<"generated OD pair number "<<ODpairVec.size()<<endl;

    ofstream OF(filename);
    OF<<ODpairVec.size()<<endl;
    for(int k=0;k<ODpairVec.size();k++){
        OF<<ODpairVec[k].first<<" "<<ODpairVec[k].second<<endl;
    }
    OF.close();
}
//function of connectivity checking
void Graph::StainingMethod(int ID){
    queue<int> Q;

    vector<bool> Stained;
    Stained.assign(node_num, false);

    Q.push(ID);
    Stained[ID]=true;
    int frontid, neiid;
    while(!Q.empty()){
        frontid=Q.front();
        Q.pop();
        for(int k=0;k<Neighbor[frontid].size();k++){
            neiid=Neighbor[frontid][k].first;
            if(!Stained[neiid]){
                Q.push(neiid);
                Stained[neiid]=true;
            }
        }
    }

    int stainNum=0;
    for(int i=0;i<node_num;i++){
        if(Stained[i])
            stainNum+=1;
    }
    //cout<<"Stained Number "<<stainNum<<endl;
    if(stainNum != node_num){
        cout<<"Incorrect!!! stain number: "<<stainNum<<" ; node number: "<<node_num<<endl;
    }

    vector<int> VertexInverted;
    VertexInverted.assign(node_num, -1);
    int j=0;
    for(int i=0;i<node_num;i++){
        if(Stained[i]){
            VertexInverted[i]=j;
            j+=1;
        }
    }
    //cout<<"Check j= "<<j<<", stainNum= "<<stainNum<<endl;

    int Orinode_num=node_num;
    node_num=stainNum;
    vector<vector<pair<vertex,int>>> Neighbor1=Neighbor;
    Neighbor.clear();
    Neighbor.assign(node_num, vector<pair<vertex,int>>());
    int InvertedID, nei, Invertednei, wei;
    for(int ID=0;ID<Orinode_num;ID++){
        if(VertexInverted[ID]!=-1){
            InvertedID=VertexInverted[ID];
            for(int k=0;k<Neighbor1[ID].size();k++){
                nei=Neighbor1[ID][k].first;
                wei=Neighbor1[ID][k].second;
                if(VertexInverted[nei]!=-1){
                    Invertednei=VertexInverted[nei];
                    Neighbor[InvertedID].push_back(make_pair(Invertednei,wei));
                }
            }
        }
    }

}
//function of checking the connectivity, set_A: the vertex set
vector<int> Graph::DFS_CC(vector<map<int,int>> & Edges, set<int> set_A, set<int> & set_LCC, int node_num) {
    /// DFS for connected component
    stack<int> stack_A;
//    set<int> set_A;//nodes waiting for visited
    unordered_set<int> set_B;//nodes visited for current component
    set_B.clear();
    int item_id,temp_id;
    vector<bool> flag_visited(node_num,false);
    bool flag_finish = false;
    int temp_num = 0;
    int component_i = 0;
    pair<unordered_set<int>,int> LCC;
    vector<int> CCs;//the vertex size of each connected component

//    for(int i=0;i<node_num;++i){
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
        if (set_B.size() > LCC.first.size()) {
            LCC.first.clear();
            LCC.first = set_B;
            LCC.second = temp_num;// /2
        }
        assert(!set_B.empty());
        CCs.push_back(set_B.size());
//        if(!set_B.empty() && set_B.size() < mcc.first.size()){
//            cout<<"Smaller connected component with vertex size "<<set_B.size()<<": ";
//            for(auto it=set_B.begin();it!=set_B.end();++it){
//                cout<<*it<<" ";
//            }
//            cout<<"; degree: ";
//            for(auto it=set_B.begin();it!=set_B.end();++it){
//                cout<<Edges[*it].size()<<" ";
//            }
//            cout<<endl;
//        }
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
        cout<<"This graph has only one connected component. ";
        cout<<"Nodes size of graph: "<< LCC.first.size() << " ; ";
        cout<<"Edges size of graph: "<< LCC.second << endl;
    }else{
        cout<<"!!! This graph has "<< component_i <<" connected component!"<<endl;
        cout<<"Nodes size of the largest connected component is: "<<LCC.first.size()<<endl;
        cout<<"Edges size of the largest connected component is: "<<LCC.second<<endl;
    }
    for(auto it=LCC.first.begin();it!=LCC.first.end();++it){
        set_LCC.insert(*it);
    }
    std::sort(CCs.begin(), CCs.end());
    return CCs;
//    return component_i;
}


