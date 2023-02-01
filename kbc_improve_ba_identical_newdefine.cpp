#include<iostream>
#include<vector>
#include<map>
#include<string>
#include<fstream>
#include<sstream>
#include<queue>
#include<stack>
#include<iomanip>
#include<set>
using namespace std;

class hetergraph {
public:
	int n_vertexclass;  
	map<char, int> vertex_class_index;  
	vector<int> vertex_num; 

	int n_edgeclass; 
	vector<vector<int>> edge_class; 
	vector<int> edge_num;  

	vector<map<int, map<int, int>>> edge_info;  

	void read_file();
	int getEdgeIndex(int x, int y); 
};


class multigraph {
public:
	int n_vertex;  //最新顶点数（加入copy的关节顶点）
	int n_vertex_org;  //原顶点数
	int m_edge;  //边数
	vector<vector<vector<int>>> edge_info;  
	vector<map<int, int> > adjlist_edge; 
;
	void getPmg(vector<char>* P);  //从异构图中获得Pmg
	void save_Mpl(vector<vector<int>>* Mpl);  //保存Mpl
	void deleteVertex(int v);
	void deleteEdge(int a, int b);
};

void CBC(int s);  //保存v后的路径实例信息
void CBC_(int s);  //直接在pred中保存边权；pred用vector  【最快】
void CBC_1(int s);  //直接在pred中保存边权；pred用map 
void CBC_ident0(int s);  //+ident，ident[s]=0
void CBC_ident1(int s);  //+ident，ident[s]>0
void CBC_reachident(int s); //+reach+ident
void save_bc();
void getbcfromfile();
void sort_bc();
void compare_result();
void Ident();
void todoType2Ident();
int judgeIdent(int i,int v,int type);
void mergeIdent2(int i, int v);
void BridgeEdgeDivide();
void b_tarjan(int index, int fa, int totalreach);
void ArticulationVertexDivide();
void a_tarjan(int index, int fa, int totalreach);
void artvertexcopy(int index, int j, int totalreach);
int getMaxComponent();
int getComponentSize(int i, vector<int>* visited);

string file = "data1";
hetergraph hg;
multigraph Pmg;
vector<double> bc;
map<int, vector<int>> IdentSet;  //顶点代表其他identical顶点的列表
map<int, map<int, int>> IdentSetEdgeindex;  //记录proxy顶点与其他代表的顶点之间边的索引
vector<int> ident;     //顶点代表其他identical顶点的个数（不包括自己）
vector<int> flag_addtoident;  //已经加入identical点集的顶点
vector<int> reach;
map<int, int> org;  //map[100]=1  表示顶点100在关节顶点划分前的原顶点是1
ofstream out_result(file + "/al5_time_result.txt");
int counta = 0, countb = 0;
int count_ident_edge = 0;
double pmid;

int main() {
	//一、根据异构图计算Pmg
	hg.read_file();

	//vector<char> P = { 'A','M','D','M','A' };
	vector<char> P = { 'A','P','V','P','A' };
	Pmg.getPmg(&P);  //得到P-multigraph
	//Pmg.show_mg();
	//二、已经保存了Pmg，从文件中读取
	/*Pmg.getPmg_read_file();*/
	/*Pmg.show_mg();*/

	/*统计最大连通分量大小*/
	int maxComponent = getMaxComponent();


	bc.resize(Pmg.n_vertex);
	reach.resize(Pmg.n_vertex,1);
	/*deta0.resize(Pmg.n_vertex);*/
	
	BridgeEdgeDivide();  //桥边分割
	int maxComponentb = getMaxComponent();

	ArticulationVertexDivide();  //关节顶点分割
	bc.resize(Pmg.n_vertex);  //新加入的copy关节顶点的bc值的存储空间
	int maxComponent1 = getMaxComponent();

	Ident();
	int maxComponent2 = getMaxComponent();


	double ps5, pt5, time5 = 0;
	double time5bfs = 0, time5back = 0;
	int count = 0;
	for (int i = 0; i < Pmg.n_vertex; i++) {
		if (flag_addtoident[i] == 0) {
			ps5 = clock();
			if (ident[i] == 0) {
				CBC_ident0(i);  //当源点i不是identical点集中的proxy顶点
			}
			else {
				CBC_ident1(i);  //当源点i是proxy顶点
			}
			pt5 = clock();
			time5 += (double)(pt5 - ps5) / CLOCKS_PER_SEC;
			time5bfs += (double)(pmid - ps5) / CLOCKS_PER_SEC;
			time5back += (double)(pt5 - pmid) / CLOCKS_PER_SEC;

		}
		count++;
	}

	//根据关节顶点的映射关系，计算关节顶点的bc值
	for (auto it : org) {
		bc[it.second] += bc[it.first];
	}

	/*for (int i = 0; i < Pmg.n_vertex; i++) {
		if (bc[i] > 0) {
			cout << "顶点" << i << "的bc值：" << bc[i] << endl;
		}
	}*/
	//compare_result();
	save_bc();
	out_result.close();
	return 0;
	
}

int getMaxComponent() {
	vector<int> visited(Pmg.n_vertex);
	int maxcsize = 0;
	for (int i = 0; i < Pmg.n_vertex; i++) {
		if (visited[i] == 0) {
			int csize = getComponentSize(i, &visited);
			maxcsize = max(maxcsize, csize);
		}
	}
	vector<int>().swap(visited);
	return maxcsize;
}
int getComponentSize(int i, vector<int>* visited) {
	int count = 0;
	queue<int> Q;
	Q.push(i);
	int totalcount = 1;
	(*visited)[i] = 1;
	while (!Q.empty()) {
		int v = Q.front();
		Q.pop();
		for (auto w : Pmg.adjlist_edge[v]) {
			if ((*visited)[w.first] == 0) {
				Q.push(w.first);
				(*visited)[w.first] = 1;
				totalcount += 1;
			}
		}
	}
	queue<int>().swap(Q);
	return totalcount;
}


void compare_result() {
	string filename = file + "/result.txt";  
	ifstream input(filename, ios::in);  //输入文件流对象input
	//判断文件是否存在
	if (!input) {
		cerr << "file error!" << endl;
		exit(1);
	}

	vector<double> bc0(Pmg.n_vertex_org);
	//读入所有顶点的bc值//保存每个顶点的bc值的文件
	string bc_info;
	while (getline(input, bc_info)) {
		istringstream sbc_info(bc_info);
		int v;
		double bc;
		sbc_info >> v >> bc;
		bc0[v] = bc;
	}

	//比较bc值
	vector<int> different_v;
	for (int i = 0; i < Pmg.n_vertex_org; i++) {
		if (abs(bc[i] - bc0[i]) > 0.000001) different_v.push_back(i);
	}
	cout << "bc值不同的顶点数：" << different_v.size() << endl;
	for (int v : different_v) {

		cout << v << ":" << bc[v] << " " << bc0[v] << endl;
	}
}

/*划分桥边、关节顶点*/
int countTotalReach(int i) {
	int totalreach = 0;
	queue<int> Q;
	vector<int> visit(Pmg.n_vertex);
	Q.push(i);  visit[i] = 1;  totalreach += reach[i];
	while (!Q.empty()) {
		int v = Q.front();
		Q.pop();
		for (auto w : Pmg.adjlist_edge[v]) {
			if (visit[w.first] == 0) {
				Q.push(w.first);
				visit[w.first] = 1;
				totalreach += reach[w.first];
			}
		}
	}
	queue<int>().swap(Q);
	vector<int>().swap(visit);
	return totalreach;
}

vector<int> dfn;
vector<int> low;
int stamp;
stack<int> stack_bi;

/*划分桥边*/
void BridgeEdgeDivide() {
	dfn.resize(Pmg.n_vertex);  
	low.resize(Pmg.n_vertex); 
	for (int i = 0; i < Pmg.n_vertex; i++) {
		stamp = 0;
		int totalreach;
		if (!dfn[i]) {
			totalreach = countTotalReach(i);  
			b_tarjan(i, -1, totalreach);  
		}
		stack<int>().swap(stack_bi);
	}
	dfn.clear(); low.clear();
}

void b_tarjan(int index, int fa, int totalreach) {
	int child = 0;
	stamp += 1;
	dfn[index] = stamp;
	low[index] = stamp;
	stack_bi.push(index);
	map<int, int> indexadj = Pmg.adjlist_edge[index];
	for (auto adjj : indexadj) {
		int j = adjj.first;
		if (j == fa) continue;
		if (!dfn[j]) {
			b_tarjan(j, index, totalreach);
			low[index] = min(low[index], low[j]);
			if (low[j] > dfn[index]) {
				countb++;
				int repj = 0, repindex = 0;
				while (stack_bi.top() != j) {  //计算j所在的分量中，顶点的reach之和
					int temp = stack_bi.top();
					stack_bi.pop();
					repj += reach[temp];
				}
				repj += reach[j]; stack_bi.pop();
				repindex = totalreach - repj; 
				//更新index、j的reach值
				reach[index] += repj;  reach[j] += repindex;

				//累加缺少的pair-dependency
				bc[index] += ((double)repindex - 1) * repj;
				bc[j] += ((double)repj - 1) * repindex;
				//修改邻接表，邻接矩阵
				Pmg.deleteEdge(index, j);
			}
		}
		else if (dfn[j] < dfn[index]) {
			low[index] = min(low[index], dfn[j]);
		}
	}
	map<int, int>().swap(indexadj);
}


/*划分关节顶点*/
void ArticulationVertexDivide() {
	dfn.resize(Pmg.n_vertex);
	low.resize(Pmg.n_vertex);
	for (int i = 0; i < Pmg.n_vertex_org; i++) {
		stamp = 0;
		int totalreach;
		if (!dfn[i]) {
			totalreach = countTotalReach(i);  
			a_tarjan(i, -1, totalreach);  
		}
		stack<int>().swap(stack_bi);
	}
	dfn.clear(); low.clear();
}

void a_tarjan(int index, int fa, int totalreach) {
	int child = 0;
	dfn[index] = ++stamp;
	low[index] = stamp;
	stack_bi.push(index);
	map<int, int> indexadj = Pmg.adjlist_edge[index];  
	for (auto iteradjj = indexadj.begin(); iteradjj != indexadj.end(); iteradjj++) {
		int j = iteradjj->first;
		if (j == fa) continue;
		if (!dfn[j]) {
			child++;
			a_tarjan(j, index, totalreach);
			low[index] = min(low[index], low[j]);
			if (low[j] >= dfn[index] && fa >= 0) {
				counta++;
				artvertexcopy(index, j, totalreach);  
			}
			if (low[j] >= dfn[index] && fa < 0) {
				auto nowj = iteradjj;
				nowj++;
				int isbi = 0;
				for (nowj; nowj != indexadj.end(); nowj++) {  
					int nb = nowj->first;
					if (!dfn[nb]) {  
						counta++;
						isbi = 1;
						break;
					}
				}
				if (isbi == 1) {
					artvertexcopy(index, j, totalreach);
				}
			}
		}
		else if (dfn[j] < dfn[index]) {
			low[index] = min(low[index], dfn[j]);
		}
	}
	map<int, int>().swap(indexadj);
}


/*tarjan算法找到关节顶点后，复制关节顶点并更新邻接矩阵、邻接表、reach值*/
void artvertexcopy(int index, int j, int totalreach) {
	int copyindex = Pmg.n_vertex++;  
	cout << "复制关节点" << copyindex << endl;
	org[copyindex] = index;//保存copyindex的原点是index
	Pmg.adjlist_edge.resize(Pmg.n_vertex);
	int repcopyindex = 0, repindex = 0;
	while (stack_bi.top() != j) {  //将index-j这个连通分量中的顶点出栈
		int temp = stack_bi.top();
		stack_bi.pop();
		repcopyindex += reach[temp];
		auto itertemp = Pmg.adjlist_edge[index].find(temp);
		if (itertemp != Pmg.adjlist_edge[index].end()) {  

			int e_index = itertemp->second; 

			//在Pmg中删除边index-temp
			Pmg.deleteEdge(index, temp);
			Pmg.adjlist_edge[copyindex].insert(map<int, int>::value_type(temp, e_index));
			Pmg.adjlist_edge[temp].insert(map<int, int>::value_type(copyindex, e_index));

		}
	}
	stack_bi.pop();  repcopyindex += reach[j];
	int e_index1 = Pmg.adjlist_edge[index][j];

	Pmg.deleteEdge(index, j);
	Pmg.adjlist_edge[copyindex].insert(map<int, int>::value_type(j, e_index1));
	Pmg.adjlist_edge[j].insert(map<int, int>::value_type(copyindex, e_index1));

	repindex = totalreach - repcopyindex;  
	reach[index] += repcopyindex;  
	reach.push_back(repindex);  
}





void Ident() {
	ident.resize(Pmg.n_vertex);
	flag_addtoident.resize(Pmg.n_vertex);

	todoType2Ident();
	//统计identical点集和顶点
	int countset = 0, countvertex = 0;
	for (auto idset : IdentSet) {
		countset++;
		for (int i : idset.second) {
			countvertex++;
		}
	}

}

void todoType2Ident() {
	for (int i = 0; i < Pmg.n_vertex; i++) {  //顶点i=0-n_vertex，寻找其identical点集
		if (flag_addtoident[i] == 0) {
			vector<int> iwithidentset = { i };
			auto iteriadj = Pmg.adjlist_edge[i].begin();
			for (iteriadj; iteriadj != Pmg.adjlist_edge[i].end(); iteriadj++) {
				int v = iteriadj->first;
				if (v > i && flag_addtoident[v] == 0) {
					int flag_isident = judgeIdent(i, v, 2);
					if (flag_isident) {
						iwithidentset.push_back(v);
					}
				}
			}
			if (iwithidentset.size()>1) {
				int proxy = iwithidentset[0];
				for (int x = 1; x < iwithidentset.size(); x++) {
					flag_addtoident[iwithidentset[x]] = 1;
					mergeIdent2(proxy, iwithidentset[x]);
				}
			}
			vector<int>().swap(iwithidentset);
		}
	}
}

void mergeIdent2(int i, int v) {
	//当考虑reach值时
	int detav = reach[v] - 1;
	int reachiident = reach[i];
	if (ident[i] > 0) {
		for (int ii : IdentSet[i])
			reachiident += reach[ii];
	}
	bc[v] += (double)detav * reachiident;

	int detai = reach[i] - 1;
	bc[i] += (double)detai * reach[v];
	if (ident[i] > 0) {
		for (int ii : IdentSet[i]) {
			int detaii = reach[ii] - 1;
			bc[ii] += (double)detaii * reach[v];
		}
	}

	//记录i与v之间边的索引
	int eindexi_v = Pmg.adjlist_edge[i][v];
	if (ident[i] == 0) {  
		map<int, int> tempvindex;
		tempvindex.insert(map<int, int>::value_type(v, eindexi_v));
		IdentSetEdgeindex.insert(map<int,map<int,int>>::value_type(i,tempvindex));
		tempvindex.clear();
	}
	else {
		IdentSetEdgeindex[i].insert(map<int, int>::value_type(v, eindexi_v));
	}
	IdentSet[i].push_back(v);
	ident[i]++;
	count_ident_edge += Pmg.adjlist_edge[v].size();
	Pmg.deleteVertex(v);
}

int judgeIdent(int i, int v, int type) {
	if (Pmg.adjlist_edge[i].size() != Pmg.adjlist_edge[v].size()) return 0;
    
	auto iteriadj = Pmg.adjlist_edge[i].begin();
	auto itervadj = Pmg.adjlist_edge[v].begin();
	for (iteriadj, itervadj; iteriadj != Pmg.adjlist_edge[i].end()&&itervadj != Pmg.adjlist_edge[v].end(); iteriadj++, itervadj++) {
		if (iteriadj->first == v) { //从i的邻居中跳过v
			iteriadj++; 
		}  
		if (itervadj->first == i) { //从v的邻居中跳过i
			itervadj++;
		} 
		if (iteriadj != Pmg.adjlist_edge[i].end()&&itervadj != Pmg.adjlist_edge[v].end()) {
			if (iteriadj->first != itervadj->first) return 0;  
			//两个邻居点相同,开始判断邻居边的信息是否相同
			int i_neibedgeindex = iteriadj->second;
			int v_neibedgeindex = itervadj->second;
			set<vector<int>> i_neib_einfo{ Pmg.edge_info[i_neibedgeindex].begin(),Pmg.edge_info[i_neibedgeindex].end() };
			set<vector<int>> v_neib_einfo{ Pmg.edge_info[v_neibedgeindex].begin(),Pmg.edge_info[v_neibedgeindex].end() };

			if (i_neib_einfo.size() != v_neib_einfo.size()) return 0;  //如果路径实例类别数不同，则不是
			set<vector<int>> einfo_result;
			set_intersection(i_neib_einfo.begin(), i_neib_einfo.end(), v_neib_einfo.begin(), v_neib_einfo.end(), inserter(einfo_result, einfo_result.begin()));
			if (einfo_result.size() != i_neib_einfo.size()) return 0;
		}
		else {
			break;
		}
	}
	return 1;
} 



//void CBC_(int s) {  //pred用vector
//	//initialize variable
//	/*cout << "源点" << s << endl;*/
//	queue<int> Q;
//	stack<int> S;
//	/*vector<map<int,double>> pred(Pmg.n_vertex);*/
//	vector<vector<vector<double>>> pred(Pmg.n_vertex);
//	vector<double> dist(Pmg.n_vertex, DBL_MAX);
//	vector<double> pathnum(Pmg.n_vertex);
//	dist[s] = 0;
//	pathnum[s] = 1;
//	Q.push(s);
//	//bfs
//	while (!Q.empty()) {
//		int v = Q.front();
//		Q.pop();
//		S.push(v);
//		map<int, int> pathinstance;
//		for (auto w : Pmg.adjlist_edge[v]) {  //访问了v的所有邻居w[邻接顶点，邻接边]
//			if (dist[w[0]] == DBL_MAX) {
//				dist[w[0]] = dist[v] + 1;
//				Q.push(w[0]);
//			}
//			if (dist[w[0]] == dist[v] + 1) {  //通过此if条件，可以得到v的所有下方的邻居w
//				/*pathnum[w[0]] += pathnum[v]*Pmg.Mp[v][w[0]];*/
//				//统计v下方的路径实例信息
//				vector<vector<int>> einfo = Pmg.edge_info[w[1]];
//				for (auto pinfo : einfo) {
//					int pclass = pinfo[0];
//					int pnum = pinfo[1];
//					pathinstance[pclass] += pnum;
//				}
//			}
//		}
//		//遍历完v的邻居点后，才能统计完整的v后路径实例的信息，获得v与其每个后置点w之间边的权重
//		for (auto w : Pmg.adjlist_edge[v]) {
//			double weight = 0;
//			if (dist[w[0]] == dist[v] + 1) {  //计算v与w[0]之间的边权
//				vector<vector<int>> einfo = Pmg.edge_info[w[1]];  //w[1]是v与w[0]之间的边
//				for (auto pinfo : einfo) {  //根据v与w[0]之间边的路径实例信息，更新pathnum[w[0]]
//					int pclass = pinfo[0];
//					int pnum = pinfo[1];
//					weight += ((double)pnum / pathinstance[pclass]);  //累加v-w[0]之间的边权
//				}
//				pathnum[w[0]] += pathnum[v] * weight;
//				pred[w[0]].push_back({ (double)v,weight });
//			}
//		}
//	}
//	/*for (int i = 0; i < pathnum.size(); i++) {
//		cout << "到点" << i << "最短路径数：" << pathnum[i] << endl;
//	}*/
//	//back propagation
//	vector<double> deta(Pmg.n_vertex);
//	while (!S.empty()) {
//		int w = S.top();
//		S.pop();
//		for (vector<double> v : pred[w]) {
//			deta[(int)v[0]] += (pathnum[(int)v[0]] * v[1] / pathnum[w]) * (1 + deta[w]);
//		}
//		if (w != s) {
//			bc[w] += deta[w];
//		}
//	}
//	/*for (int i = 0; i < deta.size(); i++) {
//		cout << "对点" << i << "源依赖：" << deta[i] << endl;
//	}*/
//}
void CBC_ident0(int s) {  //pred用vector +ident后  s不是proxy
	//initialize variable
	queue<int> Q;
	stack<int> S;
	/*vector<map<int,double>> pred(Pmg.n_vertex);*/
	vector<vector<vector<double>>> pred(Pmg.n_vertex);
	vector<double> dist(Pmg.n_vertex, DBL_MAX);
	vector<double> pathnum(Pmg.n_vertex);
	dist[s] = 0;
	pathnum[s] = 1;
	Q.push(s);
	//bfs
	while (!Q.empty()) {
		int v = Q.front();
		Q.pop();
		S.push(v);
		map<int, int> pathinstance;
		for (auto w : Pmg.adjlist_edge[v]) {  //访问了v的所有邻居w[邻接顶点，邻接边]
			if (dist[w.first] == DBL_MAX) {
				dist[w.first] = dist[v] + 1;
				Q.push(w.first);
			}
			if (dist[w.first] == dist[v] + 1) {  //通过此if条件，可以得到v的所有下方的邻居w
				/*pathnum[w[0]] += pathnum[v]*Pmg.Mp[v][w[0]];*/
				//统计v下方的路径实例信息
				vector<vector<int>> einfo = Pmg.edge_info[w.second];
				for (auto pinfo : einfo) {
					int pclass = pinfo[0];
					int pnum = pinfo[1];
					pathinstance[pclass] += pnum*(1+ident[w.first]);  //v--w[0]与v--w[0]代表的顶点的路径实例数
				}
			}
		}
		//遍历完v的邻居点后，才能统计完整的v后路径实例的信息，获得v与其每个后置点w之间边的权重
		for (auto w : Pmg.adjlist_edge[v]) {
			double weight = 0;
			if (dist[w.first] == dist[v] + 1) {  //计算v与w[0]之间的边权
				vector<vector<int>> einfo = Pmg.edge_info[w.second];  //w[1]是v与w[0]之间的边
				for (auto pinfo : einfo) {  //根据v与w[0]之间边的路径实例信息，更新pathnum[w[0]]
					int pclass = pinfo[0];
					int pnum = pinfo[1];
					weight += ((double)pnum /pathinstance[pclass]);  //累加v-w[0]之间的边权
				}
				weight = weight * einfo.size();
				pathnum[w.first] += pathnum[v] * weight*((double)1+ident[v]);
				pred[w.first].push_back({ (double)v,weight });
				
			}
		}
	}
	pmid = clock();
	//back propagation  s不是proxy顶点
	vector<double> deta(Pmg.n_vertex);
	for (int i = 0; i < Pmg.n_vertex; i++) deta[i] = (double)reach[i] - 1;
	while (!S.empty()) {
		int w = S.top();
		S.pop();
		for (vector<double> v : pred[w]) {
			if (v[0] != s) {
				double deta_add = (pathnum[(int)v[0]] * v[1] / pathnum[w]) * (1 + deta[w]);
				if (ident[w] > 0) {
					for (int wident : IdentSet[w]) {
						deta_add += (pathnum[(int)v[0]] * v[1] / pathnum[w]) * (1 + deta[wident]);
					}
				}
				deta[(int)v[0]] += deta_add;
				if (ident[v[0]] > 0) {
					for (int vident : IdentSet[v[0]]) {
						deta[vident] += deta_add;
					}
				}
			}
		}
		if (w != s) {
			bc[w] += deta[w]*reach[s];
			if (ident[w] > 0) {
				for (int wident : IdentSet[w]) {
					bc[wident] += deta[wident]*reach[s];
				}
			}
		}
	}
	//cout << s<<"对点16源依赖：" << deta[16] << endl;
	/*deta0[s] = deta[0];*/
	/*for (int i = 0; i < deta.size(); i++) {
		cout << "对点" << i << "源依赖：" << deta[i] << endl;
	}*/
}
void CBC_ident1(int s) {  //pred用vector +ident后  s是proxy顶点  根据新定义，以s代表的每个顶为源点进行bfs得到的结果都相同
	//initialize variable
	queue<int> Q;
	stack<int> S;
	/*vector<map<int,double>> pred(Pmg.n_vertex);*/
	vector<vector<vector<double>>> pred(Pmg.n_vertex);
	vector<double> dist(Pmg.n_vertex, DBL_MAX);
	vector<double> pathnum(Pmg.n_vertex);
	dist[s] = 0;
	pathnum[s] = 1;
	Q.push(s);
	//bfs
	while (!Q.empty()) {
		int v = Q.front();
		Q.pop();
		S.push(v);
		map<int, int> pathinstance;
		for (auto w : Pmg.adjlist_edge[v]) {  //访问了v的所有邻居w[邻接顶点，邻接边]
			if (dist[w.first] == DBL_MAX) {
				dist[w.first] = dist[v] + 1;
				Q.push(w.first);
			}
			if (dist[w.first] == dist[v] + 1) {  //通过此if条件，可以得到v的所有下方的邻居w
				/*pathnum[w[0]] += pathnum[v]*Pmg.Mp[v][w[0]];*/
				//统计v下方的路径实例信息
				vector<vector<int>> einfo = Pmg.edge_info[w.second];
				for (auto pinfo : einfo) {
					int pclass = pinfo[0];
					int pnum = pinfo[1];
					pathinstance[pclass] += pnum * (1 + ident[w.first]);  //v--w[0]与v--w[0]代表的顶点的路径实例数
				}
			}
		}
		if (v == s) {  //当s为ident点集中的proxy点时，第一层bfs，除了访问邻居点w外，其实还有其identical点（也是s的邻居点）也在第二层顶点中，因为要统计完整的s下方路径实例的信息，所以也要加上s——sident的路径实例信息
			for (int sident : IdentSet[s]) {
				//int edge_index_s = Pmg.Mp_edge[s][sident];
				//!!!!!!!!int edge_index_s = Pmg.adjlist_edge[s][sident]; //??????????????????????????
				int edge_index_s = IdentSetEdgeindex[s][sident];
				vector<vector<int>> einfos = Pmg.edge_info[edge_index_s];
				for (auto pinfo : einfos) {
					int pclass = pinfo[0];
					int pnum = pinfo[1];
					pathinstance[pclass] += pnum * (1 + ident[sident]);
				}
			}
		}
		//遍历完v的邻居点后，才能统计完整的v后路径实例的信息，获得v与其每个后置点w之间边的权重
		for (auto w : Pmg.adjlist_edge[v]) {
			double weight = 0;
			if (dist[w.first] == dist[v] + 1) {  //计算v与w[0]之间的边权
				vector<vector<int>> einfo = Pmg.edge_info[w.second];  //w[1]是v与w[0]之间的边
				for (auto pinfo : einfo) {  //根据v与w[0]之间边的路径实例信息，更新pathnum[w[0]]
					int pclass = pinfo[0];
					int pnum = pinfo[1];
					weight += ((double)pnum / pathinstance[pclass]);  //累加v-w[0]之间的边权
				}
				weight = weight * einfo.size();
				pathnum[w.first] += pathnum[v] * weight * (1 + (double)ident[v]);
				pred[w.first].push_back({ (double)v,weight });
			}
		}
	}
	/*for (int i = 0; i < pathnum.size(); i++) {
		cout << "到点" << i << "最短路径数：" << pathnum[i] << endl;
	}*/
	pmid = clock();
	//back propagation  s是proxy顶点
	int total_reach = reach[s];
	for (int sident : IdentSet[s]) {
		total_reach += reach[sident];
	}
	vector<double> deta(Pmg.n_vertex);
	for (int i = 0; i < Pmg.n_vertex; i++) {
		deta[i] = (double)reach[i] - 1;
	}
	while (!S.empty()) {
		int w = S.top();
		S.pop();
		for (vector<double> v : pred[w]) {
			if (v[0] != s) {
				double deta_add = (pathnum[(int)v[0]] * v[1] / pathnum[w]) * (1 + deta[w]);
				if (ident[w] > 0) {
					for (int wident : IdentSet[w]) {
						deta_add += (pathnum[(int)v[0]] * v[1] / pathnum[w]) * (1 + deta[wident]);
					}
				}
				deta[(int)v[0]] += deta_add;
				if (deta[v[0]] > 0) {
					for (int vident : IdentSet[v[0]]) {
						deta[vident] += deta_add;
					}
				}
			}
		}
		if (w != s) {
			bc[w] += deta[w] * total_reach;
			if (ident[w] > 0) {
				for (int wident : IdentSet[w]) {
					bc[wident] += deta[wident] * total_reach;
				}
			}
		}
	}
	/*deta0[s] = deta[0];
	for (int sident:IdentSet[s]) {
		deta0[sident] = deta[0];
	}*/
	/*cout << "对点0源依赖：" << deta[0] << endl;*/

	/*for (int i = 0; i < deta.size(); i++) {
		cout << "对点" << i << "源依赖：" << deta[i] << endl;
	}*/
}

bool cmp(const pair<int, double>& a, const pair<int, double>& b) {
	return a.second > b.second;
}

void sort_bc() {
	vector<pair<int, double>> bc_;
	for (int i = 0; i < bc.size(); i++) {
		bc_.push_back({i,bc[i]});
	}
	sort(bc_.begin(), bc_.end(), cmp);
	cout << "top10:" << endl;
	for (int i = 0; i < 10; i++) {
		cout <<"顶点"<< bc_[i].first << "的bc值:" << bc_[i].second << endl;
	}
}

void getbcfromfile() {
	string filename = file + "/result.txt";  
	ifstream input(filename, ios::in);  //输入文件流对象input
	//判断文件是否存在
	if (!input) {
		cerr << "file error!" << endl;
		exit(1);
	}

	//读入所有顶点的bc值//保存每个顶点的bc值的文件
	string bc_info;
	while (getline(input, bc_info)) {
		istringstream sbc_info(bc_info);
		int v;
		double bcvalue;
		sbc_info >> v >> bcvalue;
		bc.push_back(bcvalue);
	}
}


void save_bc() {
	ofstream out(file + "/cbc_result.txt");
	for (int i = 0; i < Pmg.n_vertex; i++) {
		out << i << " ";
		out<< fixed << setprecision(8) << bc[i] << endl;
	}
	cout << "由cbc计算出的bc值已保存至" << file << "/cbc_result.txt";
}

/*读文件，获取异构图信息*/
void hetergraph::read_file() {  
	string filename = file + "/base.txt";  
	ifstream input(filename, ios::in);  //输入文件流对象input
	//判断文件是否存在
	if (!input) {
		cerr << "file error!" << endl;
		exit(1);
	}

	//第一行，顶点类别数
	string vinfo;
	getline(input, vinfo);
	istringstream svinfo(vinfo);
	svinfo >> n_vertexclass;

	for (int i = 0; i < n_vertexclass; i++) {
		string vcinfo;
		getline(input, vcinfo);
		istringstream svcinfo(vcinfo);
		char vclass;
		int vcnum;
		svcinfo >> vclass >> vcnum;
		vertex_class_index[vclass] = i;  
		vertex_num.push_back(vcnum);   
	}

	//边的类别数
	string einfo;
	getline(input, einfo);
	istringstream seinfo(einfo);
	seinfo >> n_edgeclass;

	//边的类别的信息
	for (int i = 0; i < n_edgeclass; i++) {  //第i类边
		string eeinfo;
		getline(input, eeinfo);
		istringstream seeinfo(eeinfo);
		int vx, vy, edgenum;
		seeinfo >> vx >> vy >> edgenum;
		edge_class.push_back({ vx,vy });  
		edge_num.push_back(edgenum);

		//获得保存第i类边信息的文件名并读入文件
		stringstream si;
		si << i;
		string edge_filename = file + "/edge/" + si.str() + ".txt";
		ifstream edge_input(edge_filename, ios::in);  //输入文件流对象input
		//判断文件是否存在
		if (!edge_input) {
			cerr << "file error!" << endl;
			exit(1);
		}
		string edgeinformation;
		getline(edge_input, edgeinformation);
		int nx = vertex_num[vx], ny = vertex_num[vy];
		map<int, map<int, int>> edge_matrix;
		string edgedata;
		while (getline(edge_input, edgedata)) {
			istringstream sedgedata(edgedata);
			int x, y;
			sedgedata >> x >> y;
			auto iterx = edge_matrix.find(x);  //x是否已经加入
			if (iterx == edge_matrix.end()) {  //尚未加入
				map<int, int> xadjlist;
				xadjlist.insert(map<int, int>::value_type(y, 1));  
				edge_matrix.insert(map<int, map<int, int>>::value_type(x, xadjlist)); 
				xadjlist.clear();
			}
			else {  //x已经加入
				iterx->second.insert(map<int, int>::value_type(y, 1));
			}
		}
		edge_info.push_back(edge_matrix);  
		edge_matrix.clear();
		edge_input.close();
	}
	input.close();
}

/*根据边的两个顶点的类别索引号，查找边的类别索引号*/
int hetergraph::getEdgeIndex(int x, int y) {
	int index = -1;
	for (auto it : edge_class) {
		index++;
		if (x == it[0] && y == it[1])
			break;
	}
	return index;
}

/*从异构图中，生成P-multigraph*/
void multigraph::getPmg(vector<char>* P) {  
	vector<char>::iterator it = P->begin();  
	/*Pmg顶点数*/
	n_vertex = hg.vertex_num[hg.vertex_class_index[*it]];
	n_vertex_org = n_vertex;
	adjlist_edge.resize(n_vertex);  
	//adjlist.resize(n_vertex);

	vector<int> eindex;
	for (int i = 0; i < P->size() / 2 ; i++) {  
		int x = hg.vertex_class_index[*it++];
		int y = hg.vertex_class_index[*it];
		int index = hg.getEdgeIndex(x, y);  
		if (index == -1) {
			cerr << "edge error!" << endl;
			exit(1);
		}
		eindex.push_back(index);
	}

	/*计算Mpl*/
	int p = 0;
	map<int, map<int, int>> matrixpl = hg.edge_info[eindex[p]];
	hg.edge_info[eindex[p]].clear();  //此时，第一个eindex的邻接矩阵已经保存在matrixpl中，信息不再需要
	while (++p < eindex.size()) {  //将matrix与下一个邻接矩阵edge_info[eindex[p]]相乘
		int curr_eindex = eindex[p];  //当前计算的边的索引,即当前用哪个边的邻接矩阵与MPL相乘
		cout << "计算" << eindex[p - 1] << '*' << curr_eindex << "的邻接矩阵" << endl;
		double ps11, pt11, time11;
		ps11 = clock();
		map<int, map<int, int>> tmatrix;  //临时记录计算结果
		auto iterik = matrixpl.begin();
		for (iterik; iterik != matrixpl.end(); iterik++) {  //对每个有邻接序列的i
			int i = iterik->first;
			auto iterk = iterik->second.begin();  //iterk:点i的邻接点k序列
			for (iterk; iterk != iterik->second.end(); iterk++) {
				int k = iterk->first;
				int currik = iterk->second;  //(y,1) 存在即>0 表示i与k邻接，之后找与k邻接的j
				auto iterkj = hg.edge_info[curr_eindex].find(k);
				if (iterkj != hg.edge_info[curr_eindex].end()) {  //k*j的邻接表中存在点k的邻接序列
					auto iterj = iterkj->second.begin();  //iterj：点k的邻接点j序列
					for (iterj; iterj != iterkj->second.end(); iterj++) {
						int j = iterj->first;
						int currkj = iterj->second;
						//向结果中，插入i的邻接点j，邻接多重边数为currik*currkj
						auto iterri = tmatrix.find(i);
						if (iterri == tmatrix.end()) {  //若还未更新i，则更新其第一个邻接点k
							map<int, int> iadjlistj;
							iadjlistj.insert(map<int, int>::value_type(j, currik * currkj));
							tmatrix.insert(map<int, map<int, int>>::value_type(i, iadjlistj));
							iadjlistj.clear();
						}
						else {  //若i已更新，判断j是否已插入到i的邻接表中
							auto iterrj = iterri->second.find(j);
							if (iterrj == iterri->second.end()) {  //若j未作为i的邻居被插入到i的邻接表中
								iterri->second.insert(map<int, int>::value_type(j, currik * currkj));
							}
							else {
								int oldij = iterrj->second;
								int nowij = oldij + currik * currkj;
								iterrj->second = nowij;  //更新ij位置的值
							}
						}
					}
				}
			}
		}
		matrixpl = tmatrix;
		tmatrix.clear();
		hg.edge_info[curr_eindex].clear();
		pt11 = clock();
		time11 = (double)(pt11 - ps11) / CLOCKS_PER_SEC;
		cout << " 用时：" << time11 << endl;
		out_result << "计算" << eindex[p - 1] << "*" << curr_eindex << "用时:" << time11 << endl;
	}
	hg.edge_info.clear();
	eindex.clear();

	//save_Mpl(&matrixpl);  //保存Mpl
    //得到Mpl邻接表的转置
	double ps12, pt12, time12;
	ps12 = clock();
	map<int, map<int, int>> matrixplt;
	for (auto it : matrixpl) {
		int x = it.first;
		auto ity = it.second.begin();
		for (ity; ity != it.second.end(); ity++) {
			int y = ity->first;
			int valuexy = ity->second;
			//y与x邻接，值为valuexy
			auto itfindy = matrixplt.find(y);
			if (itfindy == matrixplt.end()) {  //y未插入
				map<int, int> yadjlistx;
				yadjlistx.insert(map<int, int>::value_type(x, valuexy));
				matrixplt.insert(map<int, map<int, int>>::value_type(y, yadjlistx));
				yadjlistx.clear();
			}
			else {
				itfindy->second.insert(map<int, int>::value_type(x, valuexy));
			}
		}
	}
	cout << "转置计算完毕" << endl;
	double pt14 = clock();
	double time14 = (double)(pt14 - ps12) / CLOCKS_PER_SEC;
	out_result << "计算mpl的转置用时：" << time14 << endl;


	/*计算出Mpl后，计算Mp？？,得到edge_info、adjlist_edge、Mp_edge??*/
	cout << "计算Mp邻接表" << endl;
	int edge_index = 0;
	for (int i = 0; i < n_vertex; i++) {  //计算第i个点的邻接序列
	/*cout << "计算到第" << i;*/
		auto iterik = matrixpl.find(i);  //找第i个点的k邻接序列
		if (iterik != matrixpl.end()) {  //找到了
			map<int, vector<vector<int>>> ij_einfo;  //记录i的所有邻居j，与之间边的信息
			auto iterk = iterik->second.begin();
			for (iterk; iterk != iterik->second.end(); iterk++) {  //对于i的每个邻接点k
				int k = iterk->first;
				int currik = iterk->second;
				auto iterkj = matrixplt.find(k);//找k的邻接点j序列
				if (iterkj != matrixplt.end()) {
					auto iterj = iterkj->second.begin();
					for (iterj; iterj != iterkj->second.end(); iterj++) {
						int j = iterj->first;
						int currkj = iterj->second;
						if (j > i) {
							int ij_kclasspathnum = currik * currkj;  //第ij点之间，k类路径实例数
							auto iterjeinfo = ij_einfo.find(j);  //是否已经更新过i与邻居j的边的信息
							if (iterjeinfo == ij_einfo.end()) {  //没更新过i与邻居j的边的信息
								vector<vector<int>> temp_ij_einfo;
								temp_ij_einfo.push_back({ k,ij_kclasspathnum });
								ij_einfo.insert(map<int, vector<vector<int>>>::value_type(j, temp_ij_einfo));
								temp_ij_einfo.clear();
							}
							else {  //已经更新过i与邻居j的边的信息
								iterjeinfo->second.push_back({ k,ij_kclasspathnum });
							}
						}
					}
				}
			}
			for (auto jj : ij_einfo) {
				int j = jj.first;
				adjlist_edge[i].insert(map<int, int>::value_type(j, edge_index));
				adjlist_edge[j].insert(map<int, int>::value_type(i, edge_index));
				edge_info.push_back(jj.second);
				edge_index++;
			}
			map<int, vector<vector<int>>>().swap(ij_einfo);
		}
	}
	m_edge = edge_index;
	matrixpl.clear(); matrixplt.clear();
	pt12 = clock();
	time12 = (double)(pt12 - ps12) / CLOCKS_PER_SEC;
	cout << "由mpl得到mp用时：" << time12 << endl;
	out_result << "由mpl得到mp用时：" << time12 << endl;

	/*save_Pmg();*/
}
void multigraph::save_Mpl(vector<vector<int>>* Mpl) {  
	ofstream out(file + "/mpl.txt");
	int row = Mpl->size();
	int col = (*Mpl)[0].size();
	out << row << " " << col << endl;
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			out << (*Mpl)[i][j] << " ";
		}
		out << endl;
	}
	out.close();
}

void multigraph::deleteVertex(int a) {

	//对a的每个邻接点
	for (auto an : adjlist_edge[a]) {
		//修改邻接矩阵
		/*Mp[a][an] = 0; Mp[an][a] = 0;*/
		//在an的邻接表中删除a
		map<int, int>::iterator it = adjlist_edge[an.first].begin();
		for (it; it != adjlist_edge[an.first].end(); it++) {
			if (it->first == a) {
				adjlist_edge[an.first].erase(it);
				break;
			}
		}
	}
	//清空a的邻接表（此时a为P-multigraph中一个独立的点，但通过index遍历仍然能遍历到a）
	map<int, int>().swap(adjlist_edge[a]);
}

/*在P-multigraph中删除边*/
void multigraph::deleteEdge(int a, int b) {
	/*m_edge--;*/
	map<int, int>::iterator it = adjlist_edge[a].begin();
	for (it; it != adjlist_edge[a].end(); it++) {
		if (it->first == b) {
			adjlist_edge[a].erase(it);
			break;
		}
	}

	map<int, int>::iterator itt = adjlist_edge[b].begin();
	for (itt; itt != adjlist_edge[b].end(); itt++) {
		if (itt->first == a) {
			adjlist_edge[b].erase(itt);
			break;
		}
	}
}
