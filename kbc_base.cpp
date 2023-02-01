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
#include<math.h>
using namespace std;

class hetergraph {
public:
	int n_vertexclass;  
	map<char, int> vertex_class_index;  
	vector<int> vertex_num;  

	int n_edgeclass; 
	vector<vector<int>> edge_class; 
	vector<int> edge_num;  //每类边数

	vector<map<int, map<int, int>>> edge_info;  

	void read_file();
	void show_hg();
	int getEdgeIndex(int x, int y);  
};


class multigraph {
public:
	int n_vertex;  //最新顶点数（加入copy的关节顶点）
	int n_vertex_org;  //原顶点数
	int m_edge;  //边数
	vector<vector<vector<int>>> edge_info;  
	vector<map<int,int> > adjlist_edge;  

	void getPmg(vector<char>* P);  //从异构图中获得Pmg
	void save_Mpl(vector<vector<int>>* Mpl);  //保存Mpl
};

void CBC(int s);  //保存v后的路径实例信息
void CBC_(int s);  //直接在pred中保存边权；pred用vector  【最快】
void CBC_1(int s);  //直接在pred中保存边权；pred用map 
void save_bc();
void getbcfromfile();
void compare_result();
void compare_bc();
int getMaxComponent();
int getComponentSize(int i, vector<int>* visited);

string file = "data_intro";
hetergraph hg;
multigraph Pmg;
vector<double> bc;
//vector<double> deta0;
vector<double> cbc;
vector<double> bbc;
ofstream out_result(file + "/al1_time_result.txt");
double timebfs = 0, timeback = 0;
double pmid;

int main() {
	//一、根据异构图计算Pmg
	hg.read_file();

	/*hg.show_hg();*/
	vector<char> P = { 'A','M','D','M','A' };
	Pmg.getPmg(&P);  //得到P-multigraph


	/*统计最大连通分量大小*/
	int maxComponent = getMaxComponent();

	bc.resize(Pmg.n_vertex);
	int count = 0;
	for (int i = 0; i < Pmg.n_vertex; i++) {
		CBC_(i);
		count++;
	}

	for (int i = 0; i < Pmg.n_vertex; i++) {
		/*if (bc[i] > 0) {
			cout << "顶点" << i << "的bc值：" << bc[i] << endl;
		}*/
		cout << "顶点" << i << "的bc值：" << bc[i] << endl;
	}
	/*compare_result();*/
	out_result.close();
	save_bc();
	return 0;
	
}


bool cmp(const pair<int, double>& a, const pair<int, double>& b) {
	return a.second > b.second;
}

void compare_bc() {
	vector<pair<int, double>> cbc_;
	for (int i = 0; i < cbc.size(); i++) {
		cbc_.push_back({ i,cbc[i] });
	}
	sort(cbc_.begin(), cbc_.end(), cmp);

	vector<pair<int, double>> bbc_;
	for (int i = 0; i < bbc.size(); i++) {
		bbc_.push_back({ i,bbc[i] });
	}
	sort(bbc_.begin(), bbc_.end(), cmp);
	
	/*输出排序后的结果*/
	ofstream outcbc(file + "/cbcorder.txt");
	if (!outcbc) {
		cerr << "file cbcorder.txt error" << endl;
		exit(1);
	}
	for (int i = 0; i < cbc_.size(); i++) {
		outcbc << cbc_[i].first << " " << cbc_[i].second << endl;
	}
	outcbc.close();
	ofstream outbbc(file + "/bbcorder.txt");
	if (!outbbc) {
		cerr << "file bbcorder.txt error" << endl;
		exit(1);
	}
	for (int i = 0; i < bbc_.size(); i++) {
		outbbc << bbc_[i].first << " " << bbc_[i].second << endl;
	}
	outbbc.close();
	/*输出前topn个的交差集合*/
	int topn = 10;
	set<int> bbctop;
	set<int> cbctop;

	for (int i = 0; i < topn; i++) {
		bbctop.insert(bbc_[i].first);
		cbctop.insert(cbc_[i].first);
	}
	set<int> result_inter;
	set<int> result_diffc_b;
	set<int> result_diffb_c;
	set_intersection(bbctop.begin(), bbctop.end(), cbctop.begin(), cbctop.end(), inserter(result_inter, result_inter.begin()));
	set_difference(bbctop.begin(), bbctop.end(), cbctop.begin(), cbctop.end(), inserter(result_diffb_c, result_diffb_c.begin()));
	set_difference(cbctop.begin(), cbctop.end(), bbctop.begin(), bbctop.end(), inserter(result_diffc_b, result_diffc_b.begin()));
	ofstream outcompare(file + "/comparecbcbbc.txt");
	if (!outcompare) {
		cerr << "file compare.txt error" << endl;
		exit(1);
	}
	outcompare << "交集：" << endl;
	for (int i : result_inter) {
		outcompare << i << " ";
	}
	outcompare << endl;
	outcompare << "c-b差集：" << endl;
	for (int i : result_diffc_b) {
		outcompare << i << " ";
	}
	outcompare << endl;
	outcompare << "b-c差集：" << endl;
	for (int i : result_diffb_c) {
		outcompare << i << " ";
	}
	outcompare << endl;
	outcompare.close();
}

void getbcfromfile() {
	string filename = file + "/cbc_result2.txt";  
	ifstream inputcbc(filename, ios::in);  //输入文件流对象input
	//判断文件是否存在
	if (!inputcbc) {
		cerr << "file error!" << endl;
		exit(1);
	}

	//读入所有顶点的cbc值//保存每个顶点的bc值的文件
	string cbc_info;
	while (getline(inputcbc, cbc_info)) {
		istringstream scbc_info(cbc_info);
		int v;
		double bcvalue;
		scbc_info >> v >> bcvalue;
		cbc.push_back(bcvalue);
	}
	inputcbc.close();

	string filename2 = file + "/bbc_result.txt";  
	ifstream inputbbc(filename2, ios::in);  //输入文件流对象input
	//判断文件是否存在
	if (!inputbbc) {
		cerr << "file error!" << endl;
		exit(1);
	}

	//读入所有顶点的cbc值//保存每个顶点的bc值的文件
	string bbc_info;
	while (getline(inputbbc, bbc_info)) {
		istringstream sbbc_info(bbc_info);
		int v;
		double bcvalue;
		sbbc_info >> v >> bcvalue;
		bbc.push_back(bcvalue);
	}
	inputbbc.close();
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

int getMaxComponent() {
	vector<int> visited(Pmg.n_vertex);
	int maxcsize = 0;
	for (int i = 0; i < Pmg.n_vertex; i++) {
		if (visited[i] == 0) {
			int csize = getComponentSize(i, &visited);
			maxcsize = max(maxcsize, csize);
		}
	}
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
	return totalcount;
}
//void CBC(int s) {
//	//initialize variable
//	/*cout << "源点" << s << endl;*/
//	queue<int> Q;
//	stack<int> S;
//	vector<vector<int>> pred(Pmg.n_vertex);
//	vector<double> dist(Pmg.n_vertex,DBL_MAX);
//	vector<double> pathnum(Pmg.n_vertex);
//	vector<map<int, int>> pathinstance(Pmg.n_vertex);  //每个点其下方路径实例的类型-数目信息，用于计算分母
//	dist[s] = 0;
//	pathnum[s] = 1;
//	Q.push(s);
//	//bfs
//	while (!Q.empty()) {
//		int v = Q.front();
//		Q.pop();
//		S.push(v);
//		for (auto w : Pmg.adjlist_edge[v]) {  //访问了v的所有邻居w[邻接顶点，邻接边]
//			if (dist[w[0]] == DBL_MAX) {  
//				dist[w[0]] = dist[v] + 1;
//				Q.push(w[0]);
//			}
//			if (dist[w[0]] == dist[v] + 1) {  //通过此if条件，可以得到v的所有下方的邻居w
//				/*pathnum[w[0]] += pathnum[v]*Pmg.Mp[v][w[0]];*/
//				pred[w[0]].push_back(v);    //更新w[0]的前置顶点，加入v
//				//统计v下方的路径实例信息
//				vector<vector<int>> einfo = Pmg.edge_info[w[1]];
//				for (auto pinfo : einfo) {
//					int pclass = pinfo[0];
//					int pnum = pinfo[1];
//					pathinstance[v][pclass] += pnum;
//				}
//			}
//		}
//		//遍历完v的邻居点后，才能统计完整的v后路径实例的信息，此时才能更新pathnum
//		for (auto w : Pmg.adjlist_edge[v]) {
//			if (dist[w[0]] == dist[v] + 1) { 
//				vector<vector<int>> einfo = Pmg.edge_info[w[1]];  //w[1]是v与w[0]之间的边
//				for (auto pinfo : einfo) {  //根据v与w[0]之间边的路径实例信息，更新pathnum[w[0]]
//					int pclass = pinfo[0];
//					int pnum = pinfo[1];
//					pathnum[w[0]] += pathnum[v] * ((double)pnum / pathinstance[v][pclass]);
//				}
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
//		for (int v : pred[w]) {
//			//根据v与w之间边的路径实例数信息，每类路径实例分别进行反向累积
//			int eindex = Pmg.Mp_edge[v][w];  //v与w之间边的索引号
//			vector<vector<int>> einfo = Pmg.edge_info[eindex];  //该边的信息
//			for (auto pinfo : einfo) {
//				int pclass = pinfo[0];
//				int pnum = pinfo[1];
//				deta[v] += (pathnum[v] * ((double)pnum / pathinstance[v][pclass]) / pathnum[w]) * (1 + deta[w]);
//			}
//		}
//		if (w != s) {
//			bc[w] += deta[w];
//		}
//	}
//	/*for (int i = 0; i < deta.size(); i++) {
//		cout << "对点" << i << "源依赖：" << deta[i] << endl;
//	}*/
//}

void CBC_(int s) {  //pred用vector
	//initialize variable
	/*cout << "源点" << s << endl;*/
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
		map<int, int> pathinstance;  //w下方所有路径实例按类统计
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
					pathinstance[pclass] += pnum;
				}
			}
		}
		//遍历完v的邻居点后，才能统计完整的v后路径实例的信息，获得v与其每个后置点w之间边的权重
		for (auto w : Pmg.adjlist_edge[v]) {
			double weight;
			weight = 0;
			if (dist[w.first] == dist[v] + 1) {  //计算v与w[0]之间的边权
				vector<vector<int>> einfo = Pmg.edge_info[w.second];  //w[1]是v与w[0]之间的边
				for (auto pinfo : einfo) {  //根据v与w[0]之间边的路径实例信息，更新pathnum[w[0]]
					int pclass = pinfo[0];
					int pnum = pinfo[1];
					weight += ((double)pnum / pathinstance[pclass]);  //累加v-w[0]之间的边权
					/*if (s == 1) {
						if (v == 2640 && pclass == 470) {
							cout << "2640下方到" << w[0] << "有" << pnum << "条" << pclass << "类路径实例" << endl;
							cout << pathinstance[pclass] << "pathinstance[pclass] +=" << pnum << endl;
						}
					}*/
				}
				//weight = weight * einfo.size();  //类别数
				weight = weight * einfo.size();
				pathnum[w.first] += pathnum[v] * weight;
				/*cout << weight << " ";*/
				pred[w.first].push_back({ (double)v,weight });
			}
		}
	}
	pmid = clock();
	//back propagation
	vector<double> deta(Pmg.n_vertex);
	while (!S.empty()) {
		int w = S.top();
		S.pop();
		for (vector<double> v : pred[w]) {
			deta[(int)v[0]] += (pathnum[(int)v[0]] * v[1] / pathnum[w]) * (1 + deta[w]);
			/*if (v[0] == 19) {
				cout << "w=" << w << "时，" << "deta[(int)" << 19 << "] += (" << pathnum[(int)v[0]] << " * " << v[1] << " / " << pathnum[w] << ") * (1 + " << deta[w] << ") "<< endl;
			}*/
			/*if (s == 1) {
				if ((int)v[0] == 0) {
					cout << "由" << w << "更新deta[0]:";
					cout << "deta[" << (int)v[0] << "] += (" << pathnum[(int)v[0]] << " * " << v[1] << " / " << pathnum[w] << ") * (1 + " << deta[w] << ") " << endl;
				}
			}*/
		}
		if (w != s) {
			bc[w] += deta[w];
		}
	}
	/*cout << "点"<<s<<"对点19的源依赖：" << deta[19] << endl;*/
	/*for (int i = 0; i < deta.size(); i++) {s
		cout << "对点" << i << "源依赖：" << deta[i] << endl;
	}*/
	/*deta0[s] = deta[0];*/
}

//void CBC_1(int s) {  //pred用map
//	//initialize variable
//	/*cout << "源点" << s << endl;*/
//	queue<int> Q;
//	stack<int> S;
//	vector<map<int,double>> pred(Pmg.n_vertex);
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
//					weight += ((double)pnum / pathinstance[pclass]);
//				}
//				pathnum[w[0]] += pathnum[v] * weight;
//				pred[w[0]][v]=weight;    //更新w[0]的前置顶点，加入v和v-w[0]之间的边权
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
//		for (auto v : pred[w]) {  //v表示 w的前置顶点v.first-w与v之间的边权v.second
//			deta[v.first]+= (pathnum[v.first] * v.second / pathnum[w]) * (1 + deta[w]);
//		}
//		if (w != s) {
//			bc[w] += deta[w];
//		}
//	}
//	/*for (int i = 0; i < deta.size(); i++) {
//		cout << "对点" << i << "源依赖：" << deta[i] << endl;
//	}*/
//
//}




void save_bc() {
	ofstream out(file + "/cbc_result.txt");
	for (int i = 0; i < Pmg.n_vertex; i++) {
		out << i << " ";
		out<< fixed << setprecision(8) << bc[i] << endl;
	}
	cout << "bc值已保存至" << file << "/result.txt";
}

/*读文件，获取异构图信息*/
void hetergraph::read_file() {
	string filename = file + "/base.txt";  
	ifstream input(filename, ios::in);  
	if (!input) {
		cerr << "file error!" << endl;
		exit(1);
	}

	string vinfo;
	getline(input, vinfo);
	istringstream svinfo(vinfo);
	svinfo >> n_vertexclass;

	//顶点类别vclass-顶点索引(i)  i类顶点数vcnum
	for (int i = 0; i < n_vertexclass; i++) {
		string vcinfo;
		getline(input, vcinfo);
		istringstream svcinfo(vcinfo);
		char vclass;
		int vcnum;
		svcinfo >> vclass >> vcnum;
		vertex_class_index[vclass] = i;  //顶点类别-索引号
		vertex_num.push_back(vcnum);    //每类顶点数
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
		edge_class.push_back({ vx,vy });  //第i类边连接的顶点类别索引
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
		//读保存第i类边信息的文件，并得到该类边的邻接矩阵
		string edgeinformation;
		getline(edge_input, edgeinformation);
		int nx = vertex_num[vx], ny = vertex_num[vy];
		map<int, map<int, int>> edge_matrix;
		string edgedata;
		while (getline(edge_input, edgedata)) {
			istringstream sedgedata(edgedata);
			int x, y;
			sedgedata >> x >> y;
			auto iterx = edge_matrix.find(x);  
			if (iterx == edge_matrix.end()) {  
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

	/*根据P查找P中边的索引，保存在eindex中*/
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
	cout << "找到边索引：";
	for (int i : eindex) cout << i << "-";
	cout << endl;
	/*计算Mpl*/
	int p = 0;
	map<int, map<int, int>> matrixpl = hg.edge_info[eindex[p]];
	hg.edge_info[eindex[p]].clear(); 
	while (++p < eindex.size()) {  
		int curr_eindex = eindex[p];  
		double ps11, pt11, time11;
		ps11 = clock();
		map<int, map<int, int>> tmatrix;  
		auto iterik = matrixpl.begin();
		for (iterik; iterik != matrixpl.end(); iterik++) {  
			int i = iterik->first;
			auto iterk = iterik->second.begin();  
			for (iterk; iterk != iterik->second.end(); iterk++) {
				int k = iterk->first;
				int currik = iterk->second;  
				auto iterkj = hg.edge_info[curr_eindex].find(k);
				if (iterkj != hg.edge_info[curr_eindex].end()) {  
					auto iterj = iterkj->second.begin();  
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

	int edge_index = 0;
	for (int i = 0; i < n_vertex; i++) {  //
	/*cout << "计算到第" << i;*/
		auto iterik = matrixpl.find(i);  
		if (iterik != matrixpl.end()) {  //
			map<int,vector<vector<int>>> ij_einfo;
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
							auto iterjeinfo = ij_einfo.find(j);  
							if (iterjeinfo == ij_einfo.end()) { 
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

