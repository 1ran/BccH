//9.20update
//9.21update 在读入图hg时通过邻接表形式；在计算pmg时，通过邻接表计算；pmg中删除邻接矩阵，用邻接表存储边权
#include<iostream>
#include<vector>
#include<map>
#include<string>
#include<fstream>
#include<sstream>
#include<queue>
#include<stack>
#include<iomanip>
using namespace std;

class hetergraph {
public:
	int n_vertexclass;  //顶点类别数
	map<char, int> vertex_class_index;  //顶点类别-顶点编号
	vector<int> vertex_num;  //按顶点索引顺序存储每类顶点数

	int n_edgeclass;  //边的类别数
	vector<vector<int>> edge_class;  //每类边连接的顶点类别：顶点编号x—顶点编号y
	vector<int> edge_num;  //每类边数

	//vector<vector<vector<int>>> edge_info;  //3维数组，第一维edge_num，表示每类边；第二维和第三维表示边edge的邻接矩阵
	vector<map<int,map<int,int>>> edge_info;  //3维数组，第一维edge_num，表示每类边；第二维和第三维表示边edge的邻接表

	void read_file();
	void show_hg();
	int getEdgeIndex(int x, int y);
};


class multigraph {
public:
	int n_vertex;  //最新顶点数（加入copy的关节顶点）
	int n_vertex_org;  //原顶点数
	int m_edge;  //边数？
	//vector<vector<int> > Mp;
	vector<map<int,int> > adjlist;
	void show_mg();
	void getPmg(vector<char>* P);  //从异构图中获得Pmg
	void getPmg_read_file();      //从保存的文件中读取Mpl和Pmg
	void save_Mpl(vector<vector<int>>* Mpl);  //保存Mpl
	void save_Pmg();  //保存整个Pmg
};

void BBC(int s);
void save_bc();
void compare_result();
int getMaxComponent();
int getComponentSize(int i, vector<int>* visited);

string file = "data_intro";
hetergraph hg;
multigraph Pmg;
vector<double> bc;
ofstream out_result(file + "/al1_time_result.txt");  //表示算法1的时间结果
double timebfs = 0, timeback = 0;
double pmid;


int main() {
	//一、根据异构图计算Pmg
	double ps00, pt00, time00;
	ps00 = clock();
	hg.read_file();
	pt00 = clock();
	time00 = (double(pt00 - ps00)) / CLOCKS_PER_SEC;
	out_result << "读入异构图用时：" << time00 << endl;

	/*hg.show_hg();*/
	double ps0, pt0, time0;
	ps0 = clock();
	vector<char> P = { 'A','M','D','M','A' };
	Pmg.getPmg(&P);  //得到P-multigraph
	pt0 = clock();
	time0 = (double)(pt0 - ps0) / CLOCKS_PER_SEC;
	cout << "得到Pmg用时：" << time0 << endl;
	out_result << "得到Pmg用时:" << time0 << endl;
	cout << "点数：" << Pmg.n_vertex << endl;
	out_result << "Pmg点数：" << Pmg.n_vertex << endl;
	cout << "边数：" << Pmg.m_edge << endl;
	out_result << "Pmg边数：" << Pmg.m_edge << endl;
	//二、已经保存了Pmg，从文件中读取
	/*Pmg.getPmg_read_file();*/

	/*统计最大连通分量大小*/
	double pcmps, pcmpt, timecmp;
	pcmps = clock();
	int maxComponent=getMaxComponent();
	cout << "最大连通分量：" << maxComponent << endl;
	out_result << "Pmg最大连通分量：" << maxComponent << endl;
	pcmpt = clock();
	timecmp = (double)(pcmpt - pcmps) / CLOCKS_PER_SEC;
	out_result << "计算最大连通分量大小用时：" << timecmp << endl;


	bc.resize(Pmg.n_vertex);

	double ps1, pt1, time1=0;
	int count = 0;
	for (int i = 0; i < Pmg.n_vertex; i++) {
		count++;
		ps1 = clock();
		BBC(i);
		pt1 = clock();
		timebfs += (double)(pmid - ps1) / CLOCKS_PER_SEC;
		/*cout << "bfs:" << (double)(pmid - ps1) / CLOCKS_PER_SEC<<endl;*/
		timeback += (double)(pt1 - pmid) / CLOCKS_PER_SEC;
		/*cout<<"back:"<< (double)(pt1 - pmid) / CLOCKS_PER_SEC<<endl;*/
		time1 += (double)(pt1 - ps1) / CLOCKS_PER_SEC;
		if(count%100==0) cout << count << "个点计算需要" << time1 << endl;
	}
	cout << "前向bfs的时间：" << timebfs << endl;
	cout << "反向累积的时间：" << timeback << endl;
	cout << "计算bc值的时间："<<time1 << endl;
	out_result << "前向bfs的时间：" << timebfs << endl;
	out_result << "反向累积的时间：" << timeback << endl;
	out_result << "计算bc值的时间：" << time1 << endl;

	double time_all = (double)(pt1 - ps00) / CLOCKS_PER_SEC;
	cout << "总用时：" << time_all << endl;
	out_result << "总用时：" << time_all << endl;
	cout << "运行时间已存入" << file + "/time_result.txt" << endl;
	save_bc();
	compare_result();
	out_result.close();
	return 0;
}

void compare_result() {
	string filename = file + "/result.txt";  //base中存储文件的基本信息，根据边类别的行对应的边的索引，找到边的文件"x.txt"
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
void save_bc() {
	ofstream out(file + "/bbc_result.txt");
	for (int i = 0; i < Pmg.n_vertex; i++) {
		out << i << " ";
		out << fixed << setprecision(8) << bc[i] << endl;
	}
	cout << "由bbc计算出的bc值已保存至" << file << "/result.txt"<<endl;
}

void BBC(int s) {
	//initialize variable
	/*cout << "源点" << s << endl;*/
	queue<int> Q;
	stack<int> S;
	vector<vector<pair<int,int>>> pred(Pmg.n_vertex);  //或map
	vector<double> dist(Pmg.n_vertex,DBL_MAX);
	vector<double> pathnum(Pmg.n_vertex);
	dist[s] = 0;
	pathnum[s] = 1;
	Q.push(s);
	//bfs
	while (!Q.empty()) {
		int v = Q.front();
		Q.pop();
		S.push(v);
		for (auto w : Pmg.adjlist[v]){  //访问了v的所有邻居
			if (dist[w.first] == DBL_MAX) {  
				dist[w.first] = dist[v] + 1;
				Q.push(w.first);
			}
			if (dist[w.first] == dist[v] + 1) {  //通过此if条件，可以得到v的所有下方的邻居w，在此统计下方w反向累积到v时的权值
				pathnum[w.first] +=pathnum[v] * w.second;
				pred[w.first].push_back({ v,w.second });
			}

		}
	}
	queue<int>().swap(Q);
	vector<double>().swap(dist);
	pmid = clock();
	//back propagation
	vector<double> deta(Pmg.n_vertex);
	while (!S.empty()) {
		int w = S.top();
		S.pop();
		for (auto v : pred[w]) {
			deta[v.first] += (pathnum[v.first] * v.second / pathnum[w])*(1 + deta[w]);
		}
		if (w != s) {
			bc[w] += deta[w];
		}
	}
	stack<int>().swap(S);
	vector<vector<pair<int, int>>>().swap(pred);
	vector<double>().swap(pathnum);
	vector<double>().swap(deta);
}

int getMaxComponent() {
	vector<int> visited(Pmg.n_vertex);
	int maxcsize = 0;
	for (int i = 0; i < Pmg.n_vertex; i++) {
		if (visited[i] == 0) {
			int csize = getComponentSize(i,&visited);
			maxcsize = max(maxcsize, csize);
		}
	}
	return maxcsize;
}
int getComponentSize(int i,vector<int> *visited) {
	int count = 0;
	queue<int> Q;
	Q.push(i); 
	int totalcount=1;
	(*visited)[i] = 1;
	while (!Q.empty()) {
		int v = Q.front();
		Q.pop();
		for (auto w : Pmg.adjlist[v]) {
			if ((*visited)[w.first] == 0) {
				Q.push(w.first);
				(*visited)[w.first] = 1;
				totalcount += 1;
			}
		}
	}
	return totalcount;
}


/*读文件，获取异构图信息*/
void hetergraph::read_file() {
	string filename = file + "/base.txt";  //base中存储文件的基本信息，根据边类别的行对应的边的索引，找到边的文件"x.txt"
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
		cout << "开始读入文件" << edge_filename << endl;
		ifstream edge_input(edge_filename, ios::in);  //输入文件流对象input
		//判断文件是否存在
		if (!edge_input) {
			cerr << "file error!" << endl;
			exit(1);
		}
		//读保存第i类边信息的文件，并得到该类边的邻接矩阵
		string edgeinformation;
		getline(edge_input, edgeinformation);/***edge文件第一行是否有表示边的信息***/
		int nx = vertex_num[vx], ny = vertex_num[vy];  //nx：行数 ny：列数
		map<int, map<int, int>> edge_matrix;
		string edgedata;
		while (getline(edge_input, edgedata)) {
			istringstream sedgedata(edgedata);
			int x, y;  
			sedgedata >> x >> y;
			//第x行 第y列为1
			auto iterx = edge_matrix.find(x);  //x是否已经加入
			if (iterx == edge_matrix.end()) {  //尚未加入
				map<int, int> xadjlist;
				xadjlist.insert(map<int, int>::value_type(y,1));  //将y加入x的邻接list
				edge_matrix.insert(map<int, map<int, int>>::value_type(x, xadjlist)); //将x的邻接list加入edge_matrix
				xadjlist.clear();
			}
			else {  //x已经加入
				iterx->second.insert(map<int, int>::value_type(y, 1));
			}
		}
		edge_info.push_back(edge_matrix);  //第i类边的邻接表
		edge_matrix.clear();
		edge_input.close();
	}
	input.close();
}
/*根据边的两个顶点的索引号，查找边的索引号*/
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
	vector<char>::iterator it = P->begin();  //it迭代遍历元路径P中的顶点类别
	/*Pmg顶点数*/
	n_vertex = hg.vertex_num[hg.vertex_class_index[*it]];
	n_vertex_org = n_vertex;
	adjlist.resize(n_vertex);  
	
	/*根据P查找P中边的索引，保存在eindex中*/
	vector<int> eindex;
	for (int i = 0; i < P->size() / 2 ; i++) {  //一共遍历P->size()/2条边
		int x = hg.vertex_class_index[*it++];
		int y = hg.vertex_class_index[*it];
		int index = hg.getEdgeIndex(x, y);  //根据起点和终点类别，找到该类边对应的索引，从而得到该类边的邻接矩阵
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
	map<int,map<int,int>> matrixpl = hg.edge_info[eindex[p]];  //第1类边的邻接矩阵
	hg.edge_info[eindex[p]].clear();  //此时，第一个eindex的邻接矩阵已经保存在matrixpl中，信息不再需要
	while (++p < eindex.size()) {  //将matrix与下一个邻接矩阵edge_info[eindex[p]]相乘
		int curr_eindex = eindex[p];  //当前计算的边的索引,即当前用哪个边的邻接矩阵与MPL相乘
		cout << "计算"<<eindex[p-1]<<'*' << curr_eindex << "的邻接矩阵" << endl;
		double ps11, pt11, time11;
		ps11 = clock();
		map<int,map<int,int>> tmatrix;  //临时记录计算结果
		auto iterik = matrixpl.begin();
		for(iterik;iterik!=matrixpl.end();iterik++) {  //对每个有邻接序列的i
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
		cout << " 用时：" << time11;
		cout << endl;
		out_result << "计算" << eindex[p - 1] << '*' << curr_eindex << "的邻接矩阵" << " 用时：" << time11 << endl;
	}
	//save_Mpl(&matrixpl);  //保存Mpl
	eindex.clear();
	/*for (auto it : matrixpl) {
		cout << it.first << ":";
		for (auto itt : it.second) {
			cout << itt.first << "(" << itt.second << ")  ";
		}
		cout << endl;
	}*/

	/*计算出Mpl后，计算Mp+side*/
	cout << "计算adjlist" << endl;
	double ps12, pt12, time12;
	ps12 = clock();
	m_edge = 0;

	//得到Mpl邻接表的转置
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
	double pt13 = clock();
	double time13 = (double)(pt13 - ps12) / CLOCKS_PER_SEC;
	out_result << "计算mpl的转置用时：" << time13 << endl;

	for (int i = 0; i < n_vertex; i++) {  //计算第i个点的邻接序列
		cout << "计算到第" << i;
		auto iterik = matrixpl.find(i);  //找第i个点的k邻接序列
		if (iterik != matrixpl.end()) {  //找到了
			auto iterk = iterik->second.begin();
			for (iterk; iterk != iterik->second.end(); iterk++) {  //对于i的每个邻接点k
				int k = iterk->first;
				int currik = iterk->second;
				auto iterkj = matrixplt.find(k);//找k的邻接点j
				if (iterkj != matrixplt.end()) {
					auto iterj = iterkj->second.begin();
					for (iterj; iterj != iterkj->second.end(); iterj++) {
						int j = iterj->first;
						int currkj = iterj->second;
						if (j > i) {
							//向i的adjlist中，插入j（可能已经插入过，也可能没插入过）;同时向j的adjlist中插入i
							if (adjlist[i].find(j) == adjlist[i].end()) {
								m_edge++;
								adjlist[i].insert(map<int, int>::value_type(j, currik * currkj));
							}
							else {
								int oldij = adjlist[i][j];
								int nowij = oldij + currik * currkj;
								adjlist[i][j] = nowij;
							}
							if (adjlist[j].find(i) == adjlist[j].end()) {
								adjlist[j].insert(map<int, int>::value_type(i, currik * currkj));
							}
							else {
								int oldij = adjlist[j][i];
								int nowij = oldij + currik * currkj;
								adjlist[j][i] = nowij;
							}
						}	
					}
				}
			}
		}
	}
	matrixpl.clear(); matrixplt.clear();
	pt12 = clock();
	time12 = (double)(pt12 - ps12) / CLOCKS_PER_SEC;
	cout << "由mpl得到mp用时：" << time12 << endl;
	out_result<< "由mpl得到mp用时：" << time12 << endl;
	/*cout << "side vertex:";
	for (int i = 0; i < side.size(); i++) cout << side[i] << " ";
	cout << endl;*/
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

