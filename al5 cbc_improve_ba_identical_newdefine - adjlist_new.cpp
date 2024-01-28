//识别identical顶点（最新定义，同一identical点集中，任意identical点之间的edgeinfo相同），Pmg中有adjlist
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
	int n_vertexclass;  //顶点类别数
	map<char, int> vertex_class_index;  //顶点类别-顶点编号
	vector<int> vertex_num;  //按顶点索引顺序存储每类顶点数

	int n_edgeclass;  //边的类别数
	vector<vector<int>> edge_class;  //每类边连接的顶点类别：顶点编号x—顶点编号y
	vector<int> edge_num;  //每类边数

	//vector<vector<vector<int>>> edge_info;  //3维数组，第一维edge_num，表示每类边；第二维和第三维表示边edge的邻接矩阵
	vector<map<int, map<int, int>>> edge_info;  //3维数组，第一维edge_num，表示每类边；第二维和第三维表示边edge的邻接表

	void read_file();
	void show_hg();
	int getEdgeIndex(int x, int y);  //根据输入的元路径P来匹配异构schema图上的边
};


class multigraph {
public:
	int n_vertex;  //最新顶点数（加入copy的关节顶点）
	int n_vertex_org;  //原顶点数
	int m_edge;  //边数
	vector<vector<vector<int>>> edge_info;  //存储每条边的信息[[class1index,num1],[class2index,num2]...]class即路径实例类别，即路径实例与哪个E类顶点相关
	//vector<vector<int>> Mp_edge;  //存储ai与aj之间边的索引??
	vector<map<int, int> > adjlist_edge;  //在每一行存储顶点的相邻顶点的同时，存储通过哪条边与相邻顶点相连
	//vector<vector<int> > Mp;  //??
	//vector<vector<int> > adjlist;

	void show_mg();
	void getPmg(vector<char>* P);  //从异构图中获得Pmg
	void getPmg_read_file();      //从保存的文件中读取Mpl和Pmg
	void save_Mpl(vector<vector<int>>* Mpl);  //保存Mpl
	void save_Pmg();  //保存整个Pmg
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
//vector<double> deta0;  //验证
vector<int> reach;
map<int, int> org;  //map[100]=1  表示顶点100在关节顶点划分前的原顶点是1
ofstream out_result(file + "/al1_time_result_amdma.txt");
int counta = 0, countb = 0;
int count_ident_edge = 0;
double pmid;
double weight_pathinstance_min = 100000;
double weight_pathinstance_max = 0;
double weight_D_min = 100000;
double weight_D_max = 0;
double weight_min = 100000;
double weight_max = 0;

int main() {
	//一、根据异构图计算Pmg
	double ps00, pt00, time00;
	ps00 = clock();
	hg.read_file();
	pt00 = clock();
	time00 = (double)(pt00 - ps00) / CLOCKS_PER_SEC;
	out_result << "读入异构图用时：" << time00 << endl;
	cout << "读入异构图用时：" << time00 << endl;
	/*hg.show_hg();*/
	double ps0, pt0, time0;
	ps0 = clock();
	vector<char> P = { 'A','M','W','M','A' };
	//vector<char> P = { 'A','P','V','P','A' };
	//vector<char> P = { 'B','U','B' };
	Pmg.getPmg(&P);  //得到P-multigraph
	pt0 = clock();
	time0 = (double)(pt0 - ps0) / CLOCKS_PER_SEC;
	cout << "得到Pmg用时：" << time0 << endl;
	out_result << "得到Pmg用时:" << time0 << endl;
	cout << "点数：" << Pmg.n_vertex << endl;
	out_result << "Pmg点数：" << Pmg.n_vertex << endl;
	cout << "边数：" << Pmg.m_edge << endl;
	out_result << "Pmg边数：" << Pmg.m_edge << endl;
	//Pmg.show_mg();
	//二、已经保存了Pmg，从文件中读取
	/*Pmg.getPmg_read_file();*/
	/*Pmg.show_mg();*/

	/*统计最大连通分量大小*/
	double pcmps, pcmpt, timecmp;
	pcmps = clock();
	int maxComponent = getMaxComponent();
	cout << "最大连通分量：" << maxComponent << endl;
	out_result << "Pmg最大连通分量：" << maxComponent << endl;
	pcmpt = clock();
	timecmp = (double)(pcmpt - pcmps) / CLOCKS_PER_SEC;
	out_result << "计算最大连通分量大小用时：" << timecmp << endl;


	bc.resize(Pmg.n_vertex);
	reach.resize(Pmg.n_vertex,1);
	/*deta0.resize(Pmg.n_vertex);*/
	
	double ps1, pt1, time1;
	ps1 = clock();
	BridgeEdgeDivide();  //桥边分割
	pt1 = clock();
	time1 = (double)(pt1 - ps1) / CLOCKS_PER_SEC;
	cout << "桥边分割完毕,共划分" << countb << "个桥边,总用时：" << time1 << endl;
	out_result << "桥边分割完毕,共划分(删除)" << countb << "个桥边,总用时：" << time1 << endl;

	double pcmpsb, pcmptb, timecmpb;
	pcmpsb = clock();
	int maxComponentb = getMaxComponent();
	cout << "桥边分割后，最大连通分量：" << maxComponentb << endl;
	out_result << "桥边分割后，最大连通分量：" << maxComponentb << endl;
	pcmptb = clock();
	timecmpb = (double)(pcmptb - pcmpsb) / CLOCKS_PER_SEC;
	out_result << "计算最大连通分量大小用时：" << timecmpb << endl;


	double ps2, pt2, time2;
	ps2 = clock();
	ArticulationVertexDivide();  //关节顶点分割
	bc.resize(Pmg.n_vertex);  //新加入的copy关节顶点的bc值的存储空间
	pt2 = clock();
	time2 = (double)(pt2 - ps2) / CLOCKS_PER_SEC;
	cout << "关节顶点分割完毕，共划分" << counta << "个关节顶点,总用时：" << time2 << endl;
	out_result << "关节顶点分割完毕，共划分" << counta << "个关节顶点,总用时：" << time2 << endl;

	double pcmps1, pcmpt1, timecmp1;
	pcmps1 = clock();
	int maxComponent1 = getMaxComponent();
	cout << "桥边、关节顶点分割后，最大连通分量：" << maxComponent1 << endl;
	out_result << "桥边、关节顶点分割后，最大连通分量：" << maxComponent1 << endl;
	pcmpt1 = clock();
	timecmp1 = (double)(pcmpt1 - pcmps1) / CLOCKS_PER_SEC;
	out_result << "计算最大连通分量大小用时：" << timecmp1 << endl;


	/*for (int i = 0; i < Pmg.n_vertex; i++) cout << i << ":" << bc[i] << endl;*/
	double ps4, pt4, time4;
	ps4 = clock();
	Ident();
	pt4 = clock();
	time4 = (double)(pt4 - ps4) / CLOCKS_PER_SEC;
	cout << "identical顶点识别与合并用时：" << time4 << endl;
	out_result << "identical顶点识别与合并用时：" << time4 << endl;

	double pcmps2, pcmpt2, timecmp2;
	pcmps2 = clock();
	int maxComponent2 = getMaxComponent();
	cout << "identical顶点合并完毕后，最大连通分量：" << maxComponent2 << endl;
	out_result << "identical顶点合并完毕后，最大连通分量：" << maxComponent2 << endl;
	pcmpt2 = clock();
	timecmp2 = (double)(pcmpt2 - pcmps2) / CLOCKS_PER_SEC;
	out_result << "计算最大连通分量大小用时：" << timecmp2 << endl;


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
		if(count%100==0)
		    cout << count<<"个顶点为源点遍历时间：" << time5<< endl;
	}

	//根据关节顶点的映射关系，计算关节顶点的bc值
	for (auto it : org) {
		bc[it.second] += bc[it.first];
	}
	double ptt = clock();

	cout << "前向bfs的时间：" << time5bfs << endl;
	cout << "反向累积的时间：" << time5back << endl;
	cout << "计算bc值的时间：" << time5 << endl;
	out_result << "前向bfs的时间：" << time5bfs << endl;
	out_result << "反向累积的时间：" << time5back << endl;
	out_result << "计算bc值的时间：" << time5 << endl;

	double time_all = (double)(ptt - ps00) / CLOCKS_PER_SEC;
	cout << "总用时：" << time_all << endl;
	out_result << "总用时：" << time_all << endl;

	/*for (int i = 0; i < Pmg.n_vertex; i++) {
		if (bc[i] > 0) {
			cout << "顶点" << i << "的bc值：" << bc[i] << endl;
		}
	}*/
	//compare_result();

	cout << "权重最大：" << weight_max << endl;
	cout << "权重最小：" << weight_min << endl;
	cout << "权重D最大：" << weight_D_max << endl;
	cout << "权重D最小：" << weight_D_min << endl;
	cout << "权重pathinstance最大：" << weight_pathinstance_max << endl;
	cout << "权重pathinstance最小：" << weight_pathinstance_min << endl;
	out_result << "权重最大：" << fixed << setprecision(8) << weight_max << endl;
	out_result << "权重最小：" << fixed << setprecision(8) << weight_min << endl;
	out_result << "权重D最大：" << fixed << setprecision(8) << weight_D_max << endl;
	out_result << "权重D最小：" << fixed << setprecision(8) << weight_D_min << endl;
	out_result << "权重pathinstance最大：" << fixed << setprecision(8) << weight_pathinstance_max << endl;
	out_result << "权重pathinstance最小：" << fixed << setprecision(8) << weight_pathinstance_min << endl;

	out_result.close();
	save_bc();
	return 0;
	
	//string filename = file + "/deta0_.txt";  //base中存储文件的基本信息，根据边类别的行对应的边的索引，找到边的文件"x.txt"
	//ifstream inputt(filename, ios::in);  //输入文件流对象input
	////判断文件是否存在
	//if (!inputt) {
	//	cerr << "file error!" << endl;
	//	exit(1);
	//}
	//vector<double> deta0_(Pmg.n_vertex_org);
	////读入所有顶点的deta0值
	//string deta0_info;
	//while (getline(inputt, deta0_info)) {
	//	istringstream sdeta0_info(deta0_info);
	//	int v;
	//	double deta0;
	//	sdeta0_info >> v >> deta0;
	//	deta0_[v] = deta0;
	//}
	//比较deta0值
	/*vector<int> different_v;
	for (int i = 0; i < Pmg.n_vertex_org; i++) {
		if (abs(deta0[i] - deta0_[i]) > 0.000001) different_v.push_back(i);
	}
	cout << "deta0值不同的顶点数：" << different_v.size() << endl;
	for (int v : different_v) {
		cout << v << ":" << deta0[v] << " " << deta0_[v] << endl;
	}
	ofstream outt(file + "/deta0_ident.txt");

	for (int i = 0; i < Pmg.n_vertex; i++) {
		cout << i << ":" << deta0[i] << " ";
		outt << i << " ";
		outt << fixed << setprecision(8) << deta0[i] << endl;
	}*/
	
}

//int main() {
//	getbcfromfile();
//	sort_bc();
//}

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

/*划分桥边、关节顶点*/
/*计算顶点i所在的连通分量，所有顶点的reach值之和*/
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
	dfn.resize(Pmg.n_vertex);  //dfn[i]记录了调用b_tarjan函数进行dfs时的顺序
	low.resize(Pmg.n_vertex);  //low[i]记录了点i的邻居中，最小的dfn值
	for (int i = 0; i < Pmg.n_vertex; i++) {
		stamp = 0;
		int totalreach;
		if (!dfn[i]) {
			totalreach = countTotalReach(i);  //计算i所在的连通分量中，顶点的总数，会有原图就是多个非连通分量的情况
			/*cout << "totalreach:" << totalreach << endl;*/
			/*outt << "！！！以" << i << "为源点,total reach:" << totalreach << endl;*/
			b_tarjan(i, -1, totalreach);  //对i所在的连通分量执行b_tarjan，划分桥边
		}
		/*cout << "stacklefgsize:" << stack_bi.size() << endl;*/
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
		/*outt  << index << "的邻居" << j << endl;*/
		/*int j = j_[0];*/
		if (j == fa) continue;
		if (!dfn[j]) {
			b_tarjan(j, index, totalreach);
			/*outt << "对" << j << "访问完毕" << endl;
			outt << "更新low[" << index << "]=min(" << low[index] << "," << low[j] << ")";*/
			low[index] = min(low[index], low[j]);
			/*outt << "low[" << index << "]=" << low[index] << endl;*/
			/*outt << "判断" << index << "-" << j << "是否是桥边" << endl;*/
			if (low[j] > dfn[index]) {
				//边index-j是桥边
				cout << "找到桥边：" << index << "-" << j << endl;
				countb++;
				/*outt << index << "-" << j << endl;*/
				int repj = 0, repindex = 0;
				/*cout << "栈：" << endl;*/
				while (stack_bi.top() != j) {  //计算j所在的分量中，顶点的reach之和
					int temp = stack_bi.top();
					/*cout << temp << endl;*/
					stack_bi.pop();
					repj += reach[temp];
				}
				/*cout << j << endl;*/
				repj += reach[j]; stack_bi.pop();
				repindex = totalreach - repj;  //index所在分量的顶点的reach之和，等于总reach-j所在分量的reach之和
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
			/*outt << j<<"在之前被访问过，通过dfn[" << j << "]更新low[" << index << "]=(min(" << low[index] << ", " << dfn[j] << "))  ";*/
			low[index] = min(low[index], dfn[j]);
			/*outt << "low[" << index << "]=" << low[index] << endl;*/
		}
	}
	map<int, int>().swap(indexadj);
}


//识别关节顶点，注意要在bc后增加新增顶点vi的bc值；在邻接矩阵中增加新增vi的行和列；在邻接表中用vi代替部分v；更新Pmg的n_vertex
/*划分关节顶点*/
void ArticulationVertexDivide() {
	dfn.resize(Pmg.n_vertex);
	low.resize(Pmg.n_vertex);
	for (int i = 0; i < Pmg.n_vertex_org; i++) {
		stamp = 0;
		int totalreach;
		if (!dfn[i]) {
			totalreach = countTotalReach(i);  //计算i所在的连通分量中，顶点的总数
			/*cout << "totalreach:" << totalreach << endl;*/
			a_tarjan(i, -1, totalreach);  //对i所在的连通分量执行b_tarjan，划分桥边
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
	map<int, int> indexadj = Pmg.adjlist_edge[index];  //因为下方在artvertexcopy中，可能会修改Pmg.adjlist[index]，所以先将其保存下来，在修改前的index的邻居上遍历
	for (auto iteradjj = indexadj.begin(); iteradjj != indexadj.end(); iteradjj++) {
		int j = iteradjj->first;
		if (j == fa) continue;
		if (!dfn[j]) {
			child++;
			a_tarjan(j, index, totalreach);
			low[index] = min(low[index], low[j]);
			if (low[j] >= dfn[index] && fa >= 0) {
				cout << "找到关节顶点：" << index << endl;
				counta++;
				artvertexcopy(index, j, totalreach);  //以index为割点，分割j所在的连通分量，将index-->copyindex
			}
			//当index为根节点时，遍历index的邻居，判断是否还会有另一个孩子，若没有，则index不是割点，若有，则是index-j所在分量的割点
			if (low[j] >= dfn[index] && fa < 0) {
				/*cout << j << " " << index << endl;*/
				auto nowj = iteradjj;
				nowj++;
				int isbi = 0;
				for (nowj; nowj != indexadj.end(); nowj++) {  //判断index的邻居是否都访问到了，即判断index是否还有其他邻居节点
					int nb = nowj->first;
					if (!dfn[nb]) {  //有未访问到的邻居，说明index是割点，当前栈中的所有顶点构成一个连通分量
						cout << "找到关节顶点(宽搜源点)：" << index << endl;
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
	//根节点为割点时
	//if (fa < 0 && child>1) {  //上方计算了当index不为根节点时，index为割点的连通分量，但当index为根节点，不能通过low[j]>=dfn[index]来判断index是否为割点
	//	
	//}
}


/*tarjan算法找到关节顶点后，复制关节顶点并更新邻接矩阵、邻接表、reach值*/
void artvertexcopy(int index, int j, int totalreach) {
	//顶点index是割点
	int copyindex = Pmg.n_vertex++;  //在连通分量中，copy一个index点
	cout << "复制关节点" << copyindex << endl;
	org[copyindex] = index;//保存copyindex的原点是index
	//更新邻接矩阵，插入新的顶点
	//Pmg.Mp.resize(Pmg.n_vertex);
	//Pmg.Mp_edge.resize(Pmg.n_vertex);
	//for (int i = 0; i < Pmg.n_vertex; i++) {
	//	Pmg.Mp[i].resize(Pmg.n_vertex);
	//	Pmg.Mp_edge[i].resize(Pmg.n_vertex);
	//}
	Pmg.adjlist_edge.resize(Pmg.n_vertex);
	//Pmg.adjlist.resize(Pmg.n_vertex);
	int repcopyindex = 0, repindex = 0;
	while (stack_bi.top() != j) {  //将index-j这个连通分量中的顶点出栈
		int temp = stack_bi.top();
		stack_bi.pop();
		//更新j代表的顶点数
		repcopyindex += reach[temp];
		/*更新邻接矩阵和邻接表，将连通分量中的顶点与index邻接，修改为与copyindex邻接*/
		auto itertemp = Pmg.adjlist_edge[index].find(temp);
		if (itertemp != Pmg.adjlist_edge[index].end()) {  //如果temp顶点在Pmultigraph中与index相连，则将其更新为与copyindex相连
			//更新index-temp的邻接矩阵
			int e_index = itertemp->second;  //???????????????
			//int e_index = Pmg.adjlist_edge[index][temp];

			//在Pmg中删除边index-temp
			Pmg.deleteEdge(index, temp);
			Pmg.adjlist_edge[copyindex].insert(map<int, int>::value_type(temp, e_index));
			Pmg.adjlist_edge[temp].insert(map<int, int>::value_type(copyindex, e_index));
			//Pmg.adjlist[copyindex].push_back(temp);
			//Pmg.adjlist[temp].push_back(copyindex);
		}
	}
	//最后将j出栈并更新j的代表的顶点数、邻接矩阵和邻接表
	stack_bi.pop();  repcopyindex += reach[j];
	/*Pmg.Mp[copyindex][j] = Pmg.Mp[index][j];
	Pmg.Mp[j][copyindex] = Pmg.Mp[j][index];*/
	int e_index1 = Pmg.adjlist_edge[index][j];
	/*Pmg.Mp_edge[copyindex][j] = e_index1;
	Pmg.Mp_edge[j][copyindex] = e_index1;*/
	Pmg.deleteEdge(index, j);
	Pmg.adjlist_edge[copyindex].insert(map<int, int>::value_type(j, e_index1));
	Pmg.adjlist_edge[j].insert(map<int, int>::value_type(copyindex, e_index1));
	/*Pmg.adjlist[copyindex].push_back(j);
	Pmg.adjlist[j].push_back(copyindex);*/

	repindex = totalreach - repcopyindex;  //repcopyindex是j所在连通分量上，顶点的reach值之和（不包含index顶点的reach）；repindex是分割连通分量后，index所在分量的reach值之和（包含index的reach）
	reach[index] += repcopyindex;  //更新顶点index的reach值，新加入j所在的分类上顶点的reach值之和，不同加index的reach，因为本来就有了
	reach.push_back(repindex);  //copyindex是新加入的顶点，所以要在reach中push_back；copyindex的reach值时index所在连通分量中顶点的reach之和(包括index顶点的reach)
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
	cout << "找到identical点集" << countset << "个，加入identical点集的点(不包括proxy)" << countvertex << "个" << endl;
	out_result << "找到identical点集" << countset << "个，加入identical点集的点(不包括proxy)" << countvertex << "个" << endl;
	cout << "合并identical点，删除边数：" << count_ident_edge << endl;
	out_result << "合并identical点，删除边数：" << count_ident_edge << endl;

	//identical点集中是否有任意两点对之间的边不相同
	//int countdiff = 0;
	//vector<int> different;
	//for (auto idset : IdentSet) {  //判断iset中每个顶点之间的edgeinfo是否相等
	//	vector<int> iset = { idset.first };
	//	for (int i : idset.second) iset.push_back(i);
	//	vector<vector<int>> einfo = Pmg.edge_info[Pmg.Mp_edge[iset[0]][iset[1]]];
	//	int diff = 0;
	//	for (int i = 0; i < iset.size(); i++) {
	//		for (int j = i + 1; j < iset.size(); j++) {
	//			int index1 = iset[i], index2 = iset[j];
	//			vector<vector<int>> compare_einfo = Pmg.edge_info[Pmg.Mp_edge[index1][index2]];
	//			if (einfo != compare_einfo) {
	//				diff++;
	//				break;
	//			}
	//		}
	//		if (diff > 0) break;
	//	}
	//	if (diff > 0) {
	//		different.push_back(iset[0]);
	//		countdiff++;
	//	}
	//}
	//cout << "countdiff:" << countdiff << endl;
	//for (int d : different) cout << d << " ";
	
}



//void todoType2Ident() {
//	for (int i = 0; i < Pmg.n_vertex; i++) {  //顶点i=0-n_vertex，寻找其identical点集
//		if (flag_addtoident[i] == 0) {
//			vector<int> mayidentset = { i };
//			//vector<int> iadjlist = Pmg.adjlist_edge[i];
//			//map<int, int> iadjlist = Pmg.adjlist_edge[i];  //不需要先单独存储，因为是先将i的ident点v存在mayidentset中，不修改adjli_edge[i]
//			auto iteriadj = Pmg.adjlist_edge[i].begin();
//			for (iteriadj; iteriadj != Pmg.adjlist_edge[i].end(); iteriadj++) {
//				int v = iteriadj->first;
//				if (v > i && flag_addtoident[v] == 0) {
//					int flag_isident = judgeIdent(i, v, 2);
//					if (flag_isident) {
//						/*flag_addtoident[v] = 1;
//						mergeIdent2(i, v);*/
//						mayidentset.push_back(v);
//					}
//				}
//			}
//			int flag_isident = 1;
//			if (mayidentset.size() > 2) {
//				auto iter00 = Pmg.adjlist_edge[mayidentset[0]].find(mayidentset[1]);
//				int edge_index0 = iter00->second;
//				vector<vector<int>> edgeinfo = Pmg.edge_info[edge_index0];
//				for (int x = 0; x < mayidentset.size(); x++) {
//					for (int y = x + 1; y < mayidentset.size(); y++) {
//						int indexx = mayidentset[x], indexy = mayidentset[y];
//						auto iter11 = Pmg.adjlist_edge[mayidentset[x]].find(mayidentset[y]);
//						int edge_index1 = iter11->second;
//						vector<vector<int>> compare_einfo = Pmg.edge_info[edge_index1];
//						if (edgeinfo != compare_einfo) {
//							flag_isident = 0;
//							break;
//						}
//					}
//					if (flag_isident == 0) break;
//				}
//			}
//			if (flag_isident) {
//				int proxy = mayidentset[0];
//				for (int x = 1; x < mayidentset.size(); x++) {
//					flag_addtoident[mayidentset[x]] = 1;
//					mergeIdent2(proxy, mayidentset[x]);
//				}
//			}
//			vector<int>().swap(mayidentset);
//		}
//	}
//}

void todoType2Ident() {
	for (int i = 0; i < Pmg.n_vertex; i++) {  //顶点i=0-n_vertex，寻找其identical点集
		if (flag_addtoident[i] == 0) {
			vector<int> iwithidentset = { i };
			//vector<int> iadjlist = Pmg.adjlist_edge[i];
			//map<int, int> iadjlist = Pmg.adjlist_edge[i];  //不需要先单独存储，因为是先将i的ident点v存在mayidentset中，不修改adjli_edge[i]
			auto iteriadj = Pmg.adjlist_edge[i].begin();
			for (iteriadj; iteriadj != Pmg.adjlist_edge[i].end(); iteriadj++) {
				int v = iteriadj->first;
				if (v > i && flag_addtoident[v] == 0) {
					int flag_isident = judgeIdent(i, v, 2);
					if (flag_isident) {
						/*flag_addtoident[v] = 1;
						mergeIdent2(i, v);*/
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
	if (ident[i] == 0) {  //说明v是第一个加到i为proxy的identical set的点,所以IdentSetEdgeindex中目前没有i
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
			/*cout << "i与v相同的邻居路径实例类别数：" << einfo_result.size() << endl;*/
			if (einfo_result.size() != i_neib_einfo.size()) return 0;
		}
		else {
			break;
		}
	}
	return 1;
	//set<int> i_adj{ Pmg.adjlist[i].begin(),Pmg.adjlist[i].end() };
	//set<int> v_adj{ Pmg.adjlist[v].begin(),Pmg.adjlist[v].end() };
	//if (type == 2) {  //type2类型的，要将i、v先从邻居中删除再判断
	//	i_adj.erase(v);
	//	v_adj.erase(i);
	//}
	//set<int> result;  //判断邻居是否相同
	//set_intersection(i_adj.begin(), i_adj.end(), v_adj.begin(), v_adj.end(), inserter(result, result.begin()));
	//if (result.size() != i_adj.size()) return 0;

	//for (int neib : result) {  //判断与邻居相连的边的信息是否都相同
	//	set<vector<int>> i_neib_einfo{ Pmg.edge_info[Pmg.Mp_edge[i][neib]].begin(),Pmg.edge_info[Pmg.Mp_edge[i][neib]].end() };
	//	set<vector<int>> v_neib_einfo{ Pmg.edge_info[Pmg.Mp_edge[v][neib]].begin(),Pmg.edge_info[Pmg.Mp_edge[v][neib]].end() };
	//	
	//	if (i_neib_einfo.size() != v_neib_einfo.size()) return 0;  //如果路径实例类别数不同，则不是
	//	set<vector<int>> einfo_result;
	//	set_intersection(i_neib_einfo.begin(), i_neib_einfo.end(), v_neib_einfo.begin(), v_neib_einfo.end(), inserter(einfo_result, einfo_result.begin()));
	//	/*cout << "i与v相同的邻居路径实例类别数：" << einfo_result.size() << endl;*/
	//	if (einfo_result.size() != i_neib_einfo.size()) return 0;

	//}
	//return 1;
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
	queue<int> temp_Q; //记录访问过的level层顶点

	stack<int> S;
	/*vector<map<int,double>> pred(Pmg.n_vertex);*/
	vector<vector<vector<double>>> pred(Pmg.n_vertex);
	vector<double> dist(Pmg.n_vertex, DBL_MAX);
	vector<double> pathnum(Pmg.n_vertex);
	dist[s] = 0;
	pathnum[s] = 1;
	int level = 0;

	Q.push(s);
	//bfs
	map<int, int> pathinstance;  //w下方所有路径实例按类统计

	while (!Q.empty()) {
		int v = Q.front();
		if (dist[v] == level) { //若v是level层的顶点，则访问v，将该层v点下方的路径实例信息都统计到pathinstance里面
			Q.pop();
			S.push(v);
			//统计level层路径实例类别
			temp_Q.push(v);   //记录第level层访问的顶点
			for (auto w : Pmg.adjlist_edge[v]) {  //访问了v的所有邻居w[邻接顶点，邻接边]
				if (dist[w.first] == DBL_MAX) {
					dist[w.first] = dist[v] + 1;
					Q.push(w.first);
				}
				if (dist[w.first] == dist[v] + 1) {  //通过此if条件，可以得到v的所有下方的邻居w
					//统计v下方的路径实例信息
					vector<vector<int>> einfo = Pmg.edge_info[w.second];
					for (auto pinfo : einfo) {
						int pclass = pinfo[0];
						int pnum = pinfo[1];
						pathinstance[pclass] += pnum * (1 + ident[w.first]);  //v--w[0]与v--w[0]代表的顶点的路径实例数
					}
				}
			}
		}
		else {  //若v是level+1层的，说明level层顶点已经访问完毕了；重新再访问level层的顶点rv，计算边权，更新最短路径数
			while (!temp_Q.empty()) {
				int rv = temp_Q.front();
				temp_Q.pop();
				//获得rv与其每个后置点w之间边的权重，并更新最短路径数
				for (auto w : Pmg.adjlist_edge[rv]) {
					double weight;
					weight = 0;
					if (dist[w.first] == dist[rv] + 1) {  //计算v与w[0]之间的边权
						vector<vector<int>> einfo = Pmg.edge_info[w.second];  //w[1]是v与w[0]之间的边
						for (auto pinfo : einfo) {  //根据v与w[0]之间边的路径实例信息，更新pathnum[w[0]]
							int pclass = pinfo[0];
							int pnum = pinfo[1];
							weight += ((double)pnum / (double)pathinstance[pclass]);
						}
						//weight = weight * (double)einfo.size();//1
						//weight = (double)einfo.size();//1
						if (weight_pathinstance_max < weight) weight_pathinstance_max = weight;
						if (weight_pathinstance_min > weight) weight_pathinstance_min = weight;
						if (weight_D_max < einfo.size()) weight_D_max = einfo.size();
						if (weight_D_min > einfo.size()) weight_D_min = einfo.size();
						weight = weight + (double)einfo.size();//1
						if (weight_max < weight) weight_max = weight;
						if (weight_min > weight) weight_min = weight;
						pathnum[w.first] += pathnum[rv] * weight * ((double)1 + ident[rv]);
						pred[w.first].push_back({ (double)rv,weight });
					}
				}
			}
			level = level + 1;  //接下来访问第level+1层
			pathinstance.clear(); //重新统计第level+1层的pathinstance信息
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
	queue<int> temp_Q; //记录访问过的level层顶点

	stack<int> S;
	/*vector<map<int,double>> pred(Pmg.n_vertex);*/
	vector<vector<vector<double>>> pred(Pmg.n_vertex);
	vector<double> dist(Pmg.n_vertex, DBL_MAX);
	vector<double> pathnum(Pmg.n_vertex);
	dist[s] = 0;
	pathnum[s] = 1;
	Q.push(s);
	int level = 0;

	//bfs
	map<int, int> pathinstance;  //w下方所有路径实例按类统计

	while (!Q.empty()) {
		int v = Q.front();
		if (dist[v] == level) { //若v是level层的顶点，则访问v，将该层v点下方的路径实例信息都统计到pathinstance里面
			Q.pop();
			S.push(v);
			//统计level层路径实例类别
			temp_Q.push(v);   //记录第level层访问的顶点
			for (auto w : Pmg.adjlist_edge[v]) {  //访问了v的所有邻居w[邻接顶点，邻接边]
				if (dist[w.first] == DBL_MAX) {
					dist[w.first] = dist[v] + 1;
					Q.push(w.first);
				}
				if (dist[w.first] == dist[v] + 1) {  //通过此if条件，可以得到v的所有下方的邻居w
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
		}
		else {  //若v是level+1层的，说明level层顶点已经访问完毕了；重新再访问level层的顶点rv，计算边权，更新最短路径数
			while (!temp_Q.empty()) {
				int rv = temp_Q.front();
				temp_Q.pop();
				//获得rv与其每个后置点w之间边的权重，并更新最短路径数
				for (auto w : Pmg.adjlist_edge[rv]) {
					double weight;
					weight = 0;
					if (dist[w.first] == dist[rv] + 1) {  //计算v与w[0]之间的边权
						vector<vector<int>> einfo = Pmg.edge_info[w.second];  //w[1]是v与w[0]之间的边
						for (auto pinfo : einfo) {  //根据v与w[0]之间边的路径实例信息，更新pathnum[w[0]]
							int pclass = pinfo[0];
							int pnum = pinfo[1];
							weight += ((double)pnum / (double)pathinstance[pclass]);
						}
						//weight = weight * (double)einfo.size();//1
						//weight = (double)einfo.size();//1
						if (weight_pathinstance_max < weight) weight_pathinstance_max = weight;
						if (weight_pathinstance_min > weight) weight_pathinstance_min = weight;
						if (weight_D_max < einfo.size()) weight_D_max = einfo.size();
						if (weight_D_min > einfo.size()) weight_D_min = einfo.size();
						weight = weight + (double)einfo.size();//1
						if (weight_max < weight) weight_max = weight;
						if (weight_min > weight) weight_min = weight;

						pathnum[w.first] += pathnum[rv] * weight * (1 + (double)ident[rv]);

						pred[w.first].push_back({ (double)rv,weight });
					}
				}
			}
			level = level + 1;  //接下来访问第level+1层
			pathinstance.clear(); //重新统计第level+1层的pathinstance信息
		}
		//int v = Q.front();
		//Q.pop();
		//S.push(v);
		//map<int, int> pathinstance;
		//for (auto w : Pmg.adjlist_edge[v]) {  //访问了v的所有邻居w[邻接顶点，邻接边]
		//	if (dist[w.first] == DBL_MAX) {
		//		dist[w.first] = dist[v] + 1;
		//		Q.push(w.first);
		//	}
		//	if (dist[w.first] == dist[v] + 1) {  //通过此if条件，可以得到v的所有下方的邻居w
		//		/*pathnum[w[0]] += pathnum[v]*Pmg.Mp[v][w[0]];*/
		//		//统计v下方的路径实例信息
		//		vector<vector<int>> einfo = Pmg.edge_info[w.second];
		//		for (auto pinfo : einfo) {
		//			int pclass = pinfo[0];
		//			int pnum = pinfo[1];
		//			pathinstance[pclass] += pnum * (1 + ident[w.first]);  //v--w[0]与v--w[0]代表的顶点的路径实例数
		//		}
		//	}
		//}
		//if (v == s) {  //当s为ident点集中的proxy点时，第一层bfs，除了访问邻居点w外，其实还有其identical点（也是s的邻居点）也在第二层顶点中，因为要统计完整的s下方路径实例的信息，所以也要加上s——sident的路径实例信息
		//	for (int sident : IdentSet[s]) {
		//		//int edge_index_s = Pmg.Mp_edge[s][sident];
		//		//!!!!!!!!int edge_index_s = Pmg.adjlist_edge[s][sident]; //??????????????????????????
		//		int edge_index_s = IdentSetEdgeindex[s][sident];
		//		vector<vector<int>> einfos = Pmg.edge_info[edge_index_s];
		//		for (auto pinfo : einfos) {
		//			int pclass = pinfo[0];
		//			int pnum = pinfo[1];
		//			pathinstance[pclass] += pnum * (1 + ident[sident]);
		//		}
		//	}
		//}
		////遍历完v的邻居点后，才能统计完整的v后路径实例的信息，获得v与其每个后置点w之间边的权重
		//for (auto w : Pmg.adjlist_edge[v]) {
		//	double weight = 0;
		//	if (dist[w.first] == dist[v] + 1) {  //计算v与w[0]之间的边权
		//		vector<vector<int>> einfo = Pmg.edge_info[w.second];  //w[1]是v与w[0]之间的边
		//		for (auto pinfo : einfo) {  //根据v与w[0]之间边的路径实例信息，更新pathnum[w[0]]
		//			int pclass = pinfo[0];
		//			int pnum = pinfo[1];
		//			weight += ((double)pnum / pathinstance[pclass]);  //累加v-w[0]之间的边权
		//		}
		//		weight = weight * einfo.size();
		//		pathnum[w.first] += pathnum[v] * weight * (1 + (double)ident[v]);
		//		pred[w.first].push_back({ (double)v,weight });
		//	}
		//}
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
//void CBC_ident1(int s) {  //pred用vector +ident后  s分为多个bfs树（pathnum）处理（identical顶点集中的点，之间的边的信息允许不同）
//	//initialize variable
//	/*cout << "源点" << s << endl;*/
//	vector<int> swithidentset = { s };
//	for (int sident : IdentSet[s]) swithidentset.push_back(sident);
//	queue<int> Q;
//	stack<int> S;
//	/*vector<map<int,double>> pred(Pmg.n_vertex);*/
//	vector<vector<vector<double>>> pred(Pmg.n_vertex);
//	vector<double> dist(Pmg.n_vertex, DBL_MAX);
//	vector<vector<double>> pathnum(swithidentset.size(),vector<double>(Pmg.n_vertex));
//	dist[s] = 0;
//	for (int si = 0; si < swithidentset.size(); si++) {
//		pathnum[si][swithidentset[si]] = 1;
//	}
//	//bfs0
//	int tempv = s;
//	S.push(tempv);
//	vector<map<int, int>> pathinstance_s(swithidentset.size());
//	for (auto w : Pmg.adjlist_edge[tempv]) {
//		if (dist[w[0]] == DBL_MAX) {
//			dist[w[0]] = dist[tempv] + 1;
//			Q.push(w[0]);
//		}
//		if (dist[w[0]] == dist[tempv] + 1) {
//			vector<vector<int>> einfo = Pmg.edge_info[w[1]];
//			for (auto pinfo : einfo) {
//				int pclass = pinfo[0];
//				int pnum = pinfo[1];
//				for (int i = 0; i < ident[s] + 1; i++) {
//					pathinstance_s[i][pclass] += pnum * (1 + ident[w[0]]);  //v--w[0]与v--w[0]代表的顶点的路径实例数
//				}
//			}
//		}
//	}
//	for (int si = 0; si < swithidentset.size(); si++) {  //根据idents中不同的点为源点，得到的temps下方的路径实例数，更新pathinstance_v
//		int temps = swithidentset[si];
//		for (int sident : swithidentset) {
//			if (sident != temps) {
//				int eindex = Pmg.Mp_edge[temps][sident];
//				vector<vector<int>> einfo = Pmg.edge_info[eindex];
//				for (auto pinfo : einfo) {
//					int pclass = pinfo[0];
//					int pnum = pinfo[1];
//					pathinstance_s[si][pclass] += pnum;  //因为sident一定只代表它自己，所以不需要*ident[sident](ident[w[0]])
//				}
//			}
//		}
//	}
//	for (auto w : Pmg.adjlist_edge[tempv]) {
//		if (dist[w[0]] == dist[tempv] + 1) {
//			vector<double> pred_weightset = { (double)tempv };
//			for (int si = 0; si < swithidentset.size(); si++) {
//				double weight = 0;
//				int temps = swithidentset[si];
//				int e_index = Pmg.Mp_edge[tempv][w[0]];
//				vector<vector<int>> einfo = Pmg.edge_info[e_index];
//				for (auto pinfo : einfo) {
//					int pclass = pinfo[0];
//					int pnum = pinfo[1];
//					weight += ((double)pnum / pathinstance_s[si][pclass]);
//				}
//				pathnum[si][w[0]] += pathnum[si][temps] * weight;  //不需要*ident[temps]（ident[v]），因为将不同的idents分成不同的pathnum数组
//			}
//			pred[w[0]].push_back(pred_weightset);
//		}
//	}
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
//					pathinstance[pclass] += pnum * (1 + ident[w[0]]);  //v--w[0]与v--w[0]代表的顶点的路径实例数
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
//				for (int si = 0; si < swithidentset.size(); si++) {
//					pathnum[si][w[0]] += pathnum[si][v] * weight * (1 + (double)ident[v]);
//				}
//				pred[w[0]].push_back({ (double)v,weight });
//			}
//		}
//	}
//	/*for (int i = 0; i < pathnum.size(); i++) {
//		cout << "到点" << i << "最短路径数：" << pathnum[i] << endl;
//	}*/
//	//back propagation
//	vector<vector<double>> deta(swithidentset.size(),vector<double>(Pmg.n_vertex));
//	while (!S.empty()) {
//		int w = S.top();
//		S.pop();
//		for (vector<double> v : pred[w]) {
//			if (v[0] != s) {
//				for (int si = 0; si < swithidentset.size(); si++) {
//					deta[si][(int)v[0]] += (pathnum[si][(int)v[0]] * v[1] / pathnum[si][w]) * (1 + deta[si][w]) * ((double)1 + ident[w]);
//				}
//			}
//			/*if (v[0] == 19) {
//				cout<<"w="<<w<<"时，"<< "deta[(int)"<<19<<"] += ("<<pathnum[(int)v[0]]<<" * "<<v[1]<<" / "<<pathnum[w]<<") * (1 + "<<deta[w]<<") * ((double)1 + "<<ident[w]<<")"<<endl;
//			}*/
//		}
//		if (w != s) {
//			for (int si = 0; si < swithidentset.size(); si++) {
//				bc[w] += deta[si][w];
//			}
//			if (ident[w] > 0) {
//				for (int wident : IdentSet[w]) {
//					for (int si = 0; si < swithidentset.size(); si++) {
//						bc[wident] += deta[si][w];
//					}
//				}
//			}
//			
//		}
//	}
//	/*for (int si = 0; si < swithidentset.size(); si++) {
//		cout << swithidentset[si]<<"对点0源依赖：" << deta[si][0] << endl;
//	}*/
//	for (int si = 0; si < swithidentset.size(); si++) {
//		deta0[si] = deta[si][0];
//	}
//
//	/*for (int i = 0; i < deta.size(); i++) {
//		cout << "对点" << i << "源依赖：" << deta[i] << endl;
//	}*/
//}
//void CBC_ri(int s) {  //pred用vector +reach+ident后
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
//				pathnum[w[0]] += pathnum[v] * weight * (1 + (double)ident[v]);
//				pred[w[0]].push_back({ (double)v,weight });
//			}
//		}
//	}
//	/*for (int i = 0; i < pathnum.size(); i++) {
//		cout << "到点" << i << "最短路径数：" << pathnum[i] << endl;
//	}*/
//	//back propagation
//	int reachs_total = reach[s];
//	if (ident[s] > 0) {
//		for (int sident : IdentSet[s]) {
//			reachs_total += reach[sident];
//		}
//	}
//	vector<double> deta(Pmg.n_vertex);
//	for (int i = 0; i < Pmg.n_vertex; i++) deta[i] = (double)reach[i] - 1;
//	while (!S.empty()) {
//		int w = S.top();
//		S.pop();
//		for (vector<double> v : pred[w]) {
//			double deta_add = (pathnum[(int)v[0]] * v[1] / pathnum[w]) * (1 + deta[w]);
//			if (ident[w] > 0) {
//				for (int wident : IdentSet[w]) {
//					deta_add+= (pathnum[(int)v[0]] * v[1] / pathnum[w]) * (1 + deta[wident]);
//				}
//			}
//			deta[(int)v[0]] += deta_add;
//			if (ident[(int)v[0]] > 0) {
//				for (int vident : IdentSet[(int)v[0]]) {
//					deta[vident] += deta_add;
//				}
//			}
//		}
//		if (w != s) {
//			bc[w] += deta[w]*reachs_total;
//			if (ident[w] > 0) {
//				for (int wident : IdentSet[w]) {
//					bc[wident] += deta[w]*reachs_total;
//				}
//			}
//		}
//	}
//	/*for (int i = 0; i < deta.size(); i++) {
//		cout << "对点" << i << "源依赖：" << deta[i] << endl;
//	}*/
//}

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
	ofstream out(file + "/cbc_result_amwma.txt");
	for (int i = 0; i < Pmg.n_vertex; i++) {
		out << i << " ";
		out<< fixed << setprecision(8) << bc[i] << endl;
	}
	cout << "由cbc计算出的bc值已保存至" << file << "/cbc_result_amwma.txt";
}

/*读文件，获取异构图信息*/
void hetergraph::read_file() {  //未加adjlist
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
		ifstream edge_input(edge_filename, ios::in);  //输入文件流对象input
		//判断文件是否存在
		if (!edge_input) {
			cerr << "file error!" << endl;
			exit(1);
		}
		//读保存第i类边信息的文件，并得到该类边的邻接矩阵
		string edgeinformation;
		getline(edge_input, edgeinformation);/***edge文件第一行是否有表示边的信息***/
		int nx = vertex_num[vx], ny = vertex_num[vy];
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
				xadjlist.insert(map<int, int>::value_type(y, 1));  //将y加入x的邻接list
				edge_matrix.insert(map<int, map<int, int>>::value_type(x, xadjlist)); //将x的邻接list加入edge_matrix
				xadjlist.clear();
			}
			else {  //x已经加入
				iterx->second.insert(map<int, int>::value_type(y, 1));
			}
		}
		edge_info.push_back(edge_matrix);  //第i类边的邻接矩阵
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
void multigraph::getPmg(vector<char>* P) {  //加adjlist
	vector<char>::iterator it = P->begin();  //it迭代遍历元路径P中的顶点类别
	/*Pmg顶点数*/
	n_vertex = hg.vertex_num[hg.vertex_class_index[*it]];
	n_vertex_org = n_vertex;
	adjlist_edge.resize(n_vertex);  //n个ai顶点的邻接顶点+边的编号
	//adjlist.resize(n_vertex);

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