//桥边+关节顶点+side②+identical①(桥边分割用的是j_index--的方式；关节顶点识别用的是提前保存邻居顶点的方式)
//side点的识别：位于计算Mp的过程中，遍历Mpl矩阵寻找side点；identical点的识别位于桥边、关节分割、side点移除后
//改动：1、void todoType1Ident() set<int> i2adjlist通过set得到非重复的二跳邻居序列；2、mergeIdent1()中，有些v有其他identical type2的顶点，在合并时，也加入到identSet[i]中
#include<iostream>
#include<vector>
#include<map>
#include<string>
#include<fstream>
#include<sstream>
#include<queue>
#include<stack>
#include<algorithm>
#include<time.h>
#include<ctime>
#include<set>
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
	vector<map<int, map<int, int>>> edge_info;  //3维数组，第一维edge_num，表示每类边；第二维和第三维表示边edge的邻接表

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
	vector<map<int, int> > adjlist;
	void show_mg();
	void getPmg(vector<char>* P);  //从异构图中获得Pmg
	void getPmg_read_file();      //从保存的文件中读取Mpl和Pmg
	void deleteEdge(int a, int b);  //将边ab从邻接表中删除
	void deleteVertex(int a);
	void save_Mpl(vector<vector<int>>* Mpl);  //保存Mpl
	void save_Pmg();  //保存整个Pmg
};

void BBC1(int s);
void BBC2(int s);
int countTotalReach(int i);
void BridgeEdgeDivide();
void b_tarjan(int index, int fa, int totalreach);
void ArticulationVertexDivide();
void a_tarjan(int index, int fa, int totalreach);
void artvertexcopy(int index, int j, int totalreach);
void SideAndIdent();
void deleteSide0();
void BCforSide0(vector<int> samesideset);
void todoType2Ident();
void todoType1Ident();
int judgeIdent(int i, int v, int type);
void mergeIdent1(int i, int v);
void mergeIdent2(int i, int v);
void find_side0(map<int, map<int,int>>* Mpl);
void compare_result();
void getbcfromfile();
void sort_bc();
int getMaxComponent();
int getComponentSize(int i, vector<int>* visited);
void save_bc();

string file = "data6";
hetergraph hg;
multigraph Pmg;
vector<double> bc;
vector<int> reach;
map<int, int> org;  //map[100]=1  表示顶点100在关节顶点划分前的原顶点是1
map<int, vector<int>> SideSet0;  //在计算Mp时得到的定义②的side顶点，即未对P-multigraph进行分割前，得到的side；
							  //side0[a]={Aset},a为P上中间的类别的顶点，Aset为A1类的顶点集，表示哪些A1类顶点只与顶点a相连
vector<int> flag_side;  //标志顶点是否作为side顶点被删除
map<int, vector<int>> IdentSet;  //顶点代表其他identical顶点的列表
vector<int> ident;     //顶点代表其他identical顶点的个数（不包括自己）
vector<int> flag_addtoident;  //已经加入identical点集的顶点
ofstream out_result(file + "/al7_time_result.txt");
double timebfs = 0, timeback = 0;
double pmid;
int countb = 0, counta = 0;
int count_side_edge = 0, count_ident_edge = 0;
int count_ident2 = 0, count_ident1 = 0;

int main() {
	//一、根据异构图计算Pmg
	double ps00, pt00, time00;
	ps00 = clock();
	hg.read_file();
	pt00 = clock();
	time00 = (double)(pt00 - ps00) / CLOCKS_PER_SEC;
	out_result << "读入异构图用时：" << time00 << endl;
	cout<< "读入异构图用时：" << time00 << endl;

	/*hg.show_hg();*/
	double ps0, pt0, time0;
	ps0 = clock();
	vector<char> P = { 'A','P','V','P','A' };
	Pmg.getPmg(&P);  //得到P-multigraph
	pt0 = clock();
	time0 = (double)(pt0 - ps0) / CLOCKS_PER_SEC;
	cout << "得到Pmg+side点集用时：" << time0 << endl;
	out_result << "得到Pmg+side点集用时:" << time0 << endl;
	cout << "点数：" << Pmg.n_vertex << endl;
	out_result << "Pmg点数：" << Pmg.n_vertex << endl;
	cout << "边数：" << Pmg.m_edge << endl;
	out_result << "Pmg边数：" << Pmg.m_edge << endl;
	//二、已经保存了Pmg，从文件中读取
	/*double ps0, pt0, time0;
	ps0 = clock();

	Pmg.getPmg_read_file();
	pt0 = clock();
	time0 = (double)(pt0 - ps0) / CLOCKS_PER_SEC;
	cout << "读入Pmg用时：" << time0 << endl;*/

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
	reach.resize(Pmg.n_vertex, 1);

	double ps1, pt1, time1;
	ps1 = clock();
	BridgeEdgeDivide();  //桥边分割
	pt1 = clock();
	time1 = (double)(pt1 - ps1) / CLOCKS_PER_SEC;
	cout << "桥边分割完毕,共划分"<<countb<<"个桥边,总用时：" << time1 << endl;
	out_result << "桥边分割完毕,共划分" << countb << "个桥边,总用时：" << time1 << endl;

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
	cout << "关节顶点分割完毕，共划分"<<counta<<"个关节顶点,总用时：" << time2 << endl;
	out_result << "关节顶点分割完毕，共划分" << counta << "个关节顶点,总用时：" << time2 << endl;

	double pcmps1, pcmpt1, timecmp1;
	pcmps1 = clock();
	int maxComponent1 = getMaxComponent();
	cout << "桥边、关节顶点分割后，最大连通分量：" << maxComponent1 << endl;
	out_result << "桥边、关节顶点分割后，最大连通分量：" << maxComponent1 << endl;
	pcmpt1 = clock();
	timecmp1 = (double)(pcmpt1 - pcmps1) / CLOCKS_PER_SEC;
	out_result << "计算最大连通分量大小用时：" << timecmp1 << endl;

	SideAndIdent();

	double pcmps2, pcmpt2, timecmp2;
	pcmps2 = clock();
	int maxComponent2 = getMaxComponent();
	cout << "移除side，合并identical后，最大连通分量：" << maxComponent2 << endl;
	out_result << "移除side，合并identical后，最大连通分量：" << maxComponent2 << endl;
	pcmpt2 = clock();
	timecmp2 = (double)(pcmpt2 - pcmps2) / CLOCKS_PER_SEC;
	out_result << "计算最大连通分量大小用时：" << timecmp2 << endl;


	double ps4, pt4, time4=0;
	int count = 0;
	for (int i = 0; i < Pmg.n_vertex; i++) {
		count++;
		if (flag_side[i] == 0 && flag_addtoident[i] == 0) {  //side+identical
			ps4 = clock();
			BBC2(i);
			pt4 = clock();
			timebfs += (double)(pmid - ps4) / CLOCKS_PER_SEC;
			timeback += (double)(pt4 - pmid) / CLOCKS_PER_SEC;
			time4 += (double)(pt4 - ps4) / CLOCKS_PER_SEC;
		}
		
		//if (flag_side[i] == 0) {  //side
		//	ps4 = clock();
		//	/*cout << "以" << i << "为源点宽搜" << endl;*/
		//	BBC1(i);
		//	pt4 = clock();
		//	timebfs += (double)(pmid - ps4) / CLOCKS_PER_SEC;
		//	timeback += (double)(pt4 - pmid) / CLOCKS_PER_SEC;
		//	time4 += (double)(pt4 - ps4) / CLOCKS_PER_SEC;
		//}

		//ps4 = clock();
		//BBC1(i);  //关节顶点+桥边
		//pt4 = clock();
		//timebfs += (double)(pmid - ps4) / CLOCKS_PER_SEC;
		//timeback += (double)(pt4 - pmid) / CLOCKS_PER_SEC;
		//time4 += (double)(pt4 - ps4) / CLOCKS_PER_SEC;

		if (count % 10 == 0) {
			cout << count << "个点为源点用时：" << time4 << endl;
		}
	}

	//根据关节顶点的映射关系，计算关节顶点的bc值
	for (auto it : org) {
		bc[it.second] += bc[it.first];
	}
	double ptt = clock();
	cout << "前向bfs的时间：" << timebfs << endl;
	cout << "反向累积的时间：" << timeback << endl;
	cout << "计算bc值的时间：" << time4 << endl;
	out_result << "前向bfs的时间：" << timebfs << endl;
	out_result << "反向累积的时间：" << timeback << endl;
	out_result << "计算bc值的时间：" << time4 << endl;

	double time_all = (double)(ptt - ps00) / CLOCKS_PER_SEC;
	cout << "总用时：" << time_all << endl;
	out_result << "总用时：" << time_all << endl;
	cout << "运行时间已存入" << file + "/improve_time_result.txt" << endl;

	/*for (int i = 0; i < Pmg.n_vertex_org; i++) {
		if (bc[i] > 0) {
			cout << "顶点" << i << "的bc值：" << bc[i] << endl;
		}
	}*/
	/*compare_result();*/
	save_bc();
	out_result.close();
	return 0;
}

//int main() {
//	getbcfromfile();
//	sort_bc();
//}

void save_bc() {
	ofstream out(file + "/bbc_result.txt");
	for (int i = 0; i < Pmg.n_vertex; i++) {
		out << i << " ";
		out << fixed << setprecision(8) << bc[i] << endl;
	}
	cout << "由bbc计算出的bc值已保存至" << file << "/result.txt" << endl;
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
	vector<double>().swap(bc0);
	vector<int>().swap(different_v);
}

bool cmp(const pair<int, double>& a, const pair<int, double>& b) {
	return a.second > b.second;
}

void sort_bc() {
	vector<pair<int, double>> bc_;
	for (int i = 0; i < bc.size(); i++) {
		bc_.push_back({ i,bc[i] });
	}
	sort(bc_.begin(), bc_.end(), cmp);
	cout << "top10:" << endl;
	for (int i = 0; i < 10; i++) {
		cout << "顶点" << bc_[i].first << "的bc值:" << bc_[i].second << endl;
	}
	vector<pair<int, double>>().swap(bc_);
}

void getbcfromfile() {
	string filename = file + "/result.txt";  //base中存储文件的基本信息，根据边类别的行对应的边的索引，找到边的文件"x.txt"
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
		for (auto w : Pmg.adjlist[v]) {
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
		for (auto w : Pmg.adjlist[v]) {
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
			totalreach = countTotalReach(i);  //计算i所在的连通分量中，顶点的总数
			/*cout << "totalreach:" << totalreach << endl;*/
			b_tarjan(i, -1, totalreach);  //对i所在的连通分量执行b_tarjan，划分桥边
		}
		/*cout << "stacklefgsize:" << stack_bi.size() << endl;*/
		stack<int>().swap(stack_bi);
	}
	dfn.clear(); low.clear();
}

void b_tarjan(int index, int fa, int totalreach) {
	int child = 0;
	dfn[index] = ++stamp;
	low[index] = stamp;
	stack_bi.push(index);
	map<int, int> indexadj = Pmg.adjlist[index];
	for (auto adjj:indexadj)
	/*for(int j:Pmg.adjlist[index])*/{
		int j = adjj.first;
		if (j == fa) continue;
		if (!dfn[j]) {
			b_tarjan(j, index, totalreach);
			low[index] = min(low[index], low[j]);
			if (low[j] > dfn[index]) {
				//边index-j是桥边
				cout << "找到桥边：" << index << "-" << j << endl;
				countb++;
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

				/*cout << "rep" << j << ":" << repj << "  rep" << index << ":" << repindex << endl;
				cout<<"reach"<<index<<":"<<reach[index]<< "  reach" << j << ":" << reach[j] << endl;*/


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
	
	map<int, int> ().swap(indexadj);
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
	map<int,int> indexadj = Pmg.adjlist[index];  //因为下方在artvertexcopy中，可能会修改Pmg.adjlist[index]，所以先将其保存下来，在修改前的index的邻居上遍历
	/*for (int indexj = 0; indexj < indexadj.size(); indexj++)*/
	for (auto iteradjj = indexadj.begin(); iteradjj != indexadj.end();iteradjj++) {
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
				for (nowj; nowj!=indexadj.end(); nowj++) {  //判断index的邻居是否都访问到了，即判断index是否还有其他邻居节点
					int nb = nowj->first;
					if (!dfn[nb]) {  //有未访问到的邻居，说明index是割点，当前栈中的所有顶点构成一个连通分量
						cout << "找到关节顶点(宽搜源点)：" << index << endl;
						counta++;
						isbi = 1; break;
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
	//根节点为割点时
	//if (fa < 0 && child>1) {  //上方计算了当index不为根节点时，index为割点的连通分量，但当index为根节点，不能通过low[j]>=dfn[index]来判断index是否为割点
	//	
	//}
	map<int, int>().swap(indexadj);
}

/*tarjan算法找到关节顶点后，复制关节顶点并更新邻接矩阵、邻接表、reach值*/
void artvertexcopy(int index, int j, int totalreach) {
	//顶点index是割点
	int copyindex = Pmg.n_vertex++;  //在连通分量中，copy一个index点
	cout << "复制关节点" << copyindex << endl;
	org[copyindex] = index;//保存copyindex的原点是index

	Pmg.adjlist.resize(Pmg.n_vertex);  //扩展adjlist
	int repcopyindex = 0, repindex = 0;
	while (stack_bi.top() != j) {  //将index-j这个连通分量中的顶点出栈
		int temp = stack_bi.top();
		stack_bi.pop();
		//更新j代表的顶点数
		repcopyindex += reach[temp];
		/*更新邻接矩阵和邻接表，将连通分量中的顶点与index邻接，修改为与copyindex邻接*/
		auto itertemp = Pmg.adjlist[index].find(temp);
		if (itertemp != Pmg.adjlist[index].end()) {  //如果temp顶点在Pmultigraph中与index相连，则将其更新为与copyindex相连
			int weightindex_temp = Pmg.adjlist[index][temp];
			Pmg.deleteEdge(index, temp);
			Pmg.adjlist[copyindex].insert(map<int, int>::value_type(temp, weightindex_temp));
			Pmg.adjlist[temp].insert(map<int, int>::value_type(copyindex, weightindex_temp));
		}
		//if (Pmg.Mp[temp][index] != 0) {  //如果temp顶点在Pmultigraph中与index相连，则将其更新为与copyindex相连
		//	//更新index-temp的邻接矩阵
		//	Pmg.Mp[temp][copyindex] = Pmg.Mp[temp][index];
		//	Pmg.Mp[copyindex][temp] = Pmg.Mp[temp][index];

		//	//在Pmg中删除边index-temp
		//	Pmg.deleteEdge(index, temp);
		//	Pmg.adjlist[copyindex].push_back(temp);
		//	Pmg.adjlist[temp].push_back(copyindex);
		//}
	}
	//最后将j出栈并更新j的代表的顶点数、邻接矩阵和邻接表
	stack_bi.pop();  repcopyindex += reach[j];
	int weightindex_j = Pmg.adjlist[index][j];
	//Pmg.Mp[copyindex][j] = Pmg.Mp[index][j]; Pmg.Mp[j][copyindex] = Pmg.Mp[j][index];
	Pmg.deleteEdge(index, j);
	Pmg.adjlist[copyindex].insert(map<int, int>::value_type(j, weightindex_j));
	Pmg.adjlist[j].insert(map<int, int>::value_type(copyindex, weightindex_j));

	repindex = totalreach - repcopyindex;  //repcopyindex是j所在连通分量上，顶点的reach值之和（不包含index顶点的reach）；repindex是分割连通分量后，index所在分量的reach值之和（包含index的reach）
	reach[index] += repcopyindex;  //更新顶点index的reach值，新加入j所在的分类上顶点的reach值之和，不同加index的reach，因为本来就有了
	reach.push_back(repindex);  //copyindex是新加入的顶点，所以要在reach中push_back；copyindex的reach值时index所在连通分量中顶点的reach之和(包括index顶点的reach)
}


/*side顶点和identical顶点*/
void SideAndIdent() {
	/*输出side顶点集*/
	int count_sideset = 0, count_sidenum = 0;
	for (auto it : SideSet0) {
		count_sideset++;
		/*cout << it.first << ":";*/
		for (int itt : it.second) {
			count_sidenum++;
			/*cout << itt << " ";*/
		}
		/*cout << endl;*/
	}
	cout << "side集合数" << count_sideset << ",side顶点数" << count_sidenum << endl;
	out_result<< "side集合数" << count_sideset << ",side顶点数" << count_sidenum << endl;
	/*移除side*/
	flag_side.resize(Pmg.n_vertex);//作为顶点是否为side点的标志，在最后对每个点执行bbc，不用再以side点为源点计算
	double psides, psidet, timeside;
	psides = clock();
	deleteSide0();  //delete
	SideSet0.clear();
	psidet = clock();
	timeside = (double)(psidet - psides) / CLOCKS_PER_SEC;
	cout << "移除side顶点，用时" << timeside << endl;
	out_result << "移除side顶点，用时" << timeside << endl;
	cout << "移除side顶点所删除的边数：" << count_side_edge << endl;
	out_result << "移除side顶点所删除的边数：" << count_side_edge << endl;

	double pidents, pidentt, timeident;
	double pidentmid;
	pidents = clock();
	ident.resize(Pmg.n_vertex);
	flag_addtoident.resize(Pmg.n_vertex);//作为顶点是否已被proxy顶点代表的标志，在最后对每个点执行bbc，不用再以该点为源点计算
	cout << "开始识别identical2" << endl;
	todoType2Ident();
	pidentmid = clock();
	double timemid = (double)(pidentmid - pidents) / CLOCKS_PER_SEC;
	cout << "其中，合并识别identical2用时：" << timemid << endl;
	out_result << "其中，合并识别identical2用时：" << timemid << endl;
	cout << "开始识别identical1" << endl;
	todoType1Ident();
	pidentt = clock();
	timeident = (double)(pidentt - pidents) / CLOCKS_PER_SEC;
	cout << "识别、合并identical顶点，用时" << timeident << endl;
	out_result<< "识别、合并identical顶点，用时" << timeident << endl;
	cout << "合并identical顶点所删除的边数：" << count_ident_edge << endl;
	out_result << "合并identical顶点所删除的边数：" << count_ident_edge << endl;

	//输出identical顶点集
	int count_identset = 0, count_identnum = 0;
	for (auto it : IdentSet) {
		count_identset++; 
		/*cout << it.first << ":";
		cout << "(" << reach[it.first] << ")";*/
		for (int ident : it.second)
		{
			/*cout << ident << " ";
			cout << "(" << reach[ident] << ")";*/
			count_identnum++;
		}
		/*cout << endl;*/
	}
	cout << "ident集合数" << count_identset << ",合并了的ident顶点数（不包含proxy点）" << count_identnum << endl;
	out_result << "ident集合数" << count_identset << ",合并了的ident顶点数（不包含proxy点）" << count_identnum << endl;
	cout << "ident1:" << count_ident1 << "  ident2:" << count_ident2 << endl;
	out_result << "ident1:" << count_ident1 << "  ident2:" << count_ident2 << endl;

}

//删除最初识别的side顶点
void deleteSide0() {
	double pdss, pdst, timeds=0;
	int sidecount = 0;
	for (auto it : SideSet0) {  //it为map中每个same_side点集
		//对same_side点集中的每个点对的合并更新bc值
		pdss = clock();
		sidecount++;
		for (int i = 0, len = it.second.size(); i < len; i++) {
			if (len == 2 && Pmg.adjlist[it.second[0]].find(it.second[1])==Pmg.adjlist[it.second[0]].end()) break;  //当最初识别的side点集中，有一个导演指导两个演员的孤立团，此时两个演员是side点，
																			//但在桥边删除时，已经将演员之间的边删除，此时两演员不再是同一个点集中的side点，而是两个独立的side点，合并不再有pair-dependency损失
			if (len == 1 && Pmg.adjlist[it.second[0]].size() == 0) break;  //对孤立的点（一，本来就孤立；二，因为桥边分割degree-1点导致孤立）
			int u = it.second[i];
			for (int j = i + 1; j < len; j++) {
				int v = it.second[j];
				bc[v] += ((double)reach[v] - 1) * reach[u];
				bc[u] += ((double)reach[u] - 1) * reach[v];
			}
		}
		//以same_side点集为源点进行宽搜，累加点集中的点对其他点的源依赖
		BCforSide0(it.second);
		//删除it代表的same_side点集
		for (int s : it.second) {
			count_side_edge += Pmg.adjlist[s].size();
			Pmg.deleteVertex(s);
		}
		pdst = clock();
		timeds += (double)(pdst - pdss) / CLOCKS_PER_SEC;
		cout << "移除1个side集合用时：" << timeds << endl;
		
	}
}

void BCforSide0(vector<int> samesideset) {
	int sidesetsize = samesideset.size();
	/*for (int s : samesideset) flag_side[s] = 1;*/
	queue<int> Q;
	stack<int> S;
	vector<vector<pair<int,int>>> pred(Pmg.n_vertex);  //samesideset中每个顶点为源点宽搜得到的pred都相同，因为宽搜树相同
	vector<double> dist(Pmg.n_vertex, DBL_MAX);  //samesideset中每个顶点为源点宽搜得到的dist也相同，因为宽搜树相同
	vector<vector<double>> pathnum(sidesetsize, vector<double>(Pmg.n_vertex));  //samesideset中每个顶点为源点宽搜得到的到其他顶点的最短路径数不同，因为边的权值不同
	for (int s : samesideset) dist[s] = 0;
	for (int i = 0; i < sidesetsize; i++) {
		for (int s : samesideset)
			pathnum[i][s] = 1;
	}

	//bfs0
	//对samesideset中的顶点为起点宽搜时，选择samesideset[0]为代表点，且在寻找其邻居顶点构造宽搜树时，不能是samesideset中的顶点，因为已经合并
	int tempv = samesideset[0];
	S.push(tempv);
	for (auto w : Pmg.adjlist[samesideset[0]]) {
		int flag_wnotin = 1;
		for (int i = 0; i < sidesetsize; i++) {
			if (samesideset[i] == w.first) {
				flag_wnotin = 0; break;
			}
		}
		//当w不在samesideset里面时
		if (flag_wnotin) {
			if (dist[w.first] == DBL_MAX) {
				dist[w.first] = dist[tempv] + 1;
				Q.push(w.first);
			}
			if (dist[w.first] == dist[tempv] + 1) {
				for (int i = 0; i < sidesetsize; i++) {  //以sideset中每个顶点为源点进行宽搜得到的最短路径数
					pathnum[i][w.first] += pathnum[i][samesideset[i]] * Pmg.adjlist[samesideset[i]][w.first];  //这里的权重是samesideset[i]与w.first之间的权重
				}
				pred[w.first].push_back({ tempv,w.second });
			}
		}
	}
	//bfs
	while (!Q.empty()) {
		int v = Q.front();
		Q.pop();
		S.push(v);
		for (auto w : Pmg.adjlist[v]) {
			if (dist[w.first] == DBL_MAX) {
				dist[w.first] = dist[v] + 1;
				Q.push(w.first);
			}
			if (dist[w.first] == dist[v] + 1) {
				for (int i = 0; i < sidesetsize; i++) {  //以sideset中每个顶点为源点进行宽搜得到的最短路径数
					pathnum[i][w.first] += pathnum[i][v] * w.second;
				}
				pred[w.first].push_back({ v,w.second });
			}
		}
	}
	/*for (int i = 0; i < pathnum.size(); i++) {
		cout << "到点" << i << "最短路径数：" << pathnum[i] << endl;
	}*/
	queue<int>().swap(Q);
	vector<double>().swap(dist);
	//back propagation
	vector<vector<double>> deta;
	vector<double> deta_reach(Pmg.n_vertex);
	for (int i = 0; i < Pmg.n_vertex; i++) deta_reach[i] = (double)reach[i] - 1;
	for (int i = 0; i < sidesetsize; i++) deta.push_back(deta_reach);
	vector<double>().swap(deta_reach);
	while (!S.empty()) {
		int w = S.top();
		S.pop();
		for (auto v : pred[w]) {  //当v是samesideset[0]时，对deta[i][v]的更新其实不同，因为v代表了samesideset中所有顶点，但顶点的Mp[v][w]不同；但不用单独讨论，因为当w是samesideset[0]时，不累加bc值
			if (v.first == samesideset[0]) {
				//第i组的pathnum和deta，用来更新samesideset[i]的deta值；本来(else中)是对每个i，都更新deta[i][v]，但当v是samesideset时，要更新deta[i][samesetside[i]]
				for (int i = 0; i < sidesetsize; i++) {  //根据sameside[i]这一组pathnum和deta对v进行更新
					deta[i][samesideset[i]] += (pathnum[i][samesideset[i]] * Pmg.adjlist[samesideset[i]][w] / pathnum[i][w]) * (1 + deta[i][w]);
				}
			}
			else {
				for (int i = 0; i < sidesetsize; i++) {  //根据sameside[i]这一组pathnum和deta对v进行更新
					deta[i][v.first] += (pathnum[i][v.first] * v.second / pathnum[i][w]) * (1 + deta[i][w]);
				}
			}

			/*deta[v] += (pathnum[v] * Pmg.Mp[v][w] / pathnum[w]) * (1 + deta[w]);*/
		}
		if (w != samesideset[0]) {
			for (int i = 0; i < sidesetsize; i++) {
				bc[w] += deta[i][w] * reach[samesideset[i]] + reach[samesideset[i]] * (deta[i][w] - ((double)reach[w] - 1));
				/*if (w == 17) {
					cout << "源点" << samesideset[i] << "对17的源依赖：" << deta[i][w] * reach[samesideset[i]] + reach[samesideset[i]] * (deta[i][w] - ((double)reach[w] - 1)) << endl;
				}*/
			}

			/*bc[w] += deta[w] * reach[s];*/
		}
	}
	for (int i = 0; i < sidesetsize; i++) {
		int ss = samesideset[i];
		/*if (ss == 17) {
			cout << "顶点17删除后，对17的bc值的影响：" << ((double)reach[ss] - 1) * (deta[i][ss] - ((double)reach[ss] - 1)) <<endl;
			cout << bc[ss] << endl;
		}*/
		bc[ss] += ((double)reach[ss] - 1) * (deta[i][ss] - ((double)reach[ss] - 1));  //以ss为源点，其他顶点的reach值之和等于deta-reach[ss]+1

	}
	stack<int>().swap(S);
	vector<vector<pair<int, int>>>().swap(pred);
	vector<vector<double>>().swap(pathnum);
	vector<vector<double>>().swap(deta);
	/*for (int i = 0; i < sidesetsize; i++) {
		cout << "源点" << samesideset[i] << endl;
		for (int j = 0; j < Pmg.n_vertex; j++) {
			if (j == 17) {
				cout << "对点" << j << "源依赖" << deta[i][j] << endl;
			}
		}
	}*/
}

void todoType2Ident() {
	for (int i = 0; i < Pmg.n_vertex; i++) {  //从index0——n遍历，所以对index=k的点，不用再去判断是否与0--k-1的点是否identical，以为在之前已经判断过了，若是，则flag[k]应该是1了
		if (flag_side[i] == 0 && flag_addtoident[i] == 0) {  //side顶点不去寻找identical；已经找到proxy的顶点不用在寻找identical
			//因为当识别到一个i的ident后，就在Pmg.adjlist[i]中删除，此时x会跳过一些点
			map<int,int> iadjlist = Pmg.adjlist[i];
			for (auto v : iadjlist) {
				if (v.first > i&&flag_addtoident[v.first]==0) {                //判断>i的顶点v是否与i是identical的（<i的v已经在对v的邻居点遍历时判断过了，与i都不相同）
					int flag_isident = judgeIdent(i, v.first, 2);      //判断v是否与i是type2 ident的
					if (flag_isident) {
						//add dependency and modify Mp、adjlist and update ident set、ident value of i
						flag_addtoident[v.first] = 1;
						count_ident2++;
						mergeIdent2(i, v.first);
					}
				}
			}
			map<int, int>().swap(iadjlist); 
		}
	}
}

void todoType1Ident() {
	for (int i = 0; i < Pmg.n_vertex; i++) {
		if (i % 100 == 0) cout << i << " ";
		if (flag_side[i] == 0 ) {
			set<int> i2adjlist;
			for (auto neibor : Pmg.adjlist[i]) {
				for (auto v : Pmg.adjlist[neibor.first]) {  //i的二跳邻居
					if (v.first > i && Pmg.adjlist[i].find(v.first)==Pmg.adjlist[i].end()&&flag_addtoident[v.first]==0) {  //非邻居，且未加入到任何顶点的identset中。（有些顶点v可能会代表一些与其identical type2的顶点）
						i2adjlist.insert(v.first);
					}
				}
			}
			for (int v : i2adjlist) {
				int flag_isident = judgeIdent(i, v, 1);
				if (flag_isident) {
					flag_addtoident[v] = 1;
					mergeIdent1(i, v);
					count_ident1++;
					cout << v << "ident1" << i << endl;
				}
				
			}
			set<int>().swap(i2adjlist);
		}
	}
}

int judgeIdent(int i, int v, int type) {
	if (Pmg.adjlist[i].size() != Pmg.adjlist[v].size()) return 0;

	/*method1*/
	auto iteriadj = Pmg.adjlist[i].begin();
	auto itervadj = Pmg.adjlist[v].begin();
	if (type == 2) {
		while (iteriadj != Pmg.adjlist[i].end()&&itervadj!=Pmg.adjlist[v].end()) {
			if (iteriadj->first == v) iteriadj++;
			if (itervadj->first == i) itervadj++;
			if (iteriadj == Pmg.adjlist[i].end() || itervadj == Pmg.adjlist[v].end()) break;
			if (iteriadj->first != itervadj->first) return 0;
			if (iteriadj->second != itervadj->second) return 0;
			iteriadj++; itervadj++;
		}
	}
	else {
		while (iteriadj != Pmg.adjlist[i].end() && itervadj != Pmg.adjlist[v].end()) {
			if (iteriadj->first != itervadj->first) return 0;
			if (iteriadj->second != itervadj->second) return 0;
			iteriadj++; itervadj++;
		}
	}
	return 1;
	
	/*method2*/
	//set<pair<int, int>> i_adj;
	//set<pair<int, int>> v_adj;
	//if (type == 1) {
	//	for (auto ii : Pmg.adjlist[i]) {
	//		i_adj.insert({ ii.first,ii.second });
	//	}
	//	for (auto vv : Pmg.adjlist[v]) {
	//		v_adj.insert({ vv.first,vv.second });
	//	}
	//}
	//if (type == 2) {  //type2类型
	//	for (auto ii : Pmg.adjlist[i]) {
	//		if (ii.first == v) continue;
	//		i_adj.insert({ ii.first,ii.second });
	//	}
	//	for (auto vv : Pmg.adjlist[v]) {
	//		if (vv.first == i) continue;
	//		v_adj.insert({ vv.first,vv.second });
	//	}
	//}
	//set<pair<int,int>> result;  //判断邻居、权重是否相同
	//set_intersection(i_adj.begin(), i_adj.end(), v_adj.begin(), v_adj.end(), inserter(result, result.begin()));
	//if (result.size() != i_adj.size()) return 0;
	//return 1;
}

void mergeIdent1(int i, int v) {  //注意i、v有reach值+！！且有ident值！！
	//add pair dependency
	//对i、v的pair dependency
	int reachvident = reach[v], reachiident = reach[i];
	if (ident[v] > 0) {
		for (int vident : IdentSet[v])
			reachvident += reach[vident];
	}
	if (ident[i] > 0) {
		for (int iident : IdentSet[i])
			reachiident += reach[iident];
	}
	bc[v] += ((double)reach[v] - 1) * reachiident;
	if (ident[v] > 0) {
		for (int vident : IdentSet[v]) {
			bc[vident] += ((double)reach[vident] - 1) * reachiident;
		}
	}
	bc[i] += ((double)reach[i] - 1) * reachvident;
	if (ident[i] > 0) {
		for (int iident : IdentSet[i]) {
			bc[iident] += ((double)reach[iident] - 1) * reachvident;
		}
	}
	//对邻居顶点的pair dependency
	int numpathu_i = 0;
	map<int,int> neighborlist = Pmg.adjlist[i];
	int neighbornum = neighborlist.size();
	/*vector<int> weighti, weightv;*/  //记录i、v与邻居之间的权重:i\v与任意邻居w的权重相等
	/*for (auto w : neighborlist) weighti.push_back(w.second); */   //???????????????????????????????????????????
	/*for (int w : neighborlist) weightv.push_back(Pmg.Mp[v][w]);*/
	for (auto w:neighborlist) {
		numpathu_i += w.second*w.second;
	}
	for (auto w:neighborlist) {
		bc[w.first] += reachiident * ((double)w.second * w.second / numpathu_i) * reachvident;  //因为是对称的，所以numpath_w_v=numpath_i_w=weightx[wk]
		bc[w.first] += reachvident * ((double)w.second * w.second / numpathu_i) * reachiident;  //??????????????????
	}
	//update IdentSet and ident
	IdentSet[i].push_back(v);
	ident[i]+=1;
	if (ident[v] > 0) {
		for (int identv : IdentSet[v]) {
			IdentSet[i].push_back(identv);
			ident[i] += 1;
		}
	}
	//update adjlist and Mp
	count_ident_edge += Pmg.adjlist[v].size();
	Pmg.deleteVertex(v);
}

void mergeIdent2(int i, int v) {  //注意i、v有reach值
	//add pair dependency v与i+identSet[i]
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

	//update IdentSet and ident
	IdentSet[i].push_back(v);
	ident[i]++;
	//update adjlist and Mp
	count_ident_edge += Pmg.adjlist[v].size();
	Pmg.deleteVertex(v);
}

/*以s为源点进行宽搜+反向累积，计算源点s的源依赖*/
void BBC1(int s) {  //+reach值后
	//initialize variable

	queue<int> Q;
	stack<int> S;
	vector<vector<pair<int,int>>> pred(Pmg.n_vertex);
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
		for (auto w : Pmg.adjlist[v]) {
			if (dist[w.first] == DBL_MAX) {
				dist[w.first] = dist[v] + 1;
				Q.push(w.first);
			}
			if (dist[w.first] == dist[v] + 1) {
				pathnum[w.first] += pathnum[v] * w.second;
				pred[w.first].push_back({ v,w.second });
			}
		}
	}
	/*for (int i = 0; i < pathnum.size(); i++) {
		cout << "到点" << i << "最短路径数：" << pathnum[i] << endl;
	}*/
	pmid = clock();
	//back propagation
	vector<double> deta(Pmg.n_vertex);
	for (int i = 0; i < Pmg.n_vertex; i++) deta[i] = (double)reach[i] - 1;
	while (!S.empty()) {
		int w = S.top();
		S.pop();
		for (auto v : pred[w]) {
			deta[v.first] += (pathnum[v.first] * v.second / pathnum[w]) * (1 + deta[w]);
		}
		if (w != s) {
			bc[w] += deta[w] * reach[s];

			/*if (w == 17) {
				cout << "对点" << w << "源依赖:" << deta[w] * reach[s] << endl;
			}*/
		}
	}

}

void BBC2(int s) {  //+reach值+ident后
	//initialize variable
	/*cout << "源点" << s << endl;*/
	queue<int> Q;
	stack<int> S;
	vector<vector<pair<int,int>>> pred(Pmg.n_vertex);
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
		for (auto w : Pmg.adjlist[v]) {
			if (dist[w.first] == DBL_MAX) {
				dist[w.first] = dist[v] + 1;
				Q.push(w.first);
			}
			if (dist[w.first] == dist[v] + 1) {
				pathnum[w.first] += pathnum[v] * w.second * (1 + (double)ident[v]);
				pred[w.first].push_back({ v,w.second });
			}
		}
	}
	/*for (int i = 0; i < pathnum.size(); i++) {
		cout << "到点" << i << "最短路径数：" << pathnum[i] << endl;
	}*/
	queue<int>().swap(Q);
	vector<double>().swap(dist);
	pmid = clock();
	//back propagation
	int reachs_total = reach[s];
	if (ident[s] > 0) {
		for (int sident : IdentSet[s]) {
			reachs_total += reach[sident];
		}
	}
	vector<double> deta(Pmg.n_vertex);
	for (int i = 0; i < Pmg.n_vertex; i++) deta[i] = (double)reach[i] - 1;
	while (!S.empty()) {
		int w = S.top();
		S.pop();
		for (auto v : pred[w]) {
			//当w代表其他ident顶点时，w和w代表的所有顶点对前置顶点的贡献
			double deta_add = (pathnum[v.first] * v.second / pathnum[w]) * (1 + deta[w]);
			if (ident[w] > 0) {
				for (int wident : IdentSet[w]) {
					deta_add += (pathnum[v.first] * v.second / pathnum[w]) * (1 + deta[wident]);  //deta[w]与deta[wident]不同，因为w与wident可能有不同的reach值
				}
			}
			deta[v.first] += deta_add;
			if (ident[v.first] > 0) {
				for (int vident : IdentSet[v.first]) {
					deta[vident] += deta_add;
				}
			}
		}
		if (w != s) {
			
			bc[w] += deta[w] * reachs_total;
			if (ident[w] > 0) {
				for (int wident : IdentSet[w]) {
					bc[wident] += deta[wident] * reachs_total;
				}
			}

			/*if (w == 17) {
				cout << "对点" << w << "源依赖:" << deta[w] * reach[s] << endl;
			}*/
		}
	}
	stack<int>().swap(S);
	vector<vector<pair<int, int>>>().swap(pred);
	vector<double>().swap(pathnum);
	vector<double>().swap(deta);

}

void find_side0(map<int, map<int,int>>* Mpl) {  //Mpl的邻接表，若第i的顶点的邻接表序列size为1，则该点为side
	flag_side.resize(Pmg.n_vertex);
	auto iteri = Mpl->begin();
	for (iteri; iteri != Mpl->end(); iteri++) {
		if (iteri->second.size() == 1) {  //如果i的邻接点只有一个
			int i = iteri->first;
			auto iterj = iteri->second.begin();
			int irelatedj = iterj->first;
			SideSet0[irelatedj].push_back(i);
			flag_side[i] = 1;
		}
	}
	//for (int i = 0, lenr = Mpl->size(); i < lenr; i++) {  //A1类顶点i
	//	
	//	int counti = 0;
	//	int irelatedj;
	//	int j, lenc;
	//	for (j = 0, lenc = (*Mpl)[0].size(); j < lenc; j++) {  //AL类顶点j
	//		if ((*Mpl)[i][j] > 0)
	//		{
	//			counti++;
	//			irelatedj = j;
	//		}
	//	}
	//	if (counti == 1) {
	//		SideSet0[irelatedj].push_back(i);
	//		flag_side[i] = 1;
	//	}
	//}
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
	for (int i = 0; i < P->size() / 2; i++) {  //一共遍历P->size()/2条边
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

	/*遍历Mpl，找定义②的side顶点*/
	double ps13, pt13, time13;
	ps13 = clock();
	find_side0(&matrixpl);
	pt13 = clock();
	time13 = (double)(pt13 - ps13) / CLOCKS_PER_SEC;
	cout << "识别side顶点集，用时：" << time13 << endl;
	out_result << "识别side顶点集，用时：" << time13 << endl;

	/*计算Mp+邻接表*/
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
	double pt14 = clock();
	double time14 = (double)(pt14 - ps12) / CLOCKS_PER_SEC;
	out_result << "计算mpl的转置用时：" << time14 << endl;

	for (int i = 0; i < n_vertex; i++) {  //计算第i个点的邻接序列
		/*cout << "计算到第" << i;*/
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
	out_result << "由mpl得到mp用时：" << time12 << endl;
	/*save_Pmg();*/

}

/*在P-multigraph中删除边*/
void multigraph::deleteEdge(int a, int b) {
	map<int,int>::iterator it = adjlist[a].begin();
	for (it; it != adjlist[a].end(); it++) {
		if (it->first == b) {
			adjlist[a].erase(it);
			break;
		}
	}

	map<int,int>::iterator itt = adjlist[b].begin();
	for (itt; itt != adjlist[b].end(); itt++) {
		if (itt->first == a) {
			adjlist[b].erase(itt);
			break;
		}
	}
}
/*在P-multigraph中删除点*/
void multigraph::deleteVertex(int a) {
	//对a的每个邻接点
	for (auto an : adjlist[a]) {
		//修改邻接矩阵
		/*Mp[a][an] = 0; Mp[an][a] = 0;*/
		//在an的邻接表中删除a
		map<int, int>::iterator it = adjlist[an.first].begin();
		for (it; it != adjlist[an.first].end(); it++) {
			if (it->first == a) {
				adjlist[an.first].erase(it);
				break;
			}
		}
	}
	//清空a的邻接表（此时a为P-multigraph中一个独立的点，但通过index遍历仍然能遍历到a）
	map<int, int>().swap(adjlist[a]);
}