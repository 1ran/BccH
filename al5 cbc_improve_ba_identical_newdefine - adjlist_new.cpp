//ʶ��identical���㣨���¶��壬ͬһidentical�㼯�У�����identical��֮���edgeinfo��ͬ����Pmg����adjlist
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
	int n_vertexclass;  //���������
	map<char, int> vertex_class_index;  //�������-������
	vector<int> vertex_num;  //����������˳��洢ÿ�ඥ����

	int n_edgeclass;  //�ߵ������
	vector<vector<int>> edge_class;  //ÿ������ӵĶ�����𣺶�����x��������y
	vector<int> edge_num;  //ÿ�����

	//vector<vector<vector<int>>> edge_info;  //3ά���飬��һάedge_num����ʾÿ��ߣ��ڶ�ά�͵���ά��ʾ��edge���ڽӾ���
	vector<map<int, map<int, int>>> edge_info;  //3ά���飬��һάedge_num����ʾÿ��ߣ��ڶ�ά�͵���ά��ʾ��edge���ڽӱ�

	void read_file();
	void show_hg();
	int getEdgeIndex(int x, int y);  //���������Ԫ·��P��ƥ���칹schemaͼ�ϵı�
};


class multigraph {
public:
	int n_vertex;  //���¶�����������copy�Ĺؽڶ��㣩
	int n_vertex_org;  //ԭ������
	int m_edge;  //����
	vector<vector<vector<int>>> edge_info;  //�洢ÿ���ߵ���Ϣ[[class1index,num1],[class2index,num2]...]class��·��ʵ����𣬼�·��ʵ�����ĸ�E�ඥ�����
	//vector<vector<int>> Mp_edge;  //�洢ai��aj֮��ߵ�����??
	vector<map<int, int> > adjlist_edge;  //��ÿһ�д洢��������ڶ����ͬʱ���洢ͨ�������������ڶ�������
	//vector<vector<int> > Mp;  //??
	//vector<vector<int> > adjlist;

	void show_mg();
	void getPmg(vector<char>* P);  //���칹ͼ�л��Pmg
	void getPmg_read_file();      //�ӱ�����ļ��ж�ȡMpl��Pmg
	void save_Mpl(vector<vector<int>>* Mpl);  //����Mpl
	void save_Pmg();  //��������Pmg
	void deleteVertex(int v);
	void deleteEdge(int a, int b);
};

void CBC(int s);  //����v���·��ʵ����Ϣ
void CBC_(int s);  //ֱ����pred�б����Ȩ��pred��vector  ����졿
void CBC_1(int s);  //ֱ����pred�б����Ȩ��pred��map 
void CBC_ident0(int s);  //+ident��ident[s]=0
void CBC_ident1(int s);  //+ident��ident[s]>0
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
map<int, vector<int>> IdentSet;  //�����������identical������б�
map<int, map<int, int>> IdentSetEdgeindex;  //��¼proxy��������������Ķ���֮��ߵ�����
vector<int> ident;     //�����������identical����ĸ������������Լ���
vector<int> flag_addtoident;  //�Ѿ�����identical�㼯�Ķ���
//vector<double> deta0;  //��֤
vector<int> reach;
map<int, int> org;  //map[100]=1  ��ʾ����100�ڹؽڶ��㻮��ǰ��ԭ������1
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
	//һ�������칹ͼ����Pmg
	double ps00, pt00, time00;
	ps00 = clock();
	hg.read_file();
	pt00 = clock();
	time00 = (double)(pt00 - ps00) / CLOCKS_PER_SEC;
	out_result << "�����칹ͼ��ʱ��" << time00 << endl;
	cout << "�����칹ͼ��ʱ��" << time00 << endl;
	/*hg.show_hg();*/
	double ps0, pt0, time0;
	ps0 = clock();
	vector<char> P = { 'A','M','W','M','A' };
	//vector<char> P = { 'A','P','V','P','A' };
	//vector<char> P = { 'B','U','B' };
	Pmg.getPmg(&P);  //�õ�P-multigraph
	pt0 = clock();
	time0 = (double)(pt0 - ps0) / CLOCKS_PER_SEC;
	cout << "�õ�Pmg��ʱ��" << time0 << endl;
	out_result << "�õ�Pmg��ʱ:" << time0 << endl;
	cout << "������" << Pmg.n_vertex << endl;
	out_result << "Pmg������" << Pmg.n_vertex << endl;
	cout << "������" << Pmg.m_edge << endl;
	out_result << "Pmg������" << Pmg.m_edge << endl;
	//Pmg.show_mg();
	//�����Ѿ�������Pmg�����ļ��ж�ȡ
	/*Pmg.getPmg_read_file();*/
	/*Pmg.show_mg();*/

	/*ͳ�������ͨ������С*/
	double pcmps, pcmpt, timecmp;
	pcmps = clock();
	int maxComponent = getMaxComponent();
	cout << "�����ͨ������" << maxComponent << endl;
	out_result << "Pmg�����ͨ������" << maxComponent << endl;
	pcmpt = clock();
	timecmp = (double)(pcmpt - pcmps) / CLOCKS_PER_SEC;
	out_result << "���������ͨ������С��ʱ��" << timecmp << endl;


	bc.resize(Pmg.n_vertex);
	reach.resize(Pmg.n_vertex,1);
	/*deta0.resize(Pmg.n_vertex);*/
	
	double ps1, pt1, time1;
	ps1 = clock();
	BridgeEdgeDivide();  //�ű߷ָ�
	pt1 = clock();
	time1 = (double)(pt1 - ps1) / CLOCKS_PER_SEC;
	cout << "�ű߷ָ����,������" << countb << "���ű�,����ʱ��" << time1 << endl;
	out_result << "�ű߷ָ����,������(ɾ��)" << countb << "���ű�,����ʱ��" << time1 << endl;

	double pcmpsb, pcmptb, timecmpb;
	pcmpsb = clock();
	int maxComponentb = getMaxComponent();
	cout << "�ű߷ָ�������ͨ������" << maxComponentb << endl;
	out_result << "�ű߷ָ�������ͨ������" << maxComponentb << endl;
	pcmptb = clock();
	timecmpb = (double)(pcmptb - pcmpsb) / CLOCKS_PER_SEC;
	out_result << "���������ͨ������С��ʱ��" << timecmpb << endl;


	double ps2, pt2, time2;
	ps2 = clock();
	ArticulationVertexDivide();  //�ؽڶ���ָ�
	bc.resize(Pmg.n_vertex);  //�¼����copy�ؽڶ����bcֵ�Ĵ洢�ռ�
	pt2 = clock();
	time2 = (double)(pt2 - ps2) / CLOCKS_PER_SEC;
	cout << "�ؽڶ���ָ���ϣ�������" << counta << "���ؽڶ���,����ʱ��" << time2 << endl;
	out_result << "�ؽڶ���ָ���ϣ�������" << counta << "���ؽڶ���,����ʱ��" << time2 << endl;

	double pcmps1, pcmpt1, timecmp1;
	pcmps1 = clock();
	int maxComponent1 = getMaxComponent();
	cout << "�űߡ��ؽڶ���ָ�������ͨ������" << maxComponent1 << endl;
	out_result << "�űߡ��ؽڶ���ָ�������ͨ������" << maxComponent1 << endl;
	pcmpt1 = clock();
	timecmp1 = (double)(pcmpt1 - pcmps1) / CLOCKS_PER_SEC;
	out_result << "���������ͨ������С��ʱ��" << timecmp1 << endl;


	/*for (int i = 0; i < Pmg.n_vertex; i++) cout << i << ":" << bc[i] << endl;*/
	double ps4, pt4, time4;
	ps4 = clock();
	Ident();
	pt4 = clock();
	time4 = (double)(pt4 - ps4) / CLOCKS_PER_SEC;
	cout << "identical����ʶ����ϲ���ʱ��" << time4 << endl;
	out_result << "identical����ʶ����ϲ���ʱ��" << time4 << endl;

	double pcmps2, pcmpt2, timecmp2;
	pcmps2 = clock();
	int maxComponent2 = getMaxComponent();
	cout << "identical����ϲ���Ϻ������ͨ������" << maxComponent2 << endl;
	out_result << "identical����ϲ���Ϻ������ͨ������" << maxComponent2 << endl;
	pcmpt2 = clock();
	timecmp2 = (double)(pcmpt2 - pcmps2) / CLOCKS_PER_SEC;
	out_result << "���������ͨ������С��ʱ��" << timecmp2 << endl;


	double ps5, pt5, time5 = 0;
	double time5bfs = 0, time5back = 0;
	int count = 0;
	for (int i = 0; i < Pmg.n_vertex; i++) {
		if (flag_addtoident[i] == 0) {
			ps5 = clock();
			if (ident[i] == 0) {
				CBC_ident0(i);  //��Դ��i����identical�㼯�е�proxy����
			}
			else {
				CBC_ident1(i);  //��Դ��i��proxy����
			}
			pt5 = clock();
			time5 += (double)(pt5 - ps5) / CLOCKS_PER_SEC;
			time5bfs += (double)(pmid - ps5) / CLOCKS_PER_SEC;
			time5back += (double)(pt5 - pmid) / CLOCKS_PER_SEC;

		}
		count++;
		if(count%100==0)
		    cout << count<<"������ΪԴ�����ʱ�䣺" << time5<< endl;
	}

	//���ݹؽڶ����ӳ���ϵ������ؽڶ����bcֵ
	for (auto it : org) {
		bc[it.second] += bc[it.first];
	}
	double ptt = clock();

	cout << "ǰ��bfs��ʱ�䣺" << time5bfs << endl;
	cout << "�����ۻ���ʱ�䣺" << time5back << endl;
	cout << "����bcֵ��ʱ�䣺" << time5 << endl;
	out_result << "ǰ��bfs��ʱ�䣺" << time5bfs << endl;
	out_result << "�����ۻ���ʱ�䣺" << time5back << endl;
	out_result << "����bcֵ��ʱ�䣺" << time5 << endl;

	double time_all = (double)(ptt - ps00) / CLOCKS_PER_SEC;
	cout << "����ʱ��" << time_all << endl;
	out_result << "����ʱ��" << time_all << endl;

	/*for (int i = 0; i < Pmg.n_vertex; i++) {
		if (bc[i] > 0) {
			cout << "����" << i << "��bcֵ��" << bc[i] << endl;
		}
	}*/
	//compare_result();

	cout << "Ȩ�����" << weight_max << endl;
	cout << "Ȩ����С��" << weight_min << endl;
	cout << "Ȩ��D���" << weight_D_max << endl;
	cout << "Ȩ��D��С��" << weight_D_min << endl;
	cout << "Ȩ��pathinstance���" << weight_pathinstance_max << endl;
	cout << "Ȩ��pathinstance��С��" << weight_pathinstance_min << endl;
	out_result << "Ȩ�����" << fixed << setprecision(8) << weight_max << endl;
	out_result << "Ȩ����С��" << fixed << setprecision(8) << weight_min << endl;
	out_result << "Ȩ��D���" << fixed << setprecision(8) << weight_D_max << endl;
	out_result << "Ȩ��D��С��" << fixed << setprecision(8) << weight_D_min << endl;
	out_result << "Ȩ��pathinstance���" << fixed << setprecision(8) << weight_pathinstance_max << endl;
	out_result << "Ȩ��pathinstance��С��" << fixed << setprecision(8) << weight_pathinstance_min << endl;

	out_result.close();
	save_bc();
	return 0;
	
	//string filename = file + "/deta0_.txt";  //base�д洢�ļ��Ļ�����Ϣ�����ݱ������ж�Ӧ�ıߵ��������ҵ��ߵ��ļ�"x.txt"
	//ifstream inputt(filename, ios::in);  //�����ļ�������input
	////�ж��ļ��Ƿ����
	//if (!inputt) {
	//	cerr << "file error!" << endl;
	//	exit(1);
	//}
	//vector<double> deta0_(Pmg.n_vertex_org);
	////�������ж����deta0ֵ
	//string deta0_info;
	//while (getline(inputt, deta0_info)) {
	//	istringstream sdeta0_info(deta0_info);
	//	int v;
	//	double deta0;
	//	sdeta0_info >> v >> deta0;
	//	deta0_[v] = deta0;
	//}
	//�Ƚ�deta0ֵ
	/*vector<int> different_v;
	for (int i = 0; i < Pmg.n_vertex_org; i++) {
		if (abs(deta0[i] - deta0_[i]) > 0.000001) different_v.push_back(i);
	}
	cout << "deta0ֵ��ͬ�Ķ�������" << different_v.size() << endl;
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
	string filename = file + "/result.txt";  //base�д洢�ļ��Ļ�����Ϣ�����ݱ������ж�Ӧ�ıߵ��������ҵ��ߵ��ļ�"x.txt"
	ifstream input(filename, ios::in);  //�����ļ�������input
	//�ж��ļ��Ƿ����
	if (!input) {
		cerr << "file error!" << endl;
		exit(1);
	}

	vector<double> bc0(Pmg.n_vertex_org);
	//�������ж����bcֵ//����ÿ�������bcֵ���ļ�
	string bc_info;
	while (getline(input, bc_info)) {
		istringstream sbc_info(bc_info);
		int v;
		double bc;
		sbc_info >> v >> bc;
		bc0[v] = bc;
	}

	//�Ƚ�bcֵ
	vector<int> different_v;
	for (int i = 0; i < Pmg.n_vertex_org; i++) {
		if (abs(bc[i] - bc0[i]) > 0.000001) different_v.push_back(i);
	}
	cout << "bcֵ��ͬ�Ķ�������" << different_v.size() << endl;
	for (int v : different_v) {

		cout << v << ":" << bc[v] << " " << bc0[v] << endl;
	}
}

/*�����űߡ��ؽڶ���*/
/*���㶥��i���ڵ���ͨ���������ж����reachֵ֮��*/
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

/*�����ű�*/
void BridgeEdgeDivide() {
	dfn.resize(Pmg.n_vertex);  //dfn[i]��¼�˵���b_tarjan��������dfsʱ��˳��
	low.resize(Pmg.n_vertex);  //low[i]��¼�˵�i���ھ��У���С��dfnֵ
	for (int i = 0; i < Pmg.n_vertex; i++) {
		stamp = 0;
		int totalreach;
		if (!dfn[i]) {
			totalreach = countTotalReach(i);  //����i���ڵ���ͨ�����У����������������ԭͼ���Ƕ������ͨ���������
			/*cout << "totalreach:" << totalreach << endl;*/
			/*outt << "��������" << i << "ΪԴ��,total reach:" << totalreach << endl;*/
			b_tarjan(i, -1, totalreach);  //��i���ڵ���ͨ����ִ��b_tarjan�������ű�
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
		/*outt  << index << "���ھ�" << j << endl;*/
		/*int j = j_[0];*/
		if (j == fa) continue;
		if (!dfn[j]) {
			b_tarjan(j, index, totalreach);
			/*outt << "��" << j << "�������" << endl;
			outt << "����low[" << index << "]=min(" << low[index] << "," << low[j] << ")";*/
			low[index] = min(low[index], low[j]);
			/*outt << "low[" << index << "]=" << low[index] << endl;*/
			/*outt << "�ж�" << index << "-" << j << "�Ƿ����ű�" << endl;*/
			if (low[j] > dfn[index]) {
				//��index-j���ű�
				cout << "�ҵ��űߣ�" << index << "-" << j << endl;
				countb++;
				/*outt << index << "-" << j << endl;*/
				int repj = 0, repindex = 0;
				/*cout << "ջ��" << endl;*/
				while (stack_bi.top() != j) {  //����j���ڵķ����У������reach֮��
					int temp = stack_bi.top();
					/*cout << temp << endl;*/
					stack_bi.pop();
					repj += reach[temp];
				}
				/*cout << j << endl;*/
				repj += reach[j]; stack_bi.pop();
				repindex = totalreach - repj;  //index���ڷ����Ķ����reach֮�ͣ�������reach-j���ڷ�����reach֮��
				//����index��j��reachֵ
				reach[index] += repj;  reach[j] += repindex;

				//�ۼ�ȱ�ٵ�pair-dependency
				bc[index] += ((double)repindex - 1) * repj;
				bc[j] += ((double)repj - 1) * repindex;
				//�޸��ڽӱ��ڽӾ���
				Pmg.deleteEdge(index, j);
			}
		}
		else if (dfn[j] < dfn[index]) {
			/*outt << j<<"��֮ǰ�����ʹ���ͨ��dfn[" << j << "]����low[" << index << "]=(min(" << low[index] << ", " << dfn[j] << "))  ";*/
			low[index] = min(low[index], dfn[j]);
			/*outt << "low[" << index << "]=" << low[index] << endl;*/
		}
	}
	map<int, int>().swap(indexadj);
}


//ʶ��ؽڶ��㣬ע��Ҫ��bc��������������vi��bcֵ�����ڽӾ�������������vi���к��У����ڽӱ�����vi���沿��v������Pmg��n_vertex
/*���ֹؽڶ���*/
void ArticulationVertexDivide() {
	dfn.resize(Pmg.n_vertex);
	low.resize(Pmg.n_vertex);
	for (int i = 0; i < Pmg.n_vertex_org; i++) {
		stamp = 0;
		int totalreach;
		if (!dfn[i]) {
			totalreach = countTotalReach(i);  //����i���ڵ���ͨ�����У����������
			/*cout << "totalreach:" << totalreach << endl;*/
			a_tarjan(i, -1, totalreach);  //��i���ڵ���ͨ����ִ��b_tarjan�������ű�
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
	map<int, int> indexadj = Pmg.adjlist_edge[index];  //��Ϊ�·���artvertexcopy�У����ܻ��޸�Pmg.adjlist[index]�������Ƚ��䱣�����������޸�ǰ��index���ھ��ϱ���
	for (auto iteradjj = indexadj.begin(); iteradjj != indexadj.end(); iteradjj++) {
		int j = iteradjj->first;
		if (j == fa) continue;
		if (!dfn[j]) {
			child++;
			a_tarjan(j, index, totalreach);
			low[index] = min(low[index], low[j]);
			if (low[j] >= dfn[index] && fa >= 0) {
				cout << "�ҵ��ؽڶ��㣺" << index << endl;
				counta++;
				artvertexcopy(index, j, totalreach);  //��indexΪ��㣬�ָ�j���ڵ���ͨ��������index-->copyindex
			}
			//��indexΪ���ڵ�ʱ������index���ھӣ��ж��Ƿ񻹻�����һ�����ӣ���û�У���index���Ǹ�㣬���У�����index-j���ڷ����ĸ��
			if (low[j] >= dfn[index] && fa < 0) {
				/*cout << j << " " << index << endl;*/
				auto nowj = iteradjj;
				nowj++;
				int isbi = 0;
				for (nowj; nowj != indexadj.end(); nowj++) {  //�ж�index���ھ��Ƿ񶼷��ʵ��ˣ����ж�index�Ƿ��������ھӽڵ�
					int nb = nowj->first;
					if (!dfn[nb]) {  //��δ���ʵ����ھӣ�˵��index�Ǹ�㣬��ǰջ�е����ж��㹹��һ����ͨ����
						cout << "�ҵ��ؽڶ���(����Դ��)��" << index << endl;
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
	//���ڵ�Ϊ���ʱ
	//if (fa < 0 && child>1) {  //�Ϸ������˵�index��Ϊ���ڵ�ʱ��indexΪ������ͨ����������indexΪ���ڵ㣬����ͨ��low[j]>=dfn[index]���ж�index�Ƿ�Ϊ���
	//	
	//}
}


/*tarjan�㷨�ҵ��ؽڶ���󣬸��ƹؽڶ��㲢�����ڽӾ����ڽӱ�reachֵ*/
void artvertexcopy(int index, int j, int totalreach) {
	//����index�Ǹ��
	int copyindex = Pmg.n_vertex++;  //����ͨ�����У�copyһ��index��
	cout << "���ƹؽڵ�" << copyindex << endl;
	org[copyindex] = index;//����copyindex��ԭ����index
	//�����ڽӾ��󣬲����µĶ���
	//Pmg.Mp.resize(Pmg.n_vertex);
	//Pmg.Mp_edge.resize(Pmg.n_vertex);
	//for (int i = 0; i < Pmg.n_vertex; i++) {
	//	Pmg.Mp[i].resize(Pmg.n_vertex);
	//	Pmg.Mp_edge[i].resize(Pmg.n_vertex);
	//}
	Pmg.adjlist_edge.resize(Pmg.n_vertex);
	//Pmg.adjlist.resize(Pmg.n_vertex);
	int repcopyindex = 0, repindex = 0;
	while (stack_bi.top() != j) {  //��index-j�����ͨ�����еĶ����ջ
		int temp = stack_bi.top();
		stack_bi.pop();
		//����j����Ķ�����
		repcopyindex += reach[temp];
		/*�����ڽӾ�����ڽӱ�����ͨ�����еĶ�����index�ڽӣ��޸�Ϊ��copyindex�ڽ�*/
		auto itertemp = Pmg.adjlist_edge[index].find(temp);
		if (itertemp != Pmg.adjlist_edge[index].end()) {  //���temp������Pmultigraph����index�������������Ϊ��copyindex����
			//����index-temp���ڽӾ���
			int e_index = itertemp->second;  //???????????????
			//int e_index = Pmg.adjlist_edge[index][temp];

			//��Pmg��ɾ����index-temp
			Pmg.deleteEdge(index, temp);
			Pmg.adjlist_edge[copyindex].insert(map<int, int>::value_type(temp, e_index));
			Pmg.adjlist_edge[temp].insert(map<int, int>::value_type(copyindex, e_index));
			//Pmg.adjlist[copyindex].push_back(temp);
			//Pmg.adjlist[temp].push_back(copyindex);
		}
	}
	//���j��ջ������j�Ĵ���Ķ��������ڽӾ�����ڽӱ�
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

	repindex = totalreach - repcopyindex;  //repcopyindex��j������ͨ�����ϣ������reachֵ֮�ͣ�������index�����reach����repindex�Ƿָ���ͨ������index���ڷ�����reachֵ֮�ͣ�����index��reach��
	reach[index] += repcopyindex;  //���¶���index��reachֵ���¼���j���ڵķ����϶����reachֵ֮�ͣ���ͬ��index��reach����Ϊ����������
	reach.push_back(repindex);  //copyindex���¼���Ķ��㣬����Ҫ��reach��push_back��copyindex��reachֵʱindex������ͨ�����ж����reach֮��(����index�����reach)
}





void Ident() {
	ident.resize(Pmg.n_vertex);
	flag_addtoident.resize(Pmg.n_vertex);

	todoType2Ident();
	//ͳ��identical�㼯�Ͷ���
	int countset = 0, countvertex = 0;
	for (auto idset : IdentSet) {
		countset++;
		for (int i : idset.second) {
			countvertex++;
		}
	}
	cout << "�ҵ�identical�㼯" << countset << "��������identical�㼯�ĵ�(������proxy)" << countvertex << "��" << endl;
	out_result << "�ҵ�identical�㼯" << countset << "��������identical�㼯�ĵ�(������proxy)" << countvertex << "��" << endl;
	cout << "�ϲ�identical�㣬ɾ��������" << count_ident_edge << endl;
	out_result << "�ϲ�identical�㣬ɾ��������" << count_ident_edge << endl;

	//identical�㼯���Ƿ������������֮��ı߲���ͬ
	//int countdiff = 0;
	//vector<int> different;
	//for (auto idset : IdentSet) {  //�ж�iset��ÿ������֮���edgeinfo�Ƿ����
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
//	for (int i = 0; i < Pmg.n_vertex; i++) {  //����i=0-n_vertex��Ѱ����identical�㼯
//		if (flag_addtoident[i] == 0) {
//			vector<int> mayidentset = { i };
//			//vector<int> iadjlist = Pmg.adjlist_edge[i];
//			//map<int, int> iadjlist = Pmg.adjlist_edge[i];  //����Ҫ�ȵ����洢����Ϊ���Ƚ�i��ident��v����mayidentset�У����޸�adjli_edge[i]
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
	for (int i = 0; i < Pmg.n_vertex; i++) {  //����i=0-n_vertex��Ѱ����identical�㼯
		if (flag_addtoident[i] == 0) {
			vector<int> iwithidentset = { i };
			//vector<int> iadjlist = Pmg.adjlist_edge[i];
			//map<int, int> iadjlist = Pmg.adjlist_edge[i];  //����Ҫ�ȵ����洢����Ϊ���Ƚ�i��ident��v����mayidentset�У����޸�adjli_edge[i]
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
	//������reachֵʱ
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

	//��¼i��v֮��ߵ�����
	int eindexi_v = Pmg.adjlist_edge[i][v];
	if (ident[i] == 0) {  //˵��v�ǵ�һ���ӵ�iΪproxy��identical set�ĵ�,����IdentSetEdgeindex��Ŀǰû��i
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
		if (iteriadj->first == v) { //��i���ھ�������v
			iteriadj++; 
		}  
		if (itervadj->first == i) { //��v���ھ�������i
			itervadj++;
		} 
		if (iteriadj != Pmg.adjlist_edge[i].end()&&itervadj != Pmg.adjlist_edge[v].end()) {
			if (iteriadj->first != itervadj->first) return 0;  
			//�����ھӵ���ͬ,��ʼ�ж��ھӱߵ���Ϣ�Ƿ���ͬ
			int i_neibedgeindex = iteriadj->second;
			int v_neibedgeindex = itervadj->second;
			set<vector<int>> i_neib_einfo{ Pmg.edge_info[i_neibedgeindex].begin(),Pmg.edge_info[i_neibedgeindex].end() };
			set<vector<int>> v_neib_einfo{ Pmg.edge_info[v_neibedgeindex].begin(),Pmg.edge_info[v_neibedgeindex].end() };

			if (i_neib_einfo.size() != v_neib_einfo.size()) return 0;  //���·��ʵ���������ͬ������
			set<vector<int>> einfo_result;
			set_intersection(i_neib_einfo.begin(), i_neib_einfo.end(), v_neib_einfo.begin(), v_neib_einfo.end(), inserter(einfo_result, einfo_result.begin()));
			/*cout << "i��v��ͬ���ھ�·��ʵ���������" << einfo_result.size() << endl;*/
			if (einfo_result.size() != i_neib_einfo.size()) return 0;
		}
		else {
			break;
		}
	}
	return 1;
	//set<int> i_adj{ Pmg.adjlist[i].begin(),Pmg.adjlist[i].end() };
	//set<int> v_adj{ Pmg.adjlist[v].begin(),Pmg.adjlist[v].end() };
	//if (type == 2) {  //type2���͵ģ�Ҫ��i��v�ȴ��ھ���ɾ�����ж�
	//	i_adj.erase(v);
	//	v_adj.erase(i);
	//}
	//set<int> result;  //�ж��ھ��Ƿ���ͬ
	//set_intersection(i_adj.begin(), i_adj.end(), v_adj.begin(), v_adj.end(), inserter(result, result.begin()));
	//if (result.size() != i_adj.size()) return 0;

	//for (int neib : result) {  //�ж����ھ������ıߵ���Ϣ�Ƿ���ͬ
	//	set<vector<int>> i_neib_einfo{ Pmg.edge_info[Pmg.Mp_edge[i][neib]].begin(),Pmg.edge_info[Pmg.Mp_edge[i][neib]].end() };
	//	set<vector<int>> v_neib_einfo{ Pmg.edge_info[Pmg.Mp_edge[v][neib]].begin(),Pmg.edge_info[Pmg.Mp_edge[v][neib]].end() };
	//	
	//	if (i_neib_einfo.size() != v_neib_einfo.size()) return 0;  //���·��ʵ���������ͬ������
	//	set<vector<int>> einfo_result;
	//	set_intersection(i_neib_einfo.begin(), i_neib_einfo.end(), v_neib_einfo.begin(), v_neib_einfo.end(), inserter(einfo_result, einfo_result.begin()));
	//	/*cout << "i��v��ͬ���ھ�·��ʵ���������" << einfo_result.size() << endl;*/
	//	if (einfo_result.size() != i_neib_einfo.size()) return 0;

	//}
	//return 1;
} 



//void CBC_(int s) {  //pred��vector
//	//initialize variable
//	/*cout << "Դ��" << s << endl;*/
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
//		for (auto w : Pmg.adjlist_edge[v]) {  //������v�������ھ�w[�ڽӶ��㣬�ڽӱ�]
//			if (dist[w[0]] == DBL_MAX) {
//				dist[w[0]] = dist[v] + 1;
//				Q.push(w[0]);
//			}
//			if (dist[w[0]] == dist[v] + 1) {  //ͨ����if���������Եõ�v�������·����ھ�w
//				/*pathnum[w[0]] += pathnum[v]*Pmg.Mp[v][w[0]];*/
//				//ͳ��v�·���·��ʵ����Ϣ
//				vector<vector<int>> einfo = Pmg.edge_info[w[1]];
//				for (auto pinfo : einfo) {
//					int pclass = pinfo[0];
//					int pnum = pinfo[1];
//					pathinstance[pclass] += pnum;
//				}
//			}
//		}
//		//������v���ھӵ�󣬲���ͳ��������v��·��ʵ������Ϣ�����v����ÿ�����õ�w֮��ߵ�Ȩ��
//		for (auto w : Pmg.adjlist_edge[v]) {
//			double weight = 0;
//			if (dist[w[0]] == dist[v] + 1) {  //����v��w[0]֮��ı�Ȩ
//				vector<vector<int>> einfo = Pmg.edge_info[w[1]];  //w[1]��v��w[0]֮��ı�
//				for (auto pinfo : einfo) {  //����v��w[0]֮��ߵ�·��ʵ����Ϣ������pathnum[w[0]]
//					int pclass = pinfo[0];
//					int pnum = pinfo[1];
//					weight += ((double)pnum / pathinstance[pclass]);  //�ۼ�v-w[0]֮��ı�Ȩ
//				}
//				pathnum[w[0]] += pathnum[v] * weight;
//				pred[w[0]].push_back({ (double)v,weight });
//			}
//		}
//	}
//	/*for (int i = 0; i < pathnum.size(); i++) {
//		cout << "����" << i << "���·������" << pathnum[i] << endl;
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
//		cout << "�Ե�" << i << "Դ������" << deta[i] << endl;
//	}*/
//}
void CBC_ident0(int s) {  //pred��vector +ident��  s����proxy
	//initialize variable
	queue<int> Q;
	queue<int> temp_Q; //��¼���ʹ���level�㶥��

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
	map<int, int> pathinstance;  //w�·�����·��ʵ������ͳ��

	while (!Q.empty()) {
		int v = Q.front();
		if (dist[v] == level) { //��v��level��Ķ��㣬�����v�����ò�v���·���·��ʵ����Ϣ��ͳ�Ƶ�pathinstance����
			Q.pop();
			S.push(v);
			//ͳ��level��·��ʵ�����
			temp_Q.push(v);   //��¼��level����ʵĶ���
			for (auto w : Pmg.adjlist_edge[v]) {  //������v�������ھ�w[�ڽӶ��㣬�ڽӱ�]
				if (dist[w.first] == DBL_MAX) {
					dist[w.first] = dist[v] + 1;
					Q.push(w.first);
				}
				if (dist[w.first] == dist[v] + 1) {  //ͨ����if���������Եõ�v�������·����ھ�w
					//ͳ��v�·���·��ʵ����Ϣ
					vector<vector<int>> einfo = Pmg.edge_info[w.second];
					for (auto pinfo : einfo) {
						int pclass = pinfo[0];
						int pnum = pinfo[1];
						pathinstance[pclass] += pnum * (1 + ident[w.first]);  //v--w[0]��v--w[0]����Ķ����·��ʵ����
					}
				}
			}
		}
		else {  //��v��level+1��ģ�˵��level�㶥���Ѿ���������ˣ������ٷ���level��Ķ���rv�������Ȩ���������·����
			while (!temp_Q.empty()) {
				int rv = temp_Q.front();
				temp_Q.pop();
				//���rv����ÿ�����õ�w֮��ߵ�Ȩ�أ����������·����
				for (auto w : Pmg.adjlist_edge[rv]) {
					double weight;
					weight = 0;
					if (dist[w.first] == dist[rv] + 1) {  //����v��w[0]֮��ı�Ȩ
						vector<vector<int>> einfo = Pmg.edge_info[w.second];  //w[1]��v��w[0]֮��ı�
						for (auto pinfo : einfo) {  //����v��w[0]֮��ߵ�·��ʵ����Ϣ������pathnum[w[0]]
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
			level = level + 1;  //���������ʵ�level+1��
			pathinstance.clear(); //����ͳ�Ƶ�level+1���pathinstance��Ϣ
		}
	
	}
	pmid = clock();
	//back propagation  s����proxy����
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
	//cout << s<<"�Ե�16Դ������" << deta[16] << endl;
	/*deta0[s] = deta[0];*/
	/*for (int i = 0; i < deta.size(); i++) {
		cout << "�Ե�" << i << "Դ������" << deta[i] << endl;
	}*/
}
void CBC_ident1(int s) {  //pred��vector +ident��  s��proxy����  �����¶��壬��s�����ÿ����ΪԴ�����bfs�õ��Ľ������ͬ
	//initialize variable
	queue<int> Q;
	queue<int> temp_Q; //��¼���ʹ���level�㶥��

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
	map<int, int> pathinstance;  //w�·�����·��ʵ������ͳ��

	while (!Q.empty()) {
		int v = Q.front();
		if (dist[v] == level) { //��v��level��Ķ��㣬�����v�����ò�v���·���·��ʵ����Ϣ��ͳ�Ƶ�pathinstance����
			Q.pop();
			S.push(v);
			//ͳ��level��·��ʵ�����
			temp_Q.push(v);   //��¼��level����ʵĶ���
			for (auto w : Pmg.adjlist_edge[v]) {  //������v�������ھ�w[�ڽӶ��㣬�ڽӱ�]
				if (dist[w.first] == DBL_MAX) {
					dist[w.first] = dist[v] + 1;
					Q.push(w.first);
				}
				if (dist[w.first] == dist[v] + 1) {  //ͨ����if���������Եõ�v�������·����ھ�w
					//ͳ��v�·���·��ʵ����Ϣ
					vector<vector<int>> einfo = Pmg.edge_info[w.second];
					for (auto pinfo : einfo) {
						int pclass = pinfo[0];
						int pnum = pinfo[1];
						pathinstance[pclass] += pnum * (1 + ident[w.first]);  //v--w[0]��v--w[0]����Ķ����·��ʵ����
					}
				}
			}
			if (v == s) {  //��sΪident�㼯�е�proxy��ʱ����һ��bfs�����˷����ھӵ�w�⣬��ʵ������identical�㣨Ҳ��s���ھӵ㣩Ҳ�ڵڶ��㶥���У���ΪҪͳ��������s�·�·��ʵ������Ϣ������ҲҪ����s����sident��·��ʵ����Ϣ
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
		else {  //��v��level+1��ģ�˵��level�㶥���Ѿ���������ˣ������ٷ���level��Ķ���rv�������Ȩ���������·����
			while (!temp_Q.empty()) {
				int rv = temp_Q.front();
				temp_Q.pop();
				//���rv����ÿ�����õ�w֮��ߵ�Ȩ�أ����������·����
				for (auto w : Pmg.adjlist_edge[rv]) {
					double weight;
					weight = 0;
					if (dist[w.first] == dist[rv] + 1) {  //����v��w[0]֮��ı�Ȩ
						vector<vector<int>> einfo = Pmg.edge_info[w.second];  //w[1]��v��w[0]֮��ı�
						for (auto pinfo : einfo) {  //����v��w[0]֮��ߵ�·��ʵ����Ϣ������pathnum[w[0]]
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
			level = level + 1;  //���������ʵ�level+1��
			pathinstance.clear(); //����ͳ�Ƶ�level+1���pathinstance��Ϣ
		}
		//int v = Q.front();
		//Q.pop();
		//S.push(v);
		//map<int, int> pathinstance;
		//for (auto w : Pmg.adjlist_edge[v]) {  //������v�������ھ�w[�ڽӶ��㣬�ڽӱ�]
		//	if (dist[w.first] == DBL_MAX) {
		//		dist[w.first] = dist[v] + 1;
		//		Q.push(w.first);
		//	}
		//	if (dist[w.first] == dist[v] + 1) {  //ͨ����if���������Եõ�v�������·����ھ�w
		//		/*pathnum[w[0]] += pathnum[v]*Pmg.Mp[v][w[0]];*/
		//		//ͳ��v�·���·��ʵ����Ϣ
		//		vector<vector<int>> einfo = Pmg.edge_info[w.second];
		//		for (auto pinfo : einfo) {
		//			int pclass = pinfo[0];
		//			int pnum = pinfo[1];
		//			pathinstance[pclass] += pnum * (1 + ident[w.first]);  //v--w[0]��v--w[0]����Ķ����·��ʵ����
		//		}
		//	}
		//}
		//if (v == s) {  //��sΪident�㼯�е�proxy��ʱ����һ��bfs�����˷����ھӵ�w�⣬��ʵ������identical�㣨Ҳ��s���ھӵ㣩Ҳ�ڵڶ��㶥���У���ΪҪͳ��������s�·�·��ʵ������Ϣ������ҲҪ����s����sident��·��ʵ����Ϣ
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
		////������v���ھӵ�󣬲���ͳ��������v��·��ʵ������Ϣ�����v����ÿ�����õ�w֮��ߵ�Ȩ��
		//for (auto w : Pmg.adjlist_edge[v]) {
		//	double weight = 0;
		//	if (dist[w.first] == dist[v] + 1) {  //����v��w[0]֮��ı�Ȩ
		//		vector<vector<int>> einfo = Pmg.edge_info[w.second];  //w[1]��v��w[0]֮��ı�
		//		for (auto pinfo : einfo) {  //����v��w[0]֮��ߵ�·��ʵ����Ϣ������pathnum[w[0]]
		//			int pclass = pinfo[0];
		//			int pnum = pinfo[1];
		//			weight += ((double)pnum / pathinstance[pclass]);  //�ۼ�v-w[0]֮��ı�Ȩ
		//		}
		//		weight = weight * einfo.size();
		//		pathnum[w.first] += pathnum[v] * weight * (1 + (double)ident[v]);
		//		pred[w.first].push_back({ (double)v,weight });
		//	}
		//}
	}
	/*for (int i = 0; i < pathnum.size(); i++) {
		cout << "����" << i << "���·������" << pathnum[i] << endl;
	}*/
	pmid = clock();
	//back propagation  s��proxy����
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
	/*cout << "�Ե�0Դ������" << deta[0] << endl;*/

	/*for (int i = 0; i < deta.size(); i++) {
		cout << "�Ե�" << i << "Դ������" << deta[i] << endl;
	}*/
}
//void CBC_ident1(int s) {  //pred��vector +ident��  s��Ϊ���bfs����pathnum������identical���㼯�еĵ㣬֮��ıߵ���Ϣ����ͬ��
//	//initialize variable
//	/*cout << "Դ��" << s << endl;*/
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
//					pathinstance_s[i][pclass] += pnum * (1 + ident[w[0]]);  //v--w[0]��v--w[0]����Ķ����·��ʵ����
//				}
//			}
//		}
//	}
//	for (int si = 0; si < swithidentset.size(); si++) {  //����idents�в�ͬ�ĵ�ΪԴ�㣬�õ���temps�·���·��ʵ����������pathinstance_v
//		int temps = swithidentset[si];
//		for (int sident : swithidentset) {
//			if (sident != temps) {
//				int eindex = Pmg.Mp_edge[temps][sident];
//				vector<vector<int>> einfo = Pmg.edge_info[eindex];
//				for (auto pinfo : einfo) {
//					int pclass = pinfo[0];
//					int pnum = pinfo[1];
//					pathinstance_s[si][pclass] += pnum;  //��Ϊsidentһ��ֻ�������Լ������Բ���Ҫ*ident[sident](ident[w[0]])
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
//				pathnum[si][w[0]] += pathnum[si][temps] * weight;  //����Ҫ*ident[temps]��ident[v]������Ϊ����ͬ��idents�ֳɲ�ͬ��pathnum����
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
//		for (auto w : Pmg.adjlist_edge[v]) {  //������v�������ھ�w[�ڽӶ��㣬�ڽӱ�]
//			if (dist[w[0]] == DBL_MAX) {
//				dist[w[0]] = dist[v] + 1;
//				Q.push(w[0]);
//			}
//			if (dist[w[0]] == dist[v] + 1) {  //ͨ����if���������Եõ�v�������·����ھ�w
//				/*pathnum[w[0]] += pathnum[v]*Pmg.Mp[v][w[0]];*/
//				//ͳ��v�·���·��ʵ����Ϣ
//				vector<vector<int>> einfo = Pmg.edge_info[w[1]];
//				for (auto pinfo : einfo) {
//					int pclass = pinfo[0];
//					int pnum = pinfo[1];
//					pathinstance[pclass] += pnum * (1 + ident[w[0]]);  //v--w[0]��v--w[0]����Ķ����·��ʵ����
//				}
//			}
//		}
//		//������v���ھӵ�󣬲���ͳ��������v��·��ʵ������Ϣ�����v����ÿ�����õ�w֮��ߵ�Ȩ��
//		for (auto w : Pmg.adjlist_edge[v]) {
//			double weight = 0;
//			if (dist[w[0]] == dist[v] + 1) {  //����v��w[0]֮��ı�Ȩ
//				vector<vector<int>> einfo = Pmg.edge_info[w[1]];  //w[1]��v��w[0]֮��ı�
//				for (auto pinfo : einfo) {  //����v��w[0]֮��ߵ�·��ʵ����Ϣ������pathnum[w[0]]
//					int pclass = pinfo[0];
//					int pnum = pinfo[1];
//					weight += ((double)pnum / pathinstance[pclass]);  //�ۼ�v-w[0]֮��ı�Ȩ
//				}
//				for (int si = 0; si < swithidentset.size(); si++) {
//					pathnum[si][w[0]] += pathnum[si][v] * weight * (1 + (double)ident[v]);
//				}
//				pred[w[0]].push_back({ (double)v,weight });
//			}
//		}
//	}
//	/*for (int i = 0; i < pathnum.size(); i++) {
//		cout << "����" << i << "���·������" << pathnum[i] << endl;
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
//				cout<<"w="<<w<<"ʱ��"<< "deta[(int)"<<19<<"] += ("<<pathnum[(int)v[0]]<<" * "<<v[1]<<" / "<<pathnum[w]<<") * (1 + "<<deta[w]<<") * ((double)1 + "<<ident[w]<<")"<<endl;
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
//		cout << swithidentset[si]<<"�Ե�0Դ������" << deta[si][0] << endl;
//	}*/
//	for (int si = 0; si < swithidentset.size(); si++) {
//		deta0[si] = deta[si][0];
//	}
//
//	/*for (int i = 0; i < deta.size(); i++) {
//		cout << "�Ե�" << i << "Դ������" << deta[i] << endl;
//	}*/
//}
//void CBC_ri(int s) {  //pred��vector +reach+ident��
//	//initialize variable
//	/*cout << "Դ��" << s << endl;*/
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
//		for (auto w : Pmg.adjlist_edge[v]) {  //������v�������ھ�w[�ڽӶ��㣬�ڽӱ�]
//			if (dist[w[0]] == DBL_MAX) {
//				dist[w[0]] = dist[v] + 1;
//				Q.push(w[0]);
//			}
//			if (dist[w[0]] == dist[v] + 1) {  //ͨ����if���������Եõ�v�������·����ھ�w
//				/*pathnum[w[0]] += pathnum[v]*Pmg.Mp[v][w[0]];*/
//				//ͳ��v�·���·��ʵ����Ϣ
//				vector<vector<int>> einfo = Pmg.edge_info[w[1]];
//				for (auto pinfo : einfo) {
//					int pclass = pinfo[0];
//					int pnum = pinfo[1];
//					pathinstance[pclass] += pnum;
//				}
//			}
//		}
//		//������v���ھӵ�󣬲���ͳ��������v��·��ʵ������Ϣ�����v����ÿ�����õ�w֮��ߵ�Ȩ��
//		for (auto w : Pmg.adjlist_edge[v]) {
//			double weight = 0;
//			if (dist[w[0]] == dist[v] + 1) {  //����v��w[0]֮��ı�Ȩ
//				vector<vector<int>> einfo = Pmg.edge_info[w[1]];  //w[1]��v��w[0]֮��ı�
//				for (auto pinfo : einfo) {  //����v��w[0]֮��ߵ�·��ʵ����Ϣ������pathnum[w[0]]
//					int pclass = pinfo[0];
//					int pnum = pinfo[1];
//					weight += ((double)pnum / pathinstance[pclass]);  //�ۼ�v-w[0]֮��ı�Ȩ
//				}
//				pathnum[w[0]] += pathnum[v] * weight * (1 + (double)ident[v]);
//				pred[w[0]].push_back({ (double)v,weight });
//			}
//		}
//	}
//	/*for (int i = 0; i < pathnum.size(); i++) {
//		cout << "����" << i << "���·������" << pathnum[i] << endl;
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
//		cout << "�Ե�" << i << "Դ������" << deta[i] << endl;
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
		cout <<"����"<< bc_[i].first << "��bcֵ:" << bc_[i].second << endl;
	}
}

void getbcfromfile() {
	string filename = file + "/result.txt";  
	ifstream input(filename, ios::in);  //�����ļ�������input
	//�ж��ļ��Ƿ����
	if (!input) {
		cerr << "file error!" << endl;
		exit(1);
	}

	//�������ж����bcֵ//����ÿ�������bcֵ���ļ�
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
	cout << "��cbc�������bcֵ�ѱ�����" << file << "/cbc_result_amwma.txt";
}

/*���ļ�����ȡ�칹ͼ��Ϣ*/
void hetergraph::read_file() {  //δ��adjlist
	string filename = file + "/base.txt";  //base�д洢�ļ��Ļ�����Ϣ�����ݱ������ж�Ӧ�ıߵ��������ҵ��ߵ��ļ�"x.txt"
	ifstream input(filename, ios::in);  //�����ļ�������input
	//�ж��ļ��Ƿ����
	if (!input) {
		cerr << "file error!" << endl;
		exit(1);
	}

	//��һ�У����������
	string vinfo;
	getline(input, vinfo);
	istringstream svinfo(vinfo);
	svinfo >> n_vertexclass;

	//�������vclass-��������(i)  i�ඥ����vcnum
	for (int i = 0; i < n_vertexclass; i++) {
		string vcinfo;
		getline(input, vcinfo);
		istringstream svcinfo(vcinfo);
		char vclass;
		int vcnum;
		svcinfo >> vclass >> vcnum;
		vertex_class_index[vclass] = i;  //�������-������
		vertex_num.push_back(vcnum);    //ÿ�ඥ����
	}

	//�ߵ������
	string einfo;
	getline(input, einfo);
	istringstream seinfo(einfo);
	seinfo >> n_edgeclass;

	//�ߵ�������Ϣ
	for (int i = 0; i < n_edgeclass; i++) {  //��i���
		string eeinfo;
		getline(input, eeinfo);
		istringstream seeinfo(eeinfo);
		int vx, vy, edgenum;
		seeinfo >> vx >> vy >> edgenum;
		edge_class.push_back({ vx,vy });  //��i������ӵĶ����������
		edge_num.push_back(edgenum);

		//��ñ����i�����Ϣ���ļ����������ļ�
		stringstream si;
		si << i;
		string edge_filename = file + "/edge/" + si.str() + ".txt";
		ifstream edge_input(edge_filename, ios::in);  //�����ļ�������input
		//�ж��ļ��Ƿ����
		if (!edge_input) {
			cerr << "file error!" << endl;
			exit(1);
		}
		//�������i�����Ϣ���ļ������õ�����ߵ��ڽӾ���
		string edgeinformation;
		getline(edge_input, edgeinformation);/***edge�ļ���һ���Ƿ��б�ʾ�ߵ���Ϣ***/
		int nx = vertex_num[vx], ny = vertex_num[vy];
		map<int, map<int, int>> edge_matrix;
		string edgedata;
		while (getline(edge_input, edgedata)) {
			istringstream sedgedata(edgedata);
			int x, y;
			sedgedata >> x >> y;
			//��x�� ��y��Ϊ1
			auto iterx = edge_matrix.find(x);  //x�Ƿ��Ѿ�����
			if (iterx == edge_matrix.end()) {  //��δ����
				map<int, int> xadjlist;
				xadjlist.insert(map<int, int>::value_type(y, 1));  //��y����x���ڽ�list
				edge_matrix.insert(map<int, map<int, int>>::value_type(x, xadjlist)); //��x���ڽ�list����edge_matrix
				xadjlist.clear();
			}
			else {  //x�Ѿ�����
				iterx->second.insert(map<int, int>::value_type(y, 1));
			}
		}
		edge_info.push_back(edge_matrix);  //��i��ߵ��ڽӾ���
		edge_matrix.clear();
		edge_input.close();
	}
	input.close();
}

/*���ݱߵ������������������ţ����ұߵ����������*/
int hetergraph::getEdgeIndex(int x, int y) {
	int index = -1;
	for (auto it : edge_class) {
		index++;
		if (x == it[0] && y == it[1])
			break;
	}
	return index;
}

/*���칹ͼ�У�����P-multigraph*/
void multigraph::getPmg(vector<char>* P) {  //��adjlist
	vector<char>::iterator it = P->begin();  //it��������Ԫ·��P�еĶ������
	/*Pmg������*/
	n_vertex = hg.vertex_num[hg.vertex_class_index[*it]];
	n_vertex_org = n_vertex;
	adjlist_edge.resize(n_vertex);  //n��ai������ڽӶ���+�ߵı��
	//adjlist.resize(n_vertex);

	/*����P����P�бߵ�������������eindex��*/
	vector<int> eindex;
	for (int i = 0; i < P->size() / 2 ; i++) {  //һ������P->size()/2����
		int x = hg.vertex_class_index[*it++];
		int y = hg.vertex_class_index[*it];
		int index = hg.getEdgeIndex(x, y);  //���������յ�����ҵ�����߶�Ӧ���������Ӷ��õ�����ߵ��ڽӾ���
		if (index == -1) {
			cerr << "edge error!" << endl;
			exit(1);
		}
		eindex.push_back(index);
	}
	cout << "�ҵ���������";
	for (int i : eindex) cout << i << "-";
	cout << endl;
	/*����Mpl*/
	int p = 0;
	map<int, map<int, int>> matrixpl = hg.edge_info[eindex[p]];
	hg.edge_info[eindex[p]].clear();  //��ʱ����һ��eindex���ڽӾ����Ѿ�������matrixpl�У���Ϣ������Ҫ
	while (++p < eindex.size()) {  //��matrix����һ���ڽӾ���edge_info[eindex[p]]���
		int curr_eindex = eindex[p];  //��ǰ����ıߵ�����,����ǰ���ĸ��ߵ��ڽӾ�����MPL���
		cout << "����" << eindex[p - 1] << '*' << curr_eindex << "���ڽӾ���" << endl;
		double ps11, pt11, time11;
		ps11 = clock();
		map<int, map<int, int>> tmatrix;  //��ʱ��¼������
		auto iterik = matrixpl.begin();
		for (iterik; iterik != matrixpl.end(); iterik++) {  //��ÿ�����ڽ����е�i
			int i = iterik->first;
			auto iterk = iterik->second.begin();  //iterk:��i���ڽӵ�k����
			for (iterk; iterk != iterik->second.end(); iterk++) {
				int k = iterk->first;
				int currik = iterk->second;  //(y,1) ���ڼ�>0 ��ʾi��k�ڽӣ�֮������k�ڽӵ�j
				auto iterkj = hg.edge_info[curr_eindex].find(k);
				if (iterkj != hg.edge_info[curr_eindex].end()) {  //k*j���ڽӱ��д��ڵ�k���ڽ�����
					auto iterj = iterkj->second.begin();  //iterj����k���ڽӵ�j����
					for (iterj; iterj != iterkj->second.end(); iterj++) {
						int j = iterj->first;
						int currkj = iterj->second;
						//�����У�����i���ڽӵ�j���ڽӶ��ر���Ϊcurrik*currkj
						auto iterri = tmatrix.find(i);
						if (iterri == tmatrix.end()) {  //����δ����i����������һ���ڽӵ�k
							map<int, int> iadjlistj;
							iadjlistj.insert(map<int, int>::value_type(j, currik * currkj));
							tmatrix.insert(map<int, map<int, int>>::value_type(i, iadjlistj));
							iadjlistj.clear();
						}
						else {  //��i�Ѹ��£��ж�j�Ƿ��Ѳ��뵽i���ڽӱ���
							auto iterrj = iterri->second.find(j);
							if (iterrj == iterri->second.end()) {  //��jδ��Ϊi���ھӱ����뵽i���ڽӱ���
								iterri->second.insert(map<int, int>::value_type(j, currik * currkj));
							}
							else {
								int oldij = iterrj->second;
								int nowij = oldij + currik * currkj;
								iterrj->second = nowij;  //����ijλ�õ�ֵ
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
		cout << " ��ʱ��" << time11 << endl;
		out_result << "����" << eindex[p - 1] << "*" << curr_eindex << "��ʱ:" << time11 << endl;
	}
	hg.edge_info.clear();
	eindex.clear();

	//save_Mpl(&matrixpl);  //����Mpl
    //�õ�Mpl�ڽӱ��ת��
	double ps12, pt12, time12;
	ps12 = clock();
	map<int, map<int, int>> matrixplt;
	for (auto it : matrixpl) {
		int x = it.first;
		auto ity = it.second.begin();
		for (ity; ity != it.second.end(); ity++) {
			int y = ity->first;
			int valuexy = ity->second;
			//y��x�ڽӣ�ֵΪvaluexy
			auto itfindy = matrixplt.find(y);
			if (itfindy == matrixplt.end()) {  //yδ����
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
	cout << "ת�ü������" << endl;
	double pt14 = clock();
	double time14 = (double)(pt14 - ps12) / CLOCKS_PER_SEC;
	out_result << "����mpl��ת����ʱ��" << time14 << endl;


	/*�����Mpl�󣬼���Mp����,�õ�edge_info��adjlist_edge��Mp_edge??*/
	cout << "����Mp�ڽӱ�" << endl;
	int edge_index = 0;
	for (int i = 0; i < n_vertex; i++) {  //�����i������ڽ�����
	/*cout << "���㵽��" << i;*/
		auto iterik = matrixpl.find(i);  //�ҵ�i�����k�ڽ�����
		if (iterik != matrixpl.end()) {  //�ҵ���
			map<int, vector<vector<int>>> ij_einfo;  //��¼i�������ھ�j����֮��ߵ���Ϣ
			auto iterk = iterik->second.begin();
			for (iterk; iterk != iterik->second.end(); iterk++) {  //����i��ÿ���ڽӵ�k
				int k = iterk->first;
				int currik = iterk->second;
				auto iterkj = matrixplt.find(k);//��k���ڽӵ�j����
				if (iterkj != matrixplt.end()) {
					auto iterj = iterkj->second.begin();
					for (iterj; iterj != iterkj->second.end(); iterj++) {
						int j = iterj->first;
						int currkj = iterj->second;
						if (j > i) {
							int ij_kclasspathnum = currik * currkj;  //��ij��֮�䣬k��·��ʵ����
							auto iterjeinfo = ij_einfo.find(j);  //�Ƿ��Ѿ����¹�i���ھ�j�ıߵ���Ϣ
							if (iterjeinfo == ij_einfo.end()) {  //û���¹�i���ھ�j�ıߵ���Ϣ
								vector<vector<int>> temp_ij_einfo;
								temp_ij_einfo.push_back({ k,ij_kclasspathnum });
								ij_einfo.insert(map<int, vector<vector<int>>>::value_type(j, temp_ij_einfo));
								temp_ij_einfo.clear();
							}
							else {  //�Ѿ����¹�i���ھ�j�ıߵ���Ϣ
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
	cout << "��mpl�õ�mp��ʱ��" << time12 << endl;
	out_result << "��mpl�õ�mp��ʱ��" << time12 << endl;

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

	//��a��ÿ���ڽӵ�
	for (auto an : adjlist_edge[a]) {
		//�޸��ڽӾ���
		/*Mp[a][an] = 0; Mp[an][a] = 0;*/
		//��an���ڽӱ���ɾ��a
		map<int, int>::iterator it = adjlist_edge[an.first].begin();
		for (it; it != adjlist_edge[an.first].end(); it++) {
			if (it->first == a) {
				adjlist_edge[an.first].erase(it);
				break;
			}
		}
	}
	//���a���ڽӱ���ʱaΪP-multigraph��һ�������ĵ㣬��ͨ��index������Ȼ�ܱ�����a��
	map<int, int>().swap(adjlist_edge[a]);
}

/*��P-multigraph��ɾ����*/
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