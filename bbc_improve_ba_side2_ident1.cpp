//�ű�+�ؽڶ���+side��+identical��(�ű߷ָ��õ���j_index--�ķ�ʽ���ؽڶ���ʶ���õ�����ǰ�����ھӶ���ķ�ʽ)
//side���ʶ��λ�ڼ���Mp�Ĺ����У�����Mpl����Ѱ��side�㣻identical���ʶ��λ���űߡ��ؽڷָside���Ƴ���
//�Ķ���1��void todoType1Ident() set<int> i2adjlistͨ��set�õ����ظ��Ķ����ھ����У�2��mergeIdent1()�У���Щv������identical type2�Ķ��㣬�ںϲ�ʱ��Ҳ���뵽identSet[i]��
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
	int getEdgeIndex(int x, int y);
};


class multigraph {
public:
	int n_vertex;  //���¶�����������copy�Ĺؽڶ��㣩
	int n_vertex_org;  //ԭ������
	int m_edge;  //������
	//vector<vector<int> > Mp;
	vector<map<int, int> > adjlist;
	void show_mg();
	void getPmg(vector<char>* P);  //���칹ͼ�л��Pmg
	void getPmg_read_file();      //�ӱ�����ļ��ж�ȡMpl��Pmg
	void deleteEdge(int a, int b);  //����ab���ڽӱ���ɾ��
	void deleteVertex(int a);
	void save_Mpl(vector<vector<int>>* Mpl);  //����Mpl
	void save_Pmg();  //��������Pmg
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
map<int, int> org;  //map[100]=1  ��ʾ����100�ڹؽڶ��㻮��ǰ��ԭ������1
map<int, vector<int>> SideSet0;  //�ڼ���Mpʱ�õ��Ķ���ڵ�side���㣬��δ��P-multigraph���зָ�ǰ���õ���side��
							  //side0[a]={Aset},aΪP���м�����Ķ��㣬AsetΪA1��Ķ��㼯����ʾ��ЩA1�ඥ��ֻ�붥��a����
vector<int> flag_side;  //��־�����Ƿ���Ϊside���㱻ɾ��
map<int, vector<int>> IdentSet;  //�����������identical������б�
vector<int> ident;     //�����������identical����ĸ������������Լ���
vector<int> flag_addtoident;  //�Ѿ�����identical�㼯�Ķ���
ofstream out_result(file + "/al7_time_result.txt");
double timebfs = 0, timeback = 0;
double pmid;
int countb = 0, counta = 0;
int count_side_edge = 0, count_ident_edge = 0;
int count_ident2 = 0, count_ident1 = 0;

int main() {
	//һ�������칹ͼ����Pmg
	double ps00, pt00, time00;
	ps00 = clock();
	hg.read_file();
	pt00 = clock();
	time00 = (double)(pt00 - ps00) / CLOCKS_PER_SEC;
	out_result << "�����칹ͼ��ʱ��" << time00 << endl;
	cout<< "�����칹ͼ��ʱ��" << time00 << endl;

	/*hg.show_hg();*/
	double ps0, pt0, time0;
	ps0 = clock();
	vector<char> P = { 'A','P','V','P','A' };
	Pmg.getPmg(&P);  //�õ�P-multigraph
	pt0 = clock();
	time0 = (double)(pt0 - ps0) / CLOCKS_PER_SEC;
	cout << "�õ�Pmg+side�㼯��ʱ��" << time0 << endl;
	out_result << "�õ�Pmg+side�㼯��ʱ:" << time0 << endl;
	cout << "������" << Pmg.n_vertex << endl;
	out_result << "Pmg������" << Pmg.n_vertex << endl;
	cout << "������" << Pmg.m_edge << endl;
	out_result << "Pmg������" << Pmg.m_edge << endl;
	//�����Ѿ�������Pmg�����ļ��ж�ȡ
	/*double ps0, pt0, time0;
	ps0 = clock();

	Pmg.getPmg_read_file();
	pt0 = clock();
	time0 = (double)(pt0 - ps0) / CLOCKS_PER_SEC;
	cout << "����Pmg��ʱ��" << time0 << endl;*/

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
	reach.resize(Pmg.n_vertex, 1);

	double ps1, pt1, time1;
	ps1 = clock();
	BridgeEdgeDivide();  //�ű߷ָ�
	pt1 = clock();
	time1 = (double)(pt1 - ps1) / CLOCKS_PER_SEC;
	cout << "�ű߷ָ����,������"<<countb<<"���ű�,����ʱ��" << time1 << endl;
	out_result << "�ű߷ָ����,������" << countb << "���ű�,����ʱ��" << time1 << endl;

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
	cout << "�ؽڶ���ָ���ϣ�������"<<counta<<"���ؽڶ���,����ʱ��" << time2 << endl;
	out_result << "�ؽڶ���ָ���ϣ�������" << counta << "���ؽڶ���,����ʱ��" << time2 << endl;

	double pcmps1, pcmpt1, timecmp1;
	pcmps1 = clock();
	int maxComponent1 = getMaxComponent();
	cout << "�űߡ��ؽڶ���ָ�������ͨ������" << maxComponent1 << endl;
	out_result << "�űߡ��ؽڶ���ָ�������ͨ������" << maxComponent1 << endl;
	pcmpt1 = clock();
	timecmp1 = (double)(pcmpt1 - pcmps1) / CLOCKS_PER_SEC;
	out_result << "���������ͨ������С��ʱ��" << timecmp1 << endl;

	SideAndIdent();

	double pcmps2, pcmpt2, timecmp2;
	pcmps2 = clock();
	int maxComponent2 = getMaxComponent();
	cout << "�Ƴ�side���ϲ�identical�������ͨ������" << maxComponent2 << endl;
	out_result << "�Ƴ�side���ϲ�identical�������ͨ������" << maxComponent2 << endl;
	pcmpt2 = clock();
	timecmp2 = (double)(pcmpt2 - pcmps2) / CLOCKS_PER_SEC;
	out_result << "���������ͨ������С��ʱ��" << timecmp2 << endl;


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
		//	/*cout << "��" << i << "ΪԴ�����" << endl;*/
		//	BBC1(i);
		//	pt4 = clock();
		//	timebfs += (double)(pmid - ps4) / CLOCKS_PER_SEC;
		//	timeback += (double)(pt4 - pmid) / CLOCKS_PER_SEC;
		//	time4 += (double)(pt4 - ps4) / CLOCKS_PER_SEC;
		//}

		//ps4 = clock();
		//BBC1(i);  //�ؽڶ���+�ű�
		//pt4 = clock();
		//timebfs += (double)(pmid - ps4) / CLOCKS_PER_SEC;
		//timeback += (double)(pt4 - pmid) / CLOCKS_PER_SEC;
		//time4 += (double)(pt4 - ps4) / CLOCKS_PER_SEC;

		if (count % 10 == 0) {
			cout << count << "����ΪԴ����ʱ��" << time4 << endl;
		}
	}

	//���ݹؽڶ����ӳ���ϵ������ؽڶ����bcֵ
	for (auto it : org) {
		bc[it.second] += bc[it.first];
	}
	double ptt = clock();
	cout << "ǰ��bfs��ʱ�䣺" << timebfs << endl;
	cout << "�����ۻ���ʱ�䣺" << timeback << endl;
	cout << "����bcֵ��ʱ�䣺" << time4 << endl;
	out_result << "ǰ��bfs��ʱ�䣺" << timebfs << endl;
	out_result << "�����ۻ���ʱ�䣺" << timeback << endl;
	out_result << "����bcֵ��ʱ�䣺" << time4 << endl;

	double time_all = (double)(ptt - ps00) / CLOCKS_PER_SEC;
	cout << "����ʱ��" << time_all << endl;
	out_result << "����ʱ��" << time_all << endl;
	cout << "����ʱ���Ѵ���" << file + "/improve_time_result.txt" << endl;

	/*for (int i = 0; i < Pmg.n_vertex_org; i++) {
		if (bc[i] > 0) {
			cout << "����" << i << "��bcֵ��" << bc[i] << endl;
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
	cout << "��bbc�������bcֵ�ѱ�����" << file << "/result.txt" << endl;
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
		cout << "����" << bc_[i].first << "��bcֵ:" << bc_[i].second << endl;
	}
	vector<pair<int, double>>().swap(bc_);
}

void getbcfromfile() {
	string filename = file + "/result.txt";  //base�д洢�ļ��Ļ�����Ϣ�����ݱ������ж�Ӧ�ıߵ��������ҵ��ߵ��ļ�"x.txt"
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

/*�����ű�*/
void BridgeEdgeDivide() {
	dfn.resize(Pmg.n_vertex);
	low.resize(Pmg.n_vertex);
	for (int i = 0; i < Pmg.n_vertex; i++) {
		stamp = 0;
		int totalreach;
		if (!dfn[i]) {
			totalreach = countTotalReach(i);  //����i���ڵ���ͨ�����У����������
			/*cout << "totalreach:" << totalreach << endl;*/
			b_tarjan(i, -1, totalreach);  //��i���ڵ���ͨ����ִ��b_tarjan�������ű�
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
				//��index-j���ű�
				cout << "�ҵ��űߣ�" << index << "-" << j << endl;
				countb++;
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

				/*cout << "rep" << j << ":" << repj << "  rep" << index << ":" << repindex << endl;
				cout<<"reach"<<index<<":"<<reach[index]<< "  reach" << j << ":" << reach[j] << endl;*/


				//�ۼ�ȱ�ٵ�pair-dependency
				bc[index] += ((double)repindex - 1) * repj;
				bc[j] += ((double)repj - 1) * repindex;
				//�޸��ڽӱ��ڽӾ���
				Pmg.deleteEdge(index, j);
			}
		}
		else if (dfn[j] < dfn[index]) {
			low[index] = min(low[index], dfn[j]);
		}
	}
	
	map<int, int> ().swap(indexadj);
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
	map<int,int> indexadj = Pmg.adjlist[index];  //��Ϊ�·���artvertexcopy�У����ܻ��޸�Pmg.adjlist[index]�������Ƚ��䱣�����������޸�ǰ��index���ھ��ϱ���
	/*for (int indexj = 0; indexj < indexadj.size(); indexj++)*/
	for (auto iteradjj = indexadj.begin(); iteradjj != indexadj.end();iteradjj++) {
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
				for (nowj; nowj!=indexadj.end(); nowj++) {  //�ж�index���ھ��Ƿ񶼷��ʵ��ˣ����ж�index�Ƿ��������ھӽڵ�
					int nb = nowj->first;
					if (!dfn[nb]) {  //��δ���ʵ����ھӣ�˵��index�Ǹ�㣬��ǰջ�е����ж��㹹��һ����ͨ����
						cout << "�ҵ��ؽڶ���(����Դ��)��" << index << endl;
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
	//���ڵ�Ϊ���ʱ
	//if (fa < 0 && child>1) {  //�Ϸ������˵�index��Ϊ���ڵ�ʱ��indexΪ������ͨ����������indexΪ���ڵ㣬����ͨ��low[j]>=dfn[index]���ж�index�Ƿ�Ϊ���
	//	
	//}
	map<int, int>().swap(indexadj);
}

/*tarjan�㷨�ҵ��ؽڶ���󣬸��ƹؽڶ��㲢�����ڽӾ����ڽӱ�reachֵ*/
void artvertexcopy(int index, int j, int totalreach) {
	//����index�Ǹ��
	int copyindex = Pmg.n_vertex++;  //����ͨ�����У�copyһ��index��
	cout << "���ƹؽڵ�" << copyindex << endl;
	org[copyindex] = index;//����copyindex��ԭ����index

	Pmg.adjlist.resize(Pmg.n_vertex);  //��չadjlist
	int repcopyindex = 0, repindex = 0;
	while (stack_bi.top() != j) {  //��index-j�����ͨ�����еĶ����ջ
		int temp = stack_bi.top();
		stack_bi.pop();
		//����j����Ķ�����
		repcopyindex += reach[temp];
		/*�����ڽӾ�����ڽӱ�����ͨ�����еĶ�����index�ڽӣ��޸�Ϊ��copyindex�ڽ�*/
		auto itertemp = Pmg.adjlist[index].find(temp);
		if (itertemp != Pmg.adjlist[index].end()) {  //���temp������Pmultigraph����index�������������Ϊ��copyindex����
			int weightindex_temp = Pmg.adjlist[index][temp];
			Pmg.deleteEdge(index, temp);
			Pmg.adjlist[copyindex].insert(map<int, int>::value_type(temp, weightindex_temp));
			Pmg.adjlist[temp].insert(map<int, int>::value_type(copyindex, weightindex_temp));
		}
		//if (Pmg.Mp[temp][index] != 0) {  //���temp������Pmultigraph����index�������������Ϊ��copyindex����
		//	//����index-temp���ڽӾ���
		//	Pmg.Mp[temp][copyindex] = Pmg.Mp[temp][index];
		//	Pmg.Mp[copyindex][temp] = Pmg.Mp[temp][index];

		//	//��Pmg��ɾ����index-temp
		//	Pmg.deleteEdge(index, temp);
		//	Pmg.adjlist[copyindex].push_back(temp);
		//	Pmg.adjlist[temp].push_back(copyindex);
		//}
	}
	//���j��ջ������j�Ĵ���Ķ��������ڽӾ�����ڽӱ�
	stack_bi.pop();  repcopyindex += reach[j];
	int weightindex_j = Pmg.adjlist[index][j];
	//Pmg.Mp[copyindex][j] = Pmg.Mp[index][j]; Pmg.Mp[j][copyindex] = Pmg.Mp[j][index];
	Pmg.deleteEdge(index, j);
	Pmg.adjlist[copyindex].insert(map<int, int>::value_type(j, weightindex_j));
	Pmg.adjlist[j].insert(map<int, int>::value_type(copyindex, weightindex_j));

	repindex = totalreach - repcopyindex;  //repcopyindex��j������ͨ�����ϣ������reachֵ֮�ͣ�������index�����reach����repindex�Ƿָ���ͨ������index���ڷ�����reachֵ֮�ͣ�����index��reach��
	reach[index] += repcopyindex;  //���¶���index��reachֵ���¼���j���ڵķ����϶����reachֵ֮�ͣ���ͬ��index��reach����Ϊ����������
	reach.push_back(repindex);  //copyindex���¼���Ķ��㣬����Ҫ��reach��push_back��copyindex��reachֵʱindex������ͨ�����ж����reach֮��(����index�����reach)
}


/*side�����identical����*/
void SideAndIdent() {
	/*���side���㼯*/
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
	cout << "side������" << count_sideset << ",side������" << count_sidenum << endl;
	out_result<< "side������" << count_sideset << ",side������" << count_sidenum << endl;
	/*�Ƴ�side*/
	flag_side.resize(Pmg.n_vertex);//��Ϊ�����Ƿ�Ϊside��ı�־��������ÿ����ִ��bbc����������side��ΪԴ�����
	double psides, psidet, timeside;
	psides = clock();
	deleteSide0();  //delete
	SideSet0.clear();
	psidet = clock();
	timeside = (double)(psidet - psides) / CLOCKS_PER_SEC;
	cout << "�Ƴ�side���㣬��ʱ" << timeside << endl;
	out_result << "�Ƴ�side���㣬��ʱ" << timeside << endl;
	cout << "�Ƴ�side������ɾ���ı�����" << count_side_edge << endl;
	out_result << "�Ƴ�side������ɾ���ı�����" << count_side_edge << endl;

	double pidents, pidentt, timeident;
	double pidentmid;
	pidents = clock();
	ident.resize(Pmg.n_vertex);
	flag_addtoident.resize(Pmg.n_vertex);//��Ϊ�����Ƿ��ѱ�proxy�������ı�־��������ÿ����ִ��bbc���������Ըõ�ΪԴ�����
	cout << "��ʼʶ��identical2" << endl;
	todoType2Ident();
	pidentmid = clock();
	double timemid = (double)(pidentmid - pidents) / CLOCKS_PER_SEC;
	cout << "���У��ϲ�ʶ��identical2��ʱ��" << timemid << endl;
	out_result << "���У��ϲ�ʶ��identical2��ʱ��" << timemid << endl;
	cout << "��ʼʶ��identical1" << endl;
	todoType1Ident();
	pidentt = clock();
	timeident = (double)(pidentt - pidents) / CLOCKS_PER_SEC;
	cout << "ʶ�𡢺ϲ�identical���㣬��ʱ" << timeident << endl;
	out_result<< "ʶ�𡢺ϲ�identical���㣬��ʱ" << timeident << endl;
	cout << "�ϲ�identical������ɾ���ı�����" << count_ident_edge << endl;
	out_result << "�ϲ�identical������ɾ���ı�����" << count_ident_edge << endl;

	//���identical���㼯
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
	cout << "ident������" << count_identset << ",�ϲ��˵�ident��������������proxy�㣩" << count_identnum << endl;
	out_result << "ident������" << count_identset << ",�ϲ��˵�ident��������������proxy�㣩" << count_identnum << endl;
	cout << "ident1:" << count_ident1 << "  ident2:" << count_ident2 << endl;
	out_result << "ident1:" << count_ident1 << "  ident2:" << count_ident2 << endl;

}

//ɾ�����ʶ���side����
void deleteSide0() {
	double pdss, pdst, timeds=0;
	int sidecount = 0;
	for (auto it : SideSet0) {  //itΪmap��ÿ��same_side�㼯
		//��same_side�㼯�е�ÿ����Եĺϲ�����bcֵ
		pdss = clock();
		sidecount++;
		for (int i = 0, len = it.second.size(); i < len; i++) {
			if (len == 2 && Pmg.adjlist[it.second[0]].find(it.second[1])==Pmg.adjlist[it.second[0]].end()) break;  //�����ʶ���side�㼯�У���һ������ָ��������Ա�Ĺ����ţ���ʱ������Ա��side�㣬
																			//�����ű�ɾ��ʱ���Ѿ�����Ա֮��ı�ɾ������ʱ����Ա������ͬһ���㼯�е�side�㣬��������������side�㣬�ϲ�������pair-dependency��ʧ
			if (len == 1 && Pmg.adjlist[it.second[0]].size() == 0) break;  //�Թ����ĵ㣨һ�������͹�����������Ϊ�ű߷ָ�degree-1�㵼�¹�����
			int u = it.second[i];
			for (int j = i + 1; j < len; j++) {
				int v = it.second[j];
				bc[v] += ((double)reach[v] - 1) * reach[u];
				bc[u] += ((double)reach[u] - 1) * reach[v];
			}
		}
		//��same_side�㼯ΪԴ����п��ѣ��ۼӵ㼯�еĵ���������Դ����
		BCforSide0(it.second);
		//ɾ��it�����same_side�㼯
		for (int s : it.second) {
			count_side_edge += Pmg.adjlist[s].size();
			Pmg.deleteVertex(s);
		}
		pdst = clock();
		timeds += (double)(pdst - pdss) / CLOCKS_PER_SEC;
		cout << "�Ƴ�1��side������ʱ��" << timeds << endl;
		
	}
}

void BCforSide0(vector<int> samesideset) {
	int sidesetsize = samesideset.size();
	/*for (int s : samesideset) flag_side[s] = 1;*/
	queue<int> Q;
	stack<int> S;
	vector<vector<pair<int,int>>> pred(Pmg.n_vertex);  //samesideset��ÿ������ΪԴ����ѵõ���pred����ͬ����Ϊ��������ͬ
	vector<double> dist(Pmg.n_vertex, DBL_MAX);  //samesideset��ÿ������ΪԴ����ѵõ���distҲ��ͬ����Ϊ��������ͬ
	vector<vector<double>> pathnum(sidesetsize, vector<double>(Pmg.n_vertex));  //samesideset��ÿ������ΪԴ����ѵõ��ĵ�������������·������ͬ����Ϊ�ߵ�Ȩֵ��ͬ
	for (int s : samesideset) dist[s] = 0;
	for (int i = 0; i < sidesetsize; i++) {
		for (int s : samesideset)
			pathnum[i][s] = 1;
	}

	//bfs0
	//��samesideset�еĶ���Ϊ������ʱ��ѡ��samesideset[0]Ϊ����㣬����Ѱ�����ھӶ��㹹�������ʱ��������samesideset�еĶ��㣬��Ϊ�Ѿ��ϲ�
	int tempv = samesideset[0];
	S.push(tempv);
	for (auto w : Pmg.adjlist[samesideset[0]]) {
		int flag_wnotin = 1;
		for (int i = 0; i < sidesetsize; i++) {
			if (samesideset[i] == w.first) {
				flag_wnotin = 0; break;
			}
		}
		//��w����samesideset����ʱ
		if (flag_wnotin) {
			if (dist[w.first] == DBL_MAX) {
				dist[w.first] = dist[tempv] + 1;
				Q.push(w.first);
			}
			if (dist[w.first] == dist[tempv] + 1) {
				for (int i = 0; i < sidesetsize; i++) {  //��sideset��ÿ������ΪԴ����п��ѵõ������·����
					pathnum[i][w.first] += pathnum[i][samesideset[i]] * Pmg.adjlist[samesideset[i]][w.first];  //�����Ȩ����samesideset[i]��w.first֮���Ȩ��
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
				for (int i = 0; i < sidesetsize; i++) {  //��sideset��ÿ������ΪԴ����п��ѵõ������·����
					pathnum[i][w.first] += pathnum[i][v] * w.second;
				}
				pred[w.first].push_back({ v,w.second });
			}
		}
	}
	/*for (int i = 0; i < pathnum.size(); i++) {
		cout << "����" << i << "���·������" << pathnum[i] << endl;
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
		for (auto v : pred[w]) {  //��v��samesideset[0]ʱ����deta[i][v]�ĸ�����ʵ��ͬ����Ϊv������samesideset�����ж��㣬�������Mp[v][w]��ͬ�������õ������ۣ���Ϊ��w��samesideset[0]ʱ�����ۼ�bcֵ
			if (v.first == samesideset[0]) {
				//��i���pathnum��deta����������samesideset[i]��detaֵ������(else��)�Ƕ�ÿ��i��������deta[i][v]������v��samesidesetʱ��Ҫ����deta[i][samesetside[i]]
				for (int i = 0; i < sidesetsize; i++) {  //����sameside[i]��һ��pathnum��deta��v���и���
					deta[i][samesideset[i]] += (pathnum[i][samesideset[i]] * Pmg.adjlist[samesideset[i]][w] / pathnum[i][w]) * (1 + deta[i][w]);
				}
			}
			else {
				for (int i = 0; i < sidesetsize; i++) {  //����sameside[i]��һ��pathnum��deta��v���и���
					deta[i][v.first] += (pathnum[i][v.first] * v.second / pathnum[i][w]) * (1 + deta[i][w]);
				}
			}

			/*deta[v] += (pathnum[v] * Pmg.Mp[v][w] / pathnum[w]) * (1 + deta[w]);*/
		}
		if (w != samesideset[0]) {
			for (int i = 0; i < sidesetsize; i++) {
				bc[w] += deta[i][w] * reach[samesideset[i]] + reach[samesideset[i]] * (deta[i][w] - ((double)reach[w] - 1));
				/*if (w == 17) {
					cout << "Դ��" << samesideset[i] << "��17��Դ������" << deta[i][w] * reach[samesideset[i]] + reach[samesideset[i]] * (deta[i][w] - ((double)reach[w] - 1)) << endl;
				}*/
			}

			/*bc[w] += deta[w] * reach[s];*/
		}
	}
	for (int i = 0; i < sidesetsize; i++) {
		int ss = samesideset[i];
		/*if (ss == 17) {
			cout << "����17ɾ���󣬶�17��bcֵ��Ӱ�죺" << ((double)reach[ss] - 1) * (deta[i][ss] - ((double)reach[ss] - 1)) <<endl;
			cout << bc[ss] << endl;
		}*/
		bc[ss] += ((double)reach[ss] - 1) * (deta[i][ss] - ((double)reach[ss] - 1));  //��ssΪԴ�㣬���������reachֵ֮�͵���deta-reach[ss]+1

	}
	stack<int>().swap(S);
	vector<vector<pair<int, int>>>().swap(pred);
	vector<vector<double>>().swap(pathnum);
	vector<vector<double>>().swap(deta);
	/*for (int i = 0; i < sidesetsize; i++) {
		cout << "Դ��" << samesideset[i] << endl;
		for (int j = 0; j < Pmg.n_vertex; j++) {
			if (j == 17) {
				cout << "�Ե�" << j << "Դ����" << deta[i][j] << endl;
			}
		}
	}*/
}

void todoType2Ident() {
	for (int i = 0; i < Pmg.n_vertex; i++) {  //��index0����n���������Զ�index=k�ĵ㣬������ȥ�ж��Ƿ���0--k-1�ĵ��Ƿ�identical����Ϊ��֮ǰ�Ѿ��жϹ��ˣ����ǣ���flag[k]Ӧ����1��
		if (flag_side[i] == 0 && flag_addtoident[i] == 0) {  //side���㲻ȥѰ��identical���Ѿ��ҵ�proxy�Ķ��㲻����Ѱ��identical
			//��Ϊ��ʶ��һ��i��ident�󣬾���Pmg.adjlist[i]��ɾ������ʱx������һЩ��
			map<int,int> iadjlist = Pmg.adjlist[i];
			for (auto v : iadjlist) {
				if (v.first > i&&flag_addtoident[v.first]==0) {                //�ж�>i�Ķ���v�Ƿ���i��identical�ģ�<i��v�Ѿ��ڶ�v���ھӵ����ʱ�жϹ��ˣ���i������ͬ��
					int flag_isident = judgeIdent(i, v.first, 2);      //�ж�v�Ƿ���i��type2 ident��
					if (flag_isident) {
						//add dependency and modify Mp��adjlist and update ident set��ident value of i
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
				for (auto v : Pmg.adjlist[neibor.first]) {  //i�Ķ����ھ�
					if (v.first > i && Pmg.adjlist[i].find(v.first)==Pmg.adjlist[i].end()&&flag_addtoident[v.first]==0) {  //���ھӣ���δ���뵽�κζ����identset�С�����Щ����v���ܻ����һЩ����identical type2�Ķ��㣩
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
	//if (type == 2) {  //type2����
	//	for (auto ii : Pmg.adjlist[i]) {
	//		if (ii.first == v) continue;
	//		i_adj.insert({ ii.first,ii.second });
	//	}
	//	for (auto vv : Pmg.adjlist[v]) {
	//		if (vv.first == i) continue;
	//		v_adj.insert({ vv.first,vv.second });
	//	}
	//}
	//set<pair<int,int>> result;  //�ж��ھӡ�Ȩ���Ƿ���ͬ
	//set_intersection(i_adj.begin(), i_adj.end(), v_adj.begin(), v_adj.end(), inserter(result, result.begin()));
	//if (result.size() != i_adj.size()) return 0;
	//return 1;
}

void mergeIdent1(int i, int v) {  //ע��i��v��reachֵ+��������identֵ����
	//add pair dependency
	//��i��v��pair dependency
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
	//���ھӶ����pair dependency
	int numpathu_i = 0;
	map<int,int> neighborlist = Pmg.adjlist[i];
	int neighbornum = neighborlist.size();
	/*vector<int> weighti, weightv;*/  //��¼i��v���ھ�֮���Ȩ��:i\v�������ھ�w��Ȩ�����
	/*for (auto w : neighborlist) weighti.push_back(w.second); */   //???????????????????????????????????????????
	/*for (int w : neighborlist) weightv.push_back(Pmg.Mp[v][w]);*/
	for (auto w:neighborlist) {
		numpathu_i += w.second*w.second;
	}
	for (auto w:neighborlist) {
		bc[w.first] += reachiident * ((double)w.second * w.second / numpathu_i) * reachvident;  //��Ϊ�ǶԳƵģ�����numpath_w_v=numpath_i_w=weightx[wk]
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

void mergeIdent2(int i, int v) {  //ע��i��v��reachֵ
	//add pair dependency v��i+identSet[i]
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

/*��sΪԴ����п���+�����ۻ�������Դ��s��Դ����*/
void BBC1(int s) {  //+reachֵ��
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
		cout << "����" << i << "���·������" << pathnum[i] << endl;
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
				cout << "�Ե�" << w << "Դ����:" << deta[w] * reach[s] << endl;
			}*/
		}
	}

}

void BBC2(int s) {  //+reachֵ+ident��
	//initialize variable
	/*cout << "Դ��" << s << endl;*/
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
		cout << "����" << i << "���·������" << pathnum[i] << endl;
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
			//��w��������ident����ʱ��w��w��������ж����ǰ�ö���Ĺ���
			double deta_add = (pathnum[v.first] * v.second / pathnum[w]) * (1 + deta[w]);
			if (ident[w] > 0) {
				for (int wident : IdentSet[w]) {
					deta_add += (pathnum[v.first] * v.second / pathnum[w]) * (1 + deta[wident]);  //deta[w]��deta[wident]��ͬ����Ϊw��wident�����в�ͬ��reachֵ
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
				cout << "�Ե�" << w << "Դ����:" << deta[w] * reach[s] << endl;
			}*/
		}
	}
	stack<int>().swap(S);
	vector<vector<pair<int, int>>>().swap(pred);
	vector<double>().swap(pathnum);
	vector<double>().swap(deta);

}

void find_side0(map<int, map<int,int>>* Mpl) {  //Mpl���ڽӱ�����i�Ķ�����ڽӱ�����sizeΪ1����õ�Ϊside
	flag_side.resize(Pmg.n_vertex);
	auto iteri = Mpl->begin();
	for (iteri; iteri != Mpl->end(); iteri++) {
		if (iteri->second.size() == 1) {  //���i���ڽӵ�ֻ��һ��
			int i = iteri->first;
			auto iterj = iteri->second.begin();
			int irelatedj = iterj->first;
			SideSet0[irelatedj].push_back(i);
			flag_side[i] = 1;
		}
	}
	//for (int i = 0, lenr = Mpl->size(); i < lenr; i++) {  //A1�ඥ��i
	//	
	//	int counti = 0;
	//	int irelatedj;
	//	int j, lenc;
	//	for (j = 0, lenc = (*Mpl)[0].size(); j < lenc; j++) {  //AL�ඥ��j
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

/*���ļ�����ȡ�칹ͼ��Ϣ*/
void hetergraph::read_file() {
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
		cout << "��ʼ�����ļ�" << edge_filename << endl;
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

/*���ݱߵ���������������ţ����ұߵ�������*/
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
void multigraph::getPmg(vector<char>* P) {
	vector<char>::iterator it = P->begin();  //it��������Ԫ·��P�еĶ������
	/*Pmg������*/
	n_vertex = hg.vertex_num[hg.vertex_class_index[*it]];
	n_vertex_org = n_vertex;
	adjlist.resize(n_vertex);

	/*����P����P�бߵ�������������eindex��*/
	vector<int> eindex;
	for (int i = 0; i < P->size() / 2; i++) {  //һ������P->size()/2����
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

	/*����Mpl���Ҷ���ڵ�side����*/
	double ps13, pt13, time13;
	ps13 = clock();
	find_side0(&matrixpl);
	pt13 = clock();
	time13 = (double)(pt13 - ps13) / CLOCKS_PER_SEC;
	cout << "ʶ��side���㼯����ʱ��" << time13 << endl;
	out_result << "ʶ��side���㼯����ʱ��" << time13 << endl;

	/*����Mp+�ڽӱ�*/
	cout << "����adjlist" << endl;
	double ps12, pt12, time12;
	ps12 = clock();
	m_edge = 0;

	//�õ�Mpl�ڽӱ��ת��
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

	for (int i = 0; i < n_vertex; i++) {  //�����i������ڽ�����
		/*cout << "���㵽��" << i;*/
		auto iterik = matrixpl.find(i);  //�ҵ�i�����k�ڽ�����
		if (iterik != matrixpl.end()) {  //�ҵ���
			auto iterk = iterik->second.begin();
			for (iterk; iterk != iterik->second.end(); iterk++) {  //����i��ÿ���ڽӵ�k
				int k = iterk->first;
				int currik = iterk->second;
				auto iterkj = matrixplt.find(k);//��k���ڽӵ�j
				if (iterkj != matrixplt.end()) {
					auto iterj = iterkj->second.begin();
					for (iterj; iterj != iterkj->second.end(); iterj++) {
						int j = iterj->first;
						int currkj = iterj->second;
						if (j > i) {
							//��i��adjlist�У�����j�������Ѿ��������Ҳ����û�������;ͬʱ��j��adjlist�в���i
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
	cout << "��mpl�õ�mp��ʱ��" << time12 << endl;
	out_result << "��mpl�õ�mp��ʱ��" << time12 << endl;
	/*save_Pmg();*/

}

/*��P-multigraph��ɾ����*/
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
/*��P-multigraph��ɾ����*/
void multigraph::deleteVertex(int a) {
	//��a��ÿ���ڽӵ�
	for (auto an : adjlist[a]) {
		//�޸��ڽӾ���
		/*Mp[a][an] = 0; Mp[an][a] = 0;*/
		//��an���ڽӱ���ɾ��a
		map<int, int>::iterator it = adjlist[an.first].begin();
		for (it; it != adjlist[an.first].end(); it++) {
			if (it->first == a) {
				adjlist[an.first].erase(it);
				break;
			}
		}
	}
	//���a���ڽӱ���ʱaΪP-multigraph��һ�������ĵ㣬��ͨ��index������Ȼ�ܱ�����a��
	map<int, int>().swap(adjlist[a]);
}