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
#include<cmath>
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
	vector<map<int,int> > adjlist_edge;  //��ÿһ�д洢��������ڶ����ͬʱ���洢ͨ�������������ڶ�������
	//vector<map<int,int>> adjlist_edge;  //��ÿһ�д洢��������ڶ����ͬʱ���洢ͨ�������������ڶ�������
	//vector<vector<int> > Mp;  //??
	/*vector<vector<int> > adjlist;*/

	void show_mg();
	void getPmg(vector<char>* P);  //���칹ͼ�л��Pmg
	void getPmg_read_file();      //�ӱ�����ļ��ж�ȡMpl��Pmg
	void save_Mpl(vector<vector<int>>* Mpl);  //����Mpl
	void save_Pmg();  //��������Pmg
};

void CBC(int s);  //����v���·��ʵ����Ϣ
void CBC_(int s);  //ֱ����pred�б����Ȩ��pred��vector  ����졿
void CBC_1(int s);  //ֱ����pred�б����Ȩ��pred��map 
void save_bc();
void getbcfromfile();
void compare_result();
void compare_bc();
int getMaxComponent();
int getComponentSize(int i, vector<int>* visited);
void compare_order();

string file = "data1/delete_top";
hetergraph hg;
multigraph Pmg;
vector<double> bc;
//vector<double> deta0;
vector<double> cbc;
vector<double> bbc;
ofstream out_result(file + "/al1_time_result.txt");
double timebfs = 0, timeback = 0;
double pmid;
double weight_pathinstance_min = 100000;
double weight_pathinstance_max = 0;
double weight_D_min = 100000;
double weight_D_max = 0;
double weight_min = 100000;
double weight_max = 0;
double p_type = 0;

int main() {
	//һ�������칹ͼ����Pmg
	double ps00, pt00, time00;
	ps00 = clock();
	hg.read_file();
	pt00 = clock();
	time00 = (double(pt00 - ps00)) / CLOCKS_PER_SEC;
	out_result << "�����칹ͼ��ʱ��" << time00 << endl;
	/*hg.show_hg();*/
	//vector<char> P = { 'A','P','A' };
	vector<char> P = { 'A','M','D','M','A' };
	//vector<char> P = { 'B','R','U','R','B' };
	double ps0, pt0, time0;
	ps0 = clock();
	Pmg.getPmg(&P);  //�õ�P-multigraph
	pt0 = clock();
	time0 = (double)(pt0 - ps0) / CLOCKS_PER_SEC;
	////Pmg.show_mg();
	//�����Ѿ�������Pmg�����ļ��ж�ȡ
	/*Pmg.getPmg_read_file();*/
	/*Pmg.show_mg();*/
	cout << "�õ�Pmg��ʱ��" << time0 << endl;
	out_result << "�õ�Pmg��ʱ:" << time0 << endl;
	cout << "������" << Pmg.n_vertex << endl;
	out_result << "Pmg������" << Pmg.n_vertex << endl;
	cout << "������" << Pmg.m_edge << endl;
	out_result << "Pmg������" << Pmg.m_edge << endl;

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
	/*deta0.resize(Pmg.n_vertex);*/
	double ps1, pt1, time1 = 0;
	int count = 0;
	for (int i = 0; i < Pmg.n_vertex; i++) {
		ps1 = clock();
		CBC_(i);
		pt1 = clock();
		timebfs += (double)(pmid - ps1) / CLOCKS_PER_SEC;
		/*cout << "bfs:" << (double)(pmid - ps1) / CLOCKS_PER_SEC<<endl;*/
		timeback += (double)(pt1 - pmid) / CLOCKS_PER_SEC;
		/*cout<<"back:"<< (double)(pt1 - pmid) / CLOCKS_PER_SEC<<endl;*/
		time1 += (double)(pt1 - ps1) / CLOCKS_PER_SEC;
		count++;
		if (count % 100 == 0) {
			cout << count << endl;
		}
	}
	cout << "ǰ��bfs��ʱ�䣺" << timebfs << endl;
	cout << "�����ۻ���ʱ�䣺" << timeback << endl;
	cout << "����bcֵ��ʱ�䣺" << time1 << endl;
	out_result << "ǰ��bfs��ʱ�䣺" << timebfs << endl;
	out_result << "�����ۻ���ʱ�䣺" << timeback << endl;
	out_result << "����bcֵ��ʱ�䣺" << time1 << endl;

	//for (int i = 0; i < Pmg.n_vertex; i++) {
	//	/*if (bc[i] > 0) {
	//		cout << "����" << i << "��bcֵ��" << bc[i] << endl;
	//	}*/
	//	cout << "����" << i << "��bcֵ��" << bc[i] << endl;
	//}
	double time_all = (double)(pt1 - ps00) / CLOCKS_PER_SEC;
	cout << "����ʱ��" << time_all << endl;
	out_result << "����ʱ��" << time_all << endl;
	cout << "����ʱ���ѱ���"  << endl;
	/*compare_result();*/
	
	cout << endl;
	cout << "Ȩ�����" << weight_max << endl;
	cout << "Ȩ����С��" << weight_min << endl;
	cout << "Ȩ��D���" << weight_D_max << endl;
	cout << "Ȩ��D��С��" << weight_D_min << endl;
	cout << "Ȩ��pathinstance���" << weight_pathinstance_max << endl;
	cout << "Ȩ��pathinstance��С��" << weight_pathinstance_min << endl;
	out_result << "Ȩ�����" << fixed << setprecision(8)<<weight_max << endl;
	out_result << "Ȩ����С��" << fixed << setprecision(8) << weight_min << endl;
	out_result << "Ȩ��D���" << fixed << setprecision(8) << weight_D_max << endl;
	out_result << "Ȩ��D��С��" << fixed << setprecision(8) << weight_D_min << endl;
	out_result << "Ȩ��pathinstance���" << fixed << setprecision(8) << weight_pathinstance_max << endl;
	out_result << "Ȩ��pathinstance��С��" << fixed << setprecision(8) << weight_pathinstance_min << endl;

	out_result.close();

	save_bc();

	return 0;
	
}


void save_bc() {
	ofstream out(file + "/cbc_result_amdma.txt");
	for (int i = 0; i < Pmg.n_vertex; i++) {
		out << i << " ";
		out << fixed << setprecision(8) << bc[i] << endl;
	}
	cout << "��kbc�������bcֵ�ѱ�����" << file << "/cbc_result_amdma.txt";
}

//int main() {
//	getbcfromfile();
//	compare_bc();
//	return 0;
//}


void getbcfromfile() {
	//�������ж����cbcֵ
	string filename = file + "/cbc_result_amdma.txt";  //base�д洢�ļ��Ļ�����Ϣ�����ݱ������ж�Ӧ�ıߵ��������ҵ��ߵ��ļ�"x.txt"
	ifstream inputcbc(filename, ios::in);  //�����ļ�������input
	//�ж��ļ��Ƿ����
	if (!inputcbc) {
		cerr << "file error!" << endl;
		exit(1);
	}

	//�������ж����cbcֵ//����ÿ�������bcֵ���ļ�
	string cbc_info;
	while (getline(inputcbc, cbc_info)) {
		istringstream scbc_info(cbc_info);
		int v;
		double bcvalue;
		scbc_info >> v >> bcvalue;
		cbc.push_back(bcvalue);
	}
	inputcbc.close();

	//�������ж����bbcֵ
	string filename2 = file + "/bbc_result_amdma.txt";  //base�д洢�ļ��Ļ�����Ϣ�����ݱ������ж�Ӧ�ıߵ��������ҵ��ߵ��ļ�"x.txt"
	ifstream inputbbc(filename2, ios::in);  //�����ļ�������input
	//�ж��ļ��Ƿ����
	if (!inputbbc) {
		cerr << "file error!" << endl;
		exit(1);
	}

	//�������ж����bbcֵ//����ÿ�������bcֵ���ļ�
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

bool cmp(const pair<int, double>& a, const pair<int, double>& b) {
	return a.second > b.second;
}
vector<pair<int, double>> cbc_;
vector<pair<int, double>> bbc_;
void compare_bc() {
	for (int i = 0; i < cbc.size(); i++) {
		cbc_.push_back({ i,cbc[i] });  /*�ڵ���-bcֵ*/
	}
	sort(cbc_.begin(), cbc_.end(), cmp);

	for (int i = 0; i < bbc.size(); i++) {
		bbc_.push_back({ i,bbc[i] });
	}
	sort(bbc_.begin(), bbc_.end(), cmp);
	
	/*��������Ľ��*/
	ofstream outcbc(file + "/cbc_order_amdma.txt");
	if (!outcbc) {
		cerr << "file cbcorder.txt error" << endl;
		exit(1);
	}
	for (int i = 0; i < cbc_.size(); i++) {
		outcbc << cbc_[i].first << " " << cbc_[i].second << endl;
	}
	outcbc.close();
	cout << "������cbcֵ�Ѽ�¼" << endl;
	ofstream outbbc(file + "/bbc_order_amdma.txt");
	if (!outbbc) {
		cerr << "file bbcorder.txt error" << endl;
		exit(1);
	}
	for (int i = 0; i < bbc_.size(); i++) {
		outbbc << bbc_[i].first << " " << bbc_[i].second << endl;
	}
	outbbc.close();
	cout << "������bbcֵ�Ѽ�¼" << endl;
	/*���ǰtopn���Ľ����*/
	int topn = 50;
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
	ofstream outcompare(file + "/comparebbccbc_amdma.txt");
	if (!outcompare) {
		cerr << "file compare.txt error" << endl;
		exit(1);
	}
	outcompare << "������" << result_inter.size()<< endl;
	for (int i : result_inter) {
		outcompare << i << " ";
	}
	outcompare << endl;
	outcompare << "c2-b���" << result_diffc_b.size() << endl;
	for (int i : result_diffc_b) {
		outcompare << i << " ";
	}
	outcompare << endl;
	outcompare << "b-c2���" << result_diffb_c.size()<< endl;
	for (int i : result_diffb_c) {
		outcompare << i << " "; 
	}
	outcompare << endl;
	outcompare.close();

	compare_order();
	
}

void compare_order() {   
	ofstream outorder(file + "/deta_order_bbc_cbc_amdma.txt");
	if (!outorder) {
		cerr << "file deta_order.txt error" << endl;
		exit(1);
	}
	/*��¼�����ÿ�������bbc˳��*/
	map<int, int> bbc_id_order;  //������-bbc˳��
	for (int i = 0; i < bbc_.size(); i++) {
		bbc_id_order[bbc_[i].first] = i + 1;
	}
	/*��˳�����cbc������order���:kbc��bbc������ǰ����*/
	for (int i = 0; i < cbc_.size(); i++) {
		int nodeid = cbc_[i].first;
		int order_deta = bbc_id_order[nodeid] - (i + 1);
		outorder << nodeid << " " << order_deta << endl;
	}
	outorder.close();
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
//	/*cout << "Դ��" << s << endl;*/
//	queue<int> Q;
//	stack<int> S;
//	vector<vector<int>> pred(Pmg.n_vertex);
//	vector<double> dist(Pmg.n_vertex,DBL_MAX);
//	vector<double> pathnum(Pmg.n_vertex);
//	vector<map<int, int>> pathinstance(Pmg.n_vertex);  //ÿ�������·�·��ʵ��������-��Ŀ��Ϣ�����ڼ����ĸ
//	dist[s] = 0;
//	pathnum[s] = 1;
//	Q.push(s);
//	//bfs
//	while (!Q.empty()) {
//		int v = Q.front();
//		Q.pop();
//		S.push(v);
//		for (auto w : Pmg.adjlist_edge[v]) {  //������v�������ھ�w[�ڽӶ��㣬�ڽӱ�]
//			if (dist[w[0]] == DBL_MAX) {  
//				dist[w[0]] = dist[v] + 1;
//				Q.push(w[0]);
//			}
//			if (dist[w[0]] == dist[v] + 1) {  //ͨ����if���������Եõ�v�������·����ھ�w
//				/*pathnum[w[0]] += pathnum[v]*Pmg.Mp[v][w[0]];*/
//				pred[w[0]].push_back(v);    //����w[0]��ǰ�ö��㣬����v
//				//ͳ��v�·���·��ʵ����Ϣ
//				vector<vector<int>> einfo = Pmg.edge_info[w[1]];
//				for (auto pinfo : einfo) {
//					int pclass = pinfo[0];
//					int pnum = pinfo[1];
//					pathinstance[v][pclass] += pnum;
//				}
//			}
//		}
//		//������v���ھӵ�󣬲���ͳ��������v��·��ʵ������Ϣ����ʱ���ܸ���pathnum
//		for (auto w : Pmg.adjlist_edge[v]) {
//			if (dist[w[0]] == dist[v] + 1) { 
//				vector<vector<int>> einfo = Pmg.edge_info[w[1]];  //w[1]��v��w[0]֮��ı�
//				for (auto pinfo : einfo) {  //����v��w[0]֮��ߵ�·��ʵ����Ϣ������pathnum[w[0]]
//					int pclass = pinfo[0];
//					int pnum = pinfo[1];
//					pathnum[w[0]] += pathnum[v] * ((double)pnum / pathinstance[v][pclass]);
//				}
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
//		for (int v : pred[w]) {
//			//����v��w֮��ߵ�·��ʵ������Ϣ��ÿ��·��ʵ���ֱ���з����ۻ�
//			int eindex = Pmg.Mp_edge[v][w];  //v��w֮��ߵ�������
//			vector<vector<int>> einfo = Pmg.edge_info[eindex];  //�ñߵ���Ϣ
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
//		cout << "�Ե�" << i << "Դ������" << deta[i] << endl;
//	}*/
//}

void CBC_(int s) {  //pred��vector
	//initialize variable
	/*cout << "Դ��" << s << endl;*/
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
						pathinstance[pclass] += pnum;
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
						if (weight_pathinstance_max < weight) weight_pathinstance_max = weight;
						if (weight_pathinstance_min > weight) weight_pathinstance_min = weight;
						if (weight_D_max < einfo.size()) weight_D_max = einfo.size();
						if (weight_D_min > einfo.size()) weight_D_min = einfo.size();
					    weight = weight + (double)einfo.size();//1
						if (weight_max < weight) weight_max = weight;
						if (weight_min > weight) weight_min = weight;
						//weight = (double)einfo.size();//1
						pathnum[w.first] += pathnum[rv] * weight;
						pred[w.first].push_back({ (double)rv,weight });
					}
				}
			}
			level = level + 1;  //���������ʵ�level+1��
			pathinstance.clear(); //����ͳ�Ƶ�level+1���pathinstance��Ϣ
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
		}
		if (w != s) {
			bc[w] += deta[w];
		}
	}

}

//void CBC_1(int s) {  //pred��map
//	//initialize variable
//	/*cout << "Դ��" << s << endl;*/
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
//		for (auto w : Pmg.adjlist_edge[v]) {  //������v�������ھ�w[�ڽӶ��㣬�ڽӱ�]
//			if (dist[w[0]] == DBL_MAX) {
//				dist[w[0]] = dist[v] + 1;
//				Q.push(w[0]);
//			}
//			if (dist[w[0]] == dist[v] + 1) {  //ͨ����if���������Եõ�v�������·����ھ�w
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
//					weight += ((double)pnum / pathinstance[pclass]);
//				}
//				pathnum[w[0]] += pathnum[v] * weight;
//				pred[w[0]][v]=weight;    //����w[0]��ǰ�ö��㣬����v��v-w[0]֮��ı�Ȩ
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
//		for (auto v : pred[w]) {  //v��ʾ w��ǰ�ö���v.first-w��v֮��ı�Ȩv.second
//			deta[v.first]+= (pathnum[v.first] * v.second / pathnum[w]) * (1 + deta[w]);
//		}
//		if (w != s) {
//			bc[w] += deta[w];
//		}
//	}
//	/*for (int i = 0; i < deta.size(); i++) {
//		cout << "�Ե�" << i << "Դ������" << deta[i] << endl;
//	}*/
//
//}





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
void multigraph::getPmg(vector<char>* P) {
	vector<char>::iterator it = P->begin();  //it��������Ԫ·��P�еĶ������
	/*Pmg������*/
	n_vertex = hg.vertex_num[hg.vertex_class_index[*it]];
	n_vertex_org = n_vertex;
	adjlist_edge.resize(n_vertex);  //n��ai������ڽӶ���+�ߵı��

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

	/*ofstream outnum(file + "/numberofauthor.txt");
	if (!outnum) {
		cerr << "file numberofauthor.txt error" << endl;
		exit(1);
	}
	cout << "������bbcֵ�Ѽ�¼" << endl;
	for (int j = 0; j < matrixpl[0].size(); j++) {
		int count = 0;
		for (int i = 0; i < matrixpl.size(); i++) {
			if (matrixpl[i][j] == 1) {
				count++;
			}
		}
		outnum  << j << " " << count << endl;
	}
	outnum.close();*/

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

	cout << "����Mp�ڽӱ�" << endl;
	int edge_index = 0;
	for (int i = 0; i < n_vertex; i++) {  //�����i������ڽ�����
	/*cout << "���㵽��" << i;*/
		auto iterik = matrixpl.find(i);  //�ҵ�i�����k�ڽ�����
		if (iterik != matrixpl.end()) {  //�ҵ���
			map<int,vector<vector<int>>> ij_einfo;  //��¼i�������ھ�j����֮��ߵ���Ϣ
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

