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
	vector<int> edge_num;  //ÿ�����

	vector<map<int, map<int, int>>> edge_info;  

	void read_file();
	void show_hg();
	int getEdgeIndex(int x, int y);  
};


class multigraph {
public:
	int n_vertex;  //���¶�����������copy�Ĺؽڶ��㣩
	int n_vertex_org;  //ԭ������
	int m_edge;  //����
	vector<vector<vector<int>>> edge_info;  
	vector<map<int,int> > adjlist_edge;  

	void getPmg(vector<char>* P);  //���칹ͼ�л��Pmg
	void save_Mpl(vector<vector<int>>* Mpl);  //����Mpl
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
	//һ�������칹ͼ����Pmg
	hg.read_file();

	/*hg.show_hg();*/
	vector<char> P = { 'A','M','D','M','A' };
	Pmg.getPmg(&P);  //�õ�P-multigraph


	/*ͳ�������ͨ������С*/
	int maxComponent = getMaxComponent();

	bc.resize(Pmg.n_vertex);
	int count = 0;
	for (int i = 0; i < Pmg.n_vertex; i++) {
		CBC_(i);
		count++;
	}

	for (int i = 0; i < Pmg.n_vertex; i++) {
		/*if (bc[i] > 0) {
			cout << "����" << i << "��bcֵ��" << bc[i] << endl;
		}*/
		cout << "����" << i << "��bcֵ��" << bc[i] << endl;
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
	
	/*��������Ľ��*/
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
	/*���ǰtopn���Ľ����*/
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
	outcompare << "������" << endl;
	for (int i : result_inter) {
		outcompare << i << " ";
	}
	outcompare << endl;
	outcompare << "c-b���" << endl;
	for (int i : result_diffc_b) {
		outcompare << i << " ";
	}
	outcompare << endl;
	outcompare << "b-c���" << endl;
	for (int i : result_diffb_c) {
		outcompare << i << " ";
	}
	outcompare << endl;
	outcompare.close();
}

void getbcfromfile() {
	string filename = file + "/cbc_result2.txt";  
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

	string filename2 = file + "/bbc_result.txt";  
	ifstream inputbbc(filename2, ios::in);  //�����ļ�������input
	//�ж��ļ��Ƿ����
	if (!inputbbc) {
		cerr << "file error!" << endl;
		exit(1);
	}

	//�������ж����cbcֵ//����ÿ�������bcֵ���ļ�
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
	ifstream input(filename, ios::in);  //�����ļ�������input
	//�ж��ļ��Ƿ����
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
		map<int, int> pathinstance;  //w�·�����·��ʵ������ͳ��
		for (auto w : Pmg.adjlist_edge[v]) {  //������v�������ھ�w[�ڽӶ��㣬�ڽӱ�]
			if (dist[w.first] == DBL_MAX) {
				dist[w.first] = dist[v] + 1;
				Q.push(w.first);
			}
			if (dist[w.first] == dist[v] + 1) {  //ͨ����if���������Եõ�v�������·����ھ�w
				/*pathnum[w[0]] += pathnum[v]*Pmg.Mp[v][w[0]];*/
				//ͳ��v�·���·��ʵ����Ϣ
				vector<vector<int>> einfo = Pmg.edge_info[w.second];
				for (auto pinfo : einfo) {
					int pclass = pinfo[0];
					int pnum = pinfo[1];
					pathinstance[pclass] += pnum;
				}
			}
		}
		//������v���ھӵ�󣬲���ͳ��������v��·��ʵ������Ϣ�����v����ÿ�����õ�w֮��ߵ�Ȩ��
		for (auto w : Pmg.adjlist_edge[v]) {
			double weight;
			weight = 0;
			if (dist[w.first] == dist[v] + 1) {  //����v��w[0]֮��ı�Ȩ
				vector<vector<int>> einfo = Pmg.edge_info[w.second];  //w[1]��v��w[0]֮��ı�
				for (auto pinfo : einfo) {  //����v��w[0]֮��ߵ�·��ʵ����Ϣ������pathnum[w[0]]
					int pclass = pinfo[0];
					int pnum = pinfo[1];
					weight += ((double)pnum / pathinstance[pclass]);  //�ۼ�v-w[0]֮��ı�Ȩ
					/*if (s == 1) {
						if (v == 2640 && pclass == 470) {
							cout << "2640�·���" << w[0] << "��" << pnum << "��" << pclass << "��·��ʵ��" << endl;
							cout << pathinstance[pclass] << "pathinstance[pclass] +=" << pnum << endl;
						}
					}*/
				}
				//weight = weight * einfo.size();  //�����
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
				cout << "w=" << w << "ʱ��" << "deta[(int)" << 19 << "] += (" << pathnum[(int)v[0]] << " * " << v[1] << " / " << pathnum[w] << ") * (1 + " << deta[w] << ") "<< endl;
			}*/
			/*if (s == 1) {
				if ((int)v[0] == 0) {
					cout << "��" << w << "����deta[0]:";
					cout << "deta[" << (int)v[0] << "] += (" << pathnum[(int)v[0]] << " * " << v[1] << " / " << pathnum[w] << ") * (1 + " << deta[w] << ") " << endl;
				}
			}*/
		}
		if (w != s) {
			bc[w] += deta[w];
		}
	}
	/*cout << "��"<<s<<"�Ե�19��Դ������" << deta[19] << endl;*/
	/*for (int i = 0; i < deta.size(); i++) {s
		cout << "�Ե�" << i << "Դ������" << deta[i] << endl;
	}*/
	/*deta0[s] = deta[0];*/
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




void save_bc() {
	ofstream out(file + "/cbc_result.txt");
	for (int i = 0; i < Pmg.n_vertex; i++) {
		out << i << " ";
		out<< fixed << setprecision(8) << bc[i] << endl;
	}
	cout << "bcֵ�ѱ�����" << file << "/result.txt";
}

/*���ļ�����ȡ�칹ͼ��Ϣ*/
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
			else {  //x�Ѿ�����
				iterx->second.insert(map<int, int>::value_type(y, 1));
			}
		}
		edge_info.push_back(edge_matrix);  
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
	vector<char>::iterator it = P->begin();  
	/*Pmg������*/
	n_vertex = hg.vertex_num[hg.vertex_class_index[*it]];
	n_vertex_org = n_vertex;
	adjlist_edge.resize(n_vertex);  

	/*����P����P�бߵ�������������eindex��*/
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
	cout << "�ҵ���������";
	for (int i : eindex) cout << i << "-";
	cout << endl;
	/*����Mpl*/
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

	int edge_index = 0;
	for (int i = 0; i < n_vertex; i++) {  //
	/*cout << "���㵽��" << i;*/
		auto iterik = matrixpl.find(i);  
		if (iterik != matrixpl.end()) {  //
			map<int,vector<vector<int>>> ij_einfo;
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
							auto iterjeinfo = ij_einfo.find(j);  
							if (iterjeinfo == ij_einfo.end()) { 
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

