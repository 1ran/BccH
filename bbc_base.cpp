//9.20update
//9.21update �ڶ���ͼhgʱͨ���ڽӱ���ʽ���ڼ���pmgʱ��ͨ���ڽӱ���㣻pmg��ɾ���ڽӾ������ڽӱ�洢��Ȩ
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
	int n_vertexclass;  //���������
	map<char, int> vertex_class_index;  //�������-������
	vector<int> vertex_num;  //����������˳��洢ÿ�ඥ����

	int n_edgeclass;  //�ߵ������
	vector<vector<int>> edge_class;  //ÿ������ӵĶ�����𣺶�����x��������y
	vector<int> edge_num;  //ÿ�����

	//vector<vector<vector<int>>> edge_info;  //3ά���飬��һάedge_num����ʾÿ��ߣ��ڶ�ά�͵���ά��ʾ��edge���ڽӾ���
	vector<map<int,map<int,int>>> edge_info;  //3ά���飬��һάedge_num����ʾÿ��ߣ��ڶ�ά�͵���ά��ʾ��edge���ڽӱ�

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
	vector<map<int,int> > adjlist;
	void show_mg();
	void getPmg(vector<char>* P);  //���칹ͼ�л��Pmg
	void getPmg_read_file();      //�ӱ�����ļ��ж�ȡMpl��Pmg
	void save_Mpl(vector<vector<int>>* Mpl);  //����Mpl
	void save_Pmg();  //��������Pmg
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
ofstream out_result(file + "/al1_time_result.txt");  //��ʾ�㷨1��ʱ����
double timebfs = 0, timeback = 0;
double pmid;


int main() {
	//һ�������칹ͼ����Pmg
	double ps00, pt00, time00;
	ps00 = clock();
	hg.read_file();
	pt00 = clock();
	time00 = (double(pt00 - ps00)) / CLOCKS_PER_SEC;
	out_result << "�����칹ͼ��ʱ��" << time00 << endl;

	/*hg.show_hg();*/
	double ps0, pt0, time0;
	ps0 = clock();
	vector<char> P = { 'A','M','D','M','A' };
	Pmg.getPmg(&P);  //�õ�P-multigraph
	pt0 = clock();
	time0 = (double)(pt0 - ps0) / CLOCKS_PER_SEC;
	cout << "�õ�Pmg��ʱ��" << time0 << endl;
	out_result << "�õ�Pmg��ʱ:" << time0 << endl;
	cout << "������" << Pmg.n_vertex << endl;
	out_result << "Pmg������" << Pmg.n_vertex << endl;
	cout << "������" << Pmg.m_edge << endl;
	out_result << "Pmg������" << Pmg.m_edge << endl;
	//�����Ѿ�������Pmg�����ļ��ж�ȡ
	/*Pmg.getPmg_read_file();*/

	/*ͳ�������ͨ������С*/
	double pcmps, pcmpt, timecmp;
	pcmps = clock();
	int maxComponent=getMaxComponent();
	cout << "�����ͨ������" << maxComponent << endl;
	out_result << "Pmg�����ͨ������" << maxComponent << endl;
	pcmpt = clock();
	timecmp = (double)(pcmpt - pcmps) / CLOCKS_PER_SEC;
	out_result << "���������ͨ������С��ʱ��" << timecmp << endl;


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
		if(count%100==0) cout << count << "���������Ҫ" << time1 << endl;
	}
	cout << "ǰ��bfs��ʱ�䣺" << timebfs << endl;
	cout << "�����ۻ���ʱ�䣺" << timeback << endl;
	cout << "����bcֵ��ʱ�䣺"<<time1 << endl;
	out_result << "ǰ��bfs��ʱ�䣺" << timebfs << endl;
	out_result << "�����ۻ���ʱ�䣺" << timeback << endl;
	out_result << "����bcֵ��ʱ�䣺" << time1 << endl;

	double time_all = (double)(pt1 - ps00) / CLOCKS_PER_SEC;
	cout << "����ʱ��" << time_all << endl;
	out_result << "����ʱ��" << time_all << endl;
	cout << "����ʱ���Ѵ���" << file + "/time_result.txt" << endl;
	save_bc();
	compare_result();
	out_result.close();
	return 0;
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
void save_bc() {
	ofstream out(file + "/bbc_result.txt");
	for (int i = 0; i < Pmg.n_vertex; i++) {
		out << i << " ";
		out << fixed << setprecision(8) << bc[i] << endl;
	}
	cout << "��bbc�������bcֵ�ѱ�����" << file << "/result.txt"<<endl;
}

void BBC(int s) {
	//initialize variable
	/*cout << "Դ��" << s << endl;*/
	queue<int> Q;
	stack<int> S;
	vector<vector<pair<int,int>>> pred(Pmg.n_vertex);  //��map
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
		for (auto w : Pmg.adjlist[v]){  //������v�������ھ�
			if (dist[w.first] == DBL_MAX) {  
				dist[w.first] = dist[v] + 1;
				Q.push(w.first);
			}
			if (dist[w.first] == dist[v] + 1) {  //ͨ����if���������Եõ�v�������·����ھ�w���ڴ�ͳ���·�w�����ۻ���vʱ��Ȩֵ
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
		int nx = vertex_num[vx], ny = vertex_num[vy];  //nx������ ny������
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
				xadjlist.insert(map<int, int>::value_type(y,1));  //��y����x���ڽ�list
				edge_matrix.insert(map<int, map<int, int>>::value_type(x, xadjlist)); //��x���ڽ�list����edge_matrix
				xadjlist.clear();
			}
			else {  //x�Ѿ�����
				iterx->second.insert(map<int, int>::value_type(y, 1));
			}
		}
		edge_info.push_back(edge_matrix);  //��i��ߵ��ڽӱ�
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
	map<int,map<int,int>> matrixpl = hg.edge_info[eindex[p]];  //��1��ߵ��ڽӾ���
	hg.edge_info[eindex[p]].clear();  //��ʱ����һ��eindex���ڽӾ����Ѿ�������matrixpl�У���Ϣ������Ҫ
	while (++p < eindex.size()) {  //��matrix����һ���ڽӾ���edge_info[eindex[p]]���
		int curr_eindex = eindex[p];  //��ǰ����ıߵ�����,����ǰ���ĸ��ߵ��ڽӾ�����MPL���
		cout << "����"<<eindex[p-1]<<'*' << curr_eindex << "���ڽӾ���" << endl;
		double ps11, pt11, time11;
		ps11 = clock();
		map<int,map<int,int>> tmatrix;  //��ʱ��¼������
		auto iterik = matrixpl.begin();
		for(iterik;iterik!=matrixpl.end();iterik++) {  //��ÿ�����ڽ����е�i
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
		cout << " ��ʱ��" << time11;
		cout << endl;
		out_result << "����" << eindex[p - 1] << '*' << curr_eindex << "���ڽӾ���" << " ��ʱ��" << time11 << endl;
	}
	//save_Mpl(&matrixpl);  //����Mpl
	eindex.clear();
	/*for (auto it : matrixpl) {
		cout << it.first << ":";
		for (auto itt : it.second) {
			cout << itt.first << "(" << itt.second << ")  ";
		}
		cout << endl;
	}*/

	/*�����Mpl�󣬼���Mp+side*/
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
	double pt13 = clock();
	double time13 = (double)(pt13 - ps12) / CLOCKS_PER_SEC;
	out_result << "����mpl��ת����ʱ��" << time13 << endl;

	for (int i = 0; i < n_vertex; i++) {  //�����i������ڽ�����
		cout << "���㵽��" << i;
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
	out_result<< "��mpl�õ�mp��ʱ��" << time12 << endl;
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

