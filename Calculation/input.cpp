#include "input.hpp"
#include "param.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <filesystem>
#include <vector>
#include <direct.h>
#include <sys/stat.h>
using namespace std;
namespace fs = std::filesystem;
InputData::InputData()
{
	setgoal();
	input_param();
	input_csv();
	input_mesh();
}
vector<double> InputData::getUx() {
	return Ux;
}
vector<double> InputData::getUy() {
	return Uy;
}
vector<double> InputData::getP() {
	return P;
}
vector<int> InputData::getcond() {
	return cond;
}
vector<vector<int>> InputData::getnbool1() {
	return nbool1;
}
vector<vector<int>> InputData::getnbool3() {
	return nbool3;
}
NodeP InputData::get_NodeParam() {
	return nparam;
}
Boundarycond InputData::get_BC() {
	return BC;
}
string InputData::getmodel() {
	return modelpath;
}
string InputData::getdate() {
	return datepath;
}
string InputData::getdata() {
	return datapath;
}
string InputData::getgoal() {
	return goalpath;
}
void InputData::setmodel(){
	cout << "Mesh_data_box内部" << endl;
	viewdir(start);
	cout << "選択モデル->";
	cin >> model;
	modelpath = start + "/" + model;
	struct stat statBuf;
	if (stat(modelpath.c_str(), &statBuf) != 0) {
		
	}
	else {
		cout << modelpath << "は存在しません" << endl;
		exit(-1);
	}
}
void InputData::setdate() {
	setmodel();
	cout << model << "内部" << endl;
	viewdir(modelpath);
	cout << "選択日付フォルダ->";
	cin >> date;
	datepath = modelpath + "/" + date;
	struct stat statBuf;
	if (stat(datepath.c_str(), &statBuf) != 0) {

	}
	else {
		cout << datepath << "は存在しません" << endl;
		exit(-1);
	}
}
void InputData::setdata() {
	setdate();
	cout << date << "内部" << endl;
	viewdir(datepath);
	cout << "選択フォルダ->";
	cin >> data;
	datapath = datepath + "/" + data;
	struct stat statBuf;
	if (stat(datapath.c_str(), &statBuf) != 0) {

	}
	else {
		cout << datapath << "は存在しません" << endl;
		exit(-1);
	}
}
void InputData::setgoal() {
	setdata();
	goalpath = datapath + "/" + goal;
	struct stat statBuf;
	if (stat(goalpath.c_str(), &statBuf) != 0) {

	}
	else {
		cout << goalpath << "は存在しません" << endl;
		exit(-1);
	}
}
void InputData::input_param() {
	string filename;
	filename = goalpath + "/" + "param.csv";
	ifstream inputfile(filename);
	string line;
	int count = 0;

	if (inputfile.is_open()) {
		while (getline(inputfile, line)) {
			stringstream ss(line);
			string value;
			vector<string> v;
			while (getline(ss, value, ',')) {
				v.push_back(value);
			}
			if (count == 0) {//1行目 xb,xt,yb,yt,,
				nparam.setXb(stod(v[0]));
				nparam.setXt(stod(v[1]));
				nparam.setYb(stod(v[2]));
				nparam.setYt(stod(v[3]));
			}
			if (count == 1) {//2行目 dx,dy,,,,
				nparam.setDx(stod(v[0]));
				nparam.setDy(stod(v[1]));
			}
			if (count == 2) {//3行目 xelem,yelem,,,,
				nparam.setXelem(stoi(v[0]));
				nparam.setYelem(stoi(v[1]));
			}
			if (count == 3) {//4行目 nelem,nnode,,,,
				nparam.setNelem(stoi(v[0]));
				nparam.setNnode(stoi(v[1]));
			}
			if (count == 4) {//5行目 flagL,flagR,flagU,flagD,flagC,
				BC.setBCflagL(stoi(v[0]));
				BC.setBCflagR(stoi(v[1]));
				BC.setBCflagU(stoi(v[2]));
				BC.setBCflagD(stoi(v[3]));
				BC.setBCflagC(stoi(v[4]));
			}
			count++;
		
		}
		inputfile.close();
	}
	else {
		cout << "ファイルが開けませんでした" << endl;
		exit(-1);
	}

	nparam.setLx(nparam.getXt() - nparam.getXb());
	nparam.setLy(nparam.getYt() - nparam.getYb());
	nparam.setXnode(nparam.getXelem() + 1);
	nparam.setYnode(nparam.getYelem() + 1);

}
void InputData::input_mesh() {
	string filename;
	filename = goalpath + "/" + "mesh.dat";
	ifstream inputfile(filename);
	string line;
	int count = 0;
	int nodecount = 0;

	if (inputfile.is_open()) {
		while (getline(inputfile, line)) {
			stringstream ss(line);
			string value;
			vector<string> v;
			while (getline(ss, value, ',')) {
				v.push_back(value);
			}
			if (count == 0 && nodecount == 0) {//1行目 nnode,nelem,,,,
				nparam.setNnode(stoi(v[0]));
				nparam.setNelem(stoi(v[1]));
				nodecount = 1;
			}
			if (nodecount == 1) {//2行目から座標データの終端まで
				if (stoi(v[0]) < nparam.getNnode()) {
					x.push_back(stod(v[1]));
					y.push_back(stod(v[2]));
				}
				if (stoi(v[0]) == nparam.getNnode() - 1) {
					nodecount = 2;
				}
			}
			if (nodecount == 2) {
				if (stoi(v[0]) < nparam.getNelem()) {//nbool1データの終端まで
					i1.push_back(stoi(v[1]));
					i2.push_back(stoi(v[2]));
					i3.push_back(stoi(v[3]));
					i4.push_back(stoi(v[4]));
				}
				if (stoi(v[0]) == nparam.getNelem() - 1) {
					nodecount = 3;
				}
			}
			if (nodecount == 3) {
				if (stoi(v[0]) < nparam.getNelem()) {//nbool3データの終端まで
					e1.push_back(stoi(v[1]));
					e2.push_back(stoi(v[2]));
					e3.push_back(stoi(v[3]));
					e4.push_back(stoi(v[4]));
				}
				if (stoi(v[0]) == nparam.getNelem() - 1) {
					nodecount = 4;
				}
			}
		}
		inputfile.close();
	}
	else {
		cout << "ファイルが開けませんでした" << endl;
		exit(-1);
	}

	nbool1.resize(nparam.getNelem());
	nbool3.resize(nparam.getNelem());
	for (int ie = 0; ie < nparam.getNelem(); ie++) {
		nbool1[ie].resize(4);
		nbool3[ie].resize(4);
	}
	for (int ie = 0; ie < nparam.getNelem(); ie++) {
		nbool1[ie][0] = i1[ie];
		nbool1[ie][1] = i2[ie];
		nbool1[ie][2] = i3[ie];
		nbool1[ie][3] = i4[ie];

		nbool3[ie][0] = e1[ie];
		nbool3[ie][1] = e2[ie];
		nbool3[ie][2] = e3[ie];
		nbool3[ie][3] = e4[ie];
	}
}
void InputData::input_csv() {
	
	int count = 0;
	string filename1 = goalpath + "/" + "BC.csv";
	string filename2 = goalpath + "/" + "P.csv";
	string filename3 = goalpath + "/" + "Ux.csv";
	string filename4 = goalpath + "/" + "Uy.csv";
	vector<string> filenames = { filename1,filename2,filename3,filename4 };

	for (const auto& filename : filenames) {
		ifstream file(filename);
		string line;
		if (file.is_open()) {
			vector<string> v;
			while (getline(file, line)) {
				stringstream ss(line);
				string value;
				while (getline(ss, value, ',')) {
					v.push_back(value);
				}
				
			}
			if (count == 0) {//1ファイル目BC.csv
				for (int i = 0; i < v.size(); i++) {
					cond.push_back(stoi(v[i]));
				}
			}
			if (count == 1) {//2ファイル目P.csv
				for (int i = 0; i < v.size(); i++) {
					P.push_back(stod(v[i]));
				}
			}
			if (count == 2) {//3ファイル目Ux.csv
				for (int i = 0; i < v.size(); i++) {
					Ux.push_back(stod(v[i]));
				}
			}
			if (count == 3) {//4ファイル目Uy.csv
				for (int i = 0; i < v.size(); i++) {
					Uy.push_back(stod(v[i]));
				}
			}
			file.close();
			count++;
		}
		else {
			cout << filename << "を開けませんでした" << endl;
			exit(-1);
		}
	}

}


void viewdir(string path) {
	for (const fs::directory_entry& x : fs::directory_iterator(path)) {
		cout << x.path() << endl;
	}
}


