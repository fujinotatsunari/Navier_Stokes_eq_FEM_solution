#include"output.hpp"
#include <string>
#include <time.h>
#include <direct.h>
#include <sys/stat.h>
#include <iostream>
#include <fstream>

using namespace std;

string make_directories(string str1, string str2) {

	string str;
	str = str1 + "/" + str2;
	struct stat statBuf;
	if (stat(str.c_str(), &statBuf) != 0) {
		if (_mkdir(str.c_str()) == 0) {
			// 成功
			cout << str << "が作成されました" << endl;
			return str;
		}
		else {
			// 失敗
			cout << str << "の作成に失敗しました" << endl;
			exit(-1);
		}
	}
	else {
		return str;
	}

}



OutputData::OutputData(Mesh2d& Mesh, Boundarycond& Bc, Velocity2d& v, Pressure& p)
	:mesh(Mesh), BC(Bc),V(v),P(p), scheme(-1), Filestage(0)
{
}
void OutputData::set_scheme(int Scheme) {
	scheme = Scheme;
}
void OutputData::set_Filestage(int filestage) {
	Filestage = filestage;
}
int OutputData::get_scheme() {
	return scheme;
}
int OutputData::get_Filestage() {
	return Filestage;
}
string OutputData::get_dir() {
	return dir;
}

string OutputData::directory_setup() {
	string dirname0 = "C:";
	string dirname1 = "Result";
	string dirname2 = "2d_Navier_Stokes_eq";
	string dirname3 = "Mesh_data_box";//calculation or mesh
	string dirname4 = model;
	string str;
	string str1;
	string year;
	string month;
	string day;

	time_t timer;
	struct tm local_time;
	timer = time(NULL);
	localtime_s(&local_time, &timer);
	struct stat statBuf;

	str = make_directories(make_directories(make_directories(make_directories(dirname0, dirname1), dirname2), dirname3), dirname4);

	year = to_string(local_time.tm_year + 1900);
	month = to_string(local_time.tm_mon + 1);
	day = to_string(local_time.tm_mday);

	str1 = year + month + day;
	str = make_directories(str, str1);//C:/..../day

	int check = 0;
	if (Filestage == 0) {//計算開始後初めて保存場所を作成
		for (int i = 0; check == 0; i++) {
			str1 = str + "/" + "data_" + to_string(i);//C:/../data_i
			//cout << str1 << endl;
			if (stat(str1.c_str(), &statBuf) != 0) {
				//data_iがそんざいしないとき
				str1 = "data_" + to_string(i);
				str = make_directories(str, str1);
				//cout << str << "を作成" <<  endl;
				check = 1;
				dir = str;
				Filestage = 1;
				return str;
			}
			else {
				//data_iがそんざいするときiをインクリメント
			}
		}
	}
	else if (Filestage == 1) {//保存場所が作成済み

		return dir;
	}
	else {
		exit(-1);
	}

}

void OutputData::output_csv() {
	string filename1;
	string filename2;
	string filename3;
	string filename4;
	string filename5;
	string filename6;
	string filename7;
	string filename8;
	string filename9;
	string filename10;
	string str1;
	string str2;
	string str3;
	str1 = directory_setup();//C:///.../data_(i)
	str2 = "mesh_data_box";
	str3 = make_directories(str1, str2);//C:/..../data_(i)/mesh_data_box
	filename1 = str3 + "/" + "U.csv";//x方向流速
	filename2 = str3 + "/" + "V.csv";//y方向流速
	filename3 = str3 + "/" + "P.csv";//圧力
	filename4 = str3 + "/" + "ncond.csv";//節点境界条件データ
	filename5 = str3 + "/" + "scond.csv";//要素条件データ
	filename6 = str3 + "/" + "node.csv";//節点座標配列
	filename7 = str3 + "/" + "snode.csv";//要素重心座標配列
	filename8 = str3 + "/" + "nbool1.csv";//nbool1データ
	filename9 = str3 + "/" + "nbool3.csv";//nbool3データ
	filename10 = str3 + "/" + "param.csv";//パラメーター

	ofstream outputfile1(filename1);
	ofstream outputfile2(filename2);
	ofstream outputfile3(filename3);
	ofstream outputfile4(filename4);
	ofstream outputfile5(filename5);
	ofstream outputfile6(filename6);
	ofstream outputfile7(filename7);
	ofstream outputfile8(filename8);
	ofstream outputfile9(filename9);
	ofstream outputfile10(filename10);

	
	for (int j = 0; j < mesh.ynode(); j++) {
		for (int i = 0; i < mesh.xnode(); i++) {
			int np = i + mesh.ynode() * j;
			outputfile1 << V[np][0] << ",";
		}
		outputfile1 << "\n";
	}
	outputfile1.close();

	for (int j = 0; j < mesh.ynode(); j++) {
		for (int i = 0; i < mesh.xnode(); i++) {
			int np = i + mesh.ynode() * j;
			outputfile2 << V[np][1] << ",";
		}
		outputfile2 << "\n";
	}
	outputfile2.close();


	for (int j = 0; j < mesh.yelem(); j++) {
		for (int i = 0; i < mesh.xelem(); i++) {
			int ie = i + mesh.yelem() * j;
			outputfile3 << P[ie].v() << ",";
		}
		outputfile3 << "\n";
	}
	outputfile3.close();



	for (int j = 0; j < mesh.ynode(); j++) {
		for (int i = 0; i < mesh.xnode(); i++) {
			int np = i + mesh.ynode() * j;
			outputfile4 << mesh.ncond(np) << ",";
		}
		outputfile4 << "\n";
	}
	outputfile4.close();

	for (int j = 0; j < mesh.yelem(); j++) {
		for (int i = 0; i < mesh.xelem(); i++) {
			int ie = i + mesh.yelem() * j;
			outputfile5 << mesh.scond(ie) << ",";
		}
		outputfile5 << "\n";
	}
	outputfile5.close();

	for (int i = 0; i < mesh.nnode(); i++) {
		outputfile6 << mesh.x(i) << "," << mesh.y(i) << "," << "\n";//節点座標
	}
	outputfile6.close();

	for (int ie = 0; ie < mesh.nelem(); ie++) {
		outputfile7 << mesh.eX(ie) << "," << mesh.eY(ie) << "," << "\n";//要素重心座標
	}
	outputfile7.close();

	for (int ie = 0; ie < mesh.nelem(); ie++) {
		outputfile8 << mesh.i1(ie) << "," << mesh.i2(ie) << "," << mesh.i3(ie) << "," << mesh.i4(ie) << "," << "\n";//nbool1
	}
	outputfile8.close();

	for (int ie = 0; ie < mesh.nelem(); ie++) {
		outputfile9 << mesh.e1(ie) << "," << mesh.e2(ie) << "," << mesh.e3(ie) << "," << mesh.e4(ie) << "," << "\n";//nbool3
	}
	outputfile9.close();

	outputfile10 << mesh.xb() << "," << mesh.xt() << "," << mesh.yb() << "," << mesh.yt() << "," << "," << "\n";
	outputfile10 << mesh.dx() << "," << mesh.dy() << "," << "," << "," << "," << "\n";
	outputfile10 << mesh.xelem() << "," << mesh.yelem() << "," << mesh.xnode() << "," << mesh.ynode() << "," << "," << "\n";
	outputfile10 << mesh.nelem() << "," << mesh.nnode() << "," << "," << "," << "," << "\n";
	outputfile10 << BC.getBCflagL() << "," << BC.getBCflagR() << "," << BC.getBCflagU() << "," << BC.getBCflagD() << "," << BC.getBCflagC() << "," << "\n";
	outputfile10.close();
}
void OutputData::output_dat() {
	string filename1;
	string filename2;
	string filename3;
	string filename4;
	string filename5;
	string filename6;
	string filename7;
	string filename8;
	string filename9;

	string str1;
	string str2;
	string str3;
	str1 = directory_setup();//C:///.../data_(i)
	str2 = "data_box";
	str3 = make_directories(str1, str2);//C:/..../data_(i)/mesh_data_box
	filename1 = str3 + "/" + "U.dat";//x方向流速
	filename2 = str3 + "/" + "V.dat";//y方向流速
	filename3 = str3 + "/" + "P.dat";//圧力
	filename4 = str3 + "/" + "ncond.dat";//節点境界条件データ
	filename5 = str3 + "/" + "scond.dat";//要素条件データ
	filename6 = str3 + "/" + "node.dat";//節点座標配列
	filename7 = str3 + "/" + "snode.dat";//要素重心座標配列
	filename8 = str3 + "/" + "nbool1.dat";//nbool1データ
	filename9 = str3 + "/" + "nbool3.dat";//nbool3データ
	//filename10 = str3 + "/" + "param.dat";//パラメーター

	ofstream outputfile1(filename1);
	ofstream outputfile2(filename2);
	ofstream outputfile3(filename3);
	ofstream outputfile4(filename4);
	ofstream outputfile5(filename5);
	ofstream outputfile6(filename6);
	ofstream outputfile7(filename7);
	ofstream outputfile8(filename8);
	ofstream outputfile9(filename9);
	//ofstream outputfile10(filename10);

	for (int i = 0; i < mesh.nnode(); i++) {
		outputfile1 << V[i][0] << "," << "\n";//x方向流速
	}
	outputfile1.close();

	for (int i = 0; i < mesh.nnode(); i++) {
		outputfile2 << V[i][1] << "," << "\n";//y方向流速
	}
	outputfile2.close();

	for (int ie = 0; ie < mesh.nelem(); ie++) {
		outputfile3 << P[ie][0] << "," << "\n";//圧力
	}
	outputfile3.close();

	for (int i = 0; i < mesh.nnode(); i++) {
		outputfile4 << mesh.ncond(i) << "," << "\n";//節点境界条件
	}
	outputfile4.close();

	for (int ie = 0; ie < mesh.nelem(); ie++) {
		outputfile5 << mesh.scond(ie) << "," << "\n";//要素条件
	}
	outputfile5.close();

	for (int i = 0; i < mesh.nnode(); i++) {
		outputfile6 << mesh.x(i) << "," << mesh.y(i) << "," << "\n";//節点座標
	}
	outputfile6.close();

	for (int ie = 0; ie < mesh.nelem(); ie++) {
		outputfile7 << mesh.eX(ie) << "," << mesh.eY(ie) << "," << "\n";//要素重心座標
	}
	outputfile7.close();

	for (int ie = 0; ie < mesh.nelem(); ie++) {
		outputfile8 << mesh.i1(ie) << "," << mesh.i2(ie) << "," << mesh.i3(ie) << "," << mesh.i4(ie) << "," << "\n";//nbool1
	}
	outputfile8.close();

	for (int ie = 0; ie < mesh.nelem(); ie++) {
		outputfile9 << mesh.e1(ie) << "," << mesh.e2(ie) << "," << mesh.e3(ie) << "," << mesh.e4(ie) << "," << "\n";//nbool3
	}
	outputfile9.close();

	
}

void OutputData::output_condition() {
	string str;
	string str1;
	string str2;
	string str3;
	str1 = directory_setup();//C:///.../data_(i)
	str2 = "condition";
	str3 = make_directories(str1, str2);//C:/..../data_(i)/condition
	str = str3 + "/" + "condition.txt";

	ofstream outputfile(str);
	outputfile << "xb=" << mesh.xb() << " " << "xt=" << mesh.xt() << "\n";
	outputfile << "yb=" << mesh.yb() << " " << "yt=" << mesh.yt() << "\n";
	outputfile << "dx=" << mesh.dx() << " " << "dy=" << mesh.dy() << "\n";
	outputfile << "xelem=" << mesh.xelem() << " " << "yelem=" << mesh.yelem() << "\n";
	outputfile << "nelem=" << mesh.nelem() << " " << "nnode=" << mesh.nnode() << "\n";
	outputfile << "flagU=" << BC.getBCflagU() << " " << "flagD=" << BC.getBCflagD() << " " << "flagL=" << BC.getBCflagL() << " " << "flagR=" << BC.getBCflagR() << " " << "flagC=" << BC.getBCflagC() << "\n";
	outputfile << "\n";

	outputfile.close();
}
void OutputData::set_modelname(string name) {
	model = name;
}

