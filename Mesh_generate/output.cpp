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
	nodecopy.resize(mesh.xnode());
	elemcopy.resize(mesh.xelem());
	for (int i = 0; i < nodecopy.size(); i++) {
		nodecopy[i].resize(mesh.ynode());
	}
	for (int i = 0; i < elemcopy.size(); i++) {
		elemcopy[i].resize(mesh.yelem());
	}
	//dir = directory_setup();
	/*
	x.resize(mesh.xnode());
	y.resize(mesh.ynode());
	copy.resize(mesh.xnode());
	for (int i = 0; i < copy.size(); i++) {
		copy[i].resize(mesh.ynode());
	}
	for (int j = 0; j < mesh.ynode(); j++) {
		for (int i = 0; i < mesh.xnode(); i++) {
			int np = i + mesh.xnode() * j;
			//copy[i][j] = phi[np];

			if (i == 0) {
				y[j] = mesh.y(np);
			}
			if (j == 0) {
				x[i] = mesh.x(np);
			}
		}
	}
	*/
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
	string dirname3 = "Data";//calculation or mesh
	string dirname4 = "Model";
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
void OutputData::output_time_csv(int N) {

	string str;
	string str1;
	string str2;
	string str3;

	//str1 = directories_setup(n, scheme);//C:///.../data_(i)
	str1 = directory_setup();
	str2 = "result_csv";
	str3 = make_directories(str1, str2);//C:/..../data_(i)/result_csv
	str = str3 + "/" + "step_" + to_string(N) + "_" + "phi.csv";

	ofstream outputfile(str);
	/*ファイルの中身
	outputfile << ",";
	for (int i = 0; i < x.size(); i++) {
		outputfile << x[i] << ",";
	}
	outputfile << "\n";
	for (int j = 0; j < y.size(); j++) {
		outputfile << y[j] << ",";
		for (int i = 0; i < x.size(); i++) {
			outputfile << copy[i][j] << ",";
		}
		outputfile << "\n";
	}
	*/

	outputfile.close();
}

void OutputData::output_csv() {
	string str;
	string str1;
	string str2;
	string str3;
	str1 = directory_setup();//C:///.../data_(i)
	str2 = "data_box";
	str3 = make_directories(str1, str2);//C:/..../data_(i)/condition
	str = str3 + "/" + "data.csv";
	ofstream outputfile(str);
	outputfile.close();
}
void OutputData::output_dat() {
	string str;
	string str1;
	string str2;
	string str3;
	str1 = directory_setup();//C:///.../data_(i)
	str2 = "data_box";
	str3 = make_directories(str1, str2);//C:/..../data_(i)/condition
	str = str3 + "/" + "data.dat";
	ofstream outputfile(str);
	outputfile.close();
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

	/*
	outputfile << "#calculation_condition" << "\n";
	outputfile << "##mesh_parameter" << "\n";
	outputfile << "xb=" << mesh.xb() << " " << "xt=" << mesh.xt() << "\n";
	outputfile << "yb=" << mesh.yb() << " " << "yt=" << mesh.yt() << "\n";
	outputfile << "dx=" << mesh.dx() << " " << "dy=" << mesh.dy() << "\n";
	outputfile << "xelem=" << mesh.xelem() << " " << "yelem=" << mesh.yelem() << "\n";
	outputfile << "nelem=" << mesh.nelem() << " " << "nnode=" << mesh.nnode() << "\n";
	outputfile << "\n";
	outputfile << "##time_parameter" << "\n";
	outputfile << "dt=" << t.dt() << "\n";
	outputfile << "nend=" << t.nend() << "\n";
	outputfile << "nsample=" << t.nsample() << "\n";
	outputfile << "\n";
	outputfile << "##equation_parameter" << "\n";
	outputfile << "x方向定常流速 cx=" << adp.get_cx() << " " << "y方向定常流速 cy=" << adp.get_cy() << "\n";
	outputfile << "拡散係数 alpha=" << adp.get_alpha() << "\n";
	outputfile << "x方向courant数" << adp.get_couranx() << " " << "y方向courant数" << adp.get_courany() << "\n";
	outputfile << "拡散数" << adp.get_diffusion() << " " << "Peclet数" << adp.get_Pe() << "\n";
	*/

	outputfile.close();
}

Outputcavitymesh::Outputcavitymesh(Mesh2d& Mesh, Boundarycond& Bc, Velocity2d& v, Pressure& p)
	:OutputData(Mesh, Bc, v, p)
{

}
string Outputcavitymesh::directory_setup() {
	string dirname0 = "C:";
	string dirname1 = "Result";
	string dirname2 = "2d_Navier_Stokes_eq";
	string dirname3 = "Mesh_data_box";//calculation or mesh
	string dirname4 = dir1;
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

void Outputcavitymesh::output_csv() {
	string filename1;
	string filename2;
	string filename3;
	string filename4;
	string filename5;
	string str1;
	string str2;
	string str3;
	str1 = directory_setup();//C:///.../data_(i)
	str2 = "mesh_data_box";
	str3 = make_directories(str1, str2);//C:/..../data_(i)/mesh_data_box
	filename1 = str3 + "/" + "U.csv";//x方向流速
	filename2 = str3 + "/" + "V.csv";//y方向流速
	filename3 = str3 + "/" + "P.csv";//圧力
	filename4 = str3 + "/" + "BC.csv";//境界条件
	filename5 = str3 + "/" + "param.csv";//パラメーター
	ofstream outputfile1(filename1);
	ofstream outputfile2(filename2);
	ofstream outputfile3(filename3);
	ofstream outputfile4(filename4);
	ofstream outputfile5(filename5);

	for (int j = 0; j < mesh.ynode(); j++) {
		for (int i = 0; i < mesh.xnode(); i++) {
			int np = i + mesh.xnode() * j;
			nodecopy[i][j] = V[np][0];//x方向流速
		}
	}
	for (int j = 0; j < mesh.ynode(); j++) {
		for (int i = 0; i < mesh.xnode(); i++) {
			outputfile1 << nodecopy[i][j] << ",";
		}
		outputfile1 << "\n";
	}
	outputfile1.close();

	for (int j = 0; j < mesh.ynode(); j++) {
		for (int i = 0; i < mesh.xnode(); i++) {
			int np = i + mesh.xnode() * j;
			nodecopy[i][j] = V[np][1];//y方向流速
		}
	}
	for (int j = 0; j < mesh.ynode(); j++) {
		for (int i = 0; i < mesh.xnode(); i++) {
			outputfile2 << nodecopy[i][j] << ",";
		}
		outputfile2 << "\n";
	}
	outputfile2.close();

	for (int j = 0; j < mesh.yelem(); j++) {
		for (int i = 0; i < mesh.xelem(); i++) {
			int ie = i + mesh.xelem() * j;
			elemcopy[i][j] = P[ie][0];//圧力
		}
	}
	for (int j = 0; j < mesh.yelem(); j++) {
		for (int i = 0; i < mesh.xelem(); i++) {
			outputfile3 << elemcopy[i][j] << ",";
		}
		outputfile3 << "\n";
	}
	outputfile3.close();

	for (int j = 0; j < mesh.ynode(); j++) {
		for (int i = 0; i < mesh.xnode(); i++) {
			int np = i + mesh.xnode() * j;
			nodecopy[i][j] = mesh.ncond(np);//境界条件
		}
	}
	for (int j = 0; j < mesh.ynode(); j++) {
		for (int i = 0; i < mesh.xnode(); i++) {
			outputfile4 << nodecopy[i][j] << ",";
		}
		outputfile4 << "\n";
	}
	outputfile4.close();
	/*
	outputfile << "xb=" << mesh.xb() << " " << "xt=" << mesh.xt() << "\n";
	outputfile << "yb=" << mesh.yb() << " " << "yt=" << mesh.yt() << "\n";
	outputfile << "dx=" << mesh.dx() << " " << "dy=" << mesh.dy() << "\n";
	outputfile << "xelem=" << mesh.xelem() << " " << "yelem=" << mesh.yelem() << "\n";
	outputfile << "nelem=" << mesh.nelem() << " " << "nnode=" << mesh.nnode() << "\n";
	outputfile << "flagU=" << BC.getBCflagU() << " " << "flagD=" << BC.getBCflagD() << " " << "flagL=" << BC.getBCflagL() << " " << "flagR=" << BC.getBCflagR() << "\n";
	outputfile << "\n";
	*/
	outputfile5 << mesh.xb() << "," << mesh.xt() << "," << mesh.yb() << "," << mesh.yt() << "," << "\n";
	outputfile5 << mesh.dx() << "," << mesh.dy() << "," << "," << "," << "\n";
	outputfile5 << mesh.xelem() << "," << mesh.yelem() << "," << "," << ", " << "\n";
	outputfile5 << mesh.nelem() << "," << mesh.nnode() << "," << "," << ", " << "\n";
	outputfile5 << BC.getBCflagL() << "," << BC.getBCflagR() << "," << BC.getBCflagU() << "," << BC.getBCflagD() << "," << -1 << "," << "\n";
	outputfile5.close();
}
void Outputcavitymesh::output_dat() {
	string filename1;
	string filename2;
	string filename3;

	string str1;
	string str2;
	string str3;
	str1 = directory_setup();//C:///.../data_(i)
	str2 = "mesh_data_box";
	str3 = make_directories(str1, str2);//C:/..../data_(i)/mesh_data_box
	filename1 = str3 + "/" + "mesh.dat";//meshの構成
	filename2 = str3 + "/" + "init.dat";//初期条件
	filename3 = str3 + "/" + "cond.dat";//境界条件


	ofstream outputfile1(filename1);
	ofstream outputfile2(filename2);
	ofstream outputfile3(filename3);


	outputfile1 << mesh.nnode() << "," << mesh.nelem() << "," << "," << "," << "," << "\n";
	for (int j = 0; j < mesh.ynode(); j++) {
		for (int i = 0; i < mesh.xnode(); i++) {
			int np = i + mesh.xnode() * j;
			outputfile1 << np << "," << mesh.x(np) << "," << mesh.y(np) << "," << "," << "," << "\n";
		}
	}
	for (int j = 0; j < mesh.yelem(); j++) {
		for (int i = 0; i < mesh.xelem(); i++) {
			int ie = i + mesh.xelem() * j;
			outputfile1 << ie << "," << mesh.i1(ie) << "," << mesh.i2(ie) << "," << mesh.i3(ie) << "," << mesh.i4(ie) << "," << "\n";
		}
	}
	for (int j = 0; j < mesh.yelem(); j++) {
		for (int i = 0; i < mesh.xelem(); i++) {
			int ie = i + mesh.xelem() * j;
			outputfile1 << ie << "," << mesh.e1(ie) << "," << mesh.e2(ie) << "," << mesh.e3(ie) << "," << mesh.e4(ie) << "," << "\n";
		}
	}
	outputfile1.close();

	outputfile2 << mesh.nnode() << "," << mesh.nelem() << "," << "," << "\n";
	for (int j = 0; j < mesh.ynode(); j++) {
		for (int i = 0; i < mesh.xnode(); i++) {
			int np = i + mesh.xnode() * j;
			outputfile2 << np << "," << V[np][0] << "," << V[np][1] << "," << "\n";
		}
	}
	for (int j = 0; j < mesh.yelem(); j++) {
		for (int i = 0; i < mesh.xelem(); i++) {
			int ie = i + mesh.xelem() * j;
			outputfile2 << ie << "," << P[ie][0] << "," << "," << "\n";
		}
	}
	outputfile2.close();

	outputfile3 << mesh.nnode() << "," << mesh.nelem() << "," << "," << "\n";
	for (int j = 0; j < mesh.ynode(); j++) {
		for (int i = 0; i < mesh.xnode(); i++) {
			int np = i + mesh.xnode() * j;
			outputfile3 << np << "," << mesh.ncond(np) << "," << "," << "\n";
		}
	}
	outputfile3.close();

}
void Outputcavitymesh::output_condition() {
	string filename;
	string str1;
	string str2;
	string str3;
	str1 = directory_setup();//C:///.../data_(i)
	str2 = "condition";
	str3 = make_directories(str1, str2);//C:/..../data_(i)/condition
	filename = str3 + "/" + "condition.txt";

	ofstream outputfile(filename);


	outputfile << "xb=" << mesh.xb() << " " << "xt=" << mesh.xt() << "\n";
	outputfile << "yb=" << mesh.yb() << " " << "yt=" << mesh.yt() << "\n";
	outputfile << "dx=" << mesh.dx() << " " << "dy=" << mesh.dy() << "\n";
	outputfile << "xelem=" << mesh.xelem() << " " << "yelem=" << mesh.yelem() << "\n";
	outputfile << "nelem=" << mesh.nelem() << " " << "nnode=" << mesh.nnode() << "\n";
	outputfile << "flagU=" << BC.getBCflagU() << " " << "flagD=" << BC.getBCflagD() << " " << "flagL=" << BC.getBCflagL() << " " << "flagR=" << BC.getBCflagR() << "\n";
	outputfile << "\n";

	outputfile.close();
}
