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


OutputData::OutputData(Mesh2d& Mesh, Time& T, Boundarycond& Bc, Velocity2d& v, Pressure& p, NDNSparam& Param)
	:mesh(Mesh),t(T), BC(Bc),V(v),P(p), Filestage(0),param(Param)
{

}
void OutputData::set_Filestage(int filestage) {
	Filestage = filestage;
}
void OutputData::set_model(string Model) {
	model = Model;
}
void OutputData::set_scheme(string Scheme) {
	scheme = Scheme;
}
void OutputData::set_path(string Path) {
	meshpath = Path;
}
string OutputData::get_scheme() {
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
	string dirname3 = "Calculation";//calculation or mesh
	string dirname4 = model;
	string dirname5 = scheme;
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
	str = make_directories(str, str1);//C:/..../date

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

	string dirU = "U";
	string dirV = "V";
	string dirP = "P";
	string dirMagnitude = "Magnitude";

	string filename1;//U_n.csv
	string filename2;//V_n.csv
	string filename3;//P.csv
	string filename4;//magnitude_n.csv

	string str1;
	string str2;
	string str3;
	string pathU, pathV, pathP, pathMagnitude;
	
	str1 = directory_setup();
	//str1 = directories_setup(n, scheme);//C:///.../data_(n)

	pathU = make_directories(str1, dirU);
	pathV = make_directories(str1, dirV);
	pathP = make_directories(str1, dirP);
	pathMagnitude = make_directories(str1, dirMagnitude);


	filename1 = pathU+ "/" + "U_" + to_string(N) + ".csv";
	filename2 = pathV + "/" + "V_" + to_string(N) + ".csv";
	filename3 = pathP + "/" + "P_" + to_string(N) + ".csv";
	filename4 = pathMagnitude + "/" + "magnitude_" + to_string(N) + ".csv";
	ofstream outputfile1(filename1);
	ofstream outputfile2(filename2);
	ofstream outputfile3(filename3);
	ofstream outputfile4(filename4);


	for (int j = 0; j < mesh.ynode(); j++) {
		for (int i = 0; i < mesh.xnode(); i++) {
			int np = i + mesh.xnode() * j;
			outputfile1 << V[np][0] << ",";
			outputfile2 << V[np][1] << ",";
			outputfile4 << V[np].norm() << ",";
		}
		outputfile1 << "\n";
		outputfile2 << "\n";
		outputfile4 << "\n";
	}

	for (int j = 0; j < mesh.yelem(); j++) {
		for (int i = 0; i < mesh.xelem(); i++) {
			int ie = i + mesh.xelem() * j;
			outputfile3 << P[ie].v() << ",";
		}
		outputfile3 << "\n";
	}

	outputfile1.close();
	outputfile2.close();
	outputfile3.close();
	outputfile4.close();
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
void OutputData::output_ghia(int N) {
	string dirG = "Ghia";
	string filename1;//ghiax.csv x軸上でのghiaの結果ghiav
	string filename2;//ghiay.csv y軸上でのghiaの結果ghiau
	string str1;
	string pathG;
	get_ghia();

	str1 = directory_setup();
	pathG = make_directories(str1, dirG);
	filename1 = pathG + "/" + "Ghiax_" + to_string(N) + ".csv";
	filename2 = pathG + "/" + "Ghiay_" + to_string(N) + ".csv";

	ofstream outputfile1(filename1);
	ofstream outputfile2(filename2);
	for (int i = 0; i < ghiav.size(); i++) {
		outputfile1 << ghiax[i] << "," << ghiav[i] << "," << "\n";
	}
	outputfile1.close();
	for (int j = 0; j < ghiau.size(); j++) {
		outputfile2 << ghiay[j] << "," << ghiau[j] << "," << "\n";
	}
	outputfile2.close();

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
	outputfile << "Re=" << param.get_Re() << "\n";
	outputfile << "mesh_file_path->" << meshpath << "\n";
	outputfile << "model->" << model << "\n";
	outputfile << "scheme->" << scheme << "\n";


	outputfile.close();
}
void OutputData::get_ghia() {

	for (int j = 0; j < mesh.ynode(); j++) {
		for (int i = 0; i < mesh.xnode(); i++) {
			int np = i + mesh.xnode() * j;
			if (mesh.x(np) <= mesh.Lx() / 2) {
				if (mesh.x(np + 1) > mesh.Lx() / 2) {
					ghiau.push_back(V[np][0]);//y軸上でのx方向流速
					ghiay.push_back(mesh.y(np));//対応するy軸座標
				}
			}

			if (mesh.y(np) <= mesh.Ly() / 2) {
				if (mesh.y(np + mesh.xnode()) > mesh.Ly() / 2) {
					ghiav.push_back(V[np][1]);//x軸上でのy方向流速
					ghiax.push_back(mesh.x(np));//対応するx軸座標
				}
			}
		}
	}
}
string OutputData::get_model() {
	return model;
}
void OutputData::update(Velocity2d& v, Pressure& p) {
	for (int i = 0; i < mesh.nnode(); i++) {
		V[i][0] = v[i][0];
		V[i][1] = v[i][1];
		
	}
	for (int ie = 0; ie < mesh.nelem(); ie++) {
		P[ie][0] = p[ie][0];
	}
}