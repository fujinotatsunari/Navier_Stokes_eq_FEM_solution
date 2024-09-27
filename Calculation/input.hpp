#pragma once
#include "input.hpp"
#include "param.hpp"
#include <string>
#include <vector>
using namespace std;

class InputData {
private:
	string start = "C:/Result/2d_Navier_Stokes_eq/Mesh_data_box";
	string model;//解析モデルフォルダのなまえ
	string modelpath;//modelのパス
	string date;//日付フォルダのなまえ
	string datepath;//dateのパス
	string data;//data_nフォルダのなまえ
	string datapath;//dataのパス
	string goal = "mesh_data_box";//目的フォルダ1
	string goal2 = "data_box";//目的フォルダ2
	string goalpath;//目的フォルダ1のパス
	string goalpath2;//目的フォルダ2のパス
	vector<double> Ux;
	vector<double> Uy;
	vector<double> P;
	vector<int> cond;
	vector<int> scond;
	vector<double> x;
	vector<double> y;
	vector<double> ex;
	vector<double> ey;
	vector<vector<int>> nbool1;
	vector<vector<int>> nbool3;
	vector<int> i1;
	vector<int> i2;
	vector<int> i3;
	vector<int> i4;
	vector<int> e1;
	vector<int> e2;
	vector<int> e3;
	vector<int> e4;
	NodeP nparam;
	Boundarycond BC;

	/*
	* C:/Result/2d_Navier_Stokes/Mesh_data_box/model/date/data_n/mesh_data_box
	*/
public:
	InputData();

	string getmodel();//modelディレクトリパスの取得
	string getdate();//dateディレクトリパスの取得
	string getdata();//data_nディレクトリパスの取得
	string getgoal();
	
	void setmodel();
	void setdate();
	void setdata();
	void setgoal();
	void input_param();
	void input_csv();
	void input_dat();


	double recomend_dt();//時間刻み幅の推奨
	//初期条件の流速の大きさの最大値からクーラン数が
	// C=0.01より小さくなるような時間刻み幅を提示する
	//void intput_data();
	//nboolのゲッタ
	vector<vector<int>> getnbool1();
	vector<vector<int>> getnbool3();
	//vectorのゲッタ
	vector<double> getUx();
	vector<double> getUy();
	vector<double> getx();
	vector<double> gety();
	vector<double> getex();
	vector<double> getey();
	vector<double> getP();
	vector<int> getcond();
	vector<int> getscond();
	//クラスのゲッタ
	NodeP& get_NodeParam();
	Boundarycond& get_BC();
	string get_path();
	string get_modelname();

};
void viewdir(string path);//path直下のフォルダを表示