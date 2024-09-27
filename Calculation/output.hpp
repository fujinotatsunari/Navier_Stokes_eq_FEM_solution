#pragma once
#include"output.hpp"
#include"mesh.hpp"
#include"value.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include<vector>
using namespace std;

class OutputData {
protected:
	int Filestage;//保存場所を作成したか否かのフラグ
	//0:未作成,1:作成済み
	vector<double> ghiax;//ghiavに対応するx軸座標
	vector<double> ghiay;//ghiauに対応するy軸座標
	vector<double> ghiau;//y軸上でのx方向流速
	vector<double> ghiav;//x軸上でのy方向流速
	Mesh2d& mesh;
	Time& t;
	Velocity2d& V;
	Pressure& P;
	Boundarycond& BC;
	NDNSparam& param;
	string meshpath;//計算したmeshのパス
	string model;//計算したモデルinputから持ってくる
	string scheme;//
	string dir;//ディレクトリdata_nの位置
	
public:
	OutputData(Mesh2d& Mesh, Time& T, Boundarycond& Bc, Velocity2d& v, Pressure& p, NDNSparam& Param);
	void set_Filestage(int filestage);
	void set_model(string Model);
	void set_scheme(string Scheme);
	void set_path(string Path);
	string get_scheme();
	int get_Filestage();
	string get_dir();
	string directory_setup();//下におけるdata_nを返す
	void output_setup(string Scheme, string Model, string Path);

	/*
	* C:/Result/2d_Navier_Stokes_eq/calculation/model/scheme/date/data_n
	* |data_n|-----|condtion|
	*          |---|U|//流速x方向フォルダ
	*		　 |---|V|//流速y方向フォルダ
	*		　 |---|magnitude|//流速の大きさ
	*          |---|P|//圧力
	*          |---|other|//解析対象とかで色々変わるであろうデータ
	*          |---.....//
	*/

	void output_time_csv(int N);//Nステップでのcsvの出力(U,V,P,|V|)
	void output_ghia(int N);//中心軸(x,y)上での流速
	void output_csv();//csvの出力
	//void output_dat();//datの出力

	void output_condition();//パラメータファイルの出力
	void update(Velocity2d& v, Pressure& p);
	void get_ghia();//ghiaの結果の生成(中心軸上の流速の取得)
	string get_model();
};



string make_directories(string str1, string str2);
/*ディレクトリを作成したい場所の文字列を作る
親ディレクトリの名前str1と子ディレクトリの名前str2を
結合しディレクトリstr1/str2を作りリターンする*/
