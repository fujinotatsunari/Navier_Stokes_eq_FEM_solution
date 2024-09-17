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
	int scheme;//計算スキームのフラグ
	int Filestage;//保存場所を作成したか否かのフラグ
	//0:未作成,1:作成済み
	//vector<double> x;
	//vector<double> y;
	vector<vector<double>> nodecopy;//Ux,Uy,BCのcsv用の配列
	vector<vector<double>> elemcopy;//Pのcsv用の配列
	Mesh2d& mesh;
	Velocity2d& V;
	Pressure& P;
	Boundarycond& BC;
	string dir;//ディレクトリdata_nの位置
public:
	OutputData(Mesh2d& Mesh, Boundarycond& Bc, Velocity2d& v, Pressure& p);
	void set_scheme(int Scheme);
	void set_Filestage(int filestage);
	int get_scheme();
	int get_Filestage();
	string get_dir();

	virtual string directory_setup();//下におけるdata_nを返す
	/*
	* C:-//--|result|----|data0|-----|cond|
	*           |            |
	*           |             -------|resultcsv|
	*           |-------|data1|
	*           |-------|data2|
	*/
	virtual void output_time_csv(int N);//Nステップでのcsvの出力
	virtual void output_csv();//csvの出力
	virtual void output_dat();//datの出力
	virtual void output_condition();//パラメータファイルの出力
	//出力の仕様変更
	//原則１ファイルが1個の多次元(あるいは1次元)配列に対応するようにする
	// ncond.csv scond.csv U.csv V.csv P.csv 　<=　1次元配列
	// node.csv:節点の(x,y)座標配列　 snode.csv：要素重心の(x,y)座標配列  <=2次元配列
	// nbool1.csv :nbool1のデータ nbool3.csv :nbool3のデータ <=nelem*4 配列
	// param.csv <=空間パラメータと境界条件の保存 あとモデル
	//void data_update();//データの更新
	
};

class Outputcavitymesh :public OutputData {
private:
	string dir1 = "cavity";
public:
	Outputcavitymesh(Mesh2d& Mesh, Boundarycond& Bc, Velocity2d& v, Pressure& p);
	string directory_setup() override;
	void output_csv() override;
	//void output_dat() override;
	void output_condition() override;
};
class Outputbackstepmesh :public OutputData {
private:
	string dir1 = "backstep";
public:
	Outputbackstepmesh(Mesh2d& Mesh, Boundarycond& Bc, Velocity2d& v, Pressure& p);
	string directory_setup() override;
	void output_csv() override;
	//void output_dat() override;
	void output_condition() override;
};


string make_directories(string str1, string str2);
/*ディレクトリを作成したい場所の文字列を作る
親ディレクトリの名前str1と子ディレクトリの名前str2を
結合しディレクトリstr1/str2を作りリターンする*/


//string directories_setup(int n, int scheme);
/*結果を保存するためのディレクトリを作成する*/

