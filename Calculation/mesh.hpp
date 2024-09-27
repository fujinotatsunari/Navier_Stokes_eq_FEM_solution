#pragma once
#include "input.hpp"
#include "mesh.hpp"
#include "param.hpp"
#include <vector>

class Node2d {//計算格子上の節点
private:
	int no;//節点番号
	double x;//座標x成分
	double y;//座標y成分
public:
	Node2d();
	void setNo(int no_);
	int getNo();
	void setX(double x_);
	double getX();
	void setY(double y_);
	double getY();

};


class Element2d {//計算格子上の要素
private:
	int no;//要素番号
	double x;//要素重心座標x成分
	double y;//要素重心座標y成分
	double Se;//要素面積
public:
	Element2d();
	void setNo(int no_);
	int getNo();
	void setX(double x_);
	double getX();
	void setY(double y_);
	double getY();
	void setSe(double Se_);
	double getSe();
};


class Time {//時間軸
private:
	double ntime_;//nステップでの時刻
	double dt_;//時間刻み幅
	int nend_;//終了ステップ
	int nsample_;//サンプルステップ
	int n_;//現在のステップ数
	vector<double> t;
	TimeP& tparam;
public:
	Time(TimeP& TP);
	//ゲッタ
	double ntime(int n);
	double ntime();
	void setn(int n);
	void setup();
	double dt();
	int nend();
	int nsample();
	double const& operator[](int n)const;
	double& operator[](int n);

};


class Mesh2d {//計算格子
	//今のところ四角形領域用のメッシュだけ対応(メッシュの左下から右上に向けて順番どおりに節点番号をつけられるようなmesh)
	//円柱周りメッシュとか, 円形領域メッシュとかどうしよ 
private:

	NodeP& nparam_;
	Boundarycond& Bcond_;
	vector<Node2d> node_;
	vector<Element2d> elem_;
	vector<int> ncond_;//節点境界フラグ
	vector<int> scond_;//要素フラグ
	
	//境界条件フラグの設定
	//ncond:  0:内部,	1:剛体内部(流入流出なし),
			//2:流入壁面(dirichlet)　3:流出壁面(neumann)
			//4:移動壁面条件(壁面接線方向に流速固定値1),5:滑りなし条件(壁面において(u,v)=(0,0))
			//6:滑りあり条件(dvx/dy=0,vy=0)

	//要素フラグの設定
	//scond: 0:内部 1:障害物内部
	
	vector<double> X;
	vector<double> Y;
	vector<double> EX;
	vector<double> EY;

	vector<vector<int>> nbool1_;//nbool[要素番号][要素内節点番号]=全体節点番号
	vector<vector<int>> nbool3_;//nbool3[要素番号][ローカルな要素番号(下,右,上,左)]=全体要素番号 ::ある要素に隣接する要素の番号
	double xb_, xt_, yb_, yt_, dx_, dy_, Lx_, Ly_;
	int xnode_, ynode_, xelem_, yelem_, nnode_, nelem_;
public:

	Mesh2d(InputData& input);
	Mesh2d(NodeP& NP, Boundarycond& BC);
	Mesh2d& operator=(const Mesh2d& mesh);
	Mesh2d(const Mesh2d& mesh);

	//ゲッタ
	double xb();
	double xt();
	double yb();
	double yt();
	double dx();
	double dy();
	double Lx();
	double Ly();
	double x(int i);//節点番号に対応するx座標
	double y(int i);//節点番号に対応するy座標
	double eX(int ie);//要素番号に対応するx座標(要素重心)
	double eY(int ie);//要素番号に対応するy座標(要素重心)
	double Se(int ie);//要素番号に対応する要素面積
	int xnode();
	int ynode();
	int xelem();
	int yelem();
	int nnode();
	int nelem();

	int nbool1(int ie, int np);
	int i1(int ie);//nbool1[ie][0]に対応する節点(左下)
	int i2(int ie);//nbool1[ie][1]に対応する節点(右下)
	int i3(int ie);//nbool1[ie][2]に対応する節点(右上)
	int i4(int ie);//nbool1[ie][3]に対応する節点(左上)
	int nbool3(int ie, int np);
	int e1(int ie);//nbool3[ie][0]に対応する要素(下)
	int e2(int ie);//nbool3[ie][1]に対応する要素(右)
	int e3(int ie);//nbool3[ie][2]に対応する要素(上)
	int e4(int ie);//nbool3[ie][3]に対応する要素(左)
	int ncond(int i);
	int scond(int ie);

	double length(int i1, int i2);//2点の線分の長さ
	double area(int e1, int e2);//2要素とそれを挟む辺の両端がなす面積
	double area(int i1, int i2, int i3, int i4);//4点からなる四角形の面積
	void geninputmesh();//inputしたmeshの可視化
	
};

