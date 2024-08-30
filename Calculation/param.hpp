#pragma once
#include"param.hpp"
#include<vector>
#include<iostream>
#include<cstdio>
#include<cmath>
using namespace std;

class NodeP {
private:
	double Lx, Ly;
	double xb, yb;
	double xt, yt;
	double dx, dy;

	int xnode, ynode;
	int xelem, yelem;
	int nnode, nelem;
public:
	NodeP();
	NodeP(const NodeP& Np);
	NodeP& operator=(const NodeP& Np);
	void setNodeParam();
	double getLx();
	double getLy();
	double getXb();
	double getXt();
	double getYb();
	double getYt();
	double getDx();
	double getDy();
	//double getDxDy();
	int getXnode();
	int getYnode();
	int getXelem();
	int getYelem();
	int getNnode();
	int getNelem();
};
/*
class TimeP {//時間軸パラメータクラス
private:
	double dt;
	int nend;
	int nsample;

public:
	TimeP();
	TimeP(const TimeP& Tp);
	TimeP& operator=(const TimeP& Tp);
	double getDt();
	int getNend();
	int getNsample();
	void setTparam();
};
*/
class Boundarycond {//境界条件クラス
private:
	int flagL, flagR, flagU, flagD;//壁面の境界条件
	int flagC;//障害物壁面
	//ncond:  0:内部,	1:剛体壁面内部(流入流出なし),
			//2:流入壁面(dirichlet)　3:流出壁面(neumann)
			//4:移動壁面条件(壁面接線方向に流速固定値1),5:滑りなし壁面条件(壁面において(u,v)=(0,0))
			//6:滑りあり壁面条件(dvx/dy=0,vy=0)
	//double dLx, dRx, dUx, dDx;//境界条件値x方向
	//double dLy, dRy, dUy, dDy;//境界条件値y方向

public:
	Boundarycond();
	void set_cavityBC();//キャビティ流れの境界条件
	void set_cylinderBC();//円柱流れの境界条件
	void set_backstepBC();//バックステップ流れの境界条件
	void set_userBC();//任意の境界条件
	void setBCflagL(int flag);
	void setBCflagR(int flag);
	void setBCflagD(int flag);
	void setBCflagU(int flag);
	void setBCflagC(int flag);
	int getBCflagL();
	int getBCflagR();
	int getBCflagD();
	int getBCflagU();
	int getBCflagC();
	
};


/*
class ADeq_param_2d {//一次元移流拡散方程式パラメータクラス
private:
	double cx;//定常流速x成分
	double cy;//定常流速x成分
	double alpha;//拡散係数
	double courantx;//クーラン数x方向
	double couranty;//クーラン数y方向
	double diffusion;//拡散数
	double Pe;//ペクレ数
	NodeP& nparam;
	TimeP& tparam;
public:
	//ADeq_param_2d();
	//ADeq_param_2d(TimeP& Tp);
	ADeq_param_2d(NodeP& Np, TimeP& Tp);
	void set_param();
	double get_cx();
	double get_cy();
	double get_alpha();
	double get_couranx();
	double get_courany();
	double get_diffusion();
	double get_Pe();
};
*/


