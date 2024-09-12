#pragma once
#include"solution.hpp"
#include"output.hpp"
#include"FEM.hpp"
#include"SOR.hpp"
#include"value.hpp"
#include"Mesh.hpp"
#include"param.hpp"
#include"matrix.hpp"
#include<vector>

class HSMAC_FEM {//Highly Simplified Marker and Cell methodを背景にした有限要素法によるNS方程式の求解
private:	
	Velocity2d& V;
	Pressure& P;
	Time& t;
	Mesh2d& mesh;
	NDNSparam& nsparam;
	SORparam& sorparam;
	Boundarycond& bc;
	double coutant;

public:
	HSMAC_FEM(Velocity2d& v, Pressure& p, Time& T, Mesh2d& Mesh, NDNSparam& NSP, SORparam& SRP, Boundarycond& BC);
	double Uxmax();//流速x成分最大値
	double Uymax();//流速y成分最大値
	double Pmax();//圧力P最大値
	void view_parameters();//パラメーターの表示
	//

};

class Predictor {//予測子導出
private:
	int node = 4;
public:
	Predictor();
	void euler_explicit(Velocity2d& v, Pressure& p, Time& T, Mesh2d& Mesh, NDNSparam& Param);//前進オイラー法による予測子導出

};
class SOR {//同時緩和法反復
private:
	
	SORparam& sparam;

	int node = 4;
	int nor;//反復回数 Number of repetitions

public:
	SOR(SORparam& Param);
	void do_calculation(Velocity2d& v, Pressure& p, Time& T, Mesh2d& Mesh, SORparam& param, Boundarycond& BC);
	double get_nor();

};

class Divergence :public ScalarField2d {//発散量(要素内の流速の生成消滅)
private:
	int size;//=nelem 総要素数
	int node = 4;

public:
	Divergence(Mesh2d& Mesh, Boundarycond& bc);
	void cal_divergence(Velocity2d& v);
	void cal_divergence(vector<Vector2d>& v);
	double max_div();
};
class VMPotential :public ScalarField2d {//速度修正ポテンシャル　Velocity Modification Potential
private:
	int size;//=nelem 総要素数
public:
	VMPotential(Mesh2d& Mesh, Boundarycond& bc);
	void cal_VMP(Divergence& D, SORparam& param);

};

