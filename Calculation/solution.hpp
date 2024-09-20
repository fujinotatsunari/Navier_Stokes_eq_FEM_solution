#pragma once
#include"solution.hpp"
#include"output.hpp"
#include"FEM.hpp"
#include"input.hpp"
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
	Boundarycond& BC;
	InputData& input;
	int NOR;//同時緩和法反復回数
	double max_div;//発散量の最大値

public:
	HSMAC_FEM(Velocity2d& v, Pressure& p, Time& T, Mesh2d& Mesh, NDNSparam& NSP, SORparam& SRP, Boundarycond& bc, InputData& INPUT);
	void do_solution();//HSMAC法による計算の開始
	double Uxmax();//流速x成分最大値(絶対値が一番大きい値の絶対値をつける前の値を返す)
	double Uymax();//流速y成分最大値(絶対値が一番大きい値の絶対値をつける前の値を返す)
	double Pmax();//圧力P最大値(絶対値が一番大きい値の絶対値をつける前の値を返す)
	void view_parameters(int n);//パラメーターの表示
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
	double div_max;//発散量最大値(絶対値)(これがepsより小さくなったらループを突破)
public:
	SOR(SORparam& Param);
	void do_calculation(Velocity2d& v, Pressure& p, Time& T, Mesh2d& Mesh, SORparam& param, Boundarycond& BC);
	int get_nor();
	double max_div();//発散量の最大値を返す

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

//作るやつ
//完全陰解法+直接法+人工圧縮性法 FID-ACmethod
//半陰解法(移流項の線形化)+直接法+人工圧縮性法 SID-ACmehtod
//SIMPLER法(burgers+poisson分離陰解法)
