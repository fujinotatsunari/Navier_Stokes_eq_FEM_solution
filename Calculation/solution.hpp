#pragma once
#include"solution.hpp"
#include"output.hpp"
#include"FEM.hpp"
#include"value.hpp"
#include"Mesh.hpp"
#include"param.hpp"
#include"matrix.hpp"
#include<vector>
class Explicit_FEM {//陽解法
private:
	Mesh2d& mesh;
	Boundarycond& BC;
	Time& t;
	PHI& phi;
	ADeq_param_2d& ADP;

	int node = 4;
public:
	Explicit_FEM(Mesh2d& mesh_, Time& t_, PHI& phi_, Boundarycond& BC_, ADeq_param_2d& adp_);

	void do_expcalculation();//質量集中化を用いる計算

};

class Implicit_FEM {//陰解法
private:
	Mesh2d& mesh;
	Boundarycond& BC;
	Time& t;
	PHI& phi;
	ADeq_param_2d& ADP;

	int node = 4;

	//Matrix A;//行列方程式Ax=bを解く
	//vector<double> b;
	//vector<double> x;
public:
	Implicit_FEM(Mesh2d& mesh_, Time& t_, PHI& phi_, Boundarycond& BC_, ADeq_param_2d& adp_);

	void do_impcaluculation();//陰解法

};
