#pragma once
#include "SOR.hpp"
#include "mesh.hpp"
#include "value.hpp"
#include "FEM.hpp"
#include <vector>
using namespace std;

class SORparam {//同時緩和法のパラメータ
private:
	Mesh2d& mesh;
	int nmax = 500;//同時緩和法最大反復回数
	double eps = 1.0e-5;//同時緩和法収束判定値
	vector<double> lambda;//速度修正ポテンシャル係数
	

public:
	SORparam(Mesh2d& Mesh);
	int get_nmax();
	double get_eps();
	double get_lambda(int ie);
	
};