#include "value.hpp"
#include "mesh.hpp"
#include "input.hpp"
#include <cmath>
#include <vector>
#include <iostream>
using namespace std;


Scalar2d::Scalar2d() :V(0) {}
Scalar2d::Scalar2d(double val) :V(val)
{}

double& Scalar2d::operator[](int i) {
	return V;
}
const double& Scalar2d::operator[](int i)const {
	return V;
}
double Scalar2d::v() {
	return V;
}
Scalar2d& Scalar2d::operator=(const Scalar2d& val) {

	V = val.V;
	return *this;
}
Scalar2d& Scalar2d::operator=(const double val) {
	V = val;
	return *this;
}
Scalar2d& Scalar2d::operator+=(const Scalar2d& val) {

	V += val.V;
	return *this;
}
Scalar2d& Scalar2d::operator-=(const Scalar2d& val) {

	V -= val.V;
	return *this;
}
Scalar2d Scalar2d::operator-() const {
	Scalar2d S;
	S.V = -V;
	return S;
}
Scalar2d operator+(const Scalar2d& a, const Scalar2d& b) {
	Scalar2d val;
	val[0] = a[0] + b[0];
	return val;
}
Scalar2d operator-(const Scalar2d& a, const Scalar2d& b) {
	Scalar2d val;
	val[0] = a[0] - b[0];
	return val;
}
Scalar2d operator*(const Scalar2d& a, const Scalar2d& b) {
	Scalar2d val;
	val[0] = a[0] * b[0];
	return val;
}
Scalar2d operator*(const double k, const Scalar2d& V) {
	Scalar2d val;
	val[0] = k * V[0];
	return val;
}
Scalar2d operator/(const Scalar2d& V, const double k) {
	Scalar2d val;
	val[0] = V[0] / k;
	return val;
}
ScalarField2d::ScalarField2d(Mesh2d& Mesh, Boundarycond& BC)
	:mesh(Mesh), Bcond(BC), size(Mesh.nnode()){}

const Scalar2d& ScalarField2d::operator[](int i)const {
	return scalar[i];
}
Scalar2d& ScalarField2d::operator[](int i) {
	return scalar[i];
}

ScalarField2d& ScalarField2d::operator=(const ScalarField2d& S) {
	size = S.size;
	mesh = S.mesh;
	Bcond = S.Bcond;
	scalar = S.scalar;

	return *this;
}

Vector2d::Vector2d() :x(0),y(0){}
Vector2d::Vector2d(double X,double Y):x(X),y(Y){}
Vector2d::Vector2d(const Vector2d& V) {
	x = V.x;
	y = V.y;
}
double& Vector2d::operator[](int i) {
	if (i == 0) {//x成分を返す
		return x;
	}
	else if (i == 1) {
		return y;//y成分を返す
	}
	else {
		cout << "ベクトルの要素範囲外です" << endl;
		exit(-1);
	}
}
const double& Vector2d::operator[](int i) const {
	if (i == 0) {//x成分を返す
		return x;
	}
	else if (i == 1) {
		return y;//y成分を返す
	}
	else {
		cout << "ベクトルの要素範囲外です" << endl;
		exit(-1);
	}
}
double Vector2d::norm() const {
	return sqrt(x * x + y * y);
}
Vector2d& Vector2d::operator=(const Vector2d& V) {
	x = V.x;
	y = V.y;

	return *this;
}
Vector2d& Vector2d::operator+=(const Vector2d& V) {
	x += V.x;
	y += V.y;

	return *this;
}
Vector2d& Vector2d::operator-=(const Vector2d& V) {
	x -= V.x;
	y -= V.y;

	return *this;
}
Vector2d operator+(const Vector2d& v1, const Vector2d& v2) {
	Vector2d vec;
	vec[0] = v1[0] + v2[0];
	vec[1] = v1[1] + v2[1];
	return vec;
}
Vector2d operator-(const Vector2d& v1, const Vector2d& v2) {
	Vector2d vec;
	vec[0] = v1[0] - v2[0];
	vec[1] = v1[1] - v2[1];
	return vec;
}
Vector2d operator*(const Vector2d& v1, const Vector2d& v2) {
	Vector2d vec;
	vec[0] = v1[0] * v2[0];
	vec[1] = v1[1] * v2[1];
	return vec;
}
Vector2d operator*(const double k, const Vector2d& V) {
	Vector2d vec;
	vec[0] = k * V[0];
	vec[1] = k * V[1];
	return vec;
}
Vector2d operator*(const Scalar2d& k, const Vector2d& V) {
	Vector2d vec;
	vec[0] = k[0] * V[0];
	vec[1] = k[0] * V[1];
	return vec;
}
Vector2d operator*(const Vector2d& V, const double k) {
	Vector2d vec;
	vec[0] = k * V[0];
	vec[1] = k * V[1];
	return vec;
}
Vector2d operator*(const Vector2d& V, const Scalar2d& k) {
	Vector2d vec;
	vec[0] = k[0] * V[0];
	vec[1] = k[0] * V[1];
	return vec;
}
Scalar2d operator%(const Vector2d& v1, const Vector2d& v2) {
	Scalar2d val;
	val[0] = v1[0] * v2[1] - v2[0] * v1[1];
	return val;
}
Vector2d operator/(const Vector2d& V, const double k) {
	Vector2d vec;
	vec[0] = V[0] / k;
	vec[1] = V[1] / k;
	return vec;
}

VectorField2d::VectorField2d(Mesh2d& Mesh, Boundarycond& BC) 
	:mesh(Mesh), Bcond(BC), size(Mesh.nnode())
{}
const Vector2d& VectorField2d::operator[](int i)const {
	return vector[i];
}
Vector2d& VectorField2d::operator[](int i) {
	return vector[i];
}
Pressure::Pressure(Mesh2d& Mesh, Boundarycond& BC)
	:ScalarField2d(Mesh, BC), nelem(Mesh.nelem())
{
	cout << "Object generate: Pressure " << endl;
	scalar.resize(nelem);
	
}
void Pressure::input(InputData& input) {
	for (int ie = 0; ie < nelem; ie++) {
		scalar[ie] = input.getP()[ie];
	}
}
Pressure Pressure::nodeP() {
	Pressure nodeP(mesh, Bcond);
	nodeP.scalar.resize(mesh.nnode());

	//圧力は要素内中心で定義され　要素内一定となる
	//節点値圧力を求めるために 要素ごとに節点に圧力を足し込んでいき
	//最後に平均する

	int node = 4;//要素内節点数
	for (int j = 0; j < node; j++) {//節点ループ
		for (int ie = 0; ie < mesh.nelem(); ie++) {
			int np = mesh.nbool1(ie, j);
			int i1 = mesh.i1(ie);
			int i2 = mesh.i2(ie);
			int i3 = mesh.i3(ie);
			int i4 = mesh.i4(ie);

			nodeP[np][0] = nodeP[np][0] + scalar[ie].v();//圧力の節点への足し込み
		}
	}
	for (int j = 0; j < mesh.ynode(); j++) {
		for (int i = 0; i < mesh.xnode(); i++) {
			int np = i + mesh.xnode() * j;
			if (i == 0 || i == mesh.xnode() - 1 || j == 0 || j == mesh.ynode() - 1) {
				if (i == 0) {//左壁面
					if (j != 0 && j != mesh.ynode() - 1) {//角を除く壁面
						nodeP[np][0] = nodeP[np][0] / 2.0;
					}
				}
				if (i == mesh.xnode() - 1) {//右壁面
					if (j != 0 && j != mesh.ynode() - 1) {//角を除く壁面
						nodeP[np][0] = nodeP[np][0] / 2.0;
					}
				}
				if (j == 0) {//下壁面
					if (i != 0 && i != mesh.xnode() - 1) {//角を除く壁面
						nodeP[np][0] = nodeP[np][0] / 2.0;
					}
				}
				if (j == mesh.ynode() - 1) {//上壁面
					if (i != 0 && i != mesh.xnode() - 1) {//角を除く壁面
						nodeP[np][0] = nodeP[np][0] / 2.0;
					}
				}
			}
			else {//内部
				nodeP[np][0] = nodeP[np][0] / 4.0;
			}

		}
	}
	return nodeP;
}
Velocity2d::Velocity2d(Mesh2d& Mesh, Boundarycond& BC) 
	:VectorField2d(Mesh, BC), nnode(Mesh.nnode())
{
	cout << "Object generate: Velocity2d " << endl;
	vector.resize(nnode);
}

void Velocity2d::input(InputData& input) {
	for (int np = 0; np < nnode; np++) {
		Vector2d V(input.getUx()[np], input.getUy()[np]);
		vector[np] = V;
	}
}
