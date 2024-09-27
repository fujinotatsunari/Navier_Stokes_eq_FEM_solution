#include"value.hpp"
#include"mesh.hpp"
#include<cmath>
#include<vector>
#include<iostream>
using namespace std;
# define PI 3.14159265359

Scalar2d::Scalar2d() :V(0) {}
Scalar2d::Scalar2d(double val) :V(val) {}

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
}
const double& Vector2d::operator[](int i) const {
	if (i == 0) {//x成分を返す
		return x;
	}
	else if (i == 1) {
		return y;//y成分を返す
	}
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
Vector2d operator*(const Vector2d& V, const double k) {
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
	scalar.resize(nelem);
	
}
void Pressure::init() {
	
	for (int j = 0; j < mesh.yelem(); j++) {
		for (int i = 0; i < mesh.xelem(); i++) {
			int ie = i + mesh.xelem() * j;
			scalar[ie] = 0.0;

		}
	}
}
void Pressure::cavity_init() {
	init();

}
void Pressure::backstep_init() {
	init();
}
void Pressure::Pillar_init() {
	init();
}

Velocity2d::Velocity2d(Mesh2d& Mesh, Boundarycond& BC) 
	:VectorField2d(Mesh, BC), nnode(Mesh.nnode())
{
	vector.resize(nnode);
}
void Velocity2d::init() {
	Vector2d V(0.0, 0.0);
	for (int j = 0; j < mesh.ynode(); j++) {
		for (int i = 0; i < mesh.xnode(); i++) {
			int np = i + mesh.xnode() * j;
			vector[np] = V;

		}
	}

}
void Velocity2d::cavity_init() {
	init();
	Vector2d V(1.0, 0.0);
	for (int j = 0; j < mesh.ynode(); j++) {
		for (int i = 0; i < mesh.xnode(); i++) {
			int np = i + mesh.xnode() * j;
			if (mesh.ncond(np) == 4) {//4:移動壁面条件(壁面接線方向に流速固定値1)
				vector[np] = V;
			}
		}
	}
}
void Velocity2d::backstep_init() {
	init();
	/*
	std::vector<double> inlet;//流入領域
	double y0b, y0t, y0L, y0;//流入下端,上端,流入領域長さ,領域中点

	for (int j = 0; j < mesh.ynode(); j++) {
		for (int i = 0; i < mesh.xnode(); i++) {
			int np = i + mesh.xnode() * j;
			if (mesh.ncond(np) == 2) {//2:流入壁面
				inlet.push_back(mesh.y(np));
			}
		}
	}
	y0b = inlet[0];
	y0t = inlet.back();
	y0L = y0t - y0b;
	y0 = (y0b + y0t) / 2;
	double C = y0b * y0b - 2 * y0b * y0 + y0 * y0;
	//ポアズイユ流れの流入条件
	//u = -0.5/C (y-y0)^2 + 0.5
	for (int j = 0; j < mesh.ynode(); j++) {
		for (int i = 0; i < mesh.xnode(); i++) {
			int np = i + mesh.xnode() * j;
			if (mesh.ncond(np) == 2) {//2:流入壁面(壁面法線方向に流速固定値1)
				double u = (-0.5 / C) * (mesh.y(np) - y0) * (mesh.y(np) - y0) + 0.5;
				Vector2d V(u, 0);
				vector[np] = V;
			}
		}
	}
	*/
}
void Velocity2d::Pillar_init() {
	init();
	/*
	std::vector<double> inlet;//流入領域
	double y0b, y0t, y0L, y0;//流入下端,上端,流入領域長さ,領域中点

	for (int j = 0; j < mesh.ynode(); j++) {
		for (int i = 0; i < mesh.xnode(); i++) {
			int np = i + mesh.xnode() * j;
			if (mesh.ncond(np) == 2) {//2:流入壁面
				inlet.push_back(mesh.y(np));
			}
		}
	}
	y0b = inlet[0];
	y0t = inlet.back();
	y0L = y0t - y0b;
	y0 = (y0b + y0t) / 2;
	double C = y0b * y0b - 2 * y0b * y0 + y0 * y0;
	//ポアズイユ流れの流入条件
	//u = -0.5/C (y-y0)^2 + 0.5
	for (int j = 0; j < mesh.ynode(); j++) {
		for (int i = 0; i < mesh.xnode(); i++) {
			int np = i + mesh.xnode() * j;
			if (mesh.ncond(np) == 2) {//2:流入壁面(壁面法線方向に流速固定値1)
				double u = (-0.5 / C) * (mesh.y(np) - y0) * (mesh.y(np) - y0) + 0.5;
				Vector2d V(u, 0);
				vector[np] = V;
			}
		}
	}
	*/
}

