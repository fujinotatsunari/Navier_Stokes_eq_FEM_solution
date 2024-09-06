#include "value.hpp"
#include "mesh.hpp"
#include "input.hpp"
#include <cmath>
#include <vector>
#include <iostream>
using namespace std;
# define PI 3.14159265359

Scalar2d::Scalar2d() :v(0) {}
Scalar2d::Scalar2d(double V) :v(V) {}

double& Scalar2d::operator[](int i) {
	return v;
}
const double& Scalar2d::operator[](int i)const {
	return v;
}
Scalar2d& Scalar2d::operator=(const Scalar2d& V) {

	v = V.v;
	return *this;
}
Scalar2d& Scalar2d::operator=(const double V) {
	v = V;
	return *this;
}
Scalar2d& Scalar2d::operator+=(const Scalar2d& V) {

	v += V.v;
	return *this;
}
Scalar2d& Scalar2d::operator-=(const Scalar2d& V) {

	v -= V.v;
	return *this;
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
Vector2d::Vector2d() :x_(0),y_(0){}
Vector2d::Vector2d(double X,double Y):x_(X),y_(Y){}
Vector2d::Vector2d(const Vector2d& V) {
	x_ = V.x_;
	y_ = V.y_;
}
double& Vector2d::operator[](int i) {
	if (i == 0) {//x成分を返す
		return x_;
	}
	else if (i == 1) {
		return y_;//y成分を返す
	}
	else {
		cout << "ベクトルの要素範囲外です" << endl;
		exit(-1);
	}
}
const double& Vector2d::operator[](int i) const {
	if (i == 0) {//x成分を返す
		return x_;
	}
	else if (i == 1) {
		return y_;//y成分を返す
	}
	else {
		cout << "ベクトルの要素範囲外です" << endl;
		exit(-1);
	}
}
Vector2d& Vector2d::operator=(const Vector2d& V) {
	x_ = V.x_;
	y_ = V.y_;

	return *this;
}
Vector2d& Vector2d::operator+=(const Vector2d& V) {
	x_ += V.x_;
	y_ += V.y_;

	return *this;
}
Vector2d& Vector2d::operator-=(const Vector2d& V) {
	x_ -= V.x_;
	y_ -= V.y_;

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
void Pressure::input(InputData& input) {
	for (int ie = 0; ie < nelem; ie++) {
		scalar[ie] = input.getP()[ie];
	}
}
/*
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
*/
PHI::PHI(Mesh2d& Mesh, Boundarycond& BC)
	:ScalarField2d(Mesh, BC),nnode(Mesh.nnode())
{
	scalar.resize(nnode);
}
void PHI::init() {
	for (int j = 0; j < mesh.ynode(); j++) {
		for (int i = 0; i < mesh.xnode(); i++) {
			int np = i + mesh.xnode() * j;
			scalar[np][0] = 0.0;

		}
	}
}
Velocity2d::Velocity2d(Mesh2d& Mesh, Boundarycond& BC) 
	:VectorField2d(Mesh, BC), nnode(Mesh.nnode())
{
	vector.resize(nnode);
}

void Velocity2d::input(InputData& input) {
	for (int np = 0; np < nnode; np++) {
		Vector2d V(input.getUx()[np], input.getUy()[np]);
		vector[np] = V;
	}
}
/*
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
	Vector2d V(1.0, 0.0);
	for (int j = 0; j < mesh.ynode(); j++) {
		for (int i = 0; i < mesh.xnode(); i++) {
			int np = i + mesh.xnode() * j;
			if (mesh.ncond(np) == 2) {//2:流入壁面(壁面法線方向に流速固定値1)
				vector[np] = V;
			}
		}
	}
	

}*/

/*
PHI::PHI(Mesh2d& mesh, Boundarycond& BC)
	:ScalarField2d(mesh, BC)
{
	setup();
}
void PHI::setup() {
	scalar_.resize(mesh_.nnode());
	value.resize(mesh_.nnode());
	for (int i = 0; i < scalar_.size(); i++) {
		scalar_[i].setNo(i);
		scalar_[i].setX(mesh_.x(i));
		scalar_[i].setY(mesh_.y(i));
		scalar_[i].setV(0.0);
		value[i] = scalar_[i].getV();
	}
}
void PHI::initialize_default() {
	//初期条件の設定
	for (int j = 0; j < mesh_.ynode(); j++) {
		for (int i = 0; i < mesh_.xnode(); i++) {
			int np = i + mesh_.xnode() * j;

			if (i == 0) {//左端
				value[np] = Bcond_.getdL();;
				scalar_[np].setV(value[np]);
			}
			else if (j == 0) {//下端
				value[np] = Bcond_.getdD();
				scalar_[np].setV(value[np]);
			}
			else {
				value[np] = 0.0;
				scalar_[np].setV(value[np]);
			}
			//cout << "value[" << np << "]=" << value[np] << endl;
		}
	}
	//境界条件の設定
	for (int j = 0; j < mesh_.ynode(); j++) {
		for (int i = 0; i < mesh_.xnode(); i++) {
			int np = i + mesh_.xnode() * j;
			if (mesh_.ncond(np) == 1) {//dirichlet境界条件
				if (i == 0) {//左壁面
					value[np] = Bcond_.getdL();
					scalar_[np].setV(value[np]);
				}
				if (j == 0) {//下壁面
					value[np] = Bcond_.getdD();
					scalar_[np].setV(value[np]);
				}
				if (i == mesh_.xnode() - 1) {//右壁面
					//value[np] = Bcond_.getdR();
					//scalar_[np].setV(value[np]);
				}
				if (j == mesh_.ynode() - 1) {//上壁面
					//value[np] = Bcond_.getdU();
					//scalar_[np].setV(value[np]);
				}

			}
		}
	}
	
}
*/
