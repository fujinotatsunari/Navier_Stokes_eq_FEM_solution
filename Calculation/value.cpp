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
