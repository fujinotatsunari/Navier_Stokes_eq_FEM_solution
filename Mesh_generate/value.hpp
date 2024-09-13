#pragma once
#include"value.hpp"
#include"Mesh.hpp"
#include<vector>

class Scalar2d {//二次元スカラー量クラス
private:
	double V;//スカラー量
public:
	Scalar2d();
	Scalar2d(double V);
	double v();
	double& operator[](int i);//任意の整数値を添え字に取る
	const double& operator[](int i) const;
	Scalar2d operator-() const;
	Scalar2d& operator=(const Scalar2d& V);
	Scalar2d& operator=(const double V);
	Scalar2d& operator+=(const Scalar2d& V);
	Scalar2d& operator-=(const Scalar2d& V);
};
Scalar2d operator+(const Scalar2d&, const Scalar2d&);
Scalar2d operator-(const Scalar2d&, const Scalar2d&);
Scalar2d operator*(const Scalar2d&, const Scalar2d&);

Scalar2d operator*(const double, const Scalar2d& V);
Scalar2d operator/(const Scalar2d& V, const double k);



class ScalarField2d {//二次元スカラー場クラス
	//Mesh2d上でスカラー場を定める.
protected:
	int size;//scalar_の要素数
	Mesh2d& mesh;
	Boundarycond& Bcond;
	vector<Scalar2d> scalar;

public:
	ScalarField2d(Mesh2d& Mesh, Boundarycond& BC);

	const Scalar2d& operator[](int i)const;
	Scalar2d& operator[](int i);
};

class Vector2d {
private:
	double x, y;
public:
	Vector2d();
	Vector2d(double X, double Y);
	Vector2d(const Vector2d& V);
	double& operator[](int i);
	double const& operator[](int i) const;
	Vector2d& operator=(const Vector2d& V);
	Vector2d& operator+=(const Vector2d& V);
	Vector2d& operator-=(const Vector2d& V);
};

Vector2d operator+(const Vector2d&, const Vector2d&);
Vector2d operator-(const Vector2d&, const Vector2d&);
Vector2d operator*(const Vector2d&, const Vector2d&);//内積
Vector2d operator*(const double k, const Vector2d& V);//スカラー積
Vector2d operator*(const Vector2d& V, const double k);//スカラー積
Vector2d operator*(const Scalar2d& k, const Vector2d& V);//スカラー積
Vector2d operator*(const Vector2d& V, const Scalar2d& k);//スカラー積
Scalar2d operator%(const Vector2d&, const Vector2d&);//クロス積(２次元のクロス積は実質スカラー)
Vector2d operator/(const Vector2d& V, const double k);

class VectorField2d {
protected:
	int size;//vector_の要素数
	Mesh2d& mesh;
	Boundarycond& Bcond;
	vector<Vector2d> vector;
public:
	VectorField2d(Mesh2d& Mesh, Boundarycond& BC);
	const Vector2d& operator[](int i)const;
	Vector2d& operator[](int i);

};
class Pressure :public ScalarField2d {
private:
	int nelem;//要素数
public:
	Pressure(Mesh2d& mesh, Boundarycond& BC);
	void init();//初期化
	void cavity_init();//キャビティ流れの初期化(圧力編)
	void backstep_init();//バックステップ流れの初期化
};
/*
class PHI :public ScalarField2d {//移流拡散方程式の物理量phi
private:
	int nnode;//va配列の大きさnnodeに対応
public:
	PHI(Mesh2d& mesh, Boundarycond& BC);
	void init();
	
	//void initialize_default();//左端と下端に1を与える
	//void initialize_deltafunc();//デルタ関数型の初期条件(
	//void initialize_cosfunc();//cos関数型の初期条件

};
*/
class Velocity2d :public VectorField2d {
private:
	int nnode;
public:
	Velocity2d(Mesh2d& mesh, Boundarycond& BC);
	void init();//初期化
	void cavity_init();//キャビティ流れの初期化(流速編)
	void backstep_init();//バックステップ流れの初期化
};


