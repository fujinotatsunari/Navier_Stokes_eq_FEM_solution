#pragma once
#include"FEM.hpp"
#include"value.hpp"
#include"Mesh.hpp"
#include"matrix.hpp"
#include<vector>


class CofficientMatrix {//要素内係数行列クラス
protected:
	int node = 4;//要素内節点数
	vector<Matrix> mat;//マトリクス本体:nelem個のvector
	vector<double> xi = { -1.0,1.0,1.0,-1.0 };//計算空間座標xi座標
	vector<double> eta = { -1.0,-1.0,1.0,1.0 };//計算空間座標eta座標
	Mesh2d& mesh;

public:
	CofficientMatrix(Mesh2d& mesh);//コンストラクタ
	const Matrix& operator[](int ie)const;//添字演算子
	Matrix& operator[](int ie);//添字演算子
	void view();//行列の表示

};


class Massmatrix :public CofficientMatrix {//質量行列
public:
	Massmatrix(Mesh2d& mesh);
};
class Lumped_Massmatrix :public CofficientMatrix {//集中化質量行列
public:
	Lumped_Massmatrix(Mesh2d& mesh);
};



class xAdvecmatrix :public CofficientMatrix {//x方向移流行列クラス(線形)
public:
	xAdvecmatrix(Mesh2d& mesh);
};
class yAdvecmatrix :public CofficientMatrix {//x方向移流行列クラス(線形)
public:
	yAdvecmatrix(Mesh2d& mesh);
};
class Advecmatrix :public CofficientMatrix {//非線形移流行列クラス
private:
	Velocity2d V;
public:
	Advecmatrix(Mesh2d& mesh, Velocity2d& v);
	void renew(Velocity2d& v);//流速変化による更新
};

class Diffmatrix :public CofficientMatrix {//拡散行列
public:
	Diffmatrix(Mesh2d& mesh);
};

class BTDmatrix :public CofficientMatrix {
private:
	Velocity2d V;
public:
	BTDmatrix(Mesh2d& mesh, Velocity2d& v);
	void renew(Velocity2d& v);//流速変化に伴う更新

};

class GradientVector {//要素内勾配ベクトルクラス
private:
	int node = 4;
	vector<double> xi = { -1.0,1.0,1.0,-1.0 };//計算空間座標xi座標
	vector<double> eta = { -1.0,-1.0,1.0,1.0 };//計算空間座標eta座標
	Mesh2d& mesh;
	vector<vector<Vector2d>> vec;//ie*4個のVector2d
public:
	GradientVector(Mesh2d& Mesh);
	void view();

	Vector2d i1(int ie);//要素内1節点におけるvecを返す
	Vector2d i2(int ie);
	Vector2d i3(int ie);
	Vector2d i4(int ie);

	vector<Vector2d>& operator[](int ie);//添字演算子オーバーロード
	const vector<Vector2d>& operator[](int ie) const;

	
};