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
	Mesh2d& mesh_;

public:
	CofficientMatrix(Mesh2d& mesh);//コンストラクタ
	CofficientMatrix(Mesh2d& mesh, const vector<Matrix>& Mat);

	static CofficientMatrix generate(Mesh2d& mesh);//行列の生成
	const Matrix& operator[](int ie)const;//添字演算子
	Matrix& operator[](int ie);//添字演算子
	void view();//行列の表示

};


class Massmatrix :public CofficientMatrix {//質量行列
public:
	Massmatrix(Mesh2d& mesh);
	Massmatrix(Mesh2d& mesh, const vector<Matrix>& Mat);
	static Massmatrix generate_Mass(Mesh2d& mesh);
	//const Matrix& operator[](int ie)const override;
	//Matrix& operator[](int ie) override;
};
class Lumped_Massmatrix :public CofficientMatrix {//集中化質量行列
public:
	Lumped_Massmatrix(Mesh2d& mesh);
	Lumped_Massmatrix(Mesh2d& mesh, const vector<Matrix>& Mat);
	static Lumped_Massmatrix generate_Lmass(Mesh2d& mesh);
	//const Matrix& operator[](int ie)const override;
	//Matrix& operator[](int ie) override;
};

class Diffmatrix :public CofficientMatrix {//拡散行列
public:
	Diffmatrix(Mesh2d& mesh);
	Diffmatrix(Mesh2d& mesh, const vector<Matrix>& Mat);
	static Diffmatrix generate_Diff(Mesh2d& mesh);
	//const Matrix& operator[](int ie)const override;
	//Matrix& operator[](int ie) override;
};

class xAdvecmatrix :public CofficientMatrix {//x方向移流行列クラス(線形)
public:
	xAdvecmatrix(Mesh2d& mesh);
	xAdvecmatrix(Mesh2d& mesh, const vector<Matrix>& Mat);
	static xAdvecmatrix generate_xAd(Mesh2d& mesh);
	//const Matrix& operator[](int ie)const override;
	//Matrix& operator[](int ie) override;
};
class yAdvecmatrix :public CofficientMatrix {//x方向移流行列クラス(線形)
public:
	yAdvecmatrix(Mesh2d& mesh);
	yAdvecmatrix(Mesh2d& mesh, const vector<Matrix>& Mat);
	static yAdvecmatrix generate_yAd(Mesh2d& mesh);
	//const Matrix& operator[](int ie)const override;
	//Matrix& operator[](int ie) override;
};