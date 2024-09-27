#include "mesh.hpp"
#include "param.hpp"
#include <vector>
#include <iostream>
#include <cmath>
using namespace std;

Node2d::Node2d() :no(0), x(0), y(0){}
void Node2d::setNo(int no_) {
	no = no_;
}
int Node2d::getNo() {
	return no;
}
void Node2d::setX(double x_) {
	x = x_;
}
double Node2d::getX() {
	return x;
}
void Node2d::setY(double y_) {
	y = y_;
}
double Node2d::getY() {
	return y;
}
Element2d::Element2d() :no(0), x(0), y(0), Se(0){}
void Element2d::setNo(int no_) {
	no = no_;
}
int Element2d::getNo() {
	return no;
}
void Element2d::setX(double x_) {
	x = x_;
}
double Element2d::getX() {
	return x;
}
void Element2d::setY(double y_) {
	y = y_;
}
double Element2d::getY() {
	return y;
}
void Element2d::setSe(double Se_) {
	Se = Se_;
}
double Element2d::getSe() {
	return Se;
}

Time::Time(TimeP& Tp) :ntime_(0), dt_(0), nend_(0), nsample_(0), tparam(Tp) {
	cout << "Object generate :Time" << endl;
	setup();
}
void Time::setup() {
	dt_ = tparam.getDt();
	nend_ = tparam.getNend();
	nsample_ = tparam.getNsample();
	t.resize(nend_ + 1, 0);
	for (int n = 0; n < t.size(); n++) {
		t[n] = 0.0 + (double)n * dt_;
	}
}
double Time::ntime(int n) {
	n_ = n;
	ntime_ = 0.0 + (double)n * dt_;
	return ntime_;
}
double Time::ntime() {
	ntime_ = 0.0 + (double)n_ * dt_;
	return ntime_;
}
void Time::setn(int n) {
	n_ = n;
}
double Time::dt() {
	return dt_;
}
int Time::nend() {
	return nend_;
}
int Time::nsample() {
	return nsample_;
}
double const& Time::operator[](int n)const {
	return t[n];
}
double& Time::operator[](int n) {
	return t[n];
}
Mesh2d::Mesh2d(InputData& Input)
	:nparam_(Input.get_NodeParam()), Bcond_(Input.get_BC()), Lx_(0), Ly_(0), xb_(0), xt_(0), yb_(0), yt_(0), dx_(0), dy_(0), xnode_(0), ynode_(0), xelem_(0), yelem_(0), nnode_(0), nelem_(0)
{
	cout << "Object generate: Mesh2d(input) " << endl;
	
	//setup();
	xb_ = nparam_.getXb();
	xt_ = nparam_.getXt();
	yb_ = nparam_.getYb();
	yt_ = nparam_.getYt();
	dx_ = nparam_.getDx();
	dy_ = nparam_.getDy();
	Lx_ = nparam_.getLx();
	Ly_ = nparam_.getLy();

	xnode_ = nparam_.getXnode();
	ynode_ = nparam_.getYnode();
	nnode_ = nparam_.getNnode();

	xelem_ = nparam_.getXelem();
	yelem_ = nparam_.getYelem();
	nelem_ = nparam_.getNelem();

	nbool1_ = Input.getnbool1();
	nbool3_ = Input.getnbool3();
	ncond_ = Input.getcond();
	scond_ = Input.getscond();
	X = Input.getx();
	Y = Input.gety();
	EX = Input.getex();
	EY = Input.getey();

}
Mesh2d::Mesh2d(NodeP& NP, Boundarycond& BC)
	:nparam_(NP), Bcond_(BC), Lx_(0), Ly_(0), xb_(0), xt_(0), yb_(0), yt_(0), dx_(0), dy_(0), xnode_(0), ynode_(0), xelem_(0), yelem_(0), nnode_(0), nelem_(0)
{
	cout << "Object generate: Mesh2d " << endl;
	xb_ = nparam_.getXb();
	xt_ = nparam_.getXt();
	yb_ = nparam_.getYb();
	yt_ = nparam_.getYt();
	dx_ = nparam_.getDx();
	dy_ = nparam_.getDy();
	Lx_ = nparam_.getLx();
	Ly_ = nparam_.getLy();

	xnode_ = nparam_.getXnode();
	ynode_ = nparam_.getYnode();
	nnode_ = nparam_.getNnode();

	xelem_ = nparam_.getXelem();
	yelem_ = nparam_.getYelem();
	nelem_ = nparam_.getNelem();

}

Mesh2d& Mesh2d::operator=(const Mesh2d& mesh) {
	nparam_ = mesh.nparam_;
	Bcond_ = mesh.Bcond_;
	node_ = mesh.node_;
	elem_ = mesh.elem_;
	ncond_ = mesh.ncond_;
	scond_ = mesh.scond_;
	nbool1_ = mesh.nbool1_;
	nbool3_ = mesh.nbool3_;
	xb_ = mesh.xb_;
	xt_ = mesh.xt_;
	yb_ = mesh.yb_;
	yt_ = mesh.yt_;
	dx_ = mesh.dx_;
	dy_ = mesh.dy_;
	Lx_ = mesh.Lx_;
	Ly_ = mesh.Ly_;
	xnode_ = mesh.xnode_;
	ynode_ = mesh.ynode_;
	nnode_ = mesh.nnode_;

	xelem_ = mesh.xelem_;
	yelem_ = mesh.yelem_;
	nelem_ = mesh.nelem_;

	return *this;
}
Mesh2d::Mesh2d(const Mesh2d& mesh)
	:nparam_(mesh.nparam_), Bcond_(mesh.Bcond_)
{
	cout << "Object generate :Mesh2d " << endl;

	node_ = mesh.node_;
	elem_ = mesh.elem_;
	ncond_ = mesh.ncond_;
	scond_ = mesh.scond_;
	nbool1_ = mesh.nbool1_;
	nbool3_ = mesh.nbool3_;
	xb_ = mesh.xb_;
	xt_ = mesh.xt_;
	yb_ = mesh.yb_;
	yt_ = mesh.yt_;
	dx_ = mesh.dx_;
	dy_ = mesh.dy_;
	Lx_ = mesh.Lx_;
	Ly_ = mesh.Ly_;
	xnode_ = mesh.xnode_;
	ynode_ = mesh.ynode_;
	nnode_ = mesh.nnode_;

	xelem_ = mesh.xelem_;
	yelem_ = mesh.yelem_;
	nelem_ = mesh.nelem_;
}
void Mesh2d::geninputmesh() {
	node_.resize(nnode_);
	elem_.resize(nelem_);
	for (int np = 0; np < nnode_; np++) {
		node_[np].setNo(np);
		node_[np].setX(X[np]);
		//cout << "x[" << np << "]=" << X[np] << endl;
		node_[np].setY(Y[np]);
		//cout << "y[" << np << "]=" << Y[np] << endl;
	}
	for (int ie = 0; ie < nelem_; ie++) {
		elem_[ie].setNo(ie);
		int i1 = nbool1_[ie][0];
		int i2 = nbool1_[ie][1];
		int i3 = nbool1_[ie][2];
		int i4 = nbool1_[ie][3];

		elem_[ie].setX(EX[ie]);//要素重心のx座標
		elem_[ie].setY(EY[ie]);//要素重心のy座標
		double S;
		S = area(i1, i2, i3, i4);
		elem_[ie].setSe(S);
		//cout << "S[" << ie << "]=" << S << endl;
		//cout << "EX[" << ie << "]=" << EX[ie] << "EY[" << ie << "]=" << EY[ie] << endl;
	}

	/*
	for (int ie = 0; ie < nelem_; ie++) {
		scond_[ie] = 0;
	}
	for (int ie = 0; ie < nelem_; ie++) {
		int i1 = nbool1_[ie][0];
		int i2 = nbool1_[ie][1];
		int i3 = nbool1_[ie][2];
		int i4 = nbool1_[ie][3];

		if (ncond(i1) == 1 || ncond(i2) == 1 || ncond(i3) == 1 || ncond(i4) == 1) {
			//要素内の点に剛体内部のフラグが立つ
			scond_[ie] = 1;
		}
		if (ncond(i1) != 0 && ncond(i2) != 0 && ncond(i3) != 0 && ncond(i4) != 0) {
			//要素内のすべての点が壁面のフラグ
			scond_[ie] = 1;
		}
		
	}
	*/
}


double Mesh2d::xb() {
	return xb_;
}
double Mesh2d::xt() {
	return xt_;
}
double Mesh2d::yb() {
	return yb_;
}
double Mesh2d::yt() {
	return yt_;
}
double Mesh2d::dx() {
	return dx_;
}
double Mesh2d::dy() {
	return dy_;
}
double Mesh2d::Lx() {
	return Lx_;
}
double Mesh2d::Ly() {
	return Ly_;
}
double Mesh2d::x(int i) {
	return node_[i].getX();
}
double Mesh2d::y(int i) {
	return node_[i].getY();
}
double Mesh2d::eX(int ie) {
	return elem_[ie].getX();
}
double Mesh2d::eY(int ie) {
	return elem_[ie].getY();
}
double Mesh2d::Se(int ie) {
	return elem_[ie].getSe();
}
int Mesh2d::xnode() {
	return xnode_;
}
int Mesh2d::ynode() {
	return ynode_;
}
int Mesh2d::xelem() {
	return xelem_;
}
int Mesh2d::yelem() {
	return yelem_;
}
int Mesh2d::nnode() {
	return nnode_;
}
int Mesh2d::nelem() {
	return nelem_;
}
int Mesh2d::i1(int ie) {
	return nbool1_[ie][0];
}
int Mesh2d::i2(int ie) {
	return nbool1_[ie][1];
}
int Mesh2d::i3(int ie) {
	return nbool1_[ie][2];
}
int Mesh2d::i4(int ie) {
	return nbool1_[ie][3];
}
int Mesh2d::nbool1(int ie, int np) {
	return nbool1_[ie][np];
}
int Mesh2d::e1(int ie) {
	return nbool3_[ie][0];
}
int Mesh2d::e2(int ie) {
	return nbool3_[ie][1];
}
int Mesh2d::e3(int ie) {
	return nbool3_[ie][2];
}
int Mesh2d::e4(int ie) {
	return nbool3_[ie][3];
}
int Mesh2d::nbool3(int ie, int np) {
	return nbool3_[ie][np];
}
int Mesh2d::ncond(int i) {
	if (i <= -1 || i >= nnode()) {
		return -1;
	}
	else {
		return ncond_[i];
	}

}
int Mesh2d::scond(int i) {
	if (i <= -1 || i >= nelem()) {
		return 1;
	}
	else {
		return scond_[i];
	}
	
}
double Mesh2d::length(int i1, int i2) {
	double L = 0.0;
	double x1 = node_[i1].getX();
	double x2 = node_[i2].getX();
	double y1 = node_[i1].getY();
	double y2 = node_[i2].getY();
	L = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
	return L;
}
double Mesh2d::area(int E1, int E2) {
	int e11 = nbool3(E1, 0);//E1の下の要素番号
	int e12 = nbool3(E1, 1);//E1の右の要素番号
	int e13 = nbool3(E1, 2);//E1の上の要素番号
	int e14 = nbool3(E1, 3);//E1の左の要素番号
	//要素重心の座標
	double x1 = eX(E1);
	double x3 = eX(E2);
	double y1 = eY(E1);
	double y3 = eY(E2);

	double x2, x4, y2, y4;//辺の両端の座標
	int n2, n4;//辺の両端の節点番号
	if (e11 == E2) {//E2はE1の下側
		n2 = nbool1(E1, 0);//要素左下の点
		n4 = nbool1(E1, 1);//要素右下の点
	}
	else if (e12 == E2) {//E2はE1の右側
		n2 = nbool1(E1, 1);//要素右下の点
		n4 = nbool1(E1, 2);//要素右上の点
	}
	else if (e13 == E2) {//E2はE1の上側
		n2 = nbool1(E1, 2);//要素右上の点
		n4 = nbool1(E1, 3);//要素左上の点
	}
	else if (e14 == E2) {//E2はE1の左側
		n2 = nbool1(E1, 3);//要素左上の点
		n4 = nbool1(E1, 0);//要素左下の点
	}
	else {
		cout << "E2はE1に隣接してません" << endl;
		exit(-1);
	}
	//辺両端の座標
	x2 = x(n2);
	y2 = y(n2);
	x4 = x(n4);
	y4 = y(n4);
	
	double area;
	area = ((x3 - x1) * (y4 - y2) - (x4 - x2) * (y3 - y1)) / 2;//要素四角形の面積を求める
	return area;

}
double Mesh2d::area(int n1, int n2, int n3, int n4) {

	double x1 = x(n1);
	double x2 = x(n2);
	double x3 = x(n3);
	double x4 = x(n4);
	double y1 = y(n1);
	double y2 = y(n2);
	double y3 = y(n3);
	double y4 = y(n4);
	//位置ベクトル
	double area;
	area = ((x3 - x1) * (y4 - y2) - (x4 - x2) * (y3 - y1)) / 2;//要素四角形の面積を求める
	return area;


}