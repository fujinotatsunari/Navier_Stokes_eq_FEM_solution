#include"mesh.hpp"
#include"param.hpp"
#include"value.hpp"
#include<vector>
#include<iostream>
Node2d::Node2d() :no(0), x(0), y(0)
{


}
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
Element2d::Element2d() :no(0), x(0), y(0), Se(0)
{

}
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
/*
Time::Time(TimeP& Tp) :ntime_(0), dt_(0), nend_(0), nsample_(0), tparam(Tp) {
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
	ntime_ = 0.0 + (double)n * dt_;
	return ntime_;
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
*/
Mesh2d::Mesh2d(NodeP& NP, Boundarycond& BC)
	:nparam_(NP), Bcond_(BC), Lx_(0), Ly_(0), xb_(0), xt_(0), yb_(0), yt_(0), dx_(0), dy_(0), xnode_(0), ynode_(0), xelem_(0), yelem_(0), nnode_(0), nelem_(0)
{
	setup();

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
void Mesh2d::setup() {
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

	int N = 4;//一要素の節点数
	nbool1_.resize(nelem_);
	for (int ie = 0; ie < nbool1_.size(); ie++) {
		nbool1_[ie].resize(N, 0);
	}
	int M = 4;//1要素に隣接する要素数
	nbool3_.resize(nelem_);
	for (int ie = 0; ie < nbool3_.size(); ie++) {
		nbool3_[ie].resize(M, 0);
	}

	ncond_.resize(nnode_);
	node_.resize(nnode_);
	scond_.resize(nelem_);
	elem_.resize(nelem_);
}
void Mesh2d::generate() {
	//節点座標の計算
	for (int j = 0; j < ynode_; j++) {
		for (int i = 0; i < xnode_; i++) {
			int np = i + xnode_ * j;
			node_[np].setNo(i);//節点番号の設定
			node_[np].setX(dx_ * (double)i + xb_);//節点座標xの設定
			node_[np].setY(dy_ * (double)j + yb_);//節点座標yの設定

			if (i == 0) node_[np].setX(xb_);//x左端の補正
			if (i == xnode_ - 1) node_[np].setX(xt_);//x右端の補正
			if (j == 0) node_[np].setY(yb_);//y下端の補正
			if (j == ynode_ - 1) node_[np].setY(yt_);//y上端の補正
		}
	}
	//nbool1とnbool3の割当
	for (int j = 0; j < yelem_; j++) {
		for (int i = 0; i < xelem_; i++) {
			int ie = i + xelem_ * j;

			int i1 = i + xnode_ * j;//要素の左下点の節点番号
			int i2 = i1 + 1;//要素の右下点の節点番号
			int i4 = i1 + xnode_;//要素の左上点の節点番号
			int i3 = i4 + 1;//要素の右上点の節点番号

			nbool1_[ie][0] = i1;
			nbool1_[ie][1] = i2;
			nbool1_[ie][2] = i3;
			nbool1_[ie][3] = i4;

			nbool3_[ie][2] = ie + xelem_;//上側要素
			if (j == yelem_ - 1) nbool3_[ie][2] = -1;//領域上端

			nbool3_[ie][0] = ie - xelem_;//下側要素
			if (j == 0) nbool3_[ie][0] = -1;//領域下端

			nbool3_[ie][3] = ie - 1;//左側要素
			if (i == 0) nbool3_[ie][3] = -1;//領域左端

			nbool3_[ie][1] = ie + 1;//右側要素
			if (i == xelem_ - 1) nbool3_[ie][1] = -1;//領域右端
		}
	}

	//要素座標の設定
	for (int j = 0; j < yelem_; j++) {
		for (int i = 0; i < xelem_; i++) {
			int ie = i + xelem_ * j;
			elem_[ie].setNo(ie);
			int i1 = nbool1_[ie][0];
			int i2 = nbool1_[ie][1];
			int i3 = nbool1_[ie][2];
			int i4 = nbool1_[ie][3];

			elem_[ie].setX(x(i1) + x(i2) + x(i3) + x(i4) / 4);//要素重心のx座標
			elem_[ie].setY(y(i1) + y(i2) + y(i3) + y(i4) / 4);//要素重心のy座標
			double S;//要素面積
			S = area(i1, i2, i3, i4);//算出
			elem_[ie].setSe(S);

		}
	}

	//境界条件フラグの設定
	//境界条件フラグの設定
	//ncond:  0:内部,	1:剛体壁面(流入流出なし),
			//2:流入壁面(dirichlet)　3:流出壁面(neumann)
			//4:移動壁面条件(壁面接線方向に流速固定値1),5:滑りなし壁面条件(壁面において(u,v)=(0,0))
			//6:滑りあり壁面条件(dvx/dy=0,vy=0)

	for (int j = 0; j < ynode_; j++) {
		for (int i = 0; i < xnode_; i++) {
			int np = i + xnode_ * j;
			if (i == 0) {//左壁面
				ncond_[np] = Bcond_.getBCflagL();
			}
			else if (j == 0) {//下壁面
				ncond_[np] = Bcond_.getBCflagD();
			}
			else if (i == xnode_ - 1) {//右壁面
				ncond_[np] = Bcond_.getBCflagR();
			}
			else if (j == ynode_ - 1) {//上壁面
				ncond_[np] = Bcond_.getBCflagR();
			}
			else {//それ以外の領域
				ncond_[np] = 0;
			}
		}
	}
	for (int ie = 0; ie < nelem_; ie++) {
		scond_[ie] = 0;
	}

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
	return ncond_[i];
}
int Mesh2d::scond(int i) {
	return scond_[i];
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
	double x2 = eX(E2);
	double y1 = eY(E1);
	double y2 = eY(E2);

	double x3, x4, y3, y4;//辺の両端の座標
	int n1, n2;//辺の両端の節点番号
	if (e11 == E2) {//E2はE1の下側
		n1 = nbool1(E1, 0);//要素左下の点
		n2 = nbool1(E1, 1);//要素右下の点
	}
	else if (e12 == E2) {//E2はE1の右側
		n1 = nbool1(E1, 1);//要素右下の点
		n2 = nbool1(E1, 2);//要素右上の点
	}
	else if (e13 == E2) {//E2はE1の上側
		n1 = nbool1(E1, 2);//要素右上の点
		n2 = nbool1(E1, 3);//要素左上の点
	}
	else if (e14 == E2) {//E2はE1の左側
		n1 = nbool1(E1, 3);//要素左上の点
		n2 = nbool1(E1, 0);//要素左下の点
	}
	else {
		cout << "E2はE1に隣接してません" << endl;
		exit(-1);
	}
	//辺両端の座標
	x3 = x(n1);
	y3 = y(n1);
	x4 = x(n2);
	y4 = y(n2);
	Scalar2d S;
	//位置ベクトル
	Vector2d V1(x1, y1);
	Vector2d V2(x2, y2);
	Vector2d V3(x3, y3);
	Vector2d V4(x4, y4);

	S = 0.5 * ((V3 - V1) % (V4 - V2));
	return S.v();

}
double Mesh2d::area(int n1, int n2, int n3, int n4) {

	Scalar2d S;
	double x1 = x(n1);
	double x2 = x(n2);
	double x3 = x(n3);
	double x4 = x(n4);
	double y1 = y(n1);
	double y2 = y(n2);
	double y3 = y(n3);
	double y4 = y(n4);
	//位置ベクトル
	Vector2d V1(x1, y1);
	Vector2d V2(x2, y2);
	Vector2d V3(x3, y3);
	Vector2d V4(x4, y4);

	S = 0.5 * ((V3 - V1) % (V4 - V2));
	return S.v();

}
CavityMesh2d::CavityMesh2d(NodeP& NP, Boundarycond& BC)
	:Mesh2d(NP, BC)
{

}
void CavityMesh2d::generate()
{
	//節点座標の計算
	for (int j = 0; j < ynode_; j++) {
		for (int i = 0; i < xnode_; i++) {
			int np = i + xnode_ * j;
			node_[np].setNo(i);//節点番号の設定
			node_[np].setX(dx_ * (double)i + xb_);//節点座標xの設定
			node_[np].setY(dy_ * (double)j + yb_);//節点座標yの設定

			if (i == 0) node_[np].setX(xb_);//x左端の補正
			if (i == xnode_ - 1) node_[np].setX(xt_);//x右端の補正
			if (j == 0) node_[np].setY(yb_);//y下端の補正
			if (j == ynode_ - 1) node_[np].setY(yt_);//y上端の補正
		}
	}
	//nbool1とnbool3の割当
	for (int j = 0; j < yelem_; j++) {
		for (int i = 0; i < xelem_; i++) {
			int ie = i + xelem_ * j;

			int i1 = i + xnode_ * j;//要素の左下点の節点番号
			int i2 = i1 + 1;//要素の右下点の節点番号
			int i4 = i1 + xnode_;//要素の左上点の節点番号
			int i3 = i4 + 1;//要素の右上点の節点番号

			nbool1_[ie][0] = i1;
			nbool1_[ie][1] = i2;
			nbool1_[ie][2] = i3;
			nbool1_[ie][3] = i4;

			nbool3_[ie][2] = ie + xelem_;//上側要素
			if (j == yelem_ - 1) nbool3_[ie][2] = -1;//領域外

			nbool3_[ie][0] = ie - xelem_;//下側要素
			if (j == 0) nbool3_[ie][0] = -1;//領域外

			nbool3_[ie][3] = ie - 1;//左側要素
			if (i == 0) nbool3_[ie][3] = -1;//領域外

			nbool3_[ie][1] = ie + 1;//右側要素
			if (i == xelem_ - 1) nbool3_[ie][1] = -1;//領域外
		}
	}

	//要素座標の設定
	for (int j = 0; j < yelem_; j++) {
		for (int i = 0; i < xelem_; i++) {
			int ie = i + xelem_ * j;
			elem_[ie].setNo(ie);
			int i1 = nbool1_[ie][0];
			int i2 = nbool1_[ie][1];
			int i3 = nbool1_[ie][2];
			int i4 = nbool1_[ie][3];

			elem_[ie].setX(x(i1) + x(i2) + x(i3) + x(i4) / 4);//要素重心のx座標
			elem_[ie].setY(y(i1) + y(i2) + y(i3) + y(i4) / 4);//要素重心のy座標
			double S;//要素面積
			S = area(i1, i2, i3, i4);//算出
			elem_[ie].setSe(S);
		}
	}

	//境界条件フラグの設定
	//境界条件フラグの設定
	//ncond:  0:内部,	1:剛体壁面(流入流出なし),
			//2:流入壁面(dirichlet)　3:流出壁面(neumann)
			//4:移動壁面条件(壁面接線方向に流速固定値1),5:滑りなし壁面条件(壁面において(u,v)=(0,0))
			//6:滑りあり壁面条件(dvx/dy=0,vy=0)

	for (int j = 0; j < ynode_; j++) {
		for (int i = 0; i < xnode_; i++) {
			int np = i + xnode_ * j;
			ncond_[np] = 0;

		}
	}


	for (int j = 0; j < ynode_; j++) {
		for (int i = 0; i < xnode_; i++) {
			int np = i + xnode_ * j;
			if (i == 0) {//左壁面
				if (j != 0 || j != ynode_ - 1) {//角をのぞく
					ncond_[np] = Bcond_.getBCflagL();
				}
			}
			if (j == 0) {//下壁面
				ncond_[np] = Bcond_.getBCflagD();
				if (i == 0) {//左下角
					ncond_[np] = Bcond_.getBCflagL();
				}
				if (i == xnode_ - 1) {//右下角
					ncond_[np] = Bcond_.getBCflagR();
				}

			}
			if (i == xnode_ - 1) {//右壁面
				if (j != 0 || j != ynode_ - 1) {//角をのぞく
					ncond_[np] = Bcond_.getBCflagR();
				}
			}
			if (j == ynode_ - 1) {//上壁面
				ncond_[np] = Bcond_.getBCflagU();
	
			}

		}
		
	}
	//物体要素フラグの設定
	
	for (int ie = 0; ie < nelem_; ie++) {
		scond_[ie] = 0;
	}

}
BackstepMesh2d::BackstepMesh2d(NodeP& NP, Boundarycond& BC)
	:Mesh2d(NP, BC), hx(0), hy(0)
{
	cout << "バックステップ流れMesh生成" << endl;
	cout << "段差のステップ長さ hx->";
	cin >> hx;
	cout << "段差のステップ高さ hy->";
	cin >> hy;
}
void BackstepMesh2d::generate() {
	//節点座標の計算
	for (int j = 0; j < ynode_; j++) {
		for (int i = 0; i < xnode_; i++) {
			int np = i + xnode_ * j;
			node_[np].setNo(i);//節点番号の設定
			node_[np].setX(dx_ * (double)i + xb_);//節点座標xの設定
			node_[np].setY(dy_ * (double)j + yb_);//節点座標yの設定

			if (i == 0) node_[np].setX(xb_);//x左端の補正
			if (i == xnode_ - 1) node_[np].setX(xt_);//x右端の補正
			if (j == 0) node_[np].setY(yb_);//y下端の補正
			if (j == ynode_ - 1) node_[np].setY(yt_);//y上端の補正
		}
	}
	//nbool1とnbool3の割当
	for (int j = 0; j < yelem_; j++) {
		for (int i = 0; i < xelem_; i++) {
			int ie = i + xelem_ * j;

			int i1 = i + xnode_ * j;//要素の左下点の節点番号
			int i2 = i1 + 1;//要素の右下点の節点番号
			int i4 = i1 + xnode_;//要素の左上点の節点番号
			int i3 = i4 + 1;//要素の右上点の節点番号

			nbool1_[ie][0] = i1;
			nbool1_[ie][1] = i2;
			nbool1_[ie][2] = i3;
			nbool1_[ie][3] = i4;

			nbool3_[ie][2] = ie + xelem_;//上側要素
			if (j == yelem_ - 1) nbool3_[ie][2] = -1;//領域上端

			nbool3_[ie][0] = ie - xelem_;//下側要素
			if (j == 0) nbool3_[ie][0] = -1;//領域下端

			nbool3_[ie][3] = ie - 1;//左側要素
			if (i == 0) nbool3_[ie][3] = -1;//領域左端

			nbool3_[ie][1] = ie + 1;//右側要素
			if (i == xelem_ - 1) nbool3_[ie][1] = -1;//領域右端
		}
	}

	//要素座標の設定
	for (int j = 0; j < yelem_; j++) {
		for (int i = 0; i < xelem_; i++) {
			int ie = i + xelem_ * j;
			elem_[ie].setNo(ie);
			int i1 = nbool1_[ie][0];
			int i2 = nbool1_[ie][1];
			int i3 = nbool1_[ie][2];
			int i4 = nbool1_[ie][3];

			elem_[ie].setX(x(i1) + x(i2) + x(i3) + x(i4) / 4);//要素重心のx座標
			elem_[ie].setY(y(i1) + y(i2) + y(i3) + y(i4) / 4);//要素重心のy座標
			double S;//要素面積
			S = area(i1, i2, i3, i4);//算出
			elem_[ie].setSe(S);

		}
	}

	//境界条件フラグの設定
	//境界条件フラグの設定
	//ncond:  0:内部,	1:剛体壁面(流入流出なし),
			//2:流入壁面(dirichlet)　3:流出壁面(neumann)
			//4:移動壁面条件(壁面接線方向に流速固定値1),5:滑りなし壁面条件(壁面において(u,v)=(0,0))
			//6:滑りあり壁面条件(dvx/dy=0,vy=0)

	for (int j = 0; j < ynode_; j++) {
		for (int i = 0; i < xnode_; i++) {
			int np = i + xnode_ * j;
			ncond_[np] = 0;

		}
	}


	for (int j = 0; j < ynode_; j++) {
		for (int i = 0; i < xnode_; i++) {
			int np = i + xnode_ * j;
			if (i == 0) {//左壁面
				if (j != 0 || j != ynode_ - 1) {//角をのぞく
					ncond_[np] = Bcond_.getBCflagL();
				}
			}
			if (j == 0) {//下壁面
				ncond_[np] = Bcond_.getBCflagD();
				if (i == 0) {//左下角
					ncond_[np] = Bcond_.getBCflagL();
				}
				if (i == xnode_ - 1) {//右下角
					ncond_[np] = Bcond_.getBCflagR();
				}

			}
			if (i == xnode_ - 1) {//右壁面
				if (j != 0 || j != ynode_ - 1) {//角をのぞく
					ncond_[np] = Bcond_.getBCflagR();
				}
			}
			if (j == ynode_ - 1) {//上壁面
				ncond_[np] = Bcond_.getBCflagU();
				if (i == 0) {//左上角
					ncond_[np] = Bcond_.getBCflagL();
				}
				if (i == xnode_ - 1) {//右上角
					ncond_[np] = Bcond_.getBCflagR();
				}
			}
			
		}
	}
	for (int j = 0; j < ynode_; j++) {
		for (int i = 0; i < xnode_; i++) {
			int np = i + xnode_ * j;
			if (x(np) <= (hx + xb_)&& y(np) <= (hy + yb_)) {
				//段差の部分を探る
				ncond_[np] = 1;//剛体内部
				if (x(np + 1) > (hx + xb_)) {
					ncond_[np] = Bcond_.getBCflagC();//物体壁面
				}
				if (y(np + xnode_) > (hy + yb_)) {
					ncond_[np] = Bcond_.getBCflagC();//物体壁面
				}
			}
		}
	}
	for (int ie = 0; ie < nelem_; ie++) {
		scond_[ie] = 0;
	}
	for (int ie = 0; ie < nelem_; ie++) {
		int i1 = nbool1_[ie][0];
		int i2 = nbool1_[ie][1];
		int i3 = nbool1_[ie][2];
		int i4 = nbool1_[ie][3];

		if (ncond(i1) == 1 || ncond(i2) == 1 || ncond(i3) == 1 || ncond(i4) == 1) {
			//要素内の点に剛体内部のフラグ(1)が立つ
			scond_[ie] = 1;
		}
		if (ncond(i1) != 0 && ncond(i2) != 0 && ncond(i3) != 0 && ncond(i4) != 0) {
			//要素内のすべての点が壁面のフラグ
			scond_[ie] = 1;
		}

	}
}