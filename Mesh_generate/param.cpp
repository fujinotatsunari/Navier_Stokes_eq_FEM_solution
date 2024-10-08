#include"param.hpp"
#include<cstdio>
#include<cmath>
#include<vector>
#include<iostream>
using namespace std;

NodeP::NodeP()
	:Lx(0), Ly(0), xb(0), xt(0), yb(0), yt(0), dx(0), dy(0), xnode(0), ynode(0), xelem(0), yelem(0), nnode(0), nelem(0)
{
	
	setNodeParam();
}
NodeP::NodeP(const NodeP& Np) {
	Lx = Np.Lx;
	Ly = Np.Ly;
	xb = Np.xb;
	xt = Np.xt;
	yb = Np.yb;
	yt = Np.yt;
	dx = Np.dx;
	dy = Np.dy;
	xnode = Np.xnode;
	ynode = Np.ynode;
	xelem = Np.xelem;
	yelem = Np.yelem;
	nnode = Np.nnode;
	nelem = Np.nelem;
}
NodeP& NodeP::operator=(const NodeP& Np) {
	Lx = Np.Lx;
	Ly = Np.Ly;
	xb = Np.xb;
	xt = Np.xt;
	yb = Np.yb;
	yt = Np.yt;
	dx = Np.dx;
	dy = Np.dy;
	xnode = Np.xnode;
	ynode = Np.ynode;
	xelem = Np.xelem;
	yelem = Np.yelem;
	nnode = Np.nnode;
	nelem = Np.nelem;

	return *this;
}
void NodeP::setNodeParam() {
	
	cout << "X軸左端: xb->";
	cin >> xb;
	cout << "X軸右端: xt->";
	cin >> xt;
	cout << "Y軸下端: yb->";
	cin >> yb;
	cout << "Y軸上端: yt->";
	cin >> yt;
	Lx = xt - xb;
	Ly = yt - yb;

	cout << "x方向要素数 ->";
	cin >> xelem;
	cout << "y方向要素数 ->";
	cin >> yelem;
	xnode = xelem + 1;
	ynode = yelem + 1;
	nnode = xnode * ynode;
	nelem = xelem * yelem;
	cout << "総節点数 -> 　" << nnode << ",総要素数 ->" << nelem << endl;
	dx = Lx / (double)xelem;
	dy = Ly / (double)yelem;
	cout << "dx=" << dx << endl;
	cout << "dy=" << dy << endl;
}
double NodeP::getLx() {
	return Lx;
}
double NodeP::getLy() {
	return Ly;
}
double NodeP::getXb() {
	return xb;
}
double NodeP::getXt() {
	return xt;
}
double NodeP::getYb() {
	return yb;
}
double NodeP::getYt() {
	return yt;
}
double NodeP::getDx() {
	return dx;
}
double NodeP::getDy() {
	return dy;
}
int NodeP::getXnode() {
	return xnode;
}
int NodeP::getYnode() {
	return ynode;
}
int NodeP::getXelem() {
	return xelem;
}
int NodeP::getYelem() {
	return yelem;
}
int NodeP::getNnode() {
	return nnode;
}
int NodeP::getNelem() {
	return nelem;
}


Boundarycond::Boundarycond()
	:flagL(0), flagR(0), flagU(0), flagD(0),flagC(0)
{
	
}
void Boundarycond::set_userBC() {
	cout << "境界条件の設定" << endl;
	cout << "2:流入境界条件(diriclet)  3:流出境界条件(neumann)" << endl;
	cout << "4:移動壁面条件  5:滑りなし壁面条件　6:滑りあり壁面条件" << endl;
	
	cout << "上壁面の境界条件を決めてください" << endl;
	cin >> flagU;
	cout << "左壁面の境界条件を決めてください" << endl;
	cin >> flagL;
	cout << "右壁面の境界条件を決めてください" << endl;
	cin >> flagR;
	cout << "下壁面の境界条件を決めてください" << endl;
	cin >> flagD;

	
}
void Boundarycond::set_cavityBC() {

	cout << "cavity流れの境界条件の設定" << endl;
	cout << "2:流入境界条件(diriclet)  3:流出境界条件(neumann)" << endl;
	cout << "4:移動壁面条件  5:滑りなし壁面条件　6:滑りあり壁面条件" << endl;

	cout << "上壁面の境界条件->" << 4 << endl;
	flagU = 4;
	cout << "左壁面の境界条件->" << 5 << endl;
	flagL = 5;
	cout << "右壁面の境界条件->" << 5 << endl;
	flagR = 5;
	cout << "下壁面の境界条件->" << 5 << endl;
	flagD = 5;
	flagC = -1;
}
void Boundarycond::set_cylinderBC() {
	cout << "円柱周り流れの境界条件の設定" << endl;

	cout << "2:流入境界条件(diriclet)  3:流出境界条件(neumann)" << endl;
	cout << "4:移動壁面条件  5:滑りなし壁面条件　6:滑りあり壁面条件" << endl;

	cout << "上壁面の境界条件->(5or6)";
	cin >> flagU;
	cout << "左壁面の境界条件->" << 2 << endl;
	flagL = 2;
	cout << "右壁面の境界条件->" << 3 << endl;
	flagR = 3;
	cout << "下壁面の境界条件->(5or6)";
	cin >> flagD;
	cout << "円柱壁面の境界条件->(5or6)";
	cin >> flagC;
}
void Boundarycond::set_pillarBC() {
	cout << "角柱周り流れの境界条件の設定" << endl;

	cout << "2:流入境界条件(diriclet)  3:流出境界条件(neumann)" << endl;
	cout << "4:移動壁面条件  5:滑りなし壁面条件　6:滑りあり壁面条件" << endl;

	cout << "上壁面の境界条件->(5or6)";
	cin >> flagU;
	cout << "左壁面の境界条件->" << 2 << endl;
	flagL = 2;
	cout << "右壁面の境界条件->" << 3 << endl;
	flagR = 3;
	cout << "下壁面の境界条件->(5or6)";
	cin >> flagD;
	cout << "円柱壁面の境界条件->(5or6)";
	cin >> flagC;
}
void Boundarycond::set_backstepBC() {
	cout << "バックステップ流れの境界条件の設定" << endl;

	cout << "2:流入境界条件(diriclet)  3:流出境界条件(neumann)" << endl;
	cout << "4:移動壁面条件  5:滑りなし壁面条件　6:滑りあり壁面条件" << endl;

	cout << "上壁面の境界条件->(5or6)";
	cin >> flagU;
	cout << "左壁面の境界条件->" << 2 << endl;
	flagL = 2;
	cout << "右壁面の境界条件->" << 3 << endl;
	flagR = 3;
	cout << "下壁面の境界条件->(5or6)";
	cin >> flagD;
	cout << "角柱壁面の境界条件->(5or6)";
	cin >> flagC;
}
int Boundarycond::getBCflagL() {
	return flagL;
}
int Boundarycond::getBCflagR() {
	return flagR;
}
int Boundarycond::getBCflagU() {
	return flagU;
}
int Boundarycond::getBCflagD() {
	return flagD;
}
int Boundarycond::getBCflagC() {
	return flagC;
}
void Boundarycond::setBCflagL(int flag) {
	flagL = flag;
}
void Boundarycond::setBCflagR(int flag) {
	flagR = flag;
}
void Boundarycond::setBCflagD(int flag) {
	flagD = flag;
}
void Boundarycond::setBCflagU(int flag) {
	flagU = flag;
}
void Boundarycond::setBCflagC(int flag) {
	flagC = flag;
}
/*
ADeq_param_2d::ADeq_param_2d(NodeP& Np, TimeP& Tp)
	:nparam(Np),tparam(Tp)
{


	set_param();
}
void ADeq_param_2d::set_param() {
	cout << "x方向定常流速cx->";
	cin >> cx;
	cout << "y方向定常流速cy->";
	cin >> cy;
	cout << "拡散係数　alpha->";
	cin >> alpha;
	courantx = cx * tparam.getDt() / nparam.getDx();
	couranty = cy * tparam.getDt() / nparam.getDy();
	diffusion = alpha * tparam.getDt() / (nparam.getDx() * nparam.getDy());
	cout << "x方向 courant数 Cx=" << courantx << endl;
	cout << "y方向 courant数 Cy=" << couranty << endl;
	cout << "拡散数 D=" << diffusion << endl;
	Pe = sqrt(cx * cx + cy * cy) * sqrt(nparam.getLx() * nparam.getLx() + nparam.getLy() * nparam.getLy()) / alpha;
	cout << "Peclet数 Pe=" << Pe << endl;

}
double ADeq_param_2d::get_alpha() {
	return alpha;
}
double ADeq_param_2d::get_cx() {
	return cx;
}
double ADeq_param_2d::get_cy() {
	return cy;
}
double ADeq_param_2d::get_couranx() {
	return courantx;
}
double ADeq_param_2d::get_courany() {
	return couranty;
}
double ADeq_param_2d::get_diffusion() {
	return diffusion;
}
double ADeq_param_2d::get_Pe() {
	return Pe;
}*/