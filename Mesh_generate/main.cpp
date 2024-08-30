#include<iostream>
#include"param.hpp"
#include"mesh.hpp"
#include"value.hpp"
#include"output.hpp"

using namespace std;

int main(void) {
	cout << "2次元 Navier-Stokes方程式 Cavity流れシミュレーション Mesh作成プログラム" << endl;
	NodeP np;
	Boundarycond BC;
	BC.set_cavityBC();
	CavityMesh2d mesh(np, BC);
	mesh.generate();

	Velocity2d V(mesh, BC);
	Pressure P(mesh, BC);
	V.cavity_init();
	P.cavity_init();
	Outputcavitymesh output(mesh, BC, V, P);
	output.set_Filestage(0);
	output.output_csv();
	output.output_dat();
	output.output_condition();
	cout << "Meshファイル作成完了" << endl;
	return 0;

}

