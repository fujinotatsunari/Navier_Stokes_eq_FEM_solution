#include<iostream>
#include"param.hpp"
#include"mesh.hpp"
#include"value.hpp"
#include"output.hpp"

using namespace std;

int main(void) {
	cout << "2次元 Navier-Stokes方程式 Cavity流れシミュレーション Mesh作成プログラム" << endl;
	cout << "解析modelを選択してください" << endl;
	int model;
	cout << "0:cavity流れ , 1:back step流れ, 2:角柱周り流れ " << endl;
	cin >> model;

	NodeP np;
	Boundarycond BC;

	if (model == 0) {//cavity流れ
		BC.set_cavityBC();
		CavityMesh2d mesh(np, BC);
		mesh.generate();
		Velocity2d V(mesh, BC);
		Pressure P(mesh, BC);
		V.cavity_init();
		P.cavity_init();

		OutputData output(mesh, BC, V, P);
		output.set_modelname("cavity");
		output.output_csv();
		output.output_dat();
		output.output_condition();
		cout << "Meshファイル作成完了" << endl;
	}
	if (model == 1) {//back step流れ
		BC.set_backstepBC();
		BackstepMesh2d mesh(np, BC);
		mesh.generate();
		Velocity2d V(mesh, BC);
		Pressure P(mesh, BC);
		V.backstep_init();
		P.backstep_init();

		OutputData output(mesh, BC, V, P);
		output.set_modelname("backstep");
		output.output_csv();
		output.output_dat();
		output.output_condition();
		cout << "Meshファイル作成完了" << endl;
	}
	if (model == 2) {//角柱周り流れ
		BC.set_pillarBC();
		SquarePillarMesh2d mesh(np, BC);
		mesh.generate();
		Velocity2d V(mesh, BC);
		Pressure P(mesh, BC);
		V.Pillar_init();
		P.Pillar_init();

		OutputData output(mesh, BC, V, P);
		output.set_modelname("pillar");
		output.output_csv();
		output.output_dat();
		output.output_condition();
		cout << "Meshファイル作成完了" << endl;
	}
	
	return 0;

}

