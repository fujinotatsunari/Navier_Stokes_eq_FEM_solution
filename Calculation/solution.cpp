#include"solution.hpp"
#include"output.hpp"
#include"FEM.hpp"
#include"value.hpp"
#include"Mesh.hpp"
#include"param.hpp"
#include"matrix.hpp"
#include<vector>
#include<iostream>

using namespace std;

Explicit_FEM::Explicit_FEM(Mesh2d& mesh_, Time& t_, PHI& phi_, Boundarycond& BC_, ADeq_param_2d& adp_)
	:mesh(mesh_), t(t_), phi(phi_), BC(BC_), ADP(adp_)
{

}

void Explicit_FEM::do_expcalculation() {

	Lumped_Massmatrix Fm(mesh);//集中化質量行列
	xAdvecmatrix Amx(mesh);//移流行列x方向
	yAdvecmatrix Amy(mesh);//移流行列y方向
	Diffmatrix Dm(mesh);//拡散行列
	OutputData out(mesh, t, phi, ADP, BC);
	out.set_scheme(0);//陽解法をset
	out.set_Filestage(0);
	out.directory_setup();
	//cout << out.get_dir() << endl;
	//cout << out.get_Filestage() << endl;
	//cout << out.get_dir() << endl;
	/*
	cout << "viewFm" << endl;
	Fm.view();
	cout << "viewAmx" << endl;
	Amx.view();
	cout << "viewAmy" << endl;
	Amy.view();
	cout << "viewDm" << endl;
	Dm.view();
	*/
	vector<double> phib;//1step前の物理量
	phib.resize(mesh.nnode());
	vector<double> dd;//拡散項足し込み変数
	dd.resize(mesh.nnode());
	vector<double> uu;//移流項足し込み変数
	uu.resize(mesh.nnode());
	vector<double> nn;//境界項足し込み変数
	nn.resize(mesh.nnode());
	vector<double> ff;//集中化質量行列の節点への寄与
	ff.resize(mesh.nnode());

	vector<int> nx;//境界上の単位法線ベクトルx方向
	nx.resize(mesh.nnode());
	vector<int> ny;//境界上の単位法線ベクトルx方向
	ny.resize(mesh.nnode());
	double dphidx = 0.0;
	double dphidy = 0.0;
	for (int j = 0; j < mesh.ynode(); j++) {
		for (int i = 0; i < mesh.xnode(); i++) {
			int np = i + mesh.xnode() * j;
			nx[np] = 0;
			ny[np] = 0;
			if (i == 0) {
				nx[np] = -1;//境界法線単位ベクトル
			}
			if (i == mesh.xnode() - 1) {
				nx[np] = 1;
			}
			if (j == 0) {
				ny[np] = -1;
			}
			if (j == mesh.ynode() - 1) {
				ny[np] = 1;
			}
		}
	}

	//時間進行
	for (int n = 0; n <= t.nend(); n++) {


		cout << "time=" << t.ntime(n) << endl;
		for (int j = 0; j < mesh.ynode(); j++) {
			for (int i = 0; i < mesh.xnode(); i++) {
				int np = i + mesh.xnode() * j;

				phib[np] = phi[np];
				dd[np] = 0.0;
				uu[np] = 0.0;
				ff[np] = 0.0;
				nn[np] = 0.0;
			}
		}
		for (int j = 0; j < node; j++) {//局所節点ループ
			for (int ie = 0; ie < mesh.nelem(); ie++) {//要素ループ
				int np = mesh.nbool1(ie, j);
				int i1 = mesh.i1(ie);
				int i2 = mesh.i2(ie);
				int i3 = mesh.i3(ie);
				int i4 = mesh.i4(ie);

				dd[np] = dd[np] - ADP.get_alpha() * (Dm[ie][j][0] * phi[i1] + Dm[ie][j][1] * phi[i2] + Dm[ie][j][2] * phi[i3] + Dm[ie][j][3] * phi[i4]);//拡散項
				uu[np] = uu[np] - ADP.get_cx() * (Amx[ie][j][0] * phi[i1] + Amx[ie][j][1] * phi[i2] + Amx[ie][j][2] * phi[i3] + Amx[ie][j][3] * phi[i4])//移流項
					- ADP.get_cy() * (Amy[ie][j][0] * phi[i1] + Amy[ie][j][1] * phi[i2] + Amy[ie][j][2] * phi[i3] + Amy[ie][j][3] * phi[i4]);
				ff[np] = ff[np] + Fm[ie][j][j];

				if (mesh.ncond(np) == 2) {
					nn[np] = nn[np] + ADP.get_alpha() * (nx[np] * dphidx + ny[np] * dphidy);//今回は流出のみをあつかうので境界項＝0
				}


			}
		}


		for (int j = 0; j < mesh.ynode(); j++) {
			for (int i = 0; i < mesh.xnode(); i++) {
				int np = i + mesh.xnode() * j;
				if (mesh.ncond(np) == 0) {//内部
					phi[np] = phib[np] + t.dt() * (dd[np] + uu[np]) / ff[np];


				}
				else if (mesh.ncond(np) == 1) {//dirichlet境界条件
					phi[np] = phib[np];
				}
				else if (mesh.ncond(np) == 2) {//neumann境界条件
					phi[np] = phib[np] + t.dt() * (dd[np] + uu[np] + nn[np]) / ff[np];

				}
				else {

				}
				//cout << "phi[" << np << "]=" << phi[np] << endl;
			}
		}
		if (n == 0) {
			out.output_condition();
		}
		if (n % t.nsample() == 0) {
			out.data_update(phi);
			out.output_result_csv(n);
		}
	}
}


Implicit_FEM::Implicit_FEM(Mesh2d& mesh_, Time& t_, PHI& phi_, Boundarycond& BC_, ADeq_param_2d& adp_)
	:mesh(mesh_), t(t_), phi(phi_), BC(BC_), ADP(adp_)
{

}

void Implicit_FEM::do_impcaluculation() {
	Lumped_Massmatrix Em(mesh);//質量行列
	xAdvecmatrix Amx(mesh);//移流行列x方向
	yAdvecmatrix Amy(mesh);//移流行列y方向
	Diffmatrix Dm(mesh);//拡散行列
	int Nnode = mesh.nnode();
	OutputData out(mesh, t, phi, ADP, BC);
	out.set_scheme(1);//陰解法をset
	out.set_Filestage(0);
	out.directory_setup();


	vector<vector<double>> A;
	A.resize(Nnode);
	for (int i = 0; i < A.size(); i++) {
		A[i].resize(Nnode);
	}
	vector<double> b;
	b.resize(mesh.nnode());
	vector<double> x;
	x.resize(mesh.nnode());

	vector<double> nn;//境界項足し込み変数
	nn.resize(mesh.nnode());
	double dphidx = 0.0;
	double dphidy = 0.0;
	vector<int> nx;//境界上の単位法線ベクトルx方向
	nx.resize(mesh.nnode());
	vector<int> ny;//境界上の単位法線ベクトルx方向
	ny.resize(mesh.nnode());

	for (int j = 0; j < mesh.ynode(); j++) {
		for (int i = 0; i < mesh.xnode(); i++) {
			int np = i + mesh.xnode() * j;
			nx[np] = 0;
			ny[np] = 0;
			if (i == 0) {
				nx[np] = -1;//境界法線単位ベクトル
			}
			if (i == mesh.xnode() - 1) {
				nx[np] = 1;
			}
			if (j == 0) {
				ny[np] = -1;
			}
			if (j == mesh.ynode() - 1) {
				ny[np] = 1;
			}
		}
	}

	cout << "time=" << t.ntime(0) << endl;
	out.output_condition();
	out.output_result_csv(0);
	//時間進行
	for (int n = 1; n <= t.nend(); n++) {
		cout << "time=" << t.ntime(n) << endl;
		for (int i = 0; i < A.size(); i++) {
			for (int j = 0; j < A[i].size(); j++) {
				A[i][j] = 0.0;
			}
		}
		for (int i = 0; i < mesh.nnode(); i++) {
			b[i] = 0.0;
			x[i] = 0.0;
			nn[i] = 0.0;

		}
		for (int j = 0; j < node; j++) {
			for (int ie = 0; ie < mesh.nelem(); ie++) {
				int np = mesh.nbool1(ie, j);
				int i1 = mesh.i1(ie);
				int i2 = mesh.i2(ie);
				int i3 = mesh.i3(ie);
				int i4 = mesh.i4(ie);

				if (mesh.ncond(np) == 1) {//dirichlet
					A[np][np] = 1.0;
					b[np] = phi[np];

				}
				else {//内部
					A[np][i1] = A[np][i1] + Em[ie][j][0] + t.dt() * (ADP.get_cx() * Amx[ie][j][0] + ADP.get_cy() * Amy[ie][j][0]) + t.dt() * ADP.get_alpha() * Dm[ie][j][0];
					A[np][i2] = A[np][i2] + Em[ie][j][1] + t.dt() * (ADP.get_cx() * Amx[ie][j][1] + ADP.get_cy() * Amy[ie][j][1]) + t.dt() * ADP.get_alpha() * Dm[ie][j][1];
					A[np][i3] = A[np][i3] + Em[ie][j][2] + t.dt() * (ADP.get_cx() * Amx[ie][j][2] + ADP.get_cy() * Amy[ie][j][2]) + t.dt() * ADP.get_alpha() * Dm[ie][j][2];
					A[np][i4] = A[np][i4] + Em[ie][j][3] + t.dt() * (ADP.get_cx() * Amx[ie][j][3] + ADP.get_cy() * Amy[ie][j][3]) + t.dt() * ADP.get_alpha() * Dm[ie][j][3];
					b[np] = b[np] + Em[ie][j][0] * phi[i1] + Em[ie][j][1] * phi[i2] + Em[ie][j][2] * phi[i3] + Em[ie][j][3] * phi[i4];
				}
				if (mesh.ncond(np) == 2) {//neumann
					nn[np] = nn[np] + ADP.get_alpha() * (nx[np] * dphidx + ny[np] * dphidy);//今回は流出のみをあつかうので境界項＝0
					b[np] = b[np] + nn[np];
				}
			}
		}
		LU_solve(A, x, b);
		for (int i = 0; i < x.size(); i++) {
			phi[i] = x[i];
		}
		if (n % t.nsample() == 0) {
			out.data_update(phi);
			out.output_result_csv(n);
		}

	}

}
