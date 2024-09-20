#include"FEM.hpp"
#include"value.hpp"
#include"Mesh.hpp"
#include"matrix.hpp"
#include<vector>


CofficientMatrix::CofficientMatrix(Mesh2d& Mesh)
	:mesh(Mesh)
{
	mat.resize(mesh.nelem());
	for (int i = 0; i < mat.size(); i++) {
		mat[i] = Matrix::new_zero_matrix(node, node);
	}
}


const Matrix& CofficientMatrix::operator[](int ie)const {
	return mat.at(ie);
}
Matrix& CofficientMatrix::operator[](int ie) {
	return mat.at(ie);
}

void CofficientMatrix::view() {
	for (int ie = 0; ie < mat.size(); ie++) {
		cout << "nelem=" << ie << endl;
		mat[ie].print();
		cout << endl;
	}
}

Massmatrix::Massmatrix(Mesh2d& Mesh)
	:CofficientMatrix(Mesh)
{
	//質量行列を求める
	for (int ie = 0; ie < mesh.nelem(); ie++) {
		for (int i = 0; i < node; i++) {
			for (int j = 0; j < node; j++) {
				for (int k = 0; k < node; k++) {
					for (int l = 0; l < node; l++) {
						mat[ie][i][j] += (3.0 + xi[i] * xi[j]) * (3.0 + eta[i] * eta[j])
							* ((xi[k] * eta[l] - eta[k] * xi[l]) * mesh.x(mesh.nbool1(ie, k)) * mesh.y(mesh.nbool1(ie, l))
								+ (xi[i] + xi[j]) * (xi[k] * xi[l] * (eta[l] - eta[k]) * mesh.x(mesh.nbool1(ie, k)) * mesh.y(mesh.nbool1(ie, l))) / (3.0 + xi[i] * xi[j])
								+ (eta[i] + eta[j]) * (eta[k] * eta[l] * (xi[k] - xi[l]) * mesh.x(mesh.nbool1(ie, k)) * mesh.y(mesh.nbool1(ie, l))) / (3.0 + eta[i] * eta[j])) / 576.0;

					}
				}
			}
		}
	}


}



Lumped_Massmatrix::Lumped_Massmatrix(Mesh2d& Mesh)
	:CofficientMatrix(Mesh)
{
	//集中化質量行列を求める
	double EM1;
	double EM2;
	double EM3;
	for (int ie = 0; ie < mesh.nelem(); ie++) {
		for (int i = 0; i < node; i++) {
			EM1 = 0.0;
			EM2 = 0.0;
			EM3 = 0.0;
			for (int k = 0; k < node; k++) {
				for (int l = 0; l < node; l++) {

					EM1 = EM1 + xi[k] * eta[l] * (mesh.x(mesh.nbool1(ie, k)) * mesh.y(mesh.nbool1(ie, l)) - mesh.y(mesh.nbool1(ie, k)) * mesh.x(mesh.nbool1(ie, l)));

					EM2 = EM2 + xi[k] * xi[l] * eta[l] * (mesh.x(mesh.nbool1(ie, k)) * mesh.y(mesh.nbool1(ie, l)) - mesh.y(mesh.nbool1(ie, k)) * mesh.x(mesh.nbool1(ie, l)));

					EM3 = EM3 + xi[k] * eta[k] * eta[l] * (mesh.x(mesh.nbool1(ie, k)) * mesh.y(mesh.nbool1(ie, l)) - mesh.y(mesh.nbool1(ie, k)) * mesh.x(mesh.nbool1(ie, l)));
				}
			}
			mat[ie][i][i] = EM1 / 16.0 + xi[i] * EM2 / 48.0 + eta[i] * EM3 / 48.0;

		}
	}
}


xAdvecmatrix::xAdvecmatrix(Mesh2d& Mesh)
	:CofficientMatrix(Mesh)
{
	//x方向移流行列行列を求める
	for (int ie = 0; ie < mesh.nelem(); ie++) {
		for (int i = 0; i < node; i++) {
			for (int j = 0; j < node; j++) {
				for (int k = 0; k < node; k++) {
					mat[ie][i][j] += xi[k] * eta[k] * (xi[i] * xi[j] - eta[i] * eta[j]) * mesh.y(mesh.nbool1(ie, k)) / 48.0
						+ xi[j] * eta[k] * (1 + eta[i] * eta[j] / 3.0) * mesh.y(mesh.nbool1(ie, k)) / 16.0
						- eta[j] * xi[k] * (1 + xi[i] * xi[j] / 3.0) * mesh.y(mesh.nbool1(ie, k)) / 16.0;

				}
			}
		}
	}
}

yAdvecmatrix::yAdvecmatrix(Mesh2d& Mesh)
	:CofficientMatrix(Mesh)
{
	//y方向移流行列行列を求める
	for (int ie = 0; ie < mesh.nelem(); ie++) {
		for (int i = 0; i < node; i++) {
			for (int j = 0; j < node; j++) {
				for (int k = 0; k < node; k++) {
					mat[ie][i][j] += -xi[k] * eta[k] * (xi[i] * xi[j] - eta[i] * eta[j]) * mesh.x(mesh.nbool1(ie, k)) / 48.0
						- xi[j] * eta[k] * (1 + eta[i] * eta[j] / 3.0) * mesh.x(mesh.nbool1(ie, k)) / 16.0
						+ eta[j] * xi[k] * (1 + xi[i] * xi[j] / 3.0) * mesh.x(mesh.nbool1(ie, k)) / 16.0;

				}
			}
		}
	}
}


Advecmatrix::Advecmatrix(Mesh2d& Mesh, Velocity2d& v) 
	:CofficientMatrix(Mesh),V(v)
{
	//3点数値積分により移流行列を求める
	vector<double> xi_gauss = { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };//数値積分離散点xi方向
	vector<double> eta_gauss = { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };//数値積分離散点eta方向
	vector<double> gauss_W = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };//数値積分重み
	double n_alpha, n_beta, n_gamma1, n_gamma2, n_gamma3, n_gamma4, dndxi_beta, dndeta_beta;//形状関数とその微分
	double A11, A12, A21, A22;//余因子要素
	double tax1, tax2, tax3, tax4, tay1, tay2, tay3, tay4;
	for (int ie = 0; ie < mesh.nelem(); ie++) {
		int i1 = mesh.i1(ie);
		int i2 = mesh.i2(ie);
		int i3 = mesh.i3(ie);
		int i4 = mesh.i4(ie);
		for (int i = 0; i < node; i++) {
			for (int j = 0; j < node; j++) {
				for (int gauss_i = 0; gauss_i < 3; gauss_i++) {
					for (int gauss_j = 0; gauss_j < 3; gauss_j++) {
						n_alpha = 0.25 * (1.0 + xi[i] * xi_gauss[gauss_i]) * (1.0 + eta[i] * eta_gauss[gauss_j]);
						n_beta = 0.25 * (1.0 + xi[j] * xi_gauss[gauss_i]) * (1.0 + eta[j] * eta_gauss[gauss_j]);

						n_gamma1 = 0.25 * (1.0 + xi[0] * xi_gauss[gauss_i]) * (1.0 + eta[0] * eta_gauss[gauss_j]);
						n_gamma2 = 0.25 * (1.0 + xi[1] * xi_gauss[gauss_i]) * (1.0 + eta[1] * eta_gauss[gauss_j]);
						n_gamma3 = 0.25 * (1.0 + xi[2] * xi_gauss[gauss_i]) * (1.0 + eta[2] * eta_gauss[gauss_j]);
						n_gamma4 = 0.25 * (1.0 + xi[3] * xi_gauss[gauss_i]) * (1.0 + eta[3] * eta_gauss[gauss_j]);

						dndxi_beta = 0.25 * xi[j] * (1.0 + eta[j] * eta_gauss[gauss_j]);
						dndeta_beta = 0.25 * eta[j] * (1.0 + xi[j] * xi_gauss[gauss_i]);


						A11 = 0.25 * eta[0] * (1 + xi[0] * xi_gauss[gauss_i]) * mesh.y(i1)
							+ 0.25 * eta[1] * (1 + xi[1] * xi_gauss[gauss_i]) * mesh.y(i2)
							+ 0.25 * eta[2] * (1 + xi[2] * xi_gauss[gauss_i]) * mesh.y(i3)
							+ 0.25 * eta[3] * (1 + xi[3] * xi_gauss[gauss_i]) * mesh.y(i4);


						A12 = -0.25 * xi[0] * (1 + eta[0] * eta_gauss[gauss_j]) * mesh.y(i1)
							- 0.25 * xi[1] * (1 + eta[1] * eta_gauss[gauss_j]) * mesh.y(i2)
							- 0.25 * xi[2] * (1 + eta[2] * eta_gauss[gauss_j]) * mesh.y(i3)
							- 0.25 * xi[3] * (1 + eta[3] * eta_gauss[gauss_j]) * mesh.y(i4);

						A21 = -0.25 * eta[0] * (1 + xi[0] * xi_gauss[gauss_i]) * mesh.x(i1)
							- 0.25 * eta[1] * (1 + xi[1] * xi_gauss[gauss_i]) * mesh.x(i2)
							- 0.25 * eta[2] * (1 + xi[2] * xi_gauss[gauss_i]) * mesh.x(i3)
							- 0.25 * eta[3] * (1 + xi[3] * xi_gauss[gauss_i]) * mesh.x(i4);

						A22 = 0.25 * xi[0] * (1 + eta[0] * eta_gauss[gauss_j]) * mesh.x(i1)
							+ 0.25 * xi[1] * (1 + eta[1] * eta_gauss[gauss_j]) * mesh.x(i2)
							+ 0.25 * xi[2] * (1 + eta[2] * eta_gauss[gauss_j]) * mesh.x(i3)
							+ 0.25 * xi[3] * (1 + eta[3] * eta_gauss[gauss_j]) * mesh.x(i4);

						tax1 = V[i1][0] * A11 + V[i1][1] * A21;
						tax2 = V[i2][0] * A11 + V[i2][1] * A21;
						tax3 = V[i3][0] * A11 + V[i3][1] * A21;
						tax4 = V[i4][0] * A11 + V[i4][1] * A21;

						tay1 = V[i1][0] * A12 + V[i1][1] * A22;
						tay2 = V[i2][0] * A12 + V[i2][1] * A22;
						tay3 = V[i3][0] * A12 + V[i3][1] * A22;
						tay4 = V[i4][0] * A12 + V[i4][1] * A22;

						mat[ie][i][j] += gauss_W[gauss_i] * gauss_W[gauss_j]
							* (
								(tax1 * n_gamma1 + tax2 * n_gamma2 + tax3 * n_gamma3 + tax4 * n_gamma4) * n_alpha * dndxi_beta
							  + (tay1 * n_gamma1 + tay2 * n_gamma2 + tay3 * n_gamma3 + tay4 * n_gamma4) * n_alpha * dndeta_beta);
					}
				}
			}
		}
	}
}

void Advecmatrix::renew(Velocity2d& v) {

	//3点数値積分により移流行列を更新する
	vector<double> xi_gauss = { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };//数値積分離散点xi方向
	vector<double> eta_gauss = { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };//数値積分離散点eta方向
	vector<double> gauss_W = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };//数値積分重み
	double n_alpha, n_beta, n_gamma1, n_gamma2, n_gamma3, n_gamma4, dndxi_beta, dndeta_beta;//形状関数とその微分
	double A11, A12, A21, A22;//余因子要素
	double tax1, tax2, tax3, tax4, tay1, tay2, tay3, tay4;
	for (int ie = 0; ie < mesh.nelem(); ie++) {
		int i1 = mesh.i1(ie);
		int i2 = mesh.i2(ie);
		int i3 = mesh.i3(ie);
		int i4 = mesh.i4(ie);
		for (int i = 0; i < node; i++) {
			for (int j = 0; j < node; j++) {
				for (int gauss_i = 0; gauss_i < 3; gauss_i++) {
					for (int gauss_j = 0; gauss_j < 3; gauss_j++) {
						n_alpha = 0.25 * (1.0 + xi[i] * xi_gauss[gauss_i]) * (1.0 + eta[i] * eta_gauss[gauss_j]);
						n_beta = 0.25 * (1.0 + xi[j] * xi_gauss[gauss_i]) * (1.0 + eta[j] * eta_gauss[gauss_j]);

						n_gamma1 = 0.25 * (1.0 + xi[0] * xi_gauss[gauss_i]) * (1.0 + eta[0] * eta_gauss[gauss_j]);
						n_gamma2 = 0.25 * (1.0 + xi[1] * xi_gauss[gauss_i]) * (1.0 + eta[1] * eta_gauss[gauss_j]);
						n_gamma3 = 0.25 * (1.0 + xi[2] * xi_gauss[gauss_i]) * (1.0 + eta[2] * eta_gauss[gauss_j]);
						n_gamma4 = 0.25 * (1.0 + xi[3] * xi_gauss[gauss_i]) * (1.0 + eta[3] * eta_gauss[gauss_j]);

						dndxi_beta = 0.25 * xi[j] * (1.0 + eta[j] * eta_gauss[gauss_j]);
						dndeta_beta = 0.25 * eta[j] * (1.0 + xi[j] * xi_gauss[gauss_i]);


						A11 = 0.25 * eta[0] * (1 + xi[0] * xi_gauss[gauss_i]) * mesh.y(i1)
							+ 0.25 * eta[1] * (1 + xi[1] * xi_gauss[gauss_i]) * mesh.y(i2)
							+ 0.25 * eta[2] * (1 + xi[2] * xi_gauss[gauss_i]) * mesh.y(i3)
							+ 0.25 * eta[3] * (1 + xi[3] * xi_gauss[gauss_i]) * mesh.y(i4);


						A12 = -0.25 * xi[0] * (1 + eta[0] * eta_gauss[gauss_j]) * mesh.y(i1)
							- 0.25 * xi[1] * (1 + eta[1] * eta_gauss[gauss_j]) * mesh.y(i2)
							- 0.25 * xi[2] * (1 + eta[2] * eta_gauss[gauss_j]) * mesh.y(i3)
							- 0.25 * xi[3] * (1 + eta[3] * eta_gauss[gauss_j]) * mesh.y(i4);

						A21 = -0.25 * eta[0] * (1 + xi[0] * xi_gauss[gauss_i]) * mesh.x(i1)
							- 0.25 * eta[1] * (1 + xi[1] * xi_gauss[gauss_i]) * mesh.x(i2)
							- 0.25 * eta[2] * (1 + xi[2] * xi_gauss[gauss_i]) * mesh.x(i3)
							- 0.25 * eta[3] * (1 + xi[3] * xi_gauss[gauss_i]) * mesh.x(i4);

						A22 = 0.25 * xi[0] * (1 + eta[0] * eta_gauss[gauss_j]) * mesh.x(i1)
							+ 0.25 * xi[1] * (1 + eta[1] * eta_gauss[gauss_j]) * mesh.x(i2)
							+ 0.25 * xi[2] * (1 + eta[2] * eta_gauss[gauss_j]) * mesh.x(i3)
							+ 0.25 * xi[3] * (1 + eta[3] * eta_gauss[gauss_j]) * mesh.x(i4);

						tax1 = v[i1][0] * A11 + v[i1][1] * A21;
						tax2 = v[i2][0] * A11 + v[i2][1] * A21;
						tax3 = v[i3][0] * A11 + v[i3][1] * A21;
						tax4 = v[i4][0] * A11 + v[i4][1] * A21;

						tay1 = v[i1][0] * A12 + v[i1][1] * A22;
						tay2 = v[i2][0] * A12 + v[i2][1] * A22;
						tay3 = v[i3][0] * A12 + v[i3][1] * A22;
						tay4 = v[i4][0] * A12 + v[i4][1] * A22;

						mat[ie][i][j] += gauss_W[gauss_i] * gauss_W[gauss_j]
							* (
								(tax1 * n_gamma1 + tax2 * n_gamma2 + tax3 * n_gamma3 + tax4 * n_gamma4) * n_alpha * dndxi_beta
							  + (tay1 * n_gamma1 + tay2 * n_gamma2 + tay3 * n_gamma3 + tay4 * n_gamma4) * n_alpha * dndeta_beta);
					}
				}
			}
		}
	}
}

Diffmatrix::Diffmatrix(Mesh2d& Mesh)
	:CofficientMatrix(Mesh)
{
	//3点数値積分により拡散行列を求める
	vector<double> xi_gauss = { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };//数値積分離散点xi方向
	vector<double> eta_gauss = { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };//数値積分離散点eta方向
	vector<double> gauss_W = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };//数値積分重み
	double n_alpha, n_beta, dndxi_alpha, dndxi_beta, dndeta_alpha, dndeta_beta;//形状関数とその微分
	double A11, A12, A21, A22;//余因子要素
	double Jacobian;//ヤコビアン

	for (int ie = 0; ie < mesh.nelem(); ie++) {
		int i1 = mesh.i1(ie);
		int i2 = mesh.i2(ie);
		int i3 = mesh.i3(ie);
		int i4 = mesh.i4(ie);
		for (int i = 0; i < node; i++) {
			for (int j = 0; j < node; j++) {
				for (int gauss_i = 0; gauss_i < 3; gauss_i++) {
					for (int gauss_j = 0; gauss_j < 3; gauss_j++) {
						n_alpha = 0.25 * (1.0 + xi[i] * xi_gauss[gauss_i]) * (1.0 + eta[i] * eta_gauss[gauss_j]);
						n_beta = 0.25 * (1.0 + xi[j] * xi_gauss[gauss_i]) * (1.0 + eta[j] * eta_gauss[gauss_j]);
						dndxi_alpha = 0.25 * xi[i] * (1.0 + eta[i] * eta_gauss[gauss_j]);
						dndxi_beta = 0.25 * xi[j] * (1.0 + eta[j] * eta_gauss[gauss_j]);
						dndeta_alpha = 0.25 * eta[i] * (1.0 + xi[i] * xi_gauss[gauss_i]);
						dndeta_beta = 0.25 * eta[j] * (1.0 + xi[j] * xi_gauss[gauss_i]);

						A11 = 0.25 * eta[0] * (1 + xi[0] * xi_gauss[gauss_i]) * mesh.y(i1)
							+ 0.25 * eta[1] * (1 + xi[1] * xi_gauss[gauss_i]) * mesh.y(i2)
							+ 0.25 * eta[2] * (1 + xi[2] * xi_gauss[gauss_i]) * mesh.y(i3)
							+ 0.25 * eta[3] * (1 + xi[3] * xi_gauss[gauss_i]) * mesh.y(i4);


						A12 = -0.25 * xi[0] * (1 + eta[0] * eta_gauss[gauss_j]) * mesh.y(i1)
							- 0.25 * xi[1] * (1 + eta[1] * eta_gauss[gauss_j]) * mesh.y(i2)
							- 0.25 * xi[2] * (1 + eta[2] * eta_gauss[gauss_j]) * mesh.y(i3)
							- 0.25 * xi[3] * (1 + eta[3] * eta_gauss[gauss_j]) * mesh.y(i4);

						A21 = -0.25 * eta[0] * (1 + xi[0] * xi_gauss[gauss_i]) * mesh.x(i1)
							- 0.25 * eta[1] * (1 + xi[1] * xi_gauss[gauss_i]) * mesh.x(i2)
							- 0.25 * eta[2] * (1 + xi[2] * xi_gauss[gauss_i]) * mesh.x(i3)
							- 0.25 * eta[3] * (1 + xi[3] * xi_gauss[gauss_i]) * mesh.x(i4);

						A22 = 0.25 * xi[0] * (1 + eta[0] * eta_gauss[gauss_j]) * mesh.x(i1)
							+ 0.25 * xi[1] * (1 + eta[1] * eta_gauss[gauss_j]) * mesh.x(i2)
							+ 0.25 * xi[2] * (1 + eta[2] * eta_gauss[gauss_j]) * mesh.x(i3)
							+ 0.25 * xi[3] * (1 + eta[3] * eta_gauss[gauss_j]) * mesh.x(i4);

						Jacobian = A11 * A22 - A12 * A21;

						mat[ie][i][j] += gauss_W[gauss_i] * gauss_W[gauss_j] * ((A11 * A11 + A21 * A21) * dndxi_alpha * dndxi_beta
							+ (A12 * A12 + A22 * A22) * dndeta_alpha * dndeta_beta
							+ (A11 * A12 + A21 * A22) * (dndxi_alpha * dndeta_beta + dndeta_alpha * dndxi_beta)) / Jacobian;
					}
				}
			}
		}
	}

}



BTDmatrix::BTDmatrix(Mesh2d& Mesh, Velocity2d& v)
	:CofficientMatrix(Mesh), V(v)
{
	//BTD行列を数値積分で計算
	vector<double> xi_gauss = { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };//数値積分離散点xi方向
	vector<double> eta_gauss = { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };//数値積分離散点eta方向
	vector<double> gauss_W = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };//数値積分重み
	double n_alpha, n_beta, n_gamma1, n_gamma2, n_gamma3, n_gamma4, dndxi_alpha, dndeta_alpha, dndxi_beta, dndeta_beta;//形状関数とその微分
	double A11, A12, A21, A22;//余因子要素
	double tax1, tax2, tax3, tax4, tay1, tay2, tay3, tay4;
	double Jacobian;//ヤコビアン

	for (int ie = 0; ie < mesh.nelem(); ie++) {
		int i1 = mesh.i1(ie);
		int i2 = mesh.i2(ie);
		int i3 = mesh.i3(ie);
		int i4 = mesh.i4(ie);
		for (int i = 0; i < node; i++) {
			for (int j = 0; j < node; j++) {
				for (int gauss_i = 0; gauss_i < 3; gauss_i++) {
					for (int gauss_j = 0; gauss_j < 3; gauss_j++) {

						n_alpha = 0.25 * (1.0 + xi[i] * xi_gauss[gauss_i]) * (1.0 + eta[i] * eta_gauss[gauss_j]);
						n_beta = 0.25 * (1.0 + xi[j] * xi_gauss[gauss_i]) * (1.0 + eta[j] * eta_gauss[gauss_j]);
						n_gamma1 = 0.25 * (1.0 + xi[0] * xi_gauss[gauss_i]) * (1.0 + eta[0] * eta_gauss[gauss_j]);
						n_gamma2 = 0.25 * (1.0 + xi[1] * xi_gauss[gauss_i]) * (1.0 + eta[1] * eta_gauss[gauss_j]);
						n_gamma3 = 0.25 * (1.0 + xi[2] * xi_gauss[gauss_i]) * (1.0 + eta[2] * eta_gauss[gauss_j]);
						n_gamma4 = 0.25 * (1.0 + xi[3] * xi_gauss[gauss_i]) * (1.0 + eta[3] * eta_gauss[gauss_j]);

						dndxi_alpha = 0.25 * xi[i] * (1.0 + eta[i] * eta_gauss[gauss_j]);
						dndxi_beta = 0.25 * xi[j] * (1.0 + eta[j] * eta_gauss[gauss_j]);
						dndeta_alpha = 0.25 * eta[i] * (1.0 + xi[i] * xi_gauss[gauss_i]);
						dndeta_beta = 0.25 * eta[j] * (1.0 + xi[j] * xi_gauss[gauss_i]);

						A11 = 0.25 * eta[0] * (1 + xi[0] * xi_gauss[gauss_i]) * mesh.y(i1)
							+ 0.25 * eta[1] * (1 + xi[1] * xi_gauss[gauss_i]) * mesh.y(i2)
							+ 0.25 * eta[2] * (1 + xi[2] * xi_gauss[gauss_i]) * mesh.y(i3)
							+ 0.25 * eta[3] * (1 + xi[3] * xi_gauss[gauss_i]) * mesh.y(i4);


						A12 = -0.25 * xi[0] * (1 + eta[0] * eta_gauss[gauss_j]) * mesh.y(i1)
							- 0.25 * xi[1] * (1 + eta[1] * eta_gauss[gauss_j]) * mesh.y(i2)
							- 0.25 * xi[2] * (1 + eta[2] * eta_gauss[gauss_j]) * mesh.y(i3)
							- 0.25 * xi[3] * (1 + eta[3] * eta_gauss[gauss_j]) * mesh.y(i4);

						A21 = -0.25 * eta[0] * (1 + xi[0] * xi_gauss[gauss_i]) * mesh.x(i1)
							- 0.25 * eta[1] * (1 + xi[1] * xi_gauss[gauss_i]) * mesh.x(i2)
							- 0.25 * eta[2] * (1 + xi[2] * xi_gauss[gauss_i]) * mesh.x(i3)
							- 0.25 * eta[3] * (1 + xi[3] * xi_gauss[gauss_i]) * mesh.x(i4);

						A22 = 0.25 * xi[0] * (1 + eta[0] * eta_gauss[gauss_j]) * mesh.x(i1)
							+ 0.25 * xi[1] * (1 + eta[1] * eta_gauss[gauss_j]) * mesh.x(i2)
							+ 0.25 * xi[2] * (1 + eta[2] * eta_gauss[gauss_j]) * mesh.x(i3)
							+ 0.25 * xi[3] * (1 + eta[3] * eta_gauss[gauss_j]) * mesh.x(i4);

						tax1 = V[i1][0] * A11 + V[i1][1] * A21;
						tax2 = V[i2][0] * A11 + V[i2][1] * A21;
						tax3 = V[i3][0] * A11 + V[i3][1] * A21;
						tax4 = V[i4][0] * A11 + V[i4][1] * A21;

						tay1 = V[i1][0] * A12 + V[i1][1] * A22;
						tay2 = V[i2][0] * A12 + V[i2][1] * A22;
						tay3 = V[i3][0] * A12 + V[i3][1] * A22;
						tay4 = V[i4][0] * A12 + V[i4][1] * A22;

						Jacobian = A11 * A22 - A12 * A21;

						mat[ie][i][j] += gauss_W[gauss_i] * gauss_W[gauss_j]
							* (
								((tax1 * n_gamma1 + tax2 * n_gamma2 + tax3 * n_gamma3 + tax4 * n_gamma4) * dndxi_alpha
								+ (tay1 * n_gamma1 + tay2 * n_gamma2 + tay3 * n_gamma3 + tay4 + n_gamma4) * dndeta_alpha)
								* ((tax1 * n_gamma1 + tax2 * n_gamma2 + tax3 * n_gamma3 + tax4 * n_gamma4) * dndxi_beta
								+ (tay1 * n_gamma1 + tay2 * n_gamma2 + tay3 * n_gamma3 + tay4 + n_gamma4) * dndeta_beta)
							) / Jacobian;
					}
				}
			}
		}
	}
}

void BTDmatrix::renew(Velocity2d& v) {
	//BTD行列を数値積分で計算
	vector<double> xi_gauss = { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };//数値積分離散点xi方向
	vector<double> eta_gauss = { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };//数値積分離散点eta方向
	vector<double> gauss_W = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };//数値積分重み
	double n_alpha, n_beta, n_gamma1, n_gamma2, n_gamma3, n_gamma4, dndxi_alpha, dndeta_alpha, dndxi_beta, dndeta_beta;//形状関数とその微分
	double A11, A12, A21, A22;//余因子要素
	double tax1, tax2, tax3, tax4, tay1, tay2, tay3, tay4;
	double Jacobian;//ヤコビアン

	for (int ie = 0; ie < mesh.nelem(); ie++) {
		int i1 = mesh.i1(ie);
		int i2 = mesh.i2(ie);
		int i3 = mesh.i3(ie);
		int i4 = mesh.i4(ie);
		for (int i = 0; i < node; i++) {
			for (int j = 0; j < node; j++) {
				for (int gauss_i = 0; gauss_i < 3; gauss_i++) {
					for (int gauss_j = 0; gauss_j < 3; gauss_j++) {

						n_alpha = 0.25 * (1.0 + xi[i] * xi_gauss[gauss_i]) * (1.0 + eta[i] * eta_gauss[gauss_j]);
						n_beta = 0.25 * (1.0 + xi[j] * xi_gauss[gauss_i]) * (1.0 + eta[j] * eta_gauss[gauss_j]);
						n_gamma1 = 0.25 * (1.0 + xi[0] * xi_gauss[gauss_i]) * (1.0 + eta[0] * eta_gauss[gauss_j]);
						n_gamma2 = 0.25 * (1.0 + xi[1] * xi_gauss[gauss_i]) * (1.0 + eta[1] * eta_gauss[gauss_j]);
						n_gamma3 = 0.25 * (1.0 + xi[2] * xi_gauss[gauss_i]) * (1.0 + eta[2] * eta_gauss[gauss_j]);
						n_gamma4 = 0.25 * (1.0 + xi[3] * xi_gauss[gauss_i]) * (1.0 + eta[3] * eta_gauss[gauss_j]);

						dndxi_alpha = 0.25 * xi[i] * (1.0 + eta[i] * eta_gauss[gauss_j]);
						dndxi_beta = 0.25 * xi[j] * (1.0 + eta[j] * eta_gauss[gauss_j]);
						dndeta_alpha = 0.25 * eta[i] * (1.0 + xi[i] * xi_gauss[gauss_i]);
						dndeta_beta = 0.25 * eta[j] * (1.0 + xi[j] * xi_gauss[gauss_i]);

						A11 = 0.25 * eta[0] * (1 + xi[0] * xi_gauss[gauss_i]) * mesh.y(i1)
							+ 0.25 * eta[1] * (1 + xi[1] * xi_gauss[gauss_i]) * mesh.y(i2)
							+ 0.25 * eta[2] * (1 + xi[2] * xi_gauss[gauss_i]) * mesh.y(i3)
							+ 0.25 * eta[3] * (1 + xi[3] * xi_gauss[gauss_i]) * mesh.y(i4);


						A12 = -0.25 * xi[0] * (1 + eta[0] * eta_gauss[gauss_j]) * mesh.y(i1)
							- 0.25 * xi[1] * (1 + eta[1] * eta_gauss[gauss_j]) * mesh.y(i2)
							- 0.25 * xi[2] * (1 + eta[2] * eta_gauss[gauss_j]) * mesh.y(i3)
							- 0.25 * xi[3] * (1 + eta[3] * eta_gauss[gauss_j]) * mesh.y(i4);

						A21 = -0.25 * eta[0] * (1 + xi[0] * xi_gauss[gauss_i]) * mesh.x(i1)
							- 0.25 * eta[1] * (1 + xi[1] * xi_gauss[gauss_i]) * mesh.x(i2)
							- 0.25 * eta[2] * (1 + xi[2] * xi_gauss[gauss_i]) * mesh.x(i3)
							- 0.25 * eta[3] * (1 + xi[3] * xi_gauss[gauss_i]) * mesh.x(i4);

						A22 = 0.25 * xi[0] * (1 + eta[0] * eta_gauss[gauss_j]) * mesh.x(i1)
							+ 0.25 * xi[1] * (1 + eta[1] * eta_gauss[gauss_j]) * mesh.x(i2)
							+ 0.25 * xi[2] * (1 + eta[2] * eta_gauss[gauss_j]) * mesh.x(i3)
							+ 0.25 * xi[3] * (1 + eta[3] * eta_gauss[gauss_j]) * mesh.x(i4);

						tax1 = v[i1][0] * A11 + v[i1][1] * A21;
						tax2 = v[i2][0] * A11 + v[i2][1] * A21;
						tax3 = v[i3][0] * A11 + v[i3][1] * A21;
						tax4 = v[i4][0] * A11 + v[i4][1] * A21;

						tay1 = v[i1][0] * A12 + v[i1][1] * A22;
						tay2 = v[i2][0] * A12 + v[i2][1] * A22;
						tay3 = v[i3][0] * A12 + v[i3][1] * A22;
						tay4 = v[i4][0] * A12 + v[i4][1] * A22;

						Jacobian = A11 * A22 - A12 * A21;

						mat[ie][i][j] += gauss_W[gauss_i] * gauss_W[gauss_j]
							* (
								((tax1 * n_gamma1 + tax2 * n_gamma2 + tax3 * n_gamma3 + tax4 * n_gamma4) * dndxi_alpha
									+ (tay1 * n_gamma1 + tay2 * n_gamma2 + tay3 * n_gamma3 + tay4 + n_gamma4) * dndeta_alpha)
								* ((tax1 * n_gamma1 + tax2 * n_gamma2 + tax3 * n_gamma3 + tax4 * n_gamma4) * dndxi_beta
									+ (tay1 * n_gamma1 + tay2 * n_gamma2 + tay3 * n_gamma3 + tay4 + n_gamma4) * dndeta_beta)
								) / Jacobian;
					}
				}
			}
		}
	}
}

GradientVector::GradientVector(Mesh2d& Mesh) 
	:mesh(Mesh)
{
	vec.resize(mesh.nelem());
	for (int ie = 0; ie < vec.size(); ie++) {
		vec[ie].resize(node);
	}

	//勾配ベクトルを数値積分で計算
	vector<double> xi_gauss = { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };//数値積分離散点xi方向
	vector<double> eta_gauss = { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };//数値積分離散点eta方向
	vector<double> gauss_W = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };//数値積分重み値
	double dndxi_beta, dndeta_beta;//形状関数とその微分
	double A11, A12, A21, A22;//余因子要素


	for (int ie = 0; ie < mesh.nelem(); ie++) {
		int i1 = mesh.i1(ie);
		int i2 = mesh.i2(ie);
		int i3 = mesh.i3(ie);
		int i4 = mesh.i4(ie);
		for (int beta = 0; beta < node; beta++) {
			for (int gauss_i = 0; gauss_i < 3; gauss_i++) {
				for (int gauss_j = 0; gauss_j < 3; gauss_j++) {
					dndxi_beta = 0.25 * xi[beta] * (1.0 + eta[beta] * eta_gauss[gauss_j]);
					dndeta_beta = 0.25 * eta[beta] * (1.0 + xi[beta] * xi_gauss[gauss_i]);

					A11 = 0.25 * eta[0] * (1 + xi[0] * xi_gauss[gauss_i]) * mesh.y(i1)
						+ 0.25 * eta[1] * (1 + xi[1] * xi_gauss[gauss_i]) * mesh.y(i2)
						+ 0.25 * eta[2] * (1 + xi[2] * xi_gauss[gauss_i]) * mesh.y(i3)
						+ 0.25 * eta[3] * (1 + xi[3] * xi_gauss[gauss_i]) * mesh.y(i4);


					A12 = -0.25 * xi[0] * (1 + eta[0] * eta_gauss[gauss_j]) * mesh.y(i1)
						- 0.25 * xi[1] * (1 + eta[1] * eta_gauss[gauss_j]) * mesh.y(i2)
						- 0.25 * xi[2] * (1 + eta[2] * eta_gauss[gauss_j]) * mesh.y(i3)
						- 0.25 * xi[3] * (1 + eta[3] * eta_gauss[gauss_j]) * mesh.y(i4);

					A21 = -0.25 * eta[0] * (1 + xi[0] * xi_gauss[gauss_i]) * mesh.x(i1)
						- 0.25 * eta[1] * (1 + xi[1] * xi_gauss[gauss_i]) * mesh.x(i2)
						- 0.25 * eta[2] * (1 + xi[2] * xi_gauss[gauss_i]) * mesh.x(i3)
						- 0.25 * eta[3] * (1 + xi[3] * xi_gauss[gauss_i]) * mesh.x(i4);

					A22 = 0.25 * xi[0] * (1 + eta[0] * eta_gauss[gauss_j]) * mesh.x(i1)
						+ 0.25 * xi[1] * (1 + eta[1] * eta_gauss[gauss_j]) * mesh.x(i2)
						+ 0.25 * xi[2] * (1 + eta[2] * eta_gauss[gauss_j]) * mesh.x(i3)
						+ 0.25 * xi[3] * (1 + eta[3] * eta_gauss[gauss_j]) * mesh.x(i4);

					vec[ie][beta][0] += gauss_W[gauss_i] * gauss_W[gauss_j] * (A11 * dndxi_beta + A12 * dndeta_beta);//勾配ベクトルx成分
					vec[ie][beta][1] += gauss_W[gauss_i] * gauss_W[gauss_j] * (A21 * dndxi_beta + A22 * dndeta_beta);//勾配ベクトルy成分


				}
			}
		}
	}
}
void GradientVector::view() {

	for (int ie = 0; ie < vec.size(); ie++) {
		cout << "nelem=" << ie << endl;
		for (int i = 0; i < node; i++) {
			cout << "C[" << ie << "][" << i << "]= {" << endl;
			cout << vec[ie][i][0] << endl;
			cout << vec[ie][i][1] << endl;
			cout << "}" << endl;
			cout << endl;
		}

		cout << endl;
	}
	

}
vector<Vector2d>& GradientVector::operator[](int ie) {
	return vec.at(ie);
}
const vector<Vector2d>& GradientVector::operator[](int ie) const {
	return vec.at(ie);
}
Vector2d GradientVector::i1(int ie) {
	return vec.at(ie).at(0);
}
Vector2d GradientVector::i2(int ie) {
	return vec.at(ie).at(1);
}
Vector2d GradientVector::i3(int ie) {
	return vec.at(ie).at(2);
}
Vector2d GradientVector::i4(int ie) {
	return vec.at(ie).at(3);
}
