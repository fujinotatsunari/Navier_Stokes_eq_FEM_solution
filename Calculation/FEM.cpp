#include"FEM.hpp"
#include"value.hpp"
#include"Mesh.hpp"
#include"matrix.hpp"
#include<vector>


CofficientMatrix::CofficientMatrix(Mesh2d& mesh)
	:mesh_(mesh)
{
	mat.resize(mesh_.nelem());
	for (int i = 0; i < mat.size(); i++) {
		mat[i] = Matrix::new_zero_matrix(node, node);
	}
}
CofficientMatrix::CofficientMatrix(Mesh2d& mesh, const vector<Matrix>& Mat)
	:mesh_(mesh), mat(Mat)
{
	mat.resize(mesh_.nelem());
	for (int i = 0; i < mat.size(); i++) {
		mat[i] = Matrix::new_zero_matrix(node, node);
	}

}
const Matrix& CofficientMatrix::operator[](int ie)const {
	return mat[ie];
}
Matrix& CofficientMatrix::operator[](int ie) {
	return mat[ie];
}
CofficientMatrix CofficientMatrix::generate(Mesh2d& mesh) {
	auto Mesh = Mesh2d(mesh);
	vector<Matrix> Mat;
	return CofficientMatrix(Mesh, Mat);
}
void CofficientMatrix::view() {
	for (int ie = 0; ie < mat.size(); ie++) {
		cout << "nelem=" << ie << endl;
		mat[ie].print();
		cout << endl;
	}
}

Massmatrix::Massmatrix(Mesh2d& mesh)
	:CofficientMatrix(mesh)
{
	//質量行列を求める
	for (int ie = 0; ie < mesh_.nelem(); ie++) {
		for (int i = 0; i < node; i++) {
			for (int j = 0; j < node; j++) {
				for (int k = 0; k < node; k++) {
					for (int l = 0; l < node; l++) {
						mat[ie][i][j] += (3.0 + xi[i] * xi[j]) * (3.0 + eta[i] * eta[j])
							* ((xi[k] * eta[l] - eta[k] * xi[l]) * mesh_.x(mesh_.nbool1(ie, k)) * mesh_.y(mesh_.nbool1(ie, l))
								+ (xi[i] + xi[j]) * (xi[k] * xi[l] * (eta[l] - eta[k]) * mesh_.x(mesh_.nbool1(ie, k)) * mesh_.y(mesh_.nbool1(ie, l))) / (3.0 + xi[i] * xi[j])
								+ (eta[i] + eta[j]) * (eta[k] * eta[l] * (xi[k] - xi[l]) * mesh_.x(mesh_.nbool1(ie, k)) * mesh_.y(mesh_.nbool1(ie, l))) / (3.0 + eta[i] * eta[j])) / 576.0;

					}
				}
			}
		}
	}


}
Massmatrix::Massmatrix(Mesh2d& mesh, const vector<Matrix>& Mat)
	:CofficientMatrix(mesh, Mat)
{//質量行列を求める
	for (int ie = 0; ie < mesh_.nelem(); ie++) {
		for (int i = 0; i < node; i++) {
			for (int j = 0; j < node; j++) {
				for (int k = 0; k < node; k++) {
					for (int l = 0; l < node; l++) {
						mat[ie][i][j] += (3.0 + xi[i] * xi[j]) * (3.0 + eta[i] * eta[j])
							* ((xi[k] * eta[l] - eta[k] * xi[l]) * mesh_.x(mesh_.nbool1(ie, k)) * mesh_.y(mesh_.nbool1(ie, l))
								+ (xi[i] + xi[j]) * (xi[k] * xi[l] * (eta[l] - eta[k]) * mesh_.x(mesh_.nbool1(ie, k)) * mesh_.y(mesh_.nbool1(ie, l))) / (3.0 + xi[i] * xi[j])
								+ (eta[i] + eta[j]) * (eta[k] * eta[l] * (xi[k] - xi[l]) * mesh_.x(mesh_.nbool1(ie, k)) * mesh_.y(mesh_.nbool1(ie, l))) / (3.0 + eta[i] * eta[j])) / 576.0;

					}
				}
			}
		}
	}
}
Massmatrix Massmatrix::generate_Mass(Mesh2d& mesh) {
	auto Mesh = Mesh2d(mesh);
	vector<Matrix> Mat;
	int node = 4;
	Mat.resize(Mesh.nelem());
	for (int i = 0; i < Mat.size(); i++) {
		Mat[i] = Matrix::new_zero_matrix(node, node);
	}

	return Massmatrix(Mesh, Mat);
}
/*
const Matrix& Massmatrix::operator[](int ie) const {
	return mat[ie];
}
Matrix& Massmatrix::operator[](int ie) {
	return mat[ie];
}*/
Lumped_Massmatrix::Lumped_Massmatrix(Mesh2d& mesh)
	:CofficientMatrix(mesh)
{
	//集中化質量行列を求める
	double EM1;
	double EM2;
	double EM3;
	for (int ie = 0; ie < mesh_.nelem(); ie++) {
		for (int i = 0; i < node; i++) {
			EM1 = 0.0;
			EM2 = 0.0;
			EM3 = 0.0;
			for (int k = 0; k < node; k++) {
				for (int l = 0; l < node; l++) {

					EM1 = EM1 + xi[k] * eta[l] * (mesh_.x(mesh_.nbool1(ie, k)) * mesh_.y(mesh_.nbool1(ie, l)) - mesh_.y(mesh_.nbool1(ie, k)) * mesh_.x(mesh_.nbool1(ie, l)));

					EM2 = EM2 + xi[k] * xi[l] * eta[l] * (mesh_.x(mesh_.nbool1(ie, k)) * mesh_.y(mesh_.nbool1(ie, l)) - mesh_.y(mesh_.nbool1(ie, k)) * mesh_.x(mesh_.nbool1(ie, l)));

					EM3 = EM3 + xi[k] * eta[k] * eta[l] * (mesh_.x(mesh_.nbool1(ie, k)) * mesh_.y(mesh_.nbool1(ie, l)) - mesh_.y(mesh_.nbool1(ie, k)) * mesh_.x(mesh_.nbool1(ie, l)));
				}
			}
			mat[ie][i][i] = EM1 / 16.0 + xi[i] * EM2 / 48.0 + eta[i] * EM3 / 48.0;

		}
	}
}
Lumped_Massmatrix::Lumped_Massmatrix(Mesh2d& mesh, const vector<Matrix>& Mat)
	:CofficientMatrix(mesh, Mat)
{
	//集中化質量行列を求める
	for (int ie = 0; ie < mesh_.nelem(); ie++) {
		for (int i = 0; i < node; i++) {
			for (int k = 0; k < node; k++) {
				for (int l = 0; l < node; l++) {
					double EM1;
					EM1 = xi[k] * eta[l] * (mesh_.x(mesh_.nbool1(ie, k)) * mesh_.y(mesh_.nbool1(ie, l)) - mesh_.y(mesh_.nbool1(ie, k)) * mesh_.x(mesh_.nbool1(ie, l)));
					double EM2;
					EM2 = xi[k] * xi[l] * eta[l] * (mesh_.x(mesh_.nbool1(ie, k)) * mesh_.y(mesh_.nbool1(ie, l)) - mesh_.y(mesh_.nbool1(ie, k)) * mesh_.x(mesh_.nbool1(ie, l)));
					double EM3;
					EM3 = xi[k] * eta[k] * eta[l] * (mesh_.x(mesh_.nbool1(ie, k)) * mesh_.y(mesh_.nbool1(ie, l)) - mesh_.y(mesh_.nbool1(ie, k)) * mesh_.x(mesh_.nbool1(ie, l)));

					mat[ie][i][i] = EM1 / 16.0 + xi[i] * EM2 / 48.0 + eta[i] * EM3 / 48.0;

				}
			}

		}
	}
}
Lumped_Massmatrix Lumped_Massmatrix::generate_Lmass(Mesh2d& mesh) {
	auto Mesh = Mesh2d(mesh);
	vector<Matrix> Mat;
	int node = 4;
	Mat.resize(Mesh.nelem());
	for (int i = 0; i < Mat.size(); i++) {
		Mat[i] = Matrix::new_zero_matrix(node, node);
	}

	return Lumped_Massmatrix(Mesh, Mat);
}
/*
const Matrix& Lumped_Massmatrix::operator[](int ie) const {
	return mat[ie];
}
Matrix& Lumped_Massmatrix::operator[](int ie) {
	return mat[ie];
}
*/
xAdvecmatrix::xAdvecmatrix(Mesh2d& mesh)
	:CofficientMatrix(mesh)
{
	//x方向移流行列行列を求める
	for (int ie = 0; ie < mesh_.nelem(); ie++) {
		for (int i = 0; i < node; i++) {
			for (int j = 0; j < node; j++) {
				for (int k = 0; k < node; k++) {
					mat[ie][i][j] += xi[k] * eta[k] * (xi[i] * xi[j] - eta[i] * eta[j]) * mesh_.y(mesh_.nbool1(ie, k)) / 48.0
						+ xi[j] * eta[k] * (1 + eta[i] * eta[j] / 3.0) * mesh_.y(mesh_.nbool1(ie, k)) / 16.0
						- eta[j] * xi[k] * (1 + xi[i] * xi[j] / 3.0) * mesh_.y(mesh_.nbool1(ie, k)) / 16.0;

				}
			}
		}
	}
}
xAdvecmatrix::xAdvecmatrix(Mesh2d& mesh, const vector<Matrix>& Mat)
	:CofficientMatrix(mesh, Mat)
{
	//x方向移流行列行列を求める
	for (int ie = 0; ie < mesh_.nelem(); ie++) {
		for (int i = 0; i < node; i++) {
			for (int j = 0; j < node; j++) {
				for (int k = 0; k < node; k++) {
					mat[ie][i][j] += xi[k] * eta[k] * (xi[i] * xi[j] - eta[i] * eta[j]) * mesh_.y(mesh_.nbool1(ie, k)) / 48.0
						+ xi[j] * eta[k] * (1 + eta[i] * eta[j] / 3.0) * mesh_.y(mesh_.nbool1(ie, k)) / 16.0
						- eta[j] * xi[k] * (1 + xi[i] * xi[j] / 3.0) * mesh_.y(mesh_.nbool1(ie, k)) / 16.0;

				}
			}
		}
	}
}
xAdvecmatrix xAdvecmatrix::generate_xAd(Mesh2d& mesh) {
	auto Mesh = Mesh2d(mesh);
	vector<Matrix> Mat;
	int node = 4;
	Mat.resize(Mesh.nelem());
	for (int i = 0; i < Mat.size(); i++) {
		Mat[i] = Matrix::new_zero_matrix(node, node);
	}

	return xAdvecmatrix(Mesh, Mat);
}
/*
const Matrix& xAdvecmatrix::operator[](int ie) const {
	return mat[ie];
}
Matrix& xAdvecmatrix::operator[](int ie) {
	return mat[ie];
}
*/
yAdvecmatrix::yAdvecmatrix(Mesh2d& mesh)
	:CofficientMatrix(mesh)
{
	//y方向移流行列行列を求める
	for (int ie = 0; ie < mesh_.nelem(); ie++) {
		for (int i = 0; i < node; i++) {
			for (int j = 0; j < node; j++) {
				for (int k = 0; k < node; k++) {
					mat[ie][i][j] += -xi[k] * eta[k] * (xi[i] * xi[j] - eta[i] * eta[j]) * mesh_.x(mesh_.nbool1(ie, k)) / 48.0
						- xi[j] * eta[k] * (1 + eta[i] * eta[j] / 3.0) * mesh_.x(mesh_.nbool1(ie, k)) / 16.0
						+ eta[j] * xi[k] * (1 + xi[i] * xi[j] / 3.0) * mesh_.x(mesh_.nbool1(ie, k)) / 16.0;

				}
			}
		}
	}
}
yAdvecmatrix::yAdvecmatrix(Mesh2d& mesh, const vector<Matrix>& Mat)
	:CofficientMatrix(mesh, Mat)
{
	//y方向移流行列行列を求める
	for (int ie = 0; ie < mesh_.nelem(); ie++) {
		for (int i = 0; i < node; i++) {
			for (int j = 0; j < node; j++) {
				for (int k = 0; k < node; k++) {
					mat[ie][i][j] += -xi[k] * eta[k] * (xi[i] * xi[j] - eta[i] * eta[j]) * mesh_.x(mesh_.nbool1(ie, k)) / 48.0
						- xi[j] * eta[k] * (1 + eta[i] * eta[j] / 3.0) * mesh_.x(mesh_.nbool1(ie, k)) / 16.0
						+ eta[j] * xi[k] * (1 + xi[i] * xi[j] / 3.0) * mesh_.x(mesh_.nbool1(ie, k)) / 16.0;

				}
			}
		}
	}
}
yAdvecmatrix yAdvecmatrix::generate_yAd(Mesh2d& mesh) {
	auto Mesh = Mesh2d(mesh);
	vector<Matrix> Mat;
	int node = 4;
	Mat.resize(Mesh.nelem());
	for (int i = 0; i < Mat.size(); i++) {
		Mat[i] = Matrix::new_zero_matrix(node, node);
	}

	return yAdvecmatrix(Mesh, Mat);
}
/*
const Matrix& yAdvecmatrix::operator[](int ie) const {
	return mat[ie];
}
Matrix& yAdvecmatrix::operator[](int ie) {
	return mat[ie];
}
*/
Diffmatrix::Diffmatrix(Mesh2d& mesh)
	:CofficientMatrix(mesh)
{
	//3点数値積分により拡散行列を求める
	vector<double> xi_gauss = { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };//数値積分離散点xi方向
	vector<double> eta_gauss = { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };//数値積分離散点eta方向
	vector<double> gauss_W = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };//数値積分重み
	double n_alpha, n_beta, dndxi_alpha, dndxi_beta, dndeta_alpha, dndeta_beta;//形状関数とその微分
	double A11, A12, A21, A22;//余因子要素
	double Jacobian;//ヤコビアン

	for (int ie = 0; ie < mesh_.nelem(); ie++) {
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

						A11 = 0.25 * eta[0] * (1 + xi[0] * xi_gauss[gauss_i]) * mesh_.y(mesh_.i1(ie))
							+ 0.25 * eta[1] * (1 + xi[1] * xi_gauss[gauss_i]) * mesh_.y(mesh_.i2(ie))
							+ 0.25 * eta[2] * (1 + xi[2] * xi_gauss[gauss_i]) * mesh_.y(mesh_.i3(ie))
							+ 0.25 * eta[3] * (1 + xi[3] * xi_gauss[gauss_i]) * mesh_.y(mesh_.i4(ie));


						A12 = -0.25 * xi[0] * (1 + eta[0] * eta_gauss[gauss_j]) * mesh_.y(mesh_.i1(ie))
							- 0.25 * xi[1] * (1 + eta[1] * eta_gauss[gauss_j]) * mesh_.y(mesh_.i2(ie))
							- 0.25 * xi[2] * (1 + eta[2] * eta_gauss[gauss_j]) * mesh_.y(mesh_.i3(ie))
							- 0.25 * xi[3] * (1 + eta[3] * eta_gauss[gauss_j]) * mesh_.y(mesh_.i4(ie));

						A21 = -0.25 * eta[0] * (1 + xi[0] * xi_gauss[gauss_i]) * mesh_.x(mesh_.i1(ie))
							- 0.25 * eta[1] * (1 + xi[1] * xi_gauss[gauss_i]) * mesh_.x(mesh_.i2(ie))
							- 0.25 * eta[2] * (1 + xi[2] * xi_gauss[gauss_i]) * mesh_.x(mesh_.i3(ie))
							- 0.25 * eta[3] * (1 + xi[3] * xi_gauss[gauss_i]) * mesh_.x(mesh_.i4(ie));

						A22 = 0.25 * xi[0] * (1 + eta[0] * eta_gauss[gauss_j]) * mesh_.x(mesh_.i1(ie))
							+ 0.25 * xi[1] * (1 + eta[1] * eta_gauss[gauss_j]) * mesh_.x(mesh_.i2(ie))
							+ 0.25 * xi[2] * (1 + eta[2] * eta_gauss[gauss_j]) * mesh_.x(mesh_.i3(ie))
							+ 0.25 * xi[3] * (1 + eta[3] * eta_gauss[gauss_j]) * mesh_.x(mesh_.i4(ie));

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
Diffmatrix::Diffmatrix(Mesh2d& mesh, const vector<Matrix>& Mat)
	:CofficientMatrix(mesh, Mat)
{
	//3点数値積分により拡散行列を求める
	vector<double> xi_gauss = { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };//数値積分離散点xi方向
	vector<double> eta_gauss = { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };//数値積分離散点eta方向
	vector<double> gauss_W = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };//数値積分重み
	double n_alpha, n_beta, dndxi_alpha, dndxi_beta, dndeta_alpha, dndeta_beta;//形状関数とその微分
	double A11, A12, A21, A22;//余因子要素
	double Jacobian;//ヤコビアン

	for (int ie = 0; ie < mesh_.nelem(); ie++) {
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

						A11 = 0.25 * eta[0] * (1 + xi[0] * xi_gauss[gauss_i]) * mesh_.y(mesh_.i1(ie))
							+ 0.25 * eta[1] * (1 + xi[1] * xi_gauss[gauss_i]) * mesh_.y(mesh_.i2(ie))
							+ 0.25 * eta[2] * (1 + xi[2] * xi_gauss[gauss_i]) * mesh_.y(mesh_.i3(ie))
							+ 0.25 * eta[3] * (1 + xi[3] * xi_gauss[gauss_i]) * mesh_.y(mesh_.i4(ie));


						A12 = -0.25 * xi[0] * (1 + eta[0] * eta_gauss[gauss_j]) * mesh_.y(mesh_.i1(ie))
							- 0.25 * xi[1] * (1 + eta[1] * eta_gauss[gauss_j]) * mesh_.y(mesh_.i2(ie))
							- 0.25 * xi[2] * (1 + eta[2] * eta_gauss[gauss_j]) * mesh_.y(mesh_.i3(ie))
							- 0.25 * xi[3] * (1 + eta[3] * eta_gauss[gauss_j]) * mesh_.y(mesh_.i4(ie));

						A21 = -0.25 * eta[0] * (1 + xi[0] * xi_gauss[gauss_i]) * mesh_.x(mesh_.i1(ie))
							- 0.25 * eta[1] * (1 + xi[1] * xi_gauss[gauss_i]) * mesh_.x(mesh_.i2(ie))
							- 0.25 * eta[2] * (1 + xi[2] * xi_gauss[gauss_i]) * mesh_.x(mesh_.i3(ie))
							- 0.25 * eta[3] * (1 + xi[3] * xi_gauss[gauss_i]) * mesh_.x(mesh_.i4(ie));

						A22 = 0.25 * xi[0] * (1 + eta[0] * eta_gauss[gauss_j]) * mesh_.x(mesh_.i1(ie))
							+ 0.25 * xi[1] * (1 + eta[1] * eta_gauss[gauss_j]) * mesh_.x(mesh_.i2(ie))
							+ 0.25 * xi[2] * (1 + eta[2] * eta_gauss[gauss_j]) * mesh_.x(mesh_.i3(ie))
							+ 0.25 * xi[3] * (1 + eta[3] * eta_gauss[gauss_j]) * mesh_.x(mesh_.i4(ie));

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
Diffmatrix Diffmatrix::generate_Diff(Mesh2d& mesh) {
	auto Mesh = Mesh2d(mesh);
	vector<Matrix> Mat;
	int node = 4;
	Mat.resize(Mesh.nelem());
	for (int i = 0; i < Mat.size(); i++) {
		Mat[i] = Matrix::new_zero_matrix(node, node);
	}

	return Diffmatrix(Mesh, Mat);
}
/*
const Matrix& Diffmatrix::operator[](int ie) const {
	return mat[ie];
}
Matrix& Diffmatrix::operator[](int ie){
	return mat[ie];
}
*/