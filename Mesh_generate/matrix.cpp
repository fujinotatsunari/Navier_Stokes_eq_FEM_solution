#include"matrix.hpp"
#include<iostream>
#include<vector>
using namespace std;
Matrix::Matrix() 
	:row(0),col(0),mat(NULL)
{


}
Matrix::Matrix(int row_, int col_, const vector<vector<double>>& mat_)
	:row(row_), col(col_) 
{
	mat = mat_;
	mat.resize(row);
	for (int i = 0; i < mat.size(); i++) {
		mat[i].resize(col);
	}
}
Matrix::Matrix(int row_, int col_)
	:row(row_), col(col_) 
{
	vector<vector<double>> mat_;
	mat_.resize(row);
	for (int i = 0; i < mat_.size(); i++) {
		mat_[i].resize(col);
	}
	mat = mat_;

}
Matrix::Matrix(const Matrix& mat_) {
	row = mat_.row;
	col = mat_.col;
	mat = mat_.mat;
}
Matrix& Matrix::operator=(const Matrix& mat_) {
	
	row = mat_.row;
	col = mat_.col;
	mat = mat_.mat;
	return *this;
}
vector<double> Matrix::operator[](int i) const{
	return mat.at(i);
}
vector<double>& Matrix::operator[](int i) {
	return mat.at(i);
}
/*
double& Matrix::operator()(int i, int j) {
	return mat[i][j];
}*/
int Matrix::row_() const {
	return row;
}
int Matrix::col_() const {
	return col;
}
const double& Matrix::elem(int i, int j) const  {
	return mat[i][j];
}
double& Matrix::elem(int i, int j) {
	return mat[i][j];
}
Matrix Matrix::new_zero_matrix(int n_row, int n_col) {
	vector<vector<double>> mat_(n_row, vector<double>(n_col, 0.0));
	return Matrix(n_row, n_col, mat_);
}
void Matrix::swap_elem(int i1, int j1, int i2, int j2) {
	double tmp = elem(i1, j1);
	elem(i1, j1) = elem(i2, j2);
	elem(i2, j2) = tmp;
}
void Matrix::swap_rows(int i, int j) {
	for (int k = 0; k < row; k++) {
		swap_elem(i, k, j, k);
	}
}
void Matrix::print() const{
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			cout << "         " << mat[i][j] << "         ";
		}
		cout << endl;
	}
}
void Matrix::clean() {
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			mat[i][j] = 0.0;
		}
	}
}
void Matrix::generate(int row_, int col_) {
	row = row_;
	col = col_;
	mat.resize(row);
	for (int i = 0; i < mat.size(); i++) {
		mat[i].resize(col);
	}
	clean();
}
/*
Matrix operator+(const Matrix& A, const Matrix& B) {
	int n_row = A.row_();
	int n_col = B.col_();
	auto out = Matrix::new_zero_matrix(n_row, n_col);
	for (int i = 0; i < n_row; i++) {
		for (int j = 0; j < n_col; j++) {
			out[i][j] = A[i][j] + B[i][j];
		}
	}
	return out;
}
Matrix operator-(const Matrix& A, const Matrix& B) {
	int n_row = A.row_();
	int n_col = B.col_();
	auto out = Matrix::new_zero_matrix(n_row, n_col);
	for (int i = 0; i < n_row; i++) {
		for (int j = 0; j < n_col; j++) {
			out[i][j] = A[i][j] - B[i][j];
		}
	}
	return out;
}
Matrix operator*(const Matrix& A, const Matrix& B) {
	int n_row = A.row_();
	int n_col = B.col_();
	int inner = A.col_();
	auto out = Matrix::new_zero_matrix(n_row, n_col);
	for (int i = 0; i < n_row; i++) {
		for (int j = 0; j < n_col; j++) {
			for (int k = 0; k < inner; k++) {
				out[i][j] = out[i][j] + A[i][k] * B[k][j];
			}
		}
	}
	return out;
}*/
void LU_solve(Matrix& A, vector<double>& x, vector<double>& b) {

	/* | |1(=L_00)|  2   3   4   5  =(U1->)|
	   | ----------------------------------|
	   | |2       |  3   4   5   6         |
	   | |3       |  4   5   6   7         |
	   | |4       |  5   6   7   8         |
	   | |5       |  6   7   8   9         |
	   | |(=L1->) |              =(A1)     |*/

	//cout << "===========  in LU_solve  ==============" << endl;
	int n_r = (int)A.row_();//Aの行
	int n_c = (int)A.col_();//Aの列
	int m_r = (int)b.size();//bの行
	int m_c = 1;//bの列=1
	//auto Amat_= Matrix::new_zero_matrix(n_r, n_c);
	auto L = Matrix::new_zero_matrix(n_r, n_c);
	auto U = Matrix::new_zero_matrix(n_r, n_c);
	auto LU = Matrix::new_zero_matrix(n_r, n_c);

	//A = A;//引数のAを保存
	for (int i = 0; i < n_r; i++) {
		int n = n_r - i - 1;

		L[i][i] = A[0][0];
		double L0 = L[i][i];//L_00成分をコピー

		//L1成分をコピー
		vector<double> L1(n);
		for (int j = 0; j < n; j++) {
			L1[j] = A[j + 1][0];
			L[j + i + 1][i] = L1[j];
		}
		//U1成分をコピー
		vector<double> U1(n);
		for (int j = 0; j < n; j++) {
			U1[j] = A[0][j + 1] / L0;
			U[i][j + i + 1] = U1[j];
		}
		//LUを求める
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				LU[j][k] = L1[j] * U1[k];
			}
		}
		//A1を求める(n-1次行列)
		vector<vector<double>> A1_;
		Matrix A1(n, n, A1_);
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				A1[j][k] = A[j + 1][k + 1] - LU[j][k];
			}
		}
		//A1をAとして再帰的に解く
		A = A1;
		/*
		cout << "A_" << n << "matrix" << endl;
		A.print();
		cout << endl;
		*/
	}
	/*
	cout << "L_matrix" << endl;
	L.print();
	cout << endl;

	cout << "U_matrix" << endl;
	U.print();
	cout << endl;

	cout << "LU_matrix" << endl;
	LU.print();
	cout << endl;
	*/

	//求めたLU行列で連立方程式を解く

	vector<double> y(n_r);
	for (int i = 0; i < n_r; i++) {
		double sum = 0;
		for (int k = 0; k <= i - 1; k++) {
			sum += L[i][k] * y[k];
		}
		y[i] = (b[i] - sum) / L[i][i];
	}
	vector<double> x_(n_r);
	for (int i = n_r - 1; i >= 0; i--) {
		double sum = 0;
		for (int k = i + 1; k <= n_r - 1; k++) {
			sum += U[i][k] * x_[k];
		}
		x_[i] = y[i] - sum;
	}
	x = x_;
	//cout << "LU_solve_end" << endl;

}
void LU_solve(vector<vector<double>>& Amat, vector<double>& x, vector<double> b) {

	/* | |1(=L_00)|  2   3   4   5  =(U1->)|
	   | ----------------------------------|
	   | |2       |  3   4   5   6         |
	   | |3       |  4   5   6   7         |
	   | |4       |  5   6   7   8         |
	   | |5       |  6   7   8   9         |
	   | |(=L1->) |              =(A1)     |*/


	int n_r = (int)Amat.size();//Aの行
	int n_c = (int)Amat[0].size();//Aの列
	int m_r = (int)b.size();//bの行
	int m_c = 1;//bの列=1
	Matrix A_(n_r, n_c, Amat);
	auto L = Matrix::new_zero_matrix(n_r, n_c);
	auto U = Matrix::new_zero_matrix(n_r, n_c);
	auto LU = Matrix::new_zero_matrix(n_r, n_c);


	for (int i = 0; i < n_r; i++) {
		int n = n_r - i - 1;

		L[i][i] = A_[0][0];
		double L0 = L[i][i];//L_00成分をコピー

		//L1成分をコピー
		vector<double> L1(n);
		for (int j = 0; j < n; j++) {
			L1[j] = A_[j + 1][0];
			L[j + i + 1][i] = L1[j];
		}
		//U1成分をコピー
		vector<double> U1(n);
		for (int j = 0; j < n; j++) {
			U1[j] = A_[0][j + 1] / L0;
			U[i][j + i + 1] = U1[j];
		}
		//LUを求める
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				LU[j][k] = L1[j] * U1[k];
			}
		}
		//A1を求める(n-1次行列)
		vector<vector<double>> A1_;
		Matrix A1(n, n, A1_);
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				A1[j][k] = A_[j + 1][k + 1] - LU[j][k];
			}
		}
		//A1をAとして再帰的に解く
		A_ = A1;
	}

	//求めたLU行列で連立方程式を解く

	vector<double> y(n_r);
	for (int i = 0; i < n_r; i++) {
		double sum = 0;
		for (int k = 0; k <= i - 1; k++) {
			sum += L[i][k] * y[k];
		}
		y[i] = (b[i] - sum) / L[i][i];
	}
	vector<double> x_(n_r);
	for (int i = n_r - 1; i >= 0; i--) {
		double sum = 0;
		for (int k = i + 1; k <= n_r - 1; k++) {
			sum += U[i][k] * x_[k];
		}
		x_[i] = y[i] - sum;
	}
	x = x_;
}