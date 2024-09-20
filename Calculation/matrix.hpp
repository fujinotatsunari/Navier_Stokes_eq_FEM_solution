#pragma once
#include"matrix.hpp"
#include<iostream>
#include<vector>
using namespace std;
class Matrix {
private:
	vector<vector<double>> mat;
	int row;
	int col;
public:
	//using Index = pair<int, int >;
	Matrix();
	Matrix(int row_, int col_, const vector<vector<double>>& mat_);
	Matrix(int row_, int col);
	Matrix(const Matrix& mat_);
	Matrix& operator=(const Matrix& mat_);
	vector<double> operator[](int i) const;
	vector<double>& operator[](int i);
	//double operator[](int i) const;
	//double operator()(int i, int j) const;
	static Matrix new_zero_matrix(int n_row, int n_col);
	int row_() const;
	int col_() const;
	const double& elem(int i, int j) const;
	double& elem(int i, int j);
	void swap_elem(int i1, int j1, int i2, int j2);
	void swap_rows(int i, int j);
	void print() const;
	void clean();
	void generate(int row_, int col_);

};
//Matrix operator+(const Matrix& A, const Matrix& B);
//Matrix operator-(const Matrix& A, const Matrix& B);
//Matrix operator*(const Matrix& A, const Matrix& B);
//void gaussian_elimination(Matrix& A, vector<double>& x, vector<double>& b);
//void gaussian_elimination(vector<vector<double>>& A, vector<double>& x, vector<double> b);

void LU_solve( Matrix& A, vector<double>& x, vector<double>& b);
void LU_solve( vector<vector<double>>& A, vector<double>& x, vector<double> b);