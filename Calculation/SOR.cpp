#include "SOR.hpp"
#include "mesh.hpp"
#include "value.hpp"
#include "FEM.hpp"
#include <vector>
using namespace std;

SORparam::SORparam(Mesh2d& Mesh) 
	:mesh(Mesh),nmax(0),eps(0.0)
{
	cout << "Object generate: SORparameter" << endl;
	cout << "同時緩和法最大反復回数 nmax ->";
	cin >> nmax;
	cout << "同時緩和法収束判定値 eps ->";
	cin >> eps;
	int nelem = mesh.nelem();

	lambda.resize(nelem);
	for (int ie = 0; ie < nelem; ie++) {

		vector<double> A(4);//速度修正ポテンシャルの係数値

		int Ed = mesh.nbool3(ie, 0);//下側要素番号
		int Er = mesh.nbool3(ie, 1);//右側要素番号
		int Eu = mesh.nbool3(ie, 2);//上側要素番号
		int El = mesh.nbool3(ie, 3);//左側要素番号
		int n1 = mesh.i1(ie);//要素左下点番号
		int n2 = mesh.i2(ie);//要素右下点番号
		int n3 = mesh.i3(ie);//要素右上点番号
		int n4 = mesh.i4(ie);//要素左上点番号
		int nflagd = 0;//neumann境界(流出壁面)フラグ(下側)
		int nflagr = 0;//neumann境界(流出壁面)フラグ(右側)
		int nflagu = 0;//neumann境界(流出壁面)フラグ(上側)
		int nflagl = 0;//neumann境界(流出壁面)フラグ(左側)

	
		double Sd = 0.0;//要素重心と下の要素重心と下辺の端からなる四角形の面積
		double Sr = 0.0;//要素重心と右の要素重心と右辺の端からなる四角形の面積
		double Su = 0.0;//要素重心と上の要素重心と上辺の端からなる四角形の面積
		double Sl = 0.0;//要素重心と左の要素重心と左辺の端からなる四角形の面積
		double sd = mesh.length(n1, n2);//要素下辺長さ
		double sr = mesh.length(n2, n3);//要素右辺長さ
		double su = mesh.length(n3, n4);//要素上辺長さ
		double sl = mesh.length(n4, n1);//要素左辺長さ


		if (mesh.scond(ie) != 1) {//ieが障害物要素ではない

			if (Ed != -1) {//Ed(下側要素)が範囲外でない
				if (mesh.scond(Ed) != 1) {//下側要素が障害物要素でない
					Sd = mesh.area(ie, Ed);
					A[0] = sd * sd / (2 * Sd);
				}
				else {//下側要素が障害物要素
					A[0] = 0.0;
				}
			}
			else {//Edが範囲外
				if (mesh.ncond(n1) == 3 && mesh.ncond(n2) == 3) {//Edに接する点が流出(neumann)壁面である
					nflagd = 1;
				}
				else {//Edに接する点が流出壁面以外
					A[0] = 0.0;
				}
			}

			if (Er != -1) {//Er(右側要素)が範囲外でない
				if (mesh.scond(Er) != 1) {//右側要素が障害物要素ではない
					Sr = mesh.area(ie, Er);
					A[1] = sr * sr / (2 * Sr);
				}
				else {//右側要素は障害物要素
					A[1] = 0.0;
				}
			}
			else {//Erが範囲外
				if (mesh.ncond(n2) == 3 && mesh.ncond(n3) == 3) {//Erに接する点が流出(neumann)壁面である
					nflagr = 1;
				}
				else {//Erに接する点が流出壁面以外
					A[1] = 0.0;
				}
			}

			if (Eu != -1) {//Euが範囲外でない
				if (mesh.scond(Eu) != 1) {//上側要素が障害物要素でない
					Su = mesh.area(ie, Eu);
					A[2] = su * su / (2 * Su);
				}
				else {//上側要素が障害物要素
					A[2] = 0.0;
				}
			}
			else {//Euが範囲外
				if (mesh.ncond(n3) == 3 && mesh.ncond(n4) == 3) {//Euに接する点が流出(neumann)壁面である
					nflagu = 1;
				}
				else {//Euに接する点が流出壁面以外
					A[2] = 0.0;
				}
			}

			if (El != -1) {//Elが範囲外でない
				if (mesh.scond(El) != 1) {//左側要素が障害物要素でない
					Sl = mesh.area(ie, El);
					A[3] = sl * sl / (2 * Sl);
				}
				else {//左側要素が障害物要素
					A[3] = 0.0;
				}	
			}
			else {//Elが範囲外
				if (mesh.ncond(n4) == 3 && mesh.ncond(n1) == 3) {//Elに接する点が流出(neumann)壁面である
					nflagl = 1;
				}
				else {//Elに接する点が流出壁面以外
					A[3] = 0.0;
				}
			} 
			// 1辺がneumann境界条件のとき
			if (nflagd == 1) {
				A[0] = A[2];//対面のAをコピーする
			}
			if (nflagr == 1) {
				A[1] = A[3];
			}
			if (nflagu == 1) {
				A[2] = A[0];
			}
			if (nflagl == 1) {
				A[3] = A[1];
			}
			cout << ie << endl;

			cout << "A0=" << A[0] << " A1=" << A[1] << " A2=" << A[2] << " A3=" << A[3] << endl;
			cout << "sd=" << sd << " sr=" << sr << " su=" << su << " sl=" << sl << endl;
			cout << "Sd=" << Sd << " Sr=" << Sr << " Su=" << Su << " Sl=" << Sl << endl;
			lambda[ie] = (A[0] + A[1] + A[2] + A[3]) / mesh.Se(ie);
			cout << "lambda[" << ie << "]=" << lambda[ie] << endl;

		}
		else {
			lambda[ie] = -1;
		}
		A.clear();
		
	}
}
int SORparam::get_nmax() {
	return nmax;
}
double SORparam::get_eps() {
	return eps;
}
double SORparam::get_lambda(int ie) {
	return lambda[ie];
}