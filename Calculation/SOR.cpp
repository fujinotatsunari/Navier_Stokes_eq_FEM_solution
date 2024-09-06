#include "SOR.hpp"
#include "mesh.hpp"
#include "value.hpp"
#include "FEM.hpp"
#include <vector>
using namespace std;

SORparam::SORparam(Mesh2d& Mesh) 
	:mesh(Mesh)
{
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
			if (Ed != -1 && mesh.scond(Ed) != 1) {//Edが範囲外でないかつ障害物要素ではない
				Sd = mesh.area(ie, Ed);
				A[0] = sd * sd / (2 * Sd);
			}
			else {//Edが範囲外または障害物要素
				if (mesh.ncond(n1) != 0 || mesh.ncond(n2) != 0) {//Edに接する点が内部点ではない
					if (mesh.ncond(n1) == 3 && mesh.ncond(n2) == 3) {//Edに接する点が流出(neumann)壁面である
						nflagd = 1;
					}
					else {//Edに接する点が流出壁面以外
						A[0] = 0.0;
					}
				}
			}
			if (Er != -1 && mesh.scond(Er) != 1) {//Erが範囲外でないかつ障害物要素ではない
				Sr = mesh.area(ie, Er);
				A[1] = sr * sr / (2 * Sr);
			}
			else {//Erが範囲外または障害物要素
				if (mesh.ncond(n2) != 0 || mesh.ncond(n3) != 0) {//Erに接する点が内部店ではない
					if (mesh.ncond(n2) == 3 && mesh.ncond(n3) == 3) {//Erに接する点が流出(neumann)壁面である
						nflagr = 1;
					}
					else {//Erに接する点が流出壁面以外
						A[1] = 0.0;
					}
				}
			}
			if (Eu != -1 && mesh.scond(Eu) != 1) {//Euが範囲外でないかつ障害物要素ではない
				Su = mesh.area(ie, Eu);
				A[2] = su * su / (2 * Su);
			}
			else {//Euが範囲外または障害物要素
				if (mesh.ncond(n3) != 0 || mesh.ncond(n4) != 0) {//Euに接する点が内部店ではない
					if (mesh.ncond(n3) == 3 && mesh.ncond(n4) == 3) {//Euに接する点が流出(neumann)壁面である
						nflagu = 1;
					}
					else {//Euに接する点が流出壁面以外
						A[2] = 0.0;
					}
				}
			}
			if (El != -1 && mesh.scond(El) != 1) {//Elが範囲外でないかつ障害物要素ではない
				Sl = mesh.area(ie, El);
				A[3] = sl * sl / (2 * Sl);
			}
			else {//Elが範囲外または障害物要素
				if (mesh.ncond(n4) != 0 || mesh.ncond(n1) != 0) {//Elに接する点が内部店ではない
					if (mesh.ncond(n4) == 3 && mesh.ncond(n1) == 3) {//Elに接する点が流出(neumann)壁面である
						nflagl = 1;
					}
					else {//Euに接する点が流出壁面以外
						A[3] = 0.0;
					}
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

			lambda[ie] = (A[0] + A[1] + A[2] + A[3]) / mesh.Se(ie);
		}
		else {
			lambda[ie] = -1;
		}
		A.clear();
		
	}
}