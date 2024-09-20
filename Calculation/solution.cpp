#include"solution.hpp"
#include"output.hpp"
#include"FEM.hpp"
#include"SOR.hpp"
#include"value.hpp"
#include"Mesh.hpp"
#include"param.hpp"
#include"matrix.hpp"
#include<chrono>
#include<vector>
#include<iostream>

using namespace std;

Divergence::Divergence(Mesh2d& Mesh,Boundarycond& BC)//発散量
	:ScalarField2d(Mesh, BC),size(Mesh.nelem())
{

	//発散量コンストラクタ
	scalar.resize(size);
}
void Divergence::cal_divergence(Velocity2d& v) {//発散量の計算
	
	GradientVector C(mesh);

	for (int ie = 0; ie < scalar.size(); ie++) {
		int i1 = mesh.i1(ie);
		int i2 = mesh.i2(ie);
		int i3 = mesh.i3(ie);
		int i4 = mesh.i4(ie);
		double S = mesh.Se(ie);
		//cout << "area[" << ie << "]=" << S << endl;
		if (mesh.scond(ie) != 1) {//非障害物領域
			scalar[ie] = (C[ie][0][0] * v[i1][0] + C[ie][0][1] * v[i1][1]
						+ C[ie][1][0] * v[i2][0] + C[ie][1][1] * v[i2][1]
						+ C[ie][2][0] * v[i3][0] + C[ie][2][1] * v[i3][1]
						+ C[ie][3][0] * v[i4][0] + C[ie][3][1] * v[i4][1]) / S;
		}
		else {//障害物内部
			scalar[ie] = 0.0;
		}
		//cout << "Div[" << ie << "]=" << scalar[ie].v() << endl;
	}
}
void Divergence::cal_divergence(vector<Vector2d>& v) {//発散量の計算
	
	GradientVector C(mesh);

	for (int ie = 0; ie < scalar.size(); ie++) {
		int i1 = mesh.i1(ie);
		int i2 = mesh.i2(ie);
		int i3 = mesh.i3(ie);
		int i4 = mesh.i4(ie);
		double S = mesh.Se(ie);
		if (mesh.scond(ie) != 1) {//非障害物領域
			scalar[ie] = (C[ie][0][0] * v[i1][0] + C[ie][0][1] * v[i1][1]
						+ C[ie][1][0] * v[i2][0] + C[ie][1][1] * v[i2][1]
						+ C[ie][2][0] * v[i3][0] + C[ie][2][1] * v[i3][1]
						+ C[ie][3][0] * v[i4][0] + C[ie][3][1] * v[i4][1]) / S;
		
		}
		else {//障害物内部
			scalar[ie] = 0.0;
		}
		//cout << "Div[" << ie << "]=" << scalar[ie].v() << endl;
		

	}
}
double Divergence::max_div() {//発散量の絶対値の最大値を調べる
	double max = fabs(scalar[0].v());

	for (int i = 0; i < scalar.size(); i++) {
		if (max < fabs(scalar[i].v())) {
			max = fabs(scalar[i].v());
		}
	}
	//cout << max << endl;
	return max;
}
VMPotential::VMPotential(Mesh2d& Mesh, Boundarycond& bc)//速度修正ポテンシャル
	:ScalarField2d(Mesh, bc), size(Mesh.nelem())
{
	//速度修正ポテンシャルコンストラクタ
	scalar.resize(size);
}
void VMPotential::cal_VMP(Divergence& D, SORparam& param) {//速度修正ポテンシャルの計算
	for (int ie = 0; ie < scalar.size(); ie++) {
		if (mesh.scond(ie) != 1) {//障害物領域でない
			scalar[ie] = -D[ie] / param.get_lambda(ie);
		}
		else {//障害物領域
			scalar[ie] = 0.0;
		}
		
	}
}

Predictor::Predictor()
{
	//予測子コンストラクタ
	cout << "object generate: Predictor" << endl;
}
void Predictor::euler_explicit(Velocity2d& v, Pressure& p, Time& T, Mesh2d& Mesh, NDNSparam& Param) {
	//オイラー前進法による予測子計算
	//Fm*(V~-V)/dt= Cv*P - Am*V - Dm*V/Re - dt/2*Bm*V を予測子V~についてとく(Vは既知の流速)
	
	Lumped_Massmatrix Fm(Mesh);//集中化質量行列
	Diffmatrix Dm(Mesh);//拡散行列
	Advecmatrix Am(Mesh, v);//移流行列
	GradientVector Cv(Mesh);//勾配ベクトル
	BTDmatrix Bm(Mesh, v);//BTD行列

	double rRe = 1 / Param.get_Re();//レイノルズ数逆数(拡散項の係数)
	double Cbtd = T.dt() / 2;//BTD項係数

	Vector2d aa;//移流項足し込み変数
	Vector2d dd;//拡散項足し込み変数
	Vector2d bb;//btd項足し込み変数

	vector<double> ff;//集中化質量行列の節点への寄与
	ff.resize(Mesh.nnode());

	vector<Vector2d> UU;
	UU.resize(Mesh.nnode());//速度場の足し込み値(時間進行により増加する値)
	vector<Vector2d> UB;
	UB.resize(Mesh.nnode());//速度場の1step前の値

	vector<int> nx;//境界上の法線フラグx方向
	nx.resize(Mesh.nnode());
	vector<int> ny;//境界上の法線フラグy方向
	ny.resize(Mesh.nnode());
	//領域外側を指す方向の法線
	//0:内部 正:壁面が面する方向が正(x,y) 負:壁面が面する方向が負(x,y)
	//障害物を指す方向の法線
	//0:内部 正:壁面が面する方向が負(x,y) 負:壁面が面する方向が正(x,y)
	for (int ie = 0; ie < Mesh.nelem(); ie++) {

		int i1 = Mesh.i1(ie);
		int i2 = Mesh.i2(ie);
		int i3 = Mesh.i3(ie);
		int i4 = Mesh.i4(ie);

		if (Mesh.scond(ie) != 1) {//障害物要素を除く領域
			nx[i1] += -1;
			ny[i1] += -1;

			nx[i2] += 1;
			ny[i2] += -1;

			nx[i3] += 1;
			ny[i3] += 1;

			nx[i4] += -1;
			ny[i4] += 1;

		}
		else {//障害物要素内
			nx[i1] += 0;
			ny[i1] += 0;

			nx[i2] += 0;
			ny[i2] += 0;

			nx[i3] += 0;
			ny[i3] += 0;

			nx[i4] += 0;
			ny[i4] += 0;
		}
	}



	for (int i = 0; i < Mesh.nnode(); i++) {
		UB[i] = v[i];//前ステップの値の保存

	}
	//速度場の計算
	for (int j = 0; j < node; j++) {
		for (int ie = 0; ie < Mesh.nelem(); ie++) {
			int np = Mesh.nbool1(ie, j);
			int i1 = Mesh.i1(ie);
			int i2 = Mesh.i2(ie);
			int i3 = Mesh.i3(ie);
			int i4 = Mesh.i4(ie);

			ff[np] = ff[np] + Fm[ie][j][j];

			if (Mesh.scond(ie) != 1) {//ie要素が障害物内部でない
				aa = Am[ie][j][0] * v[i1] + Am[ie][j][1] * v[i2] + Am[ie][j][2] * v[i3] + Am[ie][j][3] * v[i4];
				//aa=Am*V^n
				dd = Dm[ie][j][0] * v[i1] + Dm[ie][j][1] * v[i2] + Dm[ie][j][2] * v[i3] + Dm[ie][j][3] * v[i4];
				//dd=Dm*V^n
				bb = Bm[ie][j][0] * v[i1] + Bm[ie][j][1] * v[i2] + Bm[ie][j][2] * v[i3] + Bm[ie][j][3] * v[i4];
				//bb=Bm*V^n
				
				//足し込み
				// UU = Cv*P - Am*V - Dm*V /Re - dt/2 *Bm*V 
				UU[np] += Cv[ie][j] * p[ie] - aa - dd * rRe - Cbtd * bb;
			}
			else {
				Vector2d c(0.0, 0.0);
				UU[np] += c;
			}

			//cout << "UU[" << np << "][0]=" << UU[np][0] << endl;
			//cout << "UU[" << np << "][1]=" << UU[np][1] << endl;
		}
	}
	//予測子の計算
	for (int i = 0; i < Mesh.nnode(); i++) {
		Vector2d c(0.0, 0.0);
		if (Mesh.ncond(i) == 0) {//内部
			//普通に計算
			//予測子=(前ステップ値)+dt*(時間進行による変化)/(集中化質量行列の対角成分)
			//V~= V + dt*UU/ff

			v[i] = UB[i] + T.dt() * UU[i] / ff[i];
		}
		else if (Mesh.ncond(i) == 1) {//剛体内部
			//0固定 前ステップの値を引き継ぐ
			v[i] = c;
		}
		else if (Mesh.ncond(i) == 2) {//流入壁面
			//dirichlet 前ステップの値を引き継ぐ
			v[i] = UB[i];
		}
		else if (Mesh.ncond(i) == 3) {//流出壁面
			//流出方向 neumann 条件　接線方向は既定なし
			v[i] = UB[i] + T.dt() * UU[i] / ff[i];
		}
		else if (Mesh.ncond(i) == 4) {//移動壁面
			//壁面接線(非ゼロ)・法線方向(0) dirichlet 前ステップの値を引き継ぐ
			v[i] = UB[i];
		}
		else if (Mesh.ncond(i) == 5) {//滑りなし壁面
			//壁面接線(0)・法線方向(0) dirichlet 
			v[i] = c;
		}
		else if (Mesh.ncond(i) == 6) {//滑りあり壁面
			//壁面法線方向流速(0) dirichlet
			//壁面接線方向流速勾配(法線方向)=0 neumann
			if (nx[i] != 0 && ny[i] != 0) {//壁面の角
				v[i] = c;
			}
			else if (nx[i] != 0 && ny[i] == 0) {
				v[i][0] = 0.0;
				v[i][1] = UB[i][1] + T.dt() * UU[i][1] / ff[i];
			}
			else if (nx[i] == 0 && ny[i] != 0) {
				v[i][0] = UB[i][0] + T.dt() * UU[i][0] / ff[i];
				v[i][1] = 0.0;
			}
		}
		UU[i] = c;
		//cout << "v[" << i << "][0]=" << v[i][0] << endl;
		//cout << "v[" << i << "][1]=" << v[i][1] << endl;
	}
}

SOR::SOR(SORparam& Param)//速度圧力同時緩和法による修正クラス
	:sparam(Param), nor (0),div_max(0.0)
{
	cout << "object generate: SOR" << endl;
}
void SOR::do_calculation(Velocity2d& v, Pressure& p, Time& T, Mesh2d& Mesh, SORparam& param, Boundarycond& BC) {
	Lumped_Massmatrix Fm(Mesh);//集中化質量行列
	vector<double> ff;//集中化質量行列の対角成分のみの抽出
	ff.resize(Mesh.nnode());
	for (int j = 0; j < node; j++) {
		for (int ie = 0; ie < Mesh.nelem(); ie++) {
			int np = Mesh.nbool1(ie, j);
			ff[np] = ff[np] + Fm[ie][j][j];
		}
	}
	GradientVector Cv(Mesh);
	vector<Vector2d> Utilde;//予測子ステップの値(これを修正していく)
	vector<int> nx;//境界上の法線フラグx方向
	nx.resize(Mesh.nnode());
	vector<int> ny;//境界上の法線フラグy方向
	ny.resize(Mesh.nnode());
	//領域外側を指す方向の法線
	//0:内部 正:壁面が面する方向が正(x,y) 負:壁面が面する方向が負(x,y)
	//障害物を指す方向の法線
	//0:内部 正:壁面が面する方向が負(x,y) 負:壁面が面する方向が正(x,y)
	for (int ie = 0; ie < Mesh.nelem(); ie++) {

		int i1 = Mesh.i1(ie);
		int i2 = Mesh.i2(ie);
		int i3 = Mesh.i3(ie);
		int i4 = Mesh.i4(ie);

		if (Mesh.scond(ie) != 1) {//障害物要素を除く領域
			nx[i1] += -1;
			ny[i1] += -1;

			nx[i2] += 1;
			ny[i2] += -1;

			nx[i3] += 1;
			ny[i3] += 1;

			nx[i4] += -1;
			ny[i4] += 1;

		}
		else {//障害物要素内
			nx[i1] += 0;
			ny[i1] += 0;

			nx[i2] += 0;
			ny[i2] += 0;

			nx[i3] += 0;
			ny[i3] += 0;

			nx[i4] += 0;
			ny[i4] += 0;
		}
	}



	//修正ステップ
	//step1
	Divergence div(Mesh, BC);
	VMPotential phi(Mesh, BC);
	div.cal_divergence(v);
	/*
	for (int i = 0; i < Mesh.nelem(); i++) {
		cout << "div[" << i << "]=" << div[i].v() << endl;
	}
	*/
	
	Utilde.resize(Mesh.nnode());//予測子ステップの値
	for (int i = 0; i < Mesh.nnode(); i++) {
		Utilde[i] = v[i];//予測子ステップの値の代入
	}
	vector<Vector2d> UB;//予測子ステップの値(修正終了まで保存)
	UB.resize(Mesh.nnode());
	for (int i = 0; i < Mesh.nnode(); i++) {
		UB[i] = v[i];//予測子ステップの値の保存
	}
	vector<Vector2d> UU;//流速の修正量
	UU.resize(Mesh.nnode());


	nor = 0;//反復回数初期化
	div_max = 0.0;//発散量最大値(絶対値)(これがepsより小さくなったらループを突破)

	do{
		

		//step2

		for (int ie = 0; ie < Mesh.nelem(); ie++) {
			int i1 = Mesh.i1(ie);
			int i2 = Mesh.i2(ie);
			int i3 = Mesh.i3(ie);
			int i4 = Mesh.i4(ie);
			if (Mesh.scond(ie) != 1) {//非障害物領域
				phi[ie] = -div[ie] / param.get_lambda(ie);//速度修正ポテンシャルの更新
				cout << "lamda[" << ie << "]=" << param.get_lambda(ie) << endl;
				cout << "phi[" << ie << "]=" << phi[ie].v() << endl;
				p[ie] = p[ie] + phi[ie] / T.dt();//圧力の更新
				cout << "p[" << ie << "]=" << p[ie].v() << endl;
			}
			else {//障害物領域
				phi[ie] = 0.0;
				p[ie] = 0.0;
			}
		}
		for (int j = 0; j < node; j++) {
			for (int ie = 0; ie < Mesh.nelem(); ie++) {
				int np = Mesh.nbool1(ie, j);
				if (Mesh.scond(ie) != 1) {//非障害物領域
					UU[np] += Cv[ie][j] * phi[ie];//修正量の計算
				}
				else {//障害物領域
					Vector2d c(0, 0);
					UU[np] += c;
				}
			}
		}
		//速度の修正
		for (int i = 0; i < Mesh.nnode(); i++) {
			Vector2d O(0.0, 0.0);
			if (Mesh.ncond(i) == 0) {//内部
				//普通に計算
				Utilde[i] = Utilde[i] + UU[i] / ff[i];
				//v^(k+1)=v^k + UU/ff
				//UU=Cv*phi^k
			}
			else if (Mesh.ncond(i) == 1) {//剛体内部
				//0固定 
				Utilde[i] = O;
			}
			else if (Mesh.ncond(i) == 2) {//流入壁面
				//dirichlet 前ステップの値を引き継ぐ
				Utilde[i] = UB[i];
			}
			else if (Mesh.ncond(i) == 3) {//流出壁面
				//流出方向 neumann 条件　接線方向は既定なし
				Utilde[i] = Utilde[i] + UU[i] / ff[i];
			}
			else if (Mesh.ncond(i) == 4) {//移動壁面
				//壁面接線(非ゼロ)・法線方向(0) dirichlet 前ステップの値を引き継ぐ
				Utilde[i] = UB[i];
			}
			else if (Mesh.ncond(i) == 5) {//滑りなし壁面
				//壁面接線(0)・法線方向(0) dirichlet 
				Utilde[i] = O;
			}
			else if (Mesh.ncond(i) == 6) {//滑りあり壁面
				//壁面法線方向流速(0) dirichlet
				//壁面接線方向流速勾配(法線方向)=0 neumann
				if (nx[i] != 0 && ny[i] != 0) {//壁面の角
					Utilde[i] = O;
				}
				else if (nx[i] != 0 && ny[i] == 0) {
					Utilde[i][0] = 0.0;
					Utilde[i][1] = Utilde[i][1] + UU[i][1] / ff[i];
				}
				else if (nx[i] == 0 && ny[i] != 0) {
					Utilde[i][0] = Utilde[i][0] + UU[i][0] / ff[i];
					Utilde[i][1] = 0.0;
				}
			}
			UU[i] = O;//修正量の初期化
		}

		for (int i = 0; i < Mesh.nnode(); i++) {
			v[i] = Utilde[i];//流速の更新
		}
		//step3
		div.cal_divergence(v);//発散量の更新

		div_max = div.max_div();//発散量の最大値を計算

		nor++;//反復回数の更新
	} while ((div_max > 1.0e-5) && (nor < 500));
}
int SOR::get_nor() {
	return nor;
}
double SOR::max_div() {
	return div_max;
}
HSMAC_FEM::HSMAC_FEM(Velocity2d& v, Pressure& p, Time& T, Mesh2d& Mesh, NDNSparam& NSP, SORparam& SRP, Boundarycond& bc, InputData& INPUT)
	:V(v), P(p), t(T), mesh(Mesh), nsparam(NSP), sorparam(SRP), BC(bc), input(INPUT), NOR(0), max_div(0.0)
{
	cout << "object generate: HSMAC_FEM" << endl;
}
void HSMAC_FEM::do_solution() {
	Predictor Vp;//流速予測子計算オブジェクト
	SOR S(sorparam);//同時緩和法計算オブジェクト
	OutputData output(mesh, t, BC, V, P, nsparam);//出力用オブジェクト
	output.set_scheme("HSMAC");
	output.set_model(input.get_modelname());
	output.set_path(input.getgoal());
	output.set_Filestage(0);


	//時間進行の開始
	cout << "End of setting up initial condition" << endl;
	cout << "scheme:  HSMAC FEM" << endl;
	cout << "model :" << output.get_model() << "_flow" << endl;
	cout << "Caluculatoin Start !" << endl;
	auto start= std::chrono::high_resolution_clock::now();//総計算時間 開始時刻取得
	cout << "Step No.0   Save Initial Result" << endl;
	cout << "\n";
	cout << "\n";

	output.output_condition();//計算条件の出力
	output.output_time_csv(0);//0ステップの出力

	if (output.get_model() == "cavity") {
		output.output_ghia(0);
	}

	for (int n = 1; n <= t.nend(); n++) {//1stepからスタート
		auto step_start = std::chrono::high_resolution_clock::now();//1step計算時間 開始時刻取得


		Vp.euler_explicit(V, P, t, mesh, nsparam);//オイラー前進法による予測子算出
		S.do_calculation(V, P, t, mesh, sorparam, BC);//同時緩和法による速度圧力修正
		
		NOR = S.get_nor();//同時緩和法反復回数の取得
		max_div = S.max_div();

		
		
		if (n % t.nsample() == 0) {//ステップ数が出力ステップ数のとき
			view_parameters(n);//計算パラメータの表示
			auto step_end = std::chrono::high_resolution_clock::now();    // 1step終了時刻を取得
			auto step_duration = std::chrono::duration_cast<std::chrono::seconds>(step_end - step_start);  // 1step計算時間を秒に変換
			auto step_duration_sum = std::chrono::duration_cast<std::chrono::seconds>(step_end - start);  // 累計計算時間を秒に変換
			cout << "Step No." << n << " :calculation time : " << step_duration.count() << " s" << endl;
			cout << "Cumulative caluculation time: " << step_duration_sum.count() << " s" << endl;
			output.update(V, P);//V,Pの更新(出力のための)
			output.output_time_csv(n);
			if (output.get_model() == "cavity") {
				output.output_ghia(n);
			}

		}
		cout << "\n";
		cout << "\n";

	}
	cout << "Calculation End!" << endl;
	auto end = std::chrono::high_resolution_clock::now();    // 終了時刻を取得
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);  // 総計算時間を秒に変換
	cout << "Total caluculation time: " << duration.count() << " s" << endl;

}
double HSMAC_FEM::Uxmax() {
	double maxUx = V[0][0];//流速最大値
	double absmaxUx = fabs(V[0][0]);//流速最大値絶対値あり
	for (int i = 0; i < mesh.nnode(); i++) {
		if (absmaxUx < fabs(V[i][0])) {
			absmaxUx = fabs(V[i][0]);
			maxUx = V[i][0];
		}
	}
	return maxUx;
}
double HSMAC_FEM::Uymax() {
	double maxUy = V[0][1];//流速最大値
	double absmaxUy = fabs(V[0][1]);//流速最大値絶対値あり
	for (int i = 0; i < mesh.nnode(); i++) {
		if (absmaxUy < fabs(V[i][1])) {
			absmaxUy = fabs(V[i][1]);
			maxUy = V[i][1];
		}
	}
	return maxUy;
}
double HSMAC_FEM::Pmax() {
	double maxP = P[0].v();//圧力最大値
	double absmaxP = fabs(P[0].v());//圧力最大値絶対値あり
	for (int ie = 0; ie < mesh.nelem(); ie++) {
		if (absmaxP < fabs(P[ie].v())) {
			absmaxP = fabs(P[ie].v());
			maxP = P[ie].v();
		}
	}
	return maxP;
}
void HSMAC_FEM::view_parameters(int n) {
	cout << "Step No." << n << "    " << "Time = " << t.ntime(n) << " s" << endl;
	cout << "UMAX: " << Uxmax() << " " << "VMAX: " << Uymax() << " " << "PMAX: " << Pmax() << endl;
	cout << "SOR method Number of iterations = " << NOR << endl;
	cout << "Maximum divergence = " << max_div << endl;

}