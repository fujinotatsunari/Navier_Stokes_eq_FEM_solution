#include<iostream>
#include"input.hpp"
#include"param.hpp"
#include"mesh.hpp"
#include"value.hpp"
#include"solution.hpp"


using namespace std;

int main(void) {
	cout << "2d Navier Stokes Equation FEM Solution" << endl;
	InputData input;
	TimeP tp;
	Time t(tp);
	NDNSparam nsp;
	Mesh2d mesh(input);
	mesh.geninputmesh();
	Velocity2d V(mesh, input.get_BC());
	Pressure P(mesh, input.get_BC());
	V.input(input);
	P.input(input);

	SORparam sorp(mesh);
	
	HSMAC_FEM solve(V, P, t, mesh, nsp, sorp, input.get_BC(), input);
	solve.do_solution();

	return 0;

}

