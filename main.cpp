#include "Polygon.h"
#include "Cell.h"
#include "Mesh.h"
#include "Function.h"

#include <chrono>
#include <ctime>

using namespace std;

int main()
{
	setlocale(LC_ALL, "rus");

	Grid mesh;
	string filename = "grid.txt";

	mesh.readStruct(filename);


	int nCells = mesh.nCells;

	Cell* cells = new Cell[nCells];

	mesh.createCell(cells);
	mesh.createFaces();

	mesh.cellFunc(cells);

	remove("cells.txt");

	parameters* p = new parameters[nCells];
	changes* du = new changes[nCells];

	int Nm = 4; //NS

	Init(p, nCells, Nm);

	mesh.setZones();

	setGran(mesh);

	centerToWall(mesh, cells, nCells);

	mesh.GradCoeffs(cells);

	int itMax = 1000;
	int it = 0;

	double resMin = 0.000001;

	double res = 1;

	for (int i = 0; i < nCells; i++) {
		du[i].dU = new double[Nm];
	}
	// ¬ыделение пам€ти под градиенты параметров
	Gradient* gr = new Gradient[nCells];
	for (int i = 0; i < nCells; i++) {
		gr[i].g = new Vector[Nm];
		for (int j = 0; j < Nm; j++)
			gr[i].g[j].cx = new double[2];
	}

	auto start = chrono::system_clock::now();




	while (it<itMax && res>resMin) {
		it++;
		for (int i = 0; i < nCells; i++) {
			for (int j = 0; j < Nm; j++) {
				p[i].U[j] = p[i].U1[j]; // переприсвоение
				du[i].dU[j] = 0;
			}
		}

		Gradients(cells, mesh, gr, p, Nm); // NS

		double dt = 0.000001;

		//ѕриращение за счет нев€зких потоков
		//convect(p, du, mesh, cells,it, dt); // дл€ модельного уравнени€
		ConvectNS(p, du, mesh, cells, dt, Nm);

		//ѕриращение за счет в€зкости
		viscid(p, du, mesh, cells, dt, gr, Nm);


		for (int i = 0; i < nCells; i++) {
			for (int j = 0; j < Nm; j++) {
				p[i].U1[j] = p[i].U[j] + du[i].dU[j];
			}
		}

		getParametrs(p, nCells, Nm);

		res = 0.;
		for (int i = 0; i < nCells; i++) {
			double res1 = abs(du[i].dU[3] / p[i].U1[3]);
			if (res < res1) res = res1;
		}
		cout << "It= " << it << ", res= " << res << endl;
	}

	auto end = chrono::system_clock::now();

	chrono::duration<double> elapsed_seconds = end - start;
	time_t end_time = chrono::system_clock::to_time_t(end);

	cout << " " << elapsed_seconds.count() << endl;
	int Nx = mesh.Nx;
	int Ny = mesh.Ny;

	return 0;
}