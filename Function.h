#pragma once

#include "Mesh.h"
#include "Polygon.h"
#include <iomanip>


struct parameters {

	double ro, p, h, H, u, v, w, T, E, e;

	double mu, la, Pr;

	double Cp, Gam, Gm;

	double* U, * U1;
	double* V;

};

struct changes {
	double* dU;
};

struct Gradient {
	Vector* g;
};

void Init(parameters* (&p), int nCells, int Nm);
void viscid(parameters* p, changes* (&du), Grid mesh, Cell* cells, double dt, Gradient* gr, int Nm);
void getParametrs(parameters* (&p), int nCells, int Nm);
void convect(parameters* p, changes* (&du), Grid mesh, Cell* cells, int it, double dt);
void centerToWall(Grid mesh, Cell* (&cells), int nCells);
double Dist(Point A, Point B, Point E);
void setGran(Grid& mesh);
void Gradients(Cell* cells, Grid mesh, Gradient* (&gr), parameters* p, int Nm);
void matrixDiag(double A[4][4], double L[4], double B[4][4]);
void matrixMatrix(double A[4][4], double B[4][4], double C[4][4]);
void matrixVector(double A[4][4], double B[4], double C[4]);
void printMatrix(double A[4][4]);
void sMatr(double A[4][4], double u, double v, double ro, double p, double h, double nx, double ny, int iMod);
void ConvectNS(parameters* p, changes* (&du), Grid mesh, Cell* cells, double dt, int Nm);