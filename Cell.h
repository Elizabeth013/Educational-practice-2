#pragma once

#include "Polygon.h"


struct Vector
{
	double* cx;
};

class Cell {
private:
	Point c; // координаты центра т€жести
	double S; //площадь €чейки
	int nFaces; //число граней
	int nNodes; //число узлов
	int* nodes; //номер узлов



public:

	double Yw; // –ассто€ние от центра €чейки

	double* wk;
	Vector* ck;

	int* faces; //номер граней
	int* fType; //типы граней внутренн€€ и гранична€
	int* cells; //номера соседних €чеек

	Cell();
	~Cell();

	Point get_c() { return c; };

	void set_c(Point c1);

	void set_S(double S1) { S = S1; };
	double get_S() { return S; };

	void set_nFaces(int nf);
	int get_nFaces() { return nFaces; };

	void set_nNodes(int nn);
	int get_nNodes() { return nNodes; };

	int get_Node(int i) { return nodes[i]; };
	int get_Face(int i);
	int get_Cell(int i) { return cells[i]; };

	void set_Face(int iFace, int fIndex) { faces[iFace] = fIndex; };
	void set_fType(int iType, int ftype) { fType[iType] = ftype; };
	void set_cells(int iCell, int cel) { cells[iCell] = cel; };

	void Print(int i);

	void set_Nodes(int* nodes, int nn);
};