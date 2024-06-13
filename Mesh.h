#pragma once

#include "Cell.h"

struct Face {
	int nodes[2]; // индексы узлов грани
	int cr; //индекс правой
	int cl; //индекс левой ячейки
	bool is_bound; //граничная ли
	Point f_centr; //координаты центра грани
	double length; //длина грани
	int zone; //номер зоны
};

struct Wall {
	int vel;
	int temp;
	double value;
};

struct Boundary {
	int* rules;
	double* vals;
};

struct zone {
	int grantype;
	Wall* wall;
	Boundary* bnd;
};

class Grid {

private:


public:
	//Переменнные
	int Nx, Ny; //число точек по x y в структурированной сетке
	int nNodes, nCells, nFaces;
	Point* nodes; //координаты узлов

	Face* faces;

	int nZones;
	zone* zones;

	//Функции
	Grid();
	~Grid();

	void readStruct(string filename);

	void createCell(Cell* (&cells));

	void createFaces();

	void cellFunc(Cell* (&cells));

	void setZones();

	void GradCoeffs(Cell* (&cells));
};
