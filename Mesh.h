#pragma once

#include "Cell.h"

struct Face {
	int nodes[2]; // ������� ����� �����
	int cr; //������ ������
	int cl; //������ ����� ������
	bool is_bound; //��������� ��
	Point f_centr; //���������� ������ �����
	double length; //����� �����
	int zone; //����� ����
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
	//�����������
	int Nx, Ny; //����� ����� �� x y � ����������������� �����
	int nNodes, nCells, nFaces;
	Point* nodes; //���������� �����

	Face* faces;

	int nZones;
	zone* zones;

	//�������
	Grid();
	~Grid();

	void readStruct(string filename);

	void createCell(Cell* (&cells));

	void createFaces();

	void cellFunc(Cell* (&cells));

	void setZones();

	void GradCoeffs(Cell* (&cells));
};
