#include "Cell.h"

Cell::Cell()
{
	c.x = 0.;// координаты центра т€жести
	c.y = 0.;
	S = 0; //площадь €чейки

	nFaces = 0; //число граней
	nNodes = 0; //число узлов

	faces = new int[nFaces]; //номер граней
	nodes = new int[nNodes]; //номер узлов

	fType = new int[nFaces]; //типы граней внутренн€€ и гранична€
	cells = new int[nFaces]; //номера соседних €чеек
}

Cell::~Cell()
{
}

void Cell::set_c(Point c1)
{
	c.x = c1.x;
	c.y = c1.y;
}

void Cell::set_nFaces(int nf)
{
	nFaces = nf;
	faces = new int[nFaces];
	fType = new int[nFaces];
	cells = new int[nFaces];
}

void Cell::set_nNodes(int nn)
{
	nNodes = nn;
	nodes = new int[nNodes];
}

int Cell::get_Face(int i)
{
	if (i < nFaces)
		return faces[i];
	else
	{
		cout << "Error NO FACES" << endl;
		return -1;
	}

}

void Cell::Print(int m)
{
	string f = "cells.txt";
	ofstream record(f, ios::out | ios::app);

	if (record) {
		record << "»ндекс €чейки " << m << endl;
		record << "центр " << c.x << ", " << c.y << endl;
		record << "площадь " << S << endl;
		record << "номер грани " << nFaces << endl;
		record << "индексы грани и ее тип" << endl;
		for (int i = 0; i < nFaces; i++) {
			record << "ind= " << faces[i] << " type= " << fType[i] << endl;
		}
		record << "номер узла " << nNodes << endl;
		record << "индекс узла " << endl;
		for (int i = 0; i < nNodes; i++) {
			record << nodes[i];
			if (i < nNodes - 1) record << ", ";
		}
		record << endl;

		record << "соседние индексы " << endl;
		for (int i = 0; i < nFaces; i++) {
			record << cells[i];
			if (i < nFaces - 1) record << ", ";
		}
		record << endl;

		record << "-----------------" << endl;
	}
	else cout << "Ќе получилось открыть файл!" << f << endl;

	record.close();
}

void Cell::set_Nodes(int* nodes, int nn)
{
	nNodes = nn;
	nodes = new int[nNodes]; //Ќомер узлов окружающих €чейку

}
