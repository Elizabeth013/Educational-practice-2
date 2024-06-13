#include "Cell.h"

Cell::Cell()
{
	c.x = 0.;// ���������� ������ �������
	c.y = 0.;
	S = 0; //������� ������

	nFaces = 0; //����� ������
	nNodes = 0; //����� �����

	faces = new int[nFaces]; //����� ������
	nodes = new int[nNodes]; //����� �����

	fType = new int[nFaces]; //���� ������ ���������� � ���������
	cells = new int[nFaces]; //������ �������� �����
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
		record << "������ ������ " << m << endl;
		record << "����� " << c.x << ", " << c.y << endl;
		record << "������� " << S << endl;
		record << "����� ����� " << nFaces << endl;
		record << "������� ����� � �� ���" << endl;
		for (int i = 0; i < nFaces; i++) {
			record << "ind= " << faces[i] << " type= " << fType[i] << endl;
		}
		record << "����� ���� " << nNodes << endl;
		record << "������ ���� " << endl;
		for (int i = 0; i < nNodes; i++) {
			record << nodes[i];
			if (i < nNodes - 1) record << ", ";
		}
		record << endl;

		record << "�������� ������� " << endl;
		for (int i = 0; i < nFaces; i++) {
			record << cells[i];
			if (i < nFaces - 1) record << ", ";
		}
		record << endl;

		record << "-----------------" << endl;
	}
	else cout << "�� ���������� ������� ����!" << f << endl;

	record.close();
}

void Cell::set_Nodes(int* nodes, int nn)
{
	nNodes = nn;
	nodes = new int[nNodes]; //����� ����� ���������� ������

}
