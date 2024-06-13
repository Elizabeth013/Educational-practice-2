#pragma once

#include <iostream>
#include <string>
#include <fstream>

using namespace std;

struct Point {
	double x, y;
};

class Polygon {
private:
	Point* p;
	int n;
public:
	Polygon();
	Polygon(Point* p1, int n1);
	~Polygon();

	void DataEntry(string filename);

	int get_n() { return n; };
	Point get_p(int i) { return p[i]; };

	Point center();
	double Square();

};