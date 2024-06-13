#include "Polygon.h"


Polygon::Polygon()
{
	n = 0;
	p = new Point[n];
}

Polygon::Polygon(Point* p1, int n1)
{
	n = n1;
	p = new Point[n];

	for (int i = 0; i < n1; i++)  p[i] = p1[i];

}

Polygon::~Polygon()
{
}

void Polygon::DataEntry(string filename)
{
	int n1;
	Point* p1;

	ifstream reading;
	reading.open(filename);

	if (reading.is_open()) {
		reading >> n;
		p = new Point[n];
		for (int i = 0; i < n; i++) reading >> p[i].x >> p[i].y;
		cout << "Work " << filename << endl;
	}
	else cout << "Error  " << filename << endl;


	reading.close();
}

Point Polygon::center()
{
	int m = n + 1;
	Point* xn = new Point[m];

	for (int i = 0; i < n; i++)
	{
		xn[i].x = p[i].x;
		xn[i].y = p[i].y;
	}
	xn[n].x = xn[0].x;
	xn[n].y = xn[0].y;

	double A = 0;
	double Cx = 0, Cy = 0;

	for (int i = 0; i < n; i++)
	{
		A = A + xn[i].x * xn[i + 1].y - xn[i + 1].x * xn[i].y;
		Cx = Cx + (xn[i].x + xn[i + 1].x) * (xn[i].x * xn[i + 1].y - xn[i + 1].x * xn[i].y);
		Cy = Cy + (xn[i].y + xn[i + 1].y) * (xn[i].x * xn[i + 1].y - xn[i + 1].x * xn[i].y);
	}

	A = 0.5 * A;
	Point p1;
	p1.x = Cx / (A * 6.);
	p1.y = Cy / (A * 6.);


	return p1;
}

double Polygon::Square()
{
	int m = n + 1;
	Point* xn = new Point[m];

	for (int i = 0; i < n; i++)
	{
		xn[i].x = p[i].x;
		xn[i].y = p[i].y;
	}
	xn[n].x = xn[0].x;
	xn[n].y = xn[0].y;

	double A = 0;

	for (int i = 0; i < n; i++)
		A = A + xn[i].x * xn[i + 1].y - xn[i + 1].x * xn[i].y;


	A = 0.5 * A;
	A = fabs(A);

	return A;
}
