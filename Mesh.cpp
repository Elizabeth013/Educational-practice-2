#include "Mesh.h"

Grid::Grid()
{
	Nx = 0;
	Ny = 0;
	nNodes = 0;
	nCells = 0;
	nFaces = 0;

	nodes = new  Point[nNodes];
}

Grid::~Grid()
{
}

void Grid::readStruct(string filename)
{
	int tmp;
	ifstream reading;
	reading.open(filename);

	if (reading.is_open()) {
		reading >> Nx >> Ny >> tmp >> tmp;

		nNodes = Nx * Ny;
		nCells = (Nx - 1) * (Ny - 1);
		nFaces = Nx * (Ny - 1) + (Nx - 1) * Ny;

		nodes = new Point[nNodes];

		int count = 0;
		for (int i = 0; i < Nx; i++) {
			for (int j = 0; j < Ny; j++) {
				reading >> tmp >> tmp >> nodes[count].x >> nodes[count].y;
				count++;
			}
		}
		cout << "Work " << filename << endl;
	}
	else cout << "Error  " << filename << endl;

}

void Grid::createCell(Cell* (&cells))
{
	int count = 0;
	for (int i = 0; i < Nx - 1; i++) {
		for (int j = 0; j < Ny - 1; j++) {
			int nc = (Ny - 1) * i + j;
			cout << "Номер ячейки " << nc << endl;
			int nNodes = 4;

			int* nodes = new int[nNodes];
			// element 0
			int n1 = Ny * i + j;
			nodes[0] = n1;
			// element 1
			int n2 = Ny * (i + 1) + j;
			nodes[1] = n2;
			// element 2
			int n3 = Ny * (i + 1) + (j + 1);
			nodes[2] = n3;
			// element 3
			int n4 = Ny * i + (j + 1);
			nodes[3] = n4;

			cells[nc].set_Nodes(nodes, nNodes);
		}
	}
}

void Grid::createFaces()
{
	faces = new Face[nFaces];

	int count = 0;

	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny - 1; j++) {
			int n1 = Ny * i + j;
			int n2 = Ny * i + j + 1;
			faces[count].nodes[0] = n1;
			faces[count].nodes[1] = n2;
			if (i == 0) {
				faces[count].is_bound = true;
				faces[count].cr = -1;
				faces[count].cl = (Ny - 1) * i + j;
				faces[count].zone = 0;
			}
			else if (i == Nx - 1) {
				faces[count].is_bound = true;
				faces[count].cr = (Ny - 1) * (i - 1) + j;
				faces[count].cl = -1;
				faces[count].zone = 2;
			}
			else {
				faces[count].is_bound = false;
				faces[count].cr = (Ny - 1) * (i - 1) + j;
				faces[count].cl = (Ny - 1) * i + j;
				faces[count].zone = -1;
			}

			faces[count].f_centr.x = 0.5 * (nodes[n1].x + nodes[n2].x);
			faces[count].f_centr.y = 0.5 * (nodes[n1].y + nodes[n2].y);

			double dx = (nodes[n1].x - nodes[n2].x);
			double dy = (nodes[n1].y - nodes[n2].y);
			faces[count].length = sqrt(dx * dx + dy * dy);

			count++;
		}
	}

	for (int j = 0; j < Ny; j++) {
		for (int i = 0; i < Nx - 1; i++) {
			int n1 = Ny * i + j;
			int n2 = Ny * (i + 1) + j;
			faces[count].nodes[0] = n1;
			faces[count].nodes[1] = n2;
			if (j == 0) {
				faces[count].is_bound = true;
				faces[count].cr = (Ny - 1) * i + j;
				faces[count].cl = -1;
				faces[count].zone = 1;
			}
			else if (j == Ny - 1) {
				faces[count].is_bound = true;
				faces[count].cr = -1;
				faces[count].cl = (Ny - 1) * i + j - 1;
				faces[count].zone = 3;
			}
			else {
				faces[count].is_bound = false;
				faces[count].cr = (Ny - 1) * i + j;
				faces[count].cl = (Ny - 1) * i + j - 1;
				faces[count].zone = -1;
			}

			faces[count].f_centr.x = 0.5 * (nodes[n1].x + nodes[n2].x);
			faces[count].f_centr.y = 0.5 * (nodes[n1].y + nodes[n2].y);

			double dx = (nodes[n1].x - nodes[n2].x);
			double dy = (nodes[n1].y - nodes[n2].y);
			faces[count].length = sqrt(dx * dx + dy * dy);

			count++;
		}
	}
}

void Grid::cellFunc(Cell* (&cells))
{
	for (int i = 0; i < nCells; i++) {
		int mNodes = cells[i].get_nNodes();

		int* nds = new int[mNodes];
		Point* pnts = new Point[mNodes];

		for (int j = 0; j < mNodes; j++) {
			int n1 = cells[i].get_Node(j);
			nds[j] = n1;

			pnts[j].x = nodes[n1].x;
			pnts[j].y = nodes[n1].y;
		}

		Polygon p1(pnts, mNodes);

		Point center = p1.center();

		double Sq = p1.Square();

		cells[i].set_c(center);
		cells[i].set_S(Sq);

		cells[i].set_nFaces(mNodes);
	}

	for (int k = 0; k < nFaces; k++) {
		int n1 = faces[k].nodes[0];
		int n2 = faces[k].nodes[1];
		int cl = faces[k].cl;
		int cr = faces[k].cr;
		int ftype = 0;
		if (faces[k].is_bound) ftype = 1;

		if (cr >= 0) {
			int nn = cells[cr].get_nNodes();
			int* nncl = new int[nn];

			for (int i = 0; i < nn; i++)
				nncl[i] = cells[cr].get_Node(i);

			int iFace;

			for (int i = 0; i < nn; i++) {
				int j = i + 1;
				if (j == nn) j = 0;
				if (nncl[i] == n1 && nncl[j] == n2) {
					iFace = i;
					cells[cr].set_Face(iFace, k);
					cells[cr].set_fType(iFace, ftype);
					cells[cr].set_cells(iFace, cl);
				}
			}
		}
		if (cl >= 0) {
			int nn = cells[cl].get_nNodes();
			int* nncl = new int[nn];

			for (int i = 0; i < nn; i++)
				nncl[i] = cells[cr].get_Node(i);

			int iFace;

			for (int i = 0; i < nn; i++) {
				int j = i + 1;
				if (j == nn) j = 0;
				if (nncl[i] == n2 && nncl[j] == n1) {
					iFace = i;
					cells[cl].set_Face(iFace, k);
					cells[cl].set_fType(iFace, ftype);
					cells[cl].set_cells(iFace, cr);
				}
			}
		}
	}
}

void Grid::setZones()
{
	nZones = 4;
	zones = new zone[nZones];
}

void Grid::GradCoeffs(Cell* (&cells))
{
	for (int i = 0; i < nCells; i++) {
		int nFaces = cells[i].get_nFaces();

		cells[i].wk = new double[nFaces];
		cells[i].ck = new Vector[nFaces];

		int iDim = 2;

		for (int k = 0; k < nFaces; k++)
			cells[i].ck[k].cx = new double[iDim];

		Point xc = cells[i].get_c();

		double* dx = new double[nFaces];
		double* dy = new double[nFaces];

		double axx = 0., axy = 0., ayy = 0.;

		for (int k = 0; k, nFaces; k++) {
			int nf = cells[i].get_Face(k);
			if (!faces[nf].is_bound) {

				int nc = cells[i].get_Cell(k); //
				Point xk = cells[nc].get_c();

				dx[k] = xk.x - xc.x;
				dy[k] = xk.y - xc.y;
			}
			else {
				Point xk = faces[nf].f_centr;
				dx[k] = xk.x - xc.x;
				dy[k] = xk.y - xc.y;
			}

			double wk = 1. / sqrt(dx[k] * dx[k] + dy[k] * dy[k]);
			cells[i].wk[k] = wk;

			axx += wk * dx[k] * dx[k];
			axy += wk * dx[k] * dy[k];
			ayy += wk * dy[k] * dy[k];
		}

		double det = axx * ayy - axy * axy;
		double Mxx, Mxy, Myy;

		Mxx = ayy / det;
		Mxy = -axy / det;
		Myy = axx / det;

		for (int k = 0; k < nFaces; k++) {
			cells[i].ck[k].cx[0] = Mxx * dx[k] + Mxy * dy[k];
			cells[i].ck[k].cx[1] = Mxy * dx[k] + Myy * dy[k];
		}
	}
}
