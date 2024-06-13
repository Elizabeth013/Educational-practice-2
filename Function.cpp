#include "Function.h"

void Init(parameters* (&p), int nCells, int Nm)
{
	double TO = 293, Cp, la = 2.5658e-2, PO, Gm, Gam, R, ro, mu = 17.863e-6;
	double UO;

	//


	UO = 10.0;//

	// Давление 
	PO = 101325.0;
	Gam = 1.4;
	Gm = 28.97;

	R = 8314.41 / Gm;
	ro = PO / (R * TO);
	Cp = Gam / (Gam - 1.) * R;

	

	for (int i = 0; i < nCells; i++) {
		p[i].ro = ro;
		p[i].p = PO;
		p[i].u = UO;
		p[i].v = 0.;
		p[i].w = 0.;
		p[i].T = TO;

		p[i].Cp = Cp;
		p[i].la = la;
		p[i].mu = mu;
		p[i].Gam = Gam;
		p[i].Gm = Gm;

		p[i].Pr = p[i].mu * p[i].Cp / p[i].la;

		p[i].h = Cp * TO;
		double q = 0.5 * (p[i].u * p[i].u + p[i].v * p[i].v);

		p[i].H = p[i].h + q;
		p[i].E = p[i].H - p[i].p / p[i].ro;

		p[i].e = p[i].h / p[i].Gam;

		// NS

		p[i].U = new double[Nm];
		p[i].U1 = new double[Nm];

		p[i].V = new double[Nm];

		p[i].U1[0] = p[i].ro;
		p[i].U1[1] = p[i].ro * p[i].u;
		p[i].U1[2] = p[i].ro * p[i].v;
		p[i].U1[3] = p[i].ro * p[i].E;

		p[i].V[0] = p[i].ro;
		p[i].V[1] = p[i].u;
		p[i].V[2] = p[i].v;
		p[i].V[3] = p[i].h;
	}
}

void viscid(parameters* p, changes* (&du), Grid mesh, Cell* cells, double dt, Gradient* gr, int Nm)
{
	int nFaces = mesh.nFaces;

	double* Fv = new double[Nm];

	for (int i = 0; i < nFaces; i++)
	{
		Face face = mesh.faces[i];

		int cr = face.cr;
		int cl = face.cl;

		double length = face.length;

		int n1 = face.nodes[0];
		int n2 = face.nodes[1];

		Point x1 = mesh.nodes[n1];
		Point x2 = mesh.nodes[n2];

		double nx = -(x2.y - x1.y) / length;
		double ny = (x2.x - x1.x) / length;

		if (face.is_bound) {
			int c = max(cr, cl);
			// Координаты центра тяжести ячейки
			Point xc = cells[c].get_c();

			int z = face.zone;
			int grantype = mesh.zones[z].grantype;

			if (grantype == 1) { // Проскальзывание
				// Расстояние до стенки
				double dl = cells[c].Yw;

				double du_dx, du_dy, dv_dx, dv_dy, dh_dx, dh_dy, u_, v_;

				if (mesh.zones[z].bnd[0].rules[0] == 0) {
					// Градиенты u и v 
					double ux = gr[c].g[1].cx[0];
					double uy = gr[c].g[1].cx[1];

					du_dx = ny * ny * ux - nx * ny * uy;
					du_dy = -nx * ny * ux + nx * nx * uy;

					double vx = gr[c].g[2].cx[0];
					double vy = gr[c].g[2].cx[1];

					dv_dx = ny * ny * vx - nx * ny * vy;
					dv_dy = -nx * ny * vx + nx * nx * vy;

					u_ = p[c].u;
					v_ = p[c].v;
				}
				if (mesh.zones[z].bnd[0].rules[0] == 1) { // Без проскальзывания

					double du_dn = p[c].u / dl;
					double dv_dn = p[c].v / dl;
					// Производные скорости вдоль стенки (перпендикулярны вектору нормали)
					double du_dl = 0.;
					double dv_dl = 0.;

					du_dx = du_dn * nx - du_dl * ny;
					dv_dx = dv_dn * nx - dv_dl * ny;

					du_dy = du_dn * ny + du_dl * nx;
					dv_dy = dv_dn * ny + dv_dl * nx;

					u_ = 0.;
					v_ = 0.;
				}
				// Энергия
				if (mesh.zones[z].bnd[0].rules[1] == 1) { // Tw
					double Tw = mesh.zones[z].bnd[0].vals[0];
					double hw = p[c].Cp * Tw;
					double dh_dn = (p[c].h - hw) / dl;
					double dh_dl = 0.;
					dh_dx = dh_dn * nx - dh_dl * ny;
					dh_dy = dh_dn * ny + dh_dl * nx;
				}
				if (mesh.zones[z].bnd[0].rules[1] == 2) { // qw=0
					double tx = gr[c].g[3].cx[0];
					double ty = gr[c].g[3].cx[1];

					dh_dx = ny * ny * tx - nx * ny * ty;
					dh_dy = -nx * ny * tx + nx * nx * ty;
				}

				double mu = p[c].mu;
				
				double div = du_dx * dv_dy;
				double txx = mu * (2. * du_dx - 2. / 3. * div);
				double tyy = mu * (2. * dv_dy - 2. / 3. * div);
				double txy = mu * (du_dy + dv_dx);

				double mu_Pr = p[c].mu / p[c].Pr;
				double qx = -mu_Pr * dh_dx;
				double qy = -mu_Pr * dh_dy;

				Fv[0] = 0.;
				Fv[1] = txx * nx + txy * ny;
				Fv[2] = txy * nx + tyy * ny;
				Fv[3] = (u_ * txx + v_ * txy - qx) * nx + (u_ * txy + v_ * tyy - qy) * ny;

			}
			if (grantype == 3) {
				double ux = gr[c].g[1].cx[0];
				double uy = gr[c].g[1].cx[1];

				double du_dx = ny * ny * ux - nx * ny * uy;
				double du_dy = -nx * ny * ux + nx * nx * uy;

				double vx = gr[c].g[2].cx[0];
				double vy = gr[c].g[2].cx[1];

				double dv_dx = ny * ny * vx - nx * ny * vy;
				double dv_dy = -nx * ny * vx + nx * nx * vy;

				double u_ = p[c].u;
				double v_ = p[c].v;

				double tx = gr[c].g[3].cx[0];
				double ty = gr[c].g[3].cx[1];

				double dh_dx = ny * ny * tx - nx * ny * ty;
				double dh_dy = -nx * ny * tx + nx * nx * ty;

				double mu = p[c].mu;
				double div = du_dx * dv_dy;
				double txx = mu * (2. * du_dx - 2. / 3. * div);
				double tyy = mu * (2. * dv_dy - 2. / 3. * div);
				double txy = mu * (du_dy + dv_dx);

				double mu_Pr = p[c].mu / p[c].Pr;
				double qx = -mu_Pr * dh_dx;
				double qy = -mu_Pr * dh_dy;

				Fv[0] = 0.;
				Fv[1] = txx * nx + txy * ny;
				Fv[2] = txy * nx + tyy * ny;
				Fv[3] = (u_ * txx + v_ * txy - qx) * nx + (u_ * txy + v_ * tyy - qy) * ny;

			}
			if (grantype == 2 || grantype == 4) {
				for (int m = 0; m < Nm; m++) Fv[m] = 0;
			}
			double length = face.length;
			if (cr >= 0) {
				double Sr = cells[cr].get_S();
				for (int m = 0; m < Nm; m++) du[cr].dU[m] += -Fv[m] * length / Sr * dt;
			}
			if (cl >= 0) {
				double Sl = cells[cr].get_S();
				for (int m = 0; m < Nm; m++) du[cl].dU[m] += Fv[m] * length / Sl * dt;
			}

		}
		else {
			// Координаты правой и левой ячейки
			Point xr = cells[cr].get_c();
			Point xl = cells[cl].get_c();
			// Расстояние между ячейками
			double dx = xr.x - xl.x;
			double dy = xr.y - xl.y;
			double dl = sqrt(dx * dx + dy * dy);

			double du_dx = 0.5 * (gr[cr].g[1].cx[0] + gr[cl].g[1].cx[0]);
			double du_dy = 0.5 * (gr[cr].g[1].cx[1] + gr[cl].g[1].cx[1]);

			double dv_dx = 0.5 * (gr[cr].g[2].cx[0] + gr[cl].g[2].cx[0]);
			double dv_dy = 0.5 * (gr[cr].g[2].cx[1] + gr[cl].g[2].cx[1]);

			double dh_dx = 0.5 * (gr[cr].g[3].cx[0] + gr[cl].g[3].cx[0]);
			double dh_dy = 0.5 * (gr[cr].g[3].cx[1] + gr[cl].g[3].cx[1]);

			double mu = 0.5 * (p[cr].mu + p[cl].mu);
			double div = du_dx + dv_dy;
			double txx = mu * (2. * du_dx - 2. / 3. * div);
			double tyy = mu * (2. * dv_dy - 2. / 3. * div);
			double txy = mu * (du_dy + dv_dx);

			double mu_Pr = 0.5 * (p[cr].mu / p[cr].Pr + p[cl].mu / p[cl].Pr);

			double qx = -mu_Pr * dh_dx;
			double qy = -mu_Pr * dh_dy;

			double u_ = 0.5 * (p[cr].u + p[cl].u);
			double v_ = 0.5 * (p[cr].v + p[cl].v);

			Fv[0] = 0.;
			Fv[1] = txx * nx + txy * ny;
			Fv[2] = txy * nx + tyy * ny;
			Fv[3] = (u_ * txx + v_ * txy - qx) * nx + (u_ * txy + v_ * tyy - qy) * ny;


			// Поток через грань
			//double Fv = mu_Pr * dh_dn;

			double length = face.length; // Длина грани
			double Sr = cells[cr].get_S(); // Площадь правой ячейки
			double Sl = cells[cl].get_S(); // Площадь левой ячейки

			for (int m = 0; m < Nm; m++) {
				du[cr].dU[m] += -Fv[m] * length / Sr * dt;
				du[cl].dU[m] += Fv[m] * length / Sl * dt;
			}
		}
	}
}

void getParametrs(parameters* (&p), int nCells, int Nm)
{
	for (int i = 0; i < nCells; i++) {

		p[i].ro = p[i].U1[0];
		p[i].u = p[i].U1[1] / p[i].ro;
		p[i].v = p[i].U1[2] / p[i].ro;
		p[i].E = p[i].U1[3] / p[i].ro;

		double q = 0.5 * (p[i].u * p[i].u + p[i].v * p[i].v);

		p[i].e = p[i].E - q;
		p[i].h = p[i].e * p[i].Gam; // Справедливо когда пост адиабата
		p[i].T = p[i].h / p[i].Cp;
		p[i].H = p[i].h + q;

		double R = 8314.41 / p[i].Gm;
		p[i].p = p[i].ro * R * p[i].T;

		p[i].V[0] = p[i].ro;
		p[i].V[1] = p[i].u;
		p[i].V[2] = p[i].v;
		p[i].V[3] = p[i].h;
	}
}

void convect(parameters* p, changes* (&du), Grid mesh, Cell* cells, int it, double dt)
{
	int nFaces = mesh.nFaces;

	for (int i = 0; i < nFaces; i++) {
		Face face = mesh.faces[i];

		int cr = face.cr;
		int cl = face.cl;

		int n1 = face.nodes[0];
		int n2 = face.nodes[1];

		double dl = face.length;

		Point x1 = mesh.nodes[n1];
		Point x2 = mesh.nodes[n2];

		double nx = -(x2.y - x1.y) / dl;
		double ny = -(x2.x - x1.x) / dl;

		double Fc;

		if (face.is_bound) {
			int c = max(cr, cl);
			double S = cells[c].get_S();

			if (face.zone == 0) {
				double uInlet, vInlet, tInlet;
				uInlet = p[c].u;
				vInlet = p[c].v;
				tInlet = 800.0;

				double un = uInlet * nx + vInlet * ny;

				double h = p[c].Cp * tInlet;
				double H = h + 0.5 * (uInlet * uInlet + vInlet * vInlet);
				double E = h / p[c].Gam + 0.5 * (uInlet * uInlet + vInlet * vInlet);

				Fc = p[c].ro * H * un;
				du[c].dU[0] += -Fc * dl / S * dt;
			}
			if (face.zone == 2) {
				double un = p[c].u * nx + p[c].v * ny;
				Fc = p[c].ro * p[c].H * un;
				du[c].dU[0] += Fc * dl / S * dt;
			}
			if (face.zone == 1 || face.zone == 3) {
				du[c].dU[0] += 0;
			}
		}
		else {
			double u1, v1, un1, H1, E1, ro1;

			u1 = 0.5 * (p[cr].u + p[cl].u);
			v1 = 0.5 * (p[cr].v + p[cl].v);
			un1 = u1 * nx + v1 + ny;
			H1 = 0.5 * (p[cr].H + p[cl].H);
			E1 = 0.5 * (p[cr].E + p[cl].E);
			ro1 = 0.5 * (p[cr].ro + p[cl].ro);

			double A = H1 / E1 * un1;
			double Apl, Amn;
			Apl = 0.5 * (A + abs(A));
			Amn = 0.5 * (A - abs(A));

			double UL, UR;
			UL = p[cl].U[0];
			UR = p[cr].U[0];

			Fc = Apl * UL + Amn * UR;

			double Sr = cells[cr].get_S();
			double Sl = cells[cl].get_S();

			du[cr].dU[0] += Fc * dl / Sr * dt;
			du[cl].dU[0] += -Fc * dl / Sl * dt;
		}
	}
}

void centerToWall(Grid mesh, Cell* (&cells), int nCells)
{
	int nFaces = mesh.nFaces;

	for (int k = 0; k < nFaces; k++) {
		double z1 = 1.e10;
		Point E = cells[k].get_c();

		for (int i = 0; i < nFaces; i++) {

			int nz = mesh.faces[i].zone;
			int grantype = mesh.zones[nz].grantype;

			if (grantype == 1) {
				int n1 = mesh.faces[i].nodes[0];
				int n2 = mesh.faces[i].nodes[1];
				Point A = mesh.nodes[n1];
				Point B = mesh.nodes[n2];

				double z2 = Dist(A, B, E);

				if (z1 > z2) z1 = z2;
			}
		}
		cells[k].Yw = z1;
	}
}

double Dist(Point A, Point B, Point E)
{
	double AB[2], BE[2], AE[2];;

	double AB_BE, AB_AE, x, y;
	double x1, y1, x2, y2, mod;

	double dist;

	AB[0] = B.x - A.x;
	AB[1] = B.y - A.y;

	BE[0] = E.x - B.x;
	BE[1] = E.y - B.y;

	AE[0] = E.x - A.x;
	AE[1] = E.y - A.y;

	AB_BE = (AB[0] * BE[0] + AB[1] * BE[1]);
	AB_AE = (AB[0] * AE[0] + AB[1] * AE[1]);

	if (AB_BE > 0) {
		x = E.x - B.x;
		y = E.y - B.y;
		dist = sqrt(x * x + y * y);
	}
	else if (AB_AE < 0) {
		x = E.x - A.x;
		y = E.y - A.y;
		dist = sqrt(x * x + y * y);
	}
	else {
		x1 = AB[0];
		y1 = AB[1];
		x2 = AE[0];
		y2 = AE[1];
		mod = sqrt(x1 * x1 + y1 * y1);
		dist = abs(x1 * y2 - y1 * x2) / mod;
	}

	return dist;
}

void setGran(Grid& mesh)
{
	// NS
	mesh.zones[0].grantype = 3;
	mesh.zones[1].grantype = 1;
	mesh.zones[2].grantype = 4;
	mesh.zones[3].grantype = 2;

	int nZones = mesh.nZones;
	for (int i = 0; i < nZones; i++) {
		int tp = mesh.zones[i].grantype;
		mesh.zones[i].bnd = new Boundary[1];
		if (tp == 1) {
			mesh.zones[i].bnd[0].rules = new int[2];
			mesh.zones[i].bnd[0].vals = new double[1];
		}
		if (tp == 2) {
			mesh.zones[i].bnd[0].rules = new int[1];
			mesh.zones[i].bnd[0].vals = new double[5];
		}
	}

	mesh.zones[1].bnd[0].rules[0] = 0;
	mesh.zones[1].bnd[0].rules[1] = 2;
	mesh.zones[1].bnd[0].vals[0] = 0;

	mesh.zones[3].bnd[0].rules[0] = 1;
	mesh.zones[3].bnd[0].vals[0] = 0.08;
	mesh.zones[3].bnd[0].vals[1] = 867.9;
	mesh.zones[3].bnd[0].vals[2] = 75.1;
	mesh.zones[3].bnd[0].vals[3] = 1931.;
	mesh.zones[3].bnd[0].vals[4] = 0.;

}

void Gradients(Cell* cells, Grid mesh, Gradient* (&gr), parameters* p, int Nm)
{
	int nCells = mesh.nCells;
	// Значение вектора Vc в центре ячейки
	double* Vc = new double[Nm];
	// Значение вектора Vk-Vc у соседа k
	double* dVk = new double[Nm];
	// Временный вектор градиентов для всех компонентов вектора V
	Vector* gr1 = new Vector[Nm];

	for (int i = 0; i < Nm; i++) gr1[i].cx = new double[2];

	for (int i = 0; i < nCells; i++) {

		int NB = cells[i].get_nFaces(); //число окружающих граней

		for (int m = 0; m < Nm; m++) {
			gr1[m].cx[0] = 0;
			gr1[m].cx[1] = 0;
		}

		for (int m = 0; m < Nm; m++)
			Vc[m] = p[i].V[m]; // Промежуточный векор

		for (int k = 0; k < NB; k++) {
			int fType = cells[i].fType[k];

			int nb = cells[i].cells[k];
			int nf = cells[i].faces[k];

			if (fType == 0) {
				for (int m = 0; m < Nm; m++)
					dVk[m] = p[nb].V[m] - Vc[m];
			}
			else { // граница NS
				int z = mesh.faces[nf].zone;
				int btype = mesh.zones[z].grantype;

				if (btype == 1) {
					int vel = mesh.zones[z].bnd[0].rules[0]; // 0-проскальзывание 1-без проскальзывания
					int temp = mesh.zones[z].bnd[0].rules[1]; // 1-темпа на стенке 2-тепловой поток
					double value = mesh.zones[z].bnd[0].vals[0]; // температруа на стенке или тепловой поток
					if (vel == 0) {
						dVk[1] = 0.;
						dVk[2] = 0.;
					}
					if (vel == 1) {
						dVk[1] = 0. - Vc[1];
						dVk[2] = 0. - Vc[2];
					}
					if (temp == 1) {
						double Tw = value; //температура на стенке
						double hw = p[i].Cp * Tw;
						dVk[3] = hw - Vc[3];
						double ro = p[i].p * p[i].Gm / (8314.41 * Tw);
						dVk[0] = ro - Vc[0];
					}
					if (temp == 2) {
						dVk[3] = 0.;
						dVk[0] = 0.;
					}
				}
				if (btype == 2) {
					double u = mesh.zones[z].bnd[0].vals[1];
					double T = mesh.zones[z].bnd[0].vals[2];
					double P = mesh.zones[z].bnd[0].vals[3];

					double ro = P * p[i].Gm / (8314.41 * T);
					double h = p[i].Cp * T;

					dVk[0] = ro - Vc[0];
					dVk[1] = u - Vc[1];
					dVk[2] = 0. - Vc[2];
					dVk[3] = h - Vc[3];
				}
				if (btype == 3) {
					dVk[0] = 0.;
					dVk[1] = 0.;
					dVk[2] = 0. - Vc[2];
					dVk[3] = 0.;
				}
				if (btype == 4) {
					dVk[0] = 0.;
					dVk[1] = 0.;
					dVk[2] = 0.;
					dVk[3] = 0.;
				}
			}
			for (int m = 0; m < Nm; m++) {
				gr1[m].cx[0] += cells[i].wk[k] * dVk[m] * cells[i].ck[k].cx[0];
				gr1[m].cx[1] += cells[i].wk[k] * dVk[m] * cells[i].ck[k].cx[1];
			}
		}

		for (int m = 0; m < Nm; m++) {
			gr[i].g[m].cx[0] = gr1[m].cx[0];
			gr[i].g[m].cx[1] = gr1[m].cx[1];
		}
	}
}

void matrixDiag(double A[4][4], double L[4], double B[4][4])
{
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			B[i][j] = A[i][j] * L[j];
		}
	}
}

void matrixMatrix(double A[4][4], double B[4][4], double C[4][4])
{
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			C[i][j] = 0;
			for (int k = 0; k < 4; k++) {
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
}

void matrixVector(double A[4][4], double B[4], double C[4])
{
	for (int i = 0; i < 4; i++) {
		C[i] = 0;
		for (int j = 0; j < 4; j++) {
			C[i] += A[i][j] * B[j];
		}
	}
}

void printMatrix(double A[4][4])
{
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << setw(15) << A[i][j];
		}
		cout << endl;
	}
}

void sMatr(double A[4][4], double u, double v, double ro, double p, double h, double nx, double ny, int iMod)
{
	double gam = h / (h * p / ro);
	double a = sqrt(gam * p / ro);

	double beta = gam - 1.;
	double alfa = 0.5 * (u * u + v * v);

	double Ht = h + alfa;
	double Et = Ht - p / ro;

	double ly = nx;
	double lx = -ny;

	double U_ = u * nx + v * ny;
	double V_ = u * lx + v * ly;

	double S[4][4], S_[4][4];

	// Матрица S
	S[0][0] = a * a - alfa * beta;
	S[0][1] = beta * u;
	S[0][2] = beta * v;
	S[0][3] = -beta;

	S[1][0] = -V_;
	S[1][1] = lx;
	S[1][2] = ly;
	S[1][3] = 0.;

	S[2][0] = alfa * beta - U_ * a;
	S[2][1] = a * nx - beta * u;
	S[2][2] = a * ny - beta * v;
	S[2][3] = beta;

	S[3][0] = alfa * beta + U_ * a;
	S[3][1] = -a * nx - beta * u;
	S[3][2] = -a * ny - beta * v;
	S[3][3] = beta;

	// Матрица S^(-1)
	double a2 = a * a;
	S_[0][0] = 1. / a2;
	S_[0][1] = 0.;
	S_[0][2] = 1. / (2. * a2);
	S_[0][3] = 1. / (2. * a2);

	S_[1][0] = u / a2;
	S_[1][1] = lx;
	S_[1][2] = (u + a * nx) / (2. * a2);
	S_[1][3] = (u - a * nx) / (2. * a2);

	S_[2][0] = v / a2;
	S_[2][1] = ly;
	S_[2][2] = (v + a * ny) / (2. * a2);
	S_[2][3] = (v - a * ny) / (2. * a2);

	S_[3][0] = alfa / a2;
	S_[3][1] = V_;
	S_[3][2] = (Ht + a * U_) / (2. * a2);
	S_[3][3] = (Ht - a * U_) / (2. * a2);


	double L[4];
	L[0] = U_;
	L[1] = U_;
	L[2] = U_ + a;
	L[3] = U_ - a;
	// поправка
	double z_static = 0.5;
	double eps = z_static * (abs(U_) + a);

	if (iMod == 1) {
		// Матрица A
	}
	if (iMod == 2) {
		// Матрица A+
		for (int i = 0; i < 4; i++) {
			L[i] = 0.5 * (L[i] + sqrt(L[i] * L[i] + eps * eps));
		}
	}
	if (iMod == 3) {
		// Матрица A-
		for (int i = 0; i < 4; i++) {
			L[i] = 0.5 * (L[i] - sqrt(L[i] * L[i] + eps * eps));
		}
	}

	double Tmp[4][4];
	matrixDiag(S_, L, Tmp);
	matrixMatrix(Tmp, S, A);

}

void ConvectNS(parameters* p, changes* (&du), Grid mesh, Cell* cells, double dt, int Nm)
{
	int nFaces = mesh.nFaces;

	double* Fc = new double[Nm];

	double* UL = new double[Nm];
	double* UR = new double[Nm];
	double* ULR = new double[Nm];

	double A[4][4];
	double* AUl = new double[Nm];
	double* AUr = new double[Nm];


	for (int i = 0; i < nFaces; i++) {

		Face face = mesh.faces[i];
		// Номера ячеек справа и слева
		int cr = face.cr;
		int cl = face.cl;
		// Номера узлов
		int n1 = face.nodes[0];
		int n2 = face.nodes[1];
		// Длина грани
		double dl = face.length;
		// Точки узлов
		Point x1 = mesh.nodes[n1];
		Point x2 = mesh.nodes[n2];
		// Проекции нормали
		double nx = -(x2.y - x1.y) / dl;
		double ny = -(x2.x - x1.x) / dl;

		if (face.is_bound) {
			int c = max(cr, cl);
			double S = cells[c].get_S();

			int z = face.zone;
			int grantype = mesh.zones[z].grantype;

			if (grantype == 1 || grantype == 3) {
				Fc[0] = 0.;
				Fc[1] = p[c].p * nx;
				Fc[2] = p[c].p * ny;
				Fc[3] = 0.;
			}

			if (grantype == 2) {
				// скорость на входе
				double VV = mesh.zones[z].bnd[0].vals[1];
				// угол потока
				double angle = mesh.zones[z].bnd[0].vals[4];
				// компоненты скорости
				double u_ = VV * cos(angle);
				double v_ = VV * sin(angle);

				double T_ = mesh.zones[z].bnd[0].vals[2];
				double p_ = mesh.zones[z].bnd[0].vals[3];

				double Gm = 28.97;
				double Gam = 1.4;

				// расчет
				double R = 8314.41 / Gm;
				double ro_ = p_ / (R * T_);

				double h_ = p_ * Gam / ((Gam - 1.) * ro_);
				double H_ = h_ + 0.5 * (u_ * u_ + v_ * v_);
				double un = u_ * nx + v_ * ny;

				Fc[0] = ro_ * un;
				Fc[1] = ro_ * un * u_ + p_ * nx;
				Fc[2] = ro_ * un * v_ + p_ * ny;
				Fc[3] = ro_ * un * H_;
			}
			if (grantype == 4) {
				//потоки на выходной границе
				double un = p[c].u * nx + p[c].v * ny;
				double ro_ = p[c].ro;

				Fc[0] = ro_ * un;
				Fc[1] = ro_ * un * p[c].u + p[c].p * nx;
				Fc[2] = ro_ * un * p[c].v + p[c].p * ny;
				Fc[3] = ro_ * un * p[c].H;
			}
			if (cr >= 0) {
				// площадь правой ячейки
				double Sr = cells[cr].get_S();
				for (int m = 0; m < Nm; m++)
					du[cr].dU[m] += Fc[m] * dl / Sr * dt;
			}
			if (cl >= 0) {
				// площадь левой ячейки
				double Sl = cells[cr].get_S();
				for (int m = 0; m < Nm; m++)
					du[cl].dU[m] += -Fc[m] * dl / Sl * dt;
			}
		}
		else {

			for (int m = 0; m < Nm; m++) {
				UL[m] = p[cl].U1[m];
				UR[m] = p[cr].U1[m];
				ULR[m] = 0.5 * (UL[m] + UR[m]);
			}
			// Среднее значение
			double u_, v_, h_, E_, ro_, p_;

			ro_ = ULR[0];
			u_ = ULR[1] / ro_;
			v_ = ULR[2] / ro_;
			E_ = ULR[3] / ro_;

			double q = 0.5 * (u_ * u_ + v_ * v_);
			double e_ = E_ - q;

			double gam_ = 0.5 * (p[cr].Gam + p[cl].Gam);
			h_ = e_ * gam_;
			p_ = ro_ * (gam_ - 1) * e_;

			// Матрица A+
			int iMod = 2;
			sMatr(A, u_, v_, ro_, p_, h_, nx, ny, iMod);
			// Матрица A+ * UL
			matrixVector(A, UL, AUl);

			// Матрица A-
			iMod = 3;
			sMatr(A, u_, v_, ro_, p_, h_, nx, ny, iMod);
			// Матрица A- * UR
			matrixVector(A, UR, AUr);

			// Вектор невязких потоков 
			for (int m = 0; m < Nm; m++) {
				Fc[m] = AUl[m] + AUr[m];
			}

			double Sr = cells[cr].get_S(); // площадь правой ячейки
			double Sl = cells[cl].get_S(); // площадь левой ячейки

			// Приращение основного вектора за счет невязких потоков
			for (int m = 0; m < Nm; m++) {
				du[cr].dU[m] += Fc[m] * dl / Sr * dt;
				du[cl].dU[m] += -Fc[m] * dl / Sl * dt;
			}
		}
	}
}
