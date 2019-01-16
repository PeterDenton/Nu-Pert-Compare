#include <cmath>

#include "AJLOS.h"
#include "Parameters.h"

#define square(x) ((x)*(x))

namespace AJLOS
{
// eq. 31
double Pmue31(double a, double L, double E)
{
	L *= -1; // maps e->mu on to mu->e
	double A, Delta, alpha;
	double P1, P2, P3;

	A = a / Dmsq31;
	Delta = eVsqkm_to_GeV * Dmsq31 * L / (4 * E);
	alpha = Dmsq21 / Dmsq31;

	double sinAD = sin(A * Delta);
	double sinA1D = sin((A - 1) * Delta);

	P1 = square(alpha * s212 * c23 * sinAD / A);
	P2 = square(2 * s13 * s23 * sinA1D / (A - 1));
	P3 = 2 * alpha * s13 * s212 * s223 * cos(Delta - delta) * sinAD * sinA1D / (A * (A - 1));

	return P1 + P2 + P3;
}
// eqs. 47-48
double Pmue48(double a, double L, double E)
{
	L *= -1; // maps e->mu on to mu->e
	double A, Delta, C13, alpha;
	double P0, P11, P12, P13, P14;

	A = a / Dmsq31;
	Delta = eVsqkm_to_GeV * Dmsq31 * L / (4 * E);
	C13 = sqrt(square(s213) + square(A - c213));
	alpha = Dmsq21 / Dmsq31;

	double sinCD = sin(C13 * Delta);
	double cosCD = cos(C13 * Delta);

	P0 = square(s23 * s213 * sinCD / C13);
	P11 = -2 * sinCD * square(s12 * s23 * s213 / C13);
	P12 = Delta * cosCD * (1 - A * c213) / C13 - A * sinCD * (c213 - A) / square(C13);
	P13 = s13 * s212 * s223 * sinCD / (A * square(C13));
	P14 = sd * (cosCD - cos((1 + A) * Delta)) * C13 + cd * (C13 * sin((1 + A) * Delta) - (1 - A * c213) * sinCD);

	return P0 + alpha * (P11 * P12 + P13 * P14);
}
}

