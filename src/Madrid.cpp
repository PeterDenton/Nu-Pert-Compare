#include "Madrid.h"
#include "Parameters.h"

#define square(x) ((x)*(x))

namespace Madrid
{
double Pmue(double a, double L, double E)
{
	L *= -1; // maps e->mu on to mu->e
	double A;
	double D12, D13;
	double Btm;
	double Jt;
	double first, second, third;

	A = eVsqkm_to_GeV * a / (2 * E);

	D12 = eVsqkm_to_GeV * Dmsq21 / (2 * E);
	D13 = eVsqkm_to_GeV * Dmsq31 / (2 * E);

	Btm = fabs(A - D13);

	Jt = c13 * s212 * s223 * s213;

	first = s23sq * square(s213) * square(D13 / Btm) * square(sin(0.5 * Btm * L));
	second = c23sq * square(s212) * square(D12 / A) * square(sin(0.5 * A * L));
	third = Jt * (D12 / A) * (D13 / Btm) * sin(0.5 * A * L) * sin(0.5 * Btm * L) * cos(delta - 0.5 * D13 * L);

	return first + second + third;
}
}

