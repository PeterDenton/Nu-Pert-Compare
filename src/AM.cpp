#include <cmath>

#include "AM.h"
#include "Parameters.h"

#define square(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))

namespace AM
{
double Pmue(double a, double L, double E, double order)
{
	L *= -1; // maps e->mu on to mu->e
	double P0, P1, P32, P2, P52;
	double rA, rDelta, Delta, DeltaL, Jrr;

	rA = a / Dmsq31;
	rDelta = Dmsq21 / Dmsq31;
	Delta = Dmsq31 / (2 * E);
	DeltaL = eVsqkm_to_GeV * Delta * L;
	Jrr = c12 * s12 * c23 * s23 * s13;

	P0 = 0;
	P1 = 0;
	P32 = 0;
	P2 = 0;
	P52 = 0;

	double sin1rADL05 = sin((1 - rA) * DeltaL * 0.5);
	double sinrADL05 = sin(rA * DeltaL * 0.5);
	double cosdDL05 = cos(delta - DeltaL * 0.5);

	if (order >= 1.)
		P1 = 4 * s23sq * s13sq * square(sin1rADL05 / (1 - rA));
	if (order >= 1.5)
		P32 = 8 * Jrr * rDelta / (rA * (1 - rA)) * cosdDL05 * sinrADL05 * sin1rADL05;
	if (order >= 2.)
		P2 = 4 * c23sq * c12sq * s12sq * square(rDelta / rA) * square(sinrADL05)
			- 4 * s23sq * (square(s13sq * (1 + rA) / square(1 - rA)) - 2 * s12sq * s13sq * rDelta * rA / cube(1 - rA)) * square(sin1rADL05)
			+ 2 * s23sq * (2 * square(s13sq) * rA / cube(1 - rA) - s12sq * s13sq * rDelta / square(1 - rA)) * DeltaL * sin((1 - rA) * DeltaL);
	if (order >= 2.5)
		P52 = 8 * Jrr * s13sq * rDelta * rA / cube(1 - rA) * cd * square(sin1rADL05)
			+ 8 * Jrr * rDelta / (rA * (1 - rA)) * (-2 * s13sq * rA / square(1 - rA) + (c12sq - s12sq) * rDelta / rA + s12sq * rDelta * rA / (1 - rA))
			* cosdDL05 * sinrADL05 * sin1rADL05
			+ 8 * Jrr * s13sq * rDelta / square(1 - rA) * DeltaL * cosdDL05 * sinrADL05 * cos((1 - rA) * DeltaL * 0.5)
			- 4 * Jrr * s12sq * square(rDelta) / (rA * (1 - rA)) * DeltaL * cos(delta - rA * DeltaL * 0.5) * sinrADL05
			- 4 * Jrr * c12sq * square(rDelta) / (rA * (1 - rA)) * DeltaL * cos(delta - (1 + rA) * DeltaL * 0.5) * sin1rADL05
			- 4 * Jrr * rDelta / (rA * (1 - rA)) * (s13sq * rA / (1 - rA) - s12sq * rDelta) * DeltaL * cos(delta - (1 - rA) * DeltaL * 0.5) * sin1rADL05;

	return P0 + P1 + P32 + P2 + P52;
}
}
