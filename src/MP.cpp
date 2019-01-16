#include <cmath>

#include "MP.h"
#include "Parameters.h"

#define square(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))

namespace MP
{
double Pmue(double a, double L, double E)
{
	L *= -1; // maps e->mu on to mu->e
	double Jr, lp, lm, l0, Dlpm, Dlp0, Dlm0;

	Jr = c12 * s12 * c23 * s23 * c13sq * s13;

	lp = ((Dmsqee + a) + mo_sign * sqrt(square(Dmsqee - a) + 4 * s13sq * a * Dmsqee)) / 2. + s12sq * Dmsq21;
	lm = ((Dmsqee + a) - mo_sign * sqrt(square(Dmsqee - a) + 4 * s13sq * a * Dmsqee)) / 2. + s12sq * Dmsq21;
	l0 = c12sq * Dmsq21;

	Dlpm = lp - lm;
	Dlp0 = lp - l0;
	Dlm0 = lm - l0;

	return (s23sq * square(s213) + 4 * eps * Jr * cd * ((Dlpm - (Dmsqee - a)) / Dlp0)) * square(Dmsqee / Dlpm) * square(sin(eVsqkm_to_GeV * Dlpm * L / (4 * E))) + 8 * eps * Jr * cube(Dmsqee) / (Dlpm * Dlp0 * Dlm0) * sin(eVsqkm_to_GeV * Dlpm * L / (4 * E)) * sin(eVsqkm_to_GeV * Dlm0 * L / (4 * E)) * cos(delta - eVsqkm_to_GeV * Dlp0 * L / (4 * E));
}
}

