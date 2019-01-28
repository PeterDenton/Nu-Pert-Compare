#include <cmath>

#include "AKT.h"
#include "Parameters.h"

#define square(x) ((x)*(x))

namespace AKT
{
double Pmue(double a, double L, double E)
{
	double lp_sqrt, lpp_sqrt;
	double lpp, lpm, lppp, lppm;
	double s12psq, c12psq, s13psq, c13psq, c212p;
	double l1, l2, l3;
	double L4E, Delta21, Delta31, Delta32, sDelta21, sDelta31, sDelta32;
	double Jrm, C31, C32, C21, D;

	// tan(2theta)=num/den
	// the sqrt in the eigenvalues in sqrt(num^2+den^2)
	double den12, num13, den13;
	den12 = Dmsq21 * c212 - c13sq * a;
	num13 = Dmsqee * s213;
	den13 = Dmsqee * c213 - a;

	lp_sqrt = sqrt(square(Dmsq21 - a * c13sq) + 4 * a * c13sq * s12sq * Dmsq21);
	s12psq = 0.5 * (1 - den12 / lp_sqrt);

	lpp = 0.5 * (Dmsq21 + a * c13sq + lp_sqrt);
	lpm = lpp - lp_sqrt;

	lpp_sqrt = sqrt(square(lpp - (Dmsq31 + a * s13sq)) + 4 * square(a) * s12psq * c13sq * s13sq);
	s13psq = 0.5 * (1 - den13 / sqrt(square(num13) + square(den13)));

	lppp = 0.5 * (lpp + Dmsq31 + a * s13sq + lpp_sqrt);
	lppm = lppp - lpp_sqrt;

	l1 = lpm;
	l2 = lppm;
	l3 = lppp;

	c12psq = 1 - s12psq;
	c13psq = 1 - s13psq;

	c212p = c12psq - s12psq;

	Jrm = s23 * c23 * c13psq * sqrt(s13psq * c12psq * s12psq);
	C31 = s23sq * s13psq * c13psq * c12psq + Jrm * cd;
	C32 = s23sq * s13psq * c13psq * s12psq - Jrm * cd;
	C21 = c13psq * s12psq * c12psq * (c23sq - s13psq * s23sq) + Jrm * cd * c212p;
	D = -Jrm * sd;

	L4E = eVsqkm_to_GeV * L / (4 * E);
	Delta21 = (l2 - l1) * L4E;
	Delta31 = (l3 - l1) * L4E;
	Delta32 = Delta31 - Delta21;

	sDelta21 = sin(Delta21);
	sDelta31 = sin(Delta31);
	sDelta32 = sin(Delta32);

	return 4 * (C31 * square(sDelta31) + C32 * square(sDelta32) + C21 * square(sDelta21)) + 8 * D * sDelta21 * sDelta31 * sDelta32;
}
}
