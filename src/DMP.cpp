#include <cmath>

#include "DMP.h"
#include "Parameters.h"

#define square(x) ((x)*(x))

namespace DMP {
double Pmue(double a, double L, double E, int order)
{
	double Dmsqeea, c2phi, a12, cphisq, sphisq, s2phi, cphi13sq, Dl21, c2psi, cpsisq, spsisq, Dl31, L4E, Delta21, Delta31, Delta32, sDelta21, sDelta31, sDelta32, Jrm, C31, C32, C21, D;

	Dmsqeea = Dmsqee * sqrt(square(c213 - a / Dmsqee) + square(s213));
	c2phi = (Dmsqee * c213 - a) / Dmsqeea;
	a12 = 0.5 * (a + Dmsqee - Dmsqeea);

	cphisq = 0.5 * (1 + c2phi);
	sphisq = 0.5 * (1 - c2phi);
	s2phi = s213 * Dmsqee / Dmsqeea;
	cphi13sq = cphisq * c13sq + sphisq * s13sq + s2phi * c13 * s13;

	Dl21 = Dmsq21 * sqrt(square(c212 - a12 / Dmsq21) + cphi13sq * square(s212));
	c2psi = (Dmsq21 * c212 - a12) / Dl21;
	cpsisq = 0.5 * (1 + c2psi);
	spsisq = 0.5 * (1 - c2psi);

	Dl31 = Dmsq31 + 0.25 * a + 0.5 * (Dl21 - Dmsq21) + 0.75 * (Dmsqeea - Dmsqee);

	L4E = eVsqkm_to_GeV * L / (4 * E);
	Delta21 = Dl21 * L4E;
	Delta31 = Dl31 * L4E;
	Delta32 = Delta31 - Delta21;

	sDelta21 = sin(Delta21);
	sDelta31 = sin(Delta31);
	sDelta32 = sin(Delta32);

	Jrm = s23 * c23 * cphisq * sqrt(sphisq * cpsisq * spsisq);
	C31 = s23sq * sphisq * cphisq * cpsisq + Jrm * cd;
	C32 = s23sq * sphisq * cphisq * spsisq - Jrm * cd;
	C21 = cphisq * spsisq * cpsisq * (c23sq - sphisq * s23sq) + Jrm * cd * c2psi;
	D = -Jrm * sd;

	if (order == 1)
	{
		double cphi, sphi, s2psi, epsp, sphi13;
		double F1, F2, G1, G2, K1, K2;
		double Dl32;

		cphi = sqrt(cphisq);
		sphi = sqrt(sphisq);
		s2psi = 2 * sqrt(cpsisq * spsisq);
		sphi13 = sqrt(1 - cphi13sq);
		if (cphi > c13) sphi13 *= -1;
		epsp = eps * sphi13 * s12 * c12;

		F1 = spsisq * (0.25 * s2phi * s2psi * (c23sq + c2phi * s23sq) - s23 * c23 * cphi * (sphisq * spsisq + c2phi * cpsisq) * cd);
		F2 = cpsisq * (-0.25 * s2phi * s2psi * (c23sq + c2phi * s23sq) - s23 * c23 * cphi * (sphisq * cpsisq + c2phi * spsisq) * cd);

		G1 = -s2phi * (0.5 * s23sq * c2phi * s2psi - s23 * c23 * sphi * spsisq * cd);
		G2 = -s2phi * (-0.5 * s23sq * c2phi * s2psi - s23 * c23 * sphi * cpsisq * cd);

		K1 = -s23 * c23 * cphi * spsisq * (cphisq * cpsisq - sphisq) * sd;
		K2 = -s23 * c23 * cphi * cpsisq * (cphisq * spsisq - sphisq) * sd;

		Dl32 = Dl31 - Dl21;

		C21 += epsp * Dmsqee * (F1 / Dl31 + F2 / Dl32);
		C31 += epsp * Dmsqee * ((F1 + G1) / Dl31 - F2 / Dl32);
		C32 += epsp * Dmsqee * (-F1 / Dl31 + (F2 + G2) / Dl32);
		D += epsp * Dmsqee * (K1 / Dl31 - K2 / Dl32) * sd;
	}

	return 4 * (C31 * square(sDelta31) + C32 * square(sDelta32) + C21 * square(sDelta21)) + 8 * D * sDelta21 * sDelta31 * sDelta32;
}
}

