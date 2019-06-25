#include <cmath>
#include <complex>

#include "DPZ.h"
#include "Parameters.h"

#define sq(x) ((x)*(x))

namespace DPZ {
double Pmue(double a, double L, double E, int order)
{
	double la, lb, lc, lp, lm, l0, l1, l2, l3, tmp, phi, cphim13;

	// calculate eigenvalues
	la = a + s13sq * Dmsqee + s12sq * Dmsq21;
	lb = c12sq * Dmsq21;
	lc = c13sq * Dmsqee + s12sq * Dmsq21;

	tmp = sqrt(sq(la - lc) + sq(2 * s13 * c13 * Dmsqee));
	lp = 0.5 * (la + lc + tmp);
	lm = 0.5 * (la + lc - tmp);
	l0 = lb;

	phi = asin(sqrt((lp - lc) / (lp - lm)));
	cphim13 = cos(phi - t13);

	tmp = sqrt(sq(l0 - lm) + sq(2 * cphim13 * s12 * c12 * Dmsq21));
	l1 = 0.5 * (l0 + lm - tmp);
	l2 = 0.5 * (l0 + lm + tmp);
	l3 = lp;

	if (order == 2)
	{
		double spsisq, cpsisq, sphim13sq, epsp_Dmsqee_sq, l12, l22, l32;

		spsisq = (l2 - l0) / (l2 - l1);
		cpsisq = 1 - spsisq;

		sphim13sq = 1 - sq(cphim13);
		epsp_Dmsqee_sq = sphim13sq * s12sq * c12sq * sq(Dmsq21);

		l12 = -epsp_Dmsqee_sq * spsisq / (l3 - l1);
		l22 = -epsp_Dmsqee_sq * cpsisq / (l3 - l2);
		l32 = -l12 - l22;

		l1 += l12;
		l2 += l22;
		l3 += l32;
	} // second order eigenvalues

	// calculate submatrix eigenvalues
	double hee, hem, het, hmm, hmt, htt, sum, prod, xie, chie, xim, chim;
	hee = a + Dmsqee * s13sq + Dmsq21 * s12sq;
	hem = c13 * s12 * c12 * Dmsq21;
	het = s13 * c13 * Dmsqee;
	hmm = Dmsq21 * c12sq;
	hmt = -s13 * s12 * c12 * Dmsq21;
	htt = Dmsqee * c13sq + Dmsq21 * s12sq;

	// e
	sum = hmm + htt;
	prod = hmm * htt - sq(hmt);
	tmp = sqrt(sq(sum) - 4 * prod);
	xie = 0.5 * (sum + tmp);
	chie = 0.5 * (sum - tmp);

	// mu
	sum = hee + c23sq * htt + s23sq * hmm - 2 * s23 * c23 * cd * hmt;
	prod = hee * (c23sq * htt + s23sq * hmm - 2 * s23 * c23 * cd * hmt) - sq(fabs(c23 * het - s23 * conj(eid) * hem));
	tmp = sqrt(sq(sum) - 4 * prod);
	xim = 0.5 * (sum + tmp);
	chim = 0.5 * (sum - tmp);

	// calculate relevant Uaisq
	double Ue2sq, Ue3sq, Um3sq;

	Ue2sq = (l2 - xie) * (l2 - chie) / ((l2 - l1) * (l2 - l3));
	Ue3sq = (l3 - xie) * (l3 - chie) / ((l3 - l2) * (l3 - l1));
	Um3sq = (l3 - xim) * (l3 - chim) / ((l3 - l2) * (l3 - l1));

	// calculate mixing angles in matter
	double s13msq, c13msq, s23msq, c23msq, s12msq, c12msq;
	s13msq = Ue3sq;
	c13msq = 1 - s13msq;
	s12msq = Ue2sq / c13msq;
	s23msq = Um3sq / c13msq;

	c12msq = 1 - s12msq;
	c23msq = 1 - s23msq;

	// get CPV phase in matter
	double sdm, cdm;
	std::complex<double> eidm, I;

	sdm = s223 * sd / (2 * sqrt(s23msq * c23msq));
	cdm = sqrt(1 - sq(sdm)); // assume this gives the right quadrant
	I = std::complex<double>(0, 1);
	eidm = cdm + I * sdm;

	// kinematic terms
	double L4E, Delta21, Delta31, Delta32, sDelta21, sDelta31;
	L4E = eVsqkm_to_GeV * L / (4 * E);
	Delta21 = (l2 - l1) * L4E;
	Delta31 = (l3 - l1) * L4E;
	Delta32 = Delta31 - Delta21;

	sDelta21 = sin(Delta21);
	sDelta31 = sin(Delta31);
/*
	// amplitudes
	std::complex<double> A21, A31, Amue;

	A21 = 2 * sqrt(s12msq * c13msq) * (sqrt(c12msq * c23msq) * eidm - sqrt(s12msq * s13msq * s23msq)) * sDelta21;
	A31 = 2 * sqrt(s13msq * c13msq * s23msq) * sDelta31;
	Amue = A31 + exp(I * Delta32) * A21;

	return std::norm(Amue);
*/

double sDelta32 = sin(Delta32);
double Jrm, C31, C32, C21, D;
	Jrm = c13msq * sqrt(s23msq * c23msq * s13msq * c12msq * s12msq);
	C31 = s23msq * s13msq * c13msq * c12msq + Jrm * cdm;
	C32 = s23msq * s13msq * c13msq * s12msq - Jrm * cdm;
	C21 = c13msq * s12msq * c12msq * (c23msq - s13msq * s23msq) + Jrm * cdm * (c12msq - s12msq);
	D = -Jrm * sdm;
	return 4 * (C31 * sq(sDelta31) + C32 * sq(sDelta32) + C21 * sq(sDelta21)) + 8 * D * sDelta21 * sDelta31 * sDelta32;
}
}

