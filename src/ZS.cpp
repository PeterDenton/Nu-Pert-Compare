#include <cmath>
#include <complex>

#include "ZS.h"
#include "Parameters.h"

#define square(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))

namespace ZS {
double Pmue(double a, double L, double E)
{
	double A = Dmsq21 + Dmsq31 + a;
	double B = Dmsq31 * Dmsq21 + a * (Dmsq31 * c13sq + Dmsq21 * (c13sq * c12sq + s13sq));
	double C = a * Dmsq31 * Dmsq21 * c13sq * c12sq;
	// typo fixed:  the denominator is the 1.5 power not the 1/3 power
	double S = cos(acos((2 * cube(A) - 9 * A * B + 27 * C) / (2 * pow(square(A) - 3 * B, 1.5))) / 3.);

	double M1sq = A / 3. - sqrt(square(A) - 3 * B) * S / 3. - sqrt(square(A) - 3 * B) * sqrt(1 - square(S)) / sqrt(3.);
	double M2sq = A / 3. - sqrt(square(A) - 3 * B) * S / 3. + sqrt(square(A) - 3 * B) * sqrt(1 - square(S)) / sqrt(3.);
	double M3sq = A / 3. + 2 * sqrt(square(A) - 3 * B) * S / 3.;

	double DMsq21 = M2sq - M1sq;
	double DMsq31 = M3sq - M1sq;
	double DMsq32 = M3sq - M2sq;

	double alpha = Dmsq31 * c13sq + Dmsq21 * (c13sq * c12sq + s13sq);
	double beta = Dmsq31 * c13sq * Dmsq21 * c12sq;
	double E_ = (Dmsq31 * (M3sq - Dmsq21) - Dmsq21 * (M3sq - Dmsq31) * s12sq) * c13 * s13; // underscore added to differentiate the E function from the article from the energy
	double F = (M3sq - Dmsq31) * Dmsq21 * c12 * s12 * c13;

	double s12msq = -(square(M2sq) - alpha * M2sq + beta) * DMsq31 / (DMsq32 * (square(M1sq) - alpha * M1sq + beta) - DMsq31 * (square(M2sq) - alpha * M2sq + beta));
	double s13msq = (square(M3sq) - alpha * M3sq + beta) / (DMsq31 * DMsq32);
	double s23msq = (square(E_) * s23sq + square(F) * c23sq + 2 * E_ * F * c23 * s23 * cd) / (square(E_) + square(F));

	double c12m = sqrt(1 - s12msq);
	double s12m = sqrt(s12msq);
	double c13m = sqrt(1 - s13msq);
	double s13m = sqrt(s13msq);
	double c23m = sqrt(1 - s23msq);
	double s23m = sqrt(s23msq);

	std::complex<double> eidm, eidm_num, eidm_den1, eidm_den2;

	// typo fixed: it should be s23*c23 not s23*c23sq
	eidm_num = (square(E_) * conj(eid) - square(F) * eid) * s23 * c23 + E_ * F * (c23sq - s23sq);
	eidm_den1 = square(E_) * s23sq + square(F) * c23sq + 2 * E_ * F * c23 * s23 * cd;
	eidm_den2 = square(E_) * c23sq + square(F) * s23sq - 2 * E_ * F * c23 * s23 * cd;
	eidm = conj(eidm_num / sqrt(eidm_den1 * eidm_den2));


	double Delta21, Delta31, Delta32;
	std::complex<double> A31, A21, Amue;

	Delta21 = eVsqkm_to_GeV * DMsq21 * L / (4 * E);
	Delta31 = eVsqkm_to_GeV * DMsq31 * L / (4 * E);
	Delta32 = Delta31 - Delta21;

	A31 = 2 * s13m * c13m * s23m * sin(Delta31);
	A21 = 2 * s12m * c13m * (c12m * c23m * conj(eidm) - s12m * s13m * s23m) * sin(Delta21);

	Amue = A31 + exp(std::complex<double>(0, -Delta32)) * A21;

	return std::norm(Amue);
}
}

