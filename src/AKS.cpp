#include <cmath>

#include "AKS.h"
#include "Parameters.h"

#define square(x) ((x)*(x))

namespace AKS
{
double Pmue(double a, double L, double E)
{
	double P1, P2, P3;
	double Delta21, Delta31;

	Delta21 = eVsqkm_to_GeV * Dmsq21 * L / (4 * E);
	Delta31 = eVsqkm_to_GeV * Dmsq31 * L / (4 * E);

	P1 = 4 * square(sin(Delta31)) * c13sq * s13sq * s23sq * (1 + 2 * (1 - 2 * s13sq) * a / Dmsq31);
	P2 = 4 * Delta31 * sin(2 * Delta31) * c13sq * s13 * s23 * (-a / Dmsq31 * s13 * s23 * (1 - 2 * s13sq) + (Dmsq21 / Dmsq31) * s12 * (-s13 * s23 * s12 + cd * c23 * c12));
	P3 = 8 * Delta21 * square(sin(Delta31)) * sd * c13sq * s13 * c23 * s23 * c12 * s12;

	return P1 + P2 - P3;
}
}

