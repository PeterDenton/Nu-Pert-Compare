#include <complex>

#include "FL.h"
#include "Parameters.h"

namespace FL
{
double Pmue(double a, double L, double E)
{
	std::complex<double> G1, G2, A;
	double Delta, Delta_sol, Delta_1, Delta_2;

	Delta = Dmsq32 / (4 * E);
	Delta_sol = Dmsq21 / (4 * E);

	Delta_1 = 2 * Delta - a / (2 * E);
	Delta_2 = -a / (2 * E);

	G1 = Delta * s213 * eid;
	G2 = -Delta_sol * s212;
	A = G1 * s23 * (exp(std::complex<double>(0, eVsqkm_to_GeV * Delta_1 * L)) - 1.) / Delta_1 - G2 * c23 * (exp(std::complex<double>(0, eVsqkm_to_GeV * Delta_2 * L)) - 1.) / Delta_2;

	return std::norm(A);
}
}

