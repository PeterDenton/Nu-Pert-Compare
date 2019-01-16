#include "MF.h"
#include "Parameters.h"

#define square(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))
#define quad(x) ((x)*(x)*(x)*(x))

namespace MF
{
double Pmue(double a, double L, double E)
{
	L *= -1; // maps e->mu on to mu->e
	double alpha, Deltah, Ch, sign, Ah;
	double P0, Psd, Pcd, P1, P2, P3;

	alpha = Dmsq21 / Dmsq31;
	Deltah = eVsqkm_to_GeV * Dmsq31 * L / (4. * E);
	Ah = a / Dmsq31;
	Ch = sqrt(square(Ah - c213) + square(s213));
	sign = (Ah < c213) ? +1 : -1;

	double sin_Ch_Deltah = sin(Ch * Deltah);

	// eqs. 36
	P0 = s23sq * square(s213 * sin_Ch_Deltah / Ch);
	Psd = alpha * sd * s212 * s13 * s223 / (Ah * Ch) * sin_Ch_Deltah * (cos(Ch * Deltah) - cos((1 + Ah) * Deltah));
	Pcd = alpha * cd * s212 * s13 * s223 / (Ah * Ch) * sin_Ch_Deltah * (sin((1 + Ah) * Deltah) - sign * sin_Ch_Deltah);
	P1 = -alpha * (1 - Ah * c213) / cube(Ch) * s12sq * square(s213) * s23sq * Deltah * sin(2 * Deltah * Ch)
		+ alpha * 2 * Ah * (-Ah + c213) / quad(Ch) * s12sq * square(s213) * s23sq * square(sin_Ch_Deltah);
	P2 = alpha * (-sign + Ch + sign * Ah * c213) / (2 * square(Ch) * Ah * c13sq) * c13 * s212 * s213 * s223 * square(sin_Ch_Deltah);
	P3 = square(alpha) * 2 * Ch * c23sq * square(s212) / (square(Ah) * c13sq * (-sign * Ah + Ch + sign * c213)) * square(sin((1 + Ah - sign * Ch) * Deltah * 0.5));

	return P0 + Psd + Pcd + P1 + P2 + P3;
}
}
