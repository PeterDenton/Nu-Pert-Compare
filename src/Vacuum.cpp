#include "Vacuum.h"
#include "Parameters.h"

#define square(x) ((x)*(x))

namespace Vacuum {
double Pmue(double a, double L, double E)
{
	double L4E, Delta21, Delta31, Delta32, sDelta21, sDelta31, sDelta32, Jr, C31, C32, C21, D;

	L4E = eVsqkm_to_GeV * L / (4 * E);
	Delta21 = Dmsq21 * L4E;
	Delta31 = Dmsq31 * L4E;
	Delta32 = Delta31 - Delta21;

	sDelta21 = sin(Delta21);
	sDelta31 = sin(Delta31);
	sDelta32 = sin(Delta32);

	Jr = c13sq * s13 * c23 * s23 * c12 * s12;

	C31 = 4 * square(sDelta31) * (s23sq * s13sq * c13sq * c12sq + Jr * cd);
	C32 = 4 * square(sDelta32) * (s23sq * s13sq * c13sq * s12sq - Jr * cd);
	C21 = 4 * square(sDelta21) * (c13sq * s12sq * c12sq * (c23sq - s13sq * s23sq) + Jr * cd * c212);
	D = -8 * sDelta21 * sDelta31 * sDelta32 * Jr * sd;

	return C31 + C32 + C21 + D;
}
}
