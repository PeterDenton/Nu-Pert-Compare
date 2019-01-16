#include "Peak.h"
#include "Probabilities.h"
#include "Minimization.h"
#include "Parameters.h"

// to find the maximum with a minimizer we need -f(x)
double Pmue_helper(double E, void *params)
{
	expression m = *(expression*)params;
	return -Pmue(1.5 * E * YerhoE2a, 1300., E, m);
}

// calls the minimizer to find the first peak, parameters are chosen to ensure that we always fall into the minimum
void First_Peak(expression m, double *E, double *P)
{
	if (m == zeroth) return;

	*E = gss_min(Pmue_helper, 1.2, 10., 1e-7, &m);
	*P = -Pmue_helper(*E, &m);
}

// at the second peak it is harder to enusre that you know you are in the minimum
// first, find the first peak
// then (coarsely) walk down in energy from there until you have passed the oscillation minimum, then check for the next peak
// then find the peak within the range known to be the second peak
void Second_Peak(expression m, double *E, double *P)
{
	if (m == zeroth) return;

	double E0, E1, P0, P1, step;
	bool min_passed = false;

	step = 0.05;

	First_Peak(m, &E1, &P1);

	E0 = E1 - step;
	P0 = -Pmue_helper(E0, &m);

	while (not min_passed or P0 > P1)
	{
		E0 -= step;
		E1 -= step;
		P1 = P0;
		P0 = -Pmue_helper(E0, &m);

		if (P0 > P1) min_passed = true;
	}

	*E = gss_min(Pmue_helper, E0, E1 + step, 1e-7, &m);
	*P = -Pmue_helper(*E, &m);
}

