#ifndef Probabilities_H
#define Probabilities_H

#include <string>

enum expression : int
{
	zs,
	madrid,
	ajlos31,
	fl,
	akt,
	mp,
	dmp0,
	dmp1,
	aks,
	mf,
	ajlos48,
	am2,
	am52,
	diag,
	vacuum,
	zeroth,
	dpz0,
	dpz2,
	LAST_EXPRESSION
};

double Pmue(double a, double L, double E, expression m);
std::string Name(expression m);

#endif
