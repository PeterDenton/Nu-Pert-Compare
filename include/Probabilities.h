#ifndef Probabilities_H
#define Probabilities_H

#include <string>

enum expression : int
{
	zs,
	madrid,
	fl,
	ajlos31,
	mp,
	dmp0,
	dmp1,
	akt,
	am2,
	am52,
	mf,
	aks,
	ajlos48,
	diag,
	vacuum,
	zeroth,
	LAST_EXPRESSION
};

double Pmue(double a, double L, double E, expression m);
std::string Name(expression m);

#endif
