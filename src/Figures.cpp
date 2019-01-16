#include <fstream>

#include <iostream>

#include "Figures.h"
#include "Parameters.h"
#include "Probabilities.h"
#include "Peak.h"

void Precision()
{
	double L, E, a, exact, nu;

	nu = +1; // +1 for neutrinos, -1 for antineutrinos
	L = 1300.; // km, for DUNE

	std::ofstream data("data/Precision.txt");

	for (int i = 0; i < LAST_EXPRESSION; i++)
		data << Name((expression)i) << " ";
	data << std::endl;

	// Modify the energy range of the figure here
	for (E = nu * 0.3; fabs(E) <= 10; E *= 1.005)
	{
		// Take Ye=0.5 and rho=3 g/cc
		a = E * 0.5 * 3.0 * YerhoE2a;
		exact = Pmue(a, L, E, zs);

		data << fabs(E) << " ";

		for (int i = 0; i < LAST_EXPRESSION; i++)
			data << fabs(exact - Pmue(a, L, E, (expression)i)) / exact << " ";

		data << std::endl;
	} // E
	data.close();
}

void Probabilities()
{
	double a, L, E;

	L = 1300.; // km, for DUNE

	std::ofstream data("data/Probabilities.txt");

	for (int i = 0; i < LAST_EXPRESSION; i++)
		data << Name((expression)i) << " ";
	data << std::endl;

	// Modify the energy range of the figure here
	for (E = 0.3; E <= 10; E *= 1.001)
	{
		data << E << " ";
		for (int i = 1; i > -2; i -= 2)
		{
			E = i * fabs(E);
			// Take Ye=0.5 and rho=3 g/cc
			a = E * 0.5 * 3.0 * YerhoE2a;

			for (int j = 0; j < LAST_EXPRESSION; j++)
				data << Pmue(a, L, E, expression(j)) << " ";

		} // i, nu/nubar
		E = fabs(E);
		data << std::endl;
	} // E
	data.close();
}

void Peak_Precision_Write(expression m, std::ofstream &data)
{
	double E1_exact, P1_exact, E2_exact, P2_exact;
	double E1, P1, E2, P2;
	double dEE1, dPP1, dEE2, dPP2;

	First_Peak(zs, &E1_exact, &P1_exact);
	Second_Peak(zs, &E2_exact, &P2_exact);

	First_Peak(m, &E1, &P1);
	Second_Peak(m, &E2, &P2);

	dEE1 = fabs(E1 - E1_exact) / E1_exact;
	dPP1 = fabs(P1 - P1_exact) / P1_exact;
	dEE2 = fabs(E2 - E2_exact) / E2_exact;
	dPP2 = fabs(P2 - P2_exact) / P2_exact;

	data << Name(m) << " ";

	data << dEE1 << " ";
	data << dPP1 << " ";
	data << dEE2 << " ";
	data << dPP2 << " ";

	data << std::endl;
}

void Peak_Precision()
{
	std::ofstream data("data/Peak_Precision.txt");

	for (int i = 0; i < LAST_EXPRESSION; i++)
	{
		Peak_Precision_Write((expression)i, data);
	}
	data.close();
}

