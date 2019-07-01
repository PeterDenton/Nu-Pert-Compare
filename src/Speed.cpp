#include <cmath>
#include <chrono>
#include <fstream>

#include "Speed.h"
#include "Probabilities.h"
#include "Parameters.h"
#include "Peak.h"

double Speed(expression m, int n)
{
	double L_min, L_scale, L_max, E_min, E_scale, E_max;
	double a;

	// initialize the grid
	double *Ls = new double[n];
	double *Es = new double[n];

	L_min = 1e1;
	L_max = 1e4;

	E_min = 1e-1;
	E_max = 1e2;

	L_scale = pow(L_max / L_min, 1. / n);
	E_scale = pow(E_max / E_min, 1. / n);

	// fill in arrays outside of timing
	for (int i = 0; i < n; i++)
	{
		Ls[i] = L_min * pow(L_scale, i);
		Es[i] = E_min * pow(E_scale, i);
	} // i, n

	// begin timing loop
	std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			a = 1.5 * Es[j] * YerhoE2a;
			Pmue(a, Ls[i], Es[j], m);
		} // j, n, E
	} // i, n, L
	std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();

	// clean up
	delete[] Ls;
	delete[] Es;

	return std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count() / pow(n, 2);
}

void Write(expression m, int n, double P_exact, std::ofstream &data)
{
	double E, P;
	First_Peak(m, &E, &P);

	data << Name(m) << " ";
	data << Speed(m, n) << " ";
	data << fabs(P - P_exact) / P_exact << " ";
	data << std::endl;
}

void Speed()
{
	double P_exact, E_exact;
	First_Peak(zs, &E_exact, &P_exact);

	std::ofstream data("data/Speed.txt");

	Write(zs, 5e3, P_exact, data);
	Write(diag, 1e4, P_exact, data);
	Write(dmp0, 5e3, P_exact, data);
	Write(dmp1, 5e3, P_exact, data);
	Write(mp, 1e4, P_exact, data);
	Write(am2, 1e4, P_exact, data);
	Write(am52, 1e4, P_exact, data);
	Write(madrid, 5e4, P_exact, data);
	Write(fl, 1e4, P_exact, data);
	Write(ajlos31, 5e4, P_exact, data);
	Write(ajlos48, 1e4, P_exact, data);
	Write(mf, 1e4, P_exact, data);
	Write(aks, 1e4, P_exact, data);
	Write(zeroth, 2e4, P_exact, data);
	Write(vacuum, 5e4, P_exact, data);
	Write(akt, 5e3, P_exact, data);
	Write(dpz0, 5e3, P_exact, data);
	Write(dpz2, 5e3, P_exact, data);

	data.close();
}

