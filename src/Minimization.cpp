#include <cmath>

#include "Minimization.h"
#include "Warning.h"

// from https://en.wikipedia.org/wiki/Golden_section_search
// a<b
// p is the parameters to pass into f
// tol is abs
template <class T>
T gss_min(T (*f)(T, void*), T a_initial, T b_initial, T tol, void *p)
{
	T gr = (sqrt(5.) - 1) / 2;
	T a = a_initial;
	T b = b_initial;
	T c = b - gr * (b - a);
	T d = a + gr * (b - a);
	while (std::abs(c - d) > tol) // abs
	{
		if (f(c, p) < f(d, p))
		{
			b = d;
			d = c;
			c = b - gr * (b - a);
		}
		else
		{
			a = c;
			c = d;
			d = a + gr * (b - a);
		}
	}
	T ret = (b + a) / 2.;

	// Warns user in case final minimum is close to an edge
	if (std::abs(ret - a_initial) < 2 * tol or std::abs(b_initial - ret) < 2 * tol)
		EdgeCaseWarning W(a_initial, b_initial, ret);
	return ret;
}

template double gss_min<double>(double (*)(double, void*), double, double, double, void*);
template long double gss_min<long double>(long double (*)(long double, void*), long double, long double, long double, void*);

double first_min(double (*f)(double x, void *p), double x0, double x_step, void *p)
{
	double x, f0, f1;
	x = x0;
	f0 = f(x, p);
	x += x_step;
	f1 = f(x, p);
	if (f1 >= f0)
		EdgeCaseWarning W(x0, x, (x + x0) / 2.);
	while (f1 < f0)
	{
		f0 = f1;
		x += x_step;
		f1 = f(x, p);
	}
	return x - x_step / 2.;
}

double minimizer(double (*f)(double x, void *p), double x0, double x_step, double tol, void *p)
{
	x0 = first_min(f, x0, x_step, p);

	x_step = fabs(x_step);
	return gss_min(f, x0 - x_step, x0 + x_step, tol, p);
}

