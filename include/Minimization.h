#ifndef Minimization_H
#define Minimization_H

template <class T>
T gss_min(T (*f)(T, void*), T a, T b, T tol, void *p);

double first_min(double (*f)(double x, void *p), double x0, double x_step, void *p);
double minimizer(double (*f)(double x, void *p), double x0, double x_step, double tol, void *p); // combination of both of the above

#endif
