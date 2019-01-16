#ifndef Parameters_H
#define Parameters_H

#include <complex>

// oscillation parameters
extern double t12, t13, t23, Dmsq21, Dmsq31, delta;
extern double Dmsq32, Dmsqee, eps;
extern double c12, s12, c12sq, s12sq, s212, c212;
extern double c13, s13, c13sq, s13sq, s213, c213;
extern double c23, s23, c23sq, s23sq, s223, c223;
extern double cd, sd;
extern double mo_sign;
extern std::complex<double> eid;

void Recalc_Parameters();

// unit conversions
extern double eVsqkm_to_GeV;
extern double YerhoE2a;

#endif
