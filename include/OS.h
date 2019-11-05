#ifndef OS_H
#define OS_H

// From https://arxiv.org/abs/hep-ph/9910546

namespace OS
{
double Pmue(double a, double L, double E); // Cayley-Hamilton
double Pmue2(double a, double L, double E); // exponentiation
double lambda(double a, double L, double E, int i);
}

#endif
