#ifndef AM_H
#define AM_H

// From http://arxiv.org/abs/1103.4387
// Order should be one of 0, 1, 1.5, 2, or 2.5
// Note that order=0 => P=0.

namespace AM
{
double Pmue(double a, double L, double E, double order);
}

#endif
