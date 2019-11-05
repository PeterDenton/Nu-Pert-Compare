#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/unsupported/Eigen/MatrixFunctions>

#include <cmath>

#include "OS.h"
#include "Parameters.h"

using Eigen::Matrix3d;
using Eigen::Matrix3cd;

#define sq(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))

namespace OS
{
double Pmue(double a, double L, double E)
{
	// convert the 2E factors
	// all of these are now in eV
	double A, E21, E31, E32;
	A = a / (2 * E * 1e9);
	E21 = Dmsq21 / (2 * E * 1e9);
	E31 = Dmsq31 / (2 * E * 1e9);
	E32 = Dmsq32 / (2 * E * 1e9);

	// real pmns matrix
	Matrix3d U;
	U(0, 0) = c12 * c13;
	U(0, 1) = s12 * c13;
	U(0, 2) = s13;
	U(1, 0) = -s12 * c23 - c12 * s23 * s13;
	U(1, 1) = c12 * c23 - s12 * s23 * s13;
	U(1, 2) = s23 * c13;
	U(2, 0) = s12 * s23 - c12 * c23 * s13;
	U(2, 1) = -c12 * s23 - s12 * c23 * s13;
	U(2, 2) = c23 * c13;

	// create the T matrix
	// in eV
	Matrix3d T;
	T(0, 0) = A * sq(U(0, 0)) + (-A - E21 - E31) / 3.;
	T(0, 1) = A * U(0, 0) * U(0, 1);
	T(0, 2) = A * U(0, 0) * U(0, 2);
	T(1, 0) = T(0, 1);
	T(1, 1) = A * sq(U(0, 1)) + (-A + E21 - E32) / 3.;
	T(1, 2) = A * U(0, 1) * U(0, 2);
	T(2, 0) = T(0, 2);
	T(2, 1) = T(1, 2);
	T(2, 2) = A * sq(U(0, 2)) + (-A + E31 + E32) / 3.;

	// rotate to Ttilde, Ttildesq
	Matrix3d Ttilde = U * T * U.adjoint(); // adjoint = dagger, Uinv = Udag
	Matrix3d Ttildesq = Ttilde * Ttilde;

	// coefficients of T matrix
	// eV^3 and eV^2 respectively
	double c0, c1;
	c0 = -T.determinant();
	c1 = T(0, 0) * T(1, 1) - sq(T(0, 1)) + T(0, 0) * T(2, 2) - sq(T(0, 2)) + T(1, 1) * T(2, 2) - sq(T(1, 2));

	// eigenvalues
	// in eV
	double tmp0, tmp1, tmp2, lambdas[3];
	tmp0 = atan2(-sqrt(-sq(c0) - 4 * cube(c1) / 27.), -c0) / 3.;
	tmp1 = sqrt(-c1 / 3.) * cos(tmp0);
	tmp2 = sqrt(-c1) * sin(tmp0);
	lambdas[0] = -tmp1 + tmp2;
	lambdas[1] = -tmp1 - tmp2;
	lambdas[2] = 2 * tmp1;

	// compute the probability with eq. 48 for mu->e
	std::complex<double> Ame, e, num;
	double den;
	int j;
	Ame = 0.;
	for (int i = 0; i < 3; i++)
	{
		if (E < 0) j = 2 - i; // flip the definition for antineutrinos
		else j = i;
		e = exp(std::complex<double>(0, -L * lambdas[j] / 1.973e-10)); // km ev -> unitless
		num = lambdas[j] * Ttilde(1, 0) + Ttildesq(1, 0);
		den = 3 * sq(lambdas[j]) + c1;
		Ame += e * num / den;
	} // i, 3 (a in the paper)

	return std::norm(Ame);
}
double Pmue2(double a, double L, double E)
{
	// real pmns matrix
	Matrix3d U;
	U(0, 0) = c12 * c13;
	U(0, 1) = s12 * c13;
	U(0, 2) = s13;
	U(1, 0) = -s12 * c23 - c12 * s23 * s13;
	U(1, 1) = c12 * c23 - s12 * s23 * s13;
	U(1, 2) = s23 * c13;
	U(2, 0) = s12 * s23 - c12 * c23 * s13;
	U(2, 1) = -c12 * s23 - s12 * c23 * s13;
	U(2, 2) = c23 * c13;

	// mass matrix
	Matrix3d Msq = Matrix3d::Zero();
	Msq(1, 1) = Dmsq21;
	Msq(2, 2) = Dmsq31;

	// matter potential
	Matrix3d A = Matrix3d::Zero();
	A(0, 0) = a;

	// Hamiltonian flavor basis
	Matrix3cd H;
	H = (U * Msq * U.adjoint() + A) / (2 * E);

	// exponent
	H *= std::complex<double>(0, -L * eVsqkm_to_GeV);
	Matrix3cd e = H.exp();

	return std::norm(e(0, 1));
}
}
