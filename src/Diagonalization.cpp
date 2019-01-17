#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

#include <cassert>
#include <cmath>
#include <complex>
#include <stdexcept>

#include "Diagonalization.h"
#include "Parameters.h"

using Eigen::Matrix3d;
using Eigen::Matrix3cd;
using Eigen::SelfAdjointEigenSolver;
using Eigen::Map;
using Eigen::RowMajor;

namespace Diagonalization
{
// indices are one indexed as we write them
Matrix3d Uij(int i, int j, double theta)
{
	assert(i < j);

	i -= 1;
	j -= 1;

	Matrix3d m = Matrix3d::Identity();

	double c, s;

	c = cos(theta);
	s = sin(theta);

	m(i, i) = c;
	m(j, j) = c;
	m(i, j) = s;
	m(j, i) = -s;

	return m;
}

Matrix3cd Uij(int i, int j, double theta, double delta)
{
	assert(i < j);

	i -= 1;
	j -= 1;

	Matrix3cd m = Matrix3cd::Identity();

	double c, s;
	std::complex<double> eid;

	c = cos(theta);
	s = sin(theta);
	eid = exp(std::complex<double>(0, delta));

	m(i, i) = c;
	m(j, j) = c;
	m(i, j) = s * conj(eid);
	m(j, i) = -s * eid;

	return m;
}

Matrix3cd UPMNS()
{
	Matrix3cd m = Matrix3cd::Identity();

	m *= Uij(2, 3, t23);
	m *= Uij(1, 3, t13, delta);
	m *= Uij(1, 2, t12);

	return m;
}

// 2E * Hamiltonian
Matrix3cd H2E(double a)
{
	Matrix3d msq, matter;
	Matrix3cd u;

	msq = Matrix3d::Zero();
	matter = Matrix3d::Zero();

	msq(1, 1) = Dmsq21;
	msq(2, 2) = Dmsq31;

	matter(0, 0) = a;

	u = UPMNS();

	return u * msq * u.adjoint() + matter;
}

double Pmue(double a, double L, double E)
{
	Matrix3cd h2e;

	h2e = H2E(a);

	std::complex<double> A = 0;

	SelfAdjointEigenSolver<Matrix3cd> es;

	es.compute(h2e);
	for (int i = 0; i < 3; i++)
	{
		A += conj(es.eigenvectors()(1, i)) * es.eigenvectors()(0, i) * exp(std::complex<double>(0, -eVsqkm_to_GeV * es.eigenvalues()(i) * L / (2 * E)));
	} // i, 3

	return std::norm(A);
}
}

