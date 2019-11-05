#include <stdexcept>
#include <string>

#include "Probabilities.h"
#include "Parameters.h"
#include "ZS.h"
#include "DMP.h"
#include "Madrid.h"
#include "MP.h"
#include "AM.h"
#include "MF.h"
#include "AKS.h"
#include "Diagonalization.h"
#include "FL.h"
#include "AJLOS.h"
#include "Zeroth.h"
#include "Vacuum.h"
#include "AKT.h"
#include "DPZ.h"
#include "OS.h"

double Pmue(double a, double L, double E, expression m)
{
	switch (m)
	{
		case zs:		return ZS::Pmue(a, L, E);
		case dmp0:		return DMP::Pmue(a, L, E, 0);
		case dmp1:		return DMP::Pmue(a, L, E, 1);
		case madrid:	return Madrid::Pmue(a, L, E);
		case mp:		return MP::Pmue(a, L, E);
		case am2:		return AM::Pmue(a, L, E, 2.);
		case am52:		return AM::Pmue(a, L, E, 2.5);
		case mf:		return MF::Pmue(a, L, E);
		case aks:		return AKS::Pmue(a, L, E);
		case diag:		return Diagonalization::Pmue(a, L, E);
		case fl:		return FL::Pmue(a, L, E);
		case ajlos31:	return AJLOS::Pmue31(a, L, E);
		case ajlos48:	return AJLOS::Pmue48(a, L, E);
		case zeroth:	return Zeroth::Pmue(a, L, E);
		case vacuum:	return Vacuum::Pmue(a, L, E);
		case akt:		return AKT::Pmue(a, L, E);
		case dpz0:		return DPZ::Pmue(a, L, E, 0);
		case dpz2:		return DPZ::Pmue(a, L, E, 2);
		case os:		return OS::Pmue(a, L, E);
		case Exp:		return OS::Pmue2(a, L, E);

		default: throw std::domain_error("Model referenced not declared in Pmue(..., expression)");
	} // switch expression
}

std::string Name(expression m)
{
	switch(m)
	{
		case zs:		return "ZS";
		case dmp0:		return "DMP^0";
		case dmp1:		return "DMP^1";
		case madrid:	return "Madrid";
		case mp:		return "MP";
		case am2:		return "AM^2";
		case am52:		return "AM^{5/2}";
		case mf:		return "MF";
		case aks:		return "AKS";
		case diag:		return "Diag";
		case fl:		return "FL";
		case ajlos31:	return "AJLOS(31)";
		case ajlos48:	return "AJLOS(48)";
		case zeroth:	return "Zeroth";
		case vacuum:	return "Vacuum";
		case akt:		return "AKT";
		case dpz0:		return "DPZ^0";
		case dpz2:		return "DPZ^2";
		case os:		return "OS";
		case Exp:		return "Exp";

		default: throw std::domain_error("Model referenced not declared in Name(expression)");
	} // switch m
}

