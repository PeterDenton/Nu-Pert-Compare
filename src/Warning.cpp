#include <iostream>
#include <sstream>
#include <string>

#include "Warning.h"

int N_Warning = 0;
int N_IndexWarning = 0;
int N_EdgeCaseWarning = 0;
int N_ProbabilityWarning = 0;
int N_MaxCountWarning = 0;

Warning::Warning() 
: name("Warning"), message("<empty message>"), warn(true)
{
	Warning_number = ++N_Warning;
}
Warning::~Warning()
{
	if ((N_Warning <= 5 or (N_Warning <= 100 and N_Warning % 5 == 0) or (N_Warning <= 1000 and N_Warning % 100 == 0) or N_Warning % 1000 == 0)
		and warn)
		std::cerr << "Warning #" << Warning_number << ": " << name << ": " << message << std::endl;
}
IndexWarning::IndexWarning(int min, int max, double val)
: min(min), max(max), val(val)
{
	name = "IndexWarning";
	std::ostringstream tmp;
	tmp << "Index range = [" << min << ", " << max << "], " << val << " was provided.";
	message = tmp.str();
	N_IndexWarning++;
}
EdgeCaseWarning::EdgeCaseWarning(double min, double max, double val)
: min(min), max(max), val(val)
{
	name = "EdgeCaseWarning";
	std::ostringstream tmp;
	tmp << "Range = [" << min << ", " << max << "], " << val << " was returned.";
	message = tmp.str();
	N_EdgeCaseWarning++;
}
ProbabilityWarning::ProbabilityWarning(double *P)
{
	if (*P >= 0 and *P <= 1)
	{
		N_Warning--;
		warn = false;		
	}
	else
	{
		name = "ProbabilityWarning";
		std::ostringstream tmp;
		tmp << "Probability = " << *P << ", is outside of [0, 1]. Setting P to " << (*P < 0 ? 0 : 1) << ".";
		*P = (*P < 0 ? 0 : 1);
		message = tmp.str();
		N_ProbabilityWarning++;
	}
}
MaxCountWarning::MaxCountWarning(int count)
: count(count)
{
	name = "MaxCountWarning";
	std::ostringstream tmp;
	tmp << "Max count of " << count << " was reached.";
	message = tmp.str();
	N_MaxCountWarning++;
}
