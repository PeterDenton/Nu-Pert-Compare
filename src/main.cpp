#include "Figures.h"
#include "Speed.h"
#include "Parameters.h"

int main()
{
	Recalc_Parameters();

	// Figures
	Precision();
	Probabilities();
	Peak_Precision();

//	Speed(); // takes O(few) minutes

	return 0;
}

