#include "tnorms.h"


// Minimum-Tnorm
double min_tnorm (double in1, double in2)
{
	return MIN(in1, in2);
}

// Produkt-Tnorm
double prod_tnorm (double in1, double in2)
{
	return in1 * in2;
}
 
// Lukasiewicz-Tnorm
double lukasiewicz_tnorm (double in1, double in2)
{
	return MAX(0, in1 + in2 - 1);
}

