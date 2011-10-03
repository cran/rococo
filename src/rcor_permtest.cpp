#include <Rcpp.h>
#include "tnorms.h"

using namespace Rcpp;

/*
	Funktionen zur Berechnung der paarweisen T-Norm einschließlich der Summation
*/
#define ABS(a) ((a) < 0 ? -(a) : (a))

void shuffle_in_place(IntegerVector idx)
{
	int i, j, s;
	int len = idx.size();
	
	for (i = len - 1; i >= 1; i--)
	{
		j = rand() % (i + 1);
		if (j < 0) j += (i + 1);

		s = idx[i];
		idx[i] = idx[j];
		idx[j] = s;
	}
}





void get_sums (NumericMatrix mat_x, NumericMatrix mat_y, IntegerVector perm, double (*tnorm_fp)(double, double), double *sum, double *sum_t)
{
	int i, j, n = mat_x.nrow();
	*sum = *sum_t = 0.0;
	
	for (i = 0; i < n; i++)
	{
		for (j = i + 1; j < n; j++)
		{
			*sum += tnorm_fp(mat_x(i, j), mat_y(perm(i), perm(j)));
			*sum += tnorm_fp(mat_x(j, i), mat_y(perm(j), perm(i)));
			
			*sum_t += tnorm_fp(mat_x(i, j), mat_y(perm(j), perm(i)));
			*sum_t += tnorm_fp(mat_x(j, i), mat_y(perm(i), perm(j)));
		}
	}
}



SEXP rcor_permtest (SEXP matx, SEXP maty, SEXP tests, SEXP ogamma, double (*tnorm_fp)(double, double), SEXP alt)
{
    int i, cnt = 0;
    double c, d, gamma;

    NumericMatrix mat_x(matx); 
    NumericMatrix mat_y(maty); 

    int num_tests = IntegerVector(tests)[0];
    int alternative = IntegerVector(alt)[0];
    // 0 == two.sided, 1 == less, 2 == greater
    double old_gamma = NumericVector(ogamma)[0];
	
    IntegerVector perm(mat_x.nrow());
	
    for (i = 0; i < mat_x.nrow(); i++)
	perm[i] = i;

    double gammaSqSum = 0;

    for (i = 0; i < num_tests; i++)
    {
	shuffle_in_place(perm);	
		
	get_sums(mat_x, mat_y, perm, tnorm_fp, &c, &d);
	gamma = (fabs(c + d) <= 0.00001) ? 0 : (c - d) / (c + d);

	gammaSqSum += (gamma * gamma);

	switch (alternative)
	{
	    case 0: if (fabs(gamma) >= fabs(old_gamma)) cnt++; break;
	    case 1: if (gamma <= old_gamma) cnt++; break;
	    case 2: if (gamma >= old_gamma) cnt++; break;
	}
    }

    return List::create(_["cnt"] = cnt, _["ogamma"] = old_gamma,
	                _["H0sd"] = sqrt(gammaSqSum / num_tests));
}


/** Schnittstellen für den Permutationstest ************************************************/
RcppExport SEXP rcor_permtest_min (SEXP matx, SEXP maty, SEXP tests, SEXP ogamma, SEXP alternative)
{
	return rcor_permtest (matx, maty, tests, ogamma, min_tnorm, alternative);
}

RcppExport SEXP rcor_permtest_prod (SEXP matx, SEXP maty, SEXP tests, SEXP ogamma, SEXP alternative)
{
	return rcor_permtest (matx, maty, tests, ogamma, prod_tnorm, alternative);
}

RcppExport SEXP rcor_permtest_lukasiewicz (SEXP matx, SEXP maty, SEXP tests, SEXP ogamma, SEXP alternative)
{
	return rcor_permtest (matx, maty, tests, ogamma, lukasiewicz_tnorm, alternative);
}
/*******************************************************************************************/



/** Schnittstellen zur Berechnung des Rankkorrelationskoeffizienten ************************/
void rcor (SEXP matx, SEXP maty, double (*tnorm_fp)(double, double), double *sum, double *sum_t)
{
	int i;

	NumericMatrix mat_x(matx); 
	NumericMatrix mat_y(maty); 
	
	IntegerVector perm(mat_x.nrow());	
	
	for (i = 0; i < mat_x.nrow(); i++)
		perm[i] = i;
		
	get_sums(mat_x, mat_y, perm, tnorm_fp, sum, sum_t);	
}

RcppExport SEXP rcor_min (SEXP matx, SEXP maty)
{
	double c, d;
	
	rcor(matx, maty, min_tnorm, &c, &d);
	
	return List::create(_["c"] = c, _["d"] = d);
}

RcppExport SEXP rcor_prod (SEXP matx, SEXP maty)
{
	double c, d;
	
	rcor(matx, maty, prod_tnorm, &c, &d);

	return List::create(_["c"] = c, _["d"] = d);
}

RcppExport SEXP rcor_lukasiewicz (SEXP matx, SEXP maty)
{
	double c, d;
	
	rcor(matx, maty, lukasiewicz_tnorm, &c, &d);
	
	return List::create(_["c"] = c, _["d"] = d);
}

