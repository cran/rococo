#include <Rcpp.h>
#include "tnorms.h"
#include "rcor.h"

using namespace Rcpp;

/* produce shuffle of 1:n */

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

SEXP rcor_permtest(SEXP matx, SEXP maty, SEXP tests, SEXP ogamma,
		   double (*tnorm_fp)(double, double), SEXP alt,
		   SEXP storeValues)
{
    int i, cnt = 0;
    double c, d;
    NumericVector gamma(1);

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

    LogicalVector store(storeValues);

    int numGamma = (store[0] ? num_tests : 0);

    NumericVector permGamma(numGamma);

    for (i = 0; i < num_tests; i++)
    {
	shuffle_in_place(perm);	
		
	get_sums(mat_x, mat_y, perm, tnorm_fp, &c, &d);
	gamma[0] = (fabs(c + d) <= DBL_EPSILON) ? 0 : (c - d) / (c + d);

	gammaSqSum += (gamma[0] * gamma[0]);

	if (numGamma) permGamma[i] = gamma[0];

	/* It would, of course, be more efficient to define 'gamma' as a
	   simple double. Unfortunately, this allows the compiler to put
	   'gamma' into a register that could possibly have a different
	   precision, thereby, leading to false results when comparing two
	   numbers that should be considered as exactly equal. */

	switch (alternative)
	{
	    case 0: if (fabs(gamma[0]) >= fabs(old_gamma)) cnt++; break;
	    case 1: if (gamma[0] <= old_gamma)             cnt++; break;
	    case 2: if (gamma[0] >= old_gamma)             cnt++; break;
	}
    }

    if (numGamma)
	return List::create(_["cnt"] = cnt, _["ogamma"] = old_gamma,
			    _["H0sd"] = sqrt(gammaSqSum / num_tests),
	                    _["values"] = permGamma);
    else
	return List::create(_["cnt"] = cnt, _["ogamma"] = old_gamma,
			    _["H0sd"] = sqrt(gammaSqSum / num_tests));
}

RcppExport SEXP rcor_permtest_min(SEXP matx, SEXP maty, SEXP tests,
				  SEXP ogamma, SEXP alternative,
				  SEXP storeValues)
{
    return rcor_permtest(matx, maty, tests, ogamma, min_tnorm,
			 alternative, storeValues);
}

RcppExport SEXP rcor_permtest_prod(SEXP matx, SEXP maty, SEXP tests,
				   SEXP ogamma, SEXP alternative,
				   SEXP storeValues)
{
    return rcor_permtest(matx, maty, tests, ogamma, prod_tnorm,
			 alternative, storeValues);
}

RcppExport SEXP rcor_permtest_lukasiewicz(SEXP matx, SEXP maty, SEXP tests,
					  SEXP ogamma, SEXP alternative,
					  SEXP storeValues)
{
    return rcor_permtest(matx, maty, tests, ogamma, lukasiewicz_tnorm,
			 alternative, storeValues);
}
