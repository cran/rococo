#include <Rcpp.h>
#include "tnorms.h"
#include "rcor.h"

using namespace Rcpp;

void get_sums (NumericMatrix mat_x, NumericMatrix mat_y, IntegerVector perm,
	       double (*tnorm_fp)(double, double), double *sum, double *sum_t)
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


void rcor (SEXP matx, SEXP maty, double (*tnorm_fp)(double, double),
	   double *sum, double *sum_t)
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
