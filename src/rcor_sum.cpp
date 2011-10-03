#include "rcor_sum.h"
#include "tnorms.h"

using namespace Rcpp;

void rcor_sum(NumericMatrix mat_x, NumericMatrix mat_y, IntegerVector perm_x, IntegerVector perm_y, double (*tnorm_fp)(double, double), double *sum, double* sumt)
{
	int i, j;
	int n = mat_x.nrow();
	
	for (i = 0; i < n; i++) 
	{
		for (j = i + 1; j < n; j++)
		{
			(*sum) += tnorm_fp(mat_x(perm_x(i), perm_x(j)), mat_y(perm_y(i), perm_y(j)));
			(*sum) += tnorm_fp(mat_x(perm_x(j), perm_x(i)), mat_y(perm_y(j), perm_y(i)));
			
			// Implizite Transponierung inkl. Zeilen/Spalten-Permutation
			(*sumt) += tnorm_fp(mat_x(perm_x(i), perm_x(j)), mat_y(perm_y(j), perm_y(i))); 
			(*sumt) += tnorm_fp(mat_x(perm_x(j), perm_x(i)), mat_y(perm_y(i), perm_y(j))); 
		}
	}

	return;
}



SEXP rcor_sum_min (SEXP matx, SEXP maty, SEXP permx, SEXP permy)
{
	NumericMatrix mat_x(matx); 
	NumericMatrix mat_y(maty); 
	
	IntegerVector perm_x(permx);
	IntegerVector perm_y(permy);
	
	double sum, sumt;
	rcor_sum(mat_x, mat_y, perm_x, perm_y, min_tnorm, &sum, &sumt);
	
	NumericVector res = NumericVector::create(Named("sum") = sum, Named("sumt") = sumt);
	
	return res;
}

SEXP rcor_sum_prod (SEXP matx, SEXP maty, SEXP permx, SEXP permy)
{
	NumericMatrix mat_x(matx); 
	NumericMatrix mat_y(maty); 
	
	IntegerVector perm_x(permx);
	IntegerVector perm_y(permy);
	
	double sum, sumt;
	rcor_sum(mat_x, mat_y, perm_x, perm_y, prod_tnorm, &sum, &sumt);
	
	NumericVector res = NumericVector::create(Named("sum") = sum, Named("sumt") = sumt);
	
	return res;
}

SEXP rcor_sum_lukasiewicz (SEXP matx, SEXP maty, SEXP permx, SEXP permy)
{
	NumericMatrix mat_x(matx); 
	NumericMatrix mat_y(maty); 
	
	IntegerVector perm_x(permx);
	IntegerVector perm_y(permy);
	
	double sum, sumt;
	rcor_sum(mat_x, mat_y, perm_x, perm_y, lukasiewicz_tnorm, &sum, &sumt);
	
	NumericVector res = NumericVector::create(Named("sum") = sum, Named("sumt") = sumt);
	
	return res;
}


