#ifndef _RCOR_SUM_H
#define _RCOR_SUM_H

#include <Rcpp.h>

RcppExport SEXP rcor_sum_min(SEXP matx, SEXP maty, SEXP permx, SEXP permy);
RcppExport SEXP rcor_sum_prod(SEXP matx, SEXP maty, SEXP permx, SEXP permy);
RcppExport SEXP rcor_sum_lukasiewicz(SEXP matx, SEXP maty, SEXP permx, SEXP permy);

#endif
