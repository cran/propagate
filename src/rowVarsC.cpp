#include <Rcpp.h>
using namespace Rcpp;

RcppExport SEXP rowVarsC(SEXP mat) {
  NumericMatrix x(mat);
  int nrow = x.nrow(), ncol = x.ncol();
  NumericVector out(nrow);
  
  /* first pass => mean */ 
    for (int i = 0; i < nrow; i++) {
      double n = 0;
      double sum1 = 0, sum2 = 0;
      double mean1;    
      
      for (int j = 0; j < ncol; j++) {
        double val = x(i, j);
        if (!ISNA(val)) {
          n++;
          sum1 += val;
        }
      }
      mean1 = sum1/n;      
      
      /* second pass => variance */
        for (int j = 0; j < ncol; j++) {
          double val = x(i, j);
          if (!ISNA(val)) sum2 = sum2 + (val - mean1) * (val - mean1);  
        }
      
      double var = sum2/(n - 1);        
      out[i] = var;
    }
  
  return (out);
}
