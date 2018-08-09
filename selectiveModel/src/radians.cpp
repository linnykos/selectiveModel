#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericVector unique_sort_native(const Rcpp::NumericVector & x) {
  Rcpp::NumericVector res = Rcpp::clone(x);
  res.sort(false);

  Rcpp::NumericVector::iterator it;
  it = std::unique (res.begin(), res.end());

  res.erase(it, res.end());

  return res;
}

// [[Rcpp::export]]
Rcpp::NumericVector construct_midpoints(const Rcpp::NumericVector & x) {
  int n = x.size();
  NumericVector res = no_init(2*n-1);

  for(int i = 0; i < n-1; i++){
    res[2*i] = x[i];
    res[2*i+1] = (x[i] + x[i+1])/2;
  }

  res[2*n-2] = x[n-1];

  return res;
}

// [[Rcpp::export]]
bool theta_in_matrix(const double x,
                     const Rcpp::NumericMatrix mat) {
  int nrow = mat.nrow();

  for(int i = 0; i < nrow; i++){
    if(mat(i,0) <= x && x <= mat(i,1)) {
      return TRUE;
    }
  }

  return FALSE;
}

// [[Rcpp::export()]]
Rcpp::List i2 (Rcpp::NumericMatrix x,
               double t
) {
  int n = x.nrow() ;
  Rcpp::NumericVector y(n) ;
  Rcpp::List ret ;
  for (int it = 0 ; it < n ; it++) {
    if (x(0,it) <= t) {
      y(it) = sqrt(pow(x(0,it) - 1.3, 4)) ;
    } else {
      y(it) = x(0,it) * 2 ;
    }
  }
  ret["x"] = x ; ret["y"] = y ;
  return(ret) ;
}
