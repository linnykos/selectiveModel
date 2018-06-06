
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
