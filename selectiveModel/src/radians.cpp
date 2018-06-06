
#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericVector unique_sort_native(Rcpp::NumericVector x) {
  Rcpp::NumericVector res = Rcpp::clone(x);
  res.sort(false);

  Rcpp::NumericVector::iterator it;
  it = std::unique (res.begin(), res.end());

  res.erase(it, res.end());

  return res;
}
