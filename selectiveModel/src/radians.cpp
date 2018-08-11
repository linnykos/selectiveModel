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
bool theta_in_matrix(const double & x,
                     const Rcpp::NumericMatrix & mat) {
  int nrow = mat.nrow();

  for(int i = 0; i < nrow; i++){
    if(mat(i,0) <= x && x <= mat(i,1)) {
      return TRUE;
    }
  }

  return FALSE;
}

// [[Rcpp::export]]
Rcpp::IntegerVector which_native(const Rcpp::LogicalVector & x) {
  int nx = x.size();
  std::vector<int> y;
  y.reserve(nx);

  for(int i = 0; i < nx; i++) {
    if (x[i]) y.push_back(i+1);
  }

  return wrap(y);
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix consecutive_true(const Rcpp::LogicalVector & vec){
  int n = vec.size();
  Rcpp::LogicalVector vec2 = vec;

  // remove singletons
  for(int i = 1; i < n-1; i++){
    vec2(i) = (vec(i) && (vec(i-1) || vec(i+1)));
  }

  Rcpp::IntegerVector idx = which_native(vec2);
  Rcpp::IntegerVector idx_plus = idx + 1;
  Rcpp::IntegerVector idx_minus = idx - 1;

  Rcpp::IntegerVector start_position = setdiff(idx, idx_plus);
  Rcpp::IntegerVector end_position = setdiff(idx, idx_minus);

  int n_idx = start_position.size();
  Rcpp::IntegerMatrix mat(n_idx, 2);

  for(int i = 0; i < n_idx; i++){
    mat(i,0) = start_position(i);
    mat(i,1) = end_position(i);
  }

  return mat;
}

