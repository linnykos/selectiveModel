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
  Rcpp::NumericVector res = no_init(2*n-1);

  for(int i = 0; i < n-1; i++){
    res[2*i] = x[i];
    res[2*i+1] = (x[i] + x[i+1])/2;
  }

  res[2*n-2] = x[n-1];

  return res;
}

// from https://stackoverflow.com/questions/23849354/equivalent-of-which-function-in-rcpp
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
  Rcpp::LogicalVector vec2 = Rcpp::clone(vec);

  // remove singletons
  for(int i = 1; i < n-1; i++){
    vec2[i] = (vec[i] && (vec[i-1] || vec[i+1]));
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

// from https://stackoverflow.com/questions/30175104/how-to-effectively-combine-a-list-of-numericvectors-into-one-large-numericvector
// [[Rcpp::export]]
Rcpp::NumericVector unlist_native(const Rcpp::List & list) {
  int n = list.size();

  // Figure out the length of the output vector
  int total_length = 0;
  for (int i = 0; i < n; ++i)
    total_length += Rf_length(list[i]);

  // Allocate the vector
  Rcpp::NumericVector output = no_init(total_length);

  // Loop and fill
  int index = 0;
  for (int i = 0; i < n; ++i)
  {
    Rcpp::NumericVector el = list[i];
    std::copy(el.begin(), el.end(), output.begin() + index);

    // Update the index
    index += el.size();
  }

  return output;
}

// [[Rcpp::export]]
Rcpp::LogicalVector theta_in_matrix(const double & x,
                                    const Rcpp::NumericMatrix & mat) {
  int nrow = mat.nrow();
  Rcpp::LogicalVector result(1);
  result[0] = FALSE;

  for(int i = 0; i < nrow; i++){
    if((double) mat(i,0) <= x && x <= (double) mat(i,1)) {
      result[0] = TRUE;
      return result;
    }
  }

  return result;
}

// [[Rcpp::export]]
Rcpp::LogicalVector theta_in_all_matrix(const double & x,
                                        const Rcpp::List & list) {

  int n = list.size();
  Rcpp::LogicalVector result(1);
  result[0] = TRUE;

  for(int i = 0; i < n; i++){
    Rcpp::NumericMatrix mat = as<Rcpp::NumericMatrix>(list[i]);
    Rcpp::LogicalVector boolean = theta_in_matrix(x, mat);
    if(boolean[0] == FALSE) {
      result[0] = FALSE;
      return result;
    }
  }

  return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix intersect_intervals(const Rcpp::List & list){

  Rcpp::NumericVector vec = unlist_native(list);
  vec = unique_sort_native(vec);

  Rcpp::NumericVector vec2 = construct_midpoints(vec);
  int n = vec2.size();
  Rcpp::LogicalVector boolean(n);

  for(int i = 0; i < n; i++){
    Rcpp::LogicalVector tmp = theta_in_all_matrix(vec2[i], list);

    boolean[i] = tmp[0];
  }

  Rcpp::IntegerMatrix idx = consecutive_true(boolean);
  int nrow = idx.nrow();
  Rcpp::NumericMatrix result(nrow, 2);

  for(int i = 0; i < nrow; i++){
    result(i,0) = vec2[idx(i,0)-1];
    result(i,1) = vec2[idx(i,1)-1];
  }

  return result;
}
