#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector func_alt(double x){
  Rcout << "yolo";
  Rcout << x;
  return wrap(x);
}

// [[Rcpp::export]]
NumericVector func2(){

  NumericVector vec = NumericVector::create(1, -1);
  Rcout << vec;
  // Rcpp::LogicalVector tmp;
  // Rcpp::LogicalVector tmp = func1(vec2, list);
  NumericVector tmp = func_alt(vec[1]);

  // Rcpp::NumericVector::iterator it;
  // for(it = vec.begin(); it < vec.end(); it++){
  //   tmp = func1(*it, list);
  // }

  return tmp;
}
