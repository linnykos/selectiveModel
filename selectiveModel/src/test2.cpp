#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::LogicalVector func1(const Rcpp::NumericVector & x,
                          const Rcpp::List & list) {

  Rcout << "hi";
  Rcpp::LogicalVector result(1);
  result[0] = TRUE;

  return result;
}

double func_alt(double x){
  Rcout << "yolo";
  Rcout << x;
  return x;
}

// [[Rcpp::export]]
double func2(){

  // NumericVector vec = NumericVector::create(1, -1);
  double vec[] = {1, -1};
  // Rcpp::LogicalVector tmp;
  // Rcpp::LogicalVector tmp = func1(vec2, list);
  double tmp = func_alt(vec[1]);

  // Rcpp::NumericVector::iterator it;
  // for(it = vec.begin(); it < vec.end(); it++){
  //   tmp = func1(*it, list);
  // }

  return tmp;
}
