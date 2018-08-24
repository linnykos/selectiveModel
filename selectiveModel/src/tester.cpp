#include "plane.h"

// [[Rcpp::export]]
Rcpp::NumericMatrix c_intersect_circle_tester(Rcpp::NumericVector a,
                                              Rcpp::NumericVector b,
                                              Rcpp::NumericVector center,
                                              double radius){
  Plane plane(a, b);
  Circle circle(center, radius);
  Rcpp::NumericMatrix mat = plane.c_intersect_circle(circle);
  return mat;
}
