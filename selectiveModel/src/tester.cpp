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

// c_intersect_circle_tester(c(1,1), 1, c(0,0), 2)
// c_intersect_circle_tester(c(1,1), c(1), c(0,0), c(2))
