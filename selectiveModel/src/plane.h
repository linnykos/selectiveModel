#ifndef _PLANE_H
#define _PLANE_H

#include <Rcpp.h>
#include "radians.h"

class Plane{
public:
  Plane(Rcpp::NumericVector, Rcpp::NumericVector);
  Rcpp::NumericVector a;
  Rcpp::NumericVector b;
  void c_normalize();
  void c_intersect_basis(const Rcpp::NumericVector, const Rcpp::NumericVector, Rcpp::NumericVector);
  Rcpp::NumericVector c_point_on_plane();
  void print();
};

#endif
