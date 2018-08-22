#ifndef _PLANE_H
#define _PLANE_H

#include <Rcpp.h>
#include "radians.h"

class Circle{
public:
  Circle(Rcpp::NumericVector, Rcpp::NumericVector);
  Rcpp::NumericVector radius;
  Rcpp::NumericVector center;
};

class Plane{
public:
  Plane(Rcpp::NumericVector, Rcpp::NumericVector);
  Rcpp::NumericVector a;
  Rcpp::NumericVector b;
  void c_normalize();
  void c_intersect_basis(const Rcpp::NumericVector, const Rcpp::NumericVector, Rcpp::NumericVector);
  Rcpp::NumericVector c_point_on_plane();
  double c_distance_point_to_plane(const Rcpp::NumericVector);
  void print();
};

#endif
