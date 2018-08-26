#ifndef _INTERVAL_H
#define _INTERVAL_H

#include <Rcpp.h>
#include <math.h>
#include "plane.h"

double c_euclidean_to_radian(const Circle & circle,
                             const Rcpp::NumericVector & point);
double c_initial_theta(const Rcpp::NumericVector & y,
                       const Rcpp::NumericVector & v,
                       const Rcpp::NumericVector & w);

#endif
