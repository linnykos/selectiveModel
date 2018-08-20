#include "plane.h"

// [[Rcpp::export]]
double c_l2norm(const Rcpp::NumericVector & vec){
  double val = sqrt(sum(pow(vec, 2)));
  return val;
}

// constructor
Plane::Plane(Rcpp::NumericVector a_, Rcpp::NumericVector b_) {
  a = a_;
  b = b_;
  c_normalize();
}

void Plane::c_normalize() {
  double l2norm = c_l2norm(a);
  a = a/l2norm;
  b = b/l2norm;
}

void Plane::c_intersect_basis(const Rcpp::NumericVector y,
                              const Rcpp::NumericVector v,
                              const Rcpp::NumericVector w){
  Rcpp::NumericVector vec = Rcpp::NumericVector::create(0, 0);
  Rcpp::NumericVector intercept = Rcpp::NumericVector::create(0);
  int len = y.length();
  for(int i = 0; i < len; i++){
    vec[0] += a[i] * v[i];
    vec[1] += a[i] * w[i];
    intercept[0] += a[i] * y[i];
  }

  a = vec;
  b = b - intercept;
  c_normalize();
}

Rcpp::NumericVector Plane::c_point_on_plane(){
  int d = a.length();
  Rcpp::NumericVector vec(d);
  std::fill(vec.begin(), vec.end(), 1);

  Rcpp::LogicalVector boolean = (a != 0);
  Rcpp::IntegerVector idx = c_which_native(boolean);
  idx = idx[0] - 1;

  double tmp = 0;
  Rcpp::NumericVector tmp2;
  double sum = 0;
  for(int i = 0; i < d; i++){
    if(Rcpp::as<int>(idx) != i){
      tmp2 = a[i]*vec[i];
      sum = Rcpp::as<double>(tmp2);
      tmp += sum;
    }
  }

  Rcpp::NumericVector a_intercept = a[idx];
  tmp = (Rcpp::as<double>(b) - tmp)/(Rcpp::as<double>(a_intercept));
  vec[idx] = tmp;

  return(vec);
}

double Plane::c_distance_point_to_plane(const Rcpp::NumericVector point){
  Rcpp::NumericVector x = c_point_on_plane();
  int len = a.length();
  double tmp = 0;
  Rcpp::NumericVector tmp2;
  double sum = 0;

  for(int i = 0; i < len; i++){
    tmp2 = a[i] * (point[i] - x[i]);
    sum = Rcpp::as<double>(tmp2);
    tmp += sum;
  }

  tmp = abs(tmp/c_l2norm(a));
  return(tmp);
}

void Plane::print(){
  Rcpp::Rcout << "a = " << a << std::endl;
  Rcpp::Rcout << "b = " << b << std::endl;
}

RCPP_MODULE(module){
  Rcpp::class_<Plane>( "Plane" )
  .constructor<Rcpp::NumericVector, Rcpp::NumericVector>("documentation for constructor")
  .field( "a", &Plane::a, "documentation for plane")
  .field( "b", &Plane::b, "documentation for plane")
  .method( "print", &Plane::print, "documentation for print")
  .method( "c_intersect_basis", &Plane::c_intersect_basis, "documentation")
  .method( "c_point_on_plane", &Plane::c_point_on_plane, "documentation")
  .method( "c_distance_point_to_plane", &Plane::c_distance_point_to_plane, "documentation")
  ;
}

// in R:
// zz = new(Plane, 1:5, 2); zz$print(); zz$c_intersect_basis(c(1:5)/2, 6:10, 11:15); zz$print(); zz$c_point_on_plane(); zz$c_distance_point_to_plane(c(0,0))
