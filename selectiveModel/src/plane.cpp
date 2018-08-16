#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double c_l2norm(const Rcpp::NumericVector & vec){
  double val = sqrt(sum(pow(vec, 2)));
  return val;
}

class Plane{
public:
  Plane(Rcpp::NumericVector, Rcpp::NumericVector);
  Rcpp::NumericVector a;
  Rcpp::NumericVector b;
  void c_intersect_basis(Rcpp::NumericVector, Rcpp::NumericVector);
  void print();
};

// constructor
Plane::Plane(Rcpp::NumericVector a_, Rcpp::NumericVector b_) {
  double l2norm = c_l2norm(a_);
  a = a_/l2norm;
  b = b_/l2norm;
}

void Plane::print(){
  Rcpp::Rcout << "a = " << a << std::endl;
  Rcpp::Rcout << "b = " << b << std::endl;
}

RCPP_MODULE(planemodule){
  Rcpp::class_<Plane>( "Plane" )
  .constructor<Rcpp::NumericVector, Rcpp::NumericVector>("documentation for constructor")
  .field( "a", &Plane::a, "documentation for plane")
  .field( "b", &Plane::b, "documentation for plane")
  .method( "print", &Plane::print, "documentation for print")
  ;
}

// in R:
// zz = new(Plane, 1:5, 2); zz$print()
