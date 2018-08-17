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
  void c_intersect_basis(const Rcpp::NumericVector, const Rcpp::NumericVector, Rcpp::NumericVector);
  void print();
};

// constructor
Plane::Plane(Rcpp::NumericVector a_, Rcpp::NumericVector b_) {
  double l2norm = c_l2norm(a_);
  a = a_/l2norm;
  b = b_/l2norm;
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
    intercept[1] += a[i] * y[i];
  }

  a = vec;
  b = b - intercept;
}

void Plane::print(){
  Rcpp::Rcout << "a = " << a << std::endl;
  Rcpp::Rcout << "b = " << b << std::endl;
}

RCPP_MODULE(module){
  using namespace Rcpp;

  Rcpp::class_<Plane>( "Plane" )
  .constructor<Rcpp::NumericVector, Rcpp::NumericVector>("documentation for constructor")
  .field( "a", &Plane::a, "documentation for plane")
  .field( "b", &Plane::b, "documentation for plane")
  .method( "print", &Plane::print, "documentation for print")
  .method( "c_intersect_basis", &Plane::c_intersect_basis, "documentation")
  ;
}

// in R:
// zz = new(Plane, 1:5, 2); zz$print(); zz$c_intersect_basis(c(1:5)/2, 6:10, 11:15)
