#include "plane.h"

// [[Rcpp::export]]
double c_l2norm(const Rcpp::NumericVector & vec){
  double val = sqrt(sum(pow(vec, 2)));
  return val;
}

// [[Rcpp::export]]
Rcpp::NumericVector c_quadratic(double a, double b, double c){
  double term = pow(b, 2) - 4*a*c;
  double tol = 1e-6;
  Rcpp::NumericVector value(2);
  std::fill(value.begin(), value.end(), Rcpp::NumericVector::get_na());

  if(fabs(term) < tol){
    term = -b/(2*a);
    value[0] = term;
  } else {
    value[0] = (-b-sqrt(term))/(2*a);
    value[1] = (-b+sqrt(term))/(2*a);
  }

  return(value);
}

//constructor
Circle::Circle(Rcpp::NumericVector center_, double radius_){
  center = center_;
  radius = radius_;
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

void Plane::c_intersect_basis(const Rcpp::NumericVector & y,
                              const Rcpp::NumericVector & v,
                              const Rcpp::NumericVector & w){
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

double Plane::c_distance_point_to_plane(const Rcpp::NumericVector & point){
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

  tmp = fabs(tmp/c_l2norm(a));
  return(tmp);
}

Rcpp::NumericMatrix Plane::c_intersect_circle(const Circle & circle){
  double dis = c_distance_point_to_plane(circle.center);
  double tol = 1e-6;

  Rcpp::NumericMatrix mat(2,2);
  std::fill(mat.begin(), mat.end(), Rcpp::NumericVector::get_na());

  Rcpp::NumericVector x(2);
  Rcpp::NumericVector y(2);
  std::fill(x.begin(), x.end(), Rcpp::NumericVector::get_na());
  std::fill(y.begin(), y.end(), Rcpp::NumericVector::get_na());

  if(dis > circle.radius + tol){ return(mat);}

  if(fabs(a[1]) < tol){
    x[0] = b[0]/a[0];
    Rcpp::NumericVector tmp1 = circle.center[0];
    double tmp2 = Rcpp::as<double>(x) - Rcpp::as<double>(tmp1);
    y[0] = circle.center[1] - sqrt(pow(circle.radius, 2) - pow(tmp2, 2));
    y[1] = circle.center[1] + sqrt(pow(circle.radius, 2) - pow(tmp2, 2));

  } else if (fabs(a[0]) < tol){
    x[0] = b[0]/a[1];
    Rcpp::NumericVector tmp1 = circle.center[1];
    double tmp2 = Rcpp::as<double>(x) - Rcpp::as<double>(tmp1);
    y[0] = circle.center[0] - sqrt(pow(circle.radius, 2) - pow(tmp2, 2));
    y[1] = circle.center[0] + sqrt(pow(circle.radius, 2) - pow(tmp2, 2));

  } else {
    Rcpp::Rcout << "here" << std::endl;
    Rcpp::Rcout << "plane.a inside = " << a << std::endl;
    Rcpp::NumericVector tmp1 = Rcpp::no_init(1);
    tmp1 = a[0];
    Rcpp::Rcout << "tmp1 = " << tmp1 << std::endl;
    double a1 = Rcpp::as<double>(tmp1);
    tmp1 = a[1];
    Rcpp::Rcout << "asdfasdfasdf" << std::endl;
    double a2 = Rcpp::as<double>(tmp1);
    tmp1 = circle.center[0];
    Rcpp::Rcout << "14123" << std::endl;
    double c1 = Rcpp::as<double>(tmp1);
    tmp1 = circle.center[1];
    double c2 = Rcpp::as<double>(tmp1);
    Rcpp::Rcout << "a1 = " << a1 << std::endl;
    Rcpp::Rcout << "a2 = " << a2 << std::endl;
    Rcpp::Rcout << "c1 = " << c1 << std::endl;
    Rcpp::Rcout << "c2 = " << c2 << std::endl;

    double a_ = 1 + pow(a1/a2, 2);
    double tmp = Rcpp::as<double>(b);
    double b_ = -2*(a1/a2)*(tmp/a2 - c2) - 2*c1;
    double c_ = -pow(circle.radius, 2) + pow(tmp/a2 - c2, 2) + pow(c1, 2);

    Rcpp::Rcout << "a_ = " << a_ << std::endl;
    Rcpp::Rcout << "b_ = " << b_ << std::endl;
    Rcpp::Rcout << "c_ = " << c_ << std::endl;

    Rcpp::NumericVector x = c_quadratic(a_, b_, c_);
    Rcpp::NumericVector y(2);
    double plane_b = Rcpp::as<double>(b);

    Rcpp::Rcout << "x = " << x << std::endl;
    Rcpp::Rcout << "y = " << y << std::endl;

    tmp1 = x[0];
    double x_double = Rcpp::as<double>(tmp1);
    y[0] = (plane_b - a1*x_double)/a2;

    Rcpp::Rcout << "x = " << x << std::endl;
    Rcpp::Rcout << "y = " << y << std::endl;

    tmp1 = x[1];
    x_double = Rcpp::as<double>(tmp1);
    Rcpp::LogicalVector tmp_bool = Rcpp::NumericVector::is_na(x[1]);
    if(!Rcpp::as<bool>(tmp_bool)){
      y[1] = (plane_b - a1*x_double)/a2;
    }

    Rcpp::Rcout << "x = " << x << std::endl;
    Rcpp::Rcout << "y = " << y << std::endl;

    double tol2 = 1e-9;
    mat(0,0) = x[0];
    mat(1,0) = y[0];

    Rcpp::LogicalVector tmp_bool_x = fabs(x[0] - x[1]) > tol2;
    Rcpp::LogicalVector tmp_bool_y = fabs(y[0] - y[1]) > tol2;

    if(!Rcpp::as<bool>(tmp_bool) & Rcpp::as<bool>(tmp_bool_x) & Rcpp::as<bool>(tmp_bool_y)){
      mat(1,0) = x[1];
      mat(1,1) = y[1];
    }
  }

  return(mat);
}

void Plane::print(){
  Rcpp::Rcout << "a = " << a << std::endl;
  Rcpp::Rcout << "b = " << b << std::endl;
}

RCPP_MODULE(module){
  Rcpp::class_<Circle>( "Circle" )
  .constructor<Rcpp::NumericVector, double>("documentation for constructor")
  .field( "center", &Circle::center, "documentation for circle")
  .field( "radius", &Circle::radius, "documentation for circle")
  ;

  Rcpp::class_<Plane>( "Plane" )
  .constructor<Rcpp::NumericVector, Rcpp::NumericVector>("documentation for constructor")
  .field( "a", &Plane::a, "documentation for plane")
  .field( "b", &Plane::b, "documentation for plane")
  .method( "print", &Plane::print, "documentation for print")
  .method( "c_intersect_basis", &Plane::c_intersect_basis, "documentation")
  .method( "c_point_on_plane", &Plane::c_point_on_plane, "documentation")
  .method( "c_distance_point_to_plane", &Plane::c_distance_point_to_plane, "documentation")
  .method( "c_intersect_circle", &Plane::c_intersect_circle, "documentation")
  ;
}
