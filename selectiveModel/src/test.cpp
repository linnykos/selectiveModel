#include <Rcpp.h>
using namespace Rcpp;

class Location{
public:
  Location();
  Location(int, int);
  int x;
  int y;
  void print();
};


//constructors
Location::Location() :x(0), y(0) { }
Location::Location(int xi, int yi) :x(xi), y(yi) { }

//print function
void Location::print(){
  Rcpp::Rcout << "x = " << x << std::endl;
  Rcpp::Rcout << "y = " << y << std::endl;
}

RCPP_MODULE(locationmodule){
  Rcpp::class_<Location>( "Location" )
  .constructor("documentation for default constructor")
  .constructor<int,int>("documentation for constructor")
  .field( "x", &Location::x, "documentation for x")
  .field( "y", &Location::y, "documentation for y")
  .method( "print", &Location::print, "documentation for print")
  ;
}
