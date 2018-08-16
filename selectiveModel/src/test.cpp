#include <Rcpp.h>
using namespace Rcpp;

class Test {
public:
  Test(int x): x_(x) {}
  int getValue() { return x_; }
  void addValue(int y) { x_ += y; }
  void merge(const Test& rhs) { x_ += rhs.x_; }
private:
  int x_;
};

using namespace Rcpp;

RCPP_EXPOSED_CLASS(Test)
  RCPP_MODULE(mod_test) {

    class_<Test>("Test")

    .constructor<int>("sets initial value")

    .method("getValue", &Test::getValue, "Returns the value")
    .method("addValue", &Test::addValue, "Adds a value")
    .method("merge", &Test::merge, "Merges another Test into this object")
    ;
  }
