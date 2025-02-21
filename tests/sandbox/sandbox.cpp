#include <algorithm>
#include <complex>
#include <iostream>
#include <map>
#include <string>
#include <tuple>

// #include "mdot/include/babel_type.hpp"

int main() {
  std::complex<double> a = {1, 2};
  std::complex<double> b = {1, 2};
  // znum_t c = a + b;
  double aj = nan(0);
  // std::cout << c.real() << "+j" << c.imag() << std::endl;

  std::string blabla;
  blabla = "hello";

  std::map<int, double> plus;
  plus[1] = 2.3;
  // op_index_t test = {2, 3};
  // dopbloc_t retest = {test,{1.,2.,3.,4.,5.}};
  // int pls = std::get<0>(test);

  std::pair<double, int> ls;
  ls = {1.1, 2};
}
