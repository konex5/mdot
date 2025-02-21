
#include <map>
#include <string>
#include <vector>

#include "mdot/include/babel_type.hpp"

dopbloc_t real_single_operator(std::string name, std::string qbasis) {
  std::map<std::string, typename std::map<std::string, dopbloc_t>>
      real_single_operators;
  std::map<std::string, typename std::map<std::string, double>> normalization;
  // C mode
  normalization["sh-id"]["sh-none"] = 1. / sqrt(2);
  real_single_operators["sh-id"]["sh-none"][{0, 0}] = {1, 0, 0, 1};
  //
  real_single_operators["sh-id"]["sh-u1"][{0, 0}] = {1};
  real_single_operators["sh-id"]["sh-u1"][{1, 1}] = {1};
  ///
  normalization["sh-sp"]["sh-none"] = 0.;
  real_single_operators["sh-sp"]["sh-none"][{0, 0}] = {0, 1, 0, 0};
  //
  real_single_operators["sh-sp"]["sh-u1"][{0, 1}] = {1};
  //

  return real_single_operators[name][qbasis];
}
