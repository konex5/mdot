
#include <map>
#include <string>
#include <vector>

#include "mdot/include/babel_type.hpp"

std::pair<dopbloc_t, dnum_t> single_operator_real(std::string name,
                                                  std::string qbasis) {
  std::map<std::string, typename std::map<std::string, dopbloc_t>>
      single_operators_real;
  std::map<std::string, typename std::map<std::string, double>> normalization;
  // C mode
  normalization["sh-id"]["sh-none"] = 1. / sqrt(2);
  single_operators_real["sh-id"]["sh-none"][{0, 0}] = {1, 0, 0, 1};
  //
  single_operators_real["sh-id"]["sh-u1"][{0, 0}] = {1};
  single_operators_real["sh-id"]["sh-u1"][{1, 1}] = {1};
  ///
  normalization["sh-sp"]["sh-none"] = 0.;
  single_operators_real["sh-sp"]["sh-none"][{0, 0}] = {0, 1, 0, 0};
  //
  single_operators_real["sh-sp"]["sh-u1"][{0, 1}] = {1};
  ///
  normalization["sh-sm"]["sh-none"] = 0.;
  single_operators_real["sh-sm"]["sh-none"][{0, 0}] = {0, 0, 1, 0};
  //
  single_operators_real["sh-sm"]["sh-u1"][{1, 0}] = {1};
  //
  normalization["sh-sz"]["sh-none"] = 1. / sqrt(2);
  single_operators_real["sh-sz"]["sh-none"][{0, 0}] = {1, 0, 0, -1};
  //
  single_operators_real["sh-sz"]["sh-u1"][{0, 0}] = {1};
  single_operators_real["sh-sz"]["sh-u1"][{1, 1}] = {-1};
  //
  normalization["sh-sx"]["sh-none"] = 1. / sqrt(2);
  single_operators_real["sh-sx"]["sh-none"][{0, 0}] = {0, 1, 1, 0};
  //
  single_operators_real["sh-sx"]["sh-u1"][{0, 1}] = {1};
  single_operators_real["sh-sx"]["sh-u1"][{1, 0}] = {1};
  /// ///
  normalization["so-id"]["sh-none"] = 1. / sqrt(3);
  single_operators_real["so-id"]["so-none"][{0, 0}] = {1, 0, 0, 0, 1,
                                                       0, 0, 0, 1};
  //
  single_operators_real["so-id"]["so-u1"][{0, 0}] = {1};
  single_operators_real["so-id"]["so-u1"][{1, 1}] = {1};
  single_operators_real["so-id"]["so-u1"][{2, 2}] = {1};
  ///
  normalization["so-sp"]["so-none"] = 0.;
  single_operators_real["so-sp"]["so-none"][{0, 0}] = {0, 1, 0, 0, 0,
                                                       1, 0, 0, 0};
  //
  single_operators_real["so-sp"]["so-u1"][{1, 0}] = {1};
  single_operators_real["so-sp"]["so-u1"][{2, 1}] = {1};
  ///
  normalization["so-sm"]["so-none"] = 0.;
  single_operators_real["so-sm"]["so-none"][{0, 0}] = {0, 0, 0, 1, 0,
                                                       0, 0, 1, 0};
  //
  single_operators_real["so-sm"]["so-u1"][{0, 1}] = {1};
  single_operators_real["so-sm"]["so-u1"][{1, 2}] = {1};
  ///
  normalization["so-sz"]["so-none"] = 1. / sqrt(2);
  single_operators_real["so-sz"]["so-none"][{0, 0}] = {1, 0, 0, 0, 0,
                                                       0, 0, 0, -1};
  //
  single_operators_real["so-sz"]["so-u1"][{0, 0}] = {1};
  single_operators_real["so-sz"]["so-u1"][{2, 2}] = {-1};
  ///
  normalization["so-sx"]["so-none"] = 1. / sqrt(2);
  single_operators_real["so-sx"]["so-none"][{0, 0}] = {0, 1, 0, 1, 0,
                                                       1, 0, 1, 0};
  //
  single_operators_real["so-sx"]["so-u1"][{1, 0}] = {1};
  single_operators_real["so-sx"]["so-u1"][{0, 1}] = {1};
  single_operators_real["so-sx"]["so-u1"][{2, 1}] = {1};
  single_operators_real["so-sx"]["so-u1"][{1, 2}] = {1};
  ///
  ///
  normalization["so-sz^2"]["so-none"] = 1. / sqrt(2);
  single_operators_real["so-sz^2"]["so-none"][{0, 0}] = {1, 0, 0, 0, 0,
                                                         0, 0, 0, 1};
  //
  single_operators_real["so-sz^2"]["so-u1"][{0, 0}] = {1};
  single_operators_real["so-sz^2"]["so-u1"][{2, 2}] = {1};

  return {single_operators_real[name][qbasis], normalization[name][qbasis]};
}

std::pair<zopbloc_t, dnum_t> single_operator_cplx(std::string name,
                                                  std::string qbasis) {
  std::map<std::string, typename std::map<std::string, zopbloc_t>>
      single_operators_cplx;
  std::map<std::string, typename std::map<std::string, double>> normalization;
  // C mode
  normalization["sh-id-cplx"]["sh-none"] = 1. / sqrt(2);
  single_operators_cplx["sh-id-cplx"]["sh-none"][{0, 0}] = {
      {1, 0}, {0, 0}, {0, 0}, {1, 0}};
  //
  single_operators_cplx["sh-id-cplx"]["sh-u1"][{0, 0}] = {{1, 0}};
  single_operators_cplx["sh-id-cplx"]["sh-u1"][{1, 1}] = {{1, 0}};
  //
  normalization["sh-sy"]["sh-none"] = 1. / sqrt(2);
  single_operators_cplx["sh-sy"]["sh-none"][{0, 0}] = {
      {0, 0}, {0, 1}, {0, -1}, {0, 0}};
  //
  single_operators_cplx["sh-sy"]["sh-u1"][{0, 1}] = {{0, 1}};
  single_operators_cplx["sh-sy"]["sh-u1"][{1, 0}] = {{0, -1}};
  //////
  normalization["so-id-cplx"]["sh-none"] = 1. / sqrt(3);
  single_operators_cplx["so-id-cplx"]["sh-none"][{0, 0}] = {
      {1, 0}, {0, 0}, {0, 0}, {0, 0}, {1, 0}, {0, 0}, {0, 0}, {0, 0}, {1, 0}};
  //
  single_operators_cplx["so-id-cplx"]["sh-u1"][{0, 0}] = {{1, 0}};
  single_operators_cplx["so-id-cplx"]["sh-u1"][{1, 1}] = {{1, 0}};
  single_operators_cplx["so-id-cplx"]["sh-u1"][{2, 2}] = {{1, 0}};
  ///
  normalization["so-sy"]["sh-none"] = 1. / sqrt(3);
  single_operators_cplx["so-sy"]["sh-none"][{0, 0}] = {
      {0, 0}, {0, 1}, {0, 0}, {0, -1}, {0, 0}, {0, 1}, {0, 0}, {0, -1}, {0, 0}};
  //
  single_operators_cplx["so-sy"]["sh-u1"][{0, 0}] = {{1, 0}};
  single_operators_cplx["so-sy"]["sh-u1"][{1, 1}] = {{1, 0}};
  single_operators_cplx["so-sy"]["sh-u1"][{2, 2}] = {{1, 0}};
  //
  return {single_operators_cplx[name][qbasis], normalization[name][qbasis]};
}