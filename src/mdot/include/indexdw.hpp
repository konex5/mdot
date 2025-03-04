#pragma once

#include "mdot/include/babel_type.hpp"
#include <algorithm>

void norm_non_optimal(const std::vector<std::vector<dnum_t>> &array_of_s,
                      dnum_t &norm_out) {
  norm_out = 0;
  for (auto &itout : array_of_s)
#pragma omp
    for (auto &val : itout)
      norm_out += pow(val, 2);
  norm_out = sqrt(norm_out);
  return;
}

void norm_optimal(const std::vector<std::vector<dnum_t>> &array_of_s,
                  dnum_t &norm_out) {
  //   norm_out = 0;
  //   for (auto &itout : array_of_s)
  //     for (auto &val : itout)
  //       norm_out += val;
  //   return;
}