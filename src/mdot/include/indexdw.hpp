#pragma once

#include "mdot/include/typedef.hpp"
#include <algorithm>
#include <iostream>

void norm_non_optimal(const array_of_s_type &array_of_s, dnum_t *norm_out) {
  *norm_out = 0;
  for (auto &itout : array_of_s)
    for (auto &itin : itout)
      // norm_out += *itin;
      std::cout << "hey";
  return;
}