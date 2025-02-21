#pragma once

#include <vector>

#include "mdot/include/babel_type.hpp"

namespace mdot {
std::vector<std::tuple<t_index_t, m_index_t, m_index_t>>
indices_dst_theta_no_gate(const std::vector<m_index_t> left_indices,
                          const std::vector<m_index_t> right_indices,
                          const bool conserve_left_right = false) {
  std::vector<std::tuple<t_index_t, m_index_t, m_index_t>>
      dst_about_indices_to_contract;
  for (auto &left_index : left_indices) {
    for (auto &right_index : right_indices) {
      if (std::get<2>(left_index) == std::get<0>(right_index)) {
        if (!conserve_left_right ||
            (conserve_left_right && std::get<0>(left_index) +
                                            std::get<1>(left_index) +
                                            std::get<1>(right_index) ==
                                        std::get<2>(right_index))) {
          dst_about_indices_to_contract.push_back(
              {{std::get<0>(left_index), std::get<1>(left_index),
                std::get<1>(right_index), std::get<2>(right_index)},
               left_index,
               right_index});
        }
      }
    }
  }
  std::sort(dst_about_indices_to_contract.begin(),dst_about_indices_to_contract.end());
  auto last = std::unique(dst_about_indices_to_contract.begin(),dst_about_indices_to_contract.end());
  dst_about_indices_to_contract.erase(last,dst_about_indices_to_contract.end());
  return dst_about_indices_to_contract;
}
} // namespace mdot
