#pragma once

#include <vector>

#include "mdot/include/babel_type.hpp"

namespace mdot {

std::pair<std::vector<std::tuple<t_index_t, m_index_t, m_index_t>>,
          std::vector<std::tuple<t_index_t, m_index_t, m_index_t>>>
split_degenerate_indices(
    std::vector<std::tuple<t_index_t, m_index_t, m_index_t>> indices) {

  std::vector<std::tuple<t_index_t, m_index_t, m_index_t>> new_indices,
      dup_indices;
  std::vector<t_index_t> all_targets;
  for (auto &it : indices)
    all_targets.push_back(std::get<0>(it));

  for (std::size_t i = 0; i < all_targets.size(); i++) {
    if (std::find(all_targets.begin(), all_targets.end(), all_targets[i]) -
            all_targets.begin() ==
        static_cast<long>(i))
      new_indices.push_back(indices[i]);
    else
      dup_indices.push_back(indices[i]);
  }

  return {new_indices, dup_indices};
}

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
  std::sort(dst_about_indices_to_contract.begin(),
            dst_about_indices_to_contract.end());
  dst_about_indices_to_contract.erase(
      std::unique(dst_about_indices_to_contract.begin(),
                  dst_about_indices_to_contract.end()),
      dst_about_indices_to_contract.end());
  return dst_about_indices_to_contract;
}

static inline constexpr index_t internal_qn_sum(const index_t lhs,
                                                const index_t rhs) {
  return lhs + rhs;
}

static inline constexpr index_t internal_qn_sub(const index_t lhs,
                                                const index_t rhs) {
  return lhs - rhs;
}

std::vector<index_t> potential_middle_indices(
    std::vector<t_index_t> theta_indices, const int direction_right = -1
) {
  std::vector<index_t> middle_indices;
  if (direction_right == -2) {
      for (auto& theta_index : theta_indices) {
        middle_indices.push_back(internal_qn_sum(std::get<0>(theta_index),std::get<1>(theta_index)));
        middle_indices.push_back(internal_qn_sub(std::get<3>(theta_index),std::get<2>(theta_index)));
      } 
    } else if (direction_right == -1) {
      for (auto& theta_index : theta_indices) {
        middle_indices.push_back(internal_qn_sum(std::get<0>(theta_index),std::get<1>(theta_index)));
        middle_indices.push_back(internal_qn_sum(std::get<2>(theta_index),std::get<3>(theta_index)));
      } 
    } else if (direction_right == 1) {
      for (auto& theta_index : theta_indices) {
        middle_indices.push_back(internal_qn_sum(std::get<0>(theta_index),std::get<1>(theta_index)));
      }
    } else if (direction_right == 2) {
      for (auto& theta_index : theta_indices) {
        middle_indices.push_back(internal_qn_sub(std::get<3>(theta_index),std::get<2>(theta_index)));
      }
    } else if (direction_right == 3) {
      for (auto& theta_index : theta_indices) {
        middle_indices.push_back(internal_qn_sum(std::get<2>(theta_index),std::get<3>(theta_index)));
      }
    }
    std::sort(middle_indices.begin(), middle_indices.end());
    middle_indices.erase(std::unique(middle_indices.begin(), middle_indices.end()), middle_indices.end());
    return middle_indices;
}

} // namespace mdot
