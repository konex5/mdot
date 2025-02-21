#pragma once

#include <vector>

#include "mdot/include/babel_type.hpp"

namespace mdot {

std::pair<std::vector<std::tuple<t_index_t, m_index_t, m_index_t>>,
          std::vector<std::tuple<t_index_t, m_index_t, m_index_t>>>
split_degenerate_indices(
    const std::vector<std::tuple<t_index_t, m_index_t, m_index_t>> indices) {

  std::vector<std::tuple<t_index_t, m_index_t, m_index_t>> new_indices,
      dup_indices;
  std::vector<t_index_t> all_targets;
  for (auto &it : indices)
    all_targets.push_back(std::get<0>(it));

  for (size_t i = 0; i < all_targets.size(); i++) {
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
        if ((!conserve_left_right) ||
            (conserve_left_right &&
             std::get<0>(left_index) +
                     static_cast<index_t>(std::get<1>(left_index)) +
                     static_cast<index_t>(std::get<1>(right_index)) ==
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

std::pair<std::vector<std::tuple<t_index_t, t_index_t, g_index_t>>,
          std::vector<std::tuple<t_index_t, t_index_t, g_index_t>>>
split_degenerate_indices_with_gate(
    const std::vector<std::tuple<t_index_t, t_index_t, g_index_t>> indices) {

  std::vector<std::tuple<t_index_t, t_index_t, g_index_t>> new_indices,
      dup_indices;
  std::vector<t_index_t> all_targets;
  for (auto &it : indices)
    all_targets.push_back(std::get<0>(it));

  for (size_t i = 0; i < all_targets.size(); i++) {
    if (std::find(all_targets.begin(), all_targets.end(), all_targets[i]) -
            all_targets.begin() ==
        static_cast<long>(i))
      new_indices.push_back(indices[i]);
    else
      dup_indices.push_back(indices[i]);
  }

  return {new_indices, dup_indices};
}

std::vector<std::tuple<t_index_t, t_index_t, g_index_t>>
indices_dst_theta_with_gate(const std::vector<t_index_t> theta_indices,
                            const std::vector<g_index_t> gate_indices,
                            const bool conserve_left_right = false) {
  std::vector<std::tuple<t_index_t, t_index_t, g_index_t>>
      dst_about_indices_to_contract;
  for (auto &theta_index : theta_indices) {
    for (auto &gate_index : gate_indices) {
      if (std::get<1>(theta_index) == std::get<0>(gate_index) &&
          std::get<2>(theta_index) == std::get<3>(gate_index)) {
        if ((!conserve_left_right) ||
            (conserve_left_right &&
             std::get<0>(theta_index) +
                     static_cast<index_t>(std::get<1>(gate_index)) +
                     static_cast<index_t>(std::get<2>(gate_index)) ==
                 std::get<3>(theta_index))) {
          dst_about_indices_to_contract.push_back(
              {{std::get<0>(theta_index), std::get<1>(gate_index),
                std::get<2>(gate_index), std::get<3>(theta_index)},
               theta_index,
               gate_index});
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

std::vector<index_t>
potential_middle_indices(const std::vector<t_index_t> theta_indices,
                         const int direction_right = -1) {
  std::vector<index_t> middle_indices;
  if (direction_right == -2) {
    for (auto &theta_index : theta_indices) {
      middle_indices.push_back(
          internal_qn_sum(std::get<0>(theta_index),
                          static_cast<index_t>(std::get<1>(theta_index))));
      middle_indices.push_back(
          internal_qn_sub(std::get<3>(theta_index),
                          static_cast<index_t>(std::get<2>(theta_index))));
    }
  } else if (direction_right == -1) {
    for (auto &theta_index : theta_indices) {
      middle_indices.push_back(
          internal_qn_sum(std::get<0>(theta_index),
                          static_cast<index_t>(std::get<1>(theta_index))));
      middle_indices.push_back(
          internal_qn_sum(static_cast<index_t>(std::get<2>(theta_index)),
                          std::get<3>(theta_index)));
    }
  } else if (direction_right == 1) {
    for (auto &theta_index : theta_indices) {
      middle_indices.push_back(
          internal_qn_sum(std::get<0>(theta_index),
                          static_cast<index_t>(std::get<1>(theta_index))));
    }
  } else if (direction_right == 2) {
    for (auto &theta_index : theta_indices) {
      middle_indices.push_back(
          internal_qn_sub(std::get<3>(theta_index),
                          static_cast<index_t>(std::get<2>(theta_index))));
    }
  } else if (direction_right == 3) {
    for (auto &theta_index : theta_indices) {
      middle_indices.push_back(
          internal_qn_sum(static_cast<index_t>(std::get<2>(theta_index)),
                          std::get<3>(theta_index)));
    }
  }
  std::sort(middle_indices.begin(), middle_indices.end());
  middle_indices.erase(
      std::unique(middle_indices.begin(), middle_indices.end()),
      middle_indices.end());
  return middle_indices;
}

std::pair<std::vector<std::pair<index_t, t_index_t>>,
          std::vector<std::pair<index_t, std::vector<t_index_t>>>>
degeneracy_in_theta(const std::vector<t_index_t> theta_indices,
                    const std::vector<index_t> middle_indices,
                    const int direction_right = -1) {
  std::vector<std::pair<index_t, t_index_t>> nondeg;
  std::vector<std::pair<index_t, std::vector<t_index_t>>> degenerate;

  index_t middle_length = static_cast<index_t>(middle_indices.size());

  if (direction_right == -2) {
    for (index_t j = 0; j < middle_length; j++) {
      std::vector<t_index_t> tmp;
      for (auto &theta_index : theta_indices) {
        if (internal_qn_sum(std::get<0>(theta_index),
                            std::get<1>(theta_index)) == middle_indices[j] &&
            internal_qn_sub(std::get<3>(theta_index),
                            std::get<2>(theta_index)) == middle_indices[j])
          tmp.push_back(theta_index);
      }
      if (tmp.size() == 1)
        nondeg.push_back({j, tmp[0]});
      else
        degenerate.push_back({j, tmp});
    }
  } else if (direction_right == -1) {
    for (index_t j = 0; j < middle_length; j++) {
      std::vector<t_index_t> tmp;
      for (auto &theta_index : theta_indices) {
        if (internal_qn_sum(std::get<0>(theta_index),
                            std::get<1>(theta_index)) == middle_indices[j] &&
            internal_qn_sum(std::get<2>(theta_index),
                            std::get<3>(theta_index)) == middle_indices[j])
          tmp.push_back(theta_index);
      }
      if (tmp.size() == 1)
        nondeg.push_back({j, tmp[0]});
      else
        degenerate.push_back({j, tmp});
    }
  } else if (direction_right == 1) {
    for (index_t j = 0; j < middle_length; j++) {
      std::vector<t_index_t> tmp;
      for (auto &theta_index : theta_indices) {
        if (internal_qn_sum(std::get<0>(theta_index),
                            std::get<1>(theta_index)) == middle_indices[j])
          tmp.push_back(theta_index);
      }
      if (tmp.size() == 1)
        nondeg.push_back({j, tmp[0]});
      else
        degenerate.push_back({j, tmp});
    }
  } else if (direction_right == 2) {
    for (index_t j = 0; j < middle_length; j++) {
      std::vector<t_index_t> tmp;
      for (auto &theta_index : theta_indices) {
        if (internal_qn_sub(std::get<3>(theta_index),
                            std::get<2>(theta_index)) == middle_indices[j])
          tmp.push_back(theta_index);
      }
      if (tmp.size() == 1)
        nondeg.push_back({j, tmp[0]});
      else
        degenerate.push_back({j, tmp});
    }
  } else if (direction_right == 3) {
    for (index_t j = 0; j < middle_length; j++) {
      std::vector<t_index_t> tmp;
      for (auto &theta_index : theta_indices) {
        if (internal_qn_sum(std::get<2>(theta_index),
                            std::get<3>(theta_index)) == middle_indices[j])
          tmp.push_back(theta_index);
      }
      if (tmp.size() == 1)
        nondeg.push_back({j, tmp[0]});
      else
        degenerate.push_back({j, tmp});
    }
  }

  return {nondeg, degenerate};
}

std::vector<std::tuple<index_t, index_t,
                       typename std::vector<std::tuple<index_t, index_small_t>>,
                       typename std::vector<index_t>,
                       typename std::vector<std::pair<index_t, index_small_t>>,
                       typename std::vector<std::tuple<index_small_t, index_t>>,
                       typename std::vector<index_t>,
                       typename std::vector<std::pair<index_small_t, index_t>>>>
slices_degenerate_blocs(
    dtbloc_t theta_blocs,
    const std::vector<std::pair<index_t, std::vector<t_index_t>>>
        degenerate_list) {
  //
  std::vector<
      std::tuple<index_t, index_t,
                 typename std::vector<std::tuple<index_t, index_small_t>>,
                 typename std::vector<index_t>,
                 typename std::vector<std::pair<index_t, index_small_t>>,
                 typename std::vector<std::tuple<index_small_t, index_t>>,
                 typename std::vector<index_t>,
                 typename std::vector<std::pair<index_small_t, index_t>>>>
      subnewsize;
  for (auto &degenerate_indices : degenerate_list) {
    std::vector<std::tuple<index_t, index_small_t>> left__loc_basis;
    std::vector<std::tuple<index_small_t, index_t>> right_loc_basis;
    std::vector<std::pair<index_t, index_small_t>> left__loc_dim;
    std::vector<std::pair<index_small_t, index_t>> right_loc_dim;
    index_t total_left__dim = 0;
    index_t total_right_dim = 0;
    std::vector<index_t> left__loc_off, right_loc_off;
    // define a local basis
    for (auto &key_index : degenerate_indices.second) {
      left__loc_basis.push_back(
          {std::get<0>(key_index), std::get<1>(key_index)});
      right_loc_basis.push_back(
          {std::get<2>(key_index), std::get<3>(key_index)});
    }
    std::sort(left__loc_basis.begin(), left__loc_basis.end());
    left__loc_basis.erase(
        std::unique(left__loc_basis.begin(), left__loc_basis.end()),
        left__loc_basis.end());
    std::sort(right_loc_basis.begin(), right_loc_basis.end());
    right_loc_basis.erase(
        std::unique(right_loc_basis.begin(), right_loc_basis.end()),
        right_loc_basis.end());
    // find the local dim corresponding to left_loc_basis and right_loc_basis
    left__loc_dim.resize(left__loc_basis.size());
    right_loc_dim.resize(right_loc_basis.size());
    for (auto &key_index : degenerate_indices.second) {
      auto dims = theta_blocs[key_index].first;
      std::tuple<index_t, index_small_t> left_val = {std::get<0>(key_index),
                                                     std::get<1>(key_index)};
      auto i =
          std::find(left__loc_basis.begin(), left__loc_basis.end(), left_val) -
          left__loc_basis.begin();
      left__loc_dim[i] = {std::get<0>(dims), std::get<1>(dims)};
      std::tuple<index_small_t, index_t> right_val = {std::get<2>(key_index),
                                                      std::get<3>(key_index)};
      auto j =
          std::find(right_loc_basis.begin(), right_loc_basis.end(), right_val) -
          right_loc_basis.begin();
      right_loc_dim[j] = {std::get<2>(dims), std::get<3>(dims)};
    }
    // totdim
    for (auto &v : left__loc_dim)
      total_left__dim += std::get<0>(v) * static_cast<index_t>(std::get<1>(v));
    for (auto &v : right_loc_dim)
      total_right_dim += static_cast<index_t>(std::get<0>(v)) * std::get<1>(v);
    // offsets
    left__loc_off.push_back(0);
    for (size_t i = 1; i < left__loc_dim.size(); i++) {
      index_t acc = 0;
      for (size_t j = 0; j < i; j++)
        acc += std::get<0>(left__loc_dim[j]) *
               static_cast<index_t>(std::get<1>(left__loc_dim[j]));
      left__loc_off.push_back(acc);
    }
    right_loc_off.push_back(0);
    for (size_t i = 1; i < right_loc_dim.size(); i++) {
      index_t acc = 0;
      for (size_t j = 0; j < i; j++)
        acc += static_cast<index_t>(std::get<0>(right_loc_dim[j])) *
               std::get<1>(right_loc_dim[j]);
      right_loc_off.push_back(acc);
    }
    //
    subnewsize.push_back({total_left__dim, total_right_dim, left__loc_basis,
                          left__loc_off, left__loc_dim, right_loc_basis,
                          right_loc_off, right_loc_dim});
  }
  return subnewsize;
}

} // namespace mdot
