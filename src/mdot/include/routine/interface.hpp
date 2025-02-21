#pragma once

#include "mdot/include/babel_type.hpp"
#include "mdot/include/routine/indices.hpp"
#include "mdot/include/routine/mul_routine.hpp"
#include "mdot/include/routine/svd_routine.hpp"

namespace mdot {

void mm_to_theta_no_gate(dtbloc_t &dst_blocs, const dmbloc_t lhs_blocs,
                         const dmbloc_t rhs_blocs,
                         bool conserve_left_right = false) {
  std::vector<m_index_t> left_indices;
  for (auto &[key, _] : lhs_blocs)
    left_indices.push_back(key);
  std::vector<m_index_t> right_indices;
  for (auto &[key, _] : rhs_blocs)
    right_indices.push_back(key);

  auto about_indices_to_contract =
      split_degenerate_indices(indices_dst_theta_no_gate(
          left_indices, right_indices, conserve_left_right));

  mul_mm_blocs_new(dst_blocs, lhs_blocs, rhs_blocs,
                   about_indices_to_contract.first);
  mul_mm_blocs_dup(dst_blocs, lhs_blocs, rhs_blocs,
                   about_indices_to_contract.second);
}

void theta_to_mm(dtbloc_t &theta_blocs, dmbloc_t &lhs_blocs,
                 dmbloc_t &rhs_blocs, dnum_t &dw, const index_t chi_max,
                 const bool normalize, const int is_um,
                 const int direction_right, const double eps) {
  std::vector<t_index_t> theta_indices;
  for (auto &[key, _] : theta_blocs)
    theta_indices.push_back(key);

  auto middle_indices =
      potential_middle_indices(theta_indices, direction_right);
  auto out_nondeg_deg =
      degeneracy_in_theta(theta_indices, middle_indices, direction_right);

  // auto subnewsize_deg = : _List[_List] = []
  //  slices_degenerate_blocs(theta_blocs, deg, subnewsize_deg)

  std::vector<std::vector<dnum_t>> array_of_U, array_of_V;
  std::vector<std::vector<dnum_t>> array_of_S;
  std::vector<t_shape_t> nondeg_dims;
  svd_nondeg(theta_blocs, out_nondeg_deg.first, nondeg_dims, array_of_U,
             array_of_S, array_of_V);
  /// svd_deg(theta_blocs, deg, subnewsize_deg, array_of_U, array_of_S,
  /// array_of_V)
  auto cut = truncation_strategy(array_of_S, chi_max, dw, eps);

  // for (auto &arr : array_of_S)
  //   for (auto &it : arr)
  //     std::cout << " " << it;

  if (normalize)
    normalize_the_array(array_of_S, cut);

  // for (auto &arr : array_of_S)
  //   for (auto &it : arr)
  //     std::cout << " " << it;

  std::vector<index_t> cut_nondeg, cut_deg;
  for (std::size_t i = 0; i < out_nondeg_deg.first.size(); i++)
    cut_nondeg.push_back(cut[i]);
  // for (std::size_t i=0;i<out_nondeg_deg.second.size();i++)
  //   cut_deg.push_back(cut[out_nondeg_deg.first.size()+i]);

  /*
    mul_usv_deg(
        array_of_U,
        array_of_S,
        cut_deg,
        array_of_V,
        deg,
        subnewsize_deg,
        lhs_blocs,
        rhs_blocs,
        is_um=is_um,
    )*/
  mul_usv_nondeg(array_of_U, array_of_S, cut_nondeg, array_of_V,
                 out_nondeg_deg.first, nondeg_dims, lhs_blocs, rhs_blocs,
                 is_um);
}

} // namespace mdot