#pragma once

#include "mdot/include/babel_type.hpp"
#include "mdot/include/routine/indices.hpp"

namespace mdot {

void mm_to_theta_no_gate(
    dtbloc_t dst_blocs,
    dmbloc_t lhs_blocs,
    dmbloc_t rhs_blocs,
    bool conserve_left_right = false
) {
    std::vector<m_index_t> left_indices;
    for (auto& [key,_]: lhs_blocs)
        left_indices.push_back(key);
    std::vector<m_index_t> right_indices;
    for (auto& [key,_]: rhs_blocs)
        right_indices.push_back(key);
    
    auto about_indices_to_contract = indices_dst_theta_no_gate(left_indices,right_indices,conserve_left_right);
    /* 
    mul_mm_blocs(dst_blocs, lhs_blocs, rhs_blocs, dest_indices)
    */
}
    
}