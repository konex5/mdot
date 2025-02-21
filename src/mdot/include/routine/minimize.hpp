#pragma once

#include "mdot/include/babel_type.hpp"
#include "mdot/include/routine/eig_lanczos.hpp"

namespace mdot {

void lanczos_on_m(dmbloc_t &dst, dmenvbloc_t env_bloc, dmbloc_t &psi,
                  const size_t max_iteration, const double tolerance) {
  for (auto &[env_key, env_value] : env_bloc) {
    if (std::get<0>(env_key) == std::get<3>(env_key) &&
        std::get<1>(env_key) == std::get<4>(env_key) &&
        std::get<2>(env_key) == std::get<5>(env_key)) {
      m_index_t key = {std::get<0>(env_key), std::get<1>(env_key),
                       std::get<2>(env_key)};
      m_shape_t shape = {std::get<0>(env_value.first),
                         std::get<1>(env_value.first),
                         std::get<2>(env_value.first)};
      size_t dim = std::get<0>(shape) * std::get<1>(shape) * std::get<2>(shape);
      //
      std::vector<dnum_t> mat_psi;
      if (psi.find(key) == psi.end()) {
        auto tmp = std::vector<dnum_t>(dim, 1 / double(dim));
        mat_psi.swap(tmp);
      } else {
        mat_psi.swap(psi[key].second);
      }
      dnum_t *a_ptr = env_value.second.data();
      dnum_t eigenvalue;
      std::vector<dnum_t> eigenvector(dim);
      bool result;
      //
      result = lanczos_ev(a_ptr, mat_psi.data(), dim, max_iteration, tolerance,
                          eigenvalue, eigenvector.data());
      //
      if (result)
        printf("succeeded %f ", eigenvalue);
      else
        printf("failed %f ", eigenvalue);

      dst[key] = {shape, eigenvector};
    }
  }
}

void lanczos_on_mm(dtbloc_t &dst, dmmenvbloc_t env_bloc, dtbloc_t &psi,
                   const size_t max_iteration, const double tolerance) {
  for (auto &[env_key, env_value] : env_bloc) {
    if (std::get<0>(env_key) == std::get<4>(env_key) &&
        std::get<1>(env_key) == std::get<5>(env_key) &&
        std::get<2>(env_key) == std::get<6>(env_key) &&
        std::get<3>(env_key) == std::get<7>(env_key)) {
      t_index_t key = {std::get<0>(env_key), std::get<1>(env_key),
                       std::get<2>(env_key), std::get<3>(env_key)};
      t_shape_t shape = {
          std::get<0>(env_value.first), std::get<1>(env_value.first),
          std::get<2>(env_value.first), std::get<3>(env_value.first)};
      size_t dim = std::get<0>(shape) * std::get<1>(shape) *
                   std::get<2>(shape) * std::get<3>(shape);
      //
      std::vector<dnum_t> mat_psi;
      if (psi.find(key) == psi.end()) {
        auto tmp = std::vector<dnum_t>(dim, 1 / double(dim));
        mat_psi.swap(tmp);
      } else {
        mat_psi.swap(psi[key].second);
      }
      dnum_t *a_ptr = env_value.second.data();
      dnum_t eigenvalue;
      std::vector<dnum_t> eigenvector(dim);
      bool result;
      //
      result = lanczos_ev(a_ptr, mat_psi.data(), dim, max_iteration, tolerance,
                          eigenvalue, eigenvector.data());
      //
      if (result)
        printf("succeeded %f ", eigenvalue);
      else
        printf("failed %f ", eigenvalue);

      dst[key] = {shape, eigenvector};
    }
  }
}

} // namespace mdot