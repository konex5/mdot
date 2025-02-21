#pragma once

#include "mdot/include/babel_type.hpp"

#include <iostream>

extern "C" {
void dgemm_(const char *transa, const char *transb, const size_t *m,
            const size_t *n, const size_t *k, const dnum_t *alpha,
            const dnum_t *a, const size_t *lda, const dnum_t *b,
            const size_t *ldb, const dnum_t *beta, dnum_t *c,
            const size_t *ldc);
}

namespace mdot {

void mul_mm_blocs_new(
    dtbloc_t &new_blocs, dmbloc_t lhs_blocs, dmbloc_t rhs_blocs,
    const std::vector<std::tuple<t_index_t, m_index_t, m_index_t>> buildtarget) {
  const double alpha = 1., beta = 0.;
  for (auto &[target, lhs_it, rhs_it] : buildtarget) {
    auto N = static_cast<std::size_t>(std::get<0>(lhs_blocs[lhs_it].first) *
                                      std::get<1>(lhs_blocs[lhs_it].first));
    auto M = static_cast<std::size_t>(std::get<1>(rhs_blocs[rhs_it].first) *
                                      std::get<2>(rhs_blocs[rhs_it].first));
    auto K = static_cast<std::size_t>(std::get<2>(lhs_blocs[lhs_it].first));
    std::vector<dnum_t> mat_out(N * M);

    dgemm_((char *)"N", (char *)"N", &M, &N, &K, &alpha,
           rhs_blocs[rhs_it].second.data(), &M, lhs_blocs[lhs_it].second.data(),
           &K, &beta, mat_out.data(), &M);

    new_blocs[target] = {{std::get<0>(lhs_blocs[lhs_it].first),
                          std::get<1>(lhs_blocs[lhs_it].first),
                          std::get<1>(rhs_blocs[rhs_it].first),
                          std::get<2>(rhs_blocs[rhs_it].first)},
                         mat_out};
  }
}

void mul_mm_blocs_dup(
    dtbloc_t &new_blocs, dmbloc_t lhs_blocs, dmbloc_t rhs_blocs,
    const std::vector<std::tuple<t_index_t, m_index_t, m_index_t>> buildtarget) {
  const double alpha = 1., beta = 0.;
  for (auto &[target, lhs_it, rhs_it] : buildtarget) {
    auto N = static_cast<std::size_t>(std::get<0>(lhs_blocs[lhs_it].first) *
                                      std::get<1>(lhs_blocs[lhs_it].first));
    auto M = static_cast<std::size_t>(std::get<1>(rhs_blocs[rhs_it].first) *
                                      std::get<2>(rhs_blocs[rhs_it].first));
    auto K = static_cast<std::size_t>(std::get<2>(lhs_blocs[lhs_it].first));
    auto NM = N * M;
    std::vector<dnum_t> mat_out(NM);

    dgemm_((char *)"N", (char *)"N", &M, &N, &K, &alpha,
           rhs_blocs[rhs_it].second.data(), &M, lhs_blocs[lhs_it].second.data(),
           &K, &beta, mat_out.data(), &M);

#pragma omp
    for (size_t i = 0; i < NM; i++) {
      new_blocs[target].second[i] += mat_out[i];
    }
    std::cout << new_blocs[target].second[0] << " ";
  }
}
} // namespace mdot