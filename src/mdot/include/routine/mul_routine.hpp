#pragma once

#include "mdot/include/babel_type.hpp"

extern "C" {
void dgemm_(const char *transa, const char *transb, const size_t *m,
            const size_t *n, const size_t *k, const dnum_t *alpha,
            const dnum_t *a, const size_t *lda, const dnum_t *b,
            const size_t *ldb, const dnum_t *beta, dnum_t *c,
            const size_t *ldc);
}

namespace mdot {

void mul_mm_blocs_new(
    dtbloc_t new_blocs, dmbloc_t lhs_blocs, dmbloc_t rhs_blocs,
    std::vector<std::tuple<t_index_t, m_index_t, m_index_t>> buildtarget) {
  const double alpha = 1., beta = 0.;
  for (auto &[target, lhs_it, rhs_it] : buildtarget) {
    auto N = static_cast<std::size_t>(std::get<0>(lhs_blocs[lhs_it].first));
    auto M = static_cast<std::size_t>(std::get<2>(rhs_blocs[rhs_it].first));
    auto K = static_cast<std::size_t>(std::get<1>(lhs_blocs[lhs_it].first) *
                                      std::get<1>(lhs_blocs[lhs_it].first));
    std::vector<dnum_t> mat_out(N * M);

    dgemm_((char *)"N", (char *)"N", &M, &N, &K, &alpha,
           rhs_blocs[rhs_it].second.data(), &M, lhs_blocs[lhs_it].second.data(),
           &K, &beta, mat_out.data(), &M);
    /*
        new_blocks[target] = _np.tensordot(
            old_blocks1[it1],
            old_blocks2[it2],
            axes=(2, 0),
        )
    */
    new_blocs[target] = {{std::get<0>(lhs_blocs[lhs_it].first),
                          std::get<1>(lhs_blocs[lhs_it].first),
                          std::get<1>(rhs_blocs[rhs_it].first),
                          std::get<2>(rhs_blocs[rhs_it].first)},
                         mat_out};
  }
}

void mul_mm_blocs_dup(
    dtbloc_t new_blocks, dmbloc_t old_blocks1, dmbloc_t old_blocks2,
    std::vector<std::tuple<t_index_t, m_index_t, m_index_t>> buildtarget) {
  /*
  for target, it1, it2 in buildtarget:
      new_blocks[target] += _np.tensordot(
          old_blocks1[it1],
          old_blocks2[it2],
          axes=(2, 0),
      )
  */
}

} // namespace mdot