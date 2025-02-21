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
    const std::vector<std::tuple<t_index_t, m_index_t, m_index_t>>
        buildtarget) {
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
    const std::vector<std::tuple<t_index_t, m_index_t, m_index_t>>
        buildtarget) {
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

#pragma omp parallel
    for (size_t i = 0; i < NM; i++) {
      new_blocs[target].second[i] += mat_out[i];
    }
    std::cout << new_blocs[target].second[0] << " ";
  }
}

void mul_usv_nondeg(std::vector<std::vector<dnum_t>> &array_U,
                    std::vector<std::vector<dnum_t>> &array_S,
                    std::vector<index_t> &cut,
                    std::vector<std::vector<dnum_t>> &array_V,
                    std::vector<std::pair<index_t, t_index_t>> &nondeg,
                    std::vector<t_shape_t> &nondeg_dims,
                    dmbloc_t &dst_lhs_blocs, dmbloc_t &dst_rhs_blocs,
                    const int is_um) {
  for (std::size_t i = 0; i < nondeg.size(); i++) {
    std::cout << "hello world!" << std::endl;
    std::cout << cut[i] << " and " << cut.size() << std::endl;

    if (cut[i] > 0) {
      auto middle_index = nondeg[i].first;
      auto theta_index = nondeg[i].second;

      /*
        for _ in range(len(nondeg)):  # reversed, and passed by pop.
            Dsi = cut.pop()
            if Dsi > 0:
                dims = nondeg_dims.pop()
                tmp_nondeg = nondeg.pop()
                if is_um is None:
                    diag_sqrt = _np.diag(_np.sqrt(array_S.pop()[:Dsi]))
                    # M
                    mat_left = _np.dot(array_U.pop()[:, :Dsi],
        diag_sqrt).reshape( dims[0], dims[1], Dsi
                    )
                    # M
                    mat_right = _np.dot(diag_sqrt, array_V.pop()[:Dsi,
        :]).reshape( Dsi, dims[2], dims[3]
                    )
                elif is_um:
                    # U
                    mat_left = array_U.pop()[:, :Dsi].reshape(dims[0], dims[1],
        Dsi) # M mat_right = _np.dot( _np.diag(array_S.pop()[:Dsi]),
        array_V.pop()[:Dsi, :]
                    ).reshape(Dsi, dims[2], dims[3])
                else:
                    # M
                    mat_left = _np.dot(
                        array_U.pop()[:, :Dsi], _np.diag(array_S.pop()[:Dsi])
                    ).reshape(dims[0], dims[1], Dsi)
                    # V
                    mat_right = array_V.pop()[:Dsi, :].reshape(Dsi, dims[2],
        dims[3])

                dst_lhs_blocs[
                    (tmp_nondeg[1][0], tmp_nondeg[1][1], tmp_nondeg[0])
                ] = mat_left
                dst_rhs_blocs[
                    (tmp_nondeg[0], tmp_nondeg[1][2], tmp_nondeg[1][3])
                ] = mat_right
            else:
                array_U.pop()
                array_V.pop()
                array_S.pop()
                nondeg_dims.pop()
                nondeg.pop()
                */
    }
  }
}

} // namespace mdot