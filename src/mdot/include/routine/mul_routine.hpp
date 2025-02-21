#pragma once

#include "mdot/include/babel_type.hpp"

#include <iostream>

extern "C" {
void dgemm_(const char *transa, const char *transb, const size_t *m,
            const size_t *n, const size_t *k, const dnum_t *alpha,
            const dnum_t *a, const size_t *lda, const dnum_t *b,
            const size_t *ldb, const dnum_t *beta, dnum_t *c,
            const size_t *ldc);
void daxpy_(const size_t *N, const dnum_t *alpha, const dnum_t *X,
            const size_t *incX, dnum_t *Y, const size_t *incY);
}

namespace mdot {

void mul_mm_blocs_new(
    dtbloc_t &new_blocs, dmbloc_t lhs_blocs, dmbloc_t rhs_blocs,
    const std::vector<std::tuple<t_index_t, m_index_t, m_index_t>>
        buildtarget) {
  const double alpha = 1., beta = 0.;
  for (auto &[target, lhs_it, rhs_it] : buildtarget) {
    auto N = static_cast<size_t>(std::get<0>(lhs_blocs[lhs_it].first) *
                                 std::get<1>(lhs_blocs[lhs_it].first));
    auto M = static_cast<size_t>(std::get<1>(rhs_blocs[rhs_it].first) *
                                 std::get<2>(rhs_blocs[rhs_it].first));
    auto K = static_cast<size_t>(std::get<2>(lhs_blocs[lhs_it].first));
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
    auto N = static_cast<size_t>(std::get<0>(lhs_blocs[lhs_it].first) *
                                 std::get<1>(lhs_blocs[lhs_it].first));
    auto M = static_cast<size_t>(std::get<1>(rhs_blocs[rhs_it].first) *
                                 std::get<2>(rhs_blocs[rhs_it].first));
    auto K = static_cast<size_t>(std::get<2>(lhs_blocs[lhs_it].first));
    auto NM = N * M;
    std::vector<dnum_t> mat_out(NM);

    dgemm_((char *)"N", (char *)"N", &M, &N, &K, &alpha,
           rhs_blocs[rhs_it].second.data(), &M, lhs_blocs[lhs_it].second.data(),
           &K, &beta, mat_out.data(), &M);

#pragma omp parallel
    for (size_t i = 0; i < NM; i++) {
      new_blocs[target].second[i] += mat_out[i];
    }
  }
}

void mul_mat_diag(std::vector<dnum_t> &destination,
                  const std::vector<dnum_t> &mat, const size_t &N,
                  const size_t &M, const std::vector<dnum_t> &diag,
                  const size_t &cut) {
  // const size_t inc = 1;
  // daxpy_(&cut, &(diag[i]), &(mat[i*N]), &inc, destination.data(), &inc);
  for (size_t i = 0; i < N; i++)
#pragma omp parallel
    for (size_t j = 0; j < cut; j++)
      destination[i * cut + j] = diag[j] * mat[i * M + j];
}

void mul_diag_mat(std::vector<dnum_t> &destination,
                  const std::vector<dnum_t> &mat, const size_t &N,
                  const size_t &M, const std::vector<dnum_t> &diag,
                  const size_t &cut) {
  // const size_t inc = 1;
  // daxpy_(&cut, &(diag[i]), &(mat[i*N]), &inc, destination.data(), &inc);
  for (size_t i = 0; i < cut; i++)
#pragma omp parallel
    for (size_t j = 0; j < M; j++)
      destination[i * M + j] = diag[i] * mat[i * M + j];
}

void mul_usv_nondeg(std::vector<std::vector<dnum_t>> &array_U,
                    std::vector<std::vector<dnum_t>> &array_S,
                    std::vector<index_t> &cut,
                    std::vector<std::vector<dnum_t>> &array_V,
                    std::vector<std::pair<index_t, t_index_t>> &nondeg,
                    std::vector<t_shape_t> &nondeg_dims,
                    dmbloc_t &dst_lhs_blocs, dmbloc_t &dst_rhs_blocs,
                    const int is_um) {
  for (size_t i = 0; i < nondeg.size(); i++) {
    if (cut[i] > 0) {
      const index_t middle_index = nondeg[i].first;
      auto theta_index = nondeg[i].second;
      m_shape_t shape_left = {std::get<0>(nondeg_dims[i]),
                              std::get<1>(nondeg_dims[i]), cut[i]};
      m_shape_t shape_right = {cut[i], std::get<2>(nondeg_dims[i]),
                               std::get<3>(nondeg_dims[i])};
      std::vector<dnum_t> mat_left, mat_right;
      //
      const size_t dim0 = static_cast<size_t>(std::get<0>(shape_left));
      const size_t dim1 = static_cast<size_t>(std::get<1>(shape_left));
      const size_t dim2 = static_cast<size_t>(std::get<1>(shape_right));
      const size_t dim3 = static_cast<size_t>(std::get<2>(shape_right));

      if (is_um == 0) {
        for (auto &s : array_S[i])
          s = sqrt(s);
        mat_left.resize(dim0 * dim1 * cut[i]);
        mat_right.resize(cut[i] * dim2 * dim3);
        mul_mat_diag(mat_left, array_U[i], dim0 * dim1, array_S[i].size(),
                     array_S[i], cut[i]);
        mul_diag_mat(mat_right, array_V[i], array_S[i].size(), dim2 * dim3,
                     array_S[i], cut[i]);
      } else if (is_um == 1) {
        mat_left.swap(array_U[i]);
        mat_right.resize(cut[i] * dim2 * dim3);
        mul_diag_mat(mat_right, array_V[i], array_S[i].size(), dim2 * dim3,
                     array_S[i], cut[i]);
      } else {
        mat_left.resize(dim0 * dim1 * cut[i]);
        mul_mat_diag(mat_left, array_U[i], dim0 * dim1, array_S[i].size(),
                     array_S[i], cut[i]);
        mat_right.swap(array_V[i]);
      }

      dst_lhs_blocs[{std::get<0>(theta_index), std::get<1>(theta_index),
                     middle_index}] = {shape_left, mat_left};
      dst_rhs_blocs[{middle_index, std::get<2>(theta_index),
                     std::get<3>(theta_index)}] = {shape_right, mat_right};
    }
  }
}

void mul_usv_deg(
    std::vector<std::vector<dnum_t>> &array_U,
    std::vector<std::vector<dnum_t>> &array_S, std::vector<index_t> &cut,
    std::vector<std::vector<dnum_t>> &array_V,
    std::vector<std::pair<index_t, std::vector<t_index_t>>> &deg,
    std::vector<
        std::tuple<index_t, index_t,
                   typename std::vector<std::tuple<index_t, index_small_t>>,
                   typename std::vector<index_t>,
                   typename std::vector<std::pair<index_t, index_small_t>>,
                   typename std::vector<std::tuple<index_small_t, index_t>>,
                   typename std::vector<index_t>,
                   typename std::vector<std::pair<index_small_t, index_t>>>>
        &new_subsize,
    dmbloc_t &dst_lhs_blocs, dmbloc_t &dst_rhs_blocs, const int is_um) {

  for (size_t i = 0; i < deg.size(); i++) {
    if (cut[i]>0) {
      std::cout << "hello";
    }
  }
/*
      if (is_um == 0) {
        for (auto &s : array_S[i])
          s = sqrt(s);
        mat_left.resize(dim0 * dim1 * cut[i]);
        mat_right.resize(cut[i] * dim2 * dim3);
        mul_mat_diag(mat_left, array_U[i], dim0 * dim1, array_S[i].size(),
                     array_S[i], cut[i]);
        mul_diag_mat(mat_right, array_V[i], array_S[i].size(), dim2 * dim3,
                     array_S[i], cut[i]);
      } else if (is_um == 1) {
        mat_left.swap(array_U[i]);
        mat_right.resize(cut[i] * dim2 * dim3);
        mul_diag_mat(mat_right, array_V[i], array_S[i].size(), dim2 * dim3,
                     array_S[i], cut[i]);
      } else {
        mat_left.resize(dim0 * dim1 * cut[i]);
        mul_mat_diag(mat_left, array_U[i], dim0 * dim1, array_S[i].size(),
                     array_S[i], cut[i]);
        mat_right.swap(array_V[i]);
      }

for _ in range(len(deg)):  # reversed, and pop each value.
        Dsi = cut.pop()
        if Dsi > 0:
            if is_um is None:
                diag_sqrt = _np.diag(_np.sqrt(array_S.pop()[:Dsi]))
                # M
                mat_left = _np.dot(array_U.pop()[:, :Dsi], diag_sqrt)
                # M
                mat_right = _np.dot(diag_sqrt, array_V.pop()[:Dsi, :])
            elif is_um:
                # U
                mat_left = array_U.pop()  # [:,:Dsi]
                # M
                mat_right = _np.dot(
                    _np.diag(array_S.pop()[:Dsi]), array_V.pop()[:Dsi, :]
                )
            else:
                # M
                mat_left = _np.dot(
                    array_U.pop()[:, :Dsi], _np.diag(array_S.pop()[:Dsi])
                )
                # V
                mat_right = array_V.pop()  # [Dsi:,:]

            tmp = subnewsize.pop()
            tmp_deg = deg.pop()
            for it in tmp_deg[1]:
                posL = tmp[2].index((it[0], it[1]))
                offL = tmp[3][posL]
                dimL = tmp[4][posL]
                posR = tmp[5].index((it[2], it[3]))
                offR = tmp[6][posR]
                dimR = tmp[7][posR]

                dst_lhs_blocs[(it[0], it[1], tmp_deg[0])] = mat_left[
                    slice(offL, offL + dimL[0] * dimL[1]), :Dsi
                ].reshape(dimL[0], dimL[1], Dsi)
                dst_rhs_blocs[(tmp_deg[0], it[2], it[3])] = mat_right[
                    :Dsi, slice(offR, offR + dimR[0] * dimR[1])
                ].reshape(Dsi, dimR[0], dimR[1])
*/

}

} // namespace mdot