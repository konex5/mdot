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

void mul_th_gate_new(
    dtbloc_t &new_blocs, dtbloc_t th_blocs, dgbloc_t gate_blocs,
    const std::vector<std::tuple<t_index_t, t_index_t, g_index_t>>
        buildtarget) {
  for (auto &[target, th_key, gate_key] : buildtarget) {
    auto dim0 = static_cast<size_t>(std::get<0>(th_blocs[th_key].first));
    auto dim1 = static_cast<size_t>(std::get<3>(th_blocs[th_key].first));
    auto N = static_cast<size_t>(std::get<1>(gate_blocs[gate_key].first));
    auto M = static_cast<size_t>(std::get<2>(gate_blocs[gate_key].first));

    auto K = static_cast<size_t>(std::get<0>(gate_blocs[gate_key].first));
    auto L = static_cast<size_t>(std::get<3>(gate_blocs[gate_key].first));

    std::vector<dnum_t> mat_out(dim0 * N * M * dim1, 0);

    auto theta = th_blocs[th_key].second;
    auto gate = gate_blocs[gate_key].second;

    for (size_t n = 0; n < N; n++)
      for (size_t m = 0; m < M; m++)
        for (size_t k = 0; k < K; k++)
          for (size_t l = 0; l < L; l++)
            for (size_t i = 0; i < dim0; i++) {
              size_t out_off = i * (N * M * dim1) + n * (M * dim1) + m * (dim1);
              size_t in_off = i * (K * L * dim1) + k * (L * dim1) + l * (dim1);
              size_t gate_off = k * (N * L * M) + n * (L * M) + m * (L) + l;
#pragma omp parallel
              for (size_t j = 0; j < dim1; j++)
                mat_out[out_off + j] += theta[in_off + j] * gate[gate_off];
            }

    new_blocs[target] = {{std::get<0>(th_blocs[th_key].first),
                          std::get<1>(gate_blocs[gate_key].first),
                          std::get<2>(gate_blocs[gate_key].first),
                          std::get<3>(th_blocs[th_key].first)},
                         mat_out};
  }
}

void mul_th_gate_dup(
    dtbloc_t &new_blocs, dtbloc_t th_blocs, dgbloc_t gate_blocs,
    const std::vector<std::tuple<t_index_t, t_index_t, g_index_t>>
        buildtarget) {
  for (auto &[target, th_key, gate_key] : buildtarget) {
    auto dim0 = static_cast<size_t>(std::get<0>(th_blocs[th_key].first));
    auto dim1 = static_cast<size_t>(std::get<3>(th_blocs[th_key].first));
    auto N = static_cast<size_t>(std::get<1>(gate_blocs[gate_key].first));
    auto M = static_cast<size_t>(std::get<2>(gate_blocs[gate_key].first));

    auto K = static_cast<size_t>(std::get<0>(gate_blocs[gate_key].first));
    auto L = static_cast<size_t>(std::get<3>(gate_blocs[gate_key].first));

    std::vector<dnum_t> mat_out(dim0 * N * M * dim1, 0);

    auto theta = th_blocs[th_key].second;
    auto gate = gate_blocs[gate_key].second;

    for (size_t n = 0; n < N; n++)
      for (size_t m = 0; m < M; m++)
        for (size_t k = 0; k < K; k++)
          for (size_t l = 0; l < L; l++)
            for (size_t i = 0; i < dim0; i++) {
              size_t out_off = i * (N * M * dim1) + n * (M * dim1) + m * (dim1);
              size_t in_off = i * (K * L * dim1) + k * (L * dim1) + l * (dim1);
              size_t gate_off = k * (N * L * M) + n * (L * M) + m * (L) + l;
#pragma omp parallel
              for (size_t j = 0; j < dim1; j++)
                mat_out[out_off + j] += theta[in_off + j] * gate[gate_off];
            }

#pragma omp parallel
    for (size_t i = 0; i < dim0 * N * M * dim1; i++) {
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

void slice_matrix(std::vector<dnum_t> &dst, const size_t dst_N,
                  const size_t dst_M, std::vector<dnum_t> &src,
                  const size_t src_off0, const size_t src_off1) {
  for (size_t i = 0; i < dst_N; i++) {
    for (size_t j = 0; j < dst_M; j++) {
      dst[i * dst_M + j] = src[(src_off0 + i) * dst_M + src_off1 + j];
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
        &subnewsize,
    dmbloc_t &dst_lhs_blocs, dmbloc_t &dst_rhs_blocs, const int is_um) {

  for (size_t i = 0; i < deg.size(); i++) {
    if (cut[i] > 0) {
      const index_t middle_index = deg[i].first;
      auto theta_indices = deg[i].second;
      auto dimtotL = static_cast<size_t>(std::get<0>(subnewsize[i]));
      auto dimtotR = static_cast<size_t>(std::get<1>(subnewsize[i]));
      //
      std::vector<dnum_t> mat_left, mat_right;
      //
      if (is_um == 0) {
        for (auto &s : array_S[i])
          s = sqrt(s);
        mat_left.resize(dimtotL * cut[i]);
        mat_right.resize(cut[i] * dimtotR);
        mul_mat_diag(mat_left, array_U[i], dimtotL, array_S[i].size(),
                     array_S[i], cut[i]);
        mul_diag_mat(mat_right, array_V[i], array_S[i].size(), dimtotR,
                     array_S[i], cut[i]);
      } else if (is_um == 1) {
        mat_left.swap(array_U[i]);
        mat_right.resize(cut[i] * dimtotR);
        mul_diag_mat(mat_right, array_V[i], array_S[i].size(), dimtotR,
                     array_S[i], cut[i]);
      } else {
        mat_left.resize(dimtotL * cut[i]);
        mul_mat_diag(mat_left, array_U[i], dimtotL, array_S[i].size(),
                     array_S[i], cut[i]);
        mat_right.swap(array_V[i]);
      }
      for (auto &theta_key : theta_indices) {
        auto tmp_for_findL = std::get<2>(subnewsize[i]);
        std::tuple<index_t, index_small_t> tmp_indexL = {
            std::get<0>(theta_key), std::get<1>(theta_key)};
        auto posL =
            std::find(tmp_for_findL.begin(), tmp_for_findL.end(), tmp_indexL) -
            tmp_for_findL.begin();
        auto offL = static_cast<size_t>(std::get<3>(subnewsize[i])[posL]);
        auto dimL = std::get<4>(subnewsize[i])[posL];
        auto muldimL = static_cast<size_t>(std::get<0>(dimL)) *
                       static_cast<size_t>(std::get<1>(dimL));
        auto tmp_for_findR = std::get<5>(subnewsize[i]);
        std::tuple<index_t, index_small_t> tmp_indexR = {
            std::get<2>(theta_key), std::get<3>(theta_key)};
        auto posR =
            std::find(tmp_for_findR.begin(), tmp_for_findR.end(), tmp_indexR) -
            tmp_for_findR.begin();
        auto offR = static_cast<size_t>(std::get<6>(subnewsize[i])[posR]);
        auto dimR = std::get<7>(subnewsize[i])[posR];
        auto muldimR = static_cast<size_t>(std::get<0>(dimR)) *
                       static_cast<size_t>(std::get<1>(dimR));

        std::vector<dnum_t> out_mat_left(muldimL * cut[i]);
        m_shape_t shape_left = {std::get<0>(dimL), std::get<1>(dimL), cut[i]};
        std::vector<dnum_t> out_mat_right(cut[i] * muldimR);
        m_shape_t shape_right = {cut[i], std::get<0>(dimR), std::get<1>(dimR)};

        slice_matrix(out_mat_left, muldimL, cut[i], mat_left, offL, 0);
        slice_matrix(out_mat_right, cut[i], muldimR, mat_right, 0, offR);

        dst_lhs_blocs[{std::get<0>(theta_key), std::get<1>(theta_key),
                       middle_index}] = {shape_left, out_mat_left};
        dst_rhs_blocs[{middle_index, std::get<2>(theta_key),
                       std::get<3>(theta_key)}] = {shape_right, out_mat_right};
      }
    }
  }
}

} // namespace mdot