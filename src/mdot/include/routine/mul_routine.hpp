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

void mul_mat_diag(std::vector<dnum_t> &destination,
                         std::vector<dnum_t> &mat, const size_t &N,
                         const size_t &M, std::vector<dnum_t> &diag,
                         const size_t &cut) {
  // const size_t inc = 1;
  // daxpy_(&cut, &(diag[i]), &(mat[i*N]), &inc, destination.data(), &inc);
  for (size_t i = 0; i < N; i++)
#pragma omp parallel
    for (size_t j = 0; j < cut; j++)
      destination[i * cut + j] = diag[j] * mat[i * M + j];
}

void mul_diag_mat(std::vector<dnum_t> &destination,
                         std::vector<dnum_t> &mat, const size_t &N,
                         const size_t &M, std::vector<dnum_t> &diag,
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
  for (std::size_t i = 0; i < nondeg.size(); i++) {
    std::cout << "hello world!" << std::endl;
    std::cout << cut[i] << " and " << cut.size() << std::endl;

    if (cut[i] > 0) {
      const index_t middle_index = nondeg[i].first;
      auto theta_index = nondeg[i].second;
      m_shape_t shape_left = {std::get<0>(nondeg_dims[i]),
                              std::get<1>(nondeg_dims[i]), cut[i]};
      m_shape_t shape_right = {cut[i], std::get<2>(nondeg_dims[i]),
                               std::get<3>(nondeg_dims[i])};
      //
      const std::size_t dim0 =
          static_cast<std::size_t>(std::get<0>(shape_left));
      const std::size_t dim1 =
          static_cast<std::size_t>(std::get<1>(shape_left));
      const std::size_t dim2 =
          static_cast<std::size_t>(std::get<1>(shape_right));
      const std::size_t dim3 =
          static_cast<std::size_t>(std::get<2>(shape_right));
      std::vector<dnum_t> mat_left(dim0 * dim1, cut[i]);
      std::vector<dnum_t> mat_right(cut[i], dim2 * dim3);

      if (is_um == 0) {

        for (auto &s : array_S[i])
          s = sqrt(s);
            std::cout << "dim0=" << dim0 << "dim1=" << dim1 << "sizeU=" << array_U[i].size() << "cut" << cut[i] << "sizeS" << array_S[i].size() << std::endl;
          /*
          mul_mat_diag(mat_left, array_U[i], dim0 * dim1, array_S.size(),
                     array_S[i], cut[i]);
          */
        /*
        mul_diag_mat(mat_right, array_V[i], array_S[i].size(), dim2 * dim3,
                     array_S[i], cut[i]);
                     */
                    
      } else if (is_um == 1) {
        //mat_left.swap(array_U[i]);
        //mul_diag_mat(mat_right, array_V[i], array_S[i].size(), dim2 * dim3,
        //             array_S[i], cut[i]);
      } else {
        //mul_mat_diag(mat_left, array_U[i], dim0 * dim1, array_S[i].size(),
        //             array_S[i], cut[i]);
        //mat_right.swap(array_V[i]);
      }
      
      dst_lhs_blocs[{std::get<0>(theta_index), std::get<1>(theta_index),
                     middle_index}] = {shape_left, mat_left};
      dst_rhs_blocs[{middle_index, std::get<2>(theta_index),
                     std::get<3>(theta_index)}] = {shape_right, mat_right};
    }
    
  }
}

} // namespace mdot