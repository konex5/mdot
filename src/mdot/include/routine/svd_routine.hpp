#pragma once

#include <algorithm>
#include <bits/stl_numeric.h>

#include <tbb/tbb.h>

#include "mdot/include/babel_type.hpp"

extern "C" {
dnum_t dnrm2_(const size_t *n, const dnum_t *x, const size_t *incX);
dnum_t ddot_(const size_t *n, const dnum_t *x, const size_t *incX,
             const dnum_t *y, const size_t *incY);
#ifdef SVD_AVOID_DIVIDE_AND_CONQUER
void dgesvd_(const char *jobu, const char *jobvt, const size_t *m,
             const size_t *n, const dnum_t *a, const size_t *lda, dnum_t *s,
             dnum_t *u, const size_t *ldu, dnum_t *vt, const size_t *ldvt,
             dnum_t *work, int *lwork, int *info);
void zgesvd_(const char *jobu, const char *jobvt, const size_t *m,
             const size_t *n, const znum_t *a, const size_t *lda, double *s,
             znum_t *u, const size_t *ldu, znum_t *vt, const size_t *ldvt,
             znum_t *work, const int *lwork, double *rwork, int *info);
#else
void dgesdd_(const char *jobz, const size_t *m, const size_t *n,
             const dnum_t *a, const size_t *lda, dnum_t *s, dnum_t *u,
             const size_t *ldu, dnum_t *vt, const size_t *ldvt, dnum_t *work,
             int *lwork, int *iwork, int *info);
#endif
}

namespace mdot {

void bloc_norm(const std::vector<std::vector<dnum_t>> &array_of_s,
               const std::vector<index_t> &cut, dnum_t &norm_out) {
  norm_out = 0;
  const size_t inc = 1;

  if (cut.size() == 0) {
    for (auto &itout : array_of_s) {
      const size_t n = itout.size();
      norm_out += ddot_(&n, itout.data(), &inc, itout.data(), &inc);
    }
    norm_out = sqrt(norm_out);
  } else {
    std::for_each(tbb::counting_iterator<size_t>(0),
                  tbb::counting_iterator<size_t>(cut.size()),
                  [array_of_s, cut, &norm_out, inc](size_t i) {
                    auto itout = array_of_s[i];
                    const size_t n =
                        cut[i] < itout.size() ? cut[i] : itout.size();
                    norm_out +=
                        ddot_(&n, itout.data(), &inc, itout.data(), &inc);
                  });
    norm_out = sqrt(norm_out);
  }
}

void normalize_the_array(std::vector<std::vector<dnum_t>> &list_of_array,
                         std::vector<index_t> cut) {
  double norm;
  bloc_norm(list_of_array, cut, norm);
  std::for_each(list_of_array.begin(), list_of_array.end(),
                [&norm](std::vector<dnum_t> &itout) {
                  // std::for_each(itout.begin(),itout.end(),[&norm](dnum_t& x)
                  // { x /=norm;});
                  std::for_each(tbb::counting_iterator<size_t>(0),
                                tbb::counting_iterator<size_t>(itout.size()),
                                [&itout, norm](size_t i) { itout[i] /= norm; });
                });
}

std::vector<index_t>
truncation_strategy(const std::vector<std::vector<dnum_t>> list_of_array,
                    const index_t chi_max, dnum_t &dw,
                    const double eps_truncation_error = 1e-8) {
  // epsilon = || forall bloc s_bloc ||_2^2
  // chi_max = max chi of bloc
  // eps_truncation_error < sum_{i>chi_max} s_all,i^2
  dnum_t norm;
  bloc_norm(list_of_array, {}, norm);
  const dnum_t stop_criterion = eps_truncation_error * pow(norm, 2);
  //
  std::vector<dnum_t> tmp;
  for (auto &arr : list_of_array)
    tmp.insert(tmp.end(), arr.begin(), arr.end());
  std::sort(tmp.begin(), tmp.end());
  // square
  std::vector<dnum_t> tmp_square;
  for (dnum_t &v : tmp)
    tmp_square.push_back(v * v);

  std::vector<dnum_t> tmp_acc(tmp_square.size());
  std::partial_sum(tmp_square.begin(), tmp_square.end(), tmp_acc.begin());

  auto lower = std::lower_bound(tmp_acc.begin(), tmp_acc.end(), stop_criterion,
                                std::less_equal<dnum_t>{});
  size_t index_to_cut = std::distance(tmp_acc.begin(), lower);
  for (auto it = tmp_square.begin(); it < tmp_square.begin() + index_to_cut;
       it++)
    dw += *it; // std::accumulate(tmp_square.begin(), tmp_square.begin() +
               // index_to_cut, 0);

  dnum_t maxcutvalue = tmp[index_to_cut];

  std::vector<index_t> cut_at_index;
  for (size_t i = 0; i < list_of_array.size(); i++) {
    auto it = std::upper_bound(list_of_array[i].begin(), list_of_array[i].end(),
                               maxcutvalue, std::greater<dnum_t>{});
    auto value =
        static_cast<index_t>(std::distance(list_of_array[i].begin(), it));
    cut_at_index.push_back(std::min(value, chi_max));
  }

  return cut_at_index;
}

void svd_nondeg(dtbloc_t &theta_bloc,
                std::vector<std::pair<index_t, t_index_t>> nondeg,
                std::vector<t_shape_t> &nondeg_dims,
                std::vector<std::vector<dnum_t>> &array_of_U,
                std::vector<std::vector<dnum_t>> &array_of_S,
                std::vector<std::vector<dnum_t>> &array_of_V) {

  for (auto &it : nondeg) {
    auto key = it.second;
    auto shape = theta_bloc[key].first;
    nondeg_dims.push_back(shape);

    auto N = static_cast<size_t>(std::get<0>(shape) * std::get<1>(shape));
    auto M = static_cast<size_t>(std::get<2>(shape) * std::get<3>(shape));
    auto K = M < N ? M : N;
    ///
    std::vector<double> Uout(N * K), Sout(K), VDout(K * M);
    size_t ldA = M, ldu = K, ldvT = M;
    double worktest;
    int info, lwork = -1;

#ifdef SVD_AVOID_DIVIDE_AND_CONQUER
    dgesvd_((char *)"S", (char *)"S", &M, &N, theta_bloc[key].second.data(),
            &ldA, Sout.data(), VDout.data(), &ldvT, Uout.data(), &ldu,
            &worktest, &lwork, &info);

    lwork = (int)worktest;
    double work[lwork];
    dgesvd_((char *)"S", (char *)"S", &M, &N, theta_bloc[key].second.data(),
            &ldA, Sout.data(), VDout.data(), &ldvT, Uout.data(), &ldu, work,
            &lwork, &info);
#else
    int iwork[8 * K];
    dgesdd_((char *)"S", &M, &N, theta_bloc[key].second.data(), &ldA,
            Sout.data(), VDout.data(), &ldvT, Uout.data(), &ldu, &worktest,
            &lwork, iwork, &info);

    lwork = (int)worktest;
    double work[lwork];
    dgesdd_((char *)"S", &M, &N, theta_bloc[key].second.data(), &ldA,
            Sout.data(), VDout.data(), &ldvT, Uout.data(), &ldu, work, &lwork,
            iwork, &info);
#endif
    array_of_U.push_back(Uout);
    array_of_S.push_back(Sout);
    array_of_V.push_back(VDout);
  }
}

void fill_matrix(std::vector<dnum_t> &dst, const size_t dst_N,
                 const size_t dst_M, const size_t dst_off0,
                 const size_t dst_off1, std::vector<dnum_t> &src,
                 const size_t src_N, const size_t src_M) {
  for (size_t i = 0; i < src_N; i++) {
    for (size_t j = 0; j < src_M; j++) {
      dst[(dst_off0 + i) * dst_M + dst_off1 + j] = src[i * src_M + j];
    }
  }
}

void svd_deg(
    dtbloc_t &theta_blocs,
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
    std::vector<std::vector<dnum_t>> &array_of_U,
    std::vector<std::vector<dnum_t>> &array_of_S,
    std::vector<std::vector<dnum_t>> &array_of_V) {
  for (size_t i = 0; i < deg.size(); i++) {
    auto dim0 = static_cast<size_t>(std::get<0>(subnewsize[i]));
    auto dim1 = static_cast<size_t>(std::get<1>(subnewsize[i]));
    std::vector<dnum_t> tmp_theta(dim0 * dim1, 0);
    for (auto &theta_key : deg[i].second) {
      auto tmp_for_findL = std::get<2>(subnewsize[i]);
      std::tuple<index_t, index_small_t> tmp_indexL = {std::get<0>(theta_key),
                                                       std::get<1>(theta_key)};
      auto posL =
          std::find(tmp_for_findL.begin(), tmp_for_findL.end(), tmp_indexL) -
          tmp_for_findL.begin();
      auto offL = static_cast<size_t>(std::get<3>(subnewsize[i])[posL]);
      auto dimL = std::get<4>(subnewsize[i])[posL];
      auto muldimL = static_cast<size_t>(std::get<0>(dimL)) *
                     static_cast<size_t>(std::get<1>(dimL));
      auto tmp_for_findR = std::get<5>(subnewsize[i]);
      std::tuple<index_t, index_small_t> tmp_indexR = {std::get<2>(theta_key),
                                                       std::get<3>(theta_key)};
      auto posR =
          std::find(tmp_for_findR.begin(), tmp_for_findR.end(), tmp_indexR) -
          tmp_for_findR.begin();
      auto offR = static_cast<size_t>(std::get<6>(subnewsize[i])[posR]);
      auto dimR = std::get<7>(subnewsize[i])[posR];
      auto muldimR = static_cast<size_t>(std::get<0>(dimR)) *
                     static_cast<size_t>(std::get<1>(dimR));
      //
      fill_matrix(tmp_theta, dim0, dim1, offL, offR,
                  theta_blocs[theta_key].second, muldimL, muldimR);
    }

    size_t N = dim0;
    size_t M = dim1;
    size_t K = M < N ? M : N;
    ///
    std::vector<double> Uout(N * K), Sout(K), VDout(K * M);
    size_t ldA = M, ldu = K, ldvT = M;
    double worktest;
    int info, lwork = -1;
#ifdef SVD_AVOID_DIVIDE_AND_CONQUER
    dgesvd_((char *)"S", (char *)"S", &M, &N, tmp_theta.data(), &ldA,
            Sout.data(), VDout.data(), &ldvT, Uout.data(), &ldu, &worktest,
            &lwork, &info);

    lwork = (int)worktest;
    double work[lwork];
    dgesvd_((char *)"S", (char *)"S", &M, &N, tmp_theta.data(), &ldA,
            Sout.data(), VDout.data(), &ldvT, Uout.data(), &ldu, work, &lwork,
            &info);
#else
    int iwork[8 * K];
    dgesdd_((char *)"S", &M, &N, tmp_theta.data(), &ldA, Sout.data(),
            VDout.data(), &ldvT, Uout.data(), &ldu, &worktest, &lwork, iwork,
            &info);

    lwork = (int)worktest;
    double work[lwork];
    dgesdd_((char *)"S", &M, &N, tmp_theta.data(), &ldA, Sout.data(),
            VDout.data(), &ldvT, Uout.data(), &ldu, work, &lwork, iwork, &info);
#endif

    array_of_U.push_back(Uout);
    array_of_S.push_back(Sout);
    array_of_V.push_back(VDout);
  }
}

} // namespace mdot