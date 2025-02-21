#pragma once

#include <algorithm>
#include <bits/stl_numeric.h>

#include <tbb/tbb.h>

#include "mdot/include/babel_type.hpp"
#include <iostream>

extern "C" {
dnum_t dnrm2_(const size_t *n, const dnum_t *x, const size_t *incX);
dnum_t ddot_(const size_t *n, const dnum_t *x, const size_t *incX,
             const dnum_t *y, const size_t *incY);
void dgesvd_(const char *jobu, const char *jobvt, const size_t *m,
             const size_t *n, const dnum_t *a, const size_t *lda, dnum_t *s,
             dnum_t *u, const size_t *ldu, dnum_t *vt, const size_t *ldvt,
             dnum_t *work, int *lwork, int *info);
void zgesvd_(const char *jobu, const char *jobvt, const size_t *m,
             const size_t *n, const znum_t *a, const size_t *lda, double *s,
             znum_t *u, const size_t *ldu, znum_t *vt, const size_t *ldvt,
             znum_t *work, const int *lwork, double *rwork, int *info);
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

    auto N = static_cast<std::size_t>(std::get<0>(shape) * std::get<1>(shape));
    auto M = static_cast<std::size_t>(std::get<2>(shape) * std::get<3>(shape));
    auto K = std::min(N, M);
    ///
    std::vector<double> Uout(N * K), Sout(K), VDout(K * M);
    std::size_t ldA = M, ldu = N, ldvT = M < N ? M : N;
    double worktest;
    int info, lwork = -1;

    dgesvd_((char *)"S", (char *)"S", &M, &N, theta_bloc[key].second.data(),
            &ldA, Sout.data(), VDout.data(), &ldu, Uout.data(), &ldvT,
            &worktest, &lwork, &info);

    lwork = (int)worktest;
    double work[lwork];
    dgesvd_((char *)"S", (char *)"S", &M, &N, theta_bloc[key].second.data(),
            &ldA, Sout.data(), VDout.data(), &ldu, Uout.data(), &ldvT, work,
            &lwork, &info);

    array_of_U.push_back(Uout);
    array_of_S.push_back(Sout);
    array_of_V.push_back(VDout);
  }
}
/*
void svd_deg(
        theta_blocs
        : _Dict[tuple, _np.ndarray], deg
        : _List[_Tuple[int, _List[_Tuple[int, int, int, int]]]], subnewsize
        : _List[_List], array_of_U
        : _List[_np.ndarray], array_of_S
        : _List[_np.array], array_of_V
        : _List[_np.ndarray], )
                                  if len(theta_blocs.keys()) == 0:
                                      datatype = None
                                  else:
                                      datatype =
                                  list(theta_blocs.values())[0].dtype for i in
                                  range(len(deg)): # construct the degenerated
                                  matrix thetaDeg = _np.zeros((subnewsize[i][0],
                                  subnewsize[i][1]), dtype=datatype) # fill it
                                      for it in deg[i][1]:
                                          posL = subnewsize[i][2].index((it[0],
                                  it[1])) offL = subnewsize[i][3][posL] dimL =
                                  subnewsize[i][4][posL] posR =
                                  subnewsize[i][5].index((it[2], it[3])) offR =
                                  subnewsize[i][6][posR] dimR =
                                  subnewsize[i][7][posR] sliceL = slice(offL,
                                  offL + dimL[0] * dimL[1]) sliceR = slice(offR,
                                  offR + dimR[0] * dimR[1]) thetaDeg[sliceL,
                                  sliceR] = theta_blocs[it].reshape( dimL[0] *
                                  dimL[1], dimR[0] * dimR[1]
                                          )
                                      try:
                                          U, S, V = _svd(
                                              thetaDeg, full_matrices=False,
                                  compute_uv=True, overwrite_a=False
                                          )
                                      except:
                                          print("!!!!!!!!matrix badly
                                  conditioned!!!!!!!!!") U, S, V = _svd(
                                              thetaDeg,
                                              full_matrices=False,
                                              compute_uv=True,
                                              overwrite_a=True,
                                              lapack_driver="gesvd",
                                          )

                                      array_of_U.append(U)
                                      array_of_S.append(S)
                                      array_of_V.append(V)
                              */
} // namespace mdot