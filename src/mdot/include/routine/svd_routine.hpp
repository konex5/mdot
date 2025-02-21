#pragma once

#include "mdot/include/babel_type.hpp"


extern "C" {
  dnum_t dnrm2_(const size_t *n,  const dnum_t *x, const size_t *incX);
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


void bloc_norm(const array_of_s_type &array_of_s, dnum_t& norm_out) {
  norm_out = 0;
  const size_t inc = 1;
  for (auto &itout : array_of_s) {
    const size_t n = itout.size();
    norm_out+=ddot_(&n,itout.data(),&inc,itout.data(),&inc);
  }
  norm_out = sqrt(norm_out);
}




void normalize_the_array(
    std::vector<darr_t> list_of_array,
    std::vector<index_t> cut) { 
        double norm;
        if (cut.size()==0) {
            bloc_norm(list_of_array,norm);
            for (auto &itout : list_of_array) {
                std::for_each(itout.begin(),itout.end(), [&norm](dnum_t& x) {x /= norm;});
            }
        } else {
            bloc_norm(list_of_array,norm);
            for (auto &itout : list_of_array) {
                std::for_each(itout.begin(),itout.end(), [&norm](dnum_t& x) {x /= norm;});
            }
        }
        /*
                if isinstance(cut, list):
                    norm = _np.sqrt(
                        _np.sum(
                            [
                                _np.linalg.norm(arr[: cut[i]]) ** 2
                                for i, arr in enumerate(list_of_array)
                            ]
                        )
                    )
                    for i in range(len(list_of_array)):
                        list_of_array[i] /= norm
                else:
                    norm = _np.sqrt(_np.sum([_np.linalg.norm(arr) ** 2 for arr
                in list_of_array])) for i in range(len(list_of_array)):
                        list_of_array[i] /= norm
              */
}

void truncation_strategy(
    std::vector<darr_t> list_of_array, std::size_t chi_max,
    double eps_truncation_error = 1e-8) { /*-> _Tuple[list, float]:
                     #
                     # epsilon = || forall bloc s_bloc ||_2^2
                     # chi_max = max chi of bloc
                     # eps_truncation_error < sum_{i>chi_max} s_all,i^2
                     #
                     norm = _np.sqrt(_np.sum([_np.linalg.norm(arr) ** 2 for
                     arr in list_of_array]))

                     A = _np.sort(_np.concatenate(list_of_array, axis=0))
                     index2cutA = _np.searchsorted(
                         _np.cumsum(A ** 2), eps_truncation_error * norm **
                     2, side="left"
                     )

                     dw = _np.sum(A[:index2cutA] ** 2)
                     maxcutvalue = A[index2cutA]

                     del A
                     cut_at_index = [
                         min(arr.size - _np.searchsorted(arr[::-1],
                     maxcutvalue, "left"), chi_max) for arr in list_of_array
                     ]

                     return cut_at_index, dw
                 */
}

template <typename T>
void svd_nondeg(std::map<t_index_t, std::vector<T>> block_dict,
                std::vector<std::tuple<int, t_index_t, t_shape_t>>
                    nondeg, //: _List[_Tuple[int, _Tuple[int, int, int, int]]],
                // std::vector<bloc_index_t> nondeg_dims,//: _List[_Tuple[int,
                // int, int, int]],
                std::vector<std::vector<T>> array_of_U,
                std::vector<darr_t> array_of_S,
                std::vector<std::vector<T>>
                    array_of_V) { /*-> None:
                                       for i in range(len(nondeg)):
                                           dims = nondeg_dims[i]
                                           try:
                                               U, S, V = _svd(
                                                   block_dict[nondeg[i][1]].reshape(dims[0]
                                     * dims[1], dims[2] * dims[3]),
                                     full_matrices=False, compute_uv=True,
                                                   overwrite_a=True,
                                               )
                                           except:
                                               print("!!!!!!!!matrix badly
                                       conditioned!!!!!!!!!") U, S, V = _svd(
                                                   block_dict[nondeg[i][1]].reshape(dims[0]
                                     * dims[1], dims[2] * dims[3]),
                                     full_matrices=False, compute_uv=True,
                                                   overwrite_a=False,
                                                   lapack_driver="gesvd",
                                               )
                                           array_of_U.append(U)
                                           array_of_S.append(S)
                                           array_of_V.append(V)
                                   */
}
/*
void svd_deg(
        theta_blocs
        : _Dict[tuple, _np.ndarray], deg
        : _List[_Tuple[int, _List[_Tuple[int, int, int, int]]]], subnewsize
        : _List[_List], array_of_U
        : _List[_np.ndarray], array_of_S
        : _List[_np.array], array_of_V
        : _List[_np.ndarray], ) /*-> None:
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
