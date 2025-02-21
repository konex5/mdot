#pragma once

#include "mdot/include/babel_type.hpp"

extern "C" {
void daxpy_(const size_t *n, const dnum_t *alpha, const dnum_t *x,
            const size_t *incx, dnum_t *y, const size_t *incy);
dnum_t ddot_(const size_t *n, const dnum_t *x, const size_t *incx,
             const dnum_t *y, const size_t *incy);
void dscal_(const size_t *n, const dnum_t *a, dnum_t *x, const size_t *incx);
dnum_t dnrm2_(const size_t *n, const dnum_t *x, const size_t *incx);
void dgemv_(const char *trans, const size_t *m, const size_t *n,
            const dnum_t *alpha, const dnum_t *a, const size_t *lda,
            const dnum_t *x, const size_t *incx, const dnum_t *beta,
            const dnum_t *y, const size_t *incy);
void dstev_(const char *jobz, const size_t *n, const dnum_t *d, const dnum_t *e,
            const dnum_t *z, const size_t *ldaz, const dnum_t *work,
            int *info);
}

namespace mdot {

inline dnum_t vector_norm(dnum_t *psi, const size_t N) {
  const size_t inc = 1;
  return dnrm2_(&N, psi, &inc);
}

inline void vector_scal(const dnum_t a, dnum_t *vec, const size_t N) {
  const size_t inc = 1;
  dscal_(&N, &a, vec, &inc);
}

bool lanczos_ev(dnum_t *A, dnum_t *psi, size_t dim, const size_t max_iter,
                const dnum_t err_tol, dnum_t &eigval, dnum_t *eigvec) {
  const size_t N = dim;
  const size_t inc = 1;

  const size_t min_iter = 2;
  const dnum_t beta_err = 1E-15;
  // if (!(max_iter > min_iter)) {
  //   std::ostringstream err;
  //   err << "Maximum iteration number should be set greater than 2.";
  //   throw std::runtime_error(err.str());
  // }
  dnum_t a = 1;
  dnum_t alpha;
  dnum_t beta = 1;
  size_t M = max_iter;
  dnum_t *Vm = (dnum_t *)malloc((M + 1) * N * sizeof(double));
  dnum_t *As = (dnum_t *)malloc(M * sizeof(double));
  dnum_t *Bs = (dnum_t *)malloc(M * sizeof(double));
  dnum_t *d = (dnum_t *)malloc(M * sizeof(double));
  dnum_t *e = (dnum_t *)malloc(M * sizeof(double));
  size_t it = 0;
  memcpy(Vm, psi, N * sizeof(double));
  vector_scal(1 / vector_norm(psi, N), Vm, N); //
  memset(&Vm[(it + 1) * N], 0, N * sizeof(double));
  memset(As, 0, M * sizeof(double));
  memset(Bs, 0, M * sizeof(double));
  dnum_t e_diff = 1;
  dnum_t e0_old = 0;
  bool converged = false;
  while ((((e_diff > err_tol) && it < max_iter) || it < min_iter) &&
         beta > beta_err) {
    dnum_t minus_beta = -beta;
    dgemv_((char *)"T", &N, &N, &a, A, &N, &Vm[it * N], &inc, &minus_beta,
           &Vm[(it + 1) * N], &inc);
    alpha = ddot_(&N, &Vm[it * N], &inc, &Vm[(it + 1) * N], &inc);
    dnum_t minus_alpha = -alpha;
    daxpy_(&N, &minus_alpha, &Vm[it * N], &inc, &Vm[(it + 1) * N], &inc);

    beta = vector_norm(&Vm[(it + 1) * N], N);
    if (it < max_iter - 1)
      memcpy(&Vm[(it + 2) * N], &Vm[it * N], N * sizeof(double));
    As[it] = alpha;
    if (beta > beta_err) {
      vector_scal(1 / beta, &Vm[(it + 1) * N], N);
      if (it < max_iter - 1)
        Bs[it] = beta;
    } else
      converged = true;
    it++;
    if (it > 1) {
      dnum_t *z = (dnum_t *)malloc(it * it * sizeof(double));
      dnum_t *work = (dnum_t *)malloc(4 * it * sizeof(double));
      int info;
      memcpy(d, As, it * sizeof(double));
      memcpy(e, Bs, it * sizeof(double));
      dstev_((char *)"N", &it, d, e, z, &it, work, &info);
      // if (info != 0) {
      //   std::ostringstream err;
      //   err << "Error in Lapack function 'dstev': Lapack INFO = " << info;
      //   throw std::runtime_error(err.str());
      // }
      dnum_t base = fabs(d[0]) > 1 ? fabs(d[0]) : 1;
      e_diff = fabs(d[0] - e0_old) / base;
      e0_old = d[0];
      if (e_diff <= err_tol)
        converged = true;
    }
  }
  if (it > 1) {
    memcpy(d, As, it * sizeof(double));
    memcpy(e, Bs, it * sizeof(double));
    dnum_t *z = (dnum_t *)malloc(it * it * sizeof(double));
    dnum_t *work = (dnum_t *)malloc(4 * it * sizeof(double));
    int info;
    dstev_((char *)"V", &it, d, e, z, &it, work, &info);
    // if (info != 0) {
    //   std::ostringstream err;
    //   err << "Error in Lapack function 'dstev': Lapack INFO = " << info;
    //   throw std::runtime_error(err.str());
    // }
    memset(eigvec, 0, N * sizeof(double));

    for (size_t k = 0; k < it; k++) {
      daxpy_(&N, &z[k], &Vm[k * N], &inc, eigvec, &inc);
    }
    // max_iter = it;
    eigval = d[0];
    free(z), free(work);
  } else {
    // max_iter = 1;
    eigval = 0;
    memcpy(eigvec, psi, N * sizeof(double));
  }
  free(Vm), free(As), free(Bs), free(d), free(e);
  return converged;
}

} // namespace mdot