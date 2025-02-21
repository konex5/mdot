#pragma once

#include "mdot/include/babel_type.hpp"

namespace mdot {

inline dnum_t vector_norm(dnum_t *psi, const size_t N) {
  const size_t inc = 1;
  return dnrm2_(&N, psi, &inc);
}

inline void vector_scal(const dnum_t a, dnum_t *vec, const size_t N) {
  const size_t inc = 1;
  dscal_(&N, &a, vec, &inc);
}


bool jacobi_davidson_ev(dnum_t *A, dnum_t *psi, size_t dim, const size_t max_iter,
                const dnum_t err_tol, dnum_t &eigval, dnum_t *eigvec) {

                }
                // 1 normalize psi
                vector_scal(1./vector_norm(psi,dim),psi,dim);


                for (int iter = 0; iter < maxIterations; ++iter) {
                // 2 Compute AV = A * psi

                // build projected matrix psi Av
                }
}