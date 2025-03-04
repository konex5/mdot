#pragma once

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <btas/QSPARSE/QSDArray.h>

extern "C" {
void daxpy_(const size_t *n, const dnum_t *alpha, const dnum_t *x,
            const size_t *incx, dnum_t *y, const size_t *incy);
dnum_t ddot_(const size_t *n, const dnum_t *x, const size_t *incx,
             const dnum_t *y, const size_t *incy);
void dscal_(const size_t *n, const dnum_t *a, dnum_t *x, const size_t *incx);
void dsymv_(const size_t *n, const dnum_t *alpha, const dnum_t *a, const size_t *lda,
            const dnum_t *x, const size_t *incx, const dnum_t *beta,
            const dnum_t *y, const size_t *incy);
}

namespace mdot {

template <size_t N>
using Functor =
    boost::function<void(const btas::QSDArray<N> &, btas::QSDArray<N> &)>;

//
// Davidson's precondition
//

template <size_t N>
void precondition(const double &eval, const btas::QSDArray<N> &diag,
                  btas::QSDArray<N> &errv) {
  for (auto ir = errv.begin(); ir != errv.end(); ++ir) {
    auto id = diag.find(ir->first);
    if (id != diag.end()) {
      auto irx = ir->second->begin();
      auto idx = id->second->begin();
      for (; irx != ir->second->end(); ++irx, ++idx) {
        double denm = eval - *idx;
        if (fabs(denm) < 1.0e-12)
          denm = 1.0e-12;
        *irx /= denm;
      }
    } else {
      btas::Dscal(1.0 / eval, (*ir->second));
    }
  }
}

//
// Davidson eigen solver
//

template <size_t N>
double diagonalize(const Functor<N> &f_contract, const btas::QSDArray<N> &diag,
                   btas::QSDArray<N> &wfnc) {
  int max_ritz = 20;

  double eval = 0.0;

  // reserve working space
  std::vector<btas::QSDArray<N>> trial(max_ritz);
  std::vector<btas::QSDArray<N>> sigma(max_ritz);

  btas::QSDcopy(wfnc, trial[0]);
  btas::Normalize(trial[0]);
  f_contract(trial[0], sigma[0]);

  int niter = 0;
  int iconv = 0;
  while (iconv < 1 && niter < 20) {
    for (int m = 1; m <= max_ritz; ++m) {
      // compute small Hamiltonian matrix
      btas::DArray<2> heff(m, m);
      btas::DArray<2> ovlp(m, m);
      for (int i = 0; i < m; ++i) {
        heff(i, i) = btas::QSDdotc(trial[i], sigma[i]);
        ovlp(i, i) = btas::QSDdotc(trial[i], trial[i]);
        for (int j = 0; j < i; ++j) {
          double hij = btas::QSDdotc(trial[i], sigma[j]);
          heff(i, j) = hij;
          heff(j, i) = hij;
          double sij = btas::QSDdotc(trial[i], trial[j]);
          ovlp(i, j) = sij;
          ovlp(j, i) = sij;
        }
      }
      // solve eigenvalue problem to obtain Ritz value & vector
      btas::DArray<2> rvec;
      btas::DArray<1> rval;
      Dsyev('V', 'U', heff, rval, rvec);
      eval = rval(0);
      // rotate trial & sigma vectors by Ritz vector
      std::vector<btas::QSDArray<N>> trial_save(m);
      std::vector<btas::QSDArray<N>> sigma_save(m);
      for (int i = 0; i < m; ++i) {
        btas::QSDcopy(trial[i], trial_save[i]);
        btas::QSDcopy(sigma[i], sigma_save[i]);
        btas::QSDscal(rvec(i, i), trial[i]);
        btas::QSDscal(rvec(i, i), sigma[i]);
      }
      for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
          if (i != j) {
            btas::QSDaxpy(rvec(i, j), trial_save[i], trial[j]);
            btas::QSDaxpy(rvec(i, j), sigma_save[i], sigma[j]);
          }
        }
      }
      // compute error vector
      btas::QSDArray<N> evec;
      btas::QSDArray<N> errv;
      btas::QSDcopy(trial[0], evec);
      btas::QSDcopy(sigma[0], errv);
      btas::QSDaxpy(-eval, evec, errv);
      double rnorm = btas::QSDdotc(errv, errv);
      if (rnorm < 1.0e-8) {
        ++iconv;
        break;
      }
      // solve correction equation
      if (m < max_ritz) {
        precondition(eval, diag, errv);
        for (int i = 0; i < m; ++i) {
          btas::Normalize(errv);
          btas::Orthogonalize(trial[i], errv);
        }
        btas::Normalize(errv);
        btas::QSDcopy(errv, trial[m]);
        sigma[m].clear();
        f_contract(trial[m], sigma[m]);
      }
    }
    ++niter;
  }
  btas::QSDcopy(trial[0], wfnc);

  return eval;
}

}; // namespace davidson
