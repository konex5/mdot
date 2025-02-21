#include <boost/test/unit_test.hpp>

#include <complex>
#include <math.h>

using namespace std;

#define size_t typename std::size_t
#define znum_t std::complex<double>

extern "C" {
void zgesdd_(const char *jobz, const size_t *m, const size_t *n,
             const znum_t *a, const size_t *lda, double *s, znum_t *u,
             const size_t *ldu, znum_t *vt, const size_t *ldvt, znum_t *work,
             const int *lwork, double *rwork, int *iwork, int *info);
}

#undef size_t

BOOST_AUTO_TEST_CASE(test_simple_zgesvd) {

  const size_t N = 2;

  { // real, row major
    const znum_t A[N * N] = {{1., +1.}, {-1., -1.}, {1., -1.}, {1., 0.}};
    const znum_t Adeepcopy[N * N] = {A[0], A[1], A[2], A[3]};

    znum_t Vd[N * N] = {{-3.65148372e-01, -0.73029674},
                        {5.47722558e-01, +0.18257419},
                        {-4.08248290e-01, +0.40824829},
                        {3.33066907e-16, 0.81649658}};
    double S[N] = {2.23606798, 1.41421356};
    znum_t U[N * N] = {{-0.81649658, 0.},
                       {-0.57735027, 0.},
                       {0.40824829, +0.40824829},
                       {-0.57735027, -0.57735027}};

    for (size_t i = 0; i < N; i++)
      for (size_t j = 0; j < N; j++) {
        znum_t sum = 0;
        for (size_t k = 0; k < N; k++)
          sum += U[i * N + k] * S[k] * Vd[k * N + j];
        BOOST_CHECK(abs(Adeepcopy[i * N + j] - sum) < 1e-5);
      };

    znum_t Uout[N * N], VDout[N * N];
    double Sout[N];
    size_t ldA = N, ldu = N, ldvT = N;
    znum_t worktest;
    double rwork[5 * min(N, N)];
    int info, lwork = -1;
    int iwork[8 * N];

    zgesdd_((char *)"S", &N, &N, A, &ldA, Sout, VDout, &ldvT, Uout, &ldu,
            &worktest, &lwork, rwork, iwork, &info);

    lwork = (int)worktest.real();
    znum_t work[lwork];

    zgesdd_((char *)"S", &N, &N, A, &ldA, Sout, VDout, &ldvT, Uout, &ldu, work,
            &lwork, rwork, iwork, &info);

    for (size_t k = 0; k < N; k++)
      BOOST_CHECK(abs(S[k] - Sout[k]) < 1e-5);
    for (size_t k = 0; k < N * N; k++)
      BOOST_CHECK(abs(U[k] - Uout[k]) < 1e-5);

    for (size_t k = 0; k < N * N; k++)
      BOOST_CHECK(abs(Vd[k] - VDout[k]) < 1e-5);

    for (size_t i = 0; i < N; i++)
      for (size_t j = 0; j < N; j++) {
        znum_t sum = 0;
        for (size_t k = 0; k < N; k++)
          sum += Uout[i * N + k] * Sout[k] * VDout[k * N + j];
        BOOST_CHECK(abs(Adeepcopy[i * N + j] - sum) < 1e-5);
      };
  }
  { // real, column major
    const znum_t A[N * N] = {{1., +1.}, {1., -1.}, {-1., -1.}, {1., 0.}};
    const znum_t Adeepcopy[N * N] = {A[0], A[1], A[2], A[3]};

    znum_t U[N * N] = {{-0.36514837, -7.30296743e-01},
                       {-0.18257419, +5.47722558e-01},
                       {-0.40824829, +4.08248290e-01},
                       {-0.81649658, +2.77555756e-16}};
    double S[N] = {2.23606798, 1.41421356};
    znum_t Vd[N * N] = {{-0.81649658, 0.},
                        {-0.57735027, 0.},
                        {0.40824829, -0.40824829},
                        {-0.57735027, +0.57735027}};

    for (size_t i = 0; i < N; i++)
      for (size_t j = 0; j < N; j++) {
        znum_t sum = 0;
        for (size_t k = 0; k < N; k++)
          sum += U[i + k * N] * S[k] * Vd[k + j * N];
        BOOST_CHECK(abs(A[i + j * N] - sum) < 1e-5);
      };

    znum_t Uout[N * N], VDout[N * N];
    double Sout[N];
    size_t ldA = N, ldu = N, ldvT = N;
    znum_t worktest;
    double rwork[5 * min(N, N)];
    int info, lwork = -1;
    int iwork[8 * N];

    zgesdd_((char *)"S", &N, &N, A, &ldA, Sout, Uout, &ldvT, VDout, &ldu,
            &worktest, &lwork, rwork, iwork, &info);

    lwork = (int)worktest.real();
    znum_t work[lwork];

    zgesdd_((char *)"S", &N, &N, A, &ldA, Sout, Uout, &ldvT, VDout, &ldu, work,
            &lwork, rwork, iwork, &info);

    for (size_t k = 0; k < N; k++)
      BOOST_CHECK(S[k] - Sout[k] < 1e-5);

    for (size_t k = 0; k < N * N; k++)
      BOOST_CHECK(abs(U[k] - Uout[k]) < 1e-5);

    for (size_t k = 0; k < N * N; k++)
      BOOST_CHECK(abs(Vd[k] - VDout[k]) < 1e-5);

    for (size_t i = 0; i < N; i++)
      for (size_t j = 0; j < N; j++) {
        znum_t sum = 0;
        for (size_t k = 0; k < N; k++)
          sum += Uout[i + k * N] * Sout[k] * VDout[k + j * N];
        BOOST_CHECK(abs(Adeepcopy[i + j * N] - sum) < 1e-5);
      };
  }
}

#undef znum_t
