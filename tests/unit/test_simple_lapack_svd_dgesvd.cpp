#include <boost/test/unit_test.hpp>

#include <math.h>

using namespace std;

#define size_t typename std::size_t
#define dnum_t double

extern "C" {
void dgesvd_(const char *jobu, const char *jobvt, const size_t *m,
             const size_t *n, const dnum_t *a, const size_t *lda, dnum_t *s,
             dnum_t *u, const size_t *ldu, dnum_t *vt, const size_t *ldvt,
             dnum_t *work, int *lwork, int *info);
}

#undef size_t
#undef dnum_t

BOOST_AUTO_TEST_CASE(test_svd) {
  BOOST_CHECK(round(sqrt(5)) - 2.236067977);
  BOOST_CHECK_CLOSE(sqrt(5), 2.236067977, 1e-5);

  const std::size_t N = 3;
  const std::size_t K = 3;
  const std::size_t M = 3;

  { // real, row major

    const double A[N * M] = {1, 0, 0, -1, 0, 1, -1, 1, 0};
    const double Adeepcopy[N * M] = {A[0], A[1], A[2], A[3], A[4],
                                     A[5], A[6], A[7], A[8]};

    double U[N * K] = {4.59700843e-01,  0,
                       -8.88073834e-01, -6.27963030e-01,
                       7.07106781e-01,  -3.25057584e-01,
                       -6.27963030e-01, -7.07106781e-01,
                       -3.25057584e-01};
    double S[K] = {1.93185165, 1., 0.51763809};
    double Vd[K * M] = {0.88807383,  -0.32505758, -0.32505758,
                        0.,          -0.70710678, 0.70710678,
                        -0.45970084, -0.62796303, -0.62796303};

    for (std::size_t i = 0; i < N; i++)
      for (std::size_t j = 0; j < N; j++) {
        double sum = 0;
        for (std::size_t k = 0; k < N; k++)
          sum += U[i * N + k] * S[k] * Vd[k * N + j];
        BOOST_CHECK(abs(Adeepcopy[i * N + j] - sum) < 1e-5);
      };

    double Uout[N * K], Sout[K], VDout[K * M];
    std::size_t ldA = M, ldu = N, ldvT = M < N ? M : N;
    double worktest;
    int info, lwork = -1;

    dgesvd_((char *)"S", (char *)"S", &M, &N, A, &ldA, Sout, VDout, &ldu, Uout,
            &ldvT, &worktest, &lwork, &info);
    lwork = (int)worktest;
    double work[lwork];
    dgesvd_((char *)"S", (char *)"S", &M, &N, A, &ldA, Sout, VDout, &ldu, Uout,
            &ldvT, work, &lwork, &info);

    for (std::size_t k = 0; k < N; k++)
      BOOST_CHECK(abs(S[k] - Sout[k]) < 1e-5);
    for (std::size_t k = 0; k < N * N; k++)
      BOOST_CHECK(abs(U[k] - Uout[k]) < 1e-5);
    for (std::size_t k = 0; k < N * N; k++)
      BOOST_CHECK(abs(Vd[k] - VDout[k]) < 1e-5);

    for (std::size_t i = 0; i < N; i++)
      for (std::size_t j = 0; j < N; j++) {
        double sum = 0;
        for (std::size_t k = 0; k < N; k++)
          sum += Uout[i * N + k] * Sout[k] * VDout[k * N + j];
        BOOST_CHECK(abs(Adeepcopy[i * N + j] - sum) < 1e-5);
      };
  }
  { // real, column major
    const double A[N * M] = {1, -1, -1, 0, 0, 1, 0, 1, 0};
    const double Adeepcopy[N * M] = {A[0], A[1], A[2], A[3], A[4],
                                     A[5], A[6], A[7], A[8]};

    double U[N * K] = {-4.59700843e-01, 6.27963030e-01,  6.27963030e-01,
                       6.99362418e-17,  -7.07106781e-01, 7.07106781e-01,
                       -8.88073834e-01, -3.25057584e-01, -3.25057584e-01};
    double S[K] = {1.93185165, 1., 0.51763809};
    double Vd[K * M] = {-0.88807383, -0.,         -0.45970084,
                        0.32505758,  0.70710678,  -0.62796303,
                        0.32505758,  -0.70710678, -0.62796303};

    for (std::size_t i = 0; i < N; i++)
      for (std::size_t j = 0; j < M; j++) {
        double sum = 0;
        for (std::size_t k = 0; k < K; k++)
          sum += U[i + k * N] * S[k] * Vd[k + j * K];
        BOOST_CHECK(abs(A[i + j * N] - sum) < 1e-5);
      };

    double Uout[N * K], Sout[K], VDout[K * M];
    std::size_t ldA = N, ldu = M, ldvT = N < M ? N : M;
    double worktest;
    int info, lwork = -1;

    dgesvd_((char *)"S", (char *)"S", &N, &M, A, &ldA, Sout, Uout, &ldu, VDout,
            &ldvT, &worktest, &lwork, &info);
    lwork = (int)worktest;
    double work[lwork];
    dgesvd_((char *)"S", (char *)"S", &N, &M, A, &ldA, Sout, Uout, &ldu, VDout,
            &ldvT, work, &lwork, &info);

    for (std::size_t k = 0; k < K; k++)
      BOOST_CHECK_CLOSE(S[k], Sout[k], 1e-5);
    for (std::size_t k = 0; k < N * N; k++)
      BOOST_CHECK(abs(U[k] - Uout[k]) < 1e-5);
    for (std::size_t k = 0; k < N * N; k++)
      BOOST_CHECK(abs(Vd[k] - VDout[k]) < 1e-5);
    for (std::size_t i = 0; i < N; i++)
      for (std::size_t j = 0; j < M; j++) {
        double sum = 0;
        for (std::size_t k = 0; k < K; k++)
          sum += Uout[i + k * N] * Sout[k] * VDout[k + j * K];
        BOOST_CHECK(abs(Adeepcopy[i + j * N] - sum) < 1e-5);
      };
  }
}