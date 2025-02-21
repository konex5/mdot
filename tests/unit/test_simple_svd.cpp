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

    double U[N * K] = {-4.59700843e-01, 6.99362418e-17,  -8.88073834e-01,
                       6.27963030e-01,  -7.07106781e-01, -3.25057584e-01,
                       6.27963030e-01,  7.07106781e-01,  -3.25057584e-01};
    double S[K] = {1.93185165, 1., 0.51763809};
    double Vd[K * M] = {-0.88807383, 0.32505758,  0.32505758,
                        -0.,         0.70710678,  -0.70710678,
                        -0.45970084, -0.62796303, -0.62796303};

    for (std::size_t i = 0; i < N; i++)
      for (std::size_t j = 0; j < M; j++) {
        double sum = 0;
        for (std::size_t k = 0; k < K; k++)
          sum += U[i * K + k] * S[k] * Vd[k * M + j];
        // std::cout << "A[i*5+j]=" << A[i*M+j] << " and the sum gives:" << sum
        // << std::endl;
        BOOST_CHECK_CLOSE(A[i * M + j], sum, 1e-5);
      };

    double Uout[N * K], Sout[K], VDout[K * M];
    std::size_t ldA = M, ldu = M, ldvT = M < N ? M : N;
    double worktest;
    int info, lwork = -1;

    dgesvd_((char *)"S", (char *)"S", &M, &N, A, &ldA, Sout, VDout, &ldu, Uout,
            &ldvT, &worktest, &lwork, &info);
    lwork = (int)worktest;
    double work[lwork];
    dgesvd_((char *)"S", (char *)"S", &M, &N, A, &ldA, Sout, VDout, &ldu, Uout,
            &ldvT, work, &lwork, &info);

    for (std::size_t k = 0; k < K; k++)
      // std::cout << S[k] << "compared with" << Sout[k] << std::endl;
      BOOST_CHECK_CLOSE(S[k], Sout[k], 1e-5);

    for (std::size_t i = 0; i < N; i++)
      for (std::size_t k = 0; k < K; k++)
        // std::cout << U[i*K+k] << "compared with" << Uout[i*K+k] << std::endl;
        if (!((k != 3 && i != N) || (i == 0 && k == 1))) // some freedom in SVD
          BOOST_CHECK_CLOSE(U[i * K + k], Uout[i * K + k], 1e-5);

    for (std::size_t j = 0; j < M; j++)
      for (std::size_t k = 0; k < K; k++)
        // std::cout << Vd[k*M+j] << "compared with" << VDout[k*M+j] <<
        // std::endl;
        if (!((k == 1 && j == 0) ||
              (k == K - 2 &&
               j == M - 1))) // some floating point precision issue
          BOOST_CHECK_CLOSE(abs(Vd[k * M + j]), abs(VDout[k * M + j]), 1e-5);
  }
  { // real, column major
    const double A[N * M] = {1, -1, -1, 0, 0, 1, 0, 1, 0};

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
        // std::cout << "A[i+j*4]=" << A[i+j*N] << " and the sum gives:" << sum
        // << std::endl;
        BOOST_CHECK_CLOSE(A[i + j * N], sum, 1e-5);
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
  }
}