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
  { // real, row major
    const std::size_t N = 4;
    const std::size_t K = 3;
    const std::size_t M = 5;
    const double A[N * M] = {1, 0, 0, 0, 2, 0, 0, 3, 0, 0,
                             0, 0, 0, 0, 0, 0, 2, 0, 0, 0};
    double U[N * K] = {0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, -1};
    double S[K] = {3, sqrt(5), 2};
    double Vd[K * M] = {0, 0,          -1, 0,  0, -sqrt(0.2), 0, 0,
                        0, -sqrt(0.8), 0,  -1, 0, 0,          0};

    for (std::size_t i = 0; i < N; i++)
      for (std::size_t j = 0; j < M; j++) {
        double sum = 0;
        for (std::size_t k = 0; k < K; k++)
          sum += U[i * K + k] * S[k] * Vd[k * M + j];
        // std::cout << "A[i*5+j]=" << A[i*M+j] << " and the sum gives:" << sum
        // << std::endl;
        BOOST_CHECK_EQUAL(A[i * M + j], sum);
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
      BOOST_CHECK_EQUAL(S[k], Sout[k]);

    for (std::size_t i = 0; i < N; i++)
      for (std::size_t k = 0; k < K; k++)
        // std::cout << U[i*K+k] << "compared with" << Uout[i*K+k] << std::endl;
        if (!((k != 3 && i != N) || (i == 0 && k == 1))) // some freedom in SVD
          BOOST_CHECK_EQUAL(U[i * K + k], Uout[i * K + k]);

    for (std::size_t j = 0; j < M; j++)
      for (std::size_t k = 0; k < K; k++)
        // std::cout << Vd[k*M+j] << "compared with" << VDout[k*M+j] <<
        // std::endl;
        if (!((k == 1 && j == 0) ||
              (k == K - 2 &&
               j == M - 1))) // some floating point precision issue
          BOOST_CHECK_EQUAL(Vd[k * M + j], VDout[k * M + j]);
  }
  { // real, column major
    const std::size_t N = 4;
    const std::size_t K = 3;
    const std::size_t M = 5;
    const double A[N * M] = {1, 0, 0, 0, 0, 0, 0, 2, 0, 3,
                             0, 0, 0, 0, 0, 0, 2, 0, 0, 0};
    double U[N * K] = {0, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0, -1};
    double S[K] = {3, sqrt(5), 2};
    double Vd[K * M] = {0, -sqrt(0.2), 0, 0, 0, -1,         -1, 0,
                        0, 0,          0, 0, 0, -sqrt(0.8), 0};

    for (std::size_t i = 0; i < N; i++)
      for (std::size_t j = 0; j < M; j++) {
        double sum = 0;
        for (std::size_t k = 0; k < K; k++)
          sum += U[i + k * N] * S[k] * Vd[k + j * K];
        // std::cout << "A[i+j*4]=" << A[i+j*N] << " and the sum gives:" << sum
        // << std::endl;
        BOOST_CHECK_EQUAL(A[i + j * N], sum);
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
      BOOST_CHECK_EQUAL(S[k], Sout[k]);
  }
}