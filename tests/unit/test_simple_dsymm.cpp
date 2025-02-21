#include <boost/test/unit_test.hpp>

#include <math.h>

using namespace std;

#define size_t typename std::size_t
#define dnum_t double

extern "C" {
void dsymm_(const char *side, const char *uplo, const size_t *m,
            const size_t *n, const dnum_t *alpha, const dnum_t *a,
            const size_t *lda, const dnum_t *b, const size_t *ldb,
            const dnum_t *beta, dnum_t *c, const size_t *ldc);
}

#undef size_t
#undef dnum_t

BOOST_AUTO_TEST_CASE(test_dgemm_simple) {

  { // real, row major
    const std::size_t N = 2;
    const std::size_t M = 2;
    const std::size_t K = 2;

    const double A[N * N] = {1, 1, 1, 1};
    const double B[N * N] = {2, -1, -1, -1};
    const double C[N * N] = {1, -2, 1, -2};

    for (std::size_t i = 0; i < N; i++)
      for (std::size_t j = 0; j < M; j++) {
        double sum = 0;
        for (std::size_t k = 0; k < K; k++)
          sum += A[i * K + k] * B[k * M + j];
        BOOST_CHECK_EQUAL(C[i * M + j], sum);
      };

    double Cout[N * (N + 1) / 2];
    double alpha = 1., beta = 0.;

    dsymm_((char *)"L", (char *)"U", &M, &N, &alpha, B, &M, A, &K, &beta, Cout,
           &M);

    for (std::size_t k = 0; k < N * M; k++)
      BOOST_CHECK_EQUAL(C[k], Cout[k]);
  }

  { // real, column major
    const std::size_t M = 2;
    const std::size_t K = 2;
    const std::size_t N = 2;

    const double A[N * N] = {1, 1, 1, 1};
    const double B[N * N] = {2, -1, -1, -1};
    const double C[N * N] = {1, 1, -2, -2};

    for (std::size_t i = 0; i < M; i++)
      for (std::size_t j = 0; j < N; j++) {
        double sum = 0;
        for (std::size_t k = 0; k < K; k++)
          sum += A[i + k * M] * B[k + j * K];
        BOOST_CHECK_EQUAL(C[i + j * M], sum);
      };

    double Cout[N * M];
    double alpha = 1., beta = 0.;

    dsymm_((char *)"U", (char *)"U", &N, &M, &alpha, A, &N, B, &K, &beta, Cout,
           &N);

    for (std::size_t k = 0; k < N * M; k++)
      // std::cout << C[k] << "compared with" << Cout[k] << std::endl;
      BOOST_CHECK_EQUAL(C[k], Cout[k]);
  }
}