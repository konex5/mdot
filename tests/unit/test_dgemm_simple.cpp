#include <boost/test/unit_test.hpp>

#include <math.h>

using namespace std;

#define size_t typename std::size_t
#define dnum_t double

extern "C" {
void dgemm_(const char *transa, const char *transb, const size_t *m,
            const size_t *n, const size_t *k, const dnum_t *alpha,
            const dnum_t *a, const size_t *lda, const dnum_t *b,
            const size_t *ldb, const dnum_t *beta, dnum_t *c,
            const size_t *ldc);
}

#undef size_t
#undef dnum_t

BOOST_AUTO_TEST_CASE(test_dgemm_simple) {

  { // real, row major
    const std::size_t N = 3;
    const std::size_t K = 2;
    const std::size_t M = 2;
    const double A[N * K] = {1, 1, 1, -1, 0, 1};
    const double B[K * M] = {0, 1, 1, -1};
    const double C[N * M] = {1, 0, -1, 2, 1, -1};

    for (std::size_t i = 0; i < N; i++)
      for (std::size_t j = 0; j < M; j++) {
        double sum = 0;
        for (std::size_t k = 0; k < K; k++)
          sum += A[i * K + k] * B[k * M + j];
        // std::cout << "A[i*5+j]=" << A[i*M+j] << " and the sum gives:" << sum
        // << std::endl;
        BOOST_CHECK_CLOSE(C[i * M + j], sum, 1e-5);
      };

    double Cout[N * M];
    double alpha = 1., beta = 0.;

    dgemm_((char *)"N", (char *)"N", &M, &N, &K, &alpha, B, &M, A, &K, &beta,
           Cout, &M);

    for (std::size_t k = 0; k < N * M; k++)
      // std::cout << C[k] << "compared with" << Cout[k] << std::endl;
      BOOST_CHECK_CLOSE(C[k], Cout[k], 1e-7);

    std::cout << std::endl << std::endl;
  }
  { // real, column major
    const std::size_t M = 2;
    const std::size_t K = 2;
    const std::size_t N = 3;
    const double A[M * K] = {0, 1, 1, -1};
    const double B[K * N] = {1, 1, 1, -1, 0, 1};
    const double C[M * N] = {1, 0, -1, 2, 1, -1};

    for (std::size_t i = 0; i < M; i++)
      for (std::size_t j = 0; j < N; j++) {
        double sum = 0;
        for (std::size_t k = 0; k < K; k++)
          sum += A[i + k * M] * B[k + j * K];
        // std::cout << "A[i+j*4]=" << C[i+j*N] << " and the sum gives:" << sum
        // << std::endl;
        BOOST_CHECK_CLOSE(C[i + j * M], sum, 1e-5);
      };

    double Cout[N * M];
    double alpha = 1., beta = 0.;

    // dgemm_((char *)"T", (char *)"T", &N, &M, &K, &alpha, B, &K,
    //          A, &M, &beta, Cout, &N); // gives C^T
    dgemm_((char *)"N", (char *)"N", &M, &N, &K, &alpha, A, &M, B, &K, &beta,
           Cout, &M);

    for (std::size_t k = 0; k < N * M; k++)
      // std::cout << C[k] << "compared with" << Cout[k] << std::endl;
      BOOST_CHECK_CLOSE(C[k], Cout[k], 1e-7);
  }
}