#include <boost/test/unit_test.hpp>

#include <complex>
#include <math.h>

using namespace std;

#define size_t typename std::size_t
#define dnum_t double
#define znum_t std::complex<double>

extern "C" {

void zhemm_(const char *side, const char *uplo, const size_t *m,
            const size_t *n, const znum_t *alpha, const znum_t *a,
            const size_t *lda, const znum_t *b, const size_t *ldb,
            const znum_t *beta, znum_t *c, const size_t *ldc);
}

#undef size_t

BOOST_AUTO_TEST_CASE(test_zhemm_simple) {

  { // real, row major
    const std::size_t N = 2;
    const std::size_t K = 2;
    const std::size_t M = 2;
    const znum_t A[N * K] = {{1, 1}, {-1, 0}, {-1, 0}, {0, 2}};
    const znum_t B[K * M] = {{-1, 0}, {1, 0}, {1, 0}, {1, 1}};
    const znum_t C[N * M] = {{-2, -1}, {0, 0}, {1, 2}, {-3, 2}};

    for (std::size_t i = 0; i < N; i++)
      for (std::size_t j = 0; j < M; j++) {
        znum_t sum = 0;
        for (std::size_t k = 0; k < K; k++)
          sum += A[i * K + k] * B[k * M + j];
        // std::cout << "A[i*5+j]=" << A[i*M+j] << " and the sum gives:" << sum
        // << std::endl;
        BOOST_CHECK_EQUAL(C[i * M + j], sum);
      };

    znum_t Cout[N * M];
    znum_t alpha = {1., 0}, beta = {0., 0};

    zhemm_((char *)"L", (char *)"U", &M, &N, &alpha, B, &M, A, &K, &beta, Cout,
           &M);

    for (std::size_t k = 0; k < N * M; k++)
      // std::cout << C[k] << "compared with" << Cout[k] << std::endl;
      BOOST_CHECK_EQUAL(C[k], Cout[k]);

    std::cout << std::endl << std::endl;
  }
  { // real, column major
    const std::size_t M = 2;
    const std::size_t K = 2;
    const std::size_t N = 3;
    const znum_t A[N * K] = {{1, 1}, {-1, 0}, {-1, 0}, {0, 2}};
    const znum_t B[K * M] = {{-1, 0}, {1, 0}, {1, 0}, {1, 1}};
    const znum_t C[N * M] = {{-2, -1}, {1, 2}, {0, 0}, {-3, 2}};

    for (std::size_t i = 0; i < M; i++)
      for (std::size_t j = 0; j < N; j++) {
        znum_t sum = 0;
        for (std::size_t k = 0; k < K; k++)
          sum += A[i + k * M] * B[k + j * K];
        // std::cout << "A[i+j*4]=" << C[i+j*N] << " and the sum gives:" << sum
        // << std::endl;
        BOOST_CHECK_EQUAL(C[i + j * M], sum);
      };

    znum_t Cout[N * M];
    znum_t alpha = 1., beta = 0.;

    // zgemm_((char *)"T", (char *)"T", &N, &M, &K, &alpha, B, &K,
    //          A, &M, &beta, Cout, &N); // gives C^T
    // zgemm_((char *)"N", (char *)"N", &M, &N, &K, &alpha, A, &M, B, &K, &beta,
    //        Cout, &M);

    // for (std::size_t k = 0; k < N * M; k++)
    //   // std::cout << C[k] << "compared with" << Cout[k] << std::endl;
    //   BOOST_CHECK(abs(C[k] - Cout[k]) < 1e-7);
  }
}
#undef znum_t
