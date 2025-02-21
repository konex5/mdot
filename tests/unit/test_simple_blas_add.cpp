#include <boost/test/unit_test.hpp>

#include <math.h>

using namespace std;

#define size_t typename size_t
#define dnum_t double

#undef size_t
#undef dnum_t

BOOST_AUTO_TEST_CASE(test_dgemm_simple) {

  { // real, row major
    const size_t N = 3;
    const size_t M = 2;
    const double A[N * M] = {1, 1, 1, -1, 0, 1};
    const double B[N * M] = {1, 0, -1, 2, 1, -1};
    const double C[N * M] = {2, 1, 0, 1, 1, 0};

    for (size_t i = 0; i < N; i++)
      for (size_t j = 0; j < M; j++) {
        double sum = 0;
        sum += A[i * M + j] + B[i * M + j];
        BOOST_CHECK_EQUAL(C[i * M + j], sum);
      };

    double Cout[N * M];

#pragma omp
    for (size_t i = 0; i < N * M; i++)
      Cout[i] = A[i] + B[i];

    for (size_t k = 0; k < N * M; k++)
      // std::cout << C[k] << "compared with" << Cout[k] << std::endl;
      BOOST_CHECK_EQUAL(C[k], Cout[k]);
  }
}