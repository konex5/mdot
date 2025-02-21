#include <boost/test/unit_test.hpp>

#include <complex>
#include <math.h>

using namespace std;

#define size_t typename std::size_t
#define znum_t std::complex<double>

extern "C" {
void zheev_(const char *jobz, const char *uplo, const size_t *n,
            const znum_t *a, const size_t *lda, double *w, znum_t *work,
            const int *lwork, double *rwork, int *info);
}

#undef size_t

BOOST_AUTO_TEST_CASE(test_simple_zgesvd) {

  const std::size_t N = 2;

  { // real, row major // eigenvalues only
    const znum_t A[N * N] = {{1., +1.}, {-1., -1.}, {-1., -1.}, {1., 0.}};

    double wout[N];
    std::size_t ldA = N;
    znum_t worktest;
    double rwork[3 * N - 2];
    int info, lwork = -1;

    zheev_((char *)"N", (char *)"U", &N, A, &ldA, wout, &worktest, &lwork,
           rwork, &info);

    lwork = (int)worktest.real();
    znum_t work[lwork];

    zheev_((char *)"N", (char *)"U", &N, A, &ldA, wout, work, &lwork, rwork,
           &info);

  /*
  }
  { // real, column major // eigenvector too
    const znum_t A[N * M] = {{1., +1.}, {-1., -1.}, {-1., -1.}, {1., 0.}};

    zheev_((char *)"V", (char *)"S", &M, &N, A, &ldA, Sout, VDout, &ldvT, Uout,
            &ldu, &worktest, &lwork, rwork, &info);


    znum_t U[N * K] = {{-0.36514837, -7.30296743e-01},
                       {-0.18257419, +5.47722558e-01},
                       {-0.40824829, +4.08248290e-01},
                       {-0.81649658, +2.77555756e-16}};
    double S[K] = {2.23606798, 1.41421356};
    znum_t Vd[K * M] = {{-0.81649658, 0.},
                        {-0.57735027, 0.},
                        {0.40824829, -0.40824829},
                        {-0.57735027, +0.57735027}};

    for (std::size_t i = 0; i < N; i++)
      for (std::size_t j = 0; j < M; j++) {
        znum_t sum = 0;
        for (std::size_t k = 0; k < K; k++)
          sum += U[i + k * N] * S[k] * Vd[k + j * K];
        // std::cout << "A[i+j*4]=" << A[i+j*N] << " and the sum gives:" << sum
        // << std::endl;
        BOOST_CHECK(abs(A[i + j * N] - sum) < 1e-5);
      };

    znum_t Uout[N * K], VDout[K * M];
    double Sout[K];
    std::size_t ldA = N, ldu = M, ldvT = N < M ? N : M;
    znum_t worktest;
    int info, lwork = -1;

    zgesvd_((char *)"S", (char *)"S", &N, &M, A, &ldA, Sout, Uout, &ldu, VDout,
            &ldvT, &worktest, &lwork, &info);
    lwork = (int)worktest.real();
    znum_t work[lwork];
    zgesvd_((char *)"S", (char *)"S", &N, &M, A, &ldA, Sout, Uout, &ldu, VDout,
            &ldvT, work, &lwork, &info);

    for (std::size_t k = 0; k < K; k++)
      BOOST_CHECK(S[k] - Sout[k] < 1e-5);
  */}
}

#undef znum_t
