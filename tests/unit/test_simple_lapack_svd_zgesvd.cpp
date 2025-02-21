#include <boost/test/unit_test.hpp>

#include <complex>
#include <math.h>

using namespace std;

#define size_t typename std::size_t
#define znum_t std::complex<double>

extern "C" {
void zgesvd_(const char *jobu, const char *jobvt, const size_t *m,
             const size_t *n, const znum_t *a, const size_t *lda, double *s,
             znum_t *u, const size_t *ldu, znum_t *vt, const size_t *ldvt,
             znum_t *work, const int *lwork, double *rwork, int *info);
}

#undef size_t

BOOST_AUTO_TEST_CASE(test_simple_zgesvd) {

  const std::size_t N = 2;

  { // real, row major
    const znum_t A[N * N] = {{1., +1.}, {-1., -1.}, {1., -1.}, {1., 0.}};
    const znum_t Adeepcopy[N * N] = {A[0],A[1],A[2],A[3]};

    znum_t U[N * N] = {{-0.36514837, -7.30296743e-01},
                       {-0.40824829, +4.08248290e-01},
                       {-0.18257419, +5.47722558e-01},
                       {-0.81649658, +2.77555756e-16}};
    double S[N] = {2.23606798, 1.41421356};
    znum_t Vd[N * N] = {{-0.81649658, 0.},
                        {0.40824829, -0.40824829},
                        {-0.57735027, 0.},
                        {-0.57735027, +0.57735027}};

    for (std::size_t i = 0; i < N; i++)
      for (std::size_t j = 0; j < N; j++) {
        znum_t sum = 0;
        for (std::size_t k = 0; k < N; k++)
          sum += U[i * N + k] * S[k] * Vd[k * N + j];
        // std::cout << "A[i*5+j]=" << A[i*N+j] << " and the sum gives:" << sum
        // << std::endl;
        BOOST_CHECK(abs(A[i * N + j] - sum) < 1e-5);
      };

    znum_t Uout[N * N], VDout[N * N];
    double Sout[N];
    std::size_t ldA = N, ldu = N, ldvT = N;
    znum_t worktest;
    double rwork[5 * min(N, N)];
    int info, lwork = -1;
    
    zgesvd_((char *)"N", (char *)"N", &N, &N, A, &ldA, Sout, VDout, &ldvT, Uout,
            &ldu, &worktest, &lwork, rwork, &info);

    lwork = (int)worktest.real();
    znum_t work[lwork];

    zgesvd_((char *)"N", (char *)"N", &N, &N, A, &ldA, Sout, VDout, &ldvT, Uout,
            &ldu, work, &lwork, rwork, &info);

    std::cout << "info=" << info << " zero is ok!" << std::endl;

    for (std::size_t k = 0; k < N; k++)
      BOOST_CHECK(abs(S[k] - Sout[k]) < 1e-5);
    /*
    for (std::size_t k = 0; k < N*N; k++)
      BOOST_CHECK(abs(U[k] - Uout[k]) < 1e-5);
    
    for (std::size_t k = 0; k < N*N; k++)
      BOOST_CHECK(abs(Vd[k] - VDout[k]) < 1e-5);

    for (std::size_t i = 0; i < N; i++)
      for (std::size_t k = 0; k < N; k++)
        if (!((k != 3 && i != N) || (i == 0 && k == 1))) // some freedom in SVD
          BOOST_CHECK(abs(U[i * N + k] - Uout[i * N + k]) < 1e-5);

    for (std::size_t j = 0; j < N; j++)
      for (std::size_t k = 0; k < N; k++) {
        BOOST_CHECK(abs(Vd[k * N + j] - VDout[k * N + j]) < 1e-5);
      }
    for (std::size_t i = 0; i < N; i++)
      for (std::size_t j = 0; j < N; j++) {
        znum_t sum = 0;
        for (std::size_t k = 0; k < N; k++)
          sum += Uout[i * N + k] * Sout[k] * Vd[k * N + j];
        BOOST_CHECK(abs(Adeepcopy[i * N + j] - sum) < 1e-5);
      };
    */
  
  }
  { // real, column major
    const znum_t A[N * N] = {{1., +1.}, {1., -1.}, {-1., -1.}, {1., 0.}};
    const znum_t Adeepcopy[N * N] = {A[0],A[1],A[2],A[3]};

    znum_t U[N * N] = {{-0.36514837, -7.30296743e-01},
                       {-0.18257419, +5.47722558e-01},
                       {-0.40824829, +4.08248290e-01},
                       {-0.81649658, +2.77555756e-16}};
    double S[N] = {2.23606798, 1.41421356};
    znum_t Vd[N * N] = {{-0.81649658, 0.},
                        {-0.57735027, 0.},
                        {0.40824829, -0.40824829},
                        {-0.57735027, +0.57735027}};

    for (std::size_t i = 0; i < N; i++)
      for (std::size_t j = 0; j < N; j++) {
        znum_t sum = 0;
        for (std::size_t k = 0; k < N; k++)
          sum += U[i + k * N] * S[k] * Vd[k + j * N];
        BOOST_CHECK(abs(A[i + j * N] - sum) < 1e-5);
      };


    znum_t Uout[N * N], VDout[N * N];
    double Sout[N];
    std::size_t ldA = N, ldu = N, ldvT = N;
    znum_t worktest;
    double rwork[5 * min(N, N)];
    int info, lwork = -1;

    zgesvd_((char *)"S", (char *)"S", &N, &N, A, &ldA, Sout, Uout, &ldvT, VDout,
            &ldu, &worktest, &lwork, rwork, &info);

    lwork = (int)worktest.real();
    znum_t work[lwork];

    zgesvd_((char *)"S", (char *)"S", &N, &N, A, &ldA, Sout, Uout, &ldvT, VDout,
            &ldu, work, &lwork, rwork, &info);

    std::cout << "info=" << info << " zero is ok!" << std::endl;

    for (std::size_t k = 0; k < N; k++)
      BOOST_CHECK(S[k] - Sout[k] < 1e-5);

    for (std::size_t k = 0; k < N*N; k++)
      BOOST_CHECK(abs(U[k] - Uout[k]) < 1e-5);
    
    for (std::size_t k = 0; k < N*N; k++)
      BOOST_CHECK(abs(Vd[k] - VDout[k]) < 1e-5);
    

  for (std::size_t i = 0; i < N; i++)
      for (std::size_t j = 0; j < N; j++) {
        znum_t sum = 0;
        for (std::size_t k = 0; k < N; k++)
          sum += Uout[i + k * N] * Sout[k] * VDout[k + j * N];
        BOOST_CHECK(abs(Adeepcopy[i + j * N] - sum) < 1e-5);
      };

  }
}


BOOST_AUTO_TEST_CASE(test_simple_zgesvd_overwrite_left) {

  const std::size_t N = 2;
  const std::size_t K = 2;
  const std::size_t M = 2;

  // { // real, row major
  //   const znum_t A[N * M] = {{1., +1.}, {-1., -1.}, {1., -1.}, {1., 0.}};

  //   znum_t U[N * K] = {{-0.36514837, -7.30296743e-01},
  //                      {-0.40824829, +4.08248290e-01},
  //                      {-0.18257419, +5.47722558e-01},
  //                      {-0.81649658, +2.77555756e-16}};
  //   double S[K] = {2.23606798, 1.41421356};
  //   znum_t Vd[K * M] = {{-0.81649658, 0.},
  //                       {0.40824829, -0.40824829},
  //                       {-0.57735027, 0.},
  //                       {-0.57735027, +0.57735027}};
  //   for (std::size_t i = 0; i < N; i++)
  //     for (std::size_t j = 0; j < M; j++) {
  //       znum_t sum = 0;
  //       for (std::size_t k = 0; k < K; k++)
  //         sum += U[i * K + k] * S[k] * Vd[k * M + j];
  //       // std::cout << "A[i*5+j]=" << A[i*M+j] << " and the sum gives:" << sum
  //       // << std::endl;
  //       BOOST_CHECK(abs(A[i * M + j] - sum) < 1e-5);
  //     };
  //   /*

  //   znum_t Uout[N * K], VDout[K * M];
  //   double Sout[K];
  //   std::size_t ldA = M, ldu = N, ldvT = M < N ? M : N;
  //   znum_t worktest;
  //   double rwork[5 * min(M, N)];
  //   int info, lwork = -1;

  //   zgesvd_((char *)"S", (char *)"S", &M, &N, A, &ldA, Sout, VDout, &ldu, Uout,
  //           &ldvT, &worktest, &lwork, rwork, &info);

  //   lwork = (int)worktest.real();
  //   znum_t work[lwork];

  //   zgesvd_((char *)"S", (char *)"S", &M, &N, A, &ldA, Sout, VDout, &ldu, Uout,
  //           &ldvT, work, &lwork, rwork, &info);

  //   // std::cout << "info=" << info << " zero is ok!" << std::endl;

  //   for (std::size_t k = 0; k < K; k++)
  //     // std::cout << S[k] << "compared with" << Sout[k] << std::endl;
  //     BOOST_CHECK(abs(S[k] - Sout[k]) < 1e-5);

  //   for (std::size_t i = 0; i < N; i++)
  //     for (std::size_t k = 0; k < K; k++)
  //       // std::cout << U[i*K+k] << "compared with" << Uout[i*K+k] << std::endl;
  //       if (!((k != 3 && i != N) || (i == 0 && k == 1))) // some freedom in SVD
  //         BOOST_CHECK(abs(U[i * K + k] - Uout[i * K + k]) < 1e-5);

    
  //   for (std::size_t i = 0; i < N; i++)
  //     for (std::size_t j = 0; j < M; j++) {
  //       znum_t sum = 0;
  //       for (std::size_t k = 0; k < K; k++)
  //         sum += Uout[i * K + k] * Sout[k] * VDout[k * M + j];
  //       std::cout << "A[i*5+j]=" << A[i*M+j] << " and the sum gives:" << sum
  //       << std::endl;
  //       BOOST_CHECK(abs(A[i * M + j] - sum) < 1e-5);
  //     };
  //   */

    
  //   // for (std::size_t j = 0; j < M; j++)
  //   //   for (std::size_t k = 0; k < K; k++) {
  //   //     std::cout << Vd[k * M + j] << "compared with" << VDout[k * M + j]
  //   //               << std::endl;
  //   //     BOOST_CHECK(abs(Vd[k * M + j] - VDout[k * M + j]) < 1e-5);
  //   //   }
  //   // for (std::size_t i = 0; i < N; i++)
  //   //   for (std::size_t j = 0; j < M; j++) {
  //   //     znum_t sum = 0;
  //   //     for (std::size_t k = 0; k < K; k++)
  //   //       sum += Uout[i * K + k] * Sout[k] * Vd[k * M + j];
  //   //     // std::cout << "A[i*5+j]=" << A[i*M+j] << " and the sum gives:" << sum
  //   //     // << std::endl;
  //   //     BOOST_CHECK(abs(A[i * M + j] - sum) < 1e-5);
  //     // };

  // }
  // { // real, column major
  //   const znum_t A[N * M] = {{1., +1.}, {1., -1.}, {-1., -1.}, {1., 0.}};

  //   znum_t U[N * K] = {{-0.36514837, -7.30296743e-01},
  //                      {-0.18257419, +5.47722558e-01},
  //                      {-0.40824829, +4.08248290e-01},
  //                      {-0.81649658, +2.77555756e-16}};
  //   double S[K] = {2.23606798, 1.41421356};
  //   znum_t Vd[K * M] = {{-0.81649658, 0.},
  //                       {-0.57735027, 0.},
  //                       {0.40824829, -0.40824829},
  //                       {-0.57735027, +0.57735027}};

  //   for (std::size_t i = 0; i < N; i++)
  //     for (std::size_t j = 0; j < M; j++) {
  //       znum_t sum = 0;
  //       for (std::size_t k = 0; k < K; k++)
  //         sum += U[i + k * N] * S[k] * Vd[k + j * K];
  //       BOOST_CHECK(abs(A[i + j * N] - sum) < 1e-5);
  //     };

  //   znum_t Uout[N * K], VDout[K * M];
  //   double Sout[K];
  //   std::size_t ldA = M, ldu = min(M,N), ldvT = K;
  //   znum_t worktest;
  //   double rwork[5 * min(M, N)];
  //   int info, lwork = -1;

  //   zgesvd_((char *)"S", (char *)"S", &M, &N, A, &ldA, Sout, Uout,
  //           &ldu,  VDout, &ldvT, &worktest, &lwork, rwork, &info);

  //   lwork = (int)worktest.real();
  //   znum_t work[lwork];
  //   std::cout << "info=" << info << " zero is ok!" << std::endl;

  //   zgesvd_((char *)"S", (char *)"S", &M, &N, A, &ldA, Sout, Uout, &ldu, VDout, &ldvT,
  //            work, &lwork, rwork, &info);
  //   std::cout << "info=" << info << " zero is ok!" << std::endl;

  //   for (std::size_t k = 0; k < K; k++)
  //     BOOST_CHECK(S[k] - Sout[k] < 1e-5);

  //   for (std::size_t i = 0; i < N; i++)
  //     for (std::size_t j = 0; j < M; j++) {
  //       znum_t sum = 0;
  //       for (std::size_t k = 0; k < K; k++)
  //         sum += Uout[i + k * N] * Sout[k] * VDout[k + j * K];
  //       std::cout << "A[i*5+j]=" << A[i*M+j] << " and the sum gives:" << sum
  //       << std::endl;
  //       BOOST_CHECK(abs(A[i + j * N] - sum) < 1e-5);
  //     };


  // }
}


BOOST_AUTO_TEST_CASE(test_simple_zgesvd_overwrite_right) {

  const std::size_t N = 2;
  const std::size_t K = 2;
  const std::size_t M = 2;

  // { // real, row major
  //   const znum_t A[N * M] = {{1., +1.}, {-1., -1.}, {1., -1.}, {1., 0.}};

  //   znum_t U[N * K] = {{-0.36514837, -7.30296743e-01},
  //                      {-0.40824829, +4.08248290e-01},
  //                      {-0.18257419, +5.47722558e-01},
  //                      {-0.81649658, +2.77555756e-16}};
  //   double S[K] = {2.23606798, 1.41421356};
  //   znum_t Vd[K * M] = {{-0.81649658, 0.},
  //                       {0.40824829, -0.40824829},
  //                       {-0.57735027, 0.},
  //                       {-0.57735027, +0.57735027}};
  //   for (std::size_t i = 0; i < N; i++)
  //     for (std::size_t j = 0; j < M; j++) {
  //       znum_t sum = 0;
  //       for (std::size_t k = 0; k < K; k++)
  //         sum += U[i * K + k] * S[k] * Vd[k * M + j];
  //       // std::cout << "A[i*5+j]=" << A[i*M+j] << " and the sum gives:" << sum
  //       // << std::endl;
  //       BOOST_CHECK(abs(A[i * M + j] - sum) < 1e-5);
  //     };
  //   /*

  //   znum_t Uout[N * K], VDout[K * M];
  //   double Sout[K];
  //   std::size_t ldA = M, ldu = N, ldvT = M < N ? M : N;
  //   znum_t worktest;
  //   double rwork[5 * min(M, N)];
  //   int info, lwork = -1;

  //   zgesvd_((char *)"S", (char *)"S", &M, &N, A, &ldA, Sout, VDout, &ldu, Uout,
  //           &ldvT, &worktest, &lwork, rwork, &info);

  //   lwork = (int)worktest.real();
  //   znum_t work[lwork];

  //   zgesvd_((char *)"S", (char *)"S", &M, &N, A, &ldA, Sout, VDout, &ldu, Uout,
  //           &ldvT, work, &lwork, rwork, &info);

  //   // std::cout << "info=" << info << " zero is ok!" << std::endl;

  //   for (std::size_t k = 0; k < K; k++)
  //     // std::cout << S[k] << "compared with" << Sout[k] << std::endl;
  //     BOOST_CHECK(abs(S[k] - Sout[k]) < 1e-5);

  //   for (std::size_t i = 0; i < N; i++)
  //     for (std::size_t k = 0; k < K; k++)
  //       // std::cout << U[i*K+k] << "compared with" << Uout[i*K+k] << std::endl;
  //       if (!((k != 3 && i != N) || (i == 0 && k == 1))) // some freedom in SVD
  //         BOOST_CHECK(abs(U[i * K + k] - Uout[i * K + k]) < 1e-5);

    
  //   for (std::size_t i = 0; i < N; i++)
  //     for (std::size_t j = 0; j < M; j++) {
  //       znum_t sum = 0;
  //       for (std::size_t k = 0; k < K; k++)
  //         sum += Uout[i * K + k] * Sout[k] * VDout[k * M + j];
  //       std::cout << "A[i*5+j]=" << A[i*M+j] << " and the sum gives:" << sum
  //       << std::endl;
  //       BOOST_CHECK(abs(A[i * M + j] - sum) < 1e-5);
  //     };
  //   */

    
  //   // for (std::size_t j = 0; j < M; j++)
  //   //   for (std::size_t k = 0; k < K; k++) {
  //   //     std::cout << Vd[k * M + j] << "compared with" << VDout[k * M + j]
  //   //               << std::endl;
  //   //     BOOST_CHECK(abs(Vd[k * M + j] - VDout[k * M + j]) < 1e-5);
  //   //   }
  //   // for (std::size_t i = 0; i < N; i++)
  //   //   for (std::size_t j = 0; j < M; j++) {
  //   //     znum_t sum = 0;
  //   //     for (std::size_t k = 0; k < K; k++)
  //   //       sum += Uout[i * K + k] * Sout[k] * Vd[k * M + j];
  //   //     // std::cout << "A[i*5+j]=" << A[i*M+j] << " and the sum gives:" << sum
  //   //     // << std::endl;
  //   //     BOOST_CHECK(abs(A[i * M + j] - sum) < 1e-5);
  //     // };

  // }
  // { // real, column major
  //   const znum_t A[N * M] = {{1., +1.}, {1., -1.}, {-1., -1.}, {1., 0.}};

  //   znum_t U[N * K] = {{-0.36514837, -7.30296743e-01},
  //                      {-0.18257419, +5.47722558e-01},
  //                      {-0.40824829, +4.08248290e-01},
  //                      {-0.81649658, +2.77555756e-16}};
  //   double S[K] = {2.23606798, 1.41421356};
  //   znum_t Vd[K * M] = {{-0.81649658, 0.},
  //                       {-0.57735027, 0.},
  //                       {0.40824829, -0.40824829},
  //                       {-0.57735027, +0.57735027}};

  //   for (std::size_t i = 0; i < N; i++)
  //     for (std::size_t j = 0; j < M; j++) {
  //       znum_t sum = 0;
  //       for (std::size_t k = 0; k < K; k++)
  //         sum += U[i + k * N] * S[k] * Vd[k + j * K];
  //       BOOST_CHECK(abs(A[i + j * N] - sum) < 1e-5);
  //     };

  //   znum_t Uout[N * K], VDout[K * M];
  //   double Sout[K];
  //   std::size_t ldA = M, ldu = min(M,N), ldvT = K;
  //   znum_t worktest;
  //   double rwork[5 * min(M, N)];
  //   int info, lwork = -1;

  //   zgesvd_((char *)"S", (char *)"S", &M, &N, A, &ldA, Sout, Uout,
  //           &ldu,  VDout, &ldvT, &worktest, &lwork, rwork, &info);

  //   lwork = (int)worktest.real();
  //   znum_t work[lwork];
  //   std::cout << "info=" << info << " zero is ok!" << std::endl;

  //   zgesvd_((char *)"S", (char *)"S", &M, &N, A, &ldA, Sout, Uout, &ldu, VDout, &ldvT,
  //            work, &lwork, rwork, &info);
  //   std::cout << "info=" << info << " zero is ok!" << std::endl;

  //   for (std::size_t k = 0; k < K; k++)
  //     BOOST_CHECK(S[k] - Sout[k] < 1e-5);

  //   for (std::size_t i = 0; i < N; i++)
  //     for (std::size_t j = 0; j < M; j++) {
  //       znum_t sum = 0;
  //       for (std::size_t k = 0; k < K; k++)
  //         sum += Uout[i + k * N] * Sout[k] * VDout[k + j * K];
  //       std::cout << "A[i*5+j]=" << A[i*M+j] << " and the sum gives:" << sum
  //       << std::endl;
  //       BOOST_CHECK(abs(A[i + j * N] - sum) < 1e-5);
  //     };


  // }
}

#undef znum_t
