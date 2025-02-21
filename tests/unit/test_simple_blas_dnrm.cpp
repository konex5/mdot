#include <boost/test/unit_test.hpp>

#include <math.h>

using namespace std;

#define size_t typename std::size_t
#define dnum_t double

extern "C" {
  dnum_t dnrm2_(const size_t* n,  const dnum_t *a, const size_t* incX);
}

#undef size_t
#undef dnum_t

BOOST_AUTO_TEST_CASE(test_dnorm_simple) {

  { // real, row major
    const std::size_t N = 3;
    const double A[N] = {0, 2, 0};
    const double B[N] = {-1.2, 1.3, 2.7};
    
    double norm_out;
    const std::size_t incX=1;

    norm_out = dnrm2_(&N,A,&incX);
    
    BOOST_CHECK(norm_out == 2);
    
    norm_out = dnrm2_(&N,B,&incX);
    BOOST_CHECK(norm_out == 3.228002478313795);
    
  }
}