#include <boost/test/unit_test.hpp>

#include "mdot/include/babel_type.hpp"
#include "mdot/include/routine/eig_lanczos.hpp"

#include <iostream>

BOOST_AUTO_TEST_CASE(test_routine_mul_routine) {
  const size_t N = 3;
  const size_t M = 3;

  std::vector<dnum_t> matrix = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  std::vector<dnum_t> psi = {0.3, 0.3, 0.3};

  size_t maxiter = 20;
  const dnum_t tolerance = 1e-1;
  dnum_t eigenvalue;
  dnum_t eigenvector[3];
  bool converged;

  converged = mdot::lanczos_ev(matrix.data(), psi.data(), N, maxiter, tolerance, eigenvalue,
                   eigenvector);
  BOOST_CHECK(converged);      
  BOOST_CHECK(eigenvalue == -42.989621598238621);
}