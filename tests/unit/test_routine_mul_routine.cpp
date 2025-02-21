#include <boost/test/unit_test.hpp>

#include "mdot/include/babel_type.hpp"
#include "mdot/include/routine/mul_routine.hpp"

#include <iostream>

BOOST_AUTO_TEST_CASE(test_routine_mul_routine) {
  const size_t N = 2;
  const size_t M = 3;
  std::vector<dnum_t> matrix = {1, 2, 3, 4, 5, 6};
  const size_t cut = 2;
  std::vector<dnum_t> diagonal = {2, -1, 3};
  {
    std::vector<dnum_t> destination(N * cut);
    mdot::mul_mat_diag(destination, matrix, N, M, diagonal, cut);
    BOOST_CHECK(destination[0] == 2);
    BOOST_CHECK(destination[1] == -2);
    BOOST_CHECK(destination[2] == 8);
    BOOST_CHECK(destination[3] == -5);
    BOOST_CHECK(destination.size() == 4);
  }
  {
    std::vector<dnum_t> destination(cut * M);
    mdot::mul_diag_mat(destination, matrix, N, M, diagonal, cut);
    BOOST_CHECK(destination[0] == 2);
    BOOST_CHECK(destination[1] == 4);
    BOOST_CHECK(destination[2] == 6);
    BOOST_CHECK(destination[3] == -4);
    BOOST_CHECK(destination[4] == -5);
    BOOST_CHECK(destination[5] == -6);
    BOOST_CHECK(destination.size() == 6);
  }
}