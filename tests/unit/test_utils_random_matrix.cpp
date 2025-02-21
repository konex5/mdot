#include <boost/test/unit_test.hpp>

#include "utils/random_matrix.hpp"

BOOST_AUTO_TEST_CASE(test_utils_random_matrix) {
  auto m = utils::random_matrix();
  BOOST_CHECK(m[0] != m[1]);
  std::cout << m[0] << m[1];
  auto ma = utils::normalize(m);
  for (std::size_t i = 0; i < ma.size(); i++)
    BOOST_CHECK(ma[i] < 1);
  std::cout << ma[0] << ma[1];
}
