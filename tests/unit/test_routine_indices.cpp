#include <boost/test/unit_test.hpp>

#include "mdot/include/babel_type.hpp"
#include "mdot/include/routine/indices.hpp"

BOOST_AUTO_TEST_CASE(test_routine_indices_indices_dst_theta_no_gate) {
  std::vector<m_index_t> left_indices = {{0, 0, 0}, {0, 1, 0}};
  std::vector<m_index_t> right_indices = {{0, 0, 0}, {0, 1, 0}};
  auto out = mdot::indices_dst_theta_no_gate(left_indices, right_indices, true);
  BOOST_CHECK(std::get<0>(std::get<0>(out[0])) == 0);
}