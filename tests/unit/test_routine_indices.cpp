#include <boost/test/unit_test.hpp>

#include "mdot/include/babel_type.hpp"
#include "mdot/include/routine/indices.hpp"

BOOST_AUTO_TEST_CASE(test_routine_indices_indices_dst_theta_no_gate) {
  std::vector<m_index_t> left_indices = {{0, 0, 0}, {0, 1, 0}};
  std::vector<m_index_t> right_indices = {{0, 0, 0}, {0, 1, 0}};
  auto out = mdot::indices_dst_theta_no_gate(left_indices, right_indices, true);
  BOOST_CHECK(std::get<0>(std::get<0>(out[0])) == 0);
  auto pair_out = mdot::split_degenerate_indices(out);
  BOOST_CHECK(std::get<0>(std::get<0>(pair_out.first[0])) == 0);
}

BOOST_AUTO_TEST_CASE(test_routine_indices_indices_sum_sub) {
  constexpr index_t lhs = 3;
  constexpr index_t rhs = 1;

  BOOST_CHECK(mdot::internal_qn_sum(lhs, rhs) == 4);
  BOOST_CHECK(mdot::internal_qn_sub(lhs, rhs) == 2);

  index_t o_lhs = 5;
  index_t o_rhs = 2;

  BOOST_CHECK(mdot::internal_qn_sum(o_lhs, o_rhs) == 7);
  BOOST_CHECK(mdot::internal_qn_sub(o_lhs, o_rhs) == 3);
}

BOOST_AUTO_TEST_CASE(test_routine_indices_potential_middle_indices) {
  std::vector<t_index_t> theta_indices = {{0, 0, 0,0}, {0, 1, 0,0}};
  auto out = mdot::potential_middle_indices(theta_indices, -1);
  BOOST_CHECK(out[0] == 0);
  BOOST_CHECK(out[1] == 1);
}