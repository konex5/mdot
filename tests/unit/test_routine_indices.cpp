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
  std::vector<t_index_t> theta_indices = {
      {0, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}};
  auto middle = mdot::potential_middle_indices(theta_indices, 1);
  BOOST_CHECK(middle[0] == 0);
  BOOST_CHECK(middle[1] == 1);
  auto out = mdot::degeneracy_in_theta(theta_indices, middle, 1);
  BOOST_CHECK(std::get<0>(out.first[0].second) == 0);
  BOOST_CHECK(std::get<1>(out.first[0].second) == 1);
  BOOST_CHECK(std::get<2>(out.first[0].second) == 0);
  BOOST_CHECK(std::get<3>(out.first[0].second) == 0);
  BOOST_CHECK(std::get<0>(out.second[0].second[0]) == 0);
  BOOST_CHECK(std::get<1>(out.second[0].second[0]) == 0);
  BOOST_CHECK(std::get<2>(out.second[0].second[0]) == 0);
  BOOST_CHECK(std::get<3>(out.second[0].second[0]) == 0);
  BOOST_CHECK(std::get<0>(out.second[0].second[1]) == 0);
  BOOST_CHECK(std::get<1>(out.second[0].second[1]) == 0);
  BOOST_CHECK(std::get<2>(out.second[0].second[1]) == 1);
  BOOST_CHECK(std::get<3>(out.second[0].second[1]) == 0);
}

BOOST_AUTO_TEST_CASE(test_routine_indices_slices_degenerate_blocs) {
  std::vector<t_index_t> theta_indices_large = {
      {0, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 0, 2}, {0, 0, 0, 4}, {0, 0, 1, 0},
      {0, 0, 1, 1}, {0, 0, 1, 2}, {0, 0, 1, 4}, {0, 1, 0, 0}, {0, 1, 0, 1},
      {0, 1, 1, 0}, {0, 1, 1, 1}, {0, 1, 1, 2}, {0, 1, 1, 4}, {1, 0, 0, 3},
      {1, 0, 1, 0}, {1, 1, 0, 0}, {1, 1, 1, 3}};

  auto middle = mdot::potential_middle_indices(theta_indices_large, 1);
  BOOST_CHECK(middle.size() == 6);
  auto out_nondeg_deg =
      mdot::degeneracy_in_theta(theta_indices_large, middle, 1);
  BOOST_CHECK(out_nondeg_deg.first.size() == 0);
  BOOST_CHECK(out_nondeg_deg.second.size() == 3);
  // start of test
  std::vector<t_shape_t> theta_shape_large = {
      {1, 1, 1, 2}, {1, 1, 1, 3}, {1, 1, 1, 4}, {1, 1, 1, 8}, {1, 1, 1, 2},
      {1, 1, 1, 3}, {1, 1, 1, 4}, {1, 1, 1, 8}, {1, 1, 1, 2}, {1, 1, 1, 3},
      {1, 1, 1, 2}, {1, 1, 1, 3}, {1, 1, 1, 4}, {1, 1, 1, 8}, {3, 1, 1, 6},
      {3, 1, 1, 2}, {3, 1, 1, 2}, {3, 1, 1, 6}};
  dtbloc_t empty_theta;
  for (size_t i = 0; i < theta_indices_large.size(); i++)
    empty_theta[theta_indices_large[i]] = {theta_shape_large[i],
                                           std::vector<dnum_t>(0)};

  auto new_subsize =
      mdot::slices_degenerate_blocs(empty_theta, out_nondeg_deg.second);
}