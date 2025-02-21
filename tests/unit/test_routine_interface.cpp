#include <boost/test/unit_test.hpp>

#include "mdot/include/babel_type.hpp"
#include "mdot/include/routine/interface.hpp"

BOOST_AUTO_TEST_CASE(test_routine_indices_indices_dst_theta_no_gate) {
  dtbloc_t dst_blocs;
  dmbloc_t lhs_blocs;
  dmbloc_t rhs_blocs;
  mdot::mm_to_theta_no_gate(dst_blocs,lhs_blocs, rhs_blocs, true);
  //BOOST_CHECK(std::get<0>(std::get<0>(out[0])) == 0);
}