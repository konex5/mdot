#include <boost/test/unit_test.hpp>

#include "mdot/include/babel_type.hpp"
#include "mdot/include/routine/svd_routine.hpp"

BOOST_AUTO_TEST_CASE(test_routine_svd_routine_bloc_norm) {
  {
    array_of_dtype a = {{1, 2, 3}, {2, 1}, {1, 1, 2}};
    double norm_out;
    bloc_norm(a, {}, norm_out);
    BOOST_CHECK(norm_out == 5.);
    normalize_the_array(a, {});
    BOOST_CHECK(a[0][0] == 0.2);
  }
  //
  {
    array_of_dtype a = {{1, 2, 3}, {2, 1}, {1, 1, 2}};
    double norm_out;
    bloc_norm(a, {2, 2, 0}, norm_out);
    BOOST_CHECK(norm_out == sqrt(10.));
  }
}