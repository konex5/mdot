
#include <boost/test/unit_test.hpp>

#include "mdot/utils/babel_type.hpp"

BOOST_AUTO_TEST_CASE(unittest_babel_type) {
  using d = real_type<double>;
  bool a = std::is_same<d, double>::value;
  BOOST_CHECK(true);
  BOOST_CHECK(!a);
  // https://stackoverflow.com/questions/67116949/c-sfinae-partial-specialization
}
