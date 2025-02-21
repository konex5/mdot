#include <boost/test/unit_test.hpp>

#include "mdot/include/babel_type.hpp"
#include "mdot/include/operators.hpp"

inline darr_t mul_twotwo(darr_t a, darr_t b) {
  darr_t c;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        c[i + 2 * j] = a[i + 2 * k] * b[k + 2 * j];
  return c;
}

BOOST_AUTO_TEST_CASE(test_indexdw) {

  auto shid = real_single_operator("sh-id", "sh-none");

  BOOST_CHECK((shid[{0, 0}])[0 + 0 * 2] == 1);
  BOOST_CHECK((shid[{0, 0}])[1 + 0 * 2] == 0);
  BOOST_CHECK((shid[{0, 0}])[1 + 1 * 2] == 1);
  BOOST_CHECK((shid[{0, 0}])[0 + 1 * 2] == 0);

  auto shid_square = mul_twotwo(shid[{0, 0}], shid[{0, 0}]);
  BOOST_CHECK(shid_square[0 + 0 * 2] == 1);
  BOOST_CHECK(shid_square[1 + 0 * 2] == 0);
  BOOST_CHECK(shid_square[1 + 1 * 2] == 1);
  BOOST_CHECK(shid_square[0 + 1 * 2] == 0);
}
