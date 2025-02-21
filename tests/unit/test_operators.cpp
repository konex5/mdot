#include <boost/test/unit_test.hpp>

#include "mdot/include/babel_type.hpp"
#include "mdot/include/operators.hpp"

inline darr_t mul_twotwo(darr_t a, darr_t b) {
  darr_t c(4);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++) {
      double sum = 0;
      for (int k = 0; k < 2; k++)
        sum += a[i * 2 + k] * b[k * 2 + j];
      c[i + 2 * j] = sum;
    }
  return c;
}

BOOST_AUTO_TEST_CASE(test_operators) {

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
  //
  auto shsp = real_single_operator("sh-sp", "sh-none");
  BOOST_CHECK((shsp[{0, 0}])[1 + 0 * 2] == 1);
  auto shsm = real_single_operator("sh-sm", "sh-none");
  BOOST_CHECK((shsm[{0, 0}])[0 + 1 * 2] == 1);
  //
  auto shspsm = mul_twotwo(shsp[{0, 0}], shsm[{0, 0}]);
  BOOST_CHECK(shspsm[0 + 0 * 2] == 1);
}
