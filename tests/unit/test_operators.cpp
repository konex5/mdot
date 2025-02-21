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

BOOST_AUTO_TEST_CASE(test_operators_sh_none_real) {

  auto shid = std::get<0>(single_operator_real("sh-id", "sh-none"));
  auto shid_norm = std::get<1>(single_operator_real("sh-id", "sh-none"));

  BOOST_CHECK((shid[{0, 0}])[0 + 0 * 2] == 1);
  BOOST_CHECK((shid[{0, 0}])[1 + 0 * 2] == 0);
  BOOST_CHECK((shid[{0, 0}])[1 + 1 * 2] == 1);
  BOOST_CHECK((shid[{0, 0}])[0 + 1 * 2] == 0);

  auto shid_square = mul_twotwo(shid[{0, 0}], shid[{0, 0}]);
  BOOST_CHECK(shid_square[0 + 0 * 2] == 1);
  BOOST_CHECK(shid_square[1 + 0 * 2] == 0);
  BOOST_CHECK(shid_square[1 + 1 * 2] == 1);
  BOOST_CHECK(shid_square[0 + 1 * 2] == 0);

  BOOST_CHECK(std::abs((shid_square[0 + 0 * 2] + shid_square[1 + 1 * 2]) *
                           pow(shid_norm, 2) -
                       1) < 1e-7);
  //
  auto shsp = std::get<0>(single_operator_real("sh-sp", "sh-none"));
  BOOST_CHECK((shsp[{0, 0}])[1 + 0 * 2] == 1);
  auto shsm = std::get<0>(single_operator_real("sh-sm", "sh-none"));
  BOOST_CHECK((shsm[{0, 0}])[0 + 1 * 2] == 1);
  //
  auto shspsm = mul_twotwo(shsp[{0, 0}], shsm[{0, 0}]);
  BOOST_CHECK(shspsm[0 + 0 * 2] == 1);
  //
  auto shsz = std::get<0>(single_operator_real("sh-sz", "sh-none"));

  BOOST_CHECK((shsz[{0, 0}])[0 + 0 * 2] == 1);
  BOOST_CHECK((shsz[{0, 0}])[1 + 0 * 2] == 0);
  BOOST_CHECK((shsz[{0, 0}])[1 + 1 * 2] == -1);
  BOOST_CHECK((shsz[{0, 0}])[0 + 1 * 2] == 0);

  auto shsz_square = mul_twotwo(shid[{0, 0}], shid[{0, 0}]);
  BOOST_CHECK(shsz_square[0 + 0 * 2] == 1);
  BOOST_CHECK(shsz_square[1 + 0 * 2] == 0);
  BOOST_CHECK(shsz_square[1 + 1 * 2] == 1);
  BOOST_CHECK(shsz_square[0 + 1 * 2] == 0);

  BOOST_CHECK(std::abs((shsz_square[0 + 0 * 2] + shsz_square[1 + 1 * 2]) *
                           pow(shid_norm, 2) -
                       1) < 1e-7);
  //
  /*
  auto shsx = std::get<0>(single_operator_real("sh-sx", "sh-none"));

  BOOST_CHECK((shsx[{0, 0}])[0 + 0 * 2] == 1);
  BOOST_CHECK((shsx[{0, 0}])[1 + 0 * 2] == 0);
  BOOST_CHECK((shsx[{0, 0}])[1 + 1 * 2] == -1);
  BOOST_CHECK((shsx[{0, 0}])[0 + 1 * 2] == 0);

  auto shsx_square = mul_twotwo(shid[{0, 0}], shid[{0, 0}]);
  BOOST_CHECK(shsx_square[0 + 0 * 2] == 0);
  BOOST_CHECK(shsx_square[1 + 0 * 2] == 1);
  BOOST_CHECK(shsx_square[0 + 1 * 2] == 1);
  BOOST_CHECK(shsx_square[1 + 1 * 2] == 0);

  BOOST_CHECK(std::abs((shsx_square[0 + 0 * 2] + shsx_square[1 + 1 * 2]) *
                           pow(shid_norm, 2) -
                       1) < 1e-7);
*/
}

BOOST_AUTO_TEST_CASE(test_operators_sh_none_cplx) {}