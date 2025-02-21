#include <boost/test/unit_test.hpp>
// #include "boost/list.hpp"

#include "mdot/utils/qnum.hpp"

using namespace mdot::qnum;

BOOST_AUTO_TEST_CASE(unittest_qnum_sh) {
  BOOST_CHECK(quantum_number_crtp<sh_none>::dim == 1);
  BOOST_CHECK(quantum_number_crtp<sh_none>::sum(2, 3) == 0);
  BOOST_CHECK(quantum_number_crtp<sh_none>::sum(7, 2) == 0);

  BOOST_CHECK(quantum_number_crtp<sh_u1>::dim == 2);
  BOOST_CHECK(quantum_number_crtp<sh_u1>::sum(0, 0) == 0);
  BOOST_CHECK(quantum_number_crtp<sh_u1>::sum(1, 7) == 8);
  BOOST_CHECK(quantum_number_crtp<sh_u1>::sum(4, 9) == 13);

  BOOST_CHECK(quantum_number_crtp<sh_su2>::dim == 1);
  BOOST_CHECK(quantum_number_crtp<sh_su2>::sum(0, 0) == 0);
  BOOST_CHECK(quantum_number_crtp<sh_su2>::sum(6, 7) == 13);
  BOOST_CHECK(quantum_number_crtp<sh_su2>::sum(4, 3) == 7);
}
