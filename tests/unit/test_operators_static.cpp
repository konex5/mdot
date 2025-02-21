#include <boost/test/unit_test.hpp>

//#include "mdot/include/babel_type.hpp"
#include "mdot/include/operators_static.hpp"

BOOST_AUTO_TEST_CASE(test_operators_static_so_none_real) {
  // sh_id_no
  BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_id_no>::size == 4);
  BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_id_no>::shape[0] == 2);
  BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_id_no>::shape[1] == 2);
  BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_id_no>::array[0] == 1);
  BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_id_no>::array[3] == 1);
  BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_id_no>::times(2)[0] == 2);
  BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_id_no>::square()[0] == 1);
  // sh_id_no cplx
  BOOST_CHECK(mdot::cplx_operators_crtp<mdot::sh_id_cplx_no>::size == 4);
  BOOST_CHECK(mdot::cplx_operators_crtp<mdot::sh_id_cplx_no>::shape[0] == 2);
  BOOST_CHECK(mdot::cplx_operators_crtp<mdot::sh_id_cplx_no>::shape[1] == 2);
  BOOST_CHECK(mdot::cplx_operators_crtp<mdot::sh_id_cplx_no>::array[0].real() == 1);
  BOOST_CHECK(mdot::cplx_operators_crtp<mdot::sh_id_cplx_no>::array[3].real() == 1);
  BOOST_CHECK(mdot::cplx_operators_crtp<mdot::sh_id_cplx_no>::times(2)[0].real() == 2);
  BOOST_CHECK(mdot::cplx_operators_crtp<mdot::sh_id_cplx_no>::square()[0].real() == 1);
  
/* 
  auto soid = std::get<0>(single_operator_real("so-id", "so-none"));
  auto soid_norm = std::get<1>(single_operator_real("so-id", "so-none"));

  BOOST_CHECK((soid[{0, 0}])[0 + 0 * 3] == 1);
  BOOST_CHECK((soid[{0, 0}])[1 + 1 * 3] == 1);
  BOOST_CHECK((soid[{0, 0}])[2 + 2 * 3] == 1);
  BOOST_CHECK((soid[{0, 0}])[0 + 1 * 3] == 0);

  auto soid_square = mul_threethree(soid[{0, 0}], soid[{0, 0}]);

  BOOST_CHECK(std::abs((soid_square[0 + 0 * 3] + soid_square[1 + 1 * 3] +
                        soid_square[2 + 2 * 3]) *
                           pow(soid_norm, 2) -
                       1) < 1e-7); */
}


BOOST_AUTO_TEST_CASE(test_operators_static_so_none_cplx) {
  BOOST_CHECK(mdot::cplx_operators_crtp<mdot::sh_sy_no>::size == 4);
  BOOST_CHECK(mdot::cplx_operators_crtp<mdot::sh_sy_no>::shape[0] == 2);
  BOOST_CHECK(mdot::cplx_operators_crtp<mdot::sh_sy_no>::shape[1] == 2);
  BOOST_CHECK(mdot::cplx_operators_crtp<mdot::sh_sy_no>::array[1].imag() == 1);
  BOOST_CHECK(mdot::cplx_operators_crtp<mdot::sh_sy_no>::array[3].real() == 0);

}