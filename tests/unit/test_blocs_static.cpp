#include <boost/test/unit_test.hpp>

#include "mdot/include/blocs_static.hpp"

BOOST_AUTO_TEST_CASE(test_blocs_static_sh_real) {
  {
    // sh_id_u1
    BOOST_CHECK(mdot::real_sh_blocs_crtp<mdot::sh_id_u1>::nb_blocs == 2);
    BOOST_CHECK(mdot::real_sh_blocs_crtp<mdot::sh_id_u1>::get_size(0) == 1);
    BOOST_CHECK(mdot::real_sh_blocs_crtp<mdot::sh_id_u1>::get_size(1) == 1);
    BOOST_CHECK(std::get<0>(mdot::real_sh_blocs_crtp<mdot::sh_id_u1>::get_shape(
                    0)) == 1);
    BOOST_CHECK(std::get<1>(mdot::real_sh_blocs_crtp<mdot::sh_id_u1>::get_shape(
                    0)) == 1);
    BOOST_CHECK(std::get<0>(mdot::real_sh_blocs_crtp<mdot::sh_id_u1>::get_shape(
                    1)) == 1);
    BOOST_CHECK(std::get<1>(mdot::real_sh_blocs_crtp<mdot::sh_id_u1>::get_shape(
                    1)) == 1);
    auto a = mdot::real_sh_blocs_crtp<mdot::sh_id_u1>::get_shapes();
    BOOST_CHECK(std::get<0>(a[0]) == 1);
    BOOST_CHECK(std::get<1>(a[0]) == 1);
    BOOST_CHECK(std::get<0>(a[1]) == 1);
    BOOST_CHECK(std::get<1>(a[1]) == 1);
    auto arrays = mdot::real_sh_blocs_crtp<mdot::sh_id_u1>::get_arrays(2.);
    BOOST_CHECK(arrays[0][0] == 2.);
    BOOST_CHECK(arrays[0][1] == 2.);
    auto array1 = mdot::real_sh_blocs_crtp<mdot::sh_id_u1>::get_array(1, -2.);
    BOOST_CHECK(array1[0] == -2.);
  }
  {
    // sh_id_u1
    BOOST_CHECK(mdot::cplx_sh_blocs_crtp<mdot::sh_id_cplx_u1>::nb_blocs == 2);
    BOOST_CHECK(mdot::cplx_sh_blocs_crtp<mdot::sh_id_cplx_u1>::get_size(0) ==
                1);
    BOOST_CHECK(mdot::cplx_sh_blocs_crtp<mdot::sh_id_cplx_u1>::get_size(1) ==
                1);
    BOOST_CHECK(
        std::get<0>(
            mdot::cplx_sh_blocs_crtp<mdot::sh_id_cplx_u1>::get_shape(0)) == 1);
    BOOST_CHECK(
        std::get<1>(
            mdot::cplx_sh_blocs_crtp<mdot::sh_id_cplx_u1>::get_shape(0)) == 1);
    BOOST_CHECK(
        std::get<0>(
            mdot::cplx_sh_blocs_crtp<mdot::sh_id_cplx_u1>::get_shape(1)) == 1);
    BOOST_CHECK(
        std::get<1>(
            mdot::cplx_sh_blocs_crtp<mdot::sh_id_cplx_u1>::get_shape(1)) == 1);
    auto a = mdot::cplx_sh_blocs_crtp<mdot::sh_id_cplx_u1>::get_shapes();
    BOOST_CHECK(std::get<0>(a[0]) == 1);
    BOOST_CHECK(std::get<1>(a[0]) == 1);
    BOOST_CHECK(std::get<0>(a[1]) == 1);
    BOOST_CHECK(std::get<1>(a[1]) == 1);
    auto arrays = mdot::cplx_sh_blocs_crtp<mdot::sh_id_cplx_u1>::get_arrays(2.);
    BOOST_CHECK(arrays[0][0].real() == 2.);
    BOOST_CHECK(arrays[0][1].real() == 2.);
    auto array1 =
        mdot::cplx_sh_blocs_crtp<mdot::sh_id_cplx_u1>::get_array(1, -2.);
    BOOST_CHECK(array1[0].real() == -2.);
  }
  {
    // sh_sp_u1
    BOOST_CHECK(mdot::real_sh_blocs_crtp<mdot::sh_sp_u1>::nb_blocs == 1);
    BOOST_CHECK(mdot::real_sh_blocs_crtp<mdot::sh_sp_u1>::get_size(0) == 1);
    auto a = mdot::real_sh_blocs_crtp<mdot::sh_sp_u1>::get_shapes();
    BOOST_CHECK(std::get<0>(a[0]) == 1);
    auto array = mdot::real_sh_blocs_crtp<mdot::sh_sp_u1>::get_array(0, 1.);
    BOOST_CHECK(array[0] == 1.);
    // auto array1 = mdot::real_sh_blocs_crtp<mdot::sh_id_u1>::get_array(0,3.);
    // BOOST_CHECK(array1[0] == 3);
  }
  // BOOST_CHECK(mdot::real_sh_blocs_crtp<mdot::sh_id_u1>::sizes[0] == 1);
  // BOOST_CHECK(mdot::real_sh_blocs_crtp<mdot::sh_id_u1>::sizes[0] == 1);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_id_no>::shape[0] == 2);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_id_no>::shape[1] == 2);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_id_no>::array[0] == 1);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_id_no>::array[3] == 1);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_id_no>::times(2)[0] == 2);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_id_no>::square()[0] == 1);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_id_no>::trace() == 2);
  // // sh_id_no cplx
  // BOOST_CHECK(mdot::cplx_operators_crtp<mdot::sh_id_cplx_no>::size == 4);
  // BOOST_CHECK(mdot::cplx_operators_crtp<mdot::sh_id_cplx_no>::shape[0] == 2);
  // BOOST_CHECK(mdot::cplx_operators_crtp<mdot::sh_id_cplx_no>::shape[1] == 2);
  // BOOST_CHECK(mdot::cplx_operators_crtp<mdot::sh_id_cplx_no>::array[0].real()
  // == 1);
  // BOOST_CHECK(mdot::cplx_operators_crtp<mdot::sh_id_cplx_no>::array[3].real()
  // == 1);
  // BOOST_CHECK(mdot::cplx_operators_crtp<mdot::sh_id_cplx_no>::times(2)[0].real()
  // == 2);
  // BOOST_CHECK(mdot::cplx_operators_crtp<mdot::sh_id_cplx_no>::square()[0].real()
  // == 1);
  // BOOST_CHECK(mdot::cplx_operators_crtp<mdot::sh_id_cplx_no>::trace().real()
  // == 2);
  // // sh_sp_no
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sp_no>::size == 4);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sp_no>::shape[0] == 2);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sp_no>::shape[1] == 2);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sp_no>::array[1] == 1);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sp_no>::array[0] == 0);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sp_no>::times(2)[1] == 2);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sp_no>::trace() == 0);
  // // sh_sm_no
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sm_no>::size == 4);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sm_no>::shape[0] == 2);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sm_no>::shape[1] == 2);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sm_no>::array[2] == 1);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sm_no>::array[0] == 0);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sm_no>::times(2)[2] == 2);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sm_no>::trace() == 0);
  // // sh_sx_no
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sx_no>::size == 4);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sx_no>::shape[0] == 2);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sx_no>::shape[1] == 2);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sx_no>::array[1] == 1);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sx_no>::array[2] == 1);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sx_no>::array[0] == 0);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sx_no>::times(2)[1] == 2);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sx_no>::trace() == 0);
  // // sh_sy_no
  // BOOST_CHECK(mdot::cplx_operators_crtp<mdot::sh_sy_no>::size == 4);
  // BOOST_CHECK(mdot::cplx_operators_crtp<mdot::sh_sy_no>::shape[0] == 2);
  // BOOST_CHECK(mdot::cplx_operators_crtp<mdot::sh_sy_no>::shape[1] == 2);
  // BOOST_CHECK(mdot::cplx_operators_crtp<mdot::sh_sy_no>::array[1].imag() ==
  // 1); BOOST_CHECK(mdot::cplx_operators_crtp<mdot::sh_sy_no>::array[3].real()
  // == 0);
  // BOOST_CHECK(mdot::cplx_operators_crtp<mdot::sh_sy_no>::trace().real() ==
  // 0);
  // // sh_sz_no
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sz_no>::size == 4);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sz_no>::shape[0] == 2);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sz_no>::shape[1] == 2);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sz_no>::array[0] == 1);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sz_no>::array[3] == -1);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sz_no>::array[1] == 0);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sz_no>::times(2)[3] == -2);
  // BOOST_CHECK(mdot::real_operators_crtp<mdot::sh_sz_no>::trace() == 0);
}
