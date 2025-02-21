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
    auto indices = mdot::real_sh_blocs_crtp<mdot::sh_id_u1>::get_indices();
    BOOST_CHECK(std::get<0>(indices[0]) == 0);
    BOOST_CHECK(std::get<1>(indices[0]) == 0);
    auto index = mdot::real_sh_blocs_crtp<mdot::sh_id_u1>::get_index(1);
    BOOST_CHECK(std::get<0>(index) == 1);
    BOOST_CHECK(std::get<1>(index) == 1);
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
    auto indices = mdot::real_sh_blocs_crtp<mdot::sh_id_cplx_u1>::get_indices();
    BOOST_CHECK(std::get<1>(indices[0]) == 0);
    BOOST_CHECK(std::get<0>(indices[1]) == 1);
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
    auto index = mdot::real_sh_blocs_crtp<mdot::sh_sp_u1>::get_index(0);
    BOOST_CHECK(std::get<0>(index) == 0);
    BOOST_CHECK(std::get<1>(index) == 1);
    auto array = mdot::real_sh_blocs_crtp<mdot::sh_sp_u1>::get_array(0, 1.);
    BOOST_CHECK(array[0] == 2.);
  }
  {
    // sh_sm_u1
    BOOST_CHECK(mdot::real_sh_blocs_crtp<mdot::sh_sm_u1>::nb_blocs == 1);
    BOOST_CHECK(mdot::real_sh_blocs_crtp<mdot::sh_sm_u1>::get_size(0) == 1);
    BOOST_CHECK(
        std::get<0>(
            mdot::real_sh_blocs_crtp<mdot::sh_sm_u1>::get_shapes()[0]) == 1);
    auto index = mdot::real_sh_blocs_crtp<mdot::sh_sm_u1>::get_index(0);
    BOOST_CHECK(std::get<0>(index) == 1);
    BOOST_CHECK(std::get<1>(index) == 0);
    auto array = mdot::real_sh_blocs_crtp<mdot::sh_sm_u1>::get_array(0);
    BOOST_CHECK(array[0] == 2.);
  }
  {
    // sh_sx_u1
    BOOST_CHECK(mdot::real_sh_blocs_crtp<mdot::sh_sx_u1>::nb_blocs == 2);
    BOOST_CHECK(mdot::real_sh_blocs_crtp<mdot::sh_sx_u1>::get_size(0) == 1);
    BOOST_CHECK(
        std::get<0>(
            mdot::real_sh_blocs_crtp<mdot::sh_sx_u1>::get_shapes()[0]) == 1);
    auto indices = mdot::real_sh_blocs_crtp<mdot::sh_sx_u1>::get_indices();
    BOOST_CHECK(std::get<0>(indices[0]) == 0);
    BOOST_CHECK(std::get<1>(indices[0]) == 1);
    BOOST_CHECK(std::get<0>(indices[1]) == 1);
    BOOST_CHECK(std::get<1>(indices[1]) == 0);
    auto arrays = mdot::real_sh_blocs_crtp<mdot::sh_sx_u1>::get_arrays();
    BOOST_CHECK(arrays[0][0] == 1.);
    BOOST_CHECK(arrays[1][0] == 1.);
  }
  {
    // sh_sy_u1
    BOOST_CHECK(mdot::cplx_sh_blocs_crtp<mdot::sh_sy_u1>::nb_blocs == 2);
    BOOST_CHECK(mdot::cplx_sh_blocs_crtp<mdot::sh_sy_u1>::get_size(0) == 1);
    BOOST_CHECK(
        std::get<0>(
            mdot::cplx_sh_blocs_crtp<mdot::sh_sy_u1>::get_shapes()[0]) == 1);
    auto indices = mdot::cplx_sh_blocs_crtp<mdot::sh_sy_u1>::get_indices();
    BOOST_CHECK(std::get<0>(indices[0]) == 0);
    BOOST_CHECK(std::get<1>(indices[0]) == 1);
    BOOST_CHECK(std::get<0>(indices[1]) == 1);
    BOOST_CHECK(std::get<1>(indices[1]) == 0);
    auto arrays = mdot::cplx_sh_blocs_crtp<mdot::sh_sy_u1>::get_arrays();
    BOOST_CHECK(arrays[0][0].imag() == -1.);
    BOOST_CHECK(arrays[1][0].imag() == 1.);
  }
  {
    // sh_sz_u1
    BOOST_CHECK(mdot::real_sh_blocs_crtp<mdot::sh_sz_u1>::nb_blocs == 2);
    BOOST_CHECK(mdot::real_sh_blocs_crtp<mdot::sh_sz_u1>::get_size(0) == 1);
    BOOST_CHECK(
        std::get<0>(
            mdot::real_sh_blocs_crtp<mdot::sh_sz_u1>::get_shapes()[0]) == 1);
    auto indices = mdot::real_sh_blocs_crtp<mdot::sh_sz_u1>::get_indices();
    BOOST_CHECK(std::get<0>(indices[0]) == 0);
    BOOST_CHECK(std::get<1>(indices[0]) == 0);
    BOOST_CHECK(std::get<0>(indices[1]) == 1);
    BOOST_CHECK(std::get<1>(indices[1]) == 1);
    auto arrays = mdot::real_sh_blocs_crtp<mdot::sh_sz_u1>::get_arrays();
    BOOST_CHECK(arrays[0][0] == 1.);
    BOOST_CHECK(arrays[1][0] == -1.);
  }
}

BOOST_AUTO_TEST_CASE(test_blocs_static_sh_real_blocs) {
  auto blocs = mdot::real_sh_blocs_crtp<mdot::sh_id_u1>::get_blocs(1.);
  BOOST_CHECK(std::get<0>(blocs[{1, 1}]) == 1);
  BOOST_CHECK(std::get<0>(std::get<1>(blocs[{1, 1}])) == 1);
  BOOST_CHECK(std::get<1>(std::get<1>(blocs[{1, 1}])) == 1);
  BOOST_CHECK(std::get<2>(blocs[{1, 1}])[0] == 1);
}
