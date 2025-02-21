#include <boost/test/unit_test.hpp>

#include "mdot/utils/utils.hpp"

#include <complex>

using namespace fhmdot::utils;

BOOST_AUTO_TEST_CASE(unittest_utils_precision) {
  BOOST_CHECK(is_float<float>());
  BOOST_CHECK(!is_float<double>());
  BOOST_CHECK(is_float<typename std::complex<float>>());
  BOOST_CHECK(!is_float<typename std::complex<double>>());
  //
  BOOST_CHECK(!is_complex<float>());
  BOOST_CHECK(!is_complex<double>());
  BOOST_CHECK(is_complex<typename std::complex<float>>());
  BOOST_CHECK(is_complex<typename std::complex<double>>());
  //
  BOOST_CHECK(num_character<float>() == 'F');
  BOOST_CHECK(num_character<double>() == 'D');
  BOOST_CHECK(num_character<typename std::complex<float>>() == 'C');
  BOOST_CHECK(num_character<typename std::complex<double>>() == 'Z');
  BOOST_CHECK(num_character<uint64_t>() == '?');
}
