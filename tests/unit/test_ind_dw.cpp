#include <boost/test/unit_test.hpp>

#include "mdot/include/typedef.hpp"

BOOST_AUTO_TEST_CASE(test_ind_dw) {
  array_of_s_type array_of_s = {{0.9,0},{0.8,0.6,0.01,0},{0.3,0}};
  
  BOOST_CHECK(array_of_s[0][0] == 0.9);
  for (auto& it: array_of_s)
    BOOST_CHECK(it.at(it.size()-1) == 0);



}
