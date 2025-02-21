#include <boost/test/unit_test.hpp>

#include "mdot/include/babel_type.hpp"
#include "mdot/include/indexdw.hpp"

BOOST_AUTO_TEST_CASE(test_indexdw) {
  const std::vector<std::vector<dnum_t>> array_of_s = {
      {0.9, 0}, {0.8, 0.6, 0.01, 0}, {0.3, 0}};

  BOOST_CHECK(array_of_s[0][0] == 0.9);
  for (auto &it : array_of_s)
    BOOST_CHECK(it.at(it.size() - 1) == 0);

  dnum_t norm_out;
  norm_non_optimal(array_of_s, norm_out);
  BOOST_CHECK(norm_out == 1.3784411485442534);
}
