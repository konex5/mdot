#include "mdot/include/babel_type.hpp"
#include "mdot/include/routine/svd_routine.hpp"
#include <boost/test/unit_test.hpp>
#include <tbb/tbb.h>

BOOST_AUTO_TEST_CASE(test_routine_svd_routine_bloc_norm) {
  {
    std::vector<std::vector<dnum_t>> a = {{1, 2, 3}, {2, 1}, {1, 1, 2}};
    double norm_out;
    mdot::bloc_norm(a, {}, norm_out);
    BOOST_CHECK(norm_out == 5.);
    mdot::normalize_the_array(a, {});
    BOOST_CHECK(a[0][0] == 0.2);
  }
  //
  {
    const std::vector<std::vector<dnum_t>> a = {{3, 2, 1}, {2, 1}, {2, 1, 1}};
    double norm_out;
    mdot::bloc_norm(a, {3, 2, 3}, norm_out);
    BOOST_CHECK(norm_out == 5);

    {
      dnum_t dw = 0;
      auto cut_at_index = mdot::truncation_strategy(a, 10, dw, 0.1);
      std::vector<index_t> cut_at_index_result = {3, 2, 3};
      for (size_t i = 0; i < cut_at_index.size(); i++)
        BOOST_CHECK(cut_at_index[i] == cut_at_index_result[i]);
      BOOST_CHECK(dw == 2);
    }
    {
      dnum_t dw = 0;
      auto cut_at_index = mdot::truncation_strategy(a, 10, dw, 0.2);
      std::vector<index_t> cut_at_index_result = {2, 1, 1};
      for (size_t i = 0; i < cut_at_index.size(); i++)
        BOOST_CHECK(cut_at_index[i] == cut_at_index_result[i]);
      BOOST_CHECK(dw == 4);
    }
    {
      dnum_t dw = 0;
      auto cut_at_index = mdot::truncation_strategy(a, 10, dw, 0.8);
      std::vector<index_t> cut_at_index_result = {1, 0, 0};
      for (size_t i = 0; i < cut_at_index.size(); i++)
        BOOST_CHECK(cut_at_index[i] == cut_at_index_result[i]);
      BOOST_CHECK(dw == 16);
    }
  }
}