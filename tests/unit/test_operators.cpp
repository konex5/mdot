#include <boost/test/unit_test.hpp>

#include "mdot/include/babel_type.hpp"
#include "mdot/include/operators.hpp"

BOOST_AUTO_TEST_CASE(test_indexdw)
{

    auto a = real_single_operator("hello", "world");
    
    BOOST_CHECK((a[{0,0}])[0] == 1);
}
