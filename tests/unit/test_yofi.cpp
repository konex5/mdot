#include <boost/test/unit_test.hpp>

#include <boost/filesystem.hpp>

#include "mdot/utils/yofi.hpp"

using namespace fhmdot::yofi;

BOOST_AUTO_TEST_CASE(unittest_yofi_tuple) {

  boost::filesystem::path temp = boost::filesystem::temp_directory_path() /
                                 boost::filesystem::unique_path();
  const std::string tempstr = temp.native(); // optional
  // std::cout << "file path : " << tempstr << std::endl;

  {
    std::tuple<char, float, int> tout{'A', 1.5f, 6};
    std::tuple<char, float, int> tin;

    std::ofstream out;
    out.open(tempstr);
    write_tuple(out, tout); // write tuple to file
    out.close();
    std::ifstream in;
    in.open(tempstr);
    read_tuple(in, tin); // read the tuple from file
    std::cout << "  - Test1, tuple: (";
    print_tuple(tin); // print it out
    std::cout << ")" << std::endl;
    in.close();

    BOOST_CHECK(tout == tin);
    boost::filesystem::remove(temp);
  }
  {
    std::tuple<uint16_t, float, uint64_t, int> tout{5, 1.5f, 87, -5};
    std::tuple<uint16_t, float, uint64_t, int> tin;

    std::ofstream out(tempstr);
    write_tuple(out, tout); // write tuple to file
    out.flush();            // flush the file so we can read it immediately

    std::ifstream in(tempstr);
    read_tuple(in, tin); // read the tuple from file
    std::cout << "  - Test2, tuple: (";
    print_tuple(tin); // print it out
    std::cout << ")" << std::endl;

    BOOST_CHECK(tout == tin);
    boost::filesystem::remove(temp);
  }
}

#include <boost/mpl/list.hpp>
#include <complex>

typedef boost::mpl::list<float,
                         double> //, typename std::complex<float>> //,
                                 // std::complex<double>>
    values;

BOOST_AUTO_TEST_CASE_TEMPLATE(unittest_yofi_array, T, values) {

  boost::filesystem::path temp = boost::filesystem::temp_directory_path() /
                                 boost::filesystem::unique_path();
  const std::string tempstr = temp.native();
  // std::cout << "file path : " << tempstr << std::endl;

  {

    std::vector<T> dout{1, 3, 2, 4, 5};
    std::vector<T> din;

    std::ofstream out;
    out.open(tempstr);
    write_array(out, dout); // write tuple to file
    out.close();
    std::ifstream in(tempstr);
    read_array(in, din); // read the tuple from file
    std::cout << "  - Test*, data: (";
    print_array(din); // print it out
    std::cout << ")" << std::endl;
    in.close();

    BOOST_CHECK(dout == din);
    boost::filesystem::remove(temp);
  }
}
