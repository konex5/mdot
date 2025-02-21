#include <boost/program_options.hpp>

#include <iostream>
#include <string>

namespace po = boost::program_options;

//#include "golden/include/gold.hpp"

int main(int ac, char *av[]) {

  try {

    po::options_description desc("Allowed options");
    desc.add_options()("help,h", "produce help message")(
        "name,n", po::value<std::string>(), "set a name");

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
      std::cout << desc << "\n";
      return 0;
    }

    if (vm.count("name")) {
      std::cout << "Welcome to you, " << vm["name"].as<std::string>()
                << ", this is a cpp project!" << std::endl;
    } else {
      std::cout << "Welcome to nokx cpp golden project!" << std::endl;
    }

  } catch (std::exception &e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  } catch (...) {
    std::cerr << "Exception of unknown type!\n";
  }

  return 0;
}
