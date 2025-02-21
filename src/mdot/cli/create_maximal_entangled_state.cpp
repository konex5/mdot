#include <boost/program_options.hpp>

#include <iostream>
#include <string>

namespace po = boost::program_options;

//#include "mdot/include/mdot.hpp"

int main(int argc, char *argv[]) {

  try {

    po::options_description desc("cli, create_maximal_entangled_state");
    desc.add_options()("help,h", "produce help message")(
        "hamiltonian,H", po::value<std::string>(), "hamiltonian path")(
        "output,o", po::value<std::string>(), "output path");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help") || vm.count("hamiltonian") || vm.count("output")) {
      std::cout << desc << "\n";
      return 0;
    }

    // boost::pathname
    if (vm.count("hamiltonian")) {
      // hamiltonian_file = vm["hamiltonian"].as<std::string>();
    } else {
      std::cout << "cli, create_maximal_entangled_state: the hamiltonian path "
                << vm["hamiltonian"].as<std::string>() << " is not a valid path"
                << std::endl;
      return 1;
    }

  } catch (std::exception &e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  } catch (...) {
    std::cerr << "Exception of unknown type!\n";
  }

  return 0;
}
