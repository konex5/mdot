#include <iostream>
#include <string>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

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

    if (vm.count("help")) {
      std::cout << desc << std::endl;
      return 0;
    }
    if (vm.count("hamiltonian") == 0 || vm.count("output") == 0) {
      std::cout << desc << std::endl;
      std::cout << "cli_create_maximal_entangled_state: error: the following "
                   "arguments are required: -H/--hamiltonian, -o/--output"
                << std::endl;
      return 1;
    }
    
    boost::filesystem::path hamiltonian_path(
        vm["hamiltonian"].as<std::string>());
    boost::filesystem::path output_path(vm["output"].as<std::string>());

    if (!boost::filesystem::exists(hamiltonian_path)) {
      std::cout << "cli, create_maximal_entangled_state: the hamiltonian path "
                << hamiltonian_path << " is not a valid path." << std::endl;
      return 1;
    }

    auto output_dir = output_path.parent_path();
    if (!boost::filesystem::is_directory(output_dir)) {
      std::cout << "cli, create_maximal_entangled_state: the output dirpath "
                << output_dir << " is not a valid directory." << std::endl;
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
