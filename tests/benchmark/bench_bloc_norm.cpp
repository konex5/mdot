#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <iostream>
#include <string>

#include <sys/resource.h>
#include <sys/time.h>
#include <tbb/tbb.h>

//#include "mdot/include/routine/svd_routine.hpp"

namespace po = boost::program_options;

using namespace boost::filesystem;

void format_csv(void) {
  // std::ofstream fo;
  // std::ifstream fi;
  // fo.open("beautify.csv", std::ofstream::out);
  // fo << "tokens"
  //    << ", "
  //    << "gpus"
  //    << ", "
  //    << "tiling"
  //    << ", "
  //    << "quantize"
  //    << ", "
  //    << "encode time [s]"
  //    << ", "
  //    << "decode time [s]"
  //    << ", "
  //    << "memory encode [MB]"
  //    << ", "
  //    << "memory decode [MB]"
  //    << "\n";
  // /*    try {
  //         fi.open("./results.csv", std::ifstream::in);
  //         std::string line_encoder, line_decoder;
  //         std::getline(fi, line_encoder); // Skip the header

  //         while (std::getline(fi, line_encoder)) {
  //             std::getline(fi, line_decoder);
  //             std::stringstream ss_encoder(line_encoder),
  //    ss_decoder(line_decoder); auto getItems = [](std::stringstream ss) {
  //                 std::string val;
  //                 std::vector<std::string> codec;
  //                 codec.reserve(7);
  //                 while (ss >> val) {
  //                     val = val.back() == ',' ? val.substr(0, val.size() - 1) :
  //    val; codec.push_back(val);
  //                 }
  //                 return codec;
  //             };
  //             auto encoder(getItems(std::move(ss_encoder)));
  //             auto decoder(getItems(std::move(ss_decoder)));

  //             for (auto i = 0; i < 4; i++) {
  //                 if (encoder.at(i) != decoder.at(i)) {
  //                     std::cout << "Error when reading results.csv  " <<
  //    encoder.at(i) << "  " << decoder.at(i)
  //                               << std::endl;
  //                     fo.close();
  //                     fi.close();
  //                     return;
  //                 }
  //             }
  //             fo << encoder.at(0) << ", " << encoder.at(1) << ", " <<
  //    encoder.at(2) << ", " << encoder.at(3) << ", "
  //                << encoder.at(5) << ", " << decoder.at(5) << ", " <<
  //    encoder.at(6) << ", " << decoder.at(6) << ", "
  //                << std::endl;
  //         }
  //         if (boost::filesystem::remove(boost::filesystem::path("results.csv")))
  //    { boost::filesystem::rename(boost::filesystem::path("beautify.csv"),
  //    boost::filesystem::path("results.csv")); } else { std::cout << "Could not
  //    beautify. Could not rename" << std::endl;
  //         }
  //     } catch (const std::exception &e) {
  //         std::cout << "Could not beautify" << std::endl;
  //         fo.close();
  //         fi.close();
  //     }
  //     */
  // fo.close();
  // fi.close();
}

int main(int argc, const char *argv[]) {
  // get number of core HT on
  const int num_cores = tbb::task_scheduler_init::default_num_threads();
  // ilc::timer t("elapse time main");

  po::variables_map vm;
  po::options_description desc("Allowed options");
  // clang-format off
    
    desc.add_options()("help,h", "produce help message")
        ("level,l", po::value<std::string>()->default_value("4"),"level of compression")
        ("input,i", po::value<std::string>()->default_value("/tmp"),"Input directory data")
        ("ngpu,n", po::value<int>()->default_value(0),"number of GPUs")
        ("tokens,t", po::value<int>()->default_value(1),"number of tokens")
        ("add-header", po::bool_switch()->default_value(false),"Add CSV header only")
        ("beautify", po::bool_switch()->default_value(false),"Converts the results to a different table format")
        ;
  // clang-format on
  po::store(po::parse_command_line(argc, argv, desc), vm);

  try {
    po::notify(vm);
  } catch (std::exception &e) {
    std::cerr << "Error: " << e.what() << "\n";
  }

  if (vm.count("help")) {
    std::cout << desc;
    exit(0);
  }

  if (vm["tokens"].as<int>() > (num_cores / 2)) {
    std::cout << "Number of tokens " << vm["tokens"].as<int>()
              << " not supported" << std::endl;
    exit(0);
  }

  std::string in = vm["input"].as<std::string>();
  const int ngpu = vm["ngpu"].as<int>();
  const int ntokens = vm["tokens"].as<int>();
  const bool addheaderonly = vm["add-header"].as<bool>();
  const bool beautify = vm["beautify"].as<bool>();

  /*
  std::string model = "generic_" + vm["level"].as<std::string>();

  ilc::codec c({ilc::helper_build_path::root_path() + (bquantize ?
  "models_quantize/" : "models/")}); c.model(model); c.gpu(ngpu);
  c.tiling(btiling);
  c.set_png_compression_lvl(lvl_zlib);
  */
  // // create and open the .csv file
  // std::ofstream fs;
  // // create a name for the file output
  // std::string filename = "results.csv";
  // fs.open(filename, std::ofstream::out | std::ofstream::app);
  // if (addheaderonly) {
  //   fs << "tokens"
  //      << ", "
  //      << "gpus"
  //      << ", "
  //      << "time [s]"
  //      << ", "
  //      << "mem used [MB],"
  //      << "\n";
  //   fs.close();
  //   exit(0);
  // }
  /*
  if (beautify) {
      format_csv();
      exit(0);
  }

  path working_dir = path(working_directory);

  c.ntokens(ntokens + ngpu);
  ilc::timer tcodec("codec");
  auto codec_action = [&in, &working_dir, &codec, &c]() {
      codec == codec_type::encoder ? c.encode(in, working_dir.string()) :
  c.decode(in, working_dir.string());
  };
  codec_action();
  float res_codec = tcodec.get();

  auto meminmb = []() {
      struct rusage thisusage;
      getrusage(RUSAGE_SELF, &thisusage);
      return thisusage.ru_maxrss >> 10;
  };

  auto codec_string = [&codec]() { return codec == codec_type::encoder ?
  "encoder" : "decoder"; };

  auto float_to_string_with_precision = [](const float number, const int p) {
      std::stringstream ss;
      ss << std::fixed << std::setprecision(p) << number;
      return ss.str();
  };
  // The reading here ru_maxrss is in line with the output from /usr/bin/time -v
  command int memory_used = meminmb();

  fs << std::to_string(ntokens + ngpu) << ", " << std::to_string(ngpu) << ", "
  << (btiling ? "true" : "false") << ", "
     << (bquantize ? "true" : "false") << ", " << codec_string() << ", "
     << float_to_string_with_precision(res_codec, 2) << ", " <<
  std::to_string(memory_used) << ","
     << "\n";

  fs.close();
  */
  return 0;
}
