#pragma once

#include <map>
#include <string>

#include <boost/filesystem.hpp>

inline static const std::string spin_name(const std::size_t num) {
  // this functionality shouldn't be used in practice.
  if (num == 0)
    return "sh";
  else if (num == 1)
    return "so";
  else if (num == 2)
    return "sf";
  else if (num == 3)
    return "ldsh";
  else
    return "";
}

inline static const std::string qn_name(const std::size_t num) {
  // this functionality shouldn't be used in practice.
  if (num == 0)
    return "no";
  else if (num == 1)
    return "u1";
  else
    return "";
}

inline std::string get_model_name(boost::filesystem::path hamiltonian_path) {
  return "description";
}

inline std::map<std::string, std::size_t>
get_model_info(boost::filesystem::path hamiltonian_path) {
  std::map<std::string, std::size_t> model_dict;
  model_dict["size"] = 10;
  model_dict["interaction_range"] = 2;
  model_dict["spin_name"] = 0; // sh
  model_dict["qn_name"] = 1; // u1
  return model_dict;
}

inline std::map<std::string, double>
get_model_parameters(boost::filesystem::path hamiltonian_path) {
  std::map<std::string, double> model_dict;
  model_dict["J"] = 1.0;
  model_dict["Jz"] = 1.5;
  return model_dict;
}

// get_details_zdmrg

// get_details_tdmrg