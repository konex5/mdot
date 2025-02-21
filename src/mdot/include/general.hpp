#pragma once

#include <map>
#include <string>

#include <boost/filesystem.hpp>

#include "mdot/include/babel_type.hpp"
#include "mdot/include/blocs_static.hpp"
#include "mdot/include/operators_static.hpp"

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
  model_dict["qn_name"] = 1;   // u1
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

std::tuple<std::vector<dnum_t>, std::vector<dtbloc_t>>
create_maximal_entangled_state(std::size_t size, std::size_t spin_name,
                               std::size_t qn_name) {
  std::vector<dnum_t> coef;
  std::vector<dtbloc_t> dmps;
  t_index_t tmp_indices;
  darr_t tmp_vec;
  for (std::size_t i = 0; i < size; i++) {
    coef.push_back(1. / sqrt(2));
    if (spin_name == 0 && qn_name == 0) {
      tmp_indices = {0, 0, 0, 0};
      auto arr =
          mdot::real_sh_operators_crtp<mdot::sh_id_no>::times(1. / sqrt(2));
      tmp_vec = std::vector<dnum_t>(arr.begin(), arr.end());
      // dmps.push_back({tmp_indices, tmp_vec});
    } else if (spin_name == 0 && qn_name == 1) {
      auto tmp_blocs = mdot::real_sh_blocs_crtp<mdot::sh_id_u1>::get_blocs();
      dtbloc_t dst_blocs;
      for (auto &[key, val] : tmp_blocs) {
        // dst_blocs[{0,std::get<0>(key),std::get<1>(key),0}] =
        // {std::get<0>(val),{1,std::get<0>(std::get<1>(val)),std::get<1>(std::get<1>(val)),1},std::get<2>(val)};
      }
      // dmps.push_back(dst_blocs)
    }
  }
  return {coef, dmps};
}