#include <iostream>

#include "mdot/include/operators_static.hpp"

int main() {
  // auto arr = mdot::real_sh_operators_crtp<mdot::sh_id_no>::array;
  // auto shape = mdot::real_sh_operators_crtp<mdot::sh_id_no>::shape;
  // auto size = mdot::real_sh_operators_crtp<mdot::sh_id_no>::size;

  auto arr = mdot::cplx_sh_operators_crtp<mdot::sh_id_cplx_no>::array;
  auto shape = mdot::cplx_sh_operators_crtp<mdot::sh_id_cplx_no>::shape;
  auto size = mdot::cplx_sh_operators_crtp<mdot::sh_id_cplx_no>::size;

  std::cout << "the matrix is..." << std::endl;
  size_t nb_col = std::get<1>(shape);
  for (size_t i = 0; i < std::get<0>(shape); i++) {
    for (size_t j = 0; j < nb_col; j++)
      std::cout << arr[j + i * nb_col] << " ";
    std::cout << std::endl;
  }
  std::cout << "yeah!" << std::endl;

  return 0;
}
