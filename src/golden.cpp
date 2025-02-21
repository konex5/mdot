#include "golden.hpp"
#include <iostream>

int main(int argc, char *argv[]) {
  auto val = golden::add(4, 6);
  std::cout << "val=" << val << std::endl;
  return 0;
}