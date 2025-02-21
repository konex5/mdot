
#include <vector>
#include <map>
#include <string>

#include "mdot/include/babel_type.hpp"

dopbloc_t real_single_operator(std::string name, std::string qbasis) {
    std::map < std::string, typename std::map<std::string, dopbloc_t>> real_single_operators;
    std::map < std::string, typename std::map<std::string, double>> normalization;

    real_single_operators["sh-none"]["sh-id"][{0,0}] = {1,0,0,1}; 
    normalization["sh-none"]["sh-id"] = 1./sqrt(2);
    return real_single_operators[name][qbasis];
}

