#ifndef INCLUDE_MDOT_HPP
#define INCLUDE_MDOT_HPP

#pragma once

#include <array>
#include <tuple>
#include <map>
#include <vector>

namespace mdot {

inline constexpr std::array<double,4> shId(double mult = 1.) {
    return {mult*1.,0.,0.,mult*1.};
}

inline constexpr std::array<double,4> shSz(double mult = 1.) {
    return {-mult*0.5,0.,0.,mult*0.5};
}

constexpr std::array<double,4> test(double mult = 1.) {
    const auto id = shId(mult);
    const auto sz = shSz(mult);
    return {id[0]+sz[0],id[1]+sz[1],id[2]+sz[2],id[3]+sz[3]};
    // compile only?
}

// constexpr std::map<std::tuple<int,int>,std::vector<double> > shNoneSz() {
//     return std::map<std::tuple<int,int>,std::vector<double>>({0,0},{1,2,3,4});
// }


// constexpr std::vector<std::tuple<int,int>,std::vector<double>> shNId(double mult = 1.) {
//     return {{0,0},std::vector<double>(shId(mult))};
// }

// constexpr std::vector<typename std::tuple< std::tuple<int,int>, std::vector<double> >> shId(char qbasis, double mult = 1.) {
//     if (qbasis == 'N') {
//         return {{std::tuple<int,int>({0,0}),{mult*1.,0.,0.,mult*1.}}};
//     } else if (qbasis == '1') {
//         //return {{{0,0},{mult*1.}},{{1,1},{mult*1.}}};
//     }
// }

inline int add(int a, int b) { return a + b; }
    std::map<int,char> c;
} // namespace mdot

#endif /* INCLUDE_MDOT_HPP */
