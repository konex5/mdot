#ifndef INCLUDE_MDOT_HPP
#define INCLUDE_MDOT_HPP

#pragma once

#include <array>
#include <map>
#include <tuple>
#include <vector>

//#include <tbb/tbb.h>

namespace mdot {

// TODO: matrix triangular, hermitian
// TODO: matrix multiplication C format real+complex
// TODO: matrix svd C format

// tbb::flow::graph g;

inline constexpr std::array<double, 4> shid(double scale = 1.) {
  return {scale * 1., 0., 0., scale * 1.};
}

inline constexpr std::array<double, 4> shSz(double scale = 1.) {
  return {-scale * 0.5, 0., 0., scale * 0.5};
}

constexpr std::array<double, 4> test(double scale = 1.) {
  const auto id = shid(scale);
  const auto sz = shSz(scale);
  return {id[0] + sz[0], id[1] + sz[1], id[2] + sz[2], id[3] + sz[3]};
  // compile only?
}

// constexpr std::map<std::tuple<int,int>,std::vector<double> > shNoneSz() {
//     return
//     std::map<std::tuple<int,int>,std::vector<double>>({0,0},{1,2,3,4});
// }

// constexpr std::vector<std::tuple<int,int>,std::vector<double>> shNId(double
// scale = 1.) {
//     return {{0,0},std::vector<double>(shId(scale))};
// }

// constexpr std::vector<typename std::tuple< std::tuple<int,int>,
// std::vector<double> >> shId(char qbasis, double scale = 1.) {
//     if (qbasis == 'N') {
//         return {{std::tuple<int,int>({0,0}),{scale*1.,0.,0.,scale*1.}}};
//     } else if (qbasis == '1') {
//         //return {{{0,0},{scale*1.}},{{1,1},{scale*1.}}};
//     }
// }

inline int add(int a, int b) { return a + b; }

} // namespace mdot

#endif /* INCLUDE_MDOT_HPP */
