#pragma once

#include <stdint.h>

#include <array>
#include <complex>
#include <map>
#include <tuple>
#include <vector>

#ifndef FLOAT_PRECISION
using data_t = double;
#else
using data_t = float;
#endif

using dnum_t = data_t;
using znum_t = std::complex<data_t>;

using _darr_t = std::vector<dnum_t>;
using _zarr_t = std::vector<znum_t>;

using index_small_t = uint8_t;
using index_t = uint16_t;

using op_index_t = std::tuple<index_small_t, index_small_t>;
using op_shape_t = std::tuple<op_index_t, op_index_t>;
using dopbloc_t = std::map<op_index_t, _darr_t>;
using zopbloc_t = std::map<op_index_t, _zarr_t>;

using m_index_t = std::tuple<index_t, index_small_t, index_t>;
using m_shape_t = std::tuple<index_t, index_small_t, index_t>;
using dmbloc_t = std::map<m_index_t, std::pair<m_shape_t, _darr_t>>;
using zmbloc_t = std::map<m_index_t, std::pair<m_shape_t, _zarr_t>>;

using t_index_t = std::tuple<index_t, index_small_t, index_small_t, index_t>;
using t_shape_t = std::tuple<index_t, index_small_t, index_small_t, index_t>;
using dtbloc_t = std::map<t_index_t, std::pair<t_shape_t, _darr_t>>;
using ztbloc_t = std::map<t_index_t, std::pair<t_shape_t, _zarr_t>>;

using g_index_t =
    std::tuple<index_small_t, index_small_t, index_small_t, index_small_t>;
using g_shape_t =
    std::tuple<index_small_t, index_small_t, index_small_t, index_small_t>;
using dgbloc_t = std::map<g_index_t, std::pair<g_shape_t, _darr_t>>;
using zgbloc_t = std::map<g_index_t, std::pair<g_shape_t, _zarr_t>>;

template <typename, typename = void> struct real_type;

template <typename T>
struct real_type<T, std::enable_if_t<std::is_arithmetic_v<T>>> {
  using type = T;
};

template <typename T> struct real_type<std::complex<T>, void> {
  using type = T;
};

//
using sh_array_t = std::array<dnum_t, 4>;
