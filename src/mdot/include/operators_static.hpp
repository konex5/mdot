#pragma once

#include "mdot/include/babel_type.hpp"
#include <cmath>
#include <array>
#include <complex>

// todo CRTP name, qname
/*
::normalization
::indices_blocs_pair

::nsize
::nshape
::array
*/

namespace mdot
{

  template <class Q>
  struct real_operators_crtp
  {
    static constexpr index_t size = Q::n_size;
    static constexpr std::array<index_t, 2> shape = Q::n_shape;
    static constexpr std::array<dnum_t, 4> array = Q::n_array;
    inline static constexpr std::array<dnum_t, 4> times(const dnum_t alpha) { return {alpha * Q::n_array[0], alpha * Q::n_array[1], alpha * Q::n_array[2], alpha * Q::n_array[3]}; }
    inline static constexpr std::array<znum_t, 4> times(const znum_t alpha) { return {alpha * Q::n_array[0], alpha * Q::n_array[1], alpha * Q::n_array[2], alpha * Q::n_array[3]}; }
    inline static constexpr dnum_t trace() { return Q::n_array[0] + Q::n_array[3]; };
    inline static constexpr std::array<dnum_t, 4> square() { return {
        Q::n_array[0] * Q::n_array[0] + Q::n_array[1] * Q::n_array[2],
        Q::n_array[0] * Q::n_array[1] + Q::n_array[1] * Q::n_array[3],
        Q::n_array[2] * Q::n_array[0] + Q::n_array[3] * Q::n_array[2],
        Q::n_array[2] * Q::n_array[1] + Q::n_array[3] * Q::n_array[3]}; };
  };

  template <class Q>
  struct cplx_operators_crtp
  {
    static constexpr index_t size = Q::n_size;
    static constexpr std::array<index_t, 2> shape = Q::n_shape;
    static constexpr std::array<znum_t, 4> array = Q::n_array;
    
    inline static constexpr std::array<znum_t, 4> times(const dnum_t alpha) { return {alpha * Q::n_array[0], alpha * Q::n_array[1], alpha * Q::n_array[2], alpha * Q::n_array[3]}; }
    inline static constexpr std::array<znum_t, 4> times(const znum_t alpha) { return {alpha * Q::n_array[0], alpha * Q::n_array[1], alpha * Q::n_array[2], alpha * Q::n_array[3]}; }
    inline static constexpr znum_t trace() { return Q::n_array[0] + Q::n_array[3]; };
    inline static constexpr std::array<znum_t, 4> square() { return {
        Q::n_array[0] * Q::n_array[0] + Q::n_array[1] * Q::n_array[2],
        Q::n_array[0] * Q::n_array[1] + Q::n_array[1] * Q::n_array[3],
        Q::n_array[2] * Q::n_array[0] + Q::n_array[3] * Q::n_array[2],
        Q::n_array[2] * Q::n_array[1] + Q::n_array[3] * Q::n_array[3]}; };
  };

  struct sh_id_no : real_operators_crtp<sh_id_no>
  {
    static constexpr index_t n_size = 4;
    static constexpr std::array<index_t, 2> n_shape = {2, 2};
    static constexpr std::array<dnum_t, 4> n_array = {1, 0, 0, 1};
  };
  struct sh_id_no : cplx_operators_crtp<sh_id_no>
  {
    static constexpr index_t n_size = 4;
    static constexpr std::array<index_t, 2> n_shape = {2, 2};
    static constexpr std::array<znum_t, 4> n_array = {1 , 0, 0, 1};
  };

  struct sh_sp_no : real_operators_crtp<sh_sp_no>
  {
    static constexpr index_t n_size = 4;
    static constexpr std::array<index_t, 2> n_shape = {2, 2};
    static constexpr std::array<dnum_t, 4> n_array = {0, 1, 0, 0};
  };

    struct sh_sm_no : real_operators_crtp<sh_sm_no>
  {
    static constexpr index_t n_size = 4;
    static constexpr std::array<index_t, 2> n_shape = {2, 2};
    static constexpr std::array<dnum_t, 4> n_array = {0, 0, 1, 0};
  };

  struct sh_sx_no : real_operators_crtp<sh_sx_no>
  {
    static constexpr index_t n_size = 4;
    static constexpr std::array<index_t, 2> n_shape = {2, 2};
    static constexpr std::array<dnum_t, 4> n_array = {0, 1, 1, 0};
  };
  
  struct sh_sy_no : cplx_operators_crtp<sh_sy_no>
  {
    static constexpr index_t n_size = 4;
    static constexpr std::array<index_t, 2> n_shape = {2, 2};
    static constexpr std::array<znum_t, 4> n_array = {0 , 0+ 1.i, 0- 1.i, 0};
  };

  struct sh_sz_no : real_operators_crtp<sh_sz_no>
  {
    static constexpr index_t n_size = 4;
    static constexpr std::array<index_t, 2> n_shape = {2, 2};
    static constexpr std::array<dnum_t, 4> n_array = {1, 0, 0, -1};
  };
  
  /*


template <class Q> struct blocs_crtp {
  //  inline static indices_blocs_t get() {
  //    std::map<indices,operators> a;
  //    a[{0,0}] = operators<name>::
  //    no idea if map can organize statically
  //
  // return Q::qsum(q1, q2); }

  // static constexpr index_t n_array_pair = Q::n_array_pair;
  // static constexpr index_t n_pairs = Q::n_pairs;
  // static constexpr index_t operators = Q::n_operators;

  static constexpr index_t dim = Q::qn_dim;
  inline static index_t sum(index_t q1, index_t q2) { return Q::qsum(q1, q2); }
  inline static index_t sub(index_t q1, index_t q2) { return Q::qsub(q1, q2); }
};

struct sh_none : quantum_number_crtp<sh_none> {
  static constexpr index_t qn_dim = 1;
  inline static index_t qsum(index_t q1, index_t q2) {
    (void)q1, (void)q2;
    return 0;
  }
  inline static index_t qsub(index_t q1, index_t q2) {
    (void)q1, (void)q2;
    return 0;
  }
};

struct sh_u1 : quantum_number_crtp<sh_u1> {
  static constexpr index_t qn_dim = 2;
  inline static index_t qsum(index_t q1, index_t q2) { return q1 + q2; }
  inline static index_t qsub(index_t q1, index_t q2) { return q1 - q2; }
};

struct sh_su2 : quantum_number_crtp<sh_su2> {
  static constexpr index_t qn_dim = 1;
  inline static index_t qsum(index_t q1, index_t q2) { return q1 + q2; }
  inline static index_t qsub(index_t q1, index_t q2) { return q1 - q2; }
};







struct sh_u1 : quantum_number_crtp<sh_u1> {
  static constexpr index_t qn_dim = 2;
  inline static index_t qsum(index_t q1, index_t q2) { return q1 + q2; }
  inline static index_t qsub(index_t q1, index_t q2) { return q1 - q2; }
};

struct sh_su2 : quantum_number_crtp<sh_su2> {
  static constexpr index_t qn_dim = 1;
  inline static index_t qsum(index_t q1, index_t q2) { return q1 + q2; }
  inline static index_t qsub(index_t q1, index_t q2) { return q1 - q2; }
};

*/

} // namespace mdot
