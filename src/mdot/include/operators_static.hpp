#pragma once

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

typedef std::size_t index_t;

namespace mdot {


template <class Q> struct real_operators_crtp {
  // static constexpr double normalization = Q::n_normalization; 
  static constexpr index_t size = Q::n_size;
  static constexpr std::array<index_t,2> shape = Q::n_shape;
  static constexpr std::array<double,4> array = Q::n_array;

  //inline static index_t sub(index_t q1, index_t q2) { return Q::qsub(q1, q2); }
};


struct sh_id_no : real_operators_crtp<sh_id_no> {
  // static constexpr double n_normalization = 1./sqrt(2); 
  static constexpr index_t n_size = 4;
  static constexpr std::array<index_t,2> n_shape = {2,2};
  static constexpr std::array<double,4> n_array = {1,0,0,1};

  // inline static index_t qsub(index_t q1, index_t q2) {
  //   (void)q1, (void)q2;
  //   return 0;
  // }
};


struct sh_id_cplx_no : real_operators_crtp<sh_id_cplx_no > {
  // static constexpr double n_normalization = 1./sqrt(2); 
  static constexpr index_t n_size = 4;
  static constexpr std::array<index_t,2> n_shape = {2,2};
  //static constexpr std::array<std::complex<double> ,4> n_array = {{1,0},{0,0},{0,0},{1,0}};

  // inline static index_t qsub(index_t q1, index_t q2) {
  //   (void)q1, (void)q2;
  //   return 0;
  // }
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
