/// yet other file I/O

#include <algorithm>
#include <boost/type_traits/is_complex.hpp>
#include <fstream>
#include <iostream>
#include <iterator>
#include <tuple>
#include <type_traits>
#include <vector>

#include "mdot/include/babel_type.hpp"

namespace mdot {

namespace yofi {

// printing tuple
template <std::size_t I = 0, typename... Tp>
typename std::enable_if<I == sizeof...(Tp), void>::type
print_tuple(std::tuple<Tp...> &) {}

template <std::size_t I = 0, typename... Tp>
    typename std::enable_if <
    I<sizeof...(Tp), void>::type print_tuple(std::tuple<Tp...> &t) {
  std::cout << std::get<I>(t) << ", ";
  print_tuple<I + 1, Tp...>(t);
}

// writing tuple to a file
template <std::size_t I = 0, typename... Tp>
typename std::enable_if<I == sizeof...(Tp), void>::type
write_tuple(std::ofstream &, std::tuple<Tp...> &) {}

template <std::size_t I = 0, typename... Tp>
    typename std::enable_if <
    I<sizeof...(Tp), void>::type write_tuple(std::ofstream &ofs,
                                             std::tuple<Tp...> &t) {
  ofs.write((char *)&std::get<I>(t),
            sizeof(typename std::tuple_element<I, std::tuple<Tp...>>::type));
  write_tuple<I + 1, Tp...>(ofs, t);
}

// reading tuple from a file
template <std::size_t I = 0, typename... Tp>
typename std::enable_if<I == sizeof...(Tp), void>::type
read_tuple(std::ifstream &, std::tuple<Tp...> &) {}

template <std::size_t I = 0, typename... Tp>
    typename std::enable_if <
    I<sizeof...(Tp), void>::type read_tuple(std::ifstream &ifs,
                                            std::tuple<Tp...> &t) {
  ifs.read((char *)&std::get<I>(t),
           sizeof(typename std::tuple_element<I, std::tuple<Tp...>>::type));
  read_tuple<I + 1, Tp...>(ifs, t);
}

// printing vector
template <typename T> void print_array(std::vector<T> &arr) {
  for (const auto &item : arr) {
    std::cout << item << ",";
  }
}

// writing vector to a file
template <typename T>
void write_array(std::ofstream &ofs, std::vector<T> &arr) {
  std::copy(arr.begin(), arr.end(), std::ostreambuf_iterator<char>(ofs));
}

// reading vector to a file
template <typename T> void read_array(std::ifstream &ifs, std::vector<T> &arr) {
  auto start = std::istreambuf_iterator<char>(ifs);
  auto stop = std::istreambuf_iterator<char>();
  std::copy(start, stop, std::back_inserter(arr));
}

template <typename M> void print_mp(M &smp) {
  auto nb_blocks = smp.size();
  std::cout << "#number of blocks : " << nb_blocks << std::endl;
  std::cout << "#data type : " << nb_blocks << std::endl;

  for (auto &mp : smp.arr_) {
    //  mp.first
    //   write_tuple(ifs, index);
    //   write_tuple(ifs, shape);
    //   write_array(ifs, arr);
  }
}

template <typename M> void write_mp(std::ifstream &ifs, M &smp) {
  // for (index, shape, arr) {
  //   write_tuple(ifs, index);
  //   write_tuple(ifs, shape);
  //   write_array(ifs, arr);
  // }
}

template <typename M> void read_mp(std::ifstream &ifs, M &smp) {
  // for (index, shape, arr) {
  //   write_tuple(ifs, index);
  //   write_tuple(ifs, shape);
  //   write_array(ifs, arr);
  // }
}

} // namespace yofi

} // namespace mdot
