
#include <sys/random.h>
#include <vector>

unsigned int SEED = 2748929928;

namespace utils {

inline std::array<double, 4> random_matrix() {
  return {double(rand_r(&SEED)) / double(RAND_MAX) - 0.5,
          double(rand_r(&SEED)) / double(RAND_MAX) - 0.5,
          double(rand_r(&SEED)) / double(RAND_MAX) - 0.5,
          double(rand_r(&SEED)) / double(RAND_MAX) - 0.5};
}

inline std::array<double, 4> normalize(std::array<double, 4> array) {
  std::array<double, 4> out;

  double norm = 0;
  for (std::size_t i = 0; i < array.size(); i++)
    norm += array[i] * array[i];

  for (std::size_t i = 0; i < array.size(); i++)
    out[i] = array[i] / norm;

  return out;
}

} // namespace utils