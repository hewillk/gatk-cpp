#pragma once

#include <array>
#include <cmath>

namespace hc {

struct QualityUtils {
  static constexpr char ASCII_OFFSET = '!';
  static double
  qual_to_error_prob(char qual) {
    return cache[qual];
  }

 private:
  static inline const auto cache = [] {
    std::array<double, 128> cache{};
    for (unsigned char i = ASCII_OFFSET; i < cache.size(); i++)
      cache[i] = std::pow(10, -(i - ASCII_OFFSET) / 10.0);
    return cache;
  }();
};

}  // namespace hc