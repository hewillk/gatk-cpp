#pragma once

#include <array>
#include <cmath>

namespace hc {

struct MathUtils {
  static double
  approximate_log10_sum_log10(double a, double b) {
    if (a > b) return approximate_log10_sum_log10(b, a);
    const auto diff = b - a;
    return b + (diff < JacobianLogTable::MAX_TOLERANCE ? JacobianLogTable::get(diff) : 0.0);
  }

 private:
  struct JacobianLogTable {
    static constexpr double MAX_TOLERANCE = 8.0;
    static double
    get(double difference) {
      return cache[std::round(difference * INV_STEP)];
    }

   private:
    static constexpr double TABLE_STEP = 0.0001;
    static constexpr double INV_STEP   = 1.0 / TABLE_STEP;
    static inline auto cache           = [] {
      std::array<double, static_cast<std::size_t>(MAX_TOLERANCE / TABLE_STEP) + 1> cache{};
      for (std::size_t k = 0; k < cache.size(); k++)
        cache[k] = std::log10(1.0 + std::pow(10.0, -TABLE_STEP * k));
      return cache;
    }();
  };
};

}  // namespace hc