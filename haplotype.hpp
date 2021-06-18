#pragma once

#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>
#include <cmath>
#include <map>
#include <string>

#include "cigar.hpp"
#include "interval.hpp"
#include "variant.hpp"

namespace hc {

struct Haplotype {
  std::string bases;
  Interval location;
  std::map<int, Variant> event_map;
  Cigar cigar;
  int alignment_begin_wrt_ref = 0;
  double score                = std::numeric_limits<double>::lowest();
  int rank;

  Haplotype() = default;
  Haplotype(std::string bases, double score) : bases(std::move(bases)), score(score) {}

  std::size_t
  size() const noexcept {
    return bases.size();
  }

  auto
  get_overlapping_events(std::size_t begin) const {
    std::vector<Variant> events;
    auto upper_bound = event_map.upper_bound(begin);
    for (auto it = event_map.begin(); it != upper_bound; ++it)
      if (it->second.location.end > begin) events.push_back(it->second);
    return events;
  }

  template<typename Archive>
  void
  serialize(Archive& ar, const unsigned int version) {
    ar& bases;
    ar& location;
    ar& event_map;
    ar& cigar;
    ar& alignment_begin_wrt_ref;
    ar& score;
  }
};

}  // namespace hc