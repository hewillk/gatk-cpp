#pragma once

#include "sam.hpp"

namespace hc {

struct MappingQualityReadFilter {
  static constexpr std::uint16_t MIN_MAPPING_QUALITY_SCORE = 20;
  bool
  operator()(const SAMRecord& record) const {
    return record.MAPQ < MIN_MAPPING_QUALITY_SCORE;
  }
};

struct DuplicateReadFilter {
  bool
  operator()(const SAMRecord& record) const {
    return record.DUPLICATE_READ();
  }
};

struct SecondaryAlignmentReadFilter {
  bool
  operator()(const SAMRecord& record) const {
    return record.SECONDARY_ALIGNMENT();
  }
};

struct MinimumLengthReadFilter {
  static constexpr std::size_t MINIMUM_READ_LENGTH_AFTER_TRIMMING = 10;
  bool
  operator()(const SAMRecord& record) const {
    return record.size() < MINIMUM_READ_LENGTH_AFTER_TRIMMING;
  }
};

struct MateOnSameContigReadFilter {
  bool
  operator()(const SAMRecord& record) const {
    return record.RNEXT != "=";
  }
};

}  // namespace hc