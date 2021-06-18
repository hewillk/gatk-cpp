#pragma once

#include <boost/serialization/string.hpp>
#include <stdexcept>
#include <string>

namespace hc {

struct Interval {
  static constexpr char CONTIG_SEPARATOR    = ':';
  static constexpr char BEGIN_END_SEPARATOR = '-';
  static constexpr char END_OF_CONTIG       = '+';
  static constexpr char DIGIT_SEPARATOR     = ',';

  static bool
  is_valid(const Interval& interval) noexcept {
    return interval.end >= interval.begin;
  }

  std::string contig;
  int begin  = 0;
  int end    = 0;

  Interval() = default;

  Interval(std::string contig, int begin, int end)
  : contig(std::move(contig)), begin(begin), end(end) {
    if (!is_valid(*this))
      throw std::invalid_argument("Interval::Interval(std::string, std::size_t, std::size_t)");
  }

  Interval(const char* string) {
    std::string str = string;
    auto colon      = str.find(CONTIG_SEPARATOR);
    if (colon == std::string::npos) {
      contig = str;
      begin  = 0;
      end    = std::numeric_limits<int>::max();
    } else {
      contig      = str.substr(0, colon);
      auto remain = str.substr(colon + 1);
      remain.erase(std::remove(remain.begin(), remain.end(), DIGIT_SEPARATOR), remain.end());
      begin     = std::stoul(remain);
      auto dash = remain.find(BEGIN_END_SEPARATOR);
      if (dash == std::string::npos) {
        if (remain.back() == END_OF_CONTIG)
          end = std::numeric_limits<int>::max();
        else
          end = begin + 1;
      } else
        end = std::stoul(remain.substr(dash + 1));
    }
    if (!is_valid(*this)) throw std::invalid_argument("Interval::Interval(const char*)");
  }

  std::size_t
  size() const noexcept {
    return end - begin;
  }

  bool
  empty() const noexcept {
    return size() == 0;
  }

  bool
  overlaps(const Interval& other) const noexcept {
    return contig == other.contig && begin < other.end && other.begin < end;
  }

  bool
  contains(const Interval& other) const noexcept {
    return contig == other.contig && begin <= other.begin && end >= other.end;
  }

  Interval
  span_with(const Interval& other) {
    if (contig != other.contig)
      throw std::invalid_argument(
        "Interval::span_with(): Cannot get span for intervals on different contigs.");
    return {contig, std::min(begin, other.begin), std::max(end, other.end)};
  }

  Interval
  expand_within_contig(int padding) {
    return {contig, begin - padding, end + padding};
  }

  std::string
  to_string() const {
    return contig + CONTIG_SEPARATOR + std::to_string(begin < 0 ? 0 : begin) + BEGIN_END_SEPARATOR +
           std::to_string(end);
  }

  friend bool
  operator<(const Interval& lhs, const Interval& rhs) {
    return std::tie(lhs.contig, lhs.begin, lhs.end) < std::tie(rhs.contig, rhs.begin, rhs.end);
  }

  friend bool
  operator==(const Interval& lhs, const Interval& rhs) {
    return std::tie(lhs.contig, lhs.begin, lhs.end) == std::tie(rhs.contig, rhs.begin, rhs.end);
  }

  template<typename Archive>
  void
  serialize(Archive& ar, const unsigned int version) {
    ar& contig;
    ar& begin;
    ar& end;
  }
};

}  // namespace hc
