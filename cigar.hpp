#pragma once

#include <boost/serialization/vector.hpp>
#include <istream>
#include <ostream>
#include <string>
#include <vector>

namespace hc {

enum class CigarOperator : char {
  /** Match or mismatch */
  M  = 'M',
  /** Insertion vs. the reference. */
  I  = 'I',
  /** Deletion vs. the reference. */
  D  = 'D',
  /** Skipped region from the reference. */
  N  = 'N',
  /** Soft clip. */
  S  = 'S',
  /** Hard clip. */
  H  = 'H',
  /** Padding. */
  P  = 'P',
  /** Matches the reference. */
  EQ = '=',
  /** Mismatches the reference. */
  X  = 'X'
};

struct CigarElement {
  int length{};
  CigarOperator op{};

  constexpr CigarElement() = default;
  constexpr CigarElement(int length, CigarOperator op) : length(length), op(op) {}

  template<typename Archive>
  void
  serialize(Archive& ar, const unsigned int version) {
    ar& length;
    ar& op;
  }

  auto
  to_string() const {
    return std::to_string(length) + char(op);
  }
};

struct Cigar {
 private:
  static auto
  to_cigar_elements(const std::string& cigar_string) {
    std::vector<CigarElement> cigar_elements;
    for (std::size_t i = 0; i < cigar_string.size(); i++) {
      auto length = static_cast<std::uint32_t>(cigar_string[i] - '0');
      for (i++; std::isdigit(cigar_string[i]); i++) length = length * 10 + cigar_string[i] - '0';
      cigar_elements.emplace_back(length, CigarOperator(cigar_string[i]));
    }
    return cigar_elements;
  }

 public:
  Cigar()             = default;
  Cigar(const Cigar&) = default;
  Cigar(Cigar&&)      = default;
  Cigar&
  operator=(const Cigar&) = default;
  Cigar&
  operator=(Cigar&&) = default;
  Cigar(const std::string& cigar_string) : cigar_elements(to_cigar_elements(cigar_string)) {}
  Cigar(size_t size, CigarElement element) : cigar_elements(size, element) {}
  Cigar&
  operator=(const std::string& cigar_string) {
    cigar_elements = to_cigar_elements(cigar_string);
    return *this;
  }

  void
  push_back(CigarElement element) {
    cigar_elements.push_back(element);
  }

  void
  emplace_back(std::size_t length, CigarOperator op) {
    cigar_elements.emplace_back(length, op);
  }

  auto
  get_reference_length() const {
    std::size_t reference_length = 0;
    for (auto [length, op] : cigar_elements) {
      switch (op) {
        case CigarOperator::M:
        case CigarOperator::D:
        case CigarOperator::N:
        case CigarOperator::EQ:
        case CigarOperator::X:
          reference_length += length;
          break;
        default:
          break;
      }
    }
    return reference_length;
  }

  auto
  get_read_length() const {
    std::size_t read_length = 0;
    for (auto [length, op] : cigar_elements) {
      switch (op) {
        case CigarOperator::M:
        case CigarOperator::I:
        case CigarOperator::S:
        case CigarOperator::EQ:
        case CigarOperator::X:
          read_length += length;
          break;
        default:
          break;
      }
    }
    return read_length;
  }

  auto
  to_string() const {
    std::string cigar_string;
    for (auto elements : cigar_elements) cigar_string += elements.to_string();
    return cigar_string;
  }

  auto
  begin() {
    return cigar_elements.begin();
  }

  auto
  begin() const {
    return cigar_elements.begin();
  }

  auto
  end() {
    return cigar_elements.end();
  }

  auto
  end() const {
    return cigar_elements.end();
  }

  auto&
  front() {
    return cigar_elements.front();
  }

  auto&
  front() const {
    return cigar_elements.front();
  }

  auto&
  back() {
    return cigar_elements.back();
  }

  auto&
  back() const {
    return cigar_elements.back();
  }

  void
  reverse() {
    std::reverse(cigar_elements.begin(), cigar_elements.end());
  }

  bool
  contains(CigarOperator key) const {
    for (auto [size, op] : cigar_elements)
      if (key == op) return true;
    return false;
  }

  template<typename Archive>
  void
  serialize(Archive& ar, const unsigned int version) {
    ar& cigar_elements;
  }

  friend auto&
  operator<<(std::ostream& os, const Cigar& cigar);

 private:
  std::vector<CigarElement> cigar_elements;
};

auto&
operator<<(std::ostream& os, const Cigar& cigar) {
  for (auto [length, op] : cigar) os << length << char(op);
  return os;
}

auto&
operator>>(std::istream& is, Cigar& cigar) {
  std::string cigar_string;
  is >> cigar_string;
  cigar = cigar_string;
  return is;
}

}  // namespace hc