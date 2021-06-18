#pragma once

#include "interval.hpp"
#include "sam.hpp"

namespace hc {

struct ReadClipper {
  static void
  hard_clip_soft_clipped_bases(SAMRecord& read) {
    auto& seq                     = read.SEQ;
    auto& qual                    = read.QUAL;
    const auto& cigar             = read.CIGAR;

    auto [front_length, front_op] = cigar.front();
    if (front_op == CigarOperator::S) {
      seq  = seq.substr(front_length);
      qual = qual.substr(front_length);
    }

    auto [back_length, back_op] = cigar.back();
    if (back_op == CigarOperator::S) {
      seq  = seq.substr(0, seq.size() - back_length);
      qual = qual.substr(0, qual.size() - back_length);
    }
  }

  static void
  revert_soft_clipped_bases(SAMRecord& read) {
    auto& seq   = read.SEQ;
    auto& qual  = read.QUAL;
    auto& cigar = read.CIGAR;

    if (read.READ_REVERSE_STRAND()) {
      auto [front_length, front_op] = cigar.front();
      if (front_op == CigarOperator::S) {
        seq  = seq.substr(front_length);
        qual = qual.substr(front_length);
      }
      auto [back_length, back_op] = cigar.back();
      if (back_op == CigarOperator::S) cigar.back() = {back_length, CigarOperator::M};
    } else {
      auto [front_length, front_op] = cigar.front();
      auto alignment_begin          = read.get_alignment_begin();
      if (front_op == CigarOperator::S && alignment_begin >= front_length) {
        cigar.front() = {front_length, CigarOperator::M};
        read.POS      = alignment_begin - front_length + 1;
      }
      auto [back_length, back_op] = cigar.back();
      if (back_op == CigarOperator::S) {
        seq  = seq.substr(0, seq.size() - back_length);
        qual = qual.substr(0, qual.size() - back_length);
      }
    }
  }

  static void
  hard_clip_to_interval(SAMRecord& read, const Interval& interval) {
    auto& seq                        = read.SEQ;
    auto& qual                       = read.QUAL;

    const auto& [contig, begin, end] = interval;
    assert(read.RNAME == contig);

    auto alignment_begin = read.get_alignment_begin();
    auto alignment_end   = read.get_alignment_end();
    if (alignment_begin < begin) {
      auto clip_size = begin - alignment_begin;
      if (clip_size > seq.size()) clip_size = seq.size();
      seq  = seq.substr(clip_size);
      qual = qual.substr(clip_size);
    }
    if (alignment_end > end) {
      auto clip_size = alignment_end - end;
      seq            = seq.substr(0, seq.size() - clip_size);
      qual           = qual.substr(0, qual.size() - clip_size);
    }
  }
};

}  // namespace hc
