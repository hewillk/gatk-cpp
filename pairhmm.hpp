#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include "haplotype.hpp"
#include "quality_utils.hpp"
#include "sam.hpp"

namespace hc {

struct PairHMM {
 private:
  enum { M_TO_M, M_TO_I, M_TO_D, I_TO_M, I_TO_I, D_TO_M, D_TO_D, TRANS_PROB_ARRAY_LENGTH };
  using TrasMatrx                             = std::array<double, TRANS_PROB_ARRAY_LENGTH>;
  static constexpr double TRISTATE_CORRECTION = 3.0;
  static constexpr double MAXIMUM_BEST_ALT_LIKELIHOOD_DIFFERENCE = -4.5;
  static constexpr double EXPECTED_ERROR_RATE_PER_BASE           = 0.02;
  static constexpr double LOG10_QUALITY_PER_BASE                 = -4.0;
  static constexpr double MAXIMUM_EXPECTED_ERROR_PER_READ        = 2.0;

 public:
  static constexpr TrasMatrx ORIGINAL_DEFAULT{0.9998, 0.0001, 0.0001, 0.9, 0.1, 0.9, 0.1};

 private:
  std::vector<std::vector<double>> M, I, D, p;
  std::size_t previous_haplotype_length = 0;

 public:
  auto
  compute_likelihoods(const std::vector<Haplotype>& haplotypes, std::vector<SAMRecord>& reads,
                      const TrasMatrx& t = ORIGINAL_DEFAULT) {
    auto find_padded_size = [](const auto& v) {
      return (*std::max_element(v.begin(), v.end(),
                                [](const auto& x, const auto& y) { return x.size() < y.size(); }))
               .size() +
             1;
    };
    M = std::vector<std::vector<double>>(find_padded_size(reads),
                                         std::vector<double>(find_padded_size(haplotypes)));
    I = M;
    D = M;
    p = M;

    for (auto& read : reads) modify_read_qualities(read);

    std::vector<std::vector<double>> log_likelihoods(reads.size(),
                                                     std::vector<double>(haplotypes.size()));
    for (std::size_t i = 0; i < reads.size(); i++)
      for (std::size_t j = 0; j < haplotypes.size(); j++)
        log_likelihoods[i][j] = sub_compute_likelihood(reads[i], haplotypes[j].bases, t);
    normalize_likelihoods_and_filter_poorly_modeled_reads(reads, log_likelihoods);
    return log_likelihoods;
  }

 private:
  static inline double INITIAL_CONDITION       = std::pow(2, 1020);
  static inline double INITIAL_CONDITION_LOG10 = std::log10(INITIAL_CONDITION);

  double
  sub_compute_likelihood(SAMRecord& read, const std::string& haplotype, const TrasMatrx& t) {
    if (previous_haplotype_length == 0 || previous_haplotype_length != haplotype.length()) {
      auto initial_value = INITIAL_CONDITION / haplotype.size();
      for (std::size_t j = 0; j <= haplotype.size(); j++) D[0][j] = initial_value;
      previous_haplotype_length = haplotype.length();
    }

    initialize_priors(read, haplotype);

    for (std::size_t i = 1; i <= read.size(); i++) {
      for (std::size_t j = 1; j <= haplotype.size(); j++) {
        M[i][j] = p[i][j] * (M[i - 1][j - 1] * t[M_TO_M] + I[i - 1][j - 1] * t[I_TO_M] +
                             D[i - 1][j - 1] * t[D_TO_M]);
        I[i][j] = M[i - 1][j] * t[M_TO_I] + I[i - 1][j] * t[I_TO_I];
        D[i][j] = M[i][j - 1] * t[M_TO_D] + D[i][j - 1] * t[D_TO_D];
      }
    }

    double final_sum_prob = 0;
    auto end_i            = read.size();
    for (std::size_t j = 1; j <= haplotype.size(); j++) final_sum_prob += M[end_i][j] + D[end_i][j];

    return std::log10(final_sum_prob) - INITIAL_CONDITION_LOG10;
  }

  void
  initialize_priors(SAMRecord& read, const std::string& haplotype) {
    for (std::size_t i = 0; i < read.size(); i++) {
      auto x    = read.SEQ[i];
      auto qual = read.QUAL[i];
      for (std::size_t j = 0; j < haplotype.size(); j++) {
        auto y          = haplotype[j];
        p[i + 1][j + 1] = (x == y || x == 'N' || y == 'N' ?
                             1 - QualityUtils::qual_to_error_prob(qual) :
                             (QualityUtils::qual_to_error_prob(qual) / TRISTATE_CORRECTION));
      }
    }
  }

  void
  modify_read_qualities(SAMRecord& read) {
    char mapq = QualityUtils::ASCII_OFFSET + (char)read.MAPQ;
    for (auto& qual : read.QUAL) qual = std::min(qual, mapq);
  }

  void
  normalize_likelihoods_and_filter_poorly_modeled_reads(
    std::vector<SAMRecord>& reads, std::vector<std::vector<double>>& log_likelihoods) {
    std::vector<std::size_t> remove_indices;
    for (std::size_t i = 0; i < log_likelihoods.size(); i++) {
      auto best_likelihood =
        *std::max_element(log_likelihoods[i].begin(), log_likelihoods[i].end());
      auto cap_likelihood = best_likelihood + MAXIMUM_BEST_ALT_LIKELIHOOD_DIFFERENCE;
      for (auto& likelihood : log_likelihoods[i])
        if (likelihood < cap_likelihood) likelihood = cap_likelihood;

      auto likelihood_threshold =
        std::min(MAXIMUM_EXPECTED_ERROR_PER_READ,
                 std::ceil(reads[i].size() * EXPECTED_ERROR_RATE_PER_BASE)) *
        LOG10_QUALITY_PER_BASE;
      if (best_likelihood < likelihood_threshold) remove_indices.push_back(i);
    }
    auto remove_by_sorted_indices = [](auto& v, const auto& indices) {
      for (auto it = indices.rbegin(); it != indices.rend(); ++it) v.erase(v.begin() + *it);
    };
    remove_by_sorted_indices(log_likelihoods, remove_indices);
    remove_by_sorted_indices(reads, remove_indices);
  }
};

}  // namespace hc
