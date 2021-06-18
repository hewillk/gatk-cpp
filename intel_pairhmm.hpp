#pragma once

#include <cstdio>
#ifdef _DEBUG
#define DBG(format, args...) printf("[%s:%d] " format "\n", __FILE__, __LINE__, ##args)
#else
#define DBG(args...)
#endif

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include "haplotype.hpp"
#include "native/avx-pairhmm.h"
#include "native/shacc_pairhmm.h"
#include "sam.hpp"

namespace hc {

struct IntelPairHMM {
 private:
  static constexpr double TRISTATE_CORRECTION                    = 3.0;
  static constexpr double MAXIMUM_BEST_ALT_LIKELIHOOD_DIFFERENCE = -4.5;
  static constexpr double EXPECTED_ERROR_RATE_PER_BASE           = 0.02;
  static constexpr double LOG10_QUALITY_PER_BASE                 = -4.0;
  static constexpr double MAXIMUM_EXPECTED_ERROR_PER_READ        = 2.0;
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

 public:
  auto
  compute_likelihoods(const std::vector<Haplotype>& haplotypeDataArray,
                      std::vector<SAMRecord>& readDataArray) {
    initNative();
    std::vector<std::vector<double>> likelihoodArray(
      readDataArray.size(), std::vector<double>(haplotypeDataArray.size()));
    computeLikelihoodsNative(readDataArray, haplotypeDataArray, likelihoodArray);
    normalize_likelihoods_and_filter_poorly_modeled_reads(readDataArray, likelihoodArray);
    return likelihoodArray;
  }

 private:
  bool g_use_double;
  int g_max_threads;

  Context<float> g_ctxf;
  Context<double> g_ctxd;

  float (*g_compute_full_prob_float)(testcase* tc);
  double (*g_compute_full_prob_double)(testcase* tc);

  std::vector<std::vector<testcase>> m_testcases;

 private:
  std::vector<std::vector<testcase>>
  getData(const std::vector<SAMRecord>& readDataArray,
          const std::vector<Haplotype>& haplotypeDataArray);
  void
  initNative(bool use_double = false, int max_threads = 64);
  void
  computeLikelihoodsNative(const std::vector<SAMRecord>& readDataArray,
                           const std::vector<Haplotype>& haplotypeDataArray,
                           std::vector<std::vector<double>>& likelihoodArray);
};

void
IntelPairHMM::initNative(bool use_double, int max_threads) {
  DBG("Enter");

  g_use_double = use_double;
#ifdef _OPENMP
  int avail_threads = omp_get_max_threads();
  int req_threads   = max_threads;
  g_max_threads     = std::min(req_threads, avail_threads);

  DBG("Available threads: %d", avail_threads);
  DBG("Requested threads: %d", req_threads);
  if (req_threads > avail_threads) {
    DBG("Using %d available threads, but %d were requested", g_max_threads, req_threads);
  } else {
    DBG("Using %d threads", g_max_threads);
  }
#else
  if (max_threads != 1) {
    DBG("Ignoring request for %d threads; not using OpenMP implementation", max_threads);
  }
#endif

  // enable FTZ
  if (_MM_GET_FLUSH_ZERO_MODE() != _MM_FLUSH_ZERO_ON) {
    DBG("Flush-to-zero (FTZ) is enabled when running PairHMM");
  }
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

  g_compute_full_prob_float  = &compute_full_prob_avxs;
  g_compute_full_prob_double = &compute_full_prob_avxd;

  // init convert char table
  ConvertChar::init();
  DBG("Exit");
}

void
IntelPairHMM::computeLikelihoodsNative(const std::vector<SAMRecord>& readDataArray,
                                       const std::vector<Haplotype>& haplotypeDataArray,
                                       std::vector<std::vector<double>>& likelihoodArray) {
  DBG("Enter");

  //==================================================================
  // get data
  auto testcases = getData(readDataArray, haplotypeDataArray);

  //==================================================================
  // calcutate pairHMM

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(g_max_threads)
#endif
  for (int i = 0; i < testcases.size(); i++) {
    for (int j = 0; j < testcases[0].size(); j++) {
      double result_final = 0;

      float result_float  = g_use_double ? 0.0f : g_compute_full_prob_float(&testcases[i][j]);

      if (result_float < MIN_ACCEPTED) {
        double result_double = g_compute_full_prob_double(&testcases[i][j]);
        result_final         = log10(result_double) - g_ctxd.LOG10_INITIAL_CONSTANT;
      } else {
        result_final = (double)(log10f(result_float) - g_ctxf.LOG10_INITIAL_CONSTANT);
      }
      likelihoodArray[i][j] = result_final;
      DBG("result = %e", result_final);
    }
  }

  //==================================================================
  // release data
  DBG("Exit");
}

std::vector<std::vector<testcase>>
IntelPairHMM::getData(const std::vector<SAMRecord>& readDataArray,
                      const std::vector<Haplotype>& haplotypeDataArray) {
  int numReads      = readDataArray.size();
  int numHaplotypes = haplotypeDataArray.size();

  std::vector<const char*> haplotypes;
  std::vector<int> haplotypeLengths;

  long total_hap_length  = 0;
  long total_read_length = 0;

  // get haplotypes
  for (int i = 0; i < numHaplotypes; i++) {
    int length = haplotypeDataArray[i].bases.length();
    haplotypes.push_back(haplotypeDataArray[i].bases.data());
    haplotypeLengths.push_back(length);
    total_hap_length += length;
  }

  // get reads and create testcases
  for (int r = 0; r < numReads; r++) {
    int length            = readDataArray[r].SEQ.length();
    const char* reads     = readDataArray[r].SEQ.data();
    int readLength        = length;
    const char* insGops   = readDataArray[r].insertionGOP().data();
    const char* delGops   = readDataArray[r].deletionGOP().data();
    const char* gapConts  = readDataArray[r].overallGCP().data();
    const char* readQuals = readDataArray[r].QUAL.data();
    total_read_length += length;

    std::vector<testcase> n_testcases;
    for (int h = 0; h < numHaplotypes; h++) {
      testcase tc;
      tc.hap    = haplotypes[h];
      tc.haplen = haplotypeLengths[h];
      tc.rs     = reads;
      tc.rslen  = readLength;
      tc.i      = insGops;
      tc.d      = delGops;
      tc.c      = gapConts;
      tc.q      = readQuals;
      n_testcases.push_back(tc);
    }
    m_testcases.emplace_back(std::move(n_testcases));
  }

  return m_testcases;
}

}  // namespace hc
