#pragma once

#include <iostream>

#include "graph_wrapper.hpp"
#include "haplotype.hpp"
#include "interval.hpp"
#include "quality_utils.hpp"
#include "sam.hpp"

namespace hc {

struct Assembler {
  static constexpr std::size_t INITIAL_KMER_SIZE                 = 25;
  static constexpr std::size_t KMER_SIZE_ITERATION_INCREASE      = 10;
  static constexpr std::size_t MAX_KMER_ITERATIONS_TO_ATTEMPT    = 6;
  static constexpr std::size_t MIN_UNIQUE_KMERS_COUNT_TO_DISCARD = 5000;

 private:
  std::vector<Haplotype>
  assemble(const std::vector<SAMRecord>& reads, std::string_view ref, std::size_t kmer_size) {
    if (ref.size() < kmer_size) return {};

    GraphWrapper graph(kmer_size);

    graph.set_ref(ref);
    for (const auto& read : reads) graph.set_read(read);

    graph.build();

    if (graph.unique_kmers_count() > MIN_UNIQUE_KMERS_COUNT_TO_DISCARD) {
      std::cout << "Not using kmer size of " << kmer_size
                << " in assembler because it contains too much unique kmers\n";
      return {};
    }

    if (graph.has_cycles()) {
      std::cout << "Not using kmer size of " << kmer_size
                << " in assembler because it contains a cycle\n";
      return {};
    }

    std::cout << "Using kmer size of " << kmer_size << " in assembler\n";

    // graph.print();

    return graph.find_paths();
  }

 public:
  auto
  assemble(const std::vector<SAMRecord>& reads, std::string_view ref) {
    std::size_t iterations = 1;
    std::size_t kmer_size  = INITIAL_KMER_SIZE;
    auto haplotypes        = assemble(reads, ref, kmer_size);
    while (haplotypes.empty() && iterations < MAX_KMER_ITERATIONS_TO_ATTEMPT) {
      iterations++;
      kmer_size += KMER_SIZE_ITERATION_INCREASE;
      haplotypes = assemble(reads, ref, kmer_size);
    }
    return haplotypes;
  }
};

}  // namespace hc
