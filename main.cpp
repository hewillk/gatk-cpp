#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "assembler.hpp"
#include "interval.hpp"
#include "read_clipper.hpp"
#include "read_filter.hpp"
#include "sam.hpp"
#include "smithwaterman.hpp"
//#include "pairhmm.hpp"
#include <random>

#include "genetyper.hpp"
#include "intel_pairhmm.hpp"
#include "parser/parser.hpp"

bool first_shrink = true;
using namespace hc;
using namespace kaiba::parser;

struct ReadBuf {
  explicit ReadBuf(std::vector<SAMRecord> reads, const Interval& range) : begin(range.begin) {
    buf.resize(range.size());
    for (auto& read : reads) {
      if (read.get_alignment_begin() < begin) continue;
      auto i = read.get_alignment_begin() - begin;
      if (buf[i].size() < pos_per_reads) buf[i].emplace_back(std::move(read));
    }
  }

  auto
  shrink() {
    for (auto i = 0; i < (first_shrink ? 200 : 300); i++) {
      buf[shrink_begin].clear();
      buf[shrink_begin].shrink_to_fit();
      shrink_begin++;
    }
    first_shrink = false;
    std::cout << shrink_begin << '\n';
  }

  auto
  load(const Interval& range) {
    std::vector<SAMRecord> reads;
    for (auto i = range.begin; i < range.end; i++)
      for (const auto& read : buf[i - begin]) reads.emplace_back(read);
    return reads;
  }

  std::size_t shrink_begin = 0;
  int begin;
  std::size_t pos_per_reads = 5;
  std::vector<std::vector<SAMRecord>> buf;
};

void
call_region(std::vector<SAMRecord>& reads, std::string_view ref, const Interval& padded_region,
            const Interval& origin_region, std::ostream& os);

void
do_work(const std::string& ref_path, const std::string& bam_path, const std::string& vcf_path,
        Interval chr, int calling = 300, int padding = 100, int load_count = 5000) {
  const auto ref = Fa{ref_path}.load(chr.contig);
  auto bam       = Bam{bam_path};
  auto ofs       = std::ofstream{vcf_path};
  chr.begin      = 0;
  chr.end        = int(ref.size());

  auto called    = Interval{chr.contig, chr.begin, chr.begin + calling};
  auto padded    = Interval{chr.contig, chr.begin, chr.begin + calling + padding};
  auto loaded    = Interval{chr.contig, chr.begin, chr.begin + calling * load_count + padding};
  auto read_buf  = ReadBuf{bam.load(loaded.to_string()), loaded};
  auto passed    = 0;

  do {
    auto sub_ref = ref.substr(size_t(padded.begin), padded.size());
    auto reads   = read_buf.load(padded);
    if (reads.empty())
      std::cout << "Ignore " << called.to_string()
                << ":    (with overlap region = " << padded.to_string() << ")\n";
    else
      call_region(reads, sub_ref, padded, called, ofs);
    std::cout << "loaded: " << loaded.to_string() << '\n';
    std::cout << "called: " << called.to_string() << ", " << padded.to_string() << '\n';

    called.begin += calling;
    called.end += calling;
    padded.begin = called.begin - padding;
    padded.end   = called.end + padding;

    read_buf.shrink();
    if (++passed == load_count) {
      loaded   = Interval{chr.contig, padded.begin, called.begin + calling * load_count + padding};
      read_buf = ReadBuf{bam.load(loaded.to_string()), loaded};
      passed   = 0;
    }

  } while (called.begin < chr.end);
  std::cout << "HaplotypeCaller done." << '\n';
}

int
main(int argc, char* argv[]) {
  do_work(argv[1], argv[2], argv[3], argv[4]);
  // do_work("chrM:15680-16571");
  // do_work("chrM:245-490");
}

void
filter_reads(std::vector<SAMRecord>& reads) {
  reads.erase(std::remove_if(reads.begin(), reads.end(), MappingQualityReadFilter{}), reads.end());
  reads.erase(std::remove_if(reads.begin(), reads.end(), DuplicateReadFilter{}), reads.end());
  //    reads.erase(std::remove_if(reads.begin(), reads.end(),
  //        SecondaryAlignmentReadFilter{}
  //    ), reads.end());
  //    reads.erase(std::remove_if(reads.begin(), reads.end(),
  //        MateOnSameContigReadFilter{}
  //    ), reads.end());
}

void
hard_clip_reads(std::vector<SAMRecord>& reads, const Interval& padded_region) {
  std::for_each(reads.begin(), reads.end(),
                [](auto& read) { ReadClipper::revert_soft_clipped_bases(read); });
  std::for_each(reads.begin(), reads.end(),
                [&](auto& read) { ReadClipper::hard_clip_to_interval(read, padded_region); });
  reads.erase(std::remove_if(reads.begin(), reads.end(), MinimumLengthReadFilter{}), reads.end());
}

void
sample_reads(std::vector<SAMRecord>& reads, std::size_t n) {
  if (reads.size() > n) {
    std::vector<SAMRecord> temp;
    std::sample(reads.begin(), reads.end(), std::back_inserter(temp), n,
                std::mt19937{std::random_device{}()});
    reads.swap(temp);
  }
}

void
call_region(std::vector<SAMRecord>& reads, std::string_view ref, const Interval& padded_region,
            const Interval& origin_region, std::ostream& os) {
  Assembler assembler;
  IntelPairHMM pairhmm;
  Genetyper genetyper;
  std::cout
    << "----------------------------------------------------------------------------------\n";
  filter_reads(reads);
  hard_clip_reads(reads, padded_region);
  std::cout
    << "----------------------------------------------------------------------------------\n";
  if (reads.empty()) return;
  std::cout
    << "----------------------------------------------------------------------------------\n";
  std::cout << "Assembling " << origin_region.to_string() << " with " << reads.size()
            << " reads:    (with overlap region = " << padded_region.to_string() << ")\n";

  auto haplotypes = assembler.assemble(reads, ref);
  if (haplotypes.size() <= 1) return;

  auto likelihoods = pairhmm.compute_likelihoods(haplotypes, reads);
  auto variants    = genetyper.assign_genotype_likelihoods(reads, haplotypes, likelihoods, ref,
                                                        padded_region, origin_region);
  for (const auto& variant : variants) {
    os << variant.to_string() << '\n';
    std::cout << variant.to_string() << '\n';
  }
}
