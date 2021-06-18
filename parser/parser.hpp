#pragma once

#include <htslib/sam.h>

#include <array>
#include <fstream>
#include <map>
#include <sstream>
#include <string>

struct Bam {
  explicit Bam(const std::string& bam_path);
  std::vector<hc::SAMRecord>
  load(const std::string& region);
  ~Bam();

 private:
  samFile* bam_in;
  bam_hdr_t* bam_header;
  bam1_t* aln;
  hts_idx_t* idx;
  static constexpr auto bam_op_table = std::array{'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'};
};

Bam::Bam(const std::string& bam_path) {
  // open BAM for reading
  bam_in = sam_open(bam_path.c_str(), "r");
  assert(bam_in);
  // load the index
  idx = sam_index_load(bam_in, bam_path.c_str());
  assert(idx);
  // get the header
  bam_header = sam_hdr_read(bam_in);
  assert(bam_header);
  // Initiate the alignment record
  aln = bam_init1();
}
std::vector<SAMRecord>
Bam::load(const std::string& region) {
  std::vector<SAMRecord> reads;

  // move the iterator to the region we are interested in
  auto iter = sam_itr_querys(idx, bam_header, region.c_str());
  assert(iter);
  while (sam_itr_next(bam_in, iter, aln) >= 0) {
    auto qname  = bam_get_qname(aln);
    auto flag   = aln->core.flag;
    auto rname  = bam_header->target_name[aln->core.tid];
    auto pos    = aln->core.pos;
    auto mapq   = static_cast<short>(aln->core.qual);
    auto quali  = bam_get_qual(aln);
    auto seqi   = bam_get_seq(aln);
    auto cigari = bam_get_cigar(aln);
    auto seq    = std::string{};
    auto qual   = std::string{};
    auto cigar  = std::string{};

    seq.reserve(static_cast<std::size_t>(aln->core.l_qseq));
    qual.reserve(static_cast<std::size_t>(aln->core.l_qseq));
    for (auto i = 0u; i < aln->core.l_qseq; i++) {
      seq += seq_nt16_str[bam_seqi(seqi, i)];
      qual += static_cast<char>(33 + quali[i]);
    }
    hc::SAMRecord read;
    read.QNAME = qname;

    for (auto i = 0u; i < aln->core.n_cigar; i++) {
      auto icigar = cigari[i];
      auto len    = bam_cigar_oplen(icigar);
      auto op     = bam_op_table[bam_cigar_op(icigar)];
      cigar += std::to_string(len) += op;
    }

    read.FLAG  = flag;
    read.RNAME = rname;
    // ****************
    read.POS   = pos + 1;
    read.MAPQ  = mapq;
    read.CIGAR = cigar;
    read.SEQ   = std::move(seq);
    read.QUAL  = std::move(qual);
    reads.emplace_back(std::move(read));
  }
  hts_itr_destroy(iter);

  return reads;
}
Bam::~Bam() {
  hts_idx_destroy(idx);
  bam_destroy1(aln);
  bam_hdr_destroy(bam_header);
  sam_close(bam_in);
}
