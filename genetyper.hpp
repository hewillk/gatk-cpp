#include <numeric>
#include <set>

#include "haplotype.hpp"
#include "interval.hpp"
#include "math_utils.hpp"
#include "sam.hpp"

namespace hc {

class Genetyper {
  static inline const std::string SPAN_DEL                = "*";
  static constexpr std::size_t ALLELE_EXTENSION           = 2;
  static constexpr std::size_t MAX_GENOTYPE_QUALITY       = 99;
  static constexpr std::size_t MIN_HETEROZYGOSITY_QUALITY = 50;
  static constexpr std::size_t MAX_ALLELE_COUNT           = 7;

 private:
  static const inline auto allele_index_cache = [] {
    std::vector<std::vector<std::pair<std::size_t, std::size_t>>> cache{};
    for (std::size_t allele_count = 0; allele_count <= MAX_ALLELE_COUNT; allele_count++) {
      std::vector<std::pair<std::size_t, std::size_t>> inner_cache;
      for (std::size_t a1 = 0; a1 < allele_count; a1++)
        for (std::size_t a2 = a1; a2 < allele_count; a2++) inner_cache.emplace_back(a1, a2);
      cache.push_back(std::move(inner_cache));
    }
    return cache;
  }();

  auto
  process_cigar_for_initial_events(Haplotype& haplotype, std::string_view ref,
                                   const Interval& padded_region) {
    const auto& cigar                              = haplotype.cigar;
    const auto& hap                                = haplotype.bases;
    const auto& [contig, padded_begin, padded_end] = padded_region;

    auto ref_pos                                   = haplotype.alignment_begin_wrt_ref;
    auto hap_pos                                   = 0;
    for (auto [length, op] : cigar) {
      switch (op) {
        case CigarOperator::M: {
          std::vector<int> mismatch_offsets;
          for (auto offset = 0; offset < length; offset++)
            if (ref[ref_pos + offset] != hap[hap_pos + offset]) mismatch_offsets.push_back(offset);

          if (!mismatch_offsets.empty()) {
            for (auto offset : mismatch_offsets) {
              Variant variant;
              variant.REF         = ref[ref_pos + offset];
              variant.ALT         = hap[hap_pos + offset];
              auto mismatch_begin = padded_begin + ref_pos + offset;
              variant.location    = {contig, mismatch_begin, mismatch_begin + 1};
              haplotype.event_map.emplace(mismatch_begin, std::move(variant));
            }
          }

          ref_pos += length;
          hap_pos += length;
          break;
        }
        case CigarOperator::I: {
          if (ref_pos > 0) {
            Variant variant;
            variant.REF          = ref[ref_pos - 1];
            variant.ALT          = variant.REF + hap.substr(hap_pos, length);
            auto insertion_begin = padded_begin + ref_pos - 1;
            variant.location     = {contig, insertion_begin, insertion_begin + 1};
            haplotype.event_map.emplace(insertion_begin, std::move(variant));
          }
          hap_pos += length;
          break;
        }
        case CigarOperator::D: {
          if (ref_pos > 0) {
            Variant variant;
            variant.REF         = ref.substr(ref_pos - 1, length + 1);
            variant.ALT         = ref[ref_pos - 1];
            auto deletion_begin = padded_begin + ref_pos - 1;
            variant.location    = {contig, deletion_begin, deletion_begin + length + 1};
            haplotype.event_map.emplace(deletion_begin, std::move(variant));
          }
          ref_pos += length;
          break;
        }
        case CigarOperator::S: {
          hap_pos += length;
          break;
        }
        default:
          throw std::invalid_argument("Unsupported cigar operator created during SW alignment");
      }
    }
  }

  auto
  set_events_for_haplotypes(std::vector<Haplotype>& haplotypes, std::string_view ref,
                            const Interval& padded_region) {
    std::set<std::size_t> events_begins;
    std::size_t rank = 0;
    for (auto& h : haplotypes) {
      h.rank = rank++;
      process_cigar_for_initial_events(h, ref, padded_region);
      for (const auto& [begin, even] : h.event_map) events_begins.insert(begin);
    }
    return events_begins;
  }

  auto
  get_events_from_haplotypes(std::size_t begin, std::vector<Haplotype>& haplotypes) {
    std::set<Variant> unique_events;
    for (const auto& h : haplotypes)
      for (auto& event : h.get_overlapping_events(begin)) unique_events.insert(std::move(event));
    std::vector<Variant> results;
    std::copy(unique_events.begin(), unique_events.end(), std::back_inserter(results));
    return results;
  }

  auto
  replace_span_dels(std::vector<Variant>& events, char ref_allele, int begin) {
    for (auto& event : events) {
      if (event.location.begin != begin) {
        Variant new_event;
        new_event.location = {event.location.contig, begin, begin + 1};
        new_event.REF      = ref_allele;
        new_event.ALT      = SPAN_DEL;
        event              = std::move(new_event);
      }
    }
  }

  auto
  determine_reference_allele(const std::vector<Variant>& events) {
    return std::max_element(
             events.begin(), events.end(),
             [](const auto& lhs, const auto& rhs) { return lhs.REF.size() < rhs.REF.size(); })
      ->REF;
  }

  auto
  get_compatible_alternate_allele(const std::string& ref_allele, const Variant& event) {
    if (event.ALT == SPAN_DEL) return SPAN_DEL;
    return event.ALT + ref_allele.substr(event.REF.size());
  }

  auto
  resolve_incompatible_alleles(const std::string& ref_allele, const Variant& event,
                               std::set<std::string>& alts) {
    if (event.REF == ref_allele)
      alts.insert(event.ALT);
    else
      alts.insert(get_compatible_alternate_allele(ref_allele, event));
  }

  auto
  get_compatible_alleles(const std::vector<Variant>& events) {
    auto longest_event = events.front();
    auto ref_allele    = determine_reference_allele(events);
    std::vector<std::string> alleles{ref_allele};
    std::set<std::string> alts;
    for (const auto& event : events) {
      if (event.size() > longest_event.size()) longest_event = event;
      resolve_incompatible_alleles(ref_allele, event, alts);
    }
    std::copy(alts.begin(), alts.end(), std::back_inserter(alleles));
    return std::make_pair(alleles, longest_event.location);
  }

  auto
  get_allele_mapper(const std::vector<std::string>& alleles, std::size_t begin,
                    const std::vector<Haplotype>& haplotypes) {
    std::map<std::size_t, std::vector<std::size_t>> result;
    result[0];
    const auto& ref_allele = alleles[0];
    auto get_index         = [&](const auto& allele) {
      return std::find(alleles.begin(), alleles.end(), allele) - alleles.begin();
    };
    for (const auto& h : haplotypes) {
      auto spanning_events = h.get_overlapping_events(begin);
      if (spanning_events.empty()) result[0].push_back(h.rank);
      for (const auto& event : spanning_events) {
        if (event.location.begin == begin) {
          if (event.REF.size() == ref_allele.size())
            result[get_index(event.ALT)].push_back(h.rank);
          else if (event.REF.size() < ref_allele.size())
            result[get_index(get_compatible_alternate_allele(ref_allele, event))].push_back(h.rank);
        } else
          result[get_index(SPAN_DEL)].push_back(h.rank);
      }
    }
    return result;
  }

  auto
  get_haplotype_mapper(const std::map<std::size_t, std::vector<std::size_t>>& allele_mapper,
                       std::size_t haplotype_count) {
    std::vector<std::size_t> haplotype_mapper(haplotype_count);
    for (const auto& [allele_index, haplotype_indices] : allele_mapper)
      for (auto haplotype_index : haplotype_indices)
        haplotype_mapper[haplotype_index] = allele_index;
    return haplotype_mapper;
  }

  auto
  get_read_indices_to_keep(const std::vector<SAMRecord>& reads, const Interval& overlap) {
    std::vector<std::size_t> read_indices_to_keep;
    read_indices_to_keep.reserve(reads.size());
    for (std::size_t i = 0; i < reads.size(); i++)
      if (reads[i].get_interval().overlaps(overlap)) read_indices_to_keep.push_back(i);
    return read_indices_to_keep;
  }

  auto
  marginal_likelihoods(std::size_t allele_count, const std::vector<std::size_t>& haplotype_mapper,
                       const std::vector<std::size_t>& read_indices_to_keep,
                       const std::vector<std::vector<double>>& haplotype_likelihoods) {
    std::vector<std::vector<double>> allele_likelihoods(
      read_indices_to_keep.size(),
      std::vector<double>(allele_count, std::numeric_limits<double>::lowest()));
    for (std::size_t r = 0; r < read_indices_to_keep.size(); r++) {
      auto old_read_index = read_indices_to_keep[r];
      for (std::size_t h = 0; h < haplotype_mapper.size(); h++) {
        auto allele_index = haplotype_mapper[h];
        auto likelihood   = haplotype_likelihoods[old_read_index][h];
        if (likelihood > allele_likelihoods[r][allele_index])
          allele_likelihoods[r][allele_index] = likelihood;
      }
    }
    return allele_likelihoods;
  }

  auto
  marginalize(const std::vector<std::size_t>& haplotype_mapper, std::size_t allele_count,
              const std::vector<SAMRecord>& reads,
              const std::vector<std::vector<double>>& haplotype_likelihoods,
              const Interval& overlap) {
    auto read_indices_to_keep = get_read_indices_to_keep(reads, overlap);
    return marginal_likelihoods(allele_count, haplotype_mapper, read_indices_to_keep,
                                haplotype_likelihoods);
  }

  void
  single_component_genotype_likelihood_by_read(
    std::vector<double>& genotype_likelihoods,
    const std::vector<std::vector<double>>& allele_likelihoods, std::size_t a) {
    const auto log10_frequency = std::log10(2);
    std::transform(allele_likelihoods.begin(), allele_likelihoods.end(),
                   std::back_inserter(genotype_likelihoods),
                   [=](const auto& likelihoods) { return likelihoods[a] + log10_frequency; });
  }

  void
  two_component_genotype_likelihood_by_read(
    std::vector<double>& genotype_likelihoods,
    const std::vector<std::vector<double>>& allele_likelihoods, std::size_t a1, std::size_t a2) {
    std::transform(allele_likelihoods.begin(), allele_likelihoods.end(),
                   std::back_inserter(genotype_likelihoods), [=](const auto& likelihoods) {
                     return MathUtils::approximate_log10_sum_log10(likelihoods[a1],
                                                                   likelihoods[a2]);
                   });
  }

  auto
  calculate_read_likelihoods_by_genotype_index(
    const std::vector<std::vector<double>>& allele_likelihoods, std::size_t allele_count) {
    std::vector<std::vector<double>> read_likelihoods_by_genotype_index((allele_count + 1) *
                                                                        allele_count / 2);
    std::size_t cur_genotype_index = 0;
    for (std::size_t a1 = 0; a1 < allele_count; a1++) {
      for (std::size_t a2 = a1; a2 < allele_count; a2++) {
        auto& read_genotype_likelihoods = read_likelihoods_by_genotype_index[cur_genotype_index++];
        read_genotype_likelihoods.reserve(allele_likelihoods.size());
        if (a1 == a2)
          single_component_genotype_likelihood_by_read(read_genotype_likelihoods,
                                                       allele_likelihoods, a1);
        else
          two_component_genotype_likelihood_by_read(read_genotype_likelihoods, allele_likelihoods,
                                                    a1, a2);
      }
    }
    return read_likelihoods_by_genotype_index;
  }

  auto
  get_genotype_likelihoods(
    const std::vector<std::vector<double>>& read_likelihoods_by_genotype_index) {
    const auto genotype_count = read_likelihoods_by_genotype_index.size();
    std::vector<double> result(genotype_count);
    const auto denominator = read_likelihoods_by_genotype_index[0].size() * std::log10(2);
    for (std::size_t genotype = 0; genotype < genotype_count; genotype++)
      result[genotype] = std::accumulate(read_likelihoods_by_genotype_index[genotype].begin(),
                                         read_likelihoods_by_genotype_index[genotype].end(), 0.0) -
                         denominator;
    return result;
  }

  auto
  calculate_genotype_likelihoods(const std::vector<std::vector<double>>& allele_likelihoods,
                                 std::size_t allele_count) {
    auto read_likelihoods_by_genotype_index =
      calculate_read_likelihoods_by_genotype_index(allele_likelihoods, allele_count);
    return get_genotype_likelihoods(read_likelihoods_by_genotype_index);
  }

  auto
  get_genotype_quality_and_max_genotype_index(const std::vector<double>& genotypes) {
    double max, second_max;
    std::size_t max_index;
    if (genotypes[0] > genotypes[1]) {
      second_max = genotypes[1];
      max        = genotypes[0];
      max_index  = 0;
    } else {
      second_max = genotypes[0];
      max        = genotypes[1];
      max_index  = 1;
    }
    for (std::size_t i = 2; i < genotypes.size(); i++) {
      if (genotypes[i] >= max) {
        second_max = max;
        max        = genotypes[i];
        max_index  = i;
      } else if (genotypes[i] > second_max) {
        second_max = genotypes[i];
      }
    }
    auto genotype_quality = static_cast<std::size_t>(std::round(-10 * (second_max - max)));
    if (genotype_quality > MAX_GENOTYPE_QUALITY) genotype_quality = MAX_GENOTYPE_QUALITY;
    return std::make_pair(max_index, genotype_quality);
  }

  auto
  get_genotype(std::size_t allele_count, std::size_t genotype_index) {
    return allele_index_cache[allele_count][genotype_index];
  }

 public:
  auto
  assign_genotype_likelihoods(const std::vector<SAMRecord>& reads,
                              std::vector<Haplotype>& haplotypes,
                              const std::vector<std::vector<double>>& haplotype_likelihoods,
                              std::string_view ref, const Interval& padded_region,
                              const Interval& origin_region) {
    auto events_begins = set_events_for_haplotypes(haplotypes, ref, padded_region);
    const auto& [contig, origin_begin, origin_end] = origin_region;
    std::vector<Variant> variants;
    for (auto begin : events_begins) {
      if (begin < origin_begin || begin >= origin_end) continue;
      auto events = get_events_from_haplotypes(begin, haplotypes);
      replace_span_dels(events, ref[begin - padded_region.begin], begin);
      auto [alleles, alleles_loc] = get_compatible_alleles(events);
      auto allele_count           = alleles.size();
      if (allele_count > MAX_ALLELE_COUNT) continue;
      auto allele_mapper    = get_allele_mapper(alleles, begin, haplotypes);
      auto haplotype_mapper = get_haplotype_mapper(allele_mapper, haplotypes.size());
      auto allele_likelihoods =
        marginalize(haplotype_mapper, allele_count, reads, haplotype_likelihoods,
                    alleles_loc.expand_within_contig(ALLELE_EXTENSION));
      auto genotype_likelihoods = calculate_genotype_likelihoods(allele_likelihoods, allele_count);
      auto [genotype_index, genotype_quality] =
        get_genotype_quality_and_max_genotype_index(genotype_likelihoods);
      if (genotype_index == 0) continue;
      auto genotype = get_genotype(allele_count, genotype_index);
      if (genotype.first == 0 && genotype_quality < MIN_HETEROZYGOSITY_QUALITY) continue;
      variants.push_back(Variant{
        .location = alleles_loc, .alleles = alleles, .GT = genotype, .GQ = genotype_quality});
    }
    return variants;
  }
};

}  // namespace hc
