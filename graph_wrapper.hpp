#pragma once

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <string_view>
#include <vector>

#include "haplotype.hpp"
#include "quality_utils.hpp"
#include "sam.hpp"
#include "smithwaterman.hpp"

namespace hc {

struct GraphWrapper {
  static constexpr std::size_t DEFAULT_NUM_PATHS = 128;
  static constexpr char MIN_BASE_QUALITY_TO_USE  = 10 + QualityUtils::ASCII_OFFSET;
  static constexpr std::size_t PRUNE_FACTOR      = 2;

 private:
  struct VertexProperty {
    std::string_view kmer;
  };

  struct EdgeProperty {
    std::size_t count = 0;
    bool is_ref       = false;
    bool is_on_path   = false;
    double score      = std::numeric_limits<double>::lowest();
  };

  using Graph  = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS,
                                      VertexProperty, EdgeProperty>;
  using Vertex = boost::graph_traits<Graph>::vertex_descriptor;
  using Edge   = boost::graph_traits<Graph>::edge_descriptor;
  using Path   = std::vector<Vertex>;

  struct cycle_detector : public boost::dfs_visitor<> {
    cycle_detector(bool& has_cycle) : _has_cycle(has_cycle) {}
    template<class Edge, class Graph>
    void
    back_edge(Edge, Graph&) {
      _has_cycle = true;
    }

   protected:
    bool& _has_cycle;
  };

  struct EdgeFilter {
    bool
    operator()(Edge e) const {
      return (*g)[e].is_ref || (*g)[e].count >= PRUNE_FACTOR ||
             boost::out_degree(boost::source(e, *g), *g) == 1;
    }
    Graph* g;
  } filter{&g};

  std::size_t kmer_size;
  Graph g;
  Vertex source{}, sink{};
  std::vector<Path> paths;
  std::set<Vertex> vertices_on_paths;

  std::string_view ref;
  std::vector<std::string_view> read_segs;

  std::set<std::string_view> dup_kmers;
  std::map<std::string_view, Vertex> unique_kmers;

  auto
  create_edge(Vertex u, Vertex v, bool is_ref) {
    auto [e, success] = boost::add_edge(u, v, g);
    g[e].count++;
    g[e].is_ref = is_ref;
  }

  auto
  create_vertex(std::string_view kmer) {
    auto v    = boost::add_vertex(g);
    g[v].kmer = kmer;
    if (dup_kmers.find(kmer) == dup_kmers.end()) unique_kmers.emplace(kmer, v);
    return v;
  }

  auto
  get_vertex(std::string_view kmer) {
    if (auto it = unique_kmers.find(kmer); it != unique_kmers.end()) return it->second;
    return create_vertex(kmer);
  }

  void
  increase_counts_backwards(Vertex v, std::string_view kmer) {
    if (kmer.empty()) return;
    if (boost::in_degree(v, g) == 1) {
      for (auto e : boost::make_iterator_range(boost::in_edges(v, g))) {
        auto u = boost::source(e, g);
        if (g[u].kmer.back() == kmer.back()) {
          g[e].count++;
          increase_counts_backwards(u, kmer.substr(0, kmer.size() - 1));
        }
      }
    }
  }

  auto
  extend_chain(Vertex u, std::string_view kmer, bool is_ref) {
    for (auto e : boost::make_iterator_range(boost::out_edges(u, g))) {
      auto v = boost::target(e, g);
      if (g[v].kmer.back() == kmer.back()) {
        g[e].count++;
        return v;
      }
    }

    auto v = get_vertex(kmer);
    create_edge(u, v, is_ref);
    return v;
  }

  void
  add_seq(std::string_view seq, bool is_ref) {
    auto v = get_vertex(seq.substr(0, kmer_size));
    // increase_counts_backwards(v, seq.substr(0, kmer_size-1));
    if (is_ref) source = v;
    for (auto i = 1; i <= seq.size() - kmer_size; i++)
      v = extend_chain(v, seq.substr(i, kmer_size), is_ref);
    if (is_ref) sink = v;
  }

  void
  path_finder(Vertex from, Vertex to, Path& path) {
    path.push_back(from);
    if (from == to) {
      paths.push_back(path);
      for (auto v : path) vertices_on_paths.insert(v);
    } else {
      for (auto e : boost::make_iterator_range(boost::out_edges(from, g))) {
        if (g[e].is_ref || g[e].count >= PRUNE_FACTOR || boost::out_degree(from, g) == 1) {
          auto v = target(e, g);
          if (std::find(path.begin(), path.end(), v) == path.end()) path_finder(v, to, path);
        }
      }
    }
    path.pop_back();
  }

  void
  find_all_paths() {
    Path path;
    path_finder(source, sink, path);
  }

  void
  mark_edges_on_paths() {
    for (const auto& path : paths) {
      auto u = path[0];
      for (std::size_t i = 1; i < path.size(); i++) {
        auto v                                   = path[i];
        g[boost::edge(u, v, g).first].is_on_path = true;
        u                                        = v;
      }
    }
  }

  void
  compute_edges_score() {
    for (auto v : vertices_on_paths) {
      std::vector<Edge> edges;
      for (auto e : boost::make_iterator_range(boost::out_edges(v, g)))
        if (g[e].is_on_path) edges.push_back(e);
      double sum = 0;
      for (auto e : edges) sum += g[e].count;
      for (auto e : edges) g[e].score = std::log10(g[e].count / sum);
    }
  }

  auto
  get_haplotypes() const {
    std::vector<Haplotype> haplotypes;

    for (const auto& path : paths) {
      auto u = path[0];
      std::string seq(g[u].kmer.data(), g[u].kmer.size());
      double score = 0;
      for (std::size_t i = 1; i < path.size(); i++) {
        auto v = path[i];
        seq += g[v].kmer.back();
        score += g[boost::edge(u, v, g).first].score;
        u = v;
      }
      haplotypes.emplace_back(std::move(seq), score);
    }

    std::sort(haplotypes.begin(), haplotypes.end(),
              [](const auto& h1, const auto& h2) { return h1.score > h2.score; });

    if (haplotypes.size() > DEFAULT_NUM_PATHS)
      haplotypes.erase(haplotypes.begin() + DEFAULT_NUM_PATHS, haplotypes.end());
    if (haplotypes.size() > 1)
      std::cout << "Found " << haplotypes.size() << " candidate haplotypes.\n";
    else
      std::cout << "Found only the reference haplotype in the assembly graph.\n";

    SWAligner aligner;
    for (auto& h : haplotypes) {
      auto [alignment_begin, cigar] = aligner.align(ref, h.bases);
      // std::cout << "Adding haplotype " << cigar << " from graph with kmer " << kmer_size << '\n';
      h.alignment_begin_wrt_ref     = alignment_begin;
      h.cigar                       = std::move(cigar);
    }

    // for (auto& h : haplotypes)
    // {
    //     std::cout << h.bases << '\n';
    //     std::cout << "> Cigar = " << h.cigar << " score " << h.score << '\n';
    // }

    return haplotypes;
  }

 public:
  static auto
  get_dup_kmers(std::string_view seq, std::size_t size) {
    std::set<std::string_view> all_kmers, dup_kmers;
    for (std::size_t i = 0; i <= seq.size() - size; i++) {
      auto kmer = seq.substr(i, size);
      if (auto [iter, success] = all_kmers.insert(kmer); !success) dup_kmers.insert(kmer);
    }
    return dup_kmers;
  }

  GraphWrapper(std::size_t kmer_size) : kmer_size(kmer_size) {}

  void
  set_ref(std::string_view ref) {
    this->ref = ref;
  }
  void
  set_read(const SAMRecord& read) {
    auto seq         = static_cast<std::string_view>(read.SEQ);
    const auto& qual = read.QUAL;

    auto start       = std::string_view::npos;
    auto is_usable   = [this](auto base, auto qual) {
      return base != 'N' && qual >= MIN_BASE_QUALITY_TO_USE;
    };
    for (std::size_t i = 0; i <= seq.size(); i++) {
      if (i == seq.size() || !is_usable(seq[i], qual[i])) {
        if (start != std::string_view::npos && i - start >= kmer_size)
          read_segs.push_back(seq.substr(start, i - start));
        start = std::string_view::npos;
      } else if (start == std::string_view::npos)
        start = i;
    }
  }

  void
  build() {
    for (auto kmer : get_dup_kmers(ref, kmer_size)) dup_kmers.insert(kmer);

    for (auto seg : read_segs)
      for (auto kmer : get_dup_kmers(seg, kmer_size)) dup_kmers.insert(kmer);

    add_seq(ref, true);
    for (auto seg : read_segs) add_seq(seg, false);
  }

  bool
  has_cycles() const {
    bool has_cycle = false;
    cycle_detector vis(has_cycle);
    boost::filtered_graph<Graph, EdgeFilter> fg(g, filter);
    boost::depth_first_search(fg, boost::visitor(vis));
    return has_cycle;
  }

  auto
  unique_kmers_count() const {
    return unique_kmers.size();
  }

  auto
  find_paths() {
    find_all_paths();
    mark_edges_on_paths();
    compute_edges_score();
    return get_haplotypes();
  }

  void
  print() const {
    std::ofstream os("graph.dot");
    os << "digraph assembly_graphs {";
    for (auto e : boost::make_iterator_range(boost::edges(g))) {
      os << boost::source(e, g) << " -> " << boost::target(e, g) << " ";
      auto count = g[e].count;
      if (g[e].is_ref)
        os << "[label=" << count << ",color=red];\n";
      else if (count < PRUNE_FACTOR)
        os << "[label=" << count << ",style=dotted,color=grey];\n";
      else
        os << "[label=" << count << "];\n";
    }

    for (auto v : boost::make_iterator_range(boost::vertices(g))) {
      os << v << " ";
      auto kmer = g[v].kmer;
      if (boost::in_degree(v, g) == 0)
        os << "[label=" << kmer << ",shape=box]\n";
      else
        os << "[label=" << kmer.back() << ",shape=box]\n";
    }
    os << "}";
  }
};

}  // namespace hc
