#pragma once

#include <vector>
using Cigar = std::vector << std::pair<int, char> >> ;
namespace hc {

struct SWAligner {
 public:
  struct SWParameters {
    int w_match;
    int w_mismatch;
    int w_open;
    int w_extend;
  };

  // match=1, mismatch = -1/3, gap=-(1+k/3)
  static constexpr SWParameters ORIGINAL_DEFAULT{3, -1, -4, -3};
  static constexpr SWParameters STANDARD_NGS{25, -50, -110, -6};
  static constexpr SWParameters NEW_SW_PARAMETERS{200, -150, -260, -11};
  static constexpr SWParameters ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS{10, -15, -30, -5};

  static constexpr std::size_t MINIMAL_MISMATCH_TO_TOLERANCE = 2;

 private:
  using State = CigarOperator;

 public:
  // offset, Cigar
  std::pair<int, Cigar>
  align(std::string_view ref, std::string_view alt,
        const SWParameters& params = NEW_SW_PARAMETERS) {
    if (ref.empty() || alt.empty())
      throw std::invalid_argument("Non-null sequences are required for the SW aligner");

    if (is_all_match(ref, alt)) return {0, {1, {int(ref.size()), CigarOperator::M}}};

    std::vector<std::vector<int>> score(ref.size() + 1, std::vector<int>(alt.size() + 1, 0));
    std::vector<std::vector<int>> trace(ref.size() + 1, std::vector<int>(alt.size() + 1, 0));
    calculate_matrix(ref, alt, score, trace, params);
    return calculate_cigar(score, trace);
  }

 private:
  bool
  is_all_match(std::string_view ref, std::string_view alt) const {
    if (alt.size() == ref.size()) {
      std::size_t mismatch = 0;
      for (std::size_t i = 0; mismatch <= MINIMAL_MISMATCH_TO_TOLERANCE && i < ref.size(); i++)
        if (alt[i] != ref[i]) mismatch++;
      return mismatch <= MINIMAL_MISMATCH_TO_TOLERANCE;
    }
    return false;
  }

  void
  calculate_matrix(std::string_view ref, std::string_view alt, std::vector<std::vector<int>>& score,
                   std::vector<std::vector<int>>& trace, const SWParameters& params) {
    std::size_t row_size = score.size();
    std::size_t col_size = score[0].size();

    std::vector<int> gap_size_down(col_size + 1, 0);
    std::vector<int> best_gap_down(col_size + 1, std::numeric_limits<int>::min() / 2);
    std::vector<int> gap_size_right(row_size + 1, 0);
    std::vector<int> best_gap_right(row_size + 1, std::numeric_limits<int>::min() / 2);

    auto [w_match, w_mismatch, w_open, w_extend] = params;
    for (std::size_t i = 1; i < row_size; i++) {
      for (std::size_t j = 1; j < col_size; j++) {
        // step diag
        int step_diag     = score[i - 1][j - 1] + (ref[i - 1] == alt[j - 1] ? w_match : w_mismatch);

        // step down
        int gap_open_down = score[i - 1][j] + w_open;
        best_gap_down[j] += w_extend;
        if (gap_open_down > best_gap_down[j]) {
          best_gap_down[j] = gap_open_down;
          gap_size_down[j] = 1;
        } else
          gap_size_down[j]++;
        int step_down      = best_gap_down[j];
        int step_down_size = gap_size_down[j];

        // step right
        int gap_open_right = score[i][j - 1] + w_open;
        best_gap_right[i] += w_extend;
        if (gap_open_right > best_gap_right[i]) {
          best_gap_right[i] = gap_open_right;
          gap_size_right[i] = 1;
        } else
          gap_size_right[i]++;
        int step_right      = best_gap_right[i];
        int step_right_size = gap_size_right[i];

        // priority here will be step diagonal, step down, step right
        if (step_diag >= step_down && step_diag >= step_right) {
          score[i][j] = step_diag;
          trace[i][j] = 0;
        } else if (step_right >= step_down) {
          score[i][j] = step_right;
          trace[i][j] = -step_right_size;
        } else {
          score[i][j] = step_down;
          trace[i][j] = step_down_size;
        }
      }
    }
  }

  std::pair<size_t, Cigar>
  calculate_cigar(std::vector<std::vector<int>>& score, std::vector<std::vector<int>>& trace) {
    std::size_t ref_size       = score.size() - 1;
    std::size_t alt_size       = score[0].size() - 1;

    int max_score              = std::numeric_limits<int>::min();
    std::size_t segment_length = 0;

    // look for the largest score on the rightmost column. we use >= combined with the traversal
    // direction to ensure that if two scores are equal, the one closer to diagonal gets picked
    std::size_t pos_i          = 0;
    for (std::size_t i = 1; i <= ref_size; i++) {
      int cur_score = score[i][alt_size];
      if (cur_score >= max_score) {
        max_score = cur_score;
        pos_i     = i;
      }
    }

    // now look for a larger score on the bottom-most row
    std::size_t pos_j = alt_size;
    auto abs_diff     = [](auto x, auto y) { return x > y ? x - y : y - x; };
    for (std::size_t j = 1; j <= alt_size; j++) {
      int cur_score = score[ref_size][j];
      if (cur_score > max_score ||
          (cur_score == max_score && abs_diff(ref_size, j) < abs_diff(pos_i, pos_j))) {
        max_score      = cur_score;
        pos_i          = ref_size;
        pos_j          = j;
        // end of alternate is overhanging; we will just record it as 'M' segment
        segment_length = alt_size - j;
      }
    }

    Cigar cigar;
    if (segment_length > 0) {
      cigar.emplace_back(segment_length, State::S);
      segment_length = 0;
    }

    auto state = State::M;
    do {
      int cur_trace = trace[pos_i][pos_j], step_size;
      State new_state;
      if (cur_trace > 0) {
        new_state = State::D;
        step_size = cur_trace;
      } else if (cur_trace < 0) {
        new_state = State::I;
        step_size = -cur_trace;
      } else {
        new_state = State::M;
        step_size = 1;
      }

      // move to next best location in the sw matrix
      switch (new_state) {
        case State::M:
          pos_i--;
          pos_j--;
          break;
        case State::I:
          pos_j -= step_size;
          break;
        case State::D:
          pos_i -= step_size;
          break;
        default:
          break;
      }

      if (new_state == state)
        segment_length += step_size;
      else {
        cigar.emplace_back(segment_length, state);
        segment_length = static_cast<std::size_t>(step_size);
        state          = new_state;
      }
    } while (pos_i > 0 && pos_j > 0);

    cigar.emplace_back(segment_length, state);
    std::size_t alignment_offset = pos_i;
    if (pos_j > 0) cigar.emplace_back(pos_j, State::S);

    cigar.reverse();
    return {alignment_offset, cigar};
  }
};

}  // namespace hc
