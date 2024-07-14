#ifndef VITERBIDECODER_H
#define VITERBIDECODER_H

#include <vector>
#include <tuple>

#include "viterbiCodec.h"

struct FeedForwardTrellis;



template<typename LDInfo>
struct CompareCombinedMetric {
  bool operator()(const LDInfo& a, const LDInfo& b) const {
    return a.combined_metric >
            b.combined_metric;  // Lower combined_metric at the top
  }
};

class ViterbiDecoder {
  public:
    ViterbiDecoder(CodeInformation encoder, std::vector<CodeInformation> code_info, int max_path_to_search) {
      encoder_ = encoder;
      max_path_to_search_ = max_path_to_search;
      for (auto& code : code_info) {
        code_info_.push_back(code);
      }
    }

    ~ViterbiDecoder(){};

  // Add any other common methods or member variables here

  protected:
    CodeInformation encoder_;
    FeedForwardTrellis* encoder_trellis_ptr_;
    std::vector<CodeInformation> code_info_;
    std::vector<FeedForwardTrellis*> trellis_ptrs_;
    int max_path_to_search_;
};

#endif // VITERBIDECODER_H