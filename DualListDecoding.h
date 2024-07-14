#ifndef DUALLISTDECODER_H
#define DUALLISTDECODER_H

#include <map>
#include <queue>
#include <vector>
#include <iostream>
#include <fstream>
#include <climits>
#include <chrono>
#include <tuple>

#include "minHeap.h"
#include "stopWatch.h"
#include "viterbiDecoder.h"



class FeedForwardTrellis;

struct DLDInfo {
  double combined_metric;
  std::vector<int> message;
  std::vector<int> list_ranks;
  std::vector<double> received_signal;
  std::tuple<double, double, double> symbol_metrics;
};

class DualListMap{
  public:
    DualListMap() {};
    ~DualListMap() {};

    void insert(const MessageInformation& mi);
    int queue_size() {return agreed_messages_.size();};
    DLDInfo pop_queue();
    DLDInfo get_top();
    
    std::map<std::vector<int>, MessageInformation> dual_list_map_; // dictionary
    std::priority_queue<DLDInfo, std::vector<DLDInfo>, CompareCombinedMetric<DLDInfo>> agreed_messages_; // priority queue
};

namespace dual_list_decoder_utils {

std::vector<int> convertIntToBits(int integer, const int& length);
int hammingDistance(const std::vector<int> x, const std::vector<int>& y);
double euclideanDistance(const std::vector<double>& x,
                         const std::vector<int>& y);
std::vector<int> xOR(const std::vector<int>& x, const std::vector<int>& y);
template <typename T>
void print(const std::vector<T>& vec) {
  for (const T& element : vec) {
    std::cout << element << " ";
  }
  std::cout << std::endl;
}

template <typename T>
void print(const std::vector<std::vector<T>>& matrix) {
  for (const std::vector<T>& row : matrix) {
    for (const T& element : row) {
      std::cout << element << " ";
    }
    std::cout << std::endl;
  }
}

template <typename T>
void output(const std::vector<T>& vec, std::ofstream& outputFile) {
  for (const T& element : vec) {
    outputFile << element << " ";
  }
  // outputFile << std::endl;
}

template <typename T>
void output(const std::vector<std::vector<T>>& matrix,
            std::ofstream& outputFile) {
  for (const std::vector<T>& row : matrix) {
    for (const T& element : row) {
      outputFile << element << " ";
    }
    outputFile << std::endl;
  }
}

template <typename T>
void outputMat(const std::vector<T>& vec, std::ofstream& outputFile) {
  outputFile << "[";
  for (const T& ele : vec) {
    outputFile << ele << ", ";
  }
  outputFile << "] with size: " << vec.size();
}

template <typename T>
bool areVectorsEqual(const std::vector<T>& vector1,
                     const std::vector<T>& vector2) {
  if (vector1.size() != vector2.size()) {
    return false;  // Vectors have different sizes, so they cannot be equal.
  }

  for (size_t i = 0; i < vector1.size(); ++i) {
    if (vector1[i] != vector2[i]) {
      return false;  // Elements at index i are different.
    }
  }

  return true;  // All elements are equal.
}

// Define a custom comparison function for the priority queue
struct CompareCombinedMetric {
  bool operator()(const DLDInfo& a, const DLDInfo& b) const {
    return a.combined_metric >
           b.combined_metric;  // Lower combined_metric at the top
  }
};

// Function to combine two vectors of MessageInformation and create a priority
// queue of DLDInfo sorted in ascending order of combined_metric
std::priority_queue<DLDInfo, std::vector<DLDInfo>, CompareCombinedMetric>
combine_maps(const std::vector<MessageInformation>& vec1,
             const std::vector<MessageInformation>& vec2);
}


class DualListDecoder : private ViterbiDecoder{
 public:
  // DualListDecoder(std::vector<CodeInformation> code_info, int max_searched_path);
  DualListDecoder(CodeInformation encoder, std::vector<CodeInformation> code_info, int max_searched_path);
  ~DualListDecoder();
  
  // ----------------------------------------------  RATE 1/1 DLD DECODERS  --------------------------------------------------------
  DLDInfo adaptiveDecode(std::vector<double> received_signal, std::vector<std::chrono::milliseconds>& timeDurations);

  // Decoding function that alternates between two dual list decoders.
  // This function does not take crc degrees into consideration
  DLDInfo AdaptiveDecode_SimpleAlternate(
      std::vector<double> received_signal, std::vector<std::chrono::milliseconds>& timeDurations);
  
  // Decoding function that alternates between two dual list decoders considering crc degrees.
  // For example, if crc degree for two decoders are 3 and 5 respectively. Then every 8 and 32 codewords
  // in both lists there exist one codeword that passes crc (?).
  DLDInfo AdaptiveDecode_CRCAlternate(std::vector<double> received_signal, std::vector<std::chrono::milliseconds>& timeDurations);
  
  // -----------------------------------------  RATE 1/2 DLD DECODERS - WITH PUNCTURING  ---------------------------------------------

  // Decoding dunction that contains two rate 1/2 decoders. It alternates between two dual list decoders.
  // This function does not take crc degrees into consideration
  DLDInfo AdaptiveDecode_SimpleAlternate_rate_1_2(std::vector<double> received_signal, std::vector<std::chrono::milliseconds>& timeDurations);

  DLDInfo LookAheadDecode_SimpleAlternate_rate_1_2(
      std::vector<double> received_signal, std::vector<std::chrono::milliseconds>& timeDurations, std::vector<double> metric_0, std::vector<double> metric_1);

  // Look ahead decoding function that stop once match is found. If the match does not yield the best metric, then return the unmatched best
  DLDInfo LookAheadDecode_SimpleAlternate_StopOnceMatchFound(std::vector<double> received_signal, std::vector<std::chrono::milliseconds>& timeDurations);

  // Look ahead decoding function that stop once match is found. If the match does not yield the best metric, then return the unmatched best.
  // This function will also return the best unmatched metric if the match is not found and list size is exceeded.
  DLDInfo DLD_BAM(std::vector<double> txSig_1, std::vector<double> txSig_2, std::vector<std::chrono::milliseconds>& timeDurations);

  // Look ahead decoding function that stop once match is found. If the match does not yield the best metric, then return the unmatched best.
  // This function will also return the best unmatched metric if the match is not found and list size is exceeded.
  // This function will also only add half of the shared metric to the path metric.
  DLDInfo DLD_BAM_Half_Metric_on_Shared(std::vector<double> received_signal, std::vector<std::chrono::milliseconds>& timeDurations);
  
  
  DLDInfo DLD_BAM_double_match(std::vector<double> received_signal, std::vector<std::chrono::milliseconds>& timeDurations);

  DLDInfo DLD_Turbo(std::vector<double> txSig_1, std::vector<double> txSig_2);

  // ------------------------------------------------  TRACEBACK FUNCTIONS  ----------------------------------------------------------
  MessageInformation TraceBack_Single(minheap* heap, const CodeInformation& code, FeedForwardTrellis* trellis_ptr,
                             const std::vector<std::vector<Cell>>& trellis_states, std::vector<std::vector<int>>& prev_paths,
                             int& num_path_searched, int num_total_stages);
  
  // Trace back according to crc degree
  std::vector<MessageInformation> TraceBack_Multiple(
                            minheap* heap, const CodeInformation& code, FeedForwardTrellis* trellis_ptr,
                            const std::vector<std::vector<Cell>>& trellis_states,
                            std::vector<std::vector<int>>& prev_paths, int& num_path_searched,
                            int num_total_stages, int num_trace_back);

  
  // RATE 1/2 DLD Utilities
  std::vector<double> HardDecode(const std::vector<double>& received_signal);

 private:
  
  int crc_ratio_;

  //// ZTCC

  // Test function
  std::vector<std::vector<Cell>> constructZTListTrellis_precompute(
      const std::vector<double>& received_signal, CodeInformation code,
      FeedForwardTrellis* trellis_ptr, std::chrono::milliseconds& ssv_time);
  // Construct a ZTCC Trellis measuring the time taken by trellis construction (for both lists, iteratively)
  // Uses a regular euclidean metric.
  std::vector<std::vector<Cell>> ConstructZTCCTrellis_WithList_EuclideanMetric(
      const std::vector<double>& received_signal, CodeInformation code,
      FeedForwardTrellis* trellis_ptr, std::chrono::milliseconds& ssv_time);

  // Construct a ZTCC Trellis measuring the time taken by trellis construction (for both lists, iteratively)
  // Uses a special metric shown by Bill Ryan. Instead of calculating euclidean distance between received point and +/- 1,
  // we compute the product of these two values.
  std::vector<std::vector<Cell>> ConstructZTCCTrellis_WithList_ProductMetric(
      const std::vector<double>& received_signal, CodeInformation code,
      FeedForwardTrellis* trellis_ptr, std::chrono::milliseconds& ssv_time);

  // Construct a ZTCC Trellis measuring the time taken by trellis construction (for both lists, iteratively)
  // Uses a regular euclidean metric. This function only adds half of the shared metric to the path metric.
  std::vector<std::vector<Cell>> ConstructZTCCTrellis_WithList_EuclideanMetric_HalfSharedMetric(
      const std::vector<double>& received_signal, CodeInformation code,
      FeedForwardTrellis* trellis_ptr, std::chrono::milliseconds& ssv_time, const int& shared_index);

  //// TBCC
  // @todo
  // Construct a TBCC Trellis measuring the time taken by trellis construction (for both lists, iteratively)
  // Uses a regular euclidean metric.
  std::vector<std::vector<Cell>> ConstructTBCCTrellis_WithList_EuclideanMetric(
      const std::vector<double>& received_signal, CodeInformation code,
      FeedForwardTrellis* trellis_ptr, std::chrono::milliseconds& ssv_time);

  // @todo
  // Construct a TBCC Trellis measuring the time taken by trellis construction (for both lists, iteratively)
  // Uses a special metric shown by Bill Ryan. Instead of calculating euclidean distance between received point and +/- 1,
  // we compute the product of these two values.
  std::vector<std::vector<Cell>> ConstructTBCCTrellis_WithList_ProductMetric(
      const std::vector<double>& received_signal, CodeInformation code,
      FeedForwardTrellis* trellis_ptr, std::chrono::milliseconds& ssv_time);
  
  std::vector<int> convertPathtoMessage(
    const std::vector<int> path, FeedForwardTrellis* trellis_ptr);
  std::vector<int> convertPathtoTrimmedMessage(
    const std::vector<int> path, CodeInformation code, FeedForwardTrellis* trellis_ptr);
  std::vector<int> deconvolveCRC(const std::vector<int>& output, CodeInformation code);
  bool CRC_Check(std::vector<int> input_data, int crc_bits_num, int crc_dec);
};

#endif