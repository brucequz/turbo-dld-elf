#ifndef VITERBICODEC_H
#define VITERBICODEC_H

#include <cassert>
#include <climits>
#include <vector>
#include <iostream>
#include <fstream>
#include <random>

#include "stopWatch.h"

struct CodeInformation {
  int k;              // input length
  int n;              // output length
  int v;              // memory elements
  int list_size = 1;  // list decoder list size
  int crc_dec;
  int crc_length;
  std::vector<int> generator_poly;
};

struct MessageInformation {
  MessageInformation(){};
  int decoder_index = -1;
  std::vector<int> message;
  std::vector<int> path;
  std::pair<int, int> begin_end_states;
  double path_metric;
  int list_rank;
  int crc_passing_rank;
  bool list_size_exceeded = false;  // added, used in DLD
  std::tuple<double, double, double> symbol_metrics;
};

struct Cell {
  bool init = false;
  double pathMetric = 3000;
  int fatherState = -1;
  double subPathMetric = 3000;
  int subFatherState = -1;
};

class FeedForwardTrellis;

namespace CodecUtils {
std::vector<int> convertIntToBits(int integer, const int& length);
int hammingDistance(const std::vector<int> x, const std::vector<int>& y);
double euclideanDistance(const std::vector<double>& x, const std::vector<int>& y);
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
void output(const std::vector<std::vector<T>>& matrix, std::ofstream& outputFile) {
  for (const std::vector<T>& row : matrix) {
      for (const T& element : row) {
          outputFile << element << " ";
      }
      outputFile << std::endl;
  }
}

template<typename T>
void outputMat(const std::vector<T>& vec, std::ofstream& outputFile) {
  outputFile << "[";
  for (const T& ele : vec) {
    outputFile << ele << ", ";
  }
  outputFile << "] with size: " << vec.size();
}

template<typename T>
void printMat(const std::vector<T>& vec) {
  std::cout << "[";
  for (size_t i = 0; i < vec.size() - 1; ++i) {
    std::cout << vec[i] << ", ";
  }
  std::cout << vec.back();
  std::cout << "] with size: " << vec.size();
}

template <typename T>
bool areVectorsEqual(const std::vector<T>& vector1, const std::vector<T>& vector2) {
  if (vector1.size() != vector2.size()) {
    std::cout << "unequal sizes" << std::endl;
    return false; // Vectors have different sizes, so they cannot be equal.
  }

  for (size_t i = 0; i < vector1.size(); ++i) {
      if (vector1[i] != vector2[i]) {
          return false; // Elements at index i are different.
      }
  }

  return true; // All elements are equal.
}

}  // namespace CodecUtils

namespace BPSK {
  
std::vector<int> modulate(std::vector<int> encoded_msg);
std::vector<int> demodulate(std::vector<double> received_signal);

} // namespace BPSK

namespace CRC {
  int binSum(const int& x, const int& y);
  std::vector<int> decToBin(int input, int bit_number);
}  // namespace CRC

class ViterbiCodec {
 public:
  ViterbiCodec(int k, int n, int v, std::vector<int> poly);
  ViterbiCodec(CodeInformation code);
  ~ViterbiCodec();
  std::vector<int> encode(const std::vector<int>& message);
  std::vector<int> encodeZTCC(std::vector<int> message);
  MessageInformation softViterbiDecoding(std::vector<double> receivedMessage, std::chrono::milliseconds& ssv_traceback_time);

  std::vector<int> calculateCRC(const std::vector<int>& input);
  std::vector<int> convolveCRC(const std::vector<int>& input);
  std::vector<int> deconvolveCRC(const std::vector<int>& output);

  // Function definition in listDecoder.cpp
  std::vector<MessageInformation> ZTCCListDecoding_fullInformation_WithConstraint(
    const std::vector<double>& received_signal);
  
  
  std::vector<MessageInformation> ZTCCListDecoding_fullInformation_NoConstraint(
    const std::vector<double>& received_signal);

 private:
  int k_;  // input message length
  int n_;  // output message length
  int v_;  // number of memory elements
  int crc_dec_; // crc poly in decimal representation
  int crc_length_; // length of crc in binary representation
  double code_rate_;
  int numStates_;
  int list_size_;
  FeedForwardTrellis* trellis_ptr_;
  
  // Trellis construction for softViterbiDecoding function only. The construction process will not keep track of
  // of the suboptimal path metric / father state
  std::vector<std::vector<Cell>> ConstructZTCCTrellis_NoList_EuclideanMetric(std::vector<double> receivedMessage);

  // Trellis construction function for decoders with list decoding functionality. It will keep track of the suboptimal
  // path metrics and suboptimal father state
  std::vector<std::vector<Cell>> ConstructZTCCTrellis_WithList_EuclideanMetric(std::vector<double> receivedMessage);

  std::vector<int> convertPathtoMessage(const std::vector<int> path);
  std::vector<int> convertPathtoTrimmedMessage(const std::vector<int> path);
  bool CRC_Check(std::vector<int> input_data, int crc_bits_num, int crc_dec);
};

#endif