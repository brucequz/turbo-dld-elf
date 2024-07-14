#include "crc_search_functions.h"

std::vector<CRC_pass_count_pair> getCRCs(int crcDegree) {
  std::vector<CRC_pass_count_pair> crcs;
  for (int i = 0; i < std::pow(2, crcDegree); i++) {
    // need to pad with zeros to make sure the crc is the right length
    // std::cout << std::bitset<4>(i) << std::endl;
    // make sure i has 1 in the first bit
    if (i % 2 == 0 || i >> (crcDegree - 1) == 0) {
      continue;
    }
    crcs.push_back(CRC_pass_count_pair{i, 0});
  }
  return crcs;
}