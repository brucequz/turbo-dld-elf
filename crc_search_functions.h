
#ifndef _CRC_SEARCH_FUNCTIONS_
#define _CRC_SEARCH_FUNCTIONS_

#include <vector>
#include <cmath>

struct CRC_pass_count_pair {
	int crc;
	int pass_count;
} typedef CRC_pass_count_pair;

std::vector<CRC_pass_count_pair> getCRCs(int crcDegree);

#endif