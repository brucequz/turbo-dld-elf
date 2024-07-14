#include "listdecoder.h"
// #include "math.h"
#include <algorithm>
#include <climits>
#include <cmath>
#include <iostream>
#include <vector>

ListDecoder::ListDecoder(DualTrellis dualTrellis, int listSize, int crcDegree,
                         int crc) {
  this->nextStates = dualTrellis.getsubStates();
  this->numStates = nextStates.size();
  this->numTrellisSegLength = nextStates[0].size();
  this->numForwardPaths = nextStates[0][0].size();
  this->listSize = listSize;
  this->crcDegree = crcDegree;
  this->crc = crc;

  int n = dualTrellis.getN();
  int v = dualTrellis.getV();
  this->ztTrailingBits = (n - 1) * std::ceil((double)v / n); // n*ceil(v/(n-1))
}

// converts a path through the tb trellis to the binary message it corresponds
// with
std::vector<int> ListDecoder::pathToMessage(std::vector<int> path) {
  std::vector<int> message;
  for (int pathIndex = 0; pathIndex < pathLength - 1; pathIndex++) {
    if (pathIndex % numTrellisSegLength == 0)
      continue;
    for (int forwardPath = 0; forwardPath < numForwardPaths; forwardPath++) {
      if (nextStates[path[pathIndex]][pathIndex % numTrellisSegLength]
                    [forwardPath] == path[pathIndex + 1])
        message.push_back(forwardPath);
    }
  }
  return message;
}

// converts a path through the ztcc trellis to the binary message it corresponds
// with
std::vector<int> ListDecoder::ztPathToMessage(std::vector<int> path) {
  std::vector<int> message;
  for (int pathIndex = 0; pathIndex < pathLength - 2 - ztTrailingBits;
       pathIndex++) {
    if (pathIndex % numTrellisSegLength == 0)
      continue;
    for (int forwardPath = 0; forwardPath < numForwardPaths; forwardPath++) {
      if (nextStates[path[pathIndex]][pathIndex % numTrellisSegLength]
                    [forwardPath] == path[pathIndex + 1])
        message.push_back(forwardPath);
    }
  }
  return message;
}