#include "listdecoder.h"
#include <algorithm>
#include <climits>
#include <cmath>
#include <iostream>
#include <vector>

std::vector<std::vector<ListDecoder::cell>>
ListDecoder::constructOneTrellis(std::vector<double> receivedMessage) {

  pathLength = receivedMessage.size() + 1;
  std::vector<std::vector<cell>> trellisInfo;
  trellisInfo =
      std::vector<std::vector<cell>>(numStates, std::vector<cell>(pathLength));

  // initializes all the valid starting states
  for (int i = 0; i < numStates / 2; i++) {
    trellisInfo[i][0].pathMetric = 0;
    trellisInfo[i][0].init = true;
  }

  // precomputing euclidean distance between the received signal and +/- 1
  std::vector<std::vector<double>> precomputedMetrics;
  precomputedMetrics = std::vector<std::vector<double>>(receivedMessage.size(),
                                                        std::vector<double>(2));

  

  for (int stage = 0; stage < receivedMessage.size(); stage++) {
    // precomputedMetrics[stage][0] = std::abs(receivedMessage[stage] - 1);
    // precomputedMetrics[stage][1] = std::abs(receivedMessage[stage] + 1);
    precomputedMetrics[stage][0] = std::pow(receivedMessage[stage] - 1, 2);
    precomputedMetrics[stage][1] = std::pow(receivedMessage[stage] + 1, 2);
  }


  // building the trellis
  for (int stage = 0; stage < receivedMessage.size(); stage++) {
    for (int currentState = 0; currentState < numStates; currentState++) {
      // if the state / stage is invalid, we move on
      if (!trellisInfo[currentState][stage].init)
        continue;

      // otherwise, we compute the relevent information
      for (int forwardPathIndex = 0; forwardPathIndex < numForwardPaths;
           forwardPathIndex++) {
        // note that the forwardPathIndex is also the bit that corresponds with
        // the trellis transition

        int nextState = nextStates[currentState][stage % numTrellisSegLength]
                                  [forwardPathIndex];

        // if the nextState is invalid, we move on
        if (nextState < 0)
          continue;

        double totalPathMetric = precomputedMetrics[stage][forwardPathIndex] +
                                 trellisInfo[currentState][stage].pathMetric;

        // dealing with cases of uninitialized states, when the transition
        // becomes the optimal father state, and suboptimal father state, in
        // order
        if (!trellisInfo[nextState][stage + 1].init) {
          trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
          trellisInfo[nextState][stage + 1].optimalFatherState = currentState;
          trellisInfo[nextState][stage + 1].init = true;
        } else if (trellisInfo[nextState][stage + 1].pathMetric >
                   totalPathMetric) {
          trellisInfo[nextState][stage + 1].suboptimalPathMetric =
              trellisInfo[nextState][stage + 1].pathMetric;
          trellisInfo[nextState][stage + 1].suboptimalFatherState =
              trellisInfo[nextState][stage + 1].optimalFatherState;
          trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
          trellisInfo[nextState][stage + 1].optimalFatherState = currentState;
        } else {
          trellisInfo[nextState][stage + 1].suboptimalPathMetric =
              totalPathMetric;
          trellisInfo[nextState][stage + 1].suboptimalFatherState =
              currentState;
        }
      }
    }
  }

  return trellisInfo;
}

ListDecoder::messageInformation
ListDecoder::oneTrellisDecoding(std::vector<double> receivedMessage) {
  pathLength = receivedMessage.size() + 1;

  // trellisInfo is indexed [state][stage]
  std::vector<std::vector<cell>> trellisInfo;
  trellisInfo = constructOneTrellis(receivedMessage);

  // start search
  messageInformation output;
  // RBTree detourTree;
  minheap detourTree;
  std::vector<std::vector<int>> previousPaths;

  // create nodes for each valid ending state with no detours
  // std::cout<< "end path metrics:" <<std::endl;
  for (int i = 0; i < numStates / 2; i++) {
    DetourObject detour;
    detour.startingState = i;
    detour.pathMetric = trellisInfo[i][pathLength - 1].pathMetric;
    detourTree.insert(detour);
    // std::cout<< detour.pathMetric <<std::endl;
  }

  int numPathsSearched = 0;

  while (numPathsSearched < listSize) {
    DetourObject detour = detourTree.pop();
    std::vector<int> path(pathLength);

    int newTracebackStage = pathLength - 1;
    double forwardPartialPathMetric = 0;
    int currentState = detour.startingState;

    // if we are taking a detour from a previous path, we skip backwards to the
    // point where we take the detour from the previous path
    if (detour.originalPathIndex != -1) {
      forwardPartialPathMetric = detour.forwardPathMetric;
      newTracebackStage = detour.detourStage;

      // while we only need to copy the path from the detour to the end, this
      // simplifies things, and we'll write over the earlier data in any case
      path = previousPaths[detour.originalPathIndex];
      currentState = path[newTracebackStage];

      double suboptimalPathMetric =
          trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

      currentState =
          trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
      newTracebackStage--;

      double prevPathMetric =
          trellisInfo[currentState][newTracebackStage].pathMetric;

      forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
    }
    path[newTracebackStage] = currentState;

    // actually tracing back
    for (int stage = newTracebackStage; stage > 0; stage--) {
      double suboptimalPathMetric =
          trellisInfo[currentState][stage].suboptimalPathMetric;
      double currPathMetric = trellisInfo[currentState][stage].pathMetric;

      // if there is a detour we add to the detourTree
      if (trellisInfo[currentState][stage].suboptimalFatherState != -1) {
        DetourObject localDetour;
        localDetour.detourStage = stage;
        localDetour.originalPathIndex = numPathsSearched;
        localDetour.pathMetric =
            suboptimalPathMetric + forwardPartialPathMetric;
        localDetour.forwardPathMetric = forwardPartialPathMetric;
        localDetour.startingState = detour.startingState;
        detourTree.insert(localDetour);
      }
      currentState = trellisInfo[currentState][stage].optimalFatherState;
      double prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
      forwardPartialPathMetric += currPathMetric - prevPathMetric;
      path[stage - 1] = currentState;
    }
    previousPaths.push_back(path);

    std::vector<int> message = pathToMessage(path);
    // std::cout<< "decoded message: " << std::endl;
    // for (int i=0; i<message.size(); i++){
    // 	std::cout<< message[i];
    // }
    // std::cout<< "end of message" << std::endl;

    print_int_vector(message);
    // one trellis decoding requires both a tb and crc check
    if (path[0] == path[pathLength - 1] && crc_check(message, crcDegree, crc)) {
      output.message = message;
      output.path = path;
      output.listSize = numPathsSearched + 1;
      // return output;

      numPathsSearched++;
    }

    // numPathsSearched++;
  }
  output.listSizeExceeded = true;
  return output;
}
