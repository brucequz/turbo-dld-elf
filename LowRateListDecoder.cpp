#include "lowratelistdecoder.h"
#include <algorithm>
#include <climits>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
using namespace std;

LowRateListDecoder::LowRateListDecoder(FeedforwardTrellis feedforwardTrellis,
                                       int listSize, int crcDegree, int crc) {
  this->lowrate_nextStates = feedforwardTrellis.getNextStates();
  this->lowrate_outputs = feedforwardTrellis.getOutputs();
  this->lowrate_numStates = feedforwardTrellis.getNumStates();
  this->lowrate_symbolLength = feedforwardTrellis.getN();
  this->numForwardPaths = lowrate_nextStates[0].size();
  // this->lowrate_pathLength = lowrate_nextStates;
  this->listSize = listSize;
  this->crcDegree = crcDegree;
  this->crc = crc;

  int v = feedforwardTrellis.getV();
}

LowRateListDecoder::LowRateListDecoder(
    FeedforwardTrellis feedforwardTrellis, int listSize, int crcDegree, int crc,
    std::vector<std::vector<int>> neighboring_list,
    std::vector<std::vector<int>> neighboring_msg,
    std::vector<std::vector<int>> path_ie_state) {
  this->lowrate_nextStates = feedforwardTrellis.getNextStates();
  this->lowrate_outputs = feedforwardTrellis.getOutputs();
  this->lowrate_numStates = feedforwardTrellis.getNumStates();
  this->lowrate_symbolLength = feedforwardTrellis.getN();
  this->numForwardPaths = lowrate_nextStates[0].size();
  // this->lowrate_pathLength = lowrate_nextStates;
  this->listSize = listSize;
  this->crcDegree = crcDegree;
  this->crc = crc;
  this->neighboring_list = neighboring_list;
  this->neighboring_msg = neighboring_msg;
  this->path_ie_state = path_ie_state;

  int v = feedforwardTrellis.getV();
}

std::vector<std::vector<LowRateListDecoder::cell>>
LowRateListDecoder::constructLowRateTrellis(
    std::vector<double> receivedMessage) {
  std::vector<std::vector<cell>> trellisInfo;
  lowrate_pathLength = (receivedMessage.size() / lowrate_symbolLength) + 1;

  trellisInfo = std::vector<std::vector<cell>>(
      lowrate_numStates, std::vector<cell>(lowrate_pathLength));

  // initializes all the valid starting states
  for (int i = 0; i < lowrate_numStates; i++) {
    trellisInfo[i][0].pathMetric = 0;
    trellisInfo[i][0].init = true;
  }

  // building the trellis
  for (int stage = 0; stage < lowrate_pathLength - 1; stage++) {
    for (int currentState = 0; currentState < lowrate_numStates;
         currentState++) {
      // if the state / stage is invalid, we move on
      if (!trellisInfo[currentState][stage].init)
        continue;

      // otherwise, we compute the relevent information
      for (int forwardPathIndex = 0; forwardPathIndex < numForwardPaths;
           forwardPathIndex++) {
        // since our transitions correspond to symbols, the forwardPathIndex has
        // no correlation beyond indexing the forward path

        int nextState = lowrate_nextStates[currentState][forwardPathIndex];

        // if the nextState is invalid, we move on
        if (nextState < 0)
          continue;

        double branchMetric = 0;
        std::vector<int> output_point =
            get_point(lowrate_outputs[currentState][forwardPathIndex],
                      lowrate_symbolLength);

        for (int i = 0; i < lowrate_symbolLength; i++) {
          branchMetric +=
              std::pow(receivedMessage[lowrate_symbolLength * stage + i] -
                           (double)output_point[i],
                       2);
          // branchMetric += std::abs(receivedMessage[lowrate_symbolLength *
          // stage + i] - (double)output_point[i]);
        }
        double totalPathMetric =
            branchMetric + trellisInfo[currentState][stage].pathMetric;

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

std::vector<std::vector<std::vector<LowRateListDecoder::cell>>>
LowRateListDecoder::constructLowRateMultiTrellis(
    std::vector<double> receivedMessage) {
  // trellisInfo is indexed [TB state][state][stage]
  std::vector<std::vector<std::vector<cell>>> trellisInfo;
  lowrate_pathLength = (receivedMessage.size() / lowrate_symbolLength) + 1;

  trellisInfo = std::vector<std::vector<std::vector<cell>>>(
      lowrate_numStates,
      std::vector<std::vector<cell>>(lowrate_numStates,
                                     std::vector<cell>(lowrate_pathLength)));
  //(numStates/2, std::vector<std::vector<cell>>(numStates,
  //std::vector<cell>(pathLength)))

  // initializes all the valid starting states
  for (int i = 0; i < lowrate_numStates; i++) {
    trellisInfo[i][i][0].pathMetric = 0;
    trellisInfo[i][i][0].init = true;
  }

  // building the trellis
  for (int TBState = 0; TBState < lowrate_numStates; TBState++) {
    for (int stage = 0; stage < lowrate_pathLength - 1; stage++) {
      for (int currentState = 0; currentState < lowrate_numStates;
           currentState++) {
        // if the state / stage is invalid, we move on
        if (!trellisInfo[TBState][currentState][stage].init)
          continue;

        // otherwise, we compute the relevent information
        for (int forwardPathIndex = 0; forwardPathIndex < numForwardPaths;
             forwardPathIndex++) {
          // since our transitions correspond to symbols, the forwardPathIndex
          // has no correlation beyond indexing the forward path

          int nextState = lowrate_nextStates[currentState][forwardPathIndex];

          // if the nextState is invalid, we move on
          if (nextState < 0)
            continue;

          double branchMetric = 0;
          std::vector<int> output_point =
              get_point(lowrate_outputs[currentState][forwardPathIndex],
                        lowrate_symbolLength);

          for (int i = 0; i < lowrate_symbolLength; i++) {
            branchMetric +=
                std::pow(receivedMessage[lowrate_symbolLength * stage + i] -
                             (double)output_point[i],
                         2);
            // branchMetric += std::abs(receivedMessage[lowrate_symbolLength *
            // stage + i] - (double)output_point[i]);
          }
          double totalPathMetric =
              branchMetric +
              trellisInfo[TBState][currentState][stage].pathMetric;

          // dealing with cases of uninitialized states, when the transition
          // becomes the optimal father state, and suboptimal father state, in
          // order
          if (!trellisInfo[TBState][nextState][stage + 1].init) {
            trellisInfo[TBState][nextState][stage + 1].pathMetric =
                totalPathMetric;
            trellisInfo[TBState][nextState][stage + 1].optimalFatherState =
                currentState;
            trellisInfo[TBState][nextState][stage + 1].init = true;
          } else if (trellisInfo[TBState][nextState][stage + 1].pathMetric >
                     totalPathMetric) {
            trellisInfo[TBState][nextState][stage + 1].suboptimalPathMetric =
                trellisInfo[TBState][nextState][stage + 1].pathMetric;
            trellisInfo[TBState][nextState][stage + 1].suboptimalFatherState =
                trellisInfo[TBState][nextState][stage + 1].optimalFatherState;
            trellisInfo[TBState][nextState][stage + 1].pathMetric =
                totalPathMetric;
            trellisInfo[TBState][nextState][stage + 1].optimalFatherState =
                currentState;
          } else {
            trellisInfo[TBState][nextState][stage + 1].suboptimalPathMetric =
                totalPathMetric;
            trellisInfo[TBState][nextState][stage + 1].suboptimalFatherState =
                currentState;
          }
        }
      }
    }
  }
  return trellisInfo;
}

LowRateListDecoder::messageInformation
LowRateListDecoder::lowRateDecoding(std::vector<double> receivedMessage) {
  // trellisInfo is indexed [state][stage]
  std::vector<std::vector<cell>> trellisInfo;
  trellisInfo = constructLowRateTrellis(receivedMessage);

  // start search
  messageInformation output;
  // RBTree detourTree;
  minheap detourTree;
  std::vector<std::vector<int>> previousPaths;

  // create nodes for each valid ending state with no detours
  // std::cout<< "end path metrics:" <<std::endl;
  for (int i = 0; i < lowrate_numStates; i++) {
    DetourObject detour;
    detour.startingState = i;
    detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
    detourTree.insert(detour);
  }

  int numPathsSearched = 0;
  while (numPathsSearched < this->listSize) {
    DetourObject detour = detourTree.pop();
    std::vector<int> path(lowrate_pathLength);

    int newTracebackStage = lowrate_pathLength - 1;
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

    // one trellis decoding requires both a tb and crc check
    if (path[0] == path[lowrate_pathLength - 1] &&
        crc_check(message, crcDegree, crc)) {
      output.message = message;
      output.path = path;
      output.listSize = numPathsSearched + 1;
      output.metric = forwardPartialPathMetric;
      return output;
    }

    numPathsSearched++;
  }
  output.listSizeExceeded = true;
  return output;
}

LowRateListDecoder::messageInformation
LowRateListDecoder::multiTrellisLowRateDecoding(
    std::vector<double> receivedMessage) {
  // trellisInfo is indexed [state][stage]
  std::vector<std::vector<std::vector<cell>>> trellisInfo;
  trellisInfo = constructLowRateMultiTrellis(receivedMessage);

  // start search
  messageInformation output;
  // RBTree detourTree;
  minheap detourTree;
  std::vector<std::vector<int>> previousPaths;

  // create nodes for each valid ending state with no detours
  // std::cout<< "end path metrics:" <<std::endl;
  for (int i = 0; i < lowrate_numStates; i++) {
    DetourObject detour;
    detour.startingState = i;
    detour.pathMetric = trellisInfo[i][i][lowrate_pathLength - 1].pathMetric;
    detourTree.insert(detour);
  }

  std::ofstream outputFile_cwd;
  std::string Cwdfilename = "multiTrellis_cwd_K64v8m8.txt";
  outputFile_cwd.open(Cwdfilename, std::fstream::app); // have appending

  int numPathsSearched = 0;
  while (numPathsSearched < listSize) {
    DetourObject detour = detourTree.pop();
    std::vector<int> path(lowrate_pathLength);

    int newTracebackStage = lowrate_pathLength - 1;
    double forwardPartialPathMetric = 0;
    int currentState = detour.startingState;
    int TBState = detour.startingState;

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
          trellisInfo[TBState][currentState][newTracebackStage]
              .suboptimalPathMetric;

      currentState = trellisInfo[TBState][currentState][newTracebackStage]
                         .suboptimalFatherState;
      newTracebackStage--;

      double prevPathMetric =
          trellisInfo[TBState][currentState][newTracebackStage].pathMetric;

      forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
    }
    path[newTracebackStage] = currentState;

    // actually tracing back
    for (int stage = newTracebackStage; stage > 0; stage--) {
      double suboptimalPathMetric =
          trellisInfo[TBState][currentState][stage].suboptimalPathMetric;
      double currPathMetric =
          trellisInfo[TBState][currentState][stage].pathMetric;

      // if there is a detour we add to the detourTree
      if (trellisInfo[TBState][currentState][stage].suboptimalFatherState !=
          -1) {
        DetourObject localDetour;
        localDetour.detourStage = stage;
        localDetour.originalPathIndex = numPathsSearched;
        localDetour.pathMetric =
            suboptimalPathMetric + forwardPartialPathMetric;
        localDetour.forwardPathMetric = forwardPartialPathMetric;
        localDetour.startingState = detour.startingState;
        detourTree.insert(localDetour);
      }
      currentState =
          trellisInfo[TBState][currentState][stage].optimalFatherState;
      double prevPathMetric =
          trellisInfo[TBState][currentState][stage - 1].pathMetric;
      forwardPartialPathMetric += currPathMetric - prevPathMetric;
      path[stage - 1] = currentState;
    }
    previousPaths.push_back(path);

    std::vector<int> message = pathToMessage(path);
    std::vector<int> codeword = pathToCodeword(path);
    for (int j = 0; j < codeword.size(); j++) {
      outputFile_cwd << std::to_string(codeword[j]) + " ";
    }
    outputFile_cwd << std::endl;

    // one trellis decoding requires both a tb and crc check
    if (crc_check(message, crcDegree, crc)) {
      output.message = message;
      output.path = path;
      output.listSize = numPathsSearched + 1;
      output.metric = forwardPartialPathMetric;
      outputFile_cwd << "/////////////" << std::endl;
      outputFile_cwd.close();
      return output;
    }

    numPathsSearched++;
  }
  output.listSizeExceeded = true;
  outputFile_cwd << "/////////////" << std::endl;
  outputFile_cwd.close();
  return output;
}

LowRateListDecoder::messageInformation
LowRateListDecoder::offsetLowRateDecoding(std::vector<double> receivedMessage) {
  double weight = 0;
  for (int i = 0; i < receivedMessage.size(); i++)
    weight += std::pow(receivedMessage[i], 2);
  for (int i = 0; i < receivedMessage.size(); i++)
    receivedMessage[i] /= weight;
  // trellisInfo is indexed [state][stage]
  std::vector<std::vector<cell>> trellisInfo;
  trellisInfo = constructLowRateTrellis(receivedMessage);

  // first, we perform a single SSV pass, then use that codeword as
  // the input of the standard low rate decoder
  std::vector<int> path(lowrate_pathLength);

  int newTracebackStage = lowrate_pathLength - 1;
  double forwardPartialPathMetric = 0;
  int currentState = -1;

  double minMetric = INT_MAX;
  for (int i = 0; i < lowrate_numStates; i++) {
    if (trellisInfo[i][lowrate_pathLength - 1].pathMetric < minMetric) {
      currentState = i;
      minMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
    }
  }

  path[newTracebackStage] = currentState;

  // actually tracing back
  for (int stage = newTracebackStage; stage > 0; stage--) {
    double suboptimalPathMetric =
        trellisInfo[currentState][stage].suboptimalPathMetric;
    double currPathMetric = trellisInfo[currentState][stage].pathMetric;

    currentState = trellisInfo[currentState][stage].optimalFatherState;
    double prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
    // std::cout << "prev path metric (1 transition): " << prevPathMetric <<
    // std::endl;
    forwardPartialPathMetric += currPathMetric - prevPathMetric;
    path[stage - 1] = currentState;
  }

  std::vector<int> offsetCodeword = pathToCodeword(path);
  /*std::cout << "path metric: " << forwardPartialPathMetric << std::endl;
  double offsetWeight = 0;
  for(int i = 0; i < offsetCodeword.size(); i++){
          offsetWeight += std::pow(offsetCodeword[i] - receivedMessage[i], 2);
  }
  std::cout << "offset codeword metric: " << offsetWeight << std::endl;*/
  std::vector<double> doubleOffset(offsetCodeword.begin(),
                                   offsetCodeword.end());

  return lowRateDecoding(doubleOffset);
}

LowRateListDecoder::messageInformation
LowRateListDecoder::offsetMultiTrellisLowRateDecoding(
    std::vector<double> receivedMessage) {
  // trellisInfo is indexed [state][stage]
  std::vector<std::vector<std::vector<cell>>> trellisInfo;
  trellisInfo = constructLowRateMultiTrellis(receivedMessage);

  // first, we perform a single SSV pass, then use that codeword as
  // the input of the standard low rate decoder
  std::vector<int> path(lowrate_pathLength);

  messageInformation output;

  int newTracebackStage = lowrate_pathLength - 1;
  double forwardPartialPathMetric = 0;
  int currentState = -1;

  double minMetric = INT_MAX;
  for (int i = 0; i < lowrate_numStates; i++) {
    if (trellisInfo[i][i][lowrate_pathLength - 1].pathMetric < minMetric) {
      currentState = i;
      minMetric = trellisInfo[i][i][lowrate_pathLength - 1].pathMetric;
    }
  }

  path[newTracebackStage] = currentState;

  int TBState = currentState;

  // actually tracing back
  for (int stage = newTracebackStage; stage > 0; stage--) {
    double suboptimalPathMetric =
        trellisInfo[TBState][currentState][stage].suboptimalPathMetric;
    double currPathMetric =
        trellisInfo[TBState][currentState][stage].pathMetric;

    currentState = trellisInfo[TBState][currentState][stage].optimalFatherState;
    double prevPathMetric =
        trellisInfo[TBState][currentState][stage - 1].pathMetric;
    forwardPartialPathMetric += currPathMetric - prevPathMetric;
    path[stage - 1] = currentState;
  }

  std::vector<int> message = pathToMessage(path);

  // one trellis decoding requires both a tb and crc check
  if (crc_check(message, crcDegree, crc)) {
    output.message = message;
    output.path = path;
    output.listSize = 1;
    output.metric = forwardPartialPathMetric;
    return output;
  }

  std::vector<int> offsetCodeword = pathToCodeword(path);
  std::vector<double> doubleOffset(offsetCodeword.begin(),
                                   offsetCodeword.end());

  return lowRateDecoding(doubleOffset);
}

std::vector<int>
LowRateListDecoder::ELFDistanceSpectra(std::vector<double> receivedMessage) {
  // trellisInfo is indexed [state][stage]
  std::vector<std::vector<std::vector<cell>>> trellisInfo;
  trellisInfo = constructLowRateMultiTrellis(receivedMessage);

  // start search
  messageInformation output;
  // RBTree detourTree;
  minheap detourTree;
  std::vector<std::vector<int>> previousPaths;

  std::vector<int> spectra;
  // create nodes for each valid ending state with no detours
  // std::cout<< "end path metrics:" <<std::endl;
  for (int i = 0; i < lowrate_numStates; i++) {
    DetourObject detour;
    detour.startingState = i;
    detour.pathMetric = trellisInfo[i][i][lowrate_pathLength - 1].pathMetric;
    detourTree.insert(detour);
  }

  int numPathsSearched = 0;
  while (numPathsSearched < listSize) {
    DetourObject detour = detourTree.pop();
    std::vector<int> path(lowrate_pathLength);

    int newTracebackStage = lowrate_pathLength - 1;
    double forwardPartialPathMetric = 0;
    int currentState = detour.startingState;
    int TBState = detour.startingState;

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
          trellisInfo[TBState][currentState][newTracebackStage]
              .suboptimalPathMetric;

      currentState = trellisInfo[TBState][currentState][newTracebackStage]
                         .suboptimalFatherState;
      newTracebackStage--;

      double prevPathMetric =
          trellisInfo[TBState][currentState][newTracebackStage].pathMetric;

      forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
    }
    path[newTracebackStage] = currentState;

    // actually tracing back
    for (int stage = newTracebackStage; stage > 0; stage--) {
      double suboptimalPathMetric =
          trellisInfo[TBState][currentState][stage].suboptimalPathMetric;
      double currPathMetric =
          trellisInfo[TBState][currentState][stage].pathMetric;

      // if there is a detour we add to the detourTree
      if (trellisInfo[TBState][currentState][stage].suboptimalFatherState !=
          -1) {
        DetourObject localDetour;
        localDetour.detourStage = stage;
        localDetour.originalPathIndex = numPathsSearched;
        localDetour.pathMetric =
            suboptimalPathMetric + forwardPartialPathMetric;
        localDetour.forwardPathMetric = forwardPartialPathMetric;
        localDetour.startingState = detour.startingState;
        detourTree.insert(localDetour);
      }
      currentState =
          trellisInfo[TBState][currentState][stage].optimalFatherState;
      double prevPathMetric =
          trellisInfo[TBState][currentState][stage - 1].pathMetric;
      forwardPartialPathMetric += currPathMetric - prevPathMetric;
      path[stage - 1] = currentState;
    }
    previousPaths.push_back(path);

    std::vector<int> message = pathToMessage(path);

    // one trellis decoding requires both a tb and crc check
    if (crc_check(message, crcDegree, crc)) {
      output.message = message;
      output.path = path;
      output.listSize = numPathsSearched + 1;
      output.metric = forwardPartialPathMetric;
      spectra.push_back(output.metric);
      if (forwardPartialPathMetric == 14 * 4) {
        print_int_vector(message);
      }
      // return output;
    }

    numPathsSearched++;
  }
  return spectra;
}

void LowRateListDecoder::Find_ELFs_TB_Improved(
    FeedforwardTrellis encTrell, std::vector<double> receivedMessage,
    codeInformation code) {

  std::vector<CRC_pass_count_pair> crcs = getCRCs(code.crcDeg);
  int distance = 0;
  int previous_distance = 0;
  int tb_paths_at_distance = 0;
  int best_non_failing_distance = 0;
  vector<int> high_performers = {};

  std::ofstream outputFile;
  std::stringstream filename;
  filename << "TBCC-v=" << code.v << "(" << code.numerators[0] << ", "
           << code.numerators[1] << ")"
           << "-m=" << code.crcDeg - 1 << "K=" << code.numInfoBits << ".txt";
  outputFile.open(filename.str(), fstream::app);

  // trellisInfo is indexed [state][stage]
  std::vector<std::vector<std::vector<cell>>> trellisInfo;
  trellisInfo = constructLowRateMultiTrellis(receivedMessage);

  // RBTree detourTree;
  minheap detourTree;
  std::vector<std::vector<int>> previousPaths;

  // create nodes for each valid ending state with no detours
  // std::cout<< "end path metrics:" <<std::endl;
  for (int i = 0; i < lowrate_numStates; i++) {
    DetourObject detour;
    detour.startingState = i;
    detour.pathMetric = trellisInfo[i][i][lowrate_pathLength - 1].pathMetric;
    detourTree.insert(detour);
  }

  int numPathsSearched = 0;
  while (1) {
    DetourObject detour = detourTree.pop();
    std::vector<int> path(lowrate_pathLength);

    int newTracebackStage = lowrate_pathLength - 1;
    double forwardPartialPathMetric = 0;
    int currentState = detour.startingState;
    int TBState = detour.startingState;

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
          trellisInfo[TBState][currentState][newTracebackStage]
              .suboptimalPathMetric;

      currentState = trellisInfo[TBState][currentState][newTracebackStage]
                         .suboptimalFatherState;
      newTracebackStage--;

      double prevPathMetric =
          trellisInfo[TBState][currentState][newTracebackStage].pathMetric;

      forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
    }
    path[newTracebackStage] = currentState;

    // actually tracing back
    for (int stage = newTracebackStage; stage > 0; stage--) {
      double suboptimalPathMetric =
          trellisInfo[TBState][currentState][stage].suboptimalPathMetric;
      double currPathMetric =
          trellisInfo[TBState][currentState][stage].pathMetric;

      // if there is a detour we add to the detourTree
      if (trellisInfo[TBState][currentState][stage].suboptimalFatherState !=
          -1) {
        DetourObject localDetour;
        localDetour.detourStage = stage;
        localDetour.originalPathIndex = numPathsSearched;
        localDetour.pathMetric =
            suboptimalPathMetric + forwardPartialPathMetric;
        localDetour.forwardPathMetric = forwardPartialPathMetric;
        localDetour.startingState = detour.startingState;
        detourTree.insert(localDetour);
      }
      currentState =
          trellisInfo[TBState][currentState][stage].optimalFatherState;
      double prevPathMetric =
          trellisInfo[TBState][currentState][stage - 1].pathMetric;
      forwardPartialPathMetric += currPathMetric - prevPathMetric;
      path[stage - 1] = currentState;
    }
    previousPaths.push_back(path);

    std::vector<int> message = pathToMessage(path);

    tb_paths_at_distance++;
    // sum up the message bits

    // message should be K + m bits long
    int new_distance = forwardPartialPathMetric;

    // check if we have a new distance for the next codeword
    if (new_distance > distance) {

      previous_distance = distance;
      distance = new_distance;

      outputFile << "\n"
                 << "distance: " << previous_distance / 4 << endl;
      outputFile << "tail-biting paths_at_distance: " << tb_paths_at_distance
                 << endl;

      for (CRC_pass_count_pair crc : crcs) {
        outputFile << "CRC: " << hex << crc.crc << dec
                   << ", pass count: " << crc.pass_count << endl;
      }

      // Find the best performing CRC at the last distance
      int lowest_pass_count = crcs[0].pass_count;
      for (CRC_pass_count_pair crc : crcs) {
        if (crc.pass_count < lowest_pass_count) {
          lowest_pass_count = crc.pass_count;
        }
      }
      // check if the lowest pass count is 0. If so, prune away all the non-zero
      // pass count CRCs if all passed, then record all of these CRCs and then
      // delete the CRCs that don't have the lowest pass count
      if (lowest_pass_count == 0) {
        vector<CRC_pass_count_pair> passed_crcs;
        for (auto crc : crcs) {
          if (crc.pass_count == 0) {
            passed_crcs.push_back(crc);
          }
        }
        crcs = passed_crcs;
      } else {
        if (high_performers.size() == 0) {
          best_non_failing_distance = previous_distance / 4;
          // copy over the current CRC vector into high_performers
          for (int i = 0; i < crcs.size(); i++) {
            high_performers.push_back(crcs[i].crc);
          }
        }
        vector<CRC_pass_count_pair> passed_crcs;
        for (auto crc : crcs) {
          if (crc.pass_count == lowest_pass_count) {
            passed_crcs.push_back(crc);
          }
        }
        crcs = passed_crcs;
      }

      outputFile << "crcs left: " << crcs.size() << endl;

      // cout << "CRCs left: ";
      // for (CRC_pass_count_pair crc: crcs) {
      // 	cout << hex << ", crc: " << crc.crc << dec << ", pass count: "
      // << dec << crc.pass_count;
      // }
      // cout << "\n" << endl;

      if (crcs.size() == 1 && high_performers.size() != 0) {
        outputFile << "\n"
                   << "Best distance: " << best_non_failing_distance << endl;
        // return high performers as well as the best CRC
        // print out high_performers and best CRC
        outputFile << "High performers: " << hex;
        for (int i = 0; i < high_performers.size(); i++) {
          outputFile << high_performers[i] << " ";
        }
        outputFile << ", Best CRC: " << crcs[0].crc << ", pass count: " << dec
                   << crcs[0].pass_count << endl;

        outputFile.close();
        return;
      }
      // reset the pass count for each CRC
      for (int i = 0; i < crcs.size(); i++) {
        crcs[i].pass_count = 0;
      }
      tb_paths_at_distance = 0;
    }
    // if (numPathsSearched <= 1) {
    // 	cout << "Distance: " << previous_distance << endl;
    // 	continue;
    // }
    // for each codeword, check if the CRC passes or not.
    // If not, then increment the pass count for that CRC
    for (int i = 0; i < crcs.size(); i++) {
      if (crc_check(message, code.crcDeg, crcs[i].crc)) {
        crcs[i].pass_count++;
        // if (distance == 0) {
        // 	vector<int> codeword = encTrell.encoder(message);
        // 	cout << "Codeword: ";
        // 	for (int j = 0; j < codeword.size(); j++) {
        // 		cout << codeword[j];
        // 	}
        // 	cout << endl;
        // }
      }
    }
  }
  return;
}

void LowRateListDecoder::FindBestCRC_TBCC(FeedforwardTrellis encTrell,
                                          std::vector<double> receivedMessage,
                                          codeInformation code) {

  std::vector<CRC_pass_count_pair> crcs = getCRCs(code.crcDeg);
  int distance = 0;
  int previous_distance = 0;
  int paths_at_distance = 0;
  int tb_paths_at_distance = 0;
  int best_non_failing_distance = 0;
  vector<int> high_performers = {};

  std::ofstream outputFile;
  std::stringstream filename;
  filename << "TBCC-v=" << code.v << "(" << code.numerators[0] << ", "
           << code.numerators[1] << ")"
           << "-m=" << code.crcDeg - 1 << "K=" << code.numInfoBits << ".txt";
  outputFile.open(filename.str(), fstream::app);

  std::vector<std::vector<cell>> trellisInfo;
  trellisInfo = constructLowRateTrellis(receivedMessage);

  // start search
  messageInformation output;
  // RBTree detourTree;
  minheap detourTree;
  std::vector<std::vector<int>> previousPaths;

  // create nodes for each valid ending state with no detours
  // std::cout<< "end path metrics:" <<std::endl;
  for (int i = 0; i < lowrate_numStates; i++) {
    DetourObject detour;
    detour.startingState = i;
    detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
    detourTree.insert(detour);
  }

  int numPathsSearched = 0;
  while (1) {
    DetourObject detour = detourTree.pop();
    std::vector<int> path(lowrate_pathLength);

    int newTracebackStage = lowrate_pathLength - 1;
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

    // one trellis decoding requires both a tb and crc check

    numPathsSearched++;
    paths_at_distance++;

    // tail-biting check
    if (path[0] != path[lowrate_pathLength - 1])
      continue;

    tb_paths_at_distance++;
    // sum up the message bits

    // message should be K + m bits long
    int new_distance = forwardPartialPathMetric;

    // check if we have a new distance for the next codeword
    if (new_distance > distance) {

      previous_distance = distance;
      distance = new_distance;

      outputFile << "\n"
                 << "distance: " << previous_distance / 4 << endl;
      outputFile << "paths_at_distance: " << paths_at_distance << endl;
      outputFile << "tail-biting paths_at_distance: " << tb_paths_at_distance
                 << endl;

      for (CRC_pass_count_pair crc : crcs) {
        outputFile << "CRC: " << hex << crc.crc << dec
                   << ", pass count: " << crc.pass_count << endl;
      }

      // Find the best performing CRC at the last distance
      int lowest_pass_count = crcs[0].pass_count;
      for (CRC_pass_count_pair crc : crcs) {
        if (crc.pass_count < lowest_pass_count) {
          lowest_pass_count = crc.pass_count;
        }
      }
      // check if the lowest pass count is 0. If so, prune away all the non-zero
      // pass count CRCs if all passed, then record all of these CRCs and then
      // delete the CRCs that don't have the lowest pass count
      if (lowest_pass_count == 0) {
        vector<CRC_pass_count_pair> passed_crcs;
        for (auto crc : crcs) {
          if (crc.pass_count == 0) {
            passed_crcs.push_back(crc);
          }
        }
        crcs = passed_crcs;
      } else {
        if (high_performers.size() == 0) {
          best_non_failing_distance = previous_distance / 4;
          // copy over the current CRC vector into high_performers
          for (int i = 0; i < crcs.size(); i++) {
            high_performers.push_back(crcs[i].crc);
          }
        }
        vector<CRC_pass_count_pair> passed_crcs;
        for (auto crc : crcs) {
          if (crc.pass_count == lowest_pass_count) {
            passed_crcs.push_back(crc);
          }
        }
        crcs = passed_crcs;
      }

      outputFile << "crcs left: " << crcs.size() << endl;

      // cout << "CRCs left: ";
      // for (CRC_pass_count_pair crc: crcs) {
      // 	cout << hex << ", crc: " << crc.crc << dec << ", pass count: "
      // << dec << crc.pass_count;
      // }
      // cout << "\n" << endl;

      if (crcs.size() == 1 && high_performers.size() != 0) {
        outputFile << "\n"
                   << "Best distance: " << best_non_failing_distance << endl;
        // return high performers as well as the best CRC
        // print out high_performers and best CRC
        outputFile << "High performers: " << hex;
        for (int i = 0; i < high_performers.size(); i++) {
          outputFile << high_performers[i] << " ";
        }
        outputFile << ", Best CRC: " << crcs[0].crc << ", pass count: " << dec
                   << crcs[0].pass_count << endl;

        outputFile.close();
        return;
      }
      // reset the pass count for each CRC
      for (int i = 0; i < crcs.size(); i++) {
        crcs[i].pass_count = 0;
      }
      tb_paths_at_distance = paths_at_distance = 0;
    }
    if (numPathsSearched <= 1) {
      cout << "Distance: " << previous_distance << endl;
      continue;
    }
    // for each codeword, check if the CRC passes or not.
    // If not, then increment the pass count for that CRC
    for (int i = 0; i < crcs.size(); i++) {
      if (crc_check(message, code.crcDeg, crcs[i].crc)) {
        crcs[i].pass_count++;
        if (distance == 0) {
          vector<int> codeword = encTrell.encoder(message);
          cout << "Codeword: ";
          for (int j = 0; j < codeword.size(); j++) {
            cout << codeword[j];
          }
          cout << endl;
        }
      }
    }
  }

  outputFile.close();
  return;
}

LowRateListDecoder::messageInformation
LowRateListDecoder::linearityLowRateListSearching(
    std::vector<double> receivedMessage) {
  // trellisInfo is indexed [state][stage]
  std::vector<std::vector<cell>> trellisInfo;
  trellisInfo = constructLowRateTrellis(receivedMessage);

  // start search
  messageInformation output(this);

  // RBTree detourTree;
  minheap detourTree;
  std::vector<std::vector<int>> previousPaths;

  // create nodes for each valid ending state with no detours
  // std::cout<< "end path metrics:" <<std::endl;
  for (int i = 0; i < lowrate_numStates; i++) {
    DetourObject detour;
    detour.startingState = i;
    detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
    detourTree.insert(detour);
  }

  int numPathsSearched = 0;
  while (numPathsSearched < (this->listSize)) {
    DetourObject detour = detourTree.pop();
    std::vector<int> path(lowrate_pathLength);

    int newTracebackStage = lowrate_pathLength - 1;
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
    std::vector<int> codeword = pathToCodeword(path);

    for (int j = 0; j < output.neighbors[0].size(); j++) {
      output.neighbors[numPathsSearched][j] = codeword[j];
    }
    for (int k = 0; k < output.neighbor_msgs[0].size(); k++) {
      output.neighbor_msgs[numPathsSearched][k] = message[k];
    }
    // print_int_vector(output.neighbors[numPathsSearched]);

    output.path_ie[numPathsSearched].push_back(path[0]);
    output.path_ie[numPathsSearched].push_back(path[lowrate_pathLength - 1]);
    output.message = message;
    output.path = path;
    output.listSize = numPathsSearched;
    output.metric = forwardPartialPathMetric;

    numPathsSearched++;
  }
  output.listSizeExceeded = true;
  return output;
}

LowRateListDecoder::messageInformation
LowRateListDecoder::linearityLowRateDecoding(
    std::vector<double> receivedMessage) {
  // trellisInfo is indexed [state][stage]
  std::vector<std::vector<cell>> trellisInfo;
  trellisInfo = constructLowRateTrellis(receivedMessage);

  // start search
  messageInformation output(this);

  // RBTree detourTree;
  minheap detourTree;
  std::vector<std::vector<int>> previousPaths;

  // create nodes for each valid ending state with no detours
  // std::cout<< "end path metrics:" <<std::endl;
  for (int i = 0; i < lowrate_numStates; i++) {
    DetourObject detour;
    detour.startingState = i;
    detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
    detourTree.insert(detour);
  }

  int numPathsSearched = 0;
  while (numPathsSearched < (this->listSize)) {
    DetourObject detour = detourTree.pop();
    std::vector<int> path(lowrate_pathLength);

    int newTracebackStage = lowrate_pathLength - 1;
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
    std::vector<int> codeword = pathToCodeword(path);

    if (path[0] == path[lowrate_pathLength - 1] &&
        crc_check(message, crcDegree, crc)) {
      // std::cout << "trivially correct - linear" << std::endl;
      output.message = message;
      output.path = path;
      output.listSize = numPathsSearched + 1;
      output.metric = forwardPartialPathMetric;
      return output;
    }

    // TB and CRC check
    std::vector<std::vector<int>> passingCodewords;
    std::vector<std::vector<int>> passingMessages;
    std::vector<std::vector<int>> xor_neighbor_cwds(
        this->neighboring_list.size(),
        std::vector<int>(this->neighboring_list[0].size(), 0));
    std::vector<std::vector<int>> xor_neighbor_msgs(
        this->neighboring_msg.size(),
        std::vector<int>(this->neighboring_msg[0].size(), 0));
    std::vector<std::vector<int>> xor_path_ie(
        this->path_ie_state.size(),
        std::vector<int>(this->path_ie_state[0].size(), 0));

    for (int i = 0; i < this->path_ie_state.size(); i++) {
      // std::cout << "neighbor message path_ie: " << this->path_ie_state[i][0]
      // << " , " << this->path_ie_state[i][1] << std::endl;
      xor_path_ie[i][0] = this->path_ie_state[i][0] ^ path[0];
      xor_path_ie[i][1] =
          this->path_ie_state[i][1] ^ path[lowrate_pathLength - 1];
      // std::cout << "updated neighbor message path_ie: " <<
      // this->path_ie_state[i][0] << " , " << this->path_ie_state[i][1] <<
      // std::endl;
      if (xor_path_ie[i][0] == xor_path_ie[i][1]) {
        // std::cout << "We are at: " << i << "th neighbor" << std::endl;
        //  XOR the codewords
        for (int j = 0; j < this->neighboring_list[0].size(); j++) {
          xor_neighbor_cwds[i][j] = this->neighboring_list[i][j] ^ codeword[j];
        }
        for (int k = 0; k < this->neighboring_msg[0].size(); k++) {
          xor_neighbor_msgs[i][k] = this->neighboring_msg[i][k] ^ message[k];
        }
        // codeword to message
        // messageInformation decoding;
        // std::vector<double> temp(this->neighboring_list[i].size());
        // for (int k=0; k<temp.size(); k++){
        // 	temp[k] = static_cast<double>(this->neighboring_list[i][k]);
        // }
        // decoding = this->lowRateDecoding(temp);

        if (crc_check(xor_neighbor_msgs[i], crcDegree, crc)) {
          // std::cout << "We are at: " << i << "th neighbor" << std::endl;
          passingCodewords.push_back(xor_neighbor_cwds[i]);
          passingMessages.push_back(xor_neighbor_msgs[i]);
        }

      } else {
        continue;
      }
    }

    // std::cout << "  ++++++++++++ Passing Messages" << std::endl;
    // for (int i=0; i < passingMessages.size(); i++){
    // 	print_int_vector(passingMessages[i]);
    // }
    // Euclidean distance with received codeword of all passing codewords
    std::vector<double> euclidean_distance(passingCodewords.size(), 0);
    for (int i = 0; i < passingCodewords.size(); i++) {
      for (int j = 0; j < passingCodewords[i].size(); j++) {
        euclidean_distance[i] +=
            std::pow(receivedMessage[j] - (double)passingCodewords[i][j], 2);
      }
    }

    // Find the index of the smallest hamming distance and use that to find the
    // best message
    auto smallestIt =
        std::min_element(euclidean_distance.begin(), euclidean_distance.end());
    int index = std::distance(euclidean_distance.begin(), smallestIt);

    // Assign output.message value
    if (passingMessages.size() != 0) {
      output.message = passingMessages[index];
      output.listSize = -3;
      return output;
    }
    output.listSizeExceeded = true;
    output.listSize = -4;
    return output;
  }
  output.listSizeExceeded = true;
  return output;
}

LowRateListDecoder::messageInformation
LowRateListDecoder::TBGuaranteedListSearching(
    std::vector<double> receivedMessage) {
  // trellisInfo is indexed [state][stage]
  std::vector<std::vector<cell>> trellisInfo;
  trellisInfo = constructLowRateTrellis(receivedMessage);

  // start search
  messageInformation output(this);

  // RBTree detourTree;
  minheap detourTree;
  std::vector<std::vector<int>> previousPaths;

  // create nodes for each valid ending state with no detours
  // std::cout<< "end path metrics:" <<std::endl;
  for (int i = 0; i < lowrate_numStates; i++) {
    DetourObject detour;
    detour.startingState = i;
    detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
    detourTree.insert(detour);
  }

  int numPathsSearched = 0;
  int curListSize = 0;
  while (curListSize < (this->listSize)) {
    DetourObject detour = detourTree.pop();
    std::vector<int> path(lowrate_pathLength);

    int newTracebackStage = lowrate_pathLength - 1;
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
    std::vector<int> codeword = pathToCodeword(path);

    // TB Guarantee neighbors finder
    if (path[0] == path[lowrate_pathLength - 1]) {

      for (int j = 0; j < output.neighbors[0].size(); j++) {
        output.neighbors[curListSize][j] = codeword[j];
      }
      for (int k = 0; k < output.neighbor_msgs[0].size(); k++) {
        output.neighbor_msgs[curListSize][k] = message[k];
      }
      output.path_ie[curListSize].push_back(path[0]);
      output.path_ie[curListSize].push_back(path[lowrate_pathLength - 1]);
      // print_int_vector(output.neighbor_msgs[curListSize]);
      // std::cout << std::endl;
      // std::cout << "Path check: ";
      // print_int_vector(output.path_ie[curListSize]);
      // std::cout << std::endl;

      output.message = message;
      output.path = path;
      output.listSize = numPathsSearched;
      output.metric = forwardPartialPathMetric;

      curListSize++;
    }
    numPathsSearched++;
  }
  output.listSizeExceeded = true;
  return output;
}

LowRateListDecoder::messageInformation
LowRateListDecoder::TBGuaranteedLinearityDecoding(
    std::vector<double> receivedMessage) {
  // trellisInfo is indexed [state][stage]
  std::vector<std::vector<cell>> trellisInfo;
  trellisInfo = constructLowRateTrellis(receivedMessage);

  // start search
  messageInformation output(this);

  // RBTree detourTree;
  minheap detourTree;
  std::vector<std::vector<int>> previousPaths;

  // create nodes for each valid ending state with no detours
  // std::cout<< "end path metrics:" <<std::endl;
  for (int i = 0; i < lowrate_numStates; i++) {
    DetourObject detour;
    detour.startingState = i;
    detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
    detourTree.insert(detour);
  }

  int numPathsSearched = 0;
  while (numPathsSearched < 512) {
    DetourObject detour = detourTree.pop();
    std::vector<int> path(lowrate_pathLength);

    int newTracebackStage = lowrate_pathLength - 1;
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
    std::vector<int> codeword = pathToCodeword(path);

    // denoise process

    // if the closest TB-pass codeword also passes crc, we are done!
    // if (path[0] == path[lowrate_pathLength - 1] && crc_check(message,
    // crcDegree, crc)){
    // 	// std::cout << "trivially correct - linear" << std::endl;
    // 	output.message = message;
    // 	output.path = path;
    // 	output.listSize = numPathsSearched + 1;
    // 	output.metric = forwardPartialPathMetric;
    // 	return output;
    // }

    if (path[0] == path[lowrate_pathLength - 1]) {
      // std::cout << "Fouuuuund Ya!" << std::endl;
      std::vector<std::vector<int>> passingCodewords;
      std::vector<std::vector<int>> passingMessages;
      std::vector<std::vector<int>> xor_neighbor_cwds(
          this->neighboring_list.size(),
          std::vector<int>(this->neighboring_list[0].size(), 0));
      std::vector<std::vector<int>> xor_neighbor_msgs(
          this->neighboring_msg.size(),
          std::vector<int>(this->neighboring_msg[0].size(), 0));

      for (int i = 0; i < this->neighboring_list.size(); i++) {
        // XOR the codewords
        for (int j = 0; j < this->neighboring_list[0].size(); j++) {
          xor_neighbor_cwds[i][j] = this->neighboring_list[i][j] ^ codeword[j];
        }
        for (int k = 0; k < this->neighboring_msg[0].size(); k++) {
          xor_neighbor_msgs[i][k] = this->neighboring_msg[i][k] ^ message[k];
        }

        if (crc_check(xor_neighbor_msgs[i], crcDegree, crc)) {
          // std::cout << "We are at: " << i << "th neighbor" << std::endl;
          passingCodewords.push_back(xor_neighbor_cwds[i]);
          passingMessages.push_back(xor_neighbor_msgs[i]);
        }
      }
      std::vector<double> euclidean_distance(passingCodewords.size(), 0);
      for (int i = 0; i < passingCodewords.size(); i++) {
        for (int j = 0; j < passingCodewords[i].size(); j++) {
          euclidean_distance[i] +=
              std::pow(receivedMessage[j] - (double)passingCodewords[i][j], 2);
        }
      }

      // Find the index of the smallest hamming distance and use that to find
      // the best message
      auto smallestIt = std::min_element(euclidean_distance.begin(),
                                         euclidean_distance.end());
      int index = std::distance(euclidean_distance.begin(), smallestIt);

      // Assign output.message value
      // std::cout << "Passing message length: " << passingMessages.size() <<
      // std::endl;
      if (passingMessages.size() != 0) {
        output.message = passingMessages[index];
        // std::cout << "decoding success!" << std::endl;
        output.listSize = -3;
        return output;
      }
      output.listSizeExceeded = true;
      output.listSize = -4;
      return output;
    }
    numPathsSearched++;
  }
  output.listSizeExceeded = true;
  return output;
}

LowRateListDecoder::messageInformation
LowRateListDecoder::lowRateDecoding_listSize(
    std::vector<double> receivedMessage) {
  // trellisInfo is indexed [state][stage]
  std::vector<std::vector<cell>> trellisInfo;
  trellisInfo = constructLowRateTrellis(receivedMessage);

  // start search
  messageInformation output;
  // RBTree detourTree;
  minheap detourTree;
  std::vector<std::vector<int>> previousPaths;

  // create nodes for each valid ending state with no detours
  // std::cout<< "end path metrics:" <<std::endl;
  for (int i = 0; i < lowrate_numStates; i++) {
    DetourObject detour;
    detour.startingState = i;
    detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
    detourTree.insert(detour);
  }

  int numPathsSearched = 0;

  std::ofstream outputFile;
  std::string filename = "codewords.txt";
  outputFile.open(filename);
  while (numPathsSearched < listSize) {
    DetourObject detour = detourTree.pop();
    std::vector<int> path(lowrate_pathLength);

    int newTracebackStage = lowrate_pathLength - 1;
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

    // std::cout << "path metric at " << numPathsSearched << ": " <<
    // forwardPartialPathMetric << std::endl;

    std::vector<int> codeword = pathToCodeword(path);
    std::vector<double> justPastMidpointCodeword;

    for (int i = 0; i < codeword.size(); i++) {
      justPastMidpointCodeword.push_back(0.49 * 1 + 0.51 * codeword[i]);
    }

    justPastMidpointCodeword[47] = 0;
    justPastMidpointCodeword[60] = 0;
    justPastMidpointCodeword[129] = 0;
    justPastMidpointCodeword[504] = 0;

    // std::cout << "List size is: "<< numPathsSearched << ", Distance metric
    // is: " << forwardPartialPathMetric << std::endl;
    // print_int_vector(codeword);

    for (int i = 0; i < justPastMidpointCodeword.size() - 1; i++) {
      outputFile << justPastMidpointCodeword[i] << " ";
    }
    outputFile << justPastMidpointCodeword[justPastMidpointCodeword.size() - 1];
    outputFile << std::endl;

    numPathsSearched++;
  }
  outputFile.close();
  return output;
}

LowRateListDecoder::messageInformation
LowRateListDecoder::lowRateDecoding_recip(std::vector<double> receivedMessage) {
  // trellisInfo is indexed [state][stage]
  std::vector<std::vector<cell>> trellisInfo;
  trellisInfo = constructLowRateTrellis(receivedMessage);

  // start search
  messageInformation output;
  // RBTree detourTree;
  minheap detourTree;
  std::vector<std::vector<int>> previousPaths;

  // create nodes for each valid ending state with no detours
  // std::cout<< "end path metrics:" <<std::endl;
  for (int i = 0; i < lowrate_numStates; i++) {
    DetourObject detour;
    detour.startingState = i;
    detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
    detourTree.insert(detour);
  }

  int numPathsSearched = 0;
  while (numPathsSearched < listSize) {
    DetourObject detour = detourTree.pop();
    std::vector<int> path(lowrate_pathLength);

    int newTracebackStage = lowrate_pathLength - 1;
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
    std::vector<int> codeword = pathToCodeword(path);

    // if we reach the all-zero codeword
    std::vector<int> original;
    for (int i = 0; i < codeword.size(); i++) {
      if ((i != 47) && (i != 60) && (i != 129) && (i != 504)) {
        original.push_back(1);
      } else {
        original.push_back(0);
      }
    }
    if (codeword == original) {
      output.message = message;
      output.path = path;
      output.listSize = numPathsSearched + 1;
      output.metric = forwardPartialPathMetric;
      return output;
    }

    numPathsSearched++;
  }
  // output.listSizeExceeded = true;
  return output;
}

LowRateListDecoder::messageInformation
LowRateListDecoder::lowRateReciprocityCheck(std::vector<double> originalPoint,
                                            std::vector<double> distantPoint) {

  int localListSize = listSize;
  bool reciprocalDecoding = false;

  // trellisInfo is indexed [state][stage]
  std::vector<std::vector<cell>> trellisInfo;

  // in this case, we are taking the distant point and decoding until we return
  // to the original point
  if (distantPoint.size() > 0) {
    trellisInfo = constructLowRateTrellis(distantPoint);
    localListSize = 1e5;
    reciprocalDecoding = true;
  } else {
    trellisInfo = constructLowRateTrellis(originalPoint);
  }

  // start search
  messageInformation output;
  // RBTree detourTree;
  minheap detourTree;
  std::vector<std::vector<int>> previousPaths;

  // create nodes for each valid ending state with no detours
  for (int i = 0; i < lowrate_numStates; i++) {
    DetourObject detour;
    detour.startingState = i;
    detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
    detourTree.insert(detour);
  }

  int numPathsSearched = 0;
  while (numPathsSearched < localListSize) {
    DetourObject detour = detourTree.pop();
    std::vector<int> path(lowrate_pathLength);

    int newTracebackStage = lowrate_pathLength - 1;
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

    std::vector<int> codeword = pathToCodeword(path);
    std::vector<double> punc_msg;
    for (int i = 0; i < codeword.size(); i++) {
      if ((i != 47) && (i != 60) && (i != 129) && (i != 504)) {
        punc_msg.push_back(codeword[i]);
      } else {
        punc_msg.push_back(0);
      }
    }

    std::vector<int> message = pathToMessage(path);

    numPathsSearched++;

    if (reciprocalDecoding) {
      if (originalPoint == punc_msg) {
        output.message = message;
        output.listSize = numPathsSearched;
        return output;
      }
    } else {
      if (numPathsSearched == listSize) {
        return lowRateReciprocityCheck(originalPoint, punc_msg);
      }
    }
  }
  output.listSizeExceeded = true;
  return output;
}

std::vector<double> LowRateListDecoder::lowRateDistanceDistribution(
    std::vector<double> receivedMessage) {
  // trellisInfo is indexed [state][stage]
  std::vector<std::vector<cell>> trellisInfo;
  trellisInfo = constructLowRateTrellis(receivedMessage);

  // start search
  messageInformation output;
  // RBTree detourTree;
  minheap detourTree;
  std::vector<std::vector<int>> previousPaths;
  std::vector<double> distances;

  // create nodes for each valid ending state with no detours
  // std::cout<< "end path metrics:" <<std::endl;
  for (int i = 0; i < lowrate_numStates; i++) {
    DetourObject detour;
    detour.startingState = i;
    detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
    detourTree.insert(detour);
  }

  int numPathsSearched = 0;
  while (numPathsSearched < listSize) {
    DetourObject detour = detourTree.pop();
    std::vector<int> path(lowrate_pathLength);

    int newTracebackStage = lowrate_pathLength - 1;
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

    numPathsSearched++;

    distances.push_back(forwardPartialPathMetric);
  }
  return distances;
}

void LowRateListDecoder::printRangeOfCodewords(
    std::vector<double> receivedMessage, int start, int end) {
  // trellisInfo is indexed [state][stage]
  std::vector<std::vector<cell>> trellisInfo;
  trellisInfo = constructLowRateTrellis(receivedMessage);

  // start search
  messageInformation output;
  // RBTree detourTree;
  minheap detourTree;
  std::vector<std::vector<int>> previousPaths;

  // create nodes for each valid ending state with no detours
  // std::cout<< "end path metrics:" <<std::endl;
  for (int i = 0; i < lowrate_numStates; i++) {
    DetourObject detour;
    detour.startingState = i;
    detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
    detourTree.insert(detour);
  }

  int numPathsSearched = 0;
  while (numPathsSearched < listSize) {
    DetourObject detour = detourTree.pop();
    std::vector<int> path(lowrate_pathLength);

    int newTracebackStage = lowrate_pathLength - 1;
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

    numPathsSearched++;

    if (start < numPathsSearched && end > numPathsSearched) {
      std::cout << "on the " << numPathsSearched << "th path: " << std::endl;
      print_int_vector(pathToCodeword(path));
    }
  }
}

// converts a path through the tb trellis to the binary message it corresponds
// with
std::vector<int> LowRateListDecoder::pathToMessage(std::vector<int> path) {
  std::vector<int> message;
  for (int pathIndex = 0; pathIndex < path.size() - 1; pathIndex++) {
    for (int forwardPath = 0; forwardPath < numForwardPaths; forwardPath++) {
      if (lowrate_nextStates[path[pathIndex]][forwardPath] ==
          path[pathIndex + 1])
        message.push_back(forwardPath);
    }
  }
  return message;
}

// converts a path through the tb trellis to the BPSK it corresponds with
std::vector<int> LowRateListDecoder::pathToCodeword(std::vector<int> path) {
  std::vector<int> nopunc_codeword;
  for (int pathIndex = 0; pathIndex < path.size() - 1; pathIndex++) {
    for (int forwardPath = 0; forwardPath < numForwardPaths; forwardPath++) {
      if (lowrate_nextStates[path[pathIndex]][forwardPath] ==
          path[pathIndex + 1]) {
        std::vector<int> output_bin;
        dec_to_binary(lowrate_outputs[path[pathIndex]][forwardPath], output_bin,
                      lowrate_symbolLength);
        for (int outbit = 0; outbit < lowrate_symbolLength; outbit++) {
          nopunc_codeword.push_back(-2 * output_bin[outbit] + 1);
        }
      }
    }
  }
  // puncture 4 bits
  std::vector<int> codeword;
  for (int i = 0; i < nopunc_codeword.size(); i++) {
    if ((i != 47) && (i != 60) && (i != 129) && (i != 504)) {
      codeword.push_back(nopunc_codeword[i]);
    } else {
      codeword.push_back(0);
    }
  }
  return codeword;
}