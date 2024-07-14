#include "listdecoder.h"
#include <algorithm>
#include <climits>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
using namespace std;

std::vector<std::vector<std::vector<ListDecoder::cell>>>
ListDecoder::constructNTrellis(std::vector<double> receivedMessage) {
  // trellisInfo is indexed [starting/ending state][state][stage]
  pathLength = receivedMessage.size() + 1;
  std::vector<std::vector<std::vector<cell>>> trellisInfo;
  trellisInfo = std::vector<std::vector<std::vector<cell>>>(
      numStates / 2,
      std::vector<std::vector<cell>>(numStates, std::vector<cell>(pathLength)));

  // initializes all the valid starting states
  for (int i = 0; i < numStates / 2; i++) {
    trellisInfo[i][i][0].pathMetric = 0;
    trellisInfo[i][i][0].init = true;
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
  for (int startingState = 0; startingState < numStates / 2; startingState++) {
    for (int stage = 0; stage < receivedMessage.size(); stage++) {
      for (int currentState = 0; currentState < numStates; currentState++) {
        // if the state / stage is invalid, we move on
        if (!trellisInfo[startingState][currentState][stage].init)
          continue;

        // otherwise, we compute the relevent information
        for (int forwardPathIndex = 0; forwardPathIndex < numForwardPaths;
             forwardPathIndex++) {
          // note that the forwardPathIndex is also the bit that corresponds
          // with the trellis transition

          int nextState = nextStates[currentState][stage % numTrellisSegLength]
                                    [forwardPathIndex];

          // if the nextState is invalid, we move on
          if (nextState < 0)
            continue;

          double totalPathMetric =
              precomputedMetrics[stage][forwardPathIndex] +
              trellisInfo[startingState][currentState][stage].pathMetric;

          // dealing with cases of uninitialized states, when the transition
          // becomes the optimal father state, and suboptimal father state, in
          // order
          if (!trellisInfo[startingState][nextState][stage + 1].init) {
            trellisInfo[startingState][nextState][stage + 1].pathMetric =
                totalPathMetric;
            trellisInfo[startingState][nextState][stage + 1]
                .optimalFatherState = currentState;
            trellisInfo[startingState][nextState][stage + 1].init = true;
          } else if (trellisInfo[startingState][nextState][stage + 1]
                         .pathMetric > totalPathMetric) {
            trellisInfo[startingState][nextState][stage + 1]
                .suboptimalPathMetric =
                trellisInfo[startingState][nextState][stage + 1].pathMetric;
            trellisInfo[startingState][nextState][stage + 1]
                .suboptimalFatherState =
                trellisInfo[startingState][nextState][stage + 1]
                    .optimalFatherState;
            trellisInfo[startingState][nextState][stage + 1].pathMetric =
                totalPathMetric;
            trellisInfo[startingState][nextState][stage + 1]
                .optimalFatherState = currentState;
          } else {
            trellisInfo[startingState][nextState][stage + 1]
                .suboptimalPathMetric = totalPathMetric;
            trellisInfo[startingState][nextState][stage + 1]
                .suboptimalFatherState = currentState;
          }
        }
      }
    }
  }
  return trellisInfo;
}

ListDecoder::messageInformation
ListDecoder::nTrellisDecoding(std::vector<double> receivedMessage) {

  pathLength = receivedMessage.size() + 1;

  // trellisInfo is indexed [starting/ending state][state][stage]
  std::vector<std::vector<std::vector<cell>>> trellisInfo;
  trellisInfo = constructNTrellis(receivedMessage);

  // start search
  messageInformation output;
  // RBTree detourTree;
  minheap detourTree;

  std::vector<std::vector<int>> previousPaths;

  // create nodes for each valid ending state with no detours
  for (int i = 0; i < numStates / 2; i++) {
    DetourObject detour;
    detour.startingState = i;
    detour.pathMetric = trellisInfo[i][i][pathLength - 1].pathMetric;
    // queue.push(detour);
    detourTree.insert(detour);
  }

  int numPathsSearched = 0;

  while (numPathsSearched < listSize) {
    // DetourObject detour = queue.top();
    // queue.pop();
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
          trellisInfo[detour.startingState][currentState][newTracebackStage]
              .suboptimalPathMetric;

      currentState =
          trellisInfo[detour.startingState][currentState][newTracebackStage]
              .suboptimalFatherState;
      newTracebackStage--;

      double prevPathMetric =
          trellisInfo[detour.startingState][currentState][newTracebackStage]
              .pathMetric;

      forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
    }
    path[newTracebackStage] = currentState;

    // actually tracing back
    for (int stage = newTracebackStage; stage > 0; stage--) {
      double suboptimalPathMetric =
          trellisInfo[detour.startingState][currentState][stage]
              .suboptimalPathMetric;
      double currPathMetric =
          trellisInfo[detour.startingState][currentState][stage].pathMetric;

      // if there is a detour we add to the detourTree
      if (trellisInfo[detour.startingState][currentState][stage]
              .suboptimalFatherState != -1) {
        DetourObject localDetour;
        localDetour.detourStage = stage;
        localDetour.originalPathIndex = numPathsSearched;
        localDetour.pathMetric =
            suboptimalPathMetric + forwardPartialPathMetric;
        localDetour.forwardPathMetric = forwardPartialPathMetric;
        localDetour.startingState = detour.startingState;
        // queue.push(localDetour);
        detourTree.insert(localDetour);
      }
      currentState = trellisInfo[detour.startingState][currentState][stage]
                         .optimalFatherState;
      double prevPathMetric =
          trellisInfo[detour.startingState][currentState][stage - 1].pathMetric;
      forwardPartialPathMetric += currPathMetric - prevPathMetric;
      path[stage - 1] = currentState;
    }
    previousPaths.push_back(path);

    std::vector<int> message = pathToMessage(path);

    if (crc_check(message, crcDegree, crc)) {
      output.message = message;
      output.path = path;
      output.listSize = numPathsSearched + 1;
      return output;
    }

    numPathsSearched++;
  }
  output.listSizeExceeded = true;
  return output;
}

std::vector<int>
ListDecoder::ELFDistanceSpectra_highrate(std::vector<double> receivedMessage) {
  pathLength = receivedMessage.size() + 1;

  // trellisInfo is indexed [state][stage]
  std::vector<std::vector<std::vector<cell>>> trellisInfo;
  trellisInfo = constructNTrellis(receivedMessage);

  // start search
  messageInformation output;
  // RBTree detourTree;
  minheap detourTree;
  std::vector<std::vector<int>> previousPaths;

  std::vector<int> spectra;
  // create nodes for each valid ending state with no detours
  // std::cout<< "end path metrics:" <<std::endl;
  for (int i = 0; i < numStates / 2; i++) {
    DetourObject detour;
    detour.startingState = i;
    detour.pathMetric = trellisInfo[i][i][pathLength - 1].pathMetric;
    detourTree.insert(detour);
  }

  int numPathsSearched = 0;
  while (numPathsSearched < listSize) {
    DetourObject detour = detourTree.pop();
    std::vector<int> path(pathLength);

    int newTracebackStage = pathLength - 1;
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

    // one trellis decoding requires only crc check
    if (crc_check(message, crcDegree, crc) &&
        (path[0] == path[path.size() - 1])) {
      output.message = message;
      output.path = path;
      output.listSize = numPathsSearched + 1;
      output.metric = forwardPartialPathMetric;
      // std::cout << output.metric<< std::endl;
      spectra.push_back(output.metric);
      // if(output.metric == 20){
      // 	print_int_vector(message);
      // }
      // return output;
    }

    numPathsSearched++;
  }
  return spectra;
}

// this function is currently a dual list decoder with multiple trellises, needs
// to be updated
void ListDecoder::FindBestCRC_turbo(DualTrellis encTrell,
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

  pathLength = receivedMessage.size() + 1;

  // trellisInfo is indexed [starting/ending state][state][stage]
  std::vector<std::vector<std::vector<cell>>> trellisInfo;
  trellisInfo = constructNTrellis(receivedMessage);

  // RBTree detourTree;
  minheap detourTree;

  std::vector<std::vector<int>> previousPaths;

  // create nodes for each valid ending state with no detours
  for (int i = 0; i < numStates / 2; i++) {
    DetourObject detour;
    detour.startingState = i;
    detour.pathMetric = trellisInfo[i][i][pathLength - 1].pathMetric;
    // queue.push(detour);
    detourTree.insert(detour);
  }

  int numPathsSearched = 0;

  while (1) {
    // DetourObject detour = queue.top();
    // queue.pop();
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
          trellisInfo[detour.startingState][currentState][newTracebackStage]
              .suboptimalPathMetric;

      currentState =
          trellisInfo[detour.startingState][currentState][newTracebackStage]
              .suboptimalFatherState;
      newTracebackStage--;

      double prevPathMetric =
          trellisInfo[detour.startingState][currentState][newTracebackStage]
              .pathMetric;
      forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
    }
    path[newTracebackStage] = currentState;

    // actually tracing back
    for (int stage = newTracebackStage; stage > 0; stage--) {
      double suboptimalPathMetric =
          trellisInfo[detour.startingState][currentState][stage]
              .suboptimalPathMetric;
      double currPathMetric =
          trellisInfo[detour.startingState][currentState][stage].pathMetric;

      // if there is a detour we add to the detourTree
      if (trellisInfo[detour.startingState][currentState][stage]
              .suboptimalFatherState != -1) {
        DetourObject localDetour;
        localDetour.detourStage = stage;
        localDetour.originalPathIndex = numPathsSearched;
        localDetour.pathMetric =
            suboptimalPathMetric + forwardPartialPathMetric;
        localDetour.forwardPathMetric = forwardPartialPathMetric;
        localDetour.startingState = detour.startingState;
        // queue.push(localDetour);
        detourTree.insert(localDetour);
      }
      currentState = trellisInfo[detour.startingState][currentState][stage]
                         .optimalFatherState;
      double prevPathMetric =
          trellisInfo[detour.startingState][currentState][stage - 1].pathMetric;
      forwardPartialPathMetric += currPathMetric - prevPathMetric;
      path[stage - 1] = currentState;
    }
    previousPaths.push_back(path);

    std::vector<int> message = pathToMessage(path);
    /*todo: add the interleaved message here*/

    numPathsSearched++;
    paths_at_distance++;
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
        // if (crcs[i].pass_count == 1) {
        // 	vector<int> codeword = encTrell.encoder(message);
        // 	cout << "Codeword: ";
        // 	for (int j = 0; j < codeword.size(); j++) {
        // 		cout << codeword[j];
        // 	}
        // 	cout << endl;
        // }
      }

      /*todo: add crc_check for interleaved message here*/
    }
  }
  return;
}

ListDecoder::messageInformation ListDecoder::punctured_highrate_turbo_decoder(
    std::vector<double> receivedMessage, std::vector<int> punc_idx) {

  pathLength = receivedMessage.size() + 1;

  // trellisInfo is indexed [starting/ending state][state][stage]
  std::vector<std::vector<std::vector<cell>>> trellisInfo;
  trellisInfo = constructNTrellis(receivedMessage);

  // start search
  messageInformation output;
  // RBTree detourTree;
  minheap detourTree;

  std::vector<std::vector<int>> previousPaths;

  // create nodes for each valid ending state with no detours
  for (int i = 0; i < numStates / 2; i++) {
    DetourObject detour;
    detour.startingState = i;
    detour.pathMetric = trellisInfo[i][i][pathLength - 1].pathMetric;
    // queue.push(detour);
    detourTree.insert(detour);
  }

  int numPathsSearched = 0;

  while (numPathsSearched < listSize) {
    // DetourObject detour = queue.top();
    // queue.pop();
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
          trellisInfo[detour.startingState][currentState][newTracebackStage]
              .suboptimalPathMetric;

      currentState =
          trellisInfo[detour.startingState][currentState][newTracebackStage]
              .suboptimalFatherState;
      newTracebackStage--;

      double prevPathMetric =
          trellisInfo[detour.startingState][currentState][newTracebackStage]
              .pathMetric;

      forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
    }
    path[newTracebackStage] = currentState;

    // actually tracing back
    for (int stage = newTracebackStage; stage > 0; stage--) {
      double suboptimalPathMetric =
          trellisInfo[detour.startingState][currentState][stage]
              .suboptimalPathMetric;
      double currPathMetric =
          trellisInfo[detour.startingState][currentState][stage].pathMetric;

      // if there is a detour we add to the detourTree
      if (trellisInfo[detour.startingState][currentState][stage]
              .suboptimalFatherState != -1) {
        DetourObject localDetour;
        localDetour.detourStage = stage;
        localDetour.originalPathIndex = numPathsSearched;
        localDetour.pathMetric =
            suboptimalPathMetric + forwardPartialPathMetric;
        localDetour.forwardPathMetric = forwardPartialPathMetric;
        localDetour.startingState = detour.startingState;
        // queue.push(localDetour);
        detourTree.insert(localDetour);
      }
      currentState = trellisInfo[detour.startingState][currentState][stage]
                         .optimalFatherState;
      double prevPathMetric =
          trellisInfo[detour.startingState][currentState][stage - 1].pathMetric;
      forwardPartialPathMetric += currPathMetric - prevPathMetric;
      path[stage - 1] = currentState;
    }
    previousPaths.push_back(path);

    std::vector<int> message = pathToMessage(path);

    if (crc_check(message, crcDegree, crc)) {
      output.message = message;
      output.path = path;
      output.listSize = numPathsSearched + 1;
      return output;
    }

    numPathsSearched++;
  }
  output.listSizeExceeded = true;
  return output;
}