#include "listdecoder.h"
#include <cassert>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>

#include "listdecoder.h"

namespace{

std::vector<double> ComputeSquaredDifferences(
    const std::vector<int>& vector1, const std::vector<double>& vector2) {
  // Check if the vectors have the same size
  if (vector1.size() != vector2.size()) {
    // You can handle this error in your preferred way, e.g., throw an exception
    throw std::invalid_argument("Vectors must have the same size");
  }

  // Calculate squared differences
  std::vector<double> squaredDifferences;
  squaredDifferences.reserve(vector1.size());  // Reserve space for efficiency

  for (std::size_t i = 0; i < vector1.size(); ++i) {
    double diff = vector1[i] - vector2[i];
    squaredDifferences.push_back(diff * diff);
  }

  return squaredDifferences;
}

}

DualListDecoder::DualListDecoder(std::vector<codeInformation> code_info, int listSize) {
  // Constructor
  // Input:
  //       - DT: DualTrellis object
  //       - listSize: list size
  //       - crcDegree: degree of the CRC
  //       - crc: CRC polynomial

  for (auto& code : code_info) {
    DualTrellis dt(code.hMatrix);
    ListDecoder ld(dt, listSize, code.crcDeg, code.crc);
    list_decoders_.push_back(ld);
    FeedbackTrellis* ptr = new FeedbackTrellis(code.k, code.n, code.v, code.numerators, code.denominator);
    trellis_ptrs_.push_back(ptr);
  }
}

DualListDecoder::~DualListDecoder() {
  // Destructor
  for (auto& ptr : trellis_ptrs_) {
    delete ptr;
  }
}

void DualListDecoder::DualListMap::insert(const ListDecoder::messageInformation& mi) {
  auto it = dual_list_map_.find(mi.message); // finding the match in the dictionary
  if (it != dual_list_map_.end()) {
    DLDInfo agreed_message;
    agreed_message.combined_metric = mi.metric;
    if (it->second.decoder_index == 0) {
      agreed_message.list_ranks = {it->second.listRank, mi.listRank};
    } else if (it->second.decoder_index == 1) {
      agreed_message.list_ranks = {mi.listRank, it->second.listRank};
    } else {
      std::cerr << "Invalid decoder order" << std::endl;
    }
    agreed_message.message = mi.message;
    agreed_messages_.push(agreed_message);
    dual_list_map_.erase(mi.message);
  } else {
    // key does not exist
    dual_list_map_[mi.message] = mi;
  }
}

DLDInfo DualListDecoder::DualListMap::pop_queue() {
  if (queue_size() == 0) {
    std::cerr << "Invalid access to empty queue" << std::endl;
  }
  DLDInfo top = agreed_messages_.top();
  agreed_messages_.pop();
  return top;
}

DLDInfo DualListDecoder::DualListMap::get_top() {
  if (queue_size() == 0) {
    std::cerr << "Invalid access to empty queue" << std::endl;
  }
  DLDInfo top = agreed_messages_.top();
  return top;
}

ListDecoder::messageInformation ListDecoder::traceback_Single(minheap* detourTree, int& numPathsSearched, std::vector<std::vector<int>>& previousPaths, std::vector<std::vector<std::vector<cell>>>& trellisInfo) {
  messageInformation mi;
  bool path_found = false;

  while (!path_found) {
    
		DetourObject detour = detourTree->pop();
		std::vector<int> path(pathLength);

		int newTracebackStage = pathLength - 1;
		double forwardPartialPathMetric = 0;
		int currentState = detour.startingState;

		// if we are taking a detour from a previous path, we skip backwards to the point where we take the
		// detour from the previous path
		if(detour.originalPathIndex != -1){
			forwardPartialPathMetric = detour.forwardPathMetric;
			newTracebackStage = detour.detourStage;

			// while we only need to copy the path from the detour to the end, this simplifies things,
			// and we'll write over the earlier data in any case
			path = previousPaths[detour.originalPathIndex];
			currentState = path[newTracebackStage];

			double suboptimalPathMetric = trellisInfo[detour.startingState][currentState][newTracebackStage].suboptimalPathMetric;

			currentState = trellisInfo[detour.startingState][currentState][newTracebackStage].suboptimalFatherState;
			newTracebackStage--;
			
			double prevPathMetric = trellisInfo[detour.startingState][currentState][newTracebackStage].pathMetric;

			forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
			
		}
		path[newTracebackStage] = currentState;

		// actually tracing back
		for(int stage = newTracebackStage; stage > 0; stage--){
			double suboptimalPathMetric = trellisInfo[detour.startingState][currentState][stage].suboptimalPathMetric;
			double currPathMetric = trellisInfo[detour.startingState][currentState][stage].pathMetric;

			// if there is a detour we add to the detourTree
			if(trellisInfo[detour.startingState][currentState][stage].suboptimalFatherState != -1){
				DetourObject localDetour;
				localDetour.detourStage = stage;
				localDetour.originalPathIndex = numPathsSearched;
				localDetour.pathMetric = suboptimalPathMetric + forwardPartialPathMetric;
				localDetour.forwardPathMetric = forwardPartialPathMetric;
				localDetour.startingState = detour.startingState;
				//queue.push(localDetour);
				detourTree->insert(localDetour);
			}
			currentState = trellisInfo[detour.startingState][currentState][stage].optimalFatherState;
			double prevPathMetric = trellisInfo[detour.startingState][currentState][stage - 1].pathMetric;
			forwardPartialPathMetric += currPathMetric - prevPathMetric;
			path[stage - 1] = currentState;
		}
		previousPaths.push_back(path);


		std::vector<int> message = pathToMessage(path);
    
    if(crc_check(message, crcDegree, crc)){
			mi.message = message;
			mi.path = path;
			mi.listSize = numPathsSearched + 1;
			path_found = true;
		}
    numPathsSearched++;
  }
  return mi;
}

ListDecoder::messageInformation ListDecoder::traceback_deinterleave_Single(minheap* detourTree, int& numPathsSearched, std::vector<std::vector<int>>& previousPaths, std::vector<std::vector<std::vector<cell>>>& trellisInfo, unsigned short int* deinterleaver_ptr) {
  messageInformation mi;
  bool path_found = false;

  while (!path_found) {
    
		DetourObject detour = detourTree->pop();
		std::vector<int> path(pathLength);

		int newTracebackStage = pathLength - 1;
		double forwardPartialPathMetric = 0;
		int currentState = detour.startingState;

		// if we are taking a detour from a previous path, we skip backwards to the point where we take the
		// detour from the previous path
		if(detour.originalPathIndex != -1){
			forwardPartialPathMetric = detour.forwardPathMetric;
			newTracebackStage = detour.detourStage;

			// while we only need to copy the path from the detour to the end, this simplifies things,
			// and we'll write over the earlier data in any case
			path = previousPaths[detour.originalPathIndex];
			currentState = path[newTracebackStage];

			double suboptimalPathMetric = trellisInfo[detour.startingState][currentState][newTracebackStage].suboptimalPathMetric;

			currentState = trellisInfo[detour.startingState][currentState][newTracebackStage].suboptimalFatherState;
			newTracebackStage--;
			
			double prevPathMetric = trellisInfo[detour.startingState][currentState][newTracebackStage].pathMetric;

			forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
			
		}
		path[newTracebackStage] = currentState;

		// actually tracing back
		for(int stage = newTracebackStage; stage > 0; stage--){
			double suboptimalPathMetric = trellisInfo[detour.startingState][currentState][stage].suboptimalPathMetric;
			double currPathMetric = trellisInfo[detour.startingState][currentState][stage].pathMetric;

			// if there is a detour we add to the detourTree
			if(trellisInfo[detour.startingState][currentState][stage].suboptimalFatherState != -1){
				DetourObject localDetour;
				localDetour.detourStage = stage;
				localDetour.originalPathIndex = numPathsSearched;
				localDetour.pathMetric = suboptimalPathMetric + forwardPartialPathMetric;
				localDetour.forwardPathMetric = forwardPartialPathMetric;
				localDetour.startingState = detour.startingState;
				//queue.push(localDetour);
				detourTree->insert(localDetour);
			}
			currentState = trellisInfo[detour.startingState][currentState][stage].optimalFatherState;
			double prevPathMetric = trellisInfo[detour.startingState][currentState][stage - 1].pathMetric;
			forwardPartialPathMetric += currPathMetric - prevPathMetric;
			path[stage - 1] = currentState;
		}
		previousPaths.push_back(path);


		std::vector<int> message = pathToMessage(path);
    std::vector<int> deinterleaved_message;
    for (int i = 0; i < message.size(); i++) {
    	deinterleaved_message.push_back(message[deinterleaver_ptr[i]]);
    }

    if(crc_check(deinterleaved_message, crcDegree, crc)){
			mi.message = deinterleaved_message;
			mi.path = path;
			mi.listSize = numPathsSearched + 1;
			path_found = true;
		}
    numPathsSearched++;
  }
  return mi;
}


DLDInfo DualListDecoder::DualListDecoding_TurboELF(std::vector<double> txSig_0, std::vector<double> txSig_1, unsigned short int* deinterleaver_ptr) {
/*
* This function is the turbo version of the dual list decoder

* Args:
*       - txSig_1: X_r1 + X_sys
*       - txSig_2: X_sys + X_r2
* Output:
*/

  DLDInfo result;
  int pathLength = txSig_1.size() + 1;


  // Construct trellis
  std::vector<std::vector<std::vector<ListDecoder::cell>>> trellisInfo_0;
  std::vector<std::vector<std::vector<ListDecoder::cell>>> trellisInfo_1;
  trellisInfo_0 = list_decoders_[0].constructNTrellis(txSig_0);
  trellisInfo_1 = list_decoders_[1].constructNTrellis(txSig_1);


  int num_total_stages_0 = trellisInfo_0[0][0].size();
  int num_total_stages_1 = trellisInfo_1[0][0].size();
  std::vector<std::vector<int>> prev_paths_list_0;
  std::vector<std::vector<int>> prev_paths_list_1;
  minheap* detourTree_0 = new minheap;
  minheap* detourTree_1 = new minheap;

  // Initialize traceback queue (Detour Tree)
  for(int i = 0; i <  list_decoders_[0].numStates / 2; i++){
		DetourObject detour;
		detour.startingState = i;
		detour.pathMetric = trellisInfo_0[i][i][pathLength - 1].pathMetric;
		detourTree_0->insert(detour);
	}
  for(int i = 0; i < list_decoders_[1].numStates / 2; i++){
		DetourObject detour;
		detour.startingState = i;
		detour.pathMetric = trellisInfo_1[i][i][pathLength - 1].pathMetric;
		detourTree_1->insert(detour);
	}

  // continue traceback until the best combined metric is found
  int num_path_searched_0 = 0;
  bool decoder_0_stop = false;
  bool decoder_0_LSE = false;

  int num_path_searched_1 = 0;
  bool decoder_1_stop = false;
  bool decoder_1_LSE = false;
  bool best_combined_found = false;
  while (!best_combined_found) {

    // check list size exceeded for both decoders 
    if (num_path_searched_0 >= list_decoders_[0].listSize) { decoder_0_LSE = true;}
    if (num_path_searched_1 >= list_decoders_[1].listSize) { decoder_1_LSE = true;}

    // check if both decoders have exceeded list size
    if (decoder_0_LSE && decoder_1_LSE) {break;}

    // decoder 0 traceback
    if (!decoder_0_stop && !decoder_0_LSE) {

      // operate a single traceback
      ListDecoder::messageInformation mi_0 = list_decoders_[0].traceback_Single(
          detourTree_0, num_path_searched_0, prev_paths_list_0, trellisInfo_0);
      mi_0.decoder_index = 0;

      std::cout << "Decoder 0 message: ";
      print_int_vector(mi_0.message);
      std::cout << std::endl;

      if (!mi_0.listSizeExceeded) {
        dual_list_map_.insert(mi_0);
      }
    }

    // decoder 1 traceback
    if (!decoder_1_stop && !decoder_1_LSE) {

      // operate a single traceback
      ListDecoder::messageInformation mi_1 = list_decoders_[1].traceback_deinterleave_Single(
          detourTree_1, num_path_searched_1, prev_paths_list_1, trellisInfo_1, deinterleaver_ptr);
      mi_1.decoder_index = 1;

      std::cout << "Decoder 1 message: ";
      print_int_vector(mi_1.message);
      std::cout << std::endl;

      if (!mi_1.listSizeExceeded) {

        dual_list_map_.insert(mi_1);
      }
    }

    if (dual_list_map_.queue_size() != 0) {
      std::cout << "common queue non zero." << std::endl;
      DLDInfo best_combined = dual_list_map_.pop_queue();      
      return best_combined;
    }
  }

  return result;
}