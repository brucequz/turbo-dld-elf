#ifndef LISTDECODER_H
#define LISTDECODER_H

#include <queue>
#include <map>
#include <vector>
#include <string>

#include "dualtrellis.h"
#include "crc_search_functions.h"
#include "MinHeap.h"
#include "helper_functions.h"
#include "feedbacktrellis.h"

#include <climits>

#ifndef _CODINFO_
#define _CODINFO_
struct codeInformation{
	int k;
	int n;
	int v;
	int crcDeg;
	int crc;
	int numInfoBits;
	std::vector<int> numerators;
	int denominator;
	std::vector<std::vector<int>> hMatrix;
};
#endif

class ListDecoder{
	friend class DualListDecoder;
public:

	struct messageInformation{
		messageInformation(): listSizeExceeded(false) {};
		std::vector<int> message;
		std::vector<int> path;
		int listSize;
		bool listSizeExceeded;
		double metric;
		double full_metric;
		int decoder_index;
	};

	ListDecoder(DualTrellis DT, int listSize, int crcDegree, int crc);

	messageInformation nTrellisDecoding(std::vector<double> receivedMessage);
	std::vector<int> ELFDistanceSpectra_highrate(std::vector<double> receivedMessage);

	messageInformation oneTrellisDecoding(std::vector<double> receivedMessage);
	messageInformation punctured_highrate_turbo_decoder(std::vector<double> receivedMessage, std::vector<int> punc_idx);
	void FindBestCRC_turbo(DualTrellis encTrell, std::vector<double> receivedMessage, codeInformation code);

	messageInformation wavaDecoding(std::vector<double> receivedMessage);
	messageInformation modifiedWAVADecoding(std::vector<double> receivedMessage);

	messageInformation ztListDecoding(std::vector<double> receivedMessage);

	

protected:
	std::vector<std::vector<std::vector<int>>> nextStates;
	int numTrellisSegLength;
	int numStates;
	int numForwardPaths;
	int pathLength;
	int listSize;
	int crcDegree;
	int crc;
	int ztTrailingBits;
	int n;


	struct cell {
		cell(): optimalFatherState(-1), suboptimalFatherState(-1), 
				pathMetric(INT_MAX), suboptimalPathMetric(INT_MAX), init(false) {};
		int optimalFatherState;
		int suboptimalFatherState;
		double pathMetric;
		double suboptimalPathMetric;
		bool init;
	};

	std::vector<int> pathToMessage(std::vector<int>); 
	std::vector<int> ztPathToMessage(std::vector<int> path);

	std::vector<std::vector<std::vector<cell>>> constructNTrellis(std::vector<double> receivedMessage);

	std::vector<std::vector<cell>> constructOneTrellis(std::vector<double> receivedMessage);

	std::vector<std::vector<cell>> constructOneTrellis_SetPunctureZero(std::vector<double> receivedMessage, std::vector<int> punc_idx);

	std::vector<std::vector<cell>> constructWAVATrellis(std::vector<double> receivedMessage, std::vector<std::vector<cell>> oneTrellis);

	std::vector<std::vector<cell>> constructZTTrellis(std::vector<double> receivedMessage);

	messageInformation traceback_Single(minheap* detourTree, int& numPathsSearched, std::vector<std::vector<int>>& previousPaths, std::vector<std::vector<cell>>& trellisInfo);

	messageInformation traceback_deinterleave_Single(minheap* detourTree, int& numPathsSearched, std::vector<std::vector<int>>& previousPaths, std::vector<std::vector<cell>>& trellisInfo, unsigned short int* deinterleaver_ptr);
};

struct DLDInfo {
	DLDInfo(): combined_metric(0.0), discovered_decoder_idx(-1), discovered_partial_metric(-1.0) {
		message = std::vector<int>();
		list_ranks = std::vector<int>();
		list_ranks.push_back(-1);
		list_ranks.push_back(-1);
	};
	int discovered_decoder_idx;
	double discovered_partial_metric;
	double undiscovered_partial_metric;
	double combined_metric;
	std::string return_type;
	std::vector<int> message;
	std::vector<int> list_ranks;
};

class DualListDecoder{
public:
	DualListDecoder(std::vector<codeInformation> code_info, int listSize, std::vector<int> punc_idx);
	~DualListDecoder();

	struct CompareCombinedMetric {
  bool operator()(const DLDInfo& a, const DLDInfo& b) const {
    return a.combined_metric >
           b.combined_metric;  // Lower combined_metric at the top
  }
	};

	// Function to combine two vectors of messageInformation and create a priority
	// queue of DLDInfo sorted in ascending order of combined_metric
	std::priority_queue<DLDInfo, std::vector<DLDInfo>, CompareCombinedMetric>
	combine_maps(const std::vector<ListDecoder::messageInformation>& vec1,
							const std::vector<ListDecoder::messageInformation>& vec2);

  class DualListMap{
  public:
    DualListMap() {};
    ~DualListMap() {};

    void insert(const ListDecoder::messageInformation& mi);
    int queue_size() {return agreed_messages_.size();};
    DLDInfo pop_queue();
    DLDInfo get_top();
		
    
    std::map<std::vector<int>, ListDecoder::messageInformation> dual_list_map_; // dictionary
    std::priority_queue<DLDInfo, std::vector<DLDInfo>, CompareCombinedMetric> agreed_messages_; // priority queue
	};

	DLDInfo DualListDecoding_TurboELF(std::vector<double> txSig_1, std::vector<double> txSig_2, unsigned short int* deinterleaver_ptr);
	DLDInfo DualListDecoding_TurboELF_BAM(std::vector<double> txSig_0, std::vector<double> txSig_1, unsigned short int* interleaver_ptr, unsigned short int* deinterleaver_ptr);

private:
	std::vector<ListDecoder> list_decoders_;
	std::vector<FeedbackTrellis*> trellis_ptrs_; 
	DualListMap dual_list_map_;
	std::vector<int> punc_idx_;
};

#endif