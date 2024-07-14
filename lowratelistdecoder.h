#ifndef LOWRATELISTDECODER_H
#define LOWRATELISTDECODER_H

#include "dualtrellis.h"
#include "crc_search_functions.h"
#include "MinHeap.h"
#include "helper_functions.h"
#include "feedforwardtrellis.h"

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

class LowRateListDecoder{
public:
	LowRateListDecoder(FeedforwardTrellis FT, int listSize, int crcDegree, int crc);
	LowRateListDecoder(FeedforwardTrellis FT, int listSize, int crcDegree, int crc, std::vector<std::vector<int>> neighboring_list, std::vector<std::vector<int>> neighboring_msg, std::vector<std::vector<int>> path_ie_state);

	struct messageInformation{
		messageInformation(LowRateListDecoder * decoder){
			neighbors = std::vector<std::vector<int>>(decoder->listSize, std::vector<int>(516,0));
			neighbor_msgs = std::vector<std::vector<int>>(decoder->listSize, std::vector<int>(43,0));
			path_ie = std::vector<std::vector<int>>(decoder->listSize, std::vector<int>());
		};
		messageInformation(){};
		std::vector<std::vector<int>> neighbors; // ${listSize} x 516 matrix
		std::vector<std::vector<int>> neighbor_msgs; // ${listSize} x 43 matrix
		std::vector<std::vector<int>> path_ie; // initial and ending states of path
		std::vector<int> message;
		std::vector<int> path;
		int listSize;
		bool listSizeExceeded = false;
		double metric = -1;
		//can potentially add more information as needed for debugging
	};

	messageInformation lowRateDecoding(std::vector<double> receivedMessage);
	messageInformation multiTrellisLowRateDecoding(std::vector<double> receivedMessage);
	messageInformation offsetLowRateDecoding(std::vector<double> receivedMessage);
	messageInformation linearityLowRateListSearching(std::vector<double> receivedMessage);
	messageInformation offsetMultiTrellisLowRateDecoding(std::vector<double> receivedMessage);
	messageInformation linearityLowRateDecoding(std::vector<double> receivedMessage);
	messageInformation TBGuaranteedListSearching(std::vector<double> receivedMessage);
	messageInformation TBGuaranteedLinearityDecoding(std::vector<double> receivedMessage);
	messageInformation lowRateDecoding_listSize(std::vector<double> receivedMessage);
	messageInformation lowRateDecoding_recip(std::vector<double> receivedMessage);
	messageInformation lowRateReciprocityCheck(std::vector<double> originalPoint, std::vector<double> distantPoint = {});
	std::vector<double> lowRateDistanceDistribution(std::vector<double> receivedMessage);
	void printRangeOfCodewords(std::vector<double> receivedMessage, int start, int end);

	void FindBestCRC_TBCC(FeedforwardTrellis, std::vector<double> receivedMessage, codeInformation code);
	void Find_ELFs_TB_Improved(FeedforwardTrellis, std::vector<double> receivedMessage, codeInformation code);
	std::vector<int> ELFDistanceSpectra(std::vector<double> receivedMessage);

private:
	int numForwardPaths;
	int listSize;
	int crcDegree;
	int crc;
	int n;

	std::vector<std::vector<int>> lowrate_nextStates;
	std::vector<std::vector<int>> lowrate_outputs;
	std::vector<std::vector<int>> neighboring_list; // ${listSize} x 516 matrix
	std::vector<std::vector<int>> neighboring_msg;  // ${listSize} x 43 matrix
	std::vector<std::vector<int>> path_ie_state;
	int lowrate_numStates;
	int lowrate_symbolLength;
	int lowrate_pathLength;

	struct cell {
		int optimalFatherState = -1;
		int suboptimalFatherState = -1;
		double pathMetric = INT_MAX;
		double suboptimalPathMetric = INT_MAX;
		bool init = false;
	};
    std::vector<int> pathToMessage(std::vector<int>); 
    std::vector<int> pathToCodeword(std::vector<int>); 
	std::vector<std::vector<cell>> constructLowRateTrellis(std::vector<double> receivedMessage);
	std::vector<std::vector<std::vector<cell>>> constructLowRateMultiTrellis(std::vector<double> receivedMessage);
	
};


#endif