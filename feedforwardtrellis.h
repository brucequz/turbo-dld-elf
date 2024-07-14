#ifndef FEEDFORWARDTRELLIS_H
#define FEEDFORWARDTRELLIS_H
#include <vector>
#include <string>

class FeedforwardTrellis {
public:
	FeedforwardTrellis(int k, int n, int v, std::vector<int> numerators);
	std::vector<int> encoder(std::vector<int> originalMessage);
	std::vector<std::vector<int>> getNextStates();
	std::vector<std::vector<int>> getOutputs();
	int getNumInputSymbols();
	int getNumOutputSymbols();
	int getNumStates();
	int getV();
	int getN();
private:
	int k;
	int n;
	int v;
	int numInputSymbols;
	int numOutputSymbols;
	int numStates;
	std::vector<int> numerators;
	std::vector<std::vector<int>> nextStates;
	std::vector<std::vector<int>> outputs;

	void computeNextStates();

	std::vector<int> dec2Bin(int decimal, int length);
    int bin2Dec(std::vector<int> binary);
};
#endif
