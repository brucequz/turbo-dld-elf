#ifndef FEEDBACKTRELLIS_H
#define FEEDBACKTRELLIS_H
#include <vector>
#include <string>

class FeedbackTrellis {
public:
	FeedbackTrellis(int k, int n, int v, std::vector<int> numerators, int denominator);
	std::vector<int> ztencoder(std::vector<int> originalMessage);
	std::vector<int> terminateMsg(std::vector<int> originalMessage);
	std::vector<int> encoder(std::vector<int> originalMessage);
	std::vector<std::vector<int> > getNextStates();
	std::vector<std::vector<int> > getOutputs();
	std::vector<std::vector<int> > getTerminations();
	int getNumInputSymbols();
	int getNumOutputSymbols();
	int getNumStates();
private:
	int k;
	int n;
	int v;
	int numInputSymbols;
	int numOutputSymbols;
	int numStates;
	std::vector<int> numerators;
	int denominator;
	std::vector<std::vector<int> > nextStates;
	std::vector<std::vector<int> > outputs;
	std::vector<std::vector<int> > terminations;

	void computeTerminations();
	void computeNextStates();

	std::vector<int> dec2Bin(int decimal, int length);
};
#endif
