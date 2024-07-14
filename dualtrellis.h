#ifndef DUALTRELLIS_H
#define DUALTRELLIS_H
#include <vector>
#include <iostream>

class DualTrellis {
public:
	DualTrellis(std::vector<std::vector<int>> H);
	std::vector<std::vector<std::vector<int>>> getsubStates();
	std::vector<std::vector<std::vector<int>>> getPrevStates();
	int getV();
	int getN();
private:
	int v;
	int c;
	int n;
	int numInputSymbols;
	int numOutputSymbols;
	int numStates;
	int nextStates;
	int Nstates;
	int half_Nstates;
	std::vector<int> cur_valid_states;
	std::vector<int> next_valid_states;
	std::vector<std::vector<std::vector<int>>> subsequentStates;
	std::vector<std::vector<std::vector<int>>> prevStates;
	std::vector<std::vector<int>> H;//parity check matrix
	int find_c(std::vector<std::vector<int>> H);
	std::vector<int> dec2Bin(int decimal, int length);
	int bin2Dec(std::vector<int> binary);
};

#endif