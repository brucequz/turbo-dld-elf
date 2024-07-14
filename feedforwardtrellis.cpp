#include "feedforwardtrellis.h"
#include <string>
// #include <math.h>
#include "GaloisField.h"
#include "GaloisFieldElement.h"
#include "GaloisFieldPolynomial.h"
#include "helper_functions.h"
#include <algorithm>
#include <cmath>
#include <queue>

static const int V = 6;

FeedforwardTrellis::FeedforwardTrellis(int k, int n, int v,
                                       std::vector<int> numerators) {
  this->k = k;
  this->n = n;
  this->v = v;
  for (int i = 0; i < numerators.size(); i++) {
    this->numerators.push_back(numerators[i]);
  }
  this->numInputSymbols = pow(2.0, k);
  this->numOutputSymbols = pow(2.0, n);
  this->numStates = pow(2.0, v);
  this->nextStates = std::vector<std::vector<int>>(
      numStates, std::vector<int>(numInputSymbols));
  this->outputs = std::vector<std::vector<int>>(
      numStates, std::vector<int>(numInputSymbols));

  if (v != V) {
    std::cout
        << "MAJOR ISSUE: CONST V DOES NOT MATCH v. EDIT IN FEEDBACKTRELLIS"
        << std::endl;
    exit(1);
  }
  computeNextStates();

  // to check nextState values
  // std::cout << "nextstates:" << std::endl;
  // for (int i=0; i<256; i++){
  // 	for (int j=0; j<2; j++){
  // 		std::cout << nextStates[i][j];
  // 	}
  // 	std::cout << std::endl;
  // }

  // to check outputs values
  // std::cout << "outputs:" << std::endl;
  // for (int i=0; i<256; i++){
  // 	for (int j=0; j<2; j++){
  // 		std::cout << outputs[i][j];
  // 	}
  // 	std::cout << std::endl;
  // }
}

void FeedforwardTrellis::computeNextStates() {
  // convert to binary numerators
  std::vector<std::vector<int>> bin_numerators(n, std::vector<int>(v + 1));
  for (int i = 0; i < n; i++) {
    int tempNum = numerators[i]; // octal number in numerator(135)
    int decIn = 0;
    std::string in = std::to_string(tempNum);
    for (int p = (in.length() - 1); p >= 0; p--)
      decIn += (int)(in[p] - '0') * pow(8, (in.length() - p - 1));
    for (int j = v; j >= 0; j--) {
      if (decIn % 2 == 0)
        bin_numerators[i][j] = 0;
      else
        bin_numerators[i][j] = 1;
      decIn = decIn / 2;
    }
  }
  // calculate next states and outputs
  for (int currentState = 0; currentState < numStates; currentState++) {
    std::vector<int> mem_elements = dec2Bin(currentState, v + 1);
    for (int input = 0; input < numInputSymbols; input++) {
      mem_elements[0] = input;
      std::vector<int> output(n);
      for (int i = 0; i < n; i++) {
        output[i] = 0;
      }
      for (int x_bit = 0; x_bit < n; x_bit++) {
        for (int m_bit = 0; m_bit < V + 1; m_bit++) {
          if (bin_numerators[x_bit][m_bit] == 1) {
            output[x_bit] ^= mem_elements[m_bit];
          }
        }
      }
      outputs[currentState][input] = bin2Dec(output);
      std::vector<int> temp(V);
      for (int i = 0; i < V; i++) {
        temp[i] = mem_elements[i];
        // std::cout << temp[i] << std::endl;
      }
      nextStates[currentState][input] = bin2Dec(temp);
    }
  }
}

std::vector<int> FeedforwardTrellis::encoder(std::vector<int> originalMessage) {
  // brute force approach, there is a better way to do this assuming
  // invertibility that allows us to precompute starting / ending states,
  // reducing complexity in each encoding from O(numStates) to O(2). revisit
  // when available

  for (int m = 0; m < numStates; m++) {
    std::vector<int> output;
    int State = m;
    for (int i = 0; i < originalMessage.size(); i += k) {
      int decimal = 0;
      for (int j = 0; j < k; j++) {
        decimal += (originalMessage[i + j] * pow(2, k - j - 1));
      }
      std::vector<int> outputBinary = get_point(outputs[State][decimal], n);
      State = nextStates[State][decimal];
      for (int j = 0; j < n; j++) {
        output.push_back(outputBinary[j]);
      }
    }
    if (m == State) {
      return output;
    }
  }
  return originalMessage;
}

std::vector<int> FeedforwardTrellis::dec2Bin(int decimal, int length) {
  std::vector<int> binary(length);
  for (int j = (length - 1); j >= 0; j--) {
    if (decimal % 2 == 0)
      binary[j] = 0;
    else
      binary[j] = 1;
    decimal = decimal / 2;
  }
  return binary;
}

int FeedforwardTrellis::bin2Dec(std::vector<int> binary) {
  int decimal = 0;
  for (int i = (binary.size() - 1); i >= 0; i--) {
    decimal += (binary[i] * pow(2, (binary.size() - i - 1)));
  }
  return decimal;
}

std::vector<std::vector<int>> FeedforwardTrellis::getNextStates() {
  return nextStates;
}

std::vector<std::vector<int>> FeedforwardTrellis::getOutputs() {
  return outputs;
}

int FeedforwardTrellis::getNumInputSymbols() { return numInputSymbols; }

int FeedforwardTrellis::getNumOutputSymbols() { return numOutputSymbols; }

int FeedforwardTrellis::getNumStates() { return numStates; }

int FeedforwardTrellis::getV() { return v; }

int FeedforwardTrellis::getN() { return n; }