#include "feedbacktrellis.h"
#include <string>
// #include <math.h>
#include "GaloisField.h"
#include "GaloisFieldElement.h"
#include "GaloisFieldPolynomial.h"
#include <algorithm>
#include <cmath>
#include <queue>

static const int V = 6;
// table of primitive polynomials for v=2-15
std::vector<std::string> primitive_polys = {"111",
                                            "1011",
                                            "11001",
                                            "111101",
                                            "1100001",
                                            "10110111",
                                            "110110001",
                                            "1011010111",
                                            "11000001011",
                                            "110111100111",
                                            "1110111110101",
                                            "11111011110011",
                                            "100000101000011",
                                            "1000010000100011"};

FeedbackTrellis::FeedbackTrellis(int k, int n, int v,
                                 std::vector<int> numerators, int denominator) {
  this->k = k;
  this->n = n;
  this->v = v;
  this->denominator = denominator;
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
  // computeTerminations();
}

void FeedbackTrellis::computeNextStates() {
  std::vector<std::vector<int>> gs(k, std::vector<int>(v + 1)); // rename plz
  unsigned int b_arr[V + 1];

  for (int i = 0; i < k; i++) {
    int tempNum = numerators[i]; // octal number in numerator(135)
    int decIn = 0;
    std::string in = std::to_string(tempNum);
    for (int p = (in.length() - 1); p >= 0; p--)
      decIn += (int)(in[p] - '0') * pow(8, (in.length() - p - 1));
    for (int j = v; j >= 0; j--) {
      if (decIn % 2 == 0)
        gs[k - i - 1][j] = 0;
      else
        gs[k - i - 1][j] = 1;
      decIn = decIn / 2;
    }
  }

  // first, construct GF outside the loop with a primitive poly
  unsigned int base_poly[V + 1];
  for (int i = 0; i < V + 1; i++) {
    base_poly[i] = primitive_polys[V - 2][i] - '0';
  }
  galois::GaloisField gf1 = galois::GaloisField(V, base_poly);
  std::cout << "GF construction successful" << std::endl;

  // convert denominator to binary
  std::string b;
  for (int j = 0; j <= v; j++) {
    b += '1';
  }
  int decIn = 0;
  std::string in = std::to_string(denominator);
  for (int p = (in.length() - 1); p >= 0; p--)
    decIn += (int)(in[p] - '0') * pow(8, (in.length() - p - 1));
  for (int j = v; j >= 0; j--) {
    if (decIn % 2 == 0)
      b[j] = '0';
    decIn = decIn / 2;
  }
  for (int i = 0; i < V + 1; i++) {
    b_arr[i] = b[i] - '0';
  }

  // computing nextStates
  for (int currentState = 0; currentState < numStates; currentState++) {
    std::vector<int> cur_sigmas = dec2Bin(currentState, v);
    for (int input = 0; input < numInputSymbols; input++) {
      std::vector<int> us = dec2Bin(input, k);
      unsigned int total_numerator[V + 1];
      for (int i = 0; i < (V + 1); i++) {
        total_numerator[i] = 0;
      }
      for (int ii = 0; ii < us.size(); ii++) {
        std::vector<int> temp(v + 1);
        for (int i = 0; i < (v + 1); i++) {
          temp[i] = gs[ii][i] * us[ii]; // GFCONV
        }
        for (int i = 0; i < (v + 1); i++) {
          total_numerator[i] = temp[i] ^ total_numerator[i];
        }
      }
      for (int i = 1; i < (v + 1); i++) {
        total_numerator[i] = cur_sigmas[i - 1] ^ total_numerator[i];
      }

      // //	DEBUGGING: USED TO DISPLAY NUMERATOR AND DENOMINATOR FOR GFDIV
      // std::cout << "_____________________________________________________" <<
      // std::endl; std::cout << "total numerator: " << std::endl; for(int i =
      // 0; i < V + 1; i++){ 	std::cout << total_numerator[i] << "\t";
      // }
      // std::cout << std::endl;
      // std::cout << "b_arr: " << std::endl;
      // for(int i = 0; i < V + 1; i++){
      // 	std::cout << b_arr[i] << "\t";
      // }
      // std::cout << "decIn: " << decIn << std::endl;

      galois::GaloisFieldElement gfe1[V + 1];
      galois::GaloisFieldElement gfe2[V + 1];

      for (int i = 0; i < (V + 1); i++) {
        gfe1[i] = galois::GaloisFieldElement(&gf1, total_numerator[i]);
      }
      for (int i = 0; i < (V + 1); i++) {
        gfe2[i] = galois::GaloisFieldElement(&gf1, b_arr[i]);
      }
      galois::GaloisFieldPolynomial poly1 =
          galois::GaloisFieldPolynomial(&gf1, V, gfe1);
      galois::GaloisFieldPolynomial poly2 =
          galois::GaloisFieldPolynomial(&gf1, V, gfe2);

      galois::GaloisFieldPolynomial polyremd = poly1 % poly2;
      galois::GaloisFieldPolynomial polyq = poly1 / poly2;

      std::vector<galois::GaloisFieldElement> remdTemp = polyremd.getPoly();
      std::vector<galois::GaloisFieldElement> qTemp = polyq.getPoly();
      std::vector<int> remd;
      std::vector<int> q;
      for (int i = 0; i < remdTemp.size(); i++) {
        remd.push_back(remdTemp[i].poly());
      }
      for (int i = 0; i < qTemp.size(); i++) {
        q.push_back(qTemp[i].poly());
      }
      if (remd.size() == 0) {
        remd.push_back(0);
      }
      if (q.size() == 0) {
        q.push_back(0);
      }

      std::vector<int> testRemd;
      std::vector<int> testQ;

      if (total_numerator[v] == 1) {
        testQ.push_back(1);
        for (int i = 0; i < v + 1; i++) {
          // for(int i = 0; i < v; i++){
          testRemd.push_back(total_numerator[i] ^ b_arr[i]);
        }
      } else {
        testQ.push_back(0);
        for (int i = 0; i < v + 1; i++) {
          testRemd.push_back(total_numerator[i]);
        }
      }
      while (testRemd.size() > 1 && testRemd[testRemd.size() - 1] == 0) {
        testRemd.pop_back();
      }

      // VERBOSE DEBUGGING: USED TO CHECK EQUALITY BETWEEN MY GFDIV AND PROPER
      // GFDIV std::cout << std::endl << "remd: " << std::endl; for(int i = 0; i
      // < remd.size(); i++){ 	std::cout << remd[i] << "\t";
      // }
      // std::cout << std::endl << "testRemd: " << std::endl;
      // for(int i = 0; i < testRemd.size(); i++){
      // 	std::cout << testRemd[i] << "\t";
      // }
      // std::cout << std::endl << "quotient: " << std::endl;
      // for(int i = 0; i < q.size(); i++){
      // 	std::cout << q[i] << std::endl;
      // }
      // std::cout << std::endl << "test quotient: " << std::endl;
      // for(int i = 0; i < testQ.size(); i++){
      // 	std::cout << testQ[i] << std::endl;
      // }

      // CHECKS EQUALITY BETWEEN GALOIS FIELD LIBRARY AND MY IMPLEMENTATION-
      // SHOULD BE INCLUDED UNTIL EQS ARE VALIDATED ACROSS ALL RELEVENT CODES
      if (V != v) {
        std::cout << "constant V does not match- check feedbacktrellis"
                  << std::endl;
      }
      if (testRemd != remd) {
        std::cout << "gfdiv machine broke" << std::endl;
        return;
      }
      if (testQ != q) {
        std::cout << "gfdiv machine broke 2" << std::endl;
        return;
      }

      std::vector<int> next_sigmas(v, 0);
      for (int i = 0; i < v && i < remd.size(); i++) {
        next_sigmas[i] = remd[i];
      }

      int nextState = 0;
      for (int i = (next_sigmas.size() - 1); i >= 0; i--) {
        nextState += next_sigmas[i] *
                     pow(2, (next_sigmas.size() - 1 -
                             i)); // convert next_sigmas back into decimal value
      }

      nextStates[currentState][input] = nextState;
      std::vector<int> nextOutput(n); // k+1
      nextOutput[0] = q[0];           // q has to be size 0 or 1
      for (int i = 1; i < n; i++) {
        nextOutput[i] = us[i - 1];
      }
      int output = 0;
      for (int i = (nextOutput.size() - 1); i >= 0; i--) {
        output += nextOutput[i] * pow(2, nextOutput.size() - 1 - i);
      }
      outputs[currentState][input] = output;
      // std::cout << output << std::endl;
    }
  }
}
std::vector<int> FeedbackTrellis::encoder(std::vector<int> originalMessage) {
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
      std::vector<int> outputBinary = dec2Bin(outputs[State][decimal], n);
      State = nextStates[State][decimal];
      for (int j = 0; j < n; j++) {
        output.push_back(-1 * ((outputBinary[j] * 2) - 1));
      }
    }
    if (m == State) {
      return output;
    }
  }
  return originalMessage;
}

std::vector<int> FeedbackTrellis::dec2Bin(int decimal, int length) {
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

std::vector<std::vector<int>> FeedbackTrellis::getNextStates() {
  return nextStates;
}

std::vector<std::vector<int>> FeedbackTrellis::getOutputs() { return outputs; }

int FeedbackTrellis::getNumInputSymbols() { return numInputSymbols; }

int FeedbackTrellis::getNumOutputSymbols() { return numOutputSymbols; }

int FeedbackTrellis::getNumStates() { return numStates; }

// termination bits for primal trellis
void FeedbackTrellis::computeTerminations() {

  int num_transitions = std::ceil((double)v / (n - 1));
  std::vector<int> tree(numStates);
  std::vector<int> vis(numStates);
  vis[0] = 1;
  std::queue<int> queue;
  queue.push(0);
  while (!queue.empty()) {
    int target_state = queue.front();
    queue.pop();
    for (int i = 0; i < numStates; i++) {
      if (std::count(nextStates[i].begin(), nextStates[i].end(), target_state) >
              0 &&
          vis[i] == 0) {
        vis[i] = 1;
        queue.push(i);
        tree[i] = target_state;
      }
    }
  }

  terminations = std::vector<std::vector<int>>(
      numStates, std::vector<int>(num_transitions));
  for (int i = 0; i < num_transitions; i++) {
    terminations[0][i] = 0;
  }
  int cur = 0;
  for (int i = 1; i < numStates; i++) {
    cur = i;
    for (int j = 0; j < num_transitions; j++) {
      int fa_state = tree[cur];
      int index = std::distance(
          nextStates[cur].begin(),
          std::find(nextStates[cur].begin(), nextStates[cur].end(), fa_state));
      terminations[i][j] = index;
      cur = fa_state;
    }
  }
  std::cout << "made it to the end of the termination computations"
            << std::endl;
}

std::vector<int>
FeedbackTrellis::terminateMsg(std::vector<int> originalMessage) {
  // first encode to find the final state
  std::vector<int> output;
  int finalState = 0;
  for (int i = 0; i < originalMessage.size(); i += k) {
    int decimal = 0;
    for (int j = 0; j < k; j++) {
      decimal += (originalMessage[i + j] * pow(2, k - j - 1));
    }
    std::vector<int> outputBinary = dec2Bin(outputs[finalState][decimal], n);
    finalState = nextStates[finalState][decimal];
  }

  for (int i = 0; i < originalMessage.size(); i++) {
    output.push_back(originalMessage[i]);
  }

  // appending the termination bits that force the zero state
  std::vector<int> localTerminations = terminations[finalState];
  for (int i = 0; i < localTerminations.size(); i++) {
    std::vector<int> outputBinary = dec2Bin(localTerminations[i], k);
    for (int j = 0; j < k; j++) {
      output.push_back(outputBinary[j]);
    }
  }
  return output;
}

std::vector<int>
FeedbackTrellis::ztencoder(std::vector<int> terminatedMessage) {
  // unlike for tbcc, we always start at the zero state, add the crc, then add
  // bits to the message to guarentee we return to the zero state when we encode

  std::vector<int> output;
  int State = 0;
  for (int i = 0; i < terminatedMessage.size(); i += k) {
    int decimal = 0;
    for (int j = 0; j < k; j++) {
      decimal += (terminatedMessage[i + j] * pow(2, k - j - 1));
    }
    std::vector<int> outputBinary = dec2Bin(outputs[State][decimal], n);
    State = nextStates[State][decimal];
    for (int j = 0; j < n; j++) {
      output.push_back(-1 * ((outputBinary[j] * 2) - 1));
    }
  }
  if (State != 0) {
    std::cout << "error: not zero terminated" << std::endl;
  }
  return output;
}

std::vector<std::vector<int>> FeedbackTrellis::getTerminations() {
  return terminations;
}