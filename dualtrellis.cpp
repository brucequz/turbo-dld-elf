#include "dualtrellis.h"
#include <cmath>

DualTrellis::DualTrellis(std::vector<std::vector<int>> H)
    : numInputSymbols(2), numOutputSymbols(1) {
  v = H[0].size(); // # of columns
  c = find_c(H);
  n = H.size(); // # of rows
  Nstates = std::pow(2, v);
  half_Nstates = Nstates / 2;
  for (int i = 0; i < n; i++) {
    std::vector<int> temp;
    for (int j = 0; j < v; j++) {
      temp.push_back(H[i][v - j - 1]);
    }
    H[i] = temp;
  }
  std::vector<int> cur_valid_states(Nstates, 0);
  std::vector<int> next_valid_states(Nstates, 0);
  this->subsequentStates = std::vector<std::vector<std::vector<int>>>(
      Nstates, std::vector<std::vector<int>>(n, std::vector<int>(2)));
  this->prevStates = std::vector<std::vector<std::vector<int>>>(
      Nstates, std::vector<std::vector<int>>(n, std::vector<int>(2, -1)));

  for (int i = 0; i < (half_Nstates); i++) {
    cur_valid_states[i] = 1;
  }
  // find nextStates
  for (int i = 0; i < n; i++) {
    if ((i != c) && (i != (n - 1))) {
      for (int cur_state = 0; cur_state < Nstates; cur_state++) {
        if (cur_valid_states[cur_state] == 1) {
          std::vector<int> cur_state_bin = dec2Bin(cur_state, v);
          for (int input = 0; input < numInputSymbols; input++) {
            std::vector<int> next_state_bin(v);
            int next_state;
            for (int p = 0; p < v; p++) {
              next_state_bin[p] = (cur_state_bin[p] + (input * H[i][p])) % 2;
            }
            next_state = bin2Dec(next_state_bin);
            next_valid_states[next_state] = 1;
            subsequentStates[cur_state][i][input] = next_state;
          }
        } else {
          for (int b = 0; b < numInputSymbols; b++) {
            subsequentStates[cur_state][i][b] = -2;
          }
        }
      }
    } else if (i == c && i != (n - 1)) {
      for (int cur_state = 0; cur_state < (Nstates); cur_state++) {
        if (cur_valid_states[cur_state] == 1) {
          std::vector<int> cur_state_bin = dec2Bin(cur_state, v);
          for (int input = 0; input < (numInputSymbols - 1); input++) {
            std::vector<int> next_state_bin(v);
            for (int p = 0; p < v; p++) {
              next_state_bin[p] =
                  (cur_state_bin[p] + (cur_state_bin[v - 1] * H[i][p])) % 2;
            }
            int next_state = bin2Dec(next_state_bin);
            next_valid_states[next_state] = 1;
            subsequentStates[cur_state][i][cur_state_bin[v - 1]] = next_state;
            subsequentStates[cur_state][i][1 - cur_state_bin[v - 1]] = -1;
          }
        } else {
          for (int b = 0; b < numInputSymbols; b++) {
            subsequentStates[cur_state][i][b] = -2;
          }
        }
      }
    } else if (i != c && i == (n - 1)) {
      for (int cur_state = 0; cur_state < Nstates; cur_state++) {
        if (cur_valid_states[cur_state] == 1) {
          std::vector<int> cur_state_bin = dec2Bin(cur_state, v);
          for (int input = 0; input < (numInputSymbols - 1); input++) {
            std::vector<int> next_state_bin(v);
            for (int p = 0; p < v; p++) {
              next_state_bin[p] = (cur_state_bin[p] + (input * H[i][p])) % 2;
            }
            for (int b = (v - 1); b > 0; b--) {
              next_state_bin[b] = next_state_bin[b - 1];
            }
            next_state_bin[0] = 0;
            int next_state = bin2Dec(next_state_bin);
            next_valid_states[next_state] = 1;
            subsequentStates[cur_state][i][input] = next_state;
          }
        } else {
          for (int b = 0; b < numInputSymbols; b++) {
            subsequentStates[cur_state][i][b] = -2;
          }
        }
      }
    } else {
      for (int cur_state = 0; cur_state < Nstates; cur_state++) {
        if (cur_valid_states[cur_state] == 1) {
          std::vector<int> cur_state_bin = dec2Bin(cur_state, v);
          std::vector<int> next_state_bin(v);
          for (int p = 0; p < v; p++) {
            next_state_bin[p] =
                (cur_state_bin[p] + (cur_state_bin[v - 1] * H[i][p])) % 2;
          }
          for (int b = (v - 1); b > 0; b--) {
            next_state_bin[b] = next_state_bin[b - 1];
          }
          next_state_bin[0] = 0;
          int next_state = bin2Dec(next_state_bin);
          next_valid_states[next_state] = 1;
          subsequentStates[cur_state][i][cur_state_bin[v - 1]] = next_state;
          subsequentStates[cur_state][i][1 - cur_state_bin[v - 1]] = -1;
        } else {
          for (int b = 0; b < numInputSymbols; b++) {
            subsequentStates[cur_state][i][b] = -2;
          }
        }
      }
    }
    cur_valid_states = next_valid_states;
    next_valid_states = std::vector<int>(Nstates, 0);
  }

  for (int i = 0; i < subsequentStates.size(); i++) {
    for (int j = 0; j < subsequentStates[0].size(); j++) {
      for (int k = 0; k < subsequentStates[0][0].size(); k++) {
        if (subsequentStates[i][j][k] >= 0)
          this->prevStates[subsequentStates[i][j][k]][j][k] = i;
      }
    }
  }
}

int DualTrellis::find_c(std::vector<std::vector<int>> H) {
  for (int i = (H.size() - 1); i >= 0; i--) {
    if (H[i][0] == 1) {
      return i;
    }
  }
  return -1;
}

std::vector<int> DualTrellis::dec2Bin(int decimal, int length) {
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
int DualTrellis::bin2Dec(std::vector<int> binary) {
  int decimal = 0;
  for (int i = (binary.size() - 1); i >= 0; i--) {
    decimal += (binary[i] * pow(2, (binary.size() - i - 1)));
  }
  return decimal;
}
std::vector<std::vector<std::vector<int>>> DualTrellis::getsubStates() {
  return subsequentStates;
}

std::vector<std::vector<std::vector<int>>> DualTrellis::getPrevStates() {
  return prevStates;
}

int DualTrellis::getV() { return v; }

int DualTrellis::getN() { return n; }