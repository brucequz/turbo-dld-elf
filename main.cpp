#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <sstream>
#include <time.h>
#include <vector>
using namespace std;
#include "GaloisField.h"
#include "GaloisFieldElement.h"
#include "GaloisFieldPolynomial.h"
#include "dualtrellis.h"
#include "feedbacktrellis.h"
#include "feedforwardtrellis.h"
#include "helper_functions.h"
#include "listdecoder.h"
#include "lowratelistdecoder.h"

#include "MinHeap.h"

int LISTSIZE = 1e3;

static default_random_engine generator;

void turbo_ELF(codeInformation code);
void find_optimal_crc_turbo(codeInformation code);

std::vector<double> addNoise(std::vector<int> encodedMsg, double SNR);
void crc_calculation(std::vector<int> &input_data, int crc_bits_num,
                     int crc_dec);
int make_file_interleaver(char interleaver_file[],
                          unsigned short int interleaver[], int n);

int comboTest(codeInformation code);
void IRWEFtest(codeInformation code);
void IRWEFtest2(codeInformation code);
void find_dmin(codeInformation code);
void IRWEFdistribution(codeInformation code);
void set_combo_lists(unsigned long long combo_list[][20], int K);
void setWeightIJ(const unsigned long long combos[][20], int informaton[], int h,
                 unsigned long long int item, int Km);

void elf_turbo_simulation(codeInformation code);

/*
int main(){
        codeInformation code;
        code.k = 2;						// numerator of
the rate code.n = 3;						// denominator
of the rate code.v = 8;						// number of
memory elements code.crcDeg = 13;      // m+1 code.crc = 0xFFF; code.numInfoBits
= 128;

        // optimal code numerators and denominator are known, and are given in
octal code.numerators = {555, 631}; code.denominator = 477;
        // hMatrix is {{denom},{numerators[last]},...,{numerators[0]}},
converted to binary then flipped lr code.hMatrix = { {1,1,1,1,1,1,0,0,1},
{1,0,0,1,1,0,0,1,1},{1,0,1,1,0,1,1,0,1}}; //477

        // find_optimal_crc_turbo(code);
        turbo_ELF(code);

        return 0;
}
*/

int main() {
  codeInformation code;
  code.k = 1;      // numerator of the rate
  code.n = 2;      // denominator of the rate
  code.v = 6;      // number of memory elements
  code.crcDeg = 8; // m+1, degree of CRC, # bits of CRC polynomial
  code.crc = 215;
  code.numInfoBits = 64;

  // optimal code numerators and denominator are known, and are given in octal
  code.numerators = {133};
  code.denominator = 171;
  code.hMatrix = {{1, 0, 0, 1, 1, 1, 1}, {1, 1, 0, 1, 1, 0, 1}};

  // comboTest(code);
  // IRWEFtest(code);
  // IRWEFdistribution(code);
  // IRWEFtest2(code);
  // find_dmin(code);
  elf_turbo_simulation(code);

  return 0;
}

void elf_turbo_simulation(codeInformation code) {
  int Km = code.numInfoBits + code.crcDeg - 1; // 71

  // puncture 43 bits to achieve overall rate-1/2 when m=7
  std::vector<int> punc_idx = {0,  1,  2,  3,  4,  6,  9,  13, 14, 15, 18,
                               19, 21, 22, 23, 25, 26, 27, 28, 29, 30, 31,
                               32, 33, 34, 35, 36, 37, 38, 39, 40, 42, 45,
                               46, 52, 53, 54, 55, 59, 60, 61, 62, 63};

  FeedbackTrellis encodingTrellis(code.k, code.n, code.v, code.numerators,
                                  code.denominator);

  // interleaver
  char interleaver_file[1024]; // make sure the size is large enough
  // only K+m because half of the bits are systematic
  snprintf(interleaver_file, sizeof(interleaver_file),
           "interleaver_%d_S5_T1.dat", Km);
  unsigned short int interleaver[Km]; // Adjust size as needed
  make_file_interleaver(interleaver_file, interleaver, Km);
  unsigned short int deinterleaver[Km]; // Adjust size as needed
  for (int i = 0; i < Km; i++) {
    deinterleaver[interleaver[i]] = i;
  }

  // std::vector<int> deinterleaved_data_vector; // deinterleave example
  // for (int i = 0; i < Km*2; i++) {
  // 	deinterleaved_data_vector.push_back(interleaved_data_vector[deinterleaver[i]]);
  // }

  vector<double> SNR = {20};
  // outer loop: SNR
  for (int s = 0; s < SNR.size(); s++) {
    double cur_SNR = SNR[s];

    int targetedErrors = 100;
    int numerror = 0;
    int numtrial = 0;
    // inner loop: MC trials
    while (numerror < targetedErrors || numtrial < 1e4) {

      std::cout << "trial number: " << numtrial << std::endl;

      std::vector<int> original_message;
      for (int i = 0; i < code.numInfoBits; i++) {
        original_message.push_back(rand() % 2);
        // original_message.push_back(0);
      }

      std::cout << "original_message message: ";
      print_int_vector(original_message);
      std::cout << std::endl;

      crc_calculation(original_message, code.crcDeg, code.crc);
      std::vector<int> encodedMessage =
          encodingTrellis.encoder(original_message);
      // get parity bits
      std::vector<int> X_R1;
      std::vector<int> X_sys;
      for (int j = 0; j < encodedMessage.size(); j += 2) {
        X_sys.push_back(encodedMessage[j + 1]);
        X_R1.push_back(encodedMessage[j]);
      }

      std::cout << "X_sys message: ";
      print_int_vector(X_sys);
      std::cout << std::endl;

      std::cout << "X_R1 message: ";
      print_int_vector(X_R1);
      std::cout << std::endl;

      // interleaved original message
      std::vector<int> pi_original_message;
      for (int ii = 0; ii < Km; ii++) {
        pi_original_message.push_back(original_message[interleaver[ii]]);
      }

      std::cout << "pi_original_message: ";
      print_int_vector(pi_original_message);
      std::cout << std::endl;

      // encode interleaved message
      std::vector<int> encodedMessage_inter =
          encodingTrellis.encoder(pi_original_message);
      // get parity bits
      std::vector<int> X_R2;
      for (int j = 0; j < encodedMessage_inter.size(); j += 2) {
        X_R2.push_back(encodedMessage_inter[j]);
      }

      std::cout << "X_R2 message: ";
      print_int_vector(X_R2);
      std::cout << std::endl;

      // add noise
      std::vector<double> Y_R1 = addNoise(X_R1, cur_SNR);
      std::vector<double> Y_R2 = addNoise(X_R2, cur_SNR);
      std::vector<double> Y_sys = addNoise(X_sys, cur_SNR);
      // puncture both parity sequences
      for (int p = 0; p < Y_R1.size(); p++) {
        for (int q = 0; q < punc_idx.size(); q++) {
          if (p == punc_idx[q]) {
            Y_R1[p] = 0;
            Y_R2[p] = 0;
          }
        }
      }
      // interleave Y_sys for bottom list decoder
      std::vector<double> pi_Y_sys;
      for (int ii = 0; ii < Km; ii++) {
        pi_Y_sys.push_back(Y_sys[interleaver[ii]]);
      }

      // Y_R1.insert(Y_R1.begin(), Y_sys.begin(), Y_sys.end());
      // Y_R2.insert(Y_R2.begin(), pi_Y_sys.begin(), pi_Y_sys.end());

      std::vector<double> DLD_R1;
      std::vector<double> DLD_R2;

      for (int i = 0; i < Y_R1.size(); i++) {
        DLD_R1.push_back(Y_R1[i]);
        DLD_R1.push_back(Y_sys[i]);
      }
      for (int i = 0; i < Y_R2.size(); i++) {
        DLD_R2.push_back(Y_R2[i]);
        DLD_R2.push_back(pi_Y_sys[i]);
      }

      std::vector<codeInformation> codeList = {code, code};
      DualListDecoder DLD(codeList, LISTSIZE);
      DLDInfo result =
          DLD.DualListDecoding_TurboELF(DLD_R1, DLD_R2, deinterleaver);

      numtrial++;
    }
  }
}

void IRWEFdistribution(codeInformation code) {
  int Km = code.numInfoBits + code.crcDeg - 1;
  //    // puncture 44 bits to achieve overall rate-1/2 for m=8
  //    std::vector<int> punc_idx = {0,  1,  2,  3,  4,  6,  9, 13, 14, 15, 18,
  //    19, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37,
  //    38, 39, 40, 42, 45, 46, 52, 53, 54, 55, 59, 60, 61, 62, 63};

  // puncture 43 bits to achieve overall rate-1/2 when m=7
  std::vector<int> punc_idx = {0,  1,  2,  3,  4,  6,  9,  13, 14, 15, 18,
                               19, 21, 22, 23, 25, 26, 27, 28, 29, 30, 31,
                               32, 33, 34, 35, 36, 37, 38, 39, 40, 42, 45,
                               46, 52, 53, 54, 55, 59, 60, 61, 62, 63};

  int INTERLEAVER_LENGTH = Km;
  int info[INTERLEAVER_LENGTH];
  unsigned long long combos[Km + 1][20];
  set_combo_lists(combos, Km);

  FeedbackTrellis encodingTrellis(code.k, code.n, code.v, code.numerators,
                                  code.denominator);

  // interleaver
  char interleaver_file[1024]; // make sure the size is large enough
  // only K+m because half of the bits are systematic
  snprintf(interleaver_file, sizeof(interleaver_file),
           "interleaver_%d_S5_T1.dat", Km);
  unsigned short int interleaver[Km]; // Adjust size as needed
  make_file_interleaver(interleaver_file, interleaver, Km);
  unsigned short int deinterleaver[Km]; // Adjust size as needed
  for (int i = 0; i < Km; i++) {
    deinterleaver[interleaver[i]] = i;
  }

  // std::vector<int> deinterleaved_data_vector; // deinterleave example
  // for (int i = 0; i < Km*2; i++) {
  // 	deinterleaved_data_vector.push_back(interleaved_data_vector[deinterleaver[i]]);
  // }

  // store CRC dist spectra in a .txt file
  std::ofstream distSpectraFile;
  string filename = "turbo_full_dist_spectra_v6m7K64_ELF147_ELF183_ELF215.txt";
  distSpectraFile.open(filename);
  // // store overall dist spectra in a .txt file
  // std::ofstream allSpectraFile;
  // string filename2 = "turbo_all_spectra_v6m8K64.txt";
  // allSpectraFile.open(filename2);

  std::vector<int> candid_CRC = {147, 183, 215};
  // outer loop: CRC
  for (int c = 0; c < candid_CRC.size(); c++) {
    int crc_poly = candid_CRC[c];
    distSpectraFile << "Current CRC: " << crc_poly << std::endl;
    // inner loop: go through all message sequences, record distance spectra for
    // each candid CRC
    std::vector<int> dist_spectra(31, 0); // initialize to all 0s. starting with
                                          // weight=0, goes up to weight 30
    for (int h = 1; h <= 5; h++) {
      for (int i = 0; i < combos[Km][h]; i++) {
        setWeightIJ(combos, info, h, i, Km);
        vector<int> info_v(info, *(&info + 1)); // int array to vector
        // print_int_vector(info_v);
        if (crc_check(info_v, code.crcDeg, crc_poly)) {
          // encode message
          std::vector<int> encodedMessage = encodingTrellis.encoder(info_v);
          // when encoding fails (ie. no matching starting/ending state), the
          // message is returned
          if (encodedMessage.size() == info_v.size()) {
            continue;
            std::cout << "failed encoding" << std::endl;
          }
          // get rid of systematic bits
          std::vector<int> parityBits;
          for (int j = 0; j < encodedMessage.size(); j += 2) {
            parityBits.push_back(encodedMessage[j]);
          }
          // interleave message bits
          std::vector<int> interleave_info;
          for (int ii = 0; ii < Km; ii++) {
            interleave_info.push_back(info_v[interleaver[ii]]);
          }
          // encode interleaved message
          std::vector<int> encodedMessage_inter =
              encodingTrellis.encoder(interleave_info);
          // get rid of systematic bits
          std::vector<int> parityBits_inter;
          for (int j = 0; j < encodedMessage_inter.size(); j += 2) {
            parityBits_inter.push_back(encodedMessage_inter[j]);
          }
          // puncture both parity seauences & record weight
          for (int p = 0; p < parityBits.size(); p++) {
            for (int q = 0; q < punc_idx.size(); q++) {
              if (p == punc_idx[q]) {
                parityBits[p] = 0;
                parityBits_inter[p] = 0;
              }
            }
          }
          int parityWeight =
              std::count(parityBits.begin(), parityBits.end(), -1);
          int interleaveWeight =
              std::count(parityBits_inter.begin(), parityBits_inter.end(), -1);
          // calculate total cwd weight
          int totalWeight = h + parityWeight + interleaveWeight;
          // store in dist_spectra (up to weight 30)
          if (totalWeight <= 30)
            dist_spectra[totalWeight]++;
        }
      }
    }
    for (int d = 0; d < dist_spectra.size(); d++) {
      distSpectraFile << "weight-" << d << ": " << dist_spectra[d] << std::endl;
    }
  }
  distSpectraFile.close();
}

void find_dmin(codeInformation code) {
  // run a straightforward sieve and find the CRC that maximizes the dmin of the
  // (100,64 code). We can do this sieve search for each of the two component
  // codes. We want to keep track of the minimum distance of all the CRCs not
  // just the best one. codeword corresponds to the all-zeros message
  std::vector<double> rxSig;
  // this code is (64,100)
  for (int i = 0; i < 144; i++) {
    rxSig.push_back(1);
  }

  DualTrellis trell(code.hMatrix);
  int L = 2;
  ListDecoder decoder = ListDecoder(trell, L, code.crcDeg, code.crc);
  ListDecoder::messageInformation decoded = decoder.oneTrellisDecoding(rxSig);
}

void IRWEFtest2(codeInformation code) {
  int Km = code.numInfoBits + code.crcDeg - 1;
  //    // puncture 44 bits to achieve overall rate-1/2 when m=8
  //    std::vector<int> punc_idx = {0,  1,  2,  3,  4,  6,  9, 13, 14, 15, 18,
  //    19, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37,
  //    38, 39, 40, 42, 45, 46, 52, 53, 54, 55, 59, 60, 61, 62, 63};

  // puncture 43 bits to achieve overall rate-1/2 when m=7
  std::vector<int> punc_idx = {0,  1,  2,  3,  4,  6,  9,  13, 14, 15, 18,
                               19, 21, 22, 23, 25, 26, 27, 28, 29, 30, 31,
                               32, 33, 34, 35, 36, 37, 38, 39, 40, 42, 45,
                               46, 52, 53, 54, 55, 59, 60, 61, 62, 63};

  int INTERLEAVER_LENGTH = Km;
  int info[code.numInfoBits];
  unsigned long long combos[code.numInfoBits + 1][20];
  set_combo_lists(combos, code.numInfoBits);

  FeedbackTrellis encodingTrellis(code.k, code.n, code.v, code.numerators,
                                  code.denominator);
  //    std::cout << "numInputSymbols:" << std::endl;
  //    std::cout << encodingTrellis.getNumInputSymbols() << std::endl;
  //    std::cout << "NumOutputSymbols:" << std::endl;
  //    std::cout << encodingTrellis.getNumOutputSymbols() << std::endl;
  //    std::cout << "n states:" << std::endl;
  // 	std::cout << encodingTrellis.getNumStates() << std::endl;
  // 	std::cout << "outputs:" << std::endl;
  // 	for (int i=0; i<encodingTrellis.getNumStates(); i++){
  // 		print_int_vector(encodingTrellis.getOutputs()[i]);
  // 	}
  // 	std::cout << "next states:" << std::endl;
  // 	for (int i=0; i<encodingTrellis.getNumStates(); i++){
  // 		print_int_vector(encodingTrellis.getNextStates()[i]);
  // 	}
  // 	return;

  // interleaver
  char interleaver_file[1024]; // make sure the size is large enough
  // only K+m because half of the bits are systematic
  snprintf(interleaver_file, sizeof(interleaver_file),
           "interleaver_%d_S5_T1.dat", INTERLEAVER_LENGTH);
  unsigned short int interleaver[Km]; // Adjust size as needed
  make_file_interleaver(interleaver_file, interleaver, Km);
  unsigned short int deinterleaver[Km]; // Adjust size as needed
  for (int i = 0; i < Km; i++) {
    deinterleaver[interleaver[i]] = i;
  }

  // std::vector<int> deinterleaved_data_vector; // deinterleave example
  // for (int i = 0; i < Km*2; i++) {
  // 	deinterleaved_data_vector.push_back(interleaved_data_vector[deinterleaver[i]]);
  // }

  // store CRC dist spectra in a .txt file
  std::ofstream distSpectraFile;
  string filename = "turbo_dist_spectra_v6m8K64_enumerateMessageOnly.txt";
  distSpectraFile.open(filename);

  std::vector<int> dist_spectra; // record weight of each CRC-passing cwd
  // outer loop: CRC
  // for(int crc = 0; crc < 128; crc++){
  // 	int crc_poly = 257 + crc*2; // actual crc polynomial: 257:2:511
  for (int crc = 0; crc < 64; crc++) {
    int crc_poly = 129 + crc * 2; // actual crc polynomial:
    distSpectraFile << "Current CRC: " << crc_poly << std::endl;
    // inner loop: all message sequences
    for (int h = 1; h < 4; h++) {
      for (int i = 0; i < combos[code.numInfoBits][h]; i++) {
        setWeightIJ(combos, info, h, i, code.numInfoBits);
        vector<int> info_v(info, *(&info + 1)); // int array to vector
        // add CRC
        crc_calculation(info_v, code.crcDeg, crc_poly); // 72 bits
        // print_int_vector(info_v);
        // encode message
        std::vector<int> encodedMessage = encodingTrellis.encoder(info_v);
        // print_int_vector(encodedMessage);
        // when encoding fails (ie. no matching starting/ending state), the
        // message is returned
        if (encodedMessage.size() == info_v.size()) {
          std::cout << "fail to encode" << std::endl;
          continue;
        }
        // get rid of systematic bits
        std::vector<int> parityBits;
        for (int j = 0; j < encodedMessage.size(); j += 2) {
          parityBits.push_back(encodedMessage[j]);
        }
        // interleave message bits
        std::vector<int> interleave_info;
        for (int ii = 0; ii < Km; ii++) {
          interleave_info.push_back(info_v[interleaver[ii]]);
        }
        // encode interleaved message
        std::vector<int> encodedMessage_inter =
            encodingTrellis.encoder(interleave_info);
        // get rid of systematic bits
        std::vector<int> parityBits_inter;
        for (int j = 0; j < encodedMessage_inter.size(); j += 2) {
          parityBits_inter.push_back(encodedMessage_inter[j]);
        }
        // puncture both parity seauences & record weight
        for (int p = 0; p < parityBits.size(); p++) {
          for (int q = 0; q < punc_idx.size(); q++) {
            if (p == punc_idx[q]) {
              parityBits[p] = 0;
              parityBits_inter[p] = 0;
            }
          }
        }
        int parityWeight = std::count(parityBits.begin(), parityBits.end(), -1);
        int interleaveWeight =
            std::count(parityBits_inter.begin(), parityBits_inter.end(), -1);
        // calculate total cwd weight
        int totalWeight = h + parityWeight + interleaveWeight;
        // store in dist spectra
        distSpectraFile << totalWeight << " ";
      }
      // distSpectraFile << std::endl;
    }
    distSpectraFile << std::endl;
  }
  distSpectraFile.close();
}

void IRWEFtest(codeInformation code) {
  int Km = code.numInfoBits + code.crcDeg - 1;
  //    // puncture 44 bits to achieve overall rate-1/2 for m=8
  //    std::vector<int> punc_idx = {0,  1,  2,  3,  4,  6,  9, 13, 14, 15, 18,
  //    19, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37,
  //    38, 39, 40, 42, 45, 46, 52, 53, 54, 55, 59, 60, 61, 62, 63};

  // puncture 43 bits to achieve overall rate-1/2 when m=7
  std::vector<int> punc_idx = {0,  1,  2,  3,  4,  6,  9,  13, 14, 15, 18,
                               19, 21, 22, 23, 25, 26, 27, 28, 29, 30, 31,
                               32, 33, 34, 35, 36, 37, 38, 39, 40, 42, 45,
                               46, 52, 53, 54, 55, 59, 60, 61, 62, 63};

  int INTERLEAVER_LENGTH = Km;
  int info[INTERLEAVER_LENGTH];
  unsigned long long combos[Km + 1][20];
  set_combo_lists(combos, Km);

  FeedbackTrellis encodingTrellis(code.k, code.n, code.v, code.numerators,
                                  code.denominator);

  // interleaver
  char interleaver_file[1024]; // make sure the size is large enough
  // only K+m because half of the bits are systematic
  snprintf(interleaver_file, sizeof(interleaver_file),
           "interleaver_%d_S5_T1.dat", Km);
  unsigned short int interleaver[Km]; // Adjust size as needed
  make_file_interleaver(interleaver_file, interleaver, Km);
  unsigned short int deinterleaver[Km]; // Adjust size as needed
  for (int i = 0; i < Km; i++) {
    deinterleaver[interleaver[i]] = i;
  }

  // std::vector<int> deinterleaved_data_vector; // deinterleave example
  // for (int i = 0; i < Km*2; i++) {
  // 	deinterleaved_data_vector.push_back(interleaved_data_vector[deinterleaver[i]]);
  // }

  // store CRC dist spectra in a .txt file
  std::ofstream distSpectraFile;
  string filename = "turbo_dist_spectra_v6m7K64.txt";
  distSpectraFile.open(filename);
  // // store overall dist spectra in a .txt file
  // std::ofstream allSpectraFile;
  // string filename2 = "turbo_all_spectra_v6m8K64.txt";
  // allSpectraFile.open(filename2);

  std::vector<int> dist_spectra; // record weight of each CRC-passing cwd
  // outer loop: CRC
  for (int crc = 0; crc < 64; crc++) {
    int crc_poly = 129 + crc * 2; // actual crc polynomial: 257:2:511
    distSpectraFile << "Current CRC: " << crc_poly << std::endl;
    // inner loop: all message sequences
    for (int h = 1; h <= 5; h++) {
      // allSpectraFile << "Current message weight: " << h << std::endl;
      int min_h = INT_MAX;
      int counter_h = 0;
      for (int i = 0; i < combos[Km][h]; i++) {
        setWeightIJ(combos, info, h, i, Km);
        vector<int> info_v(info, *(&info + 1)); // int array to vector
        // print_int_vector(info_v);
        if (crc_check(info_v, code.crcDeg, crc_poly)) {
          // counter_h ++ ;
          // encode message
          std::vector<int> encodedMessage = encodingTrellis.encoder(info_v);
          // when encoding fails (ie. no matching starting/ending state), the
          // message is returned
          if (encodedMessage.size() == info_v.size()) {
            continue;
            std::cout << "failed encoding" << std::endl;
          }
          // get rid of systematic bits
          std::vector<int> parityBits;
          for (int j = 0; j < encodedMessage.size(); j += 2) {
            parityBits.push_back(encodedMessage[j]);
          }
          // interleave message bits
          std::vector<int> interleave_info;
          for (int ii = 0; ii < Km; ii++) {
            interleave_info.push_back(info_v[interleaver[ii]]);
          }
          // encode interleaved message
          std::vector<int> encodedMessage_inter =
              encodingTrellis.encoder(interleave_info);
          // get rid of systematic bits
          std::vector<int> parityBits_inter;
          for (int j = 0; j < encodedMessage_inter.size(); j += 2) {
            parityBits_inter.push_back(encodedMessage_inter[j]);
          }
          // puncture both parity seauences & record weight
          for (int p = 0; p < parityBits.size(); p++) {
            for (int q = 0; q < punc_idx.size(); q++) {
              if (p == punc_idx[q]) {
                parityBits[p] = 0;
                parityBits_inter[p] = 0;
              }
            }
          }
          int parityWeight =
              std::count(parityBits.begin(), parityBits.end(), -1);
          int interleaveWeight =
              std::count(parityBits_inter.begin(), parityBits_inter.end(), -1);
          // calculate total cwd weight
          int totalWeight = h + parityWeight + interleaveWeight;
          if (totalWeight < min_h) {
            min_h = totalWeight;
            counter_h = 1;
          } else if (totalWeight == min_h) {
            counter_h++;
          }
          // store in dist spectra
          // distSpectraFile << totalWeight << " ";
          // allSpectraFile << totalWeight << " ";
        }
      }
      // allSpectraFile << std::endl;
      distSpectraFile << "h: " << h << " dmin: " << min_h
                      << " nearest neighbors: " << counter_h << std::endl;
    }
    // distSpectraFile << std::endl;
  }
  distSpectraFile.close();
  // allSpectraFile.close();
}

int comboTest(codeInformation code) {
  int Km = code.numInfoBits + code.crcDeg - 1; // K=64, m=8
                                               // int h = 3;
  int INTERLEAVER_LENGTH = Km;
  int info[INTERLEAVER_LENGTH];
  unsigned long long combos[Km + 1][20];
  set_combo_lists(combos, Km);

  FeedforwardTrellis encodingTrellis =
      FeedforwardTrellis(code.k, code.n, code.v, code.numerators);
  std::vector<int> puncpat = {1, 1, 1, 0}; // puncture rate-1/2 to rate-2/3
  int maxNumTrials = 1e5;
  int listSize = 1e5;
  LowRateListDecoder listDecoder(encodingTrellis, listSize, code.crcDeg,
                                 code.crc);

  /*
// interleaver
char interleaver_file[1024]; // make sure the size is large enough
// only K+m because half of the bits are systematic
  snprintf(interleaver_file, sizeof(interleaver_file),
"interleaver_%d_S7_T2.dat", Km); unsigned short int interleaver[Km]; // Adjust
size as needed make_file_interleaver(interleaver_file, interleaver, Km);
  unsigned short int deinterleaver[Km]; // Adjust size as needed
  for (int i = 0; i < Km; i++) {
          deinterleaver[interleaver[i]] = i;
  }
  */

  // std::vector<int> deinterleaved_data_vector; // deinterleave example
  // for (int i = 0; i < Km*2; i++) {
  // 	deinterleaved_data_vector.push_back(interleaved_data_vector[deinterleaver[i]]);
  // }
  int optimal_crc;
  // for(int h = 1; h < 5; h++) {
  for (int h = 1; h < 3; h++) {
    int max_dmin_prev_h = 0;
    std::vector<int> d_mins; // record for each CRC
    std::vector<int> A_dmins;
    for (int i = 0; i < combos[Km][h]; i++) {
      setWeightIJ(combos, info, h, i, Km);
      vector<int> info_v(info, *(&info + 1)); // int array to vector
      for (int crc = 0; crc < 128; crc++) {
        int crc_poly = 257 + crc * 2; // actual crc polynomial: 257:2:511
        d_mins.push_back(INT_MAX);
        A_dmins.push_back(0);
        /*
// see if any weight-1 or 2 sequences pass any CRC
        if (crc_check(info_v, code.crcDeg, crc)){
                std::cout << "weight-" << h << " pass crc: " << crc_poly <<
std::endl;
        }
        */

        // check CRC[crc]
        if (crc_check(info_v, code.crcDeg, crc)) {
          // if any passes encode with code 1 and get weight
          std::vector<int> encodedMessage = encodingTrellis.encoder(info_v);
          int cwd_weight =
              std::accumulate(encodedMessage.begin(), encodedMessage.end(), 0);

          // interleave only the nonsystematic bits
          // TODO: CHANGE THIS PART TO ONLY INTERLEAVE ON EVERY OTHER BIT
          std::vector<int> interleaved_info;
          for (int i = 0; i < Km; i++) {
            // interleaved_info.push_back(info_v[interleaver[i]]);
          }
          // check crc of interleaved  message
          if (crc_check(interleaved_info, code.crcDeg, crc)) {
            // encode with code 1 the interleaved info and get weight
            std::vector<int> interleaved_encodedMessage =
                encodingTrellis.encoder(interleaved_info);
            int interleaved_cwd_weight =
                std::accumulate(interleaved_encodedMessage.begin(),
                                interleaved_encodedMessage.end(), 0);
            // compare and update weight
            int total_weight = h + cwd_weight + interleaved_cwd_weight;
            if (total_weight < d_mins[crc]) {
              d_mins[crc] = total_weight;
              A_dmins[crc] = 1;
            } else if (total_weight == d_mins[crc]) {
              A_dmins[crc]++;
            }
            // When h + 2 * d min of main code > any dmin for any crc(j) j=0,1,
            // .. crc Count, then remove these CRCs. int max_dmin =
            // *max_element(std::begin(d_mins), std::end(d_mins));
            if (h + 2 * cwd_weight > max_dmin_prev_h) {
              d_mins[crc] = -1;
              A_dmins[crc] = -1;
            }
          }
        }
      }
    }
    //   max_dmin_prev_h = *max_element(std::begin(d_mins), std::end(d_mins));
    //   if (h==5){
    // 	optimal_crc = *find(d_mins.begin(), d_mins.end(), max_dmin_prev_h);
    //   }
  }
  //    std::cout << "optimal CRC is: " << optimal_crc << std::endl;
  exit(0);
  return 0;
}

void set_combo_lists(unsigned long long combo_list[][20], int K) {
  combo_list[0][0] = 1;
  for (unsigned int i = 1; i <= K; i++) {
    combo_list[i][0] = 1;
    if (i < 20)
      combo_list[i][i] = 1;
    for (unsigned int j = 1; j < i & j < 20; j++)
      combo_list[i][j] = combo_list[i - 1][j] + combo_list[i - 1][j - 1];
  }
}

void setWeightIJ(const unsigned long long combos[][20], int informaton[], int h,
                 unsigned long long int item, int Km) {
  if (h >= 20) {
    printf("ERROR: h is too large, h: %d, max 19\n", h);
    exit(1);
  }
  for (int i = 0; i < Km; i++) {
    informaton[i] = 0;
  }
  if (h == 0)
    return;
  unsigned int i = (0), j = (0);
  unsigned long long int current = item;
  while (j < h) {
    if (current < combos[Km - i - 1][h - j - 1]) {
      informaton[i] = 1;
      j++;
    } else {
      current = current - combos[Km - i - 1][h - j - 1];
    }
    i++;
  }
}

void turbo_ELF(codeInformation code) {
  // for m=12, puncture 12 bits
  std::vector<int> punc_idx = {39, 30, 37, 33, 32, 36,
                               27, 35, 40, 38, 34, 28}; // 477

  vector<double> SNR = {4.5};

  std::cout << "running turbo-ELF performance simulations" << endl;
  srand(time(NULL));
  string filename = "";

  filename += "Turbo-ELF rate" + to_string(code.k) + "-" + to_string(code.n) +
              ",v" + to_string(code.v) + ",crcdeg" + to_string(code.crcDeg) +
              ",crc" + to_string(code.crc) + ",infolen" +
              to_string(code.numInfoBits) + ".txt";
  std::cout << "working on code " << filename.substr(0, filename.size() - 4)
            << endl;
  if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
    std::cout << "invalid msg + crc length" << endl;
    return;
  }

  std::ofstream outputFile;
  outputFile.open(filename, fstream::app);

  // use rate-2/3 RSC from Lin & Costello
  FeedbackTrellis encodingTrellis(code.k, code.n, code.v, code.numerators,
                                  code.denominator);
  std::cout << "feedback trellis construction complete" << endl;
  DualTrellis decodingTrellis(code.hMatrix);

  int n = 140;                 // or any value you want
  char interleaver_file[1024]; // make sure the size is large enough
  snprintf(interleaver_file, sizeof(interleaver_file),
           "interleaver_%d_S7_T2.dat", n);
  unsigned short int interleaver[n]; // Adjust size as needed
  make_file_interleaver(interleaver_file, interleaver, n);

  unsigned short int deinterleaver[n]; // Adjust size as needed
  for (int i = 0; i < n; i++) {
    deinterleaver[interleaver[i]] = i;
  }

  // interleaver example
  std::vector<int> data_vector; // 140 bits
  for (int i = 1; i < n + 1; i++) {
    data_vector.push_back(i);
  }
  std::vector<int> interleaved_data_vector;
  for (int i = 0; i < n; i++) {
    interleaved_data_vector.push_back(data_vector[interleaver[i]]);
  }
  for (int i = 0; i < n; i++) {
    printf("%d\n", interleaved_data_vector[i]);
  }

  std::vector<int> deinterleaved_data_vector;
  for (int i = 0; i < n; i++) {
    deinterleaved_data_vector.push_back(
        interleaved_data_vector[deinterleaver[i]]);
  }
  for (int i = 0; i < n; i++) {
    printf("%d\n", deinterleaved_data_vector[i]);
  }

  int maxNumTrials = 1e5;
  int listSize = 1e5;
  ListDecoder listDecoder(decodingTrellis, listSize, code.crcDeg, code.crc);

  std::cout << "beginning simulations" << endl;
  // simulate the comms system
  for (int i = 0; i < SNR.size(); i++) {
    std::cout << "beginning simulations at SNR " << SNR[i]
              << "____________________" << std::endl;
    int numTrials = 0;
    int numTBErrors = 0;
    int listsize_sum = 0;
    while (numTrials < maxNumTrials || numTBErrors < 200) {
      double snr = SNR[i];

      std::vector<int> originalMessage;
      for (int i = 0; i < code.numInfoBits; i++) {
        originalMessage.push_back(rand() % 2);
      }
      // compute the CRC for the systematic bits
      crc_calculation(originalMessage, code.crcDeg, code.crc); // 128+m bits

      std::vector<int> interleavedMessage;
      for (int i = 0; i < originalMessage.size(); i++) {
        interleavedMessage.push_back(
            originalMessage[interleaver[i]]); // original message is 140 bits
      }

      // encode the message
      std::vector<int> encodedMessage =
          encodingTrellis.encoder(originalMessage);
      std::vector<int> encodedInterleavedMessage =
          encodingTrellis.encoder(interleavedMessage);
      std::cout << "finished encoding" << std::endl;

      // extract parity bits and systematic bits
      std::vector<int> parity_bits;
      std::vector<int> parity_interleaved_bits;
      std::vector<int> systematic_bits;
      for (int i = 0; i < encodedMessage.size() / code.n; i++) {
        parity_bits.push_back(
            encodedMessage[i * code.n]); // the first bit is parity
        parity_interleaved_bits.push_back(
            encodedInterleavedMessage[i * code.n]); // the first bit is parity
        systematic_bits.push_back(encodedMessage[i * code.n + 1]);
        systematic_bits.push_back(encodedMessage[i * code.n + 2]);
      }

      // add noise
      std::vector<double> SystematicOutput = addNoise(systematic_bits, snr);
      std::vector<double> ParityOutput = addNoise(parity_bits, snr);
      std::vector<double> ParityInterleavedOutput =
          addNoise(parity_interleaved_bits, snr);

      // puncture the parity bits
      for (int index = 0; index < punc_idx.size(); index++) {
        parity_bits[punc_idx[index]] = 0;
        parity_interleaved_bits[punc_idx[index]] = 0;
      } // parity_bits is (128+m)/2 bits

      // systemmatic interleaved output
      std::vector<double> ReceivedOutput;   // top rail
      std::vector<double> ReceivedPiOutput; // bottom rail
      for (int i = 0; i < code.numInfoBits / 2; i++) {
        ReceivedOutput.push_back(parity_bits[i]);          // top to top rail
        ReceivedOutput.push_back(SystematicOutput[2 * i]); // middle to top rail
        ReceivedOutput.push_back(
            SystematicOutput[2 * i + 1]); // middle to top rail
        ReceivedPiOutput.push_back(
            parity_interleaved_bits[i]); // bottom to bottom rail
        ReceivedPiOutput.push_back(
            SystematicOutput[interleaver[2 * i]]); // middle to bottom rail
        ReceivedPiOutput.push_back(
            SystematicOutput[interleaver[2 * i + 1]]); // middle to bottom rail
      }

      ListDecoder::messageInformation decodedMessage =
          listDecoder.punctured_highrate_turbo_decoder(ReceivedOutput,
                                                       punc_idx);
      // deinterleave before checking crc
      ListDecoder::messageInformation decodedPiMessage =
          listDecoder.punctured_highrate_turbo_decoder(ReceivedPiOutput,
                                                       punc_idx);

      if (decodedMessage.message != originalMessage) {
        numTBErrors++;
      }
      listsize_sum += decodedMessage.listSize;
      numTrials++;
      if (numTrials % 10000 == 0) {
        std::cout << "currently at " << numTrials << " trials with "
                  << numTBErrors << " errors with TB decoding" << std::endl;
      }
    }
    outputFile << "current SNR: " << SNR[i]
               << ". FER: " << (double)numTBErrors / (double)numTrials
               << ". list size: " << listsize_sum / (double)numTrials
               << std::endl;
    if ((double)numTBErrors / (double)numTrials <
        1E-4) { // run until the SNR point with FER<1e-4
      outputFile.close();
      return;
    }
  }
}

void find_optimal_crc_turbo(codeInformation code) {
  DualTrellis trell = DualTrellis(code.hMatrix);
  ListDecoder decoder = ListDecoder(trell, -1, code.crcDeg, code.crc);
  FeedforwardTrellis encTrell =
      FeedforwardTrellis(code.k, code.n, code.v, code.numerators);

  vector<double> allZerosMessage;
  for (int i = 0; i < (code.numInfoBits + code.crcDeg - 1) * 2; i++) {
    allZerosMessage.push_back(1);
  }

  decoder.FindBestCRC_turbo(trell, allZerosMessage, code);
}

int make_file_interleaver(char interleaver_file[],
                          unsigned short int interleaver[], int n) {
  FILE *fp;
  if ((fp = fopen(interleaver_file, "r")) == NULL) {
    printf("Error opening interleaver file for reading.\n");
    printf("%s\n", interleaver_file);
    exit(1); /* exit program with error condition  */
  }
  printf("Reading interleaver from file.\n");
  for (int i = 0; i < n; i++) {
    if (fscanf(fp, "%hu", &interleaver[i])) {
    };
    /* Translate from 1..K to 0..(K-1) */
    //       interleaver[i]--;
  }

  printf("%s\n", interleaver_file);
  fclose(fp);
  return 0;
}

// adds noise at the given snr to each bit, based on a normal distribution
std::vector<double> addNoise(std::vector<int> encodedMsg, double SNR) {
  std::vector<double> noisyMsg;

  /*
  // the below lines break the noise generator for some compilers
  std::random_device rd;
  std::default_random_engine generator;
  generator.seed( rd() ); //Now this is seeded differently each time.
  */
  double variance = pow(10.0, -SNR / 10.0);
  double sigma = sqrt(variance);
  normal_distribution<double> distribution(0.0, sigma);

  // cout << variance << endl;

  for (int i = 0; i < encodedMsg.size(); i++) {
    noisyMsg.push_back(encodedMsg[i] + distribution(generator));
  }
  return noisyMsg;
}
