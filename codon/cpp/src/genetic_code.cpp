#include "genetic_code.h"

namespace daa {
namespace genetic {
namespace {

const char* kBaseChars = "TCAG";  // internal base encoding: T=0, C=1, A=2, G=3

int baseValue(char ch) {
  switch (ch) {
    case 'T': return 0;
    case 'C': return 1;
    case 'A': return 2;
    case 'G': return 3;
    default:  return -1;
  }
}

}  // namespace

int aaIndex(char letter) {
  switch (letter) {
    case 'A': return 0; case 'R': return 1; case 'N': return 2; case 'D': return 3;
    case 'C': return 4; case 'Q': return 5; case 'E': return 6; case 'G': return 7;
    case 'H': return 8; case 'I': return 9; case 'L': return 10; case 'K': return 11;
    case 'M': return 12; case 'F': return 13; case 'P': return 14; case 'S': return 15;
    case 'T': return 16; case 'W': return 17; case 'Y': return 18; case 'V': return 19;
    default: return -1;
  }
}

int aaIndexForCodon(int codonIndex) {
  // Standard genetic code table in codon-index order (bases ordered T, C, A, G).
  // Characters are amino acid letters; '*' marks stop codons.
  static const std::string kTable =
      "FFLLLLLLIIIMVVVVSSSSPPPPTTTTAAAAYY**HHQQNNKKDDEECCWRRRRSSRRGGGG";
  const char c = kTable[codonIndex];
  if (c == '*') return -1;
  return aaIndex(c);
}

int codonIndexForString(const std::string& codon) {
  if (codon.size() != 3) return -1;
  int b[3];
  for (int i = 0; i < 3; ++i) {
    const int v = baseValue(codon[i]);
    if (v < 0) return -1;
    b[i] = v;
  }
  return ((b[0] * 4) + b[1]) * 4 + b[2];
}

void codonString(int codonIndex, char out[3]) {
  const int b1 = codonIndex / 16;
  const int rem = codonIndex % 16;
  const int b2 = rem / 4;
  const int b3 = rem % 4;
  out[0] = kBaseChars[b1];
  out[1] = kBaseChars[b2];
  out[2] = kBaseChars[b3];
}

int degeneracy(int aa) { return (int)codonsForAa(aa).size(); }

int senseCodonCount(int aa) { return (int)codonsForAa(aa).size(); }

int gcCount(int codonIndex) {
  char buf[4];
  codonString(codonIndex, buf);
  int g = 0;
  for (int i = 0; i < 3; ++i) {
    if (buf[i] == 'G' || buf[i] == 'C') ++g;
  }
  return g;
}

const std::vector<int>& codonsForAa(int aa) {
  static std::vector<std::vector<int>> cache(kNumAminoAcids);
  static bool ready = false;
  if (!ready) {
    for (int i = 0; i < 64; ++i) {
      const int a = aaIndexForCodon(i);
      if (a >= 0) cache[a].push_back(i);
    }
    ready = true;
  }
  return cache[aa];
}

const std::vector<int>& codonsForAaReduced(int aa) {
  static std::vector<std::vector<int>> cache(kNumAminoAcids);
  static bool ready = false;
  if (!ready) {
    for (int a = 0; a < kNumAminoAcids; ++a) {
      cache[a].clear();
      int taken = 0;
      for (int c = 0; c < 64; ++c) {
        const int a2 = aaIndexForCodon(c);
        if (a2 != a) continue;
        if (taken >= 2) break;
        cache[a].push_back(c);
        ++taken;
      }
    }
    ready = true;
  }
  return cache[aa];
}

}  // namespace genetic
}  // namespace daa
