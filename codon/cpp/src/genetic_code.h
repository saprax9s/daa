#pragma once

#include <string>
#include <vector>

namespace daa {
namespace genetic {

constexpr int kNumAminoAcids = 20;

// Amino acid indices:
// A=0 R=1 N=2 D=3 C=4 Q=5 E=6 G=7 H=8 I=9 L=10 K=11 M=12 F=13 P=14 S=15 T=16 W=17 Y=18 V=19

// Returns the amino acid index for a letter, or -1 for stop/invalid codes.
int aaIndex(char letter);

// Returns the amino acid index encoded by a codon index (0..63), or -1 for stop.
int aaIndexForCodon(int codonIndex);

// Maps a codon string (bases in internal order T/C/A/G) to its index, or -1 if invalid.
int codonIndexForString(const std::string& codon);

// Writes the codon string (internal base order T/C/A/G) into out[0..2].
void codonString(int codonIndex, char out[3]);

// Number of codons encoding the amino acid (1, 2, 3, 4 or 6).
int degeneracy(int aa);

// Number of sense codons for an amino acid (same as degeneracy; stops excluded).
int senseCodonCount(int aa);

// GC content of a codon as a count of G/C nucleotides (0..3).
int gcCount(int codonIndex);

// All sense codon indices encoding the amino acid, ascending by codon index.
const std::vector<int>& codonsForAa(int aa);

// Reduced code: exactly 2 codons per amino acid (1 where only 1 exists),
// fixed deterministic choice (lowest codon indices).
const std::vector<int>& codonsForAaReduced(int aa);

}  // namespace genetic
}  // namespace daa
