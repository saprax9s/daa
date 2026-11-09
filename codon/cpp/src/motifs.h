#pragma once

#include <string>
#include <vector>

namespace daa {
namespace motifs {

struct Motif {
  std::string pattern;
  std::string source;
};

struct MotifSet {
  std::vector<std::string> patterns;
  std::vector<std::string> sources;

  // Adds a concrete motif and its reverse complement (double-stranded).
  void add(const std::string& pattern, const std::string& source);

  // Adds every concrete string of an IUPAC pattern and its reverse complement.
  void addIupac(const std::string& pattern, const std::string& source);

  const std::vector<std::string>& allPatterns() const { return patterns; }

  // Complement of one nucleotide or IUPAC code, then reverse the string.
  static std::string reverseComplement(const std::string& s);

  // Expands an IUPAC pattern (A/C/G/T/R/Y/S/W/K/M/B/D/H/V/N) into concrete strings.
  static std::vector<std::string> expandIupac(const std::string& pattern);
};

// Independent checker: naive substring search for any motif (both strands kept
// in the MotifSet).
bool naiveMotifSearch(const MotifSet& ms, const std::string& dna);

}  // namespace motifs
}  // namespace daa
