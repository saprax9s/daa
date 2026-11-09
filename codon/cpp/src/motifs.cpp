#include "motifs.h"

#include <algorithm>

namespace daa {
namespace motifs {
namespace {

char complementChar(char ch) {
  switch (ch) {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';
    case 'R': return 'Y';
    case 'Y': return 'R';
    case 'S': return 'S';
    case 'W': return 'W';
    case 'K': return 'M';
    case 'M': return 'K';
    case 'B': return 'V';
    case 'V': return 'B';
    case 'D': return 'H';
    case 'H': return 'D';
    case 'N': return 'N';
    default: return 'N';
  }
}

}  // namespace

void MotifSet::add(const std::string& pattern, const std::string& source) {
  patterns.push_back(pattern);
  sources.push_back(source);
  const std::string rev = MotifSet::reverseComplement(pattern);
  if (rev != pattern) {
    patterns.push_back(rev);
    sources.push_back(source);
  }
}

void MotifSet::addIupac(const std::string& pattern, const std::string& source) {
  const std::vector<std::string> expanded = MotifSet::expandIupac(pattern);
  for (const std::string& e : expanded) {
    add(e, source);
  }
}

std::string MotifSet::reverseComplement(const std::string& s) {
  std::string out = s;
  for (char& c : out) c = complementChar(c);
  std::reverse(out.begin(), out.end());
  return out;
}

static std::vector<std::string> iupacBases(char ch) {
  switch (ch) {
    case 'A': return {"A"};
    case 'C': return {"C"};
    case 'G': return {"G"};
    case 'T': case 'U': return {"T"};
    case 'R': return {"A", "G"};
    case 'Y': return {"C", "T"};
    case 'S': return {"C", "G"};
    case 'W': return {"A", "T"};
    case 'K': return {"G", "T"};
    case 'M': return {"A", "C"};
    case 'B': return {"C", "G", "T"};
    case 'D': return {"A", "G", "T"};
    case 'H': return {"A", "C", "T"};
    case 'V': return {"A", "C", "G"};
    case 'N': return {"A", "C", "G", "T"};
    default: return {std::string(1, ch)};
  }
}

std::vector<std::string> MotifSet::expandIupac(const std::string& pattern) {
  std::vector<std::string> results = {""};
  for (char ch : pattern) {
    const std::vector<std::string> choices = iupacBases(ch);
    std::vector<std::string> next;
    next.reserve(results.size() * choices.size());
    for (const std::string& prefix : results) {
      for (const std::string& base : choices) {
        next.push_back(prefix + base);
      }
    }
    results = std::move(next);
  }
  return results;
}

bool naiveMotifSearch(const MotifSet& ms, const std::string& dna) {
  for (const std::string& p : ms.allPatterns()) {
    if (dna.find(p) != std::string::npos) return true;
  }
  return false;
}

}  // namespace motifs
}  // namespace daa
