#pragma once

#include <array>
#include <string>
#include <utility>
#include <vector>

#include "motifs.h"

namespace daa {
namespace automaton {

struct AhoCorasick {
  struct Node {
    std::array<int, 4> next;
    int fail;
    bool out;
    Node() : next{0, 0, 0, 0}, fail(0), out(false) {}
  };

  std::vector<Node> nodes;

  // Builds the trie with BFS failure links and completed goto transition table.
  // Returns the number of automaton states.
  int build(const motifs::MotifSet& ms);

  // Returns (newState, bad); bad means the new state is terminal (motif hit).
  std::pair<int, bool> step(int state, int base) const;

  // True if any motif occurs in dna (including boundary-spanning occurrences).
  bool search(const std::string& dna) const;

  int numStates() const;
};

int baseIndex(char ch);

}  // namespace automaton
}  // namespace daa
