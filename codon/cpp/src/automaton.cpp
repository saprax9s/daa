#include "automaton.h"

namespace daa {
namespace automaton {

int AhoCorasick::build(const motifs::MotifSet& ms) {
  nodes.clear();
  nodes.emplace_back();
  for (const std::string& p : ms.allPatterns()) {
    int state = 0;
    for (char ch : p) {
      const int b = baseIndex(ch);
      if (nodes[state].next[b] == 0) {
        nodes.emplace_back();
        nodes[state].next[b] = static_cast<int>(nodes.size()) - 1;
      }
      state = nodes[state].next[b];
    }
    nodes[state].out = true;
  }
  // BFS: failure links + completed goto transitions.
  std::vector<int> queue;
  queue.reserve(nodes.size());
  for (int b = 0; b < 4; ++b) {
    const int v = nodes[0].next[b];
    if (v != 0) {
      nodes[v].fail = 0;
      queue.push_back(v);
    } else {
      nodes[0].next[b] = 0;
    }
  }
  for (std::size_t head = 0; head < queue.size(); ++head) {
    const int u = queue[head];
    for (int b = 0; b < 4; ++b) {
      const int v = nodes[u].next[b];
      if (v != 0) {
        nodes[v].fail = nodes[nodes[u].fail].next[b];
        nodes[v].out = nodes[v].out || nodes[nodes[v].fail].out;
        queue.push_back(v);
      } else {
        nodes[u].next[b] = nodes[nodes[u].fail].next[b];
      }
    }
  }
  return static_cast<int>(nodes.size());
}

std::pair<int, bool> AhoCorasick::step(int state, int base) const {
  const int ns = nodes[state].next[base];
  return std::make_pair(ns, nodes[ns].out);
}

bool AhoCorasick::search(const std::string& dna) const {
  int state = 0;
  for (char ch : dna) {
    const std::pair<int, bool> pr = step(state, baseIndex(ch));
    state = pr.first;
    if (pr.second) return true;
  }
  return false;
}

int AhoCorasick::numStates() const { return static_cast<int>(nodes.size()); }

int baseIndex(char ch) {
  switch (ch) {
    case 'A': return 2;
    case 'C': return 1;
    case 'G': return 3;
    case 'T': return 0;
    default: return 0;
  }
}

}  // namespace automaton
}  // namespace daa
