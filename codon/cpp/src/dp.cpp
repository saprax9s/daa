#include "dp.h"

#include <algorithm>
#include <map>
#include <tuple>
#include <vector>

namespace daa {
namespace dp {

namespace {
using DpKey = std::tuple<int, int, std::vector<int>, int>;
struct Entry {
  DpKey key;
  int64_t score;
  int parent;
  int codon;
};
}  // namespace

DpResult solve(const daa::automaton::AhoCorasick& ac,
              const daa::scores::ScoreTable& st,
              const daa::motifs::MotifSet& ms,
              const daa::dp::DpInstance& inst,
              const daa::dp::DpFlags& flags) {
  (void)ms;
  DpResult res;
  const int n = inst.n;

  // Lower bound on remaining GC additions (for GC-bound pruning).
  std::vector<int64_t> minRem(n + 1, 0);
  for (int i = n - 1; i >= 0; --i) {
    int64_t m = 0;
    bool first = true;
    for (int c : inst.allowed[i]) {
      const int g = daa::genetic::gcCount(c);
      if (first || g < m) {
        m = g;
        first = false;
      }
    }
    minRem[i] = minRem[i + 1] + m;
  }

  std::vector<std::vector<Entry>> layers;
  std::vector<Entry> cur;
  {
    DpKey k(-1, 0, std::vector<int>(), 0);
    cur.push_back(Entry{k, 0, -1, -1});
  }
  layers.push_back(cur);

  const bool doGcPruning = flags.gcPruning || flags.dominance;
  const bool doDominance = flags.dominance;

  for (int i = 0; i < n; ++i) {
    std::vector<Entry> nxt;
    std::map<DpKey, int> index;
    for (int ei = 0; ei < static_cast<int>(cur.size()); ++ei) {
      const Entry& e = cur[ei];
      const int last = std::get<0>(e.key);
      const int q = std::get<1>(e.key);
      const std::vector<int> w = std::get<2>(e.key);
      const int gc = std::get<3>(e.key);
      const int64_t baseScore = e.score;
      for (int c : inst.allowed[i]) {
        int qq = q;
        bool bad = false;
        char buf[4];
        daa::genetic::codonString(c, buf);
        for (int b = 0; b < 3; ++b) {
          const std::pair<int, bool> pr = ac.step(qq, daa::automaton::baseIndex(buf[b]));
          qq = pr.first;
          if (pr.second) {
            bad = true;
            break;
          }
        }
        if (bad) continue;
        const int codonGc = daa::genetic::gcCount(c);
        const int newGc = gc + codonGc;
        std::vector<int> w2 = w;
        if (inst.windowK > 0) {
          w2.push_back(codonGc);
          if (static_cast<int>(w2.size()) >= inst.windowK) {
            int sum = 0;
            for (int j = static_cast<int>(w2.size()) - inst.windowK; j < static_cast<int>(w2.size()); ++j) {
              sum += w2[j];
            }
            if (sum < inst.windowLo || sum > inst.windowHi) continue;
          }
          if (static_cast<int>(w2.size()) > inst.windowK - 1) {
            w2.erase(w2.begin());
          }
        }
        if (doGcPruning) {
          if (static_cast<int64_t>(newGc) > inst.hi) continue;
          if (static_cast<int64_t>(newGc) + minRem[i + 1] < inst.lo) continue;
        }
        DpKey nk(c, qq, w2, newGc);
        const int64_t add = st.codonScore[c] + (last >= 0 ? st.pairScore[last][c] : static_cast<int64_t>(0));
        const int64_t sc = baseScore + add;
        auto it = index.find(nk);
        if (it == index.end()) {
          index[nk] = static_cast<int>(nxt.size());
          nxt.push_back(Entry{nk, sc, ei, c});
        } else {
          Entry& ex = nxt[it->second];
          if (sc > ex.score) {
            ex.score = sc;
            ex.parent = ei;
            ex.codon = c;
          }
        }
      }
    }
    if (doDominance) {
      // Dominance: within group (last codon, automaton state, window vector),
      // a state with GC <= and score >= another dominates it. Only valid when lo <= 0.
      std::vector<char> drop(nxt.size(), 0);
      std::map<std::tuple<int, int, std::vector<int>>, std::vector<int>> groups;
      for (int i = 0; i < static_cast<int>(nxt.size()); ++i) {
        const DpKey& k = nxt[i].key;
        groups[std::make_tuple(std::get<0>(k), std::get<1>(k), std::get<2>(k))].push_back(i);
      }
      for (auto& g : groups) {
        const std::vector<int>& v = g.second;
        for (int i : v) {
          if (drop[i]) continue;
          for (int j : v) {
            if (i == j || drop[j]) continue;
            const Entry& a = nxt[i];
            const Entry& b = nxt[j];
            const int gcA = std::get<3>(a.key);
            const int gcB = std::get<3>(b.key);
            if (gcB <= gcA && b.score >= a.score) {
              drop[i] = 1;
              break;
            }
          }
        }
      }
      std::vector<Entry> compacted;
      compacted.reserve(nxt.size());
      for (int i = 0; i < static_cast<int>(nxt.size()); ++i) {
        if (!drop[i]) compacted.push_back(nxt[i]);
      }
      nxt = compacted;
    }
    layers.push_back(nxt);
    cur = nxt;
  }

  uint64_t total = 0;
  uint64_t peak = 0;
  for (const auto& lyr : layers) {
    total += lyr.size();
    if (lyr.size() > peak) peak = lyr.size();
  }
  res.totalStates = total;
  res.peakStates = peak;

  const std::vector<Entry>& finalLayer = layers.back();
  int bestIdx = -1;
  for (int i = 0; i < static_cast<int>(finalLayer.size()); ++i) {
    const int gc = std::get<3>(finalLayer[i].key);
    if (gc < inst.lo || gc > inst.hi) continue;
    if (bestIdx < 0 || finalLayer[i].score > finalLayer[bestIdx].score) bestIdx = i;
  }
  if (bestIdx < 0) {
    res.feasible = false;
    res.score = 0;
    res.codons.clear();
    return res;
  }
  res.feasible = true;
  res.score = finalLayer[bestIdx].score;
  std::vector<int> codons;
  int idx = bestIdx;
  for (int l = n; l >= 1; --l) {
    const Entry& e = layers[l][idx];
    codons.push_back(e.codon);
    idx = e.parent;
  }
  std::reverse(codons.begin(), codons.end());
  res.codons = codons;
  return res;
}

}  // namespace dp
}  // namespace daa
