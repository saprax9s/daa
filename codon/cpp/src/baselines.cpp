#include "baselines.h"

#include <cmath>
#include <functional>

namespace daa {
namespace baselines {

namespace {

bool windowOk(const daa::dp::DpInstance& inst, const std::vector<int>& codons) {
  if (inst.windowK <= 0) return true;
  const int k = inst.windowK;
  const int n = static_cast<int>(codons.size());
  for (int s = 0; s + k <= n; ++s) {
    int sum = 0;
    for (int j = s; j < s + k; ++j) sum += daa::genetic::gcCount(codons[j]);
    if (sum < inst.windowLo || sum > inst.windowHi) return false;
  }
  return true;
}

bool sequenceValid(const daa::dp::DpInstance& inst,
                   const daa::motifs::MotifSet& ms,
                   const std::vector<int>& codons) {
  std::string dna;
  int gc = 0;
  for (int c : codons) {
    char buf[4];
    daa::genetic::codonString(c, buf);
    dna += buf;
    gc += daa::genetic::gcCount(c);
  }
  if (daa::motifs::naiveMotifSearch(ms, dna)) return false;
  if (gc < inst.lo || gc > inst.hi) return false;
  if (!windowOk(inst, codons)) return false;
  return true;
}

}  // namespace

BruteResult bruteForce(const daa::dp::DpInstance& inst,
                       const daa::scores::ScoreTable& st,
                       const daa::motifs::MotifSet& ms) {
  BruteResult res;
  res.feasible = false;
  res.score = 0;
  res.evaluated = 0;
  std::vector<int> codons;
  std::function<void(int)> rec = [&](int pos) {
    if (pos == inst.n) {
      std::string dna;
      int gc = 0;
      for (int c : codons) {
        char buf[4];
        daa::genetic::codonString(c, buf);
        dna += buf;
        gc += daa::genetic::gcCount(c);
      }
      if (daa::motifs::naiveMotifSearch(ms, dna)) return;
      if (gc < inst.lo || gc > inst.hi) return;
      if (!windowOk(inst, codons)) return;
      const int64_t sc = daa::scores::sequenceScore(codons, st);
      ++res.evaluated;
      if (!res.feasible || sc > res.score) {
        res.feasible = true;
        res.score = sc;
        res.codons = codons;
      }
      return;
    }
    for (int c : inst.allowed[pos]) {
      codons.push_back(c);
      rec(pos + 1);
      codons.pop_back();
    }
  };
  rec(0);
  return res;
}

GreedyResult greedy(const daa::dp::DpInstance& inst,
                    const daa::scores::ScoreTable& st,
                    const daa::automaton::AhoCorasick& ac,
                    const daa::motifs::MotifSet& ms,
                    bool gcSteered) {
  GreedyResult res;
  res.feasible = false;
  res.score = 0;
  int q = 0;
  int last = -1;
  int gc = 0;
  std::vector<int> codons;
  for (int i = 0; i < inst.n; ++i) {
    bool found = false;
    int bestC = -1;
    int64_t bestAdd = 0;
    int bestNewGc = 0;
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
      const int64_t add = st.codonScore[c] + (last >= 0 ? st.pairScore[last][c] : static_cast<int64_t>(0));
      const int newGc = gc + daa::genetic::gcCount(c);
      if (!found) {
        found = true;
        bestC = c;
        bestAdd = add;
        bestNewGc = newGc;
      } else {
        bool better = false;
        if (add > bestAdd) {
          better = true;
        } else if (add == bestAdd) {
          if (gcSteered) {
            const int mid = static_cast<int>((inst.lo + inst.hi) / 2);
            const int dNew = newGc > mid ? newGc - mid : mid - newGc;
            const int dOld = bestNewGc > mid ? bestNewGc - mid : mid - bestNewGc;
            if (dNew != dOld) better = (dNew < dOld);
            else better = (c < bestC);
          } else {
            better = (c < bestC);
          }
        }
        if (better) {
          bestC = c;
          bestAdd = add;
          bestNewGc = newGc;
        }
      }
    }
    if (!found) break;  // greedy cannot continue: infeasible prefix
    {
      int qq = q;
      char buf[4];
      daa::genetic::codonString(bestC, buf);
      for (int b = 0; b < 3; ++b) {
        const std::pair<int, bool> pr = ac.step(qq, daa::automaton::baseIndex(buf[b]));
        qq = pr.first;
      }
      q = qq;
    }
    last = bestC;
    gc = bestNewGc;
    codons.push_back(bestC);
  }
  // final validity check
  if (static_cast<int>(codons.size()) == inst.n) {
    if (sequenceValid(inst, ms, codons)) {
      res.feasible = true;
      res.score = daa::scores::sequenceScore(codons, st);
      res.codons = codons;
    }
  }
  return res;
}

SampleResult sampleAndReject(const daa::dp::DpInstance& inst,
                             const daa::scores::ScoreTable& st,
                             const daa::motifs::MotifSet& ms,
                             std::uint64_t seed,
                             long long budget,
                             double temperature) {
  SampleResult res;
  res.feasible = false;
  res.score = 0;
  res.trials = 0;
  res.wasted = 0;
  daa::prng::SplitMix64 rng(seed);
  std::vector<int> bestCodons;
  int64_t bestScore = 0;
  bool haveBest = false;
  for (long long t = 0; t < budget; ++t) {
    std::vector<int> cs(inst.n);
    int last = -1;
    int gc = 0;
    for (int i = 0; i < inst.n; ++i) {
      const std::vector<int>& opts = inst.allowed[i];
      double total = 0.0;
      std::vector<double> w(opts.size());
      for (std::size_t j = 0; j < opts.size(); ++j) {
        const int c = opts[j];
        const int64_t s = st.codonScore[c] + (last >= 0 ? st.pairScore[last][c] : static_cast<int64_t>(0));
        w[j] = std::exp(static_cast<double>(s) / temperature);
        total += w[j];
      }
      double r = rng.nextDouble() * total;
      double acc = 0.0;
      int chosen = opts.front();
      for (std::size_t j = 0; j < w.size(); ++j) {
        acc += w[j];
        if (r <= acc) {
          chosen = opts[j];
          break;
        }
      }
      cs[i] = chosen;
      last = chosen;
      gc += daa::genetic::gcCount(chosen);
    }
    std::string dna;
    int gc2 = 0;
    for (int c : cs) {
      char buf[4];
      daa::genetic::codonString(c, buf);
      dna += buf;
      gc2 += daa::genetic::gcCount(c);
    }
    ++res.trials;
    bool ok = true;
    if (daa::motifs::naiveMotifSearch(ms, dna)) ok = false;
    if (ok && (gc2 < inst.lo || gc2 > inst.hi)) ok = false;
    if (ok && !windowOk(inst, cs)) ok = false;
    if (ok) {
      const int64_t sc = daa::scores::sequenceScore(cs, st);
      if (!haveBest || sc > bestScore) {
        haveBest = true;
        bestScore = sc;
        bestCodons = cs;
      }
    } else {
      ++res.wasted;
    }
  }
  if (haveBest) {
    res.feasible = true;
    res.score = bestScore;
    res.codons = bestCodons;
  }
  return res;
}

}  // namespace baselines
}  // namespace daa
