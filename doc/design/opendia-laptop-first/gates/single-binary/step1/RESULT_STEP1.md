# PLAN_SINGLE_BINARY step 1 result: the C++ scorer

**Passed.** OpenDIA::CrossFitScorer (OpenMS b851c0504a) replaces train_gap.py + make_w.py + calib_tail.py.

**Inputs.** The PLAN_GAP arm-E scoring arrays (heuristics 11, feat12 6, featD 3) and the C++ interference features
(featC 4), i.e. the D11FDC inputs. The training rows come from train_windows.tsv: the anchor, anchor_decoy and decoy
groups, which are the ones train_gap uses. The Python references are the codon D11FDC outputs at seed offsets 0-4.

**Evaluation.** `eval_scores.py`, which uses the code of evaluate_v2.py with the arm-E library. It reproduces the Python reference
lines of `gap/arms/score/cppcheck` exactly (e.g. raw01 D11FDCWT 1464 / 1837 / 1928, q <= 0.01 1714, FDP_c 0.0070).

| | C++ seeds 0-4, N@FDP_c 1 % | median | Python D11FDC seeds 0-4 | median | ratio |
|---|---|---|---|---|---|
| raw01 | 1805, 1854, 1864, 1822, 1804 | 1822 | 1837, 1857, 1844, 1855, 1848 | 1848 | 0.986 |
| raw02 | 1826, 2021, 1905, 1866, 1930 | 1905 | 1912, 1956, 1903, 1925, 1976 | 1925 | 0.990 |

- **FDP_c at q <= 0.01:** C++ gives 0.0058-0.0082 (raw01) and 0.0025-0.0062 (raw02). Python gives 0.0058-0.0071 and 0.0012-0.0071.
- **Variability:** the C++ seed spread is wider on raw02 (1826-2021 against 1903-1976).
  - At 2 % C++ gives raw01 1950-1985 (Python 1928-2005) and raw02 2031-2116 (Python 2094-2141).
  - The C++ scorer is therefore at least as good on raw01 at 2 % and slightly weaker on raw02 at 2 %.
  - The plan did not make 2 % a criterion.
- **Acceptance** (median at least 98 % of Python, FDP_c <= 0.015, seed-deterministic, thread-independent):
  - The median and FDP_c criteria are met.
  - On raw01 seed 0, 1 and 8 threads give byte-identical output.
- **Speed:** raw01 takes 45-90 s on 8 threads with 0.57 GB peak, and 160 s on 1 thread. Python took 546 s on 6 threads.

**Not identical to Python by design.**
- The C++ scorer uses its own RNG (splitmix64, partial Fisher-Yates minibatches), double precision and hand-derived gradients. The unit test checks these gradients against finite differences, to a worst relative error below 1e-5.
- The default initialisation is uniform ±1/sqrt(fan_in) for weights and biases, as in torch's Linear bias. Torch's weight init is uniform ±1/sqrt(fan_in) too (kaiming_uniform with a = sqrt(5)).

Files: `EVAL.txt`, the per-seed logs `raw0?_s?.log`, `cpp_scorer_seeds.sh`, `eval_scores.py`.
