# RESULT_BTCPP: hard-decoy training does not transfer to the C++ network scorer (21 Nov 2026)

The gate is in `PRECHECK_BTCPP.md` (`64e28decdc`). OpenMS `f3005a4b14`, container `opendia-pd10-btopt.sif` (sha256 `a3a75146…`). The 12 runs are in `$C/btcpp-20261121/runs` and all exited 0. The evaluator is the frozen `eval_pdec.py` (sha256 `bb03387d…`), with r = 1.1225.

**Byte check:** cur `precursors.tsv` is byte-identical to `bt156-20261121/{D,H}_s{0,1}` (container rtrdef) for all 4 runs, so the option is inert at its default.

| Arm | D s0 | D s1 | D mean Δ | D lower (pooled) | H s0 | H s1 | H mean Δ | H lower (pooled) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| cur | 16,678 | 17,038 | | 0.87 % | 32,310 | 32,415 | | 0.67 % |
| bt1: `-hard_decoys 0.02` | 15,806 | 15,984 | **−5.7 %** | 0.78 % | 31,394 | 31,493 | −2.8 % | 0.66 % |
| bt2: + `-scorer_self_q "0.01 0.01 0.01"` | 15,844 | 15,825 | **−6.1 %** | 0.77 % | 32,037 | 31,407 | −1.9 % | 0.68 % |

The values are N at combined FDP 1 %.

**Verdict: FAIL.** Both arms lose on every run. The flag stays opt-in and unused.

**Why it fails here but passed in Python (RESULT_BT156, GBM +3.7 % / +6.0 %).** The C++ learner scores with a decoy-whitened readout, w = (S_N + λI)⁻¹ (μ_P − μ_N). The negatives set both the centre μ_N and the covariance S_N that defines the score metric. By default the negatives are the winners of all decoy precursors (about 150k per fold), which describes the whole decoy population. With hard decoys they shrink to about 19-23k windows from about 4,500 precursors (fold logs: `hard_neg 19190 of 4515 prec`). Those are the decoys that already look most like targets, so the whitening is fitted to a narrow, target-like cloud. Most decoys then sit outside it and score high. Fewer targets are accepted at the same FDR, and the lower bound falls (D 0.87 → 0.78 %) because the extra decoys crowd out entrapments as well. The GBM has no such coupling: its negatives only shape the decision boundary.

**What this means for the scorer lever.** The +4-6 % seen with the GBM requires a learner whose fit does not depend on a negative covariance. The options for a future precheck, none of them gated yet:

- (a) Keep the default whitening (all decoy winners) but up-weight hard decoys in the network's loss.
- (b) Add a small GBM or logistic final round on the network embedding.
- (c) Replace the readout for the last round only.

Each one is a scorer redesign, not a flag.
