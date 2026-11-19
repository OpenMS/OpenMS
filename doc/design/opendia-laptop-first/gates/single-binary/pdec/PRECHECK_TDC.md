# PRECHECK_TDC: where the residual ~1.1× sits at r ≈ 7 (diagnostic, ODP D and H outputs)

Written 19 Nov 05:05. This uses only the existing PLAN_ODP outputs: D and H, CUR and NEW, 5 seeds. Labels are used for
evaluation only, and no verdict changes. Scripts: `tdcdiag.py` (output `TDCDIAG.txt`) and `sharediag.py`
(`SHAREDIAG.txt`).

## 1. How the tool's q is formed
`pooledPairedTDC` (OpenDIAScorer.cpp) works as follows:
- It takes the best window per precursor.
- Within each target-decoy pair it takes the better member; ties go to the decoy.
- It ranks the pair winners and sets q = (D + 1) / T.

So the tool reports D/T ≈ q, and it does so exactly (D/T = 0.0098-0.0099 at q ≤ 0.01).

## 2. Accepted decoys by pair type (q ≤ 0.01, pooled over 5 seeds)
| run, arm | T | D | entrapments E | decoys of entrapment pairs D_E | decoys of HeLa pairs D_H | E / D_E |
|---|---|---|---|---|---|---|
| D CUR | 43,243 | 425 | 458 | 332 | 93 | 1.38 |
| D NEW | 52,096 | 514 | 553 | 408 | 106 | 1.36 |
| H CUR | 106,684 | 1,060 | 1,128 | 858 | 202 | 1.31 |
| H NEW | 126,332 | 1,257 | 1,391 | 1,035 | 222 | 1.34 |

**The ratio E / D_E by cut:**
- q ≤ 0.005: **1.68-1.87**
- q ≤ 0.01: 1.31-1.38
- q ≤ 0.02: 1.16-1.23

It falls towards 1 as the cut loosens.

**Reading.** Paired TDC assumes that each accepted false target is matched by an accepted decoy from the same kind of
pair. For entrapment pairs, the accepted entrapment targets outnumber their own decoys by 1.3-1.9×, most strongly at the
strictest cut.
- This is not visible in the pair win above, because P(win | best q ≤ 0.05) is dominated by pairs near the cut.
- At the top of the ranking, an entrapment target beats its own decoy more often than a decoy does. A few entrapment
  targets score like true targets.
- About 19 % of the gap is offset by HeLa-pair decoys (D_H). This is why the tool's total D/T is close to the true FDP.

## 3. Why the top entrapments score high: shared signal
For each accepted winner (D and H seed 0, CUR), the diagnostic finds the largest number of its fragments within 20 ppm of
a co-eluting identified HeLa target (|ΔRT| ≤ 10 s, |Δm/z| ≤ 12.5).

| class | D: ≥ 2 shared / ≥ 3 shared | H: ≥ 2 / ≥ 3 |
|---|---|---|
| accepted entrapment targets | 0.191 / **0.138** | 0.170 / **0.140** |
| accepted decoys, entrapment pairs | 0.190 / 0.095 | 0.122 / **0.052** |
| accepted decoys, HeLa pairs | 0.200 / 0.150 | 0.150 / 0.050 |
| accepted HeLa targets (sample) | 0.041 / 0.023 | 0.041 / 0.029 |

- Accepted entrapment targets share ≥ 3 fragments with a co-eluting identified HeLa precursor about 1.5-2.7× as often as
  accepted decoys of entrapment pairs (H: 14.0 % against 5.2 %).
- This is the single-interferer, shared-ladder effect PAIRWIN found earlier (about 13 % of the excess then). It survives
  predicted decoys.
- A foreign peptide's fragments can coincide with a real co-eluting HeLa peptide's fragments, because both are real
  b/y-ion masses near the same m/z. Its reversed decoy's fragments coincide less often, because they are shifted.
- **Caveat:** these are small counts (n 63-229 per class), from one seed.

## 4. Implication
- **The residual is not a q-estimator bug.** Plain TDC equals D/T here, so swapping the estimator changes nothing. It is
  a scoring problem: some false targets explain real signal that belongs to another, identified precursor, and decoys
  cannot mimic this.
- **Candidate fix, label-free and run-derived:** stronger shared-signal attribution (interference resolution).
  - If a window's supporting fragments are already explained by a co-eluting, higher-scoring precursor, discount them.
  - OpenDIA has an InterferenceResolver feature set. The question is whether its features reach the top-ranked winners
    or whether the scorer underweights them.
  - Next precheck: for accepted winners, does the existing interference feature separate "entrapment sharing ≥ 3
    fragments" from "decoy"? That needs one `-features_out` run (D s0, about 2 h).
- **Benchmark side:** the same effect produces genuine false IDs in real searches, for example of near-isobaric human
  peptides. So this is not only an entrapment artefact. Fixing it is real FDR work, not tuning to the benchmark.
