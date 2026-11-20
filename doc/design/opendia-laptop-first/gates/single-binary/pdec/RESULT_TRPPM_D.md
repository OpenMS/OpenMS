# RESULT_TRPPM_D: a ppm trace window gives +18 % IDs at 1 % true FDP on r ≈ 7 D. PLAN_TRPPM passes every criterion

- **Plan:** PLAN_TRPPM, frozen 0d80832844. `PLAN_TRPPM.sha256` reads OK, with sha 8c909603… and written_at 2026-11-20T15:09:37Z.
- **Runs:** codon `entlib-trppm-20261120`, jobs 38642251-67, submitted after the freeze. The binary sha256 was checked at submit time (`SHA256SUMS` OK).
- **Evaluation:** `eval_pdec.py` (sha256 bb03387d…). Output is in `trppm/EVAL_D.txt`, and wall time and RSS are in `trppm/fetch_D.txt`.

## Per seed

| Arm | Seed | N @ q ≤ 1 % | E | Lower | Combined | N @ FDP_c 1 % |
|---|---|---|---|---|---|---|
| CUR | 0 | 12,432 | 149 | 0.0120 | 0.0137 | 11,740 |
| CUR | 1 | 12,639 | 152 | 0.0120 | 0.0137 | 12,124 |
| CUR | 2 | 12,482 | 146 | 0.0117 | 0.0134 | 11,907 |
| CUR | 3 | 12,636 | 160 | 0.0127 | 0.0145 | 11,815 |
| CUR | 4 | 12,398 | 136 | 0.0110 | 0.0125 | 11,727 |
| NEW | 0 | 14,562 | 157 | 0.0108 | 0.0123 | 13,924 |
| NEW | 1 | 14,670 | 162 | 0.0110 | 0.0126 | 14,267 |
| NEW | 2 | 14,612 | 171 | 0.0117 | 0.0134 | 14,009 |
| NEW | 3 | 14,578 | 159 | 0.0109 | 0.0125 | 13,985 |
| NEW | 4 | 14,390 | 155 | 0.0108 | 0.0123 | 13,889 |

Every NEW seed is above every CUR seed. The CUR runs reproduce the current default: the seed 0 count of 11,740 matches the round-4 memory run.

## Criteria

| # | Criterion | Limit | NEW | Verdict |
|---|---|---|---|---|
| 1 | median N at FDP_c 1 % ≥ 1.05 × CUR (11,815) | ≥ 12,405.8 | **13,985** (1.184×) | pass |
| 2a | pooled lower bound ≤ CUR + 0.0030 | ≤ 0.0149 | **0.0110** (CUR 0.0119) | pass |
| 2b | pooled combined bound ≤ CUR + 0.0040 | ≤ 0.0176 | **0.0126** (CUR 0.0136) | pass |
| 3 | pair win ≤ 0.56 | ≤ 0.56 | **0.499** (n 7,287) | pass |
| 4a | median wall ≤ 1.3 × CUR (600.5 s) | ≤ 780.7 s | **569.2 s** (0.95×) | pass |
| 4b | peak RSS ≤ 64 GB | ≤ 64 GB | **34.2 GiB** max | pass |

**Verdict: PASS.**

## Reading

- **+18.4 % IDs at 1 % true FDP** (11,815 → 13,985). That is larger than at r ≈ 1.56 (+10-12 %), as expected when there are more entrapment competitors to pick up foreign ions.
- **Calibration got better here, not worse.** The pooled lower bound went 1.19 → 1.10 % and the combined bound 1.36 → 1.26 %. At r ≈ 1.56 the bound moved the other way (1.15 → 1.56 %), so that worsening did not carry over to the large library. The residual excess of about 1.1× (b9) is still there, in both arms.
- No extra cost: wall time 0.95×, the same RSS.
- **Against DIA-NN on D** (26,006 at FDP_c 1 %): 45 % → 54 %.
- **Next, per the plan:** H is submitted only now, after this verdict is committed, and is reported either way. Making 10 ppm the default is a separate commit, with ctests and a byte check against the tested arm.
