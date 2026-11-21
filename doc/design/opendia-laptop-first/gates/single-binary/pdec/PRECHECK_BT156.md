# PRECHECK_BT156: does label-free hard-decoy training (Bt) hold at r ≈ 1.56 on the current default? (21 Nov 2026, diagnostic)

**Why.** RESULT_CEIL0 rule 4 made Bt a candidate: r0 +6.5 / +6.4 %, and r7 N@FDP_c 1 % +4.3 % (D) / +5.9 % (H) at seed 0. The same idea on the older default gained only +1.2 % at r ≈ 1.56 D, and its lower bound grew each round (PRECHECK_NOCAND section 3). This precheck repeats Bt at r ≈ 1.56 on the current default before any C++ port.

## Setup (fixed before any result)

- **Library:** `/homes/sing/hela_entrap_codon.predicted_decoys.pqp`, with r = 193,749 / 172,608 = 1.1225 as `eval_pdec.py` defines it.
- **Runs:** D and H, seeds 0 and 1, current default (`opendia-pd10-rtrdef.sif`), with `-features_out`. That is 4 runs.
- **Script:** `precheck_ceil0.py` (sha256 `19955add…`), unchanged, with r = 1.1225. Arms A (tool) and Bt are read from it. The other arms are reported but not used.
- **Measures:** N at combined FDP 1 % and the lower bound at own q ≤ 0.01, Bt against A, per run and seed.

## Decision rule (written before any result)

- **Pass**, so that a C++ port of Bt is written as an opt-in flag and then goes to a frozen r7 plan, if all of these hold:
  - the mean Bt gain over the 2 seeds is ≥ +3 % at N@FDP_c 1 % on D **and** on H;
  - every one of the 4 runs gains;
  - the mean lower bound rises by ≤ +0.2 pp on each run.
- **Fail** otherwise. Bt is then recorded as an r7-only effect and dropped.
