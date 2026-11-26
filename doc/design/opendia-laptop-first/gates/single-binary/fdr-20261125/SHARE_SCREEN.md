# SHARE_SCREEN: a label-free shared-fragment competition step cannot fix the diaPASEF excess (26 Nov 2026, dev screen)

**Question (user request):** test the remaining lead from RESULT_EVSHARE and CALSEARCH. That is a label-free step aimed at targets that share 4 or more library fragments with a higher-ranked, co-eluting, co-isolated ID on diaPASEF.

**Status:** offline screen on already-used data only. The selection rule for a held-out precheck was fixed in `share/share_screen.py` before running. No rule passed it, so **no precheck was frozen and no held-out run was made**.

## Set-up

- **Data:** existing default (G0) run tables, 5 seeds each.
  - HeLa D and H: `entlib-gbt-20261122`;
  - HeLa N2: `evshare-20261126`, used in RESULT_EVSHARE;
  - Astral raw01-03: `astral-gbmfix-20261123`, 15 runs.
- **Labels:** instr `_p_target`, with the step A relabel lists. They are used only for evaluation and for the ORACLE ceilings.
- **share(p):** the most library fragments (15 ppm) that p shares with one donor.
  - A donor is an accepted target (baseline q <= 0.01, entrapments included, because a real rule cannot see labels).
  - It must have a higher score, a different modified sequence, an apex within 3 s, and a precursor m/z within the isolation half-width (12.5 Th diaPASEF, 1 Th Astral).
  - In the `_IM` variants it must also have a library 1/K0 within 0.05.
  - share is computed for every row (target or decoy) with baseline q <= 0.05.
- **RULE_K:** remove every row with share >= K, then recompute paired TDC. This is label-free.
- **ORACLE_K:** remove only the entrapment targets with share >= K. This is the ceiling for any lossless fix of that class.
- **ORACLE_ALLE:** remove every entrapment in the region. This is a sanity row.
- **Selection rule (fixed before running):** a RULE goes forward only if, on every set, the pooled combined bound is <= 1.00 % and the median N@FDP_c 1 % is >= 1.00x the baseline.

## Results (pooled combined bound at q <= 0.01 / median N@FDP_c 1 % vs baseline)

| Set | Baseline | RULE_2 | RULE_3 | RULE_4 | ORACLE_2 (ceiling) | ORACLE_4 (ceiling) |
|---|---|---|---|---|---|---|
| HeLa D | 1.07 % | 1.01 % / 0.999x | 1.02 % / 1.002x | 1.03 % / 0.999x | 0.85 % / 1.028x | 0.99 % / 1.014x |
| HeLa H | 1.28 % | 1.21 % / 0.995x | 1.22 % / 0.999x | 1.23 % / 0.999x | **1.09 %** / 1.009x | **1.17 %** / 1.006x |
| HeLa N2 | 1.27 % | 1.26 % / 0.997x | 1.27 % / 0.998x | 1.27 % / 0.998x | **1.20 %** / 1.001x | **1.24 %** / 1.000x |
| Astral (15) | 1.10 % | 1.07 % / 1.000x | 1.07 % / 1.000x | 1.07 % / 1.000x | 1.07 % / 1.000x | 1.07 % / 1.000x |

- K = 5 and 6 are weaker still. The IM gate changes almost nothing: it moves only 0-2 accepted human targets per run out of the share >= 4 class. Donors that share 4+ fragments are almost always IM neighbours already.
- **No RULE passes the selection rule.** None reaches <= 1 % on H, N2 or Astral, and most lose IDs (0.994-0.999x).
- **ORACLE_ALLE** gives 0 % and IDs of 1.06-1.30x. It is only a sanity row.

## Reading

- **The ceiling decides it.** Even removing every borrowing entrapment, using the labels, still leaves the bound above 1 %:
  - H: 1.09 % with share >= 2 and 1.17 % with share >= 4;
  - N2: 1.20 % and 1.24 %;
  - Astral: 1.07 %.
  
  So no label-free rule aimed at this class can bring H, N2 or Astral under 1 %. A held-out precheck would fail by construction, and one was not run.
- **The 4+ class is small.** Per run it holds 9-35 accepted entrapments (D 14-19, H 33-35, N2 9-12), against 190-530 accepted entrapments in all. The DIAG_STEPA excess counts in that class were real but are a small part of the total.
- **The label-free rule removes more than the oracle.** It also removes real co-isolated human targets: about 120 per run on N2 and 260-290 on D/H have share >= 4. That is why it costs IDs on every HeLa set.
- **Where the remaining excess is:** most of it sits in entrapments that share at most 1 fragment with any accepted co-eluting ID (DIAG_STEPA: on H, 1,822 E against 1,604 predicted in that class). It is not explained by co-isolation borrowing as measured here. Candidate causes still open:
  - the residual P of about 0.52-0.53 for entrapment against own decoy in the "neither borrows" class (DIAG_STEPA);
  - entrapment-specific library properties;
  - a GBM-wide calibration shift. On H, network-only is 0.98 % lower bound and the default is 1.12 % (DIAG_STEPA).

## Held-out set prepared, not used

- PXD017703 `HeLa_Evosep_diaPASEF_RAW.zip` (PRIDE SHA-1 `16410713…`, OK) was downloaded to codon `share-20261126/data`.
- It holds six never-used HeLa diaPASEF runs: 60SPD 2450/2451/2452 (200 ng) and 100SPD 2731/2732/2735 (100 ng), from the same lab and LC as D.
- No search has been run on them. They remain a clean held-out set for a future candidate.
