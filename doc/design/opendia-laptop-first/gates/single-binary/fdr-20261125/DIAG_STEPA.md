# DIAG_STEPA: HeLa's excess is half a labelling artifact, and the rest is co-isolation borrowing (25 Nov 2026, diagnostic)

- **Scope:** step A of the post-RESULT_BLEND plan. This is analysis only. No tool change and no new run.
- **Inputs:** existing run tables.
  - Astral raws 01-03 × 5 seeds: G0 = default and B5 = `-gbm_blend 0.5`.
  - HeLa D and H at r ≈ 7 × 5 seeds:
    - NET = network only, `entlib-gbt-20261122/runs/*_cur_*`;
    - G0 = default with GBM;
    - B5 = `-gbm_blend 0.5`.
- **Records:** scripts and outputs are in `stepA/`. Inputs are hashed in `stepA/inputs.sha256`.
- **Caveat:** every split and rule here is post hoc. Nothing below is a validated result.

## Result 1: the HeLa entrapment set contains 2,648 human peptides

- **How it happens:** the HeLa library build digests each foreign entrapment FASTA entry again, with up to 2 missed cleavages. A foreign entry that passed the homolog-free rule as a whole can still contain a sub-peptide that is an exact human sequence. Example: `LGANSLLDIVVFGRACANR` (A. thaliana) yields `LGANSLLDIVVFGR` (human SDHA). That precursor maps only to `_p_target` accessions, so every evaluator counts it as an entrapment.
- **Audit:** `label_audit.py`. Exact substring of the build-time human FASTA (`uniprotkb_Human…2026_02_26`, the same sha256 as on Codon), with I = L.
  - **HeLa:** 2,648 of 515,756 entrapment peptides (0.51 %, 2,912 precursors). Mostly length 7-8.
  - **Astral:** 278 of 50,316 against the same human FASTA, 472 precursors. Only 19 are exact substrings of the library's own human entries.
- **They are the hits that recur.** At own q ≤ 0.01, entrapment peptides accepted in at least 3 of 5 seeds overlap between the two independent acquisitions D and H as follows (`recur.py`, `recur_decoy.py`):

  | | recurrent D | recurrent H | overlap D × H | expected by chance |
  |---|---:|---:|---:|---:|
  | accepted entrapments | 193 | 422 | **39** | 0.16 |
  | counted decoys (score ≥ the lowest accepted target) | 182 | 365 | 9 | |

  - 23 of the 39 are exact human substrings (`fasta_near.py`): histones, IF4A1, RAB5, RPS4X and similar.
  - The control set of matched-length entrapment peptides that do not recur sits at median Hamming 3 from any human window.
- **A side note on the frozen evaluator:** its SQL `LIKE '%_p_target%'` is case-insensitive and treats `_` as a wildcard. So it also labels as entrapment the human target peptides that end in P, whose accession reads `…P_target`: 19 proteins on HeLa and 101 on Astral. Using the exact substring leaves Astral E unchanged. HeLa E at q ≤ 0.01 falls by 20-25 per arm (D default 1,222 to 1,202, H default 2,476 to 2,451), and median N @ FDP_c 1 % moves by up to 0.5 % (`eval_blend_instr.txt`). Any new evaluator should use `instr`. The "as labelled" column below already uses the exact substring.

## Result 2: with corrected labels, only the GBM default on H still shows a provable excess

Same evaluator logic, with the 2,912 human-substring precursors counted as human (`eval_pdec_relabel.py`, `eval_relabel_*.txt`). Own q ≤ 0.01, pooled over 5 seeds:

| Set | Arm | lower (as labelled) | lower (relabelled) | combined (relabelled) | median N @ FDP_c 1 % (relabelled) |
|---|---|---:|---:|---:|---:|
| HeLa D | network only | 1.08 % | **0.94 %** | 1.08 % | 20,361 |
| HeLa D | default (GBM) | 1.09 % | **0.94 %** | 1.07 % | 21,805 |
| HeLa D | blend 0.5 | 1.01 % | 0.87 % | 0.99 % | 20,985 |
| HeLa H | network only | 1.14 % | **0.98 %** | 1.13 % | 35,333 |
| HeLa H | default (GBM) | 1.27 % | **1.12 %** | 1.28 % | 37,897 |
| HeLa H | blend 0.5 | 1.15 % | 1.00 % | 1.14 % | 36,021 |
| Astral, 15 runs | default (GBM) | 0.59 % | 0.58 % | 1.16 % | 2,208 |
| Astral, 15 runs | blend 0.5 | 0.40 % | 0.40 % | 0.81 % | 2,231 |

**Reading:**
- On HeLa, about half of the reported excess was the label artifact.
  - **D:** the remaining lower bound is below 1 % for every arm.
  - **H:** only the default with GBM still has a lower bound above 1 % (1.12 %, 90 % interval 1.08-1.16 %).
- The ID ranking is unchanged:
  - on D, the default stays +7 % over network only and +4 % over the blend;
  - on H, it stays +7 % and +5 %.
- **RESULT_BLEND's verdict is unchanged.** With corrected labels, the blend's HeLa IDs are still 0.962× (D) and 0.950× (H) of the default.
- Earlier HeLa records also used the old labels: BOUNDS_R5, RESULT_GBT_D/H, RESULT_DIANN_R5 and RESULT_BLEND. Their FDR statements are therefore somewhat pessimistic for every tool. The DIA-NN comparison used the same labels, so its relative ranking is unaffected, but its absolute bounds should be re-derived the same way.
- **The Astral excess is not explained by labels.** Relabelling there is negligible (E 187 to 183).

## Result 3: the remaining excess on both instruments is co-isolation borrowing

**Test** (`pair_borrow.py`, relabelled, `pair_borrow_hela_relabel.txt`):
- For each member of an entrapment / own-decoy pair, compute "share". This is the most library fragments it shares with one accepted human ID whose apex lies within 3 s and whose precursor m/z lies within the isolation window: ±12.5 Th on diaPASEF, ±1 Th on Astral.
- **P** is the fraction of these pairs in which the entrapment scores higher than its own decoy. It should be 0.5 if the decoy is a fair stand-in for the entrapment.

| HeLa H, default (pairs with best q ≤ 0.05) | n | P(entrapment beats own decoy) |
|---|---:|---:|
| neither member shares ≥ 2 fragments with a co-eluting ID | 17,348 | 0.533 |
| entrapment shares 2-3 | 953 | 0.717 |
| entrapment shares 4+ | 309 | 0.951 |
| decoy shares 2-3 | 933 | 0.319 |
| decoy shares 4+ | 216 | 0.042 |

- **Borrowing is symmetric:** decoys borrow too, at about the same rate as entrapments. Of the counted decoys on H, 12 % share ≥ 2 fragments with a co-eluting ID, against 15 % of accepted entrapments (`borrow_hela.txt`).
- **The asymmetry sits in the top class (4+ shared).** At q ≤ 0.01, the 4+ class holds:
  - D: 68 E against the 41 that decoys predict;
  - H: 168 against 98.
  
  The 0-1 class is about at parity on D (835 against 857), and above it on H (1,822 against 1,604).
- **The GBM amplifies it.** Network only on H has 113 E in the 4+ class against 55 predicted. With the GBM it is 168 against 98.
- **Astral has the same shape at a smaller scale** (`pair_borrow_astral.txt`): 5 E against 1 predicted in the 4+ class. Most Astral near-isomers have no co-eluting confident partner (`isomer_partner_astral.txt`). Only 18 % of accepted near-isomer entrapments have one within 3 s, which is why the Astral shadow rule removed so little.
- **Why it is asymmetric.** A pseudo-reverse decoy keeps the precursor mass but reverses the fragment ladder, so its fragments line up with a co-isolated real peptide less often. When the evidence comes from 4 or more shared fragments, the false target has the stronger borrowed match. This is a property of target-decoy competition with reversed decoys, and the GBM learns to reward it.
- **Also found:** after relabelling, a residual P ≈ 0.51-0.53 remains in the "neither borrows" class on HeLa. Its source is not identified here. Library-level differences (RT, IM, intensity shape, ion types) move P by at most ±0.04 (`pair_lib_hela.txt`).

## Result 4: a label-free competition rule, simulated offline, removes some of the excess but costs IDs on D

**Rule** (`sim_compete.py`): after scoring, remove any target or decoy that shares ≥ K library fragments with a higher-scored, co-eluting (3 s), co-isolated target of a different peptide. Then recompute paired TDC. The same rule applies to every row, and no labels are used.

| Set (old labels) | K | lower before → after | combined before → after | median N @ FDP_c 1 % before → after |
|---|---:|---|---|---|
| Astral, 15 runs | 2 | 0.57 → 0.55 % | 1.13 → 1.09 % | 2,160 → 2,159 |
| HeLa D | 3 | 1.09 → 1.04 % | 1.25 → 1.19 % | 21,495 → 21,374 |
| HeLa D | 4 | 1.09 → 1.04 % | 1.25 → 1.18 % | 21,495 → 21,357 |
| HeLa H | 3 | 1.27 → 1.20 % | 1.46 → 1.37 % | 37,376 → 37,500 |
| HeLa H | 4 | 1.27 → 1.21 % | 1.46 → 1.38 % | 37,376 → 37,458 |

- **The rule is too blunt.** On HeLa it removes 1,260-1,620 accepted human targets per 5 runs, many of them real co-eluting peptides that share fragments. It removes 104-292 accepted entrapments.
- **Gains:** a +0.3 % ID gain on H, and less than 0.1 pp change on Astral.
- **Losses:** −0.6 % IDs on D.
- **Verdict:** this is not a fix as written. A better fix would share evidence rather than remove whole rows. One example: score each precursor's fragments after subtracting the part explained by a higher-ranked co-isolated ID.

## What this changes

1. **HeLa's FDR problem was mostly measurement.**
   - With correct labels, no network-only arm and not the default on D shows a provable excess: the lower bounds are 0.94-0.98 %.
   - The combined bounds stay at 1.07-1.13 %, so calibration is not proven either.
   - The remaining provable issue is H with the GBM: lower bound 1.12 %, about 0.1 pp over 1 %.
2. **The real mechanism is one shared cause, borrowing from co-isolated IDs, on both instruments.** The GBM amplifies it.
   - Astral: through near-isobars, few of which have a confident co-eluting partner.
   - HeLa: through 4+ shared fragments within the 25 Th diaPASEF windows.
3. **No new decoy type is needed** for what remains. The candidate fixes, in order:
   - **(a) Evidence sharing.** Fragment-level deconvolution against higher-ranked co-isolated IDs, or a feature computed the same way for targets and decoys. It must not be a label-dependent rule.
   - **(b) Stop the GBM rewarding borrowed evidence.** For example, exclude or decorrelate the shared-fragment features from the GBM's inputs.
   - Either needs a frozen precheck. HeLa D and H and Astral raw 03 have been used, so a new held-out set is required.
4. **Evaluation must use corrected labels from now on.** A label audit of the library (`label_audit.py`) should become part of building any entrapment library.
