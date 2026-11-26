# DIAG_NONBORROW: near-human errors and a residual HeLa shift

26 November 2026. Retrospective diagnosis on used data, not a precheck or an FDR-control success. No binary, reference table, baseline label list or held-out run was changed.

## Findings

Near-human sequence similarity explains much of the Astral excess, but does not establish a labelling error. On HeLa it explains only part of the excess. Removing similar entrapments changes the estimand and is not a valid way to declare the original gate passed.

All analyses use accession `instr('_p_target')`, the existing step A exact-human relabel, tool q <= 0.01 and the original seeds. Hamming distance treats I and L as equivalent and searches same-length windows of the human FASTA. Distance 3 means >= 3, not an exact distance. Astral uses the human peptide FASTA extracted from the benchmark FASTA, whereas HeLa uses the reviewed human protein FASTA. Thus the two distance populations are not directly comparable.

## Sequence similarity and decoy control

| Set, G0 | Accepted E | Hamming-1 share E / decoys | Original combined | Decoy-matched H1 counterfactual |
|---|---:|---:|---:|---:|
| D | 1,032 | 13.0 / 2.7 % | 1.07 % | 0.96 % |
| H | 2,144 | 15.5 / 10.2 % | 1.28 % | 1.20 % |
| N2 | 2,555 | 13.7 / 10.8 % | 1.27 % | 1.23 % |
| Astral raw01 | 66 | 53.0 / 22.8 % | 1.23 % | 0.75 % |
| Astral raw02 | 48 | 56.2 / 25.3 % | 1.03 % | 0.60 % |
| Astral raw03 | 69 | 52.2 / 26.2 % | 1.21 % | 0.78 % |

Counterfactual: E_cf = E_non-H1 / (1 - accepted-decoy H1 share), with original N and r. This is an explanatory calculation, not an estimator with validated coverage or a candidate selection rule. It is not length-, mass- or score-matched. Astral accepted E are 82 % length 7-8, versus 15.5 % of library entrapments. Both decoys and entrapments are enriched near human sequences, so similarity alone does not justify relabelling.

In HeLa D, 37/134 accepted H1 entrapments have a human neighbour in the library and 32 have an accepted neighbour. Of those 32, 96.9 % are within 3 s. In H, the corresponding counts are 117/332 and 98, with 71.4 % within 3 s. This supports a mechanism for a subset, not for every H1 error. Modification-equivalent substitutions alone leave D 1.03 %, H 1.18 % and N2 1.19 %.

## Why the apparently passing exclusion analysis is not a fix

A library-only list of every entrapment at Hamming <= 2 was computed after inspecting the distance pattern. Excluding these E from accepted N and E and from library r gives G0 combined values D 0.83 %, H 0.90 %, N2 0.91 %, Astral raw01/02/03 0.41/0.16/0.17 %. These are sensitivity results only. The cutoff was not frozen prospectively and the remaining population is different.

The stricter check keeps every accepted identification and every entrapment as an error. It divides human targets and entrapments into near (distance <= 2 to a distinct human sequence for human targets) and far strata, using each stratum's own library ratio r_s. The diagnostic combined estimate is sum_s E_s(1+1/r_s)/N. Human exact self-matches are excluded when finding a distinct neighbour. This estimate still relies on within-stratum exchangeability and is not independently validated.

| G0 set | Stratified total | Near stratum | Far stratum |
|---|---:|---:|---:|
| D | 1.08 % | 1.35 % | 1.01 % |
| H | 1.29 % | 1.75 % | 1.15 % |
| N2 | 1.28 % | 1.56 % | 1.18 % |
| Astral raw01 | 1.32 % | 5.82 % | 0.46 % |
| Astral raw02 | 1.13 % | 5.95 % | 0.18 % |
| Astral raw03 | 1.33 % | 7.00 % | 0.19 % |

The far HeLa class remains above 1 %. Exclusion retained all human IDs, including near-human ones, in its denominator. The stratum calculation exposes why that could appear to pass. Neither calculation replaces the original benchmark.

## GBM promotion and features

For D, GBM-kept N=98,770 and E=522 give 0.60 % combined. Promoted N=11,140 and E=510 give 5.24 %, but the net implied false-target increase (+85) tracks the counted-decoy increase (+83). For H, kept N=175,006 and E=854 give 0.56 %. Promoted N=17,233 and E=1,290 give 8.57 %. The net implied false-target increase is +452 versus +144 counted decoys. These are descriptive class calculations with the global r, not independently calibrated class FDRs.

In H, 68.3 % of promoted E are distance >= 3, versus 70.1 % of kept E. Thus promotion is not exclusively a near-human mechanism. Score-binned feature contrasts on D/H seeds 0-1 show modest, mixed differences rather than one decisive separator: the largest reported accepted-region standardized mean difference is 0.21 decoy SD (D IM-centroid offset), and H's leading differences are boundary-apex agreement 0.144 and mt_z_med 0.137. These exploratory contrasts do not prove exchangeability or identify a lossless feature filter.

## Decision and evidence

- No FDR fix or benchmark correction is claimed. Original bounds and the ID-preservation requirement remain binding.
- The strongest remaining HeLa lead is the GBM-promoted far-from-human class, particularly H. A next candidate should target calibration or training rather than simply remove similar sequences.
- No new option, precheck or remote job was introduced here. The six PXD017703 runs remain unused for a prospectively frozen candidate.
- Scripts and textual outputs are in `near/`, with `SHA256SUMS`. Full generated lists and intermediate rows remain in `/home/singjc/.jcode/scratch/near/`. `strict_list.py` prints generic entrapment wording even in HUMAN/DECOYS/SEQS modes. Interpret those logs according to the mode, not that generic heading.
- `cf_hd1.txt` contains the Astral counterfactual rows only. HeLa rows are reproducible with `cf_hd1.py`, N=109910/192239/230993 and r=6.9176. The table above records the captured results.
- HeLa full-library list: 45,360 H1, 92,258 H2, 375,490 far peptides. Astral: 1,160 H1, 5,310 H2, 43,568 far. Evaluator r uses precursor counts, not these peptide counts. Existing sample-based `near_audit.py` estimates of r are superseded by `eval_strict.py` full-library results.
