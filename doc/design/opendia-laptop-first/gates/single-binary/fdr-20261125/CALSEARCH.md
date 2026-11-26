# CALSEARCH: label-free calibration candidates on dev data (26 Nov 2026, diagnostic)

**Status:** exploration on dev data only. Nothing here is a gate result. The candidates and cuts were looked at with entrapment labels, so any winner needs a frozen precheck on new held-out data.

**Data:**
- Astral raw01/02 x seeds 0-4, G0 (default) and A (`-evidence_share`) runs, `features.tsv`;
- HeLa D and H seeds 0-1, default, re-run with `-features_out` on codon (`calsearch/cal_hela.sbatch`, container `opendia-evshare-20261126.sif`).
  - All 4 HeLa runs are byte-identical to the G0 references on every output file, so the feature export does not change the result.
  - The HeLa tables were reduced on codon to the scored rows and 6 columns (`reduce.sh`, sha256 in `hela_reduced.sha256`).
- Labels: instr `_p_target`. HeLa uses the step A relabel list.
- Evaluation: precursor-paired TDC as the tool; combined bound E(1 + 1/r)/N at q <= 0.01; N @ FDP_c 1 %.

## Candidates (all label-free)

| Name | Rule |
|---|---|
| GATE_t | rows whose network-only q > t (targets and decoys alike) are ranked below all others |
| GRP_X@p | split the TDC winners into two groups by feature X above the p quantile of decoy winners with q <= 0.2; compute q per group |
| GRP_ANY@p | the same, with group = any of `d_lib`, `d_max`, `i_nclaim` above its cut |
| X@3g | three groups, cuts at the 0.75 and 0.95 quantiles |

`cal_explore.py` (Astral A) motivated the features. At q <= 0.01, 2E/D is highest in the top bins of `d_lib` (5.0), `d_max` (3.8) and `i_nclaim` (3.3), against 1.2 overall. In the null region, T/D is about 1 in every bin.

## Results

**Astral dev, 10 runs** (`rob_G0.txt`, `rob_A.txt`, `sim_*.txt`):

| Base | Rule | Combined | N@FDP_c 1 % (sum) | Worst single run |
|---|---|---:|---:|---:|
| G0 | none | 1.13 % | 1.000x | |
| G0 | GATE_0.02 | 0.99 % | 1.024x | |
| G0 | GRP_d_max@0.9 | 0.94 % | 1.027x | 0.971x |
| G0 | GRP_ANY@0.75 | 0.96 % | 1.038x | 0.991x |
| A | none | 1.12 % | 1.000x (1.035x of G0) | |
| A | GRP_d_lib@0.9 | 0.90 % | 1.013x | 0.981x |
| A | GRP_ANY@0.75 | 0.95 % | 1.047x | 0.995x |

- Results depend on the cut: neighbouring cuts give 0.90-1.11 % and single-run losses down to 0.84x.
- GATE is not monotone in t (0.01: 0.84 %, 0.05: 1.17 %).

**HeLa D and H dev, 2 runs each, default** (`rob_D.txt`, `rob_H.txt`):

| Set | Rule | Combined | N@FDP_c 1 % |
|---|---|---:|---:|
| D | none | 1.12 % | 1.000x |
| D | best (GRP_ANY@0.75) | 1.06 % | 1.009x |
| H | none | 1.24 % | 1.000x |
| H | best (GRP_ANY@0.9 / 0.95) | 1.21 % | 1.001x |

## Reading

- **Exposure-split TDC is not a general fix.** On Astral it brings the bound under 1 % with a small ID gain. On HeLa it removes only 0.03-0.06 pp, and the bound stays well above 1 %.
- **Why (likely, not proven):** on Astral the excess sits in a small, extreme tail that one feature cut can isolate. On HeLa, step A found excess also in the class with no shared-fragment partner (H: 1,822 E against 1,604 predicted), which no exposure split can reach. The 4+ class (H: 168 against 98) is not what `d_lib`, `d_max` or `i_nclaim` above the decoy-90 % cut select.
- No candidate meets the goal (bound <= 1 % at no ID cost) on both instruments. Nothing is proposed for a precheck from this round.
- A label-free feature that measures 4+ fragments shared with a higher-ranked co-isolated ID within 25 Th is the remaining lead on HeLa. The current `i_*` and `e_*` features do not measure exactly that.
