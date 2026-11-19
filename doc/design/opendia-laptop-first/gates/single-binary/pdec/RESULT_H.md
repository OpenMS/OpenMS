# RESULT_PDEC_H and RESULT_ODP_H: held-out HeLa H (reported only, no criteria)

- **Run:** H is HeLa 200 ng diaPASEF (21 min) on the r ≈ 7 foreign entrapment library, seeds 0-4, CUR and NEW.
- **Setup:** both sets ran with amendment 1: container `opendia-pd3-bc4e4f8d9b.sif`, `-trace_tensor_gb 96`, 192 GB.
- **Order:** both were submitted after the D verdicts were committed. PDEC H was submitted at 01:48Z, after 1f393c2184 and
  amendment 31b10fd49e. ODP H was submitted at 02:17:41Z, after 480991c28a.
- **Evaluation:** `eval_pdec.py` (bb03387d…), on the tool's own q. Outputs are in `EVAL_PDEC_H.txt` and `EVAL_ODP_H.txt`.
- **Cost:** about 1:00-1:07 wall and 110-117 GB peak per run. H is faster than D despite more windows (node variance, not
  analysed).
- **ODP library:** all 20 runs (D and H) wrote the same predicted library (sha256 76f2a6c2…).

## Results (H, pooled over 5 seeds, q ≤ 0.01)
| quantity (CUR / NEW) | DIALibGen-predicted (PDEC) | **OpenDIA `-decoys predict`** (ODP) |
|---|---|---|
| lower bound N_E/N | 0.0111 / 0.0109 | **0.0106 / 0.0110** |
| combined bound | 0.0127 / 0.0124 | **0.0121 / 0.0126** |
| 90 % interval, lower bound (CUR) | 0.0106-0.0117 | 0.0101-0.0111 |
| median N@FDP_c 1 % | 20,932 / 24,502 | 20,745 / 24,549 |
| N at q ≤ 0.01 (median) | 21,450 / 25,091 | 21,309 / 25,300 |
| P(entrapment beats own decoy), q ≤ 0.05 | 0.536 / 0.533 | **0.505 / 0.500** (±0.008) |

**For comparison, D (5 seeds):**
- PDEC: lower 0.0109 / 0.0111, combined 0.0124 / 0.0127, pair win 0.542 / 0.547.
- ODP: lower 0.0106 / 0.0106, combined 0.0121 / 0.0121, pair win 0.534 / 0.539.

## Reading
1. **H reproduces D for both routes.** The true FDP at q ≤ 0.01 is about 1.06-1.27 %. It is the same for both decoy
   sources, both arms and both runs. The roughly 85 % reduction from pseudo-reverse holds on the held-out run.
2. **`-pg_scores` holds on H.** It gives +17 % (PDEC) and +18 % (ODP) N@FDP_c 1 %, with the same calibration as CUR.
3. **New observation: with OpenDIA's decoys on H, pairs are fair (0.500-0.505), yet the excess is unchanged (1.06-1.10 %).**
   - So the residual excess **is not** a target-decoy pair asymmetry, at least on H.
   - On D, pair win was 0.534-0.539 with the same excess, so pair win and excess do not move together.
   - This points away from decoy construction and towards the **q-value estimation** at r ≈ 7. Candidates, none tested:
     - **The TDC estimator itself.** OpenDIA's paired competition counts decoy winners. The excess of about 1.1× could
       come from how ties go to the decoy, how the +1 correction is handled, or from unpaired members (0.15 %).
     - **Entrapment-specific competition.** An entrapment target competes with its own decoy, but a HeLa target competing
       with an entrapment precursor in the same window is a different case. The many entrapments at r ≈ 7 may crowd out
       HeLa targets' decoys or the reverse. That breaks the assumption that a false target is replaced only by its own
       decoy.
     - **Multiple candidates per precursor.** The precursor winner is the best window. A false target's best-of-N may
       differ in distribution from its decoy's best-of-N when N differs (more candidates in a dense region).
4. **The residual is small but consistent:** about 1.1× at q ≤ 0.01 (q ≤ 0.005 and q ≤ 0.02 on D gave 1.0-1.4×). Before
   predicted decoys it was 1.6-1.9×.

## Proposed next step (not started)
**PRECHECK_TDC.** On existing D and H outputs only, without new runs, recompute q three ways and compare with the
entrapment bounds:
1. The tool's own paired competition.
2. Plain unpaired TDC on all precursor winners: (D + 1) / T.
3. Paired TDC restricted to complete pairs.

If (2) or (3) is calibrated where (1) is not, the fix is in OpenDIA's q-value code, not in the library. That would be a
small, label-free C++ change for a frozen plan.

This also needs the window-winner counts per precursor, which `-features_out` gives on one rerun.
