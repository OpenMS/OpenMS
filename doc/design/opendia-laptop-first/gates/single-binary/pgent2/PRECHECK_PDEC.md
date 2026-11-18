# PRECHECK_PDEC: decoys with their own predicted RT, IM and intensities remove the TDC asymmetry (b3, diagnostic)

Written 18 November 2026. HeLa D, seed 0, CUR and NEW. Container `opendia-pg-a2aa8fc56c.sif`, codon.

## Library (`entlib-pdec-20261118`, `hela_pdec.pqp`, sha256 0920a1a3…)
- **Targets and entrapments:** exactly the PGENT2 `targets.pqp`. These are the 58,998 HeLa peptides plus 58,998
  homolog-free Arabidopsis peptides, as predicted by DIALibGen.
- **Decoys:** the pseudo-reversed sequence of every target peptide, **predicted by DIALibGen like a target**: its own RT,
  its own CCS/1/K0, and its own fragment intensities and top-6 ions. Pairing is by charge and oxidation count.
  - Built by `make_decoy_fasta.py` and `merge_predicted_decoys.py`.
- **Coverage:** 298,543 of 441,910 targets have a predicted decoy. The other 143,367 targets have none, because their
  reversed precursor or fragments fall outside DIALibGen's windows.
- **Pipeline bug fixed while building this:** OpenMS's PQP reader `INNER JOIN`s `PEPTIDE_GENE_MAPPING` whenever `GENE`
  has rows. Peptides without a gene mapping are silently dropped. The first build lost every decoy because of this.

## Results
| quantity | pseudo-reverse decoys (PGENT2, s0) | predicted decoys |
|---|---|---|
| **P(entrapment beats its own decoy), best q ≤ 0.05** (`PAIRWIN_pdec.txt`) | 0.579 / 0.581 | **0.503 / 0.502** |
| P(...), best q ≤ 0.01 | 0.630 / 0.644 | 0.470 / 0.555 (n 100-119) |
| P(...), all pairs | 0.500 | 0.499 / 0.500 |
| paired-only TDC at q ≤ 0.01: lower bound N_E/N (CUR / NEW) | 0.0088 / 0.0095 (r 1.56) | **0.0042 / 0.0050** (r 1.03) |
| paired-only TDC at q ≤ 0.01: combined upper bound | 0.0145 / 0.0155 | **0.0084 / 0.0099** |
| tool output as written: combined FDP_c | 0.0145 / 0.0155 | 0.0164 / 0.0168 |
| N at q ≤ 0.01 (tool) | 11,452 / 13,416 | 11,372 / 13,012 |

## Reading
1. **The high-score asymmetry is gone with predicted decoys.**
   - False targets beat their own decoy 50.2-50.3 % of the time among pairs good enough to be accepted (q ≤ 0.05).
     With pseudo-reverse decoys it was 58 %.
   - So the cause is decoy construction: a pseudo-reverse decoy carries its target's predicted RT, IM and intensities.
     Those are mismatched to its own sequence, so it gathers co-eluting fragment evidence less easily than a
     self-consistent false target.
2. **With paired-only competition, q ≤ 0.01 is calibrated: true FDP 0.4-1.0 %.**
   - Both bounds sit at or below 1 % for CUR and NEW.
   - This is a diagnostic re-competition, not the tool's own q. It uses only pairs that have both members.
3. **The tool's own output is not yet calibrated** (combined 1.6-1.7 %), because 32 % of targets have no predicted
   decoy:
   - They compete with nothing and pass on the target score alone.
   - 514 / 594 of the accepted IDs are such unpaired targets, and 67 of each are entrapments.
   - This is an artefact of how this test library was built (DIALibGen windows), not of the method.
   - A production version must generate a decoy for every target, for example by relaxing the windows for decoys or by
     predicting the decoy at the target's m/z window.
4. **`-pg_scores` remains calibrated like CUR,** and keeps its gain: +14 % N at q ≤ 0.01 with paired-only competition.

## Consequence (b)
The diaPASEF calibration defect is caused by how decoys are made for predicted libraries, and the fix is known:
**decoys must be predicted like targets.** That is the DIA-NN approach. Options:

- **(i) Library side:** generate decoys in DIALibGen with a predict-decoys mode, one per target, and use them with
  OpenDIA unchanged.
  - The fastest route. It needs DIALibGen to emit a decoy for every target, which it does for its own
    `-generation:decoys` (but those currently copy the target's predictions; to be checked).
- **(ii) OpenDIA side:** OpenDIA builds its decoys itself from the target library, using an RT/IM/MS2 predictor at run
  time.
  - This is self-contained for users with any library, but it adds a model dependency to OpenDIA.

Which is preferable is a design decision for the user. Either way, the next frozen plan would test it with:
- the r ≈ 7 foreign entrapment (tight bounds);
- the criterion that the lower bound at q ≤ 0.01 is ≤ 1.0 % and the combined bound is ≤ 1.2 %;
- D as design and H as held out.

**Independent corroboration:** DIALibGen's own help text for `-search:intensities` describes the target-assay decoy (the
decoy in its target's fragment slots) as giving "weaker decoys". The predicted option is its default for its built-in
search.
