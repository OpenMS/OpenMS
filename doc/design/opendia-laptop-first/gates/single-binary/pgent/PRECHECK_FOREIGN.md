# PRECHECK_FOREIGN: is OpenDIA's diaPASEF q too optimistic, or are the entrapments too close to human? (b1, diagnostic)

> **Correction (PLAN_PGENT2 verdict, `pgent2/VERDICT_DESIGN.md`):** the post-hoc homolog filter below flattered both
> arms. With near-homologs removed at library construction, FDP_c at q ≤ 0.01 is 1.43 % (CUR) and 1.59 % (NEW) over 5
> seeds. The conclusion "q is close to calibrated" is withdrawn: the diaPASEF excess is real, at about 1.5 times
> nominal.

Written 18 November 2026. HeLa D, seed 0, CUR and NEW (`-pg_scores`), on codon. Labels are used here for evaluation only.

## Library (`entlib-ath-20261118`, sha256 b93f2872…)
- **Targets:** the 58,998 HeLa library peptides, the same as entlib v1.
- **Entrapments:** 58,914 *Arabidopsis thaliana* tryptic peptides (UniProt UP000006548, reviewed).
  - Drawn with the same length distribution as the targets.
  - Every peptide that is a human tryptic peptide or a substring of a human protein (I=L) was removed.
- **Prediction:** DIALibGen with the same timsTOF settings as entlib v1. Decoys are pseudo-reverse. r = 1.557.
- **Inputs and scripts:** `inputs_ath.sha256`, `make_fasta_ath.py`.

## Results (`EVAL_D_ath.txt`, `RHO_ath.txt`, `homolog.py`, `strict.py`)
| entrapment set | arm | N at q ≤ 0.01 | FDP_c at q ≤ 0.01 (90 % Poisson) | N@FDP_c 1 % |
|---|---|---|---|---|
| all Arabidopsis | CUR | 11,067 | 0.0172 | 9,114 |
| all Arabidopsis | NEW | 13,208 | 0.0148 | 12,225 |
| Hamming > 2 to human (82 % of entrapments) | CUR | 11,026 | **0.0116** (0.0095-0.0141) | 10,800 |
| Hamming > 2 to human (82 % of entrapments) | NEW | 13,164 | **0.0097** (0.0080-0.0118) | 13,197 |

**Accepted Arabidopsis entrapments are mostly near-homologs of human peptides:**

| group | within Hamming ≤ 1 of a human peptide | Hamming ≤ 2 | share a ≥ 6-residue y-suffix |
|---|---|---|---|
| accepted at q ≤ 0.01 (CUR / NEW) | 22 % / 26 % | 32 % / 34 % | 28 % / 35 % |
| random entrapment winners (null) | 3 % | 9-11 % | 9 % |

## Reading
1. **Most of the "1.7 % true FDR" is a benchmark artefact.**
   - Accepted entrapments are 3-4 times enriched for peptides one or two residues from a real human peptide. These are
     peptides that share most fragments with something present in the sample, so they are partly correct
     identifications of a homologous sequence, not random false targets.
   - The same holds for shuffled entrapments that keep cleavage sites and composition.
   - Excluding near-homologs from both the count and r, FDP_c at q ≤ 0.01 is 1.16 % for CUR and 0.97 % for NEW.
     Both are within noise of nominal.
2. **OpenDIA's q calibration on diaPASEF is close to correct.**
   - The TDC decoys model genuinely unrelated false targets.
   - They do not model near-homologs, and neither does any standard decoy scheme. DIA-NN and FDRBench-style analyses
     treat homologous entrapment the same way.
3. **`-pg_scores` is at least as well calibrated as the default on this benchmark,** and gives +22 % IDs at 1 % true FDR
   (strict set: 10,800 → 13,197).
4. **Caveats:**
   - The Hamming > 2 cut was chosen after looking at the enrichment, so it is a post-hoc diagnostic, not a frozen
     criterion.
   - This is one run and one seed.
   - The cut also removes some genuinely false near-homolog hits, so the strict FDP_c is a lower estimate of the
     all-entrapment rate.
   - The two numbers (1.0-1.2 % strict, 1.5-1.7 % all) bracket the true rate for unrelated false targets.

## Consequence for option (b)
- There is no large calibration defect to fix in the decoy model or in q.
- What remains is how OpenDIA should treat near-homolog matches: a peptide differing by one residue from a present
  peptide.
  - That is a peptide-inference and specificity question.
  - A run-derived remedy would be a "sibling" competition: a target loses if a library peptide that shares most
    fragments scores higher at the same apex. OpenDIA already has an interference graph that could carry this.
- **Proposed next frozen plan:** PLAN_PGENT2, with the benchmark fixed before any new run:
  - entrapments: Arabidopsis, homolog-filtered at Hamming > 2 on the full library (the rule above, now fixed);
  - criterion 1: FDP_c at q ≤ 0.01 ≤ 0.012;
  - criterion 2: N@FDP_c 1 % ≥ 1.05 × CUR;
  - D design, H held out, Astral restated;
  - a fresh library build is needed so the homolog filter is applied at construction.
