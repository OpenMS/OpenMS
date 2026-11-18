# PLAUS_D: are the extra IDs that `-pg_scores` adds on HeLa D, and that DIA-NN does not report, real? (diagnostic)

Written 18 November 2026, at the user's request. This is not a plan, and the default stays off.

`-pg_scores` was restored as an opt-in flag (OpenMS a2aa8fc56c). Its binary is byte-identical to the tested one, the
default output is unchanged, and 23/23 tests pass. It is there to support this investigation.

**Inputs:**
- HeLa D, seed 0, CUR and NEW.
- DIA-NN on D at all q (its report goes up to q 0.5).
- OpenSWATH SVM on D.
- DIA-NN on the other HeLa run, H.
- Support at the same RT means the other tool's peak is within ±3 s of our apex.

## 1. Support of NEW's additions, against CUR's own non-DIA-NN IDs (`PLAUS_D.txt`)
| group | n | DIA-NN q ≤ 0.05, same RT | OpenSWATH q ≤ 0.01, same RT | either | DIA-NN 1 % on run H |
|---|---|---|---|---|---|
| CUR accepted, in DIA-NN | 14,911 | 0.992 | 0.829 | 0.995 | 0.853 |
| NEW added, in DIA-NN | 2,359 | 0.975 | 0.584 | 0.983 | 0.747 |
| CUR accepted, not in DIA-NN | 1,111 | 0.225 | 0.521 | 0.636 | 0.385 |
| NEW added, not in DIA-NN | 694 | 0.353 | 0.313 | 0.552 | 0.539 |

- Support at a shifted RT (our apex ±30 s) is 0.000 on every indicator. Coincidence at the same RT is therefore not
  chance.
- DIA-NN reports 71 % of the 694 at some q. 99 % of those are within 10 s of our apex, with a median difference of
  1.0 s, so they are the same peak.
- **DIA-NN's q for the 694** (where it reports one):

  | DIA-NN q | count |
  |---|---|
  | 0.01-0.02 | 133 |
  | 0.02-0.05 | 127 |
  | 0.05-0.1 | 76 |
  | 0.1-0.2 | 54 |
  | 0.2-0.5 | 102 |

## 2. Reading (revised with the decoy null, `PLAUS_D_NULL.txt`)
The ±30 s shift is the wrong null. Decoy winners in the same score band as the 694 also look "supported" at our apex:

| group | n | DIA-NN q ≤ 0.05, same RT | OpenSWATH q ≤ 0.01, same RT | either | DIA-NN 1 % on run H |
|---|---|---|---|---|---|
| decoys in the score band of the 694 (paired target key, decoy's apex) | 269 | 0.342 | 0.204 | 0.364 | 0.636 |
| all decoy winners (3,000 random) | 3,000 | 0.180 | 0.091 | 0.193 | 0.526 |

- Half of these decoy apexes sit within 3 s of their target's peak. A decoy that scores well is often the real peptide's
  signal read through shared fragments. "Support at the same RT" therefore tests the peak, not the identity, and it is
  partly chance.
- "DIA-NN 1 % on run H" is uninformative: decoys reach 0.53-0.64, the same as the 694.
- **Using the decoys as the false-ID rate:**
  - Fraction true = (group − decoys) / (CUR accepted in DIA-NN − decoys).
  - On "OpenSWATH q ≤ 0.01" this gives about **20-30 % true for the 694**, against about 50-80 % for CUR's own 1,111
    non-DIA-NN IDs.
  - "Either" gives about 31 % against 45 %.
- **Conclusion:**
  - The 2,359 additions that DIA-NN confirms are real (+14.7 % of CUR's IDs).
  - Of the 694 that DIA-NN does not report, **most look false**: about 70-80 %, so roughly 500 of 694.
  - These are worse than CUR's own non-DIA-NN IDs.
  - About 500 extra false IDs among 18,639 is about 2.7 % of NEW's IDs. With CUR's share, NEW on D is clearly above its
    nominal 1 %.
  - This is the same direction as Astral raw 02, but stronger.
- So the q-values are too optimistic in the tail on diaPASEF. The extra true IDs are still real: about 2,360 plus a few
  hundred, set against about 500 false.
- **Caveats:**
  - These are proxies, on seed 0 only.
  - The decoy null assumes that a false target is supported as often as a well-scored decoy.

### First reading (superseded)
- The text below used the ±30 s null. It overstated how plausible the 694 are.
- **The 2,359 additions that DIA-NN confirms are essentially all real.** At our RT, 97.5 % have DIA-NN q ≤ 0.05.
- **The 694 additions that DIA-NN does not report at 1 % are a mix, not pure noise.**
  - 55 % have same-RT support from DIA-NN at q ≤ 0.05 or from OpenSWATH at 1 %.
  - 54 % are found by DIA-NN at 1 % on the other HeLa run.
  - Both rates are about as good as CUR's own 1,111 IDs that DIA-NN misses (64 % and 39 %), which we already accept.
  - So, per ID, the new non-DIA-NN IDs are about as plausible as the old ones.
- **Rough split of the 694:** about 40-65 % look real, depending on the indicator. That is roughly 300-450 true and
  250-400 false.
  - If 250-400 are false among the 18,639 NEW IDs, the extra false rate is about 1.3-2.1 % of NEW's IDs.
  - Our q ≤ 0.01 should allow about 186 false IDs in total.
  - So NEW on D is probably somewhat over its nominal FDR. This matches Astral raw 02 (true FDP 0.0100 against 0.0045
    for CUR). It is not complete garbage, and it is most of the reason agreement fell from 0.930 to 0.910.
- **Caveats:**
  - These indicators are proxies. A true peptide can be missed by both other tools, and agreement with another tool is
    not proof.
  - The null for "DIA-NN at a looser q" is not clean. DIA-NN also reports chance matches at q 0.2-0.5, and our own
    top decoys reach 0.53 on DIA-NN q ≤ 0.05 at the same RT. That points to dense co-eluting signal in this 5-min Evosep
    run.
  - Only seed 0 was examined.

## 3. What this means for the next step
- About 2,400 of NEW's ~2,600 extra IDs on D are real. The other ~500 false IDs come from q being too optimistic in the
  tail.
- **Why the decoys under-count them is the key question.**
  - A well-scored decoy often sits on a real peptide's peak.
  - A false target in the tail does the same, but the new scores can reward "a real peak is here" more than "it is this
    peptide".
  - The next check is per score, among winners at the same RT as another precursor's peak: does it separate the
    identity, or only the presence of a peak?
- **Candidate fixes**, each needing a precheck and a frozen plan:
  - a more conservative scorer as the input dimension grows;
  - identity-specific features such as library correlation over the picked peak, as an interaction;
  - a stricter cut on diaPASEF.
- The flag stays opt-in until a plan passes.
