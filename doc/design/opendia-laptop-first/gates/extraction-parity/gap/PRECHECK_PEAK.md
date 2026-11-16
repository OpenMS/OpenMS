# PRECHECK_PEAK: room in the wrong-peak / missed-peak losses (16 November 2026)

Diagnostic, written before any plan.
- **Inputs:** arm-E arrays and the C++ scorer at seeds 0-4. The seed 0-2 rescore also writes every window's raw score
  and linear part; it reproduces `cppscore/` byte for byte.
- **Reference:** DIA-NN 2.2 (same targets, own decoys), human IDs at q ≤ 0.01, used only to annotate the right peak
  (apex within 6 s of DIA-NN's).
- **Scripts:** `peak_split.py`, `wrule_whatif.py`, plus the inline script behind `CAP_DETAIL.txt`.

## 1. Where DIA-NN's peak is lost before scoring (`PEAK_SPLIT.txt`)

| cause | raw01 | raw02 |
|---|---|---|
| no search candidate at all | 0 | 0 |
| candidates, none within 6 s, DIA-NN apex **outside** our RT window | 49 | 63 |
| candidates, none within 6 s, apex **inside** the window (picker / 5-peak cap) | 47 | 37 |
| right candidate found but the pair is dropped by the population rule | 3 | 1 |
| right candidate in the population, **another window kept** (W rule) | 37-40 | 34-37 |
| right peak kept | 2205-2208 | 2171-2174 |

**Inside-window misses.**
- 96 % / 97 % of them sit at the 5-peak cap.
- DIA-NN's quantity for these IDs is 0.4× the median, so they are weak peaks.
- Raising the cap would bring in more candidates for every precursor: 83-88 % of all precursors already exceed 5 peaks.
  The extra candidates would add chance windows on both sides.

**Outside-window misses.** These are the extrapolated gradient ends. Two plans already failed on them (PLAN_RTPROFILE
and PLAN_CALREFINE2).

## 2. Is the W rule choosing badly? (`WRULE_WHATIF.txt`, seeds 0-2)
The what-if keeps the window with the highest full score instead of the largest linear part.

| | right peak kept, raw01 | right peak kept, raw02 | N@1 % raw01 | N@1 % raw02 |
|---|---|---|---|---|
| W rule (now) | 2205-2208 | 2171-2174 | 1822 / 1871 / 1860 | 1835 / 2012 / 1896 |
| max full score | 2049-2062 | 2066-2084 | 1839 / 1867 / 1839 | 1837 / 1967 / 1884 |

- The W rule picks the right peak far more often than the full score does: 146-159 more on raw01, 87-106 more on
  raw02. It is the better rule, and N moves only within seed noise.
- Of the kept right peaks, 97 % are the picker's top-ranked candidate by co-elution (`CAP_DETAIL.txt`).

## Reading
- The losses before scoring total about 100 per raw.
- Of these, 49-63 come from the RT window ends, which two earlier plans failed to fix.
- About 40-47 are weak peaks behind the 5-peak cap.
- About 35-40 are W-rule choices, and the W rule is already the better rule.
- **No cheap fix is left here.** A cap change could recover at most the 37-47 inside-window weak peaks. Even those must
  then beat the cut, which most weak right peaks fail (the 570 / 463). The expected net gain is about +10-20, below the
  seed noise.
- **Not pursued as a plan.**
