# A+B fix, step 1: extended RT calibration (search level only; no scoring yet)

**Cause** (`DIAG_B_RT.txt`, `DIAG_A_POP.txt`): OpenDIA's `RTMap` holds the end knot's run RT constant outside the calibration knots. The knots cover library RT 0.098–0.752 on raw 01, but the library spans −0.049–0.940, so 44 % of library precursors (266,799 of 600,222) get the end value. About 70 % of A+B misses are such precursors: 265 / 261 of 351 / 331.

**Change** (`extend_cal.py`, label-free): add one knot at each library RT end, on the line through the end knot and the 3rd knot from that end. Inside the knots the map, the tolerance and the ms2/im blocks are unchanged. Decoys share their target's library RT, so the change is the same for both. Search options are the same as the base pipeline (`search_rtx.sh`).

**Result** (`RTX_CHECK_raw0{1,2}.txt`; DIA-NN and entrapment used only to measure):

| Raw | Search | Population precursors T / D | Entrapment share of T | A | B | Has a window at the DIA-NN RT |
|---|---|---|---|---|---|---|
| 01 | old | 93,705 / 92,997 | 0.4913 | 185 | 166 | 2392 |
| 01 | extended | 82,934 / 82,167 | 0.4908 | 137 | 90 | **2516 (+124)** |
| 02 | old | 95,079 / 94,385 | 0.4921 | 151 | 180 | 2400 |
| 02 | extended | 87,088 / 86,339 | 0.4906 | 117 | 100 | **2514 (+114)** |

- The population shrinks by about 11 % and stays balanced between targets and decoys, with the entrapment share unchanged.
- **Not shown yet:** whether the recovered precursors are accepted. That needs the scoring pipeline rebuilt on the new search:
  - the OpenDIACNN `-stage2` features (about 13–20 min and about 25 GB RSS per raw);
  - the 6- and 12-fragment `.xipm` exports and arrays (about 45 min per raw);
  - then D11FWT under a frozen plan.
