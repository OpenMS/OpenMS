# DESIGN: OpenDIA predicts its own decoys (downstream of PLAN_PDEC)

Written 18 November 2026. This is a design note, not a frozen plan. It is the user's downstream goal: "we would want it
all in OpenDIA". DIALibGen-predicted decoys (PLAN_PDEC) are the interim route.

## Why
PRECHECK_PDEC found the cause of the diaPASEF excess:
- A pseudo-reverse decoy keeps its target's predicted RT, ion mobility and intensities, which do not fit its own
  sequence.
- So a self-consistent false target beats it 58-64 % of the time among strong candidates.
- Decoys predicted like targets bring that to 50 %. With 96.6 % of targets paired, the tool's own q ≤ 0.01 gives a true
  FDP of 0.69-1.13 % on r = 1.56.
- A user with one file and one predicted library should get this without a second tool.

## What exists already (OpenMS tree, all in a2aa8fc56c)
- **`PeptDeepLibraryPredictor`** (`ANALYSIS/OPENSWATH`, `WITH_ONNX`): native ONNX AlphaPeptDeep RT, CCS→1/K0 and MS2 for
  b/y z1-2.
  - Models: `peptdeep_{rt,ms2,ccs}_dynamic.onnx`, fetched at build time and checked by sha256.
  - DIALibGen runs AlphaPeptDeep models too (its own help says so), so the two are the same model family. Whether the
    weights are identical is still to be checked.
- **`OpenSwathDecoyGenerator -predict_decoys`** (7fc3d1a1cb, `predictDecoyAssaysLight_`):
  - It predicts each MRMDecoy decoy at its target's charge.
  - It maps predicted RT and IM onto the library scale with lines fitted on up to 2,000 targets (prediction against
    library value).
  - It keeps the decoy's own top fragments, as many as its target has, inside the targets' fragment m/z range.
  - It estimates NCE and instrument from the library itself: the best mean cosine to the library intensities over every
    PeptDeep instrument and NCE 20-40.
  - It drops decoys the models reject.
- **Gap:** the `opendia-pg-a2aa8fc56c.sif` container has the binary but not the model files, so `-predict_decoys`
  cannot run there yet.

## Proposed design
1. **In OpenDIAWorkflow, after loading the library:**
   - If the library has no decoys, or `-decoys predict` is set, OpenDIA builds pseudo-reverse decoys from the targets
     and predicts them with `PeptDeepLibraryPredictor`.
   - This reuses `predictDecoyAssaysLight_`, moved to a shared helper so OpenSwathDecoyGenerator and OpenDIA call the
     same code.
   - The output is the in-memory search library, so there is no extra PQP on disk.
2. **Parameters come from the file:**
   - NCE and instrument: the existing cosine fit to the library's own intensities.
   - RT and IM scale: the existing target lines.
   - Nothing is read from labels or from other runs.
3. **When to predict:** the same cosine fit says whether the library looks predicted. A high mean cosine at the best
   setting means a predicted library, where predicted decoys are right.
   - For an empirical library (DIA-NN or Spectronaut output), predicted decoys would be smoother than the empirical
     targets. That creates a new asymmetry in the other direction.
   - The default there stays pseudo-reverse until that case has its own entrapment test.
   - The threshold must be fixed in a plan, not tuned.
4. **Coverage:** one decoy per target precursor, at the same charge and with the modification moved with its residue
   (pairing v2 of `merge_predicted_decoys.py`).
   - A target whose decoy the models reject is reported and keeps a pseudo-reverse decoy, so no target competes alone.
   - PRECHECK_PDEC showed that unpaired targets were the whole remaining excess.
5. **Cost:**
   - DIALibGen predicts 1.39 M decoy precursors in 2:40 wall on 16 CPU threads.
   - OpenMS ONNX throughput is to be measured. The budget is that a HeLa D run (about 8-10 min on the plain library)
     grows by at most about 20 %.
   - Prediction runs once per library and can be cached next to it (sha256 of the library plus the models).

## Order of work
1. **Precheck:**
   - Run `OpenSwathDecoyGenerator -predict_decoys` (models staged into the container) on the PRECHECK_PDEC targets
     (r = 1.56).
   - On D seed 0, compare to the DIALibGen-predicted decoys on three numbers: P(entrapment beats own decoy | q ≤ 0.05),
     the lower and combined bounds, and N@FDP_c 1 %.
   - If they match within noise, OpenMS's predictor can replace DIALibGen. If not, find out why first (weights, NCE
     estimate, RT line, fragment selection).
2. **Frozen plan:** OpenDIA's built-in decoys on the r ≈ 7 library, with the same criteria as PLAN_PDEC. D is design,
   H is held out.
3. **Empirical-library entrapment check** before item 3 of the design (the predicted or empirical switch) becomes a
   default.

## Done already
The PQP reader dropped peptides without a gene mapping whenever GENE had rows. That silently removed merged decoys. It is
fixed in 2982bd647d, with unchanged row order and a class test that fails without the fix.

## Precheck step 1, build (18 Nov, 22:19)
`entlib-osdg-20261118`, using `OpenSwathDecoyGenerator -method pseudo-reverse -predict_decoys true` on the PGENT2
`targets.pqp`:
- **Models:** staged into `models/`, hashes in `models.sha256`.
- **Output:** `hela_osdg.pqp`, sha256 d1c44334…
- **Settings estimated from the library:** timsTOF at NCE 30, with mean cosine 1 to the library intensities. That is
  expected: the same model family predicted the library.
  - RT line: identity.
  - Cost: 1:48 wall and 4.5 GB on 16 threads.
- **Coverage:** 404,719 decoys predicted and 0 rejected by the models. But 36,927 (8.4 % of targets) were **dropped**,
  because they had fewer fragments than their target inside 200-1800 m/z.
  - Those targets compete without a decoy. This is the same failure mode PRECHECK_PDEC found at 32 %.
  - So the fallback in design item 4 is required, not optional: keep a pseudo-reverse decoy when the predicted one is
    dropped.
- **Search:** D seed 0 CUR/NEW were submitted (34769520/1). The results will be compared with PRECHECK_PDEC "Full
  coverage" using `eval_pdec.py`.

## Fallback implemented: `-predict_keep_short` (OpenMS fa75f8f46d, default off, output unchanged)
- **How DIALibGen handles it:** it never drops a decoy for having fewer fragments. On the PRECHECK_PDEC library, 8.0 % of
  its decoys have fewer fragments than their target and 7.4 % have more. That library calibrated at 0.69-1.13 %.
- **What OpenSwathDecoyGenerator now does with `-predict_keep_short true`:** the same thing. A short decoy keeps all of
  its usable fragments (at least one) instead of being dropped. A decoy with no usable fragment is still dropped.
- **Tests:**
  - The class test shows short decoys are dropped without the flag and kept with it, holding 1 to k−1 fragments.
  - Decoys kept before are unchanged.
  - 342 OpenSwath, PQP and OpenDIA tests pass.
- **Unrelated failure:** `TOPP_OpenDIA_write_ini_massim_defaults` fails with and without this change. It is not caused by
  it.
- **Next:** rebuild the container, then rerun precheck step 1 with `-predict_keep_short true`. The design's item 4 should
  then hold: no target competes alone, except for decoys with zero usable fragments.
