# PLAN_ASTRAL_FULL, 23 November 2026

## Purpose and status

User-requested full predicted-library capacity benchmark, not a scoring intervention. OpenDIA has strong experimental identification gains on HeLa but unestablished q calibration, unvalidated quantitative accuracy and a large compute gap. The small Astral library is not evidence of full-proteome scalability. Freeze and push this plan and scripts before any new library/search outcomes. No OpenMS edits, parameter tuning, seed selection, new default or revised earlier verdicts.

## Inputs and provenance

Codon root `/nfs/research/juan/singj/opendia`, work directory `astral-full-20261123`.

- All three native RAWs: `proteobench-astral-entrapment/transfer/raw/LFQ_Astral_DIA_15min_50ng_Human_{01,02,03}.raw`.
- Full unsampled peptide FASTA: `proteobench-astral-entrapment/transfer/full_fasta/ProteoBenchFASTA_Entrapment_Human_with_contaminants_entrapment_pep.fasta`, 2,841,736 records by metadata inspection. Do NOT substitute 100k sampled FASTA, filter by identifications, or use observed RT/spectra to refine it. This is a peptide FASTA, not protein FASTA to redigest.
- GPU generator: `entlib-hela-20261118/dialibgen-cuda-0.13.1.sif`, DIALibGen 0.13.1, pinned SHA256 in execution manifest.
- OpenDIA: `opendia-pd10-fastint2.sif`, code `3b1d90ed67`, previously byte-checked 9/9. Use unchanged default GBM/network/RT/peak selection and `-decoys library` (default) after one prepared library is frozen.
- Reference context: the supplied ProteoBench DIA-NN 2.3.1 metadata at local dataset `proteobench_diann_v231_res_20260717_123332/param_0..txt` reports 8,649,692 library precursors, peptide7-30, charge1-4, CAM C, no variable mods, no methionine excision, precursor300-1800 and fragment200-1800. Reference used MBR/reanalysis and different hardware. It is not a matched speed or single-run ID comparator, and its nominal-q IDs are not entrapment-calibrated ground truth. No new DIA-NN run is authorized by this plan.

## Frozen generation

One A100, 16 CPU threads, 384 GiB allocation, 12h limit. DIALibGen CUDA with require_gpu, Astral/Lumos model, NCE25, no cleavage, missed cleavages0 (input already peptides), peptide7-30, charge1-4, precursor300-1800, fragment200-1800, CAM C, empty variable modifications, no methionine excision, fragment charge<=2, **6-20 fragments per target precursor**. Full means every eligible precursor under this recipe, not every mathematically possible charge/modification. Actual count need not equal DIA-NN because predictors and eligibility differ.

Generate targets only. DIALibGen's conventional decoys are not own-sequence predictions and can corrupt PQP peptide flags, so do not export them as if they were the predicted decoys validated in OpenDIA. Next use unchanged OpenSwathDecoyGenerator from fastint2, CPU16, 512 GiB, 24h, pseudo-reverse with its default K/R switch, predict_decoys=true, keep_short=true, NCE25/Lumos. This decoy-preparation step is CPU, explicitly separate from GPU library generation, and its cost counts in end-to-end reporting. No claim that the entire preparation is GPU-accelerated.

Preserve master predicted library capped at20, plus search6 and scoring12 views with identical precursor IDs and same target/decoy membership, selecting by predicted intensity with deterministic ties. Scoring numeric transition IDs must encode descending intensity because current loader takes the first12. Do not silently drop unmatched targets, collision cases, or short-decoy pairs, regenerate decoys to improve results, or alter default decoy construction. Record missing partners and fragment/sibling asymmetry as library limitations. Stop for invalid schema, inconsistent IDs/flags, >20 fragments, or unsupported writer behavior. Current algorithm consumes6 search and up to12 scoring fragments: this is NOT a 20-fragment scoring experiment.

## Frozen search population and resources

RAW01/02/03 x seeds0,1,2,3,4 =15 independent executions, default scoring,16 threads each. Capacity allocation256GiB and24h per run, at most3 concurrently. This deliberately measures full-library capacity rather than artificially claiming a32GiB capped run. Report any observed RSS>32GiB as failure of the laptop memory target. Observed RSS<=32GiB in a256GiB allocation is not a proven32GiB-cap acceptance test. Do not quietly downsample, disable GBM, reduce library size, add MBR, or increase resources after a failure. Infrastructure retries require a written record and exact same scientific configuration, not selection of the best completed attempt. OOM, timeout or calibration failure is an operational failure and remains in the denominator15.

No all-window features dump in15 capacity runs because it materially changes I/O/storage and these are timing measurements. Exact candidate-window null-symmetry diagnosis remains a separate missing-artifact task, not falsely claimed from precursor winner tables. No protein identification/quantification claims from peptide-FASTA accessions that label individual peptides as proteins.

## Frozen endpoints, no post hoc pass bar

1. **IDs:** human and entrapment precursor winners and unique human modified peptide sequences at own q grid0.1,0.25,0.5,0.75,1,1.5,2,3,5%; all seeds shown, median/min/max per raw. Reuse complete tool output, never recompetition, no best-seed reporting. Distinguish total targets including entrapments from human IDs.
2. **Calibration:** E/N and E*(1+1/r)/N at every grid point, r measured on the prepared search-library target population. Undefined for no discoveries. These are descriptive entrapment estimators under assumptions, not exact realized-FDP brackets or proof of FDR control. Same-raw seeds are not independent replicates, no pooled Poisson confidence or simultaneous curve claim. Paired estimator remains blocked unless original generator map is authenticated under a separate preregistration. Entrapment labels only grade and define explicitly reported post-search quantification strata, never feed scoring, q estimation, peak selection or selection of normalization methods.
3. **Quantification:** companion QUANT_PLAN and synthetic-tested evaluator define no-imputation area completeness, overlap, across-three-run CV and log2 correlations at q<=0.5%,1%,2%, reported for human-only and all targets separately for each seed. Raw and explicit median-normalized areas, shared-feature conditioning stated. These50ng human runs offer precision/reproducibility, not established fold-change accuracy, interference freedom or absolute accuracy. If matching technical-replicate provenance cannot be authenticated, call CV cross-run dispersion only. No pooled15-run CV.
4. **Resources:** per-job hostname/CPU model, allocated CPUs/memory, /usr/bin/time wall/user/sys and MaxRSS, Slurm state/exit/elapsed/MaxRSS, output/library bytes. Separate queue delay, GPU prediction, CPU decoy generation, view construction, run-only cost and total elapsed/CPU-hours. Note caching, shared-filesystem/node variability and `%3` concurrency. Include failed attempts. Record GPU type/driver and require_gpu evidence. Timing comparisons to existing differently-sized libraries or unmatched DIA-NN runs are contextual only.

Operational completeness requires15/15 successful runs, valid structural checks, and all requested endpoint tables. No numerical ID/FDR/quant quality pass threshold is invented from these outcomes. A completed benchmark does not mean statistical validity is proven. Report failures plainly.

## Execution and gating

SHA256 all scripts, generator/tool containers, RAWs, full FASTA and all library artifacts. Run synthetic tests and bash -n first. Commit/push this plan before submitting prediction. Chain GPU prediction -> CPU preparation ->15-run array ->evaluation with afterok dependencies for preparation and afterany for reporting search failures. Save job IDs and scheduler states. No OpenMS push. Run board_sync.py --archive-stale before each docs commit. Preserve every original log and output, do not overwrite failed attempts.
