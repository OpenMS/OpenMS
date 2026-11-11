# Voxel v4 phase 1: window ranking term and MS1 isotope features for the single decoy-whitened network (written before any v4 training)

Context: `v3seq/MISSES_D11.md`. Against DIA-NN 1 % precursors (exact library match, human, reference only), D11 loses about 400 per raw because a wrong window wins (C), and about 700 because the right window scores below the cut (D). On C, BASE picks the right window for 332 / 307 of them. Below its top 1 %, the D11 score does not follow the features.

## Arms (one variable each, all against D11 = PLAN_V3_P2b scores, same population, folds, rounds, seeds and λ)

| Arm | Change from D11 |
|---|---|
| D11R | + within-precursor ranking term (below) |
| D11M | + 6 MS1 isotope features added to z (`v4/ms1_features.py`: log M+0 area, M+0–fragment-sum correlation, observed-vs-theoretical isotope cosine from the elemental composition, M+0–M+1 correlation, M+0 apex shift to the fragment apex, M+0 apex share) |
| D11MR | both |

**Ranking term.** For each positive row in a batch, the target windows of the same precursor whose apex is more than 10 s from the positive's apex are its siblings (at most 8, sampled). The score is s = h·w / √obj, with w and obj from the same batch readout (decoy-SD units). L_rank is the mean over positives with at least one sibling of logsumexp(s_pos, s_sib) − s_pos. Training loss: −obj + α · obj.detach() · L_rank, with α = 1 (fixed, not tuned). Early stopping uses log(obj_v) − α · L_rank,v on the validation rows (the same direction as the training gradient). For D11M, α = 0 and early stopping uses obj_v, as in D11.

## Pass rule (each arm, on both raws)

1. PLAN_V2 rules 1–4 against BASE (ΔN@2 % interval > 0; ΔN@1 % interval not entirely < 0; FDP_c at q ≤ 0.01 ≤ 0.015; entrapment AUC ≤ BASE).
2. N@2 % above D11 on both raws (the paired bootstrap interval against D11 is reported).
3. Null-winner AUC ≤ D11's + 0.005 (sequence or composition tell guard; the MS1 envelope depends on composition).

**Reported only:**
- Funnel against DIA-NN (`v3seq/misses.py <arm>`): stage C (wrong window) and stage D (below the cut) counts.
- N at 0.5/1/1.5 %.
- Entrapment AUC.
- Epochs and time.

An arm that passes is a candidate for the C++ `OpenDIAScoringModel` under a new plan, and the base of phase 2.

## Phase 2 (planned here, frozen in its own plan before its runs)

- The 12-fragment library (DIALibGen 0.13.1, `astral-lib12-20261111/`; at most 14 threads) adds up to 6 extra fragments per precursor.
- Search and candidate windows stay those of the 6-fragment library. The extra fragments are extracted only for the existing windows, through an `-out_xipm` export.
- Then: extra-fragment co-elution and 12-fragment library-fit features.
