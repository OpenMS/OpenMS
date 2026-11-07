> **Archived 7 November 2026.** This laptop-first WP-00 to WP-10 plan is closed and superseded by the [OpenDIA rewrite master plan](../../../README.md). It is kept unchanged for provenance. Its native-input policy, read-only inputs and gate discipline carry over.

# WP-07: Early localization experiment and family scoring

- **Track:** science.
- **Status:** proposed. **Owner / reviewer:** unassigned.
- **Parent epic:** [#10661](https://github.com/OpenMS/OpenMS/issues/10661). **Issue:** [#10669](https://github.com/OpenMS/OpenMS/issues/10669).
- **Shared policies:** [architecture](../architecture.md) and [plan](../README.md).

## Objective

Test the localization null/decoy model before freezing contracts, then implement shared-evidence family scoring informed by that experiment. The 27 October owner clarification also authorizes replacing identification/confidence scoring, including learned models, under the [scoring exploration protocol](../performance-strategy.md#scoring-exploration-is-allowed-including-replacement). This is planned research, not a validated faster model.

## Group proposal and scoring boundary, 5 November 2026

The grouped peak-picker experiment belongs upstream as a candidate/boundary proposal, not as a replacement for transition evidence or the family scorer. Score each retained candidate from its original fragment observations, including missingness, weak and discordant transitions, interference, RT/IM separation and shared/site-discriminating evidence. Test whether early group filtering reduces scoring visits while preserving exhaustive-reference candidate recall and independent identification/localization quality. A faster picker alone, fewer candidates alone, or fewer nominal IDs is not acceptance.

## Scope

- Start with an exhaustive small-family reference, IPF `PyProphet_inference_test.osw` and known-site data.
- Compare proposed localization nulls/decoys against reference inference, measuring localization error separately from identification error.
- Record competing/coexisting forms, unresolved groups and quantitative ambiguity rather than forcing a winner.
- Test shared-observation attribution, weak site-discriminating ions, optional MS1 and RT/IM evidence before optimizing sparse accumulation.
- New family science is single-run in v1. New multi-run inference and experiment-wide FDR claims are deferred.
- Compare batched/fused existing numerical scoring with one compact learned alternative before adopting a model. Distinguish expensive raw access and handcrafted feature generation from classifier fitting/inference; adding a neural network after unchanged feature computation does not remove that computation.
- Keep training/tuning/evaluation split by appropriate peptide/family, shared spectrum and run/instrument groups. Validate independent identification/localization error and quantitative quality, including model preparation cost and CPU memory. A new model need not reproduce legacy scores, but nominal q-values or higher ID totals alone do not establish acceptance.

## Dependencies

- **Early experiment / local experiment merge:** none. Existing fixtures and reference code suffice, independently of WP-00/01/06.
- **Final implementation integration:** WP-05 and WP-06 provide calibrated evidence and candidates.
- Experiment review gates later WP-01 contract freeze, not initial v0 closure. Keep experiment and final-integration milestones separate to avoid a cycle.

## Local merge criteria

- [ ] A reproducible experiment records input checksums, null/decoy definitions, seeds and exhaustive reference outputs for IPF and known-site data.
- [ ] A report compares localization and identification errors separately against predeclared thresholds, including unresolved cases and observed failures.
- [ ] Duplicating references to shared observations does not increase evidence weight, and exhaustive versus sparse scores agree within declared tolerances.
- [ ] Weak discriminating-ion, no-MS1 and overlapping RT/IM fixtures yield the expected localized, coexisting or unresolved outputs.
- [ ] Group proposal plus downstream scoring preserves candidate recall and reference quality for weak/missing fragments, a single strong interferer, delayed coelution, overlapping peaks and mobility-separated evidence; report scoring-work reduction separately from scientific quality.
- [ ] Before contract freeze, review records accepted/rejected model assumptions and required contract revisions, not an unsupported calibration claim.
- [ ] For the later implementation slice, WP-05/06 adapter tests reproduce reference scores and ambiguity across family page sizes without dropping alternatives.

**Integration gates:** WP-09 owns final public family-workflow checks, and WP-10 owns representative scientific/resource/performance acceptance.
