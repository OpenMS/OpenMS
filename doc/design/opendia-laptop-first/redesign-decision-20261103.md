# OpenDIA system redesign decision, 3 November 2026

This decision responds to the latest full-workflow A/B result and resets the implementation direction. The objective remains a fast, memory-efficient DIA search engine for large predicted libraries, with native Bruker and Thermo inputs. This is an architecture and sequencing decision, not a claim that the proposed engine is already faster or scientifically equivalent.

## Decision

Stop promoting isolated allocation, cache, or adapter changes as the main optimization path. The latest three-run A/B for the UIS score-vector reservation change produced median wall times of 399.26 s (parent) and 394.20 s (candidate), a nominal 1.3% difference with overlapping run ranges. Median user CPU moved in the wrong direction, from 1,466.62 s to 1,511.05 s. All six output databases passed integrity and had matching feature/score row counts. This is useful evidence that the tiny change is safe, but it does not establish a meaningful speedup. Keep the code change only if desired for its small allocation reduction; do not spend more cycles repeating this A/B.

Build a separate, private OpenDIA search core under `src/topp/OpenDIA/`. Retain OpenSWATH as the behavior/reference path and reuse vendor decoding and calibration where they do not force the new engine through expensive compatibility containers. Do not rewrite codecs, require converted mzML, or change shared OpenMS APIs just to enable the redesign.

## What the current path tells us

The current tool performs substantial orchestration in `OpenDIA.cpp`, then calls `OpenSwathWorkflow::performExtraction`. The latter materializes transition chromatograms and routes them through the existing `MRMFeatureFinderScoring`/peak-picking pipeline. The private calibrated-region adapter changes how spectra are read, but still services chromatogram extraction and leaves the downstream object-heavy scoring path in place. Thus copy avoidance and a micro-optimization in an existing scorer are not a change to the central search dataflow.

The code and existing profiles do not justify claiming that one particular formula, model, or decoder is the sole bottleneck. The architecture should remove repeated work by construction, and its first vertical slice should count candidate/fragment visits and materialized evidence bytes so that the algorithmic change is directly observable. Do not run another general timing campaign before that slice exists.

## Target dataflow

1. **Compile the library once.** Convert empirical or predicted input into a compact, versioned hypothesis store. Keep stable precursor, peptidoform-family, target/decoy and provenance IDs; sparse fragment masses/intensities; and RT/mobility predictions with uncertainty. Index precursor hypotheses by m/z and acquisition eligibility. Store family membership and shared versus site-discriminating fragment incidence without expanding every pair of forms into a UIS list. Page cold metadata instead of copying it into every worker.
2. **Stream native acquisitions in bounded tiles.** Reuse native readers/calibration, but make an acquisition tile the unit of work. Decode each tile once, preserve exact scan/frame, isolation, RT, calibration, discrete FAIMS and mobility identity, and pass it to all eligible candidate queries before releasing it. Batch by acquisition, not lane-major accessor clones. Use a single byte budget, bounded queues, backpressure and cancellation.
3. **Join observations to candidate hypotheses.** Use indexed precursor and fragment-m/z interval queries, with RT/isolation/IM constraints, instead of creating and retaining a dense chromatogram object for every transition as the internal representation. Emit sparse evidence only for observed candidate support, retaining missingness and exact coordinates. Query bins may accelerate lookup but never substitute for measured coordinates. Preserve a complete/recovery path for hypotheses whose predicted anchors are absent or weak.
4. **Generate and score local candidate families.** Build competing peptidoform/charge candidates from the sparse evidence, jointly accounting for shared and discriminating fragments, alternatives, nulls and completeness. Batch compact contiguous score inputs. Start with a correctness/reference engineering scorer on the new evidence representation; evaluate any scientifically different compact model as a separate, versioned lane. Do not run the old expensive feature suite and add a model on top if the new model is intended to replace it.
5. **Finalize statistics and outputs after the complete population is known.** Keep target/decoy populations, rejected candidates and local alternatives through confidence estimation. Do not compute tile-level q-values or silently truncate candidates. Emit compact results by default. Provide an explicit compatibility/rescoring adapter for full evidence when requested, and keep existing OSW output through the reference path until the new output contract is independently accepted.

## First implementation slice

Implement one private, end-to-end vertical slice that takes one native acquisition tile plus a compact candidate index and emits sparse fragment evidence for eligible candidates. The core acceptance check is that the same tile/query set yields the same extraction values, axes, missingness and candidate identities as the current extractor for supported tophat cases, including IM and isolation boundaries. Add synthetic adversarial cases for overlapping m/z intervals, duplicate centers, empty scans, missing mobility, calibration changes and RT boundaries.

The slice must report deterministic counters for tiles decoded, hypotheses visited, fragment queries tested, evidence points emitted, bytes retained and work-queue high-water. Its algorithmic acceptance condition is that unrelated library entries are not scanned for every raw point and that retained data scales with observed sparse evidence plus bounded workspace, not all candidates times all scans. Keep OpenSWATH as an executable reference during this increment. Do not make a speed claim from counters alone.

Then wire this evidence to a candidate-family scoring batch and real output path. Only after the end-to-end native route exists, perform a small, predeclared cold-run comparison on the frozen TimSim and sampled-100K Astral scopes. Include preparation, calibration, search, statistics and writing. Stop if the new path fails science/resource limits; do not broaden into an unbounded parameter or benchmark matrix.

## Milestones and gates

- **R1: Indexed sparse extraction.** Private records and index, one-tile indexed interval join, exact reference tests, mobility/isolation/boundary coverage and counters. No new public OpenMS API.
- **R2: Family candidate scoring.** Competing-family formation and batched scoring over sparse evidence, complete target/decoy and missingness semantics. Keep old scorer selectable as a reference.
- **R3: Native end-to-end execution.** Bounded scheduling/admission, cancellation, output publication and rescore compatibility. Bruker and Thermo remain native inputs, with source-specific capabilities behind private contracts.
- **R4: Scientific and scale acceptance.** Independently validate identification, localization and quantification. Measure complete cold/warm cost, process/cgroup memory and disk on frozen representative scopes. Compare external tools only with matched search scopes and properly validated error control. No competitive claim until these gates pass.

The existing WP-00 evidence and public acceptance baseline remain frozen. WP-01/02/03/05 work is reusable groundwork, not proof that the resulting engine is fast. Update package ownership and integration status separately as each milestone lands. Do not merge the private redesign to `develop` until combined correctness, speed and memory gates are satisfied.

## Explicitly deferred

No further isolated cache tier, decoder rewrite, mandatory conversion format, large neural model, GPU dependency, or repeated timing of the one-line scorer change. Revisit only if the indexed sparse path exposes a concrete need and the complete workload can pay for it.
