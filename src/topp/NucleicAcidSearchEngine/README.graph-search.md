# Graph-guided RNA candidates and route scoring

NASE keeps `candidates:method=enumeration` and `scoring:method=hyperscore`
as defaults. Graph candidate generation and route scoring are independent,
experimental options. The full-spectrum MVH and HyperScore implementations
remain in OpenMS core.

## Candidate modes

| `candidates:method` | Behavior |
| --- | --- |
| `enumeration` | Existing precursor-mass search and modification enumeration. |
| `graph_audit` | Score enumeration candidates and annotate graph support. |
| `graph` | Filter enumeration candidates by validated RNA tags. |
| `graph_generate` | Prune modification placement with conservative mass/tag bounds. |
| `graph_generate_audit` | Compare eligible placements against exhaustive graph filtering; fail on a mismatch. |
| `graph_routes` | Generate spectrum-specific placements from database-compatible anchored routes. |

Graph edges retain every isobaric residue label and ordered two-residue
alternative. Separate fragment series, orientations, terminal origins, and
charge hypotheses are validated against `NucleicAcidSpectrumGenerator` output.
One- and two-nucleotide anchors and an optional missing intermediate rung keep
missing early fragments from eliminating every route. A mass-difference path
alone never establishes a sequence assignment.

The tag/filter modes preserve enumeration fallback for absent evidence,
unsupported models, or exceeded graph/query bounds. `graph_generate` also
falls back to exhaustive placement if its per-oligo placement bound is exceeded.
`graph:verify_matches` compares the optimized matcher with the reference.

`graph_routes` instead combines routes against database prefix/suffix tries,
completes modifications, and validates exact precursor masses. Its explicit
state, completion, candidate, and support limits fail the entire search rather
than export partial results. `routes:verify_candidates` checks complete
candidate identities and scaffold witnesses against the unindexed oracle;
`routes:precursor_index=false` selects that oracle directly.

## Reuse and RNA-specific code

| Responsibility | Existing OpenMS component |
| --- | --- |
| Sequence identities, formulas, terminal groups and modification masses | `NASequence`, `RibonucleotideDB` |
| Modification-dependent RNase cleavage and parent mappings | `RNaseDigestion`, `ModifiedNASequenceGenerator` |
| Fragment masses, charges, orientations and annotations | `NucleicAcidSpectrumGenerator` |
| MS1 correction and MS2 preprocessing | Existing FLASH integration, `RNADeisotoper`, window and intensity filters |
| Native spectrum scores | `MetaboliteSpectralMatching::computeMVHScore` and `computeHyperScore` through one NASE dispatch |
| Longest connected ordinal series | `SpectrumAnnotator::computeLongestIonSeries`, shared with peptide annotation |
| Sorted mass-window access | `MSSpectrum::MZBegin` |
| Hit competition, FDR and exports | Existing NASE/IdentificationData, OMS, idXML and mzTab paths |

The peptide `Tagger` cannot directly represent RNA residue alternatives,
modification-aware terminal chemistry and bounded missing-rung paths. Closest-only
spectrum alignment also cannot supply every competing peak/charge interpretation
needed by route validation. The RNA DAG, candidate completion, diagnostics and
anchored route dynamic programming therefore remain tool-local. The shared
ordinal counter supplies an admissible bound for route scoring; anchoring,
peak reuse and scaffold consistency are still enforced by the RNA scorer.

## Preprocessing and modification evidence

MS1 correction is selectable separately from MS2 deisotoping. To use FLASH
for MS1 and legacy processing for MS2:

```sh
NucleicAcidSearchEngine -in spectra.mzML -database targets_decoys.fasta \
  -out results.mzTab -id_out results.idXML \
  -precursor:correct_with_flash -preprocessing:deisotoping legacy \
  -candidates:method graph_audit -out_graph graph.tsv
```

Raw MS2 diagnostic evaluation precedes noise filtering and deisotoping.
`diagnostic:fixed_ions` supplies polarity-specific `code=mz` rules for labile
residues; other generic diagnostics use existing single-nucleotide chemistry.
`diagnostic:filter` enables gating. `diagnostic:filter_scope=explicit` gates
only explicit rules, allowing intact graph-supported modifications without a
terminal monomer. Generic charge hypotheses above one require a raw isotope
partner. Diagnostic permission does not establish modification localization.

The default `routes:completion=bounded` permits configured variable modifications
at uncovered positions without requiring a diagnostic reporter. It limits variable
placements without intact graph evidence to `routes:max_unanchored_modifications`
(default two). `labile` instead keeps uncovered sites canonical or permits only
diagnostic-supported labile modifications. `bounded_labile` applies that limit to
labile completion only; `all` enumerates permitted completions under resource
limits. Graph-placed modifications still consume the full per-oligo budget.
The default per-spectrum limits are one million route states, five million candidates,
twenty million supports and one billion completion nodes. These are failure bounds,
not preallocated buffers; exhaustion aborts without output and names the option to
raise. Bounded completion can be substantially slower than enumeration. It still
requires an anchored route and is not guaranteed to reproduce enumeration.

The default modification-aware cleavage model uses `RNaseDigestion` to validate
each completed product against its parent and missed-cleavage budget.
Cuts created by a modification (e.g. RNase T1 after inosine) count against
that budget, so such a product at the full budget is dropped here although
enumeration keeps it. A cached digest supplies its own missed-cleavage budget.
`oligo:cleavage_model=canonical` is an explicit alternative for `graph_routes`.
Cached digests must match the database and digestion model and come from
`digest_out`; files containing search observation matches are rejected. Enumeration/tag modes
reject route digests because cleavage-sensitive placements have not been pre-expanded. Route generation
requires monoisotopic masses. By default (`routes:terminal_states=enzyme`) each
product receives the terminal groups of the configured RNase at its parent
position, as in enumeration; a parent match is retained only if the candidate's
termini agree with it. Deferred scaffolds are held outside the output molecule
store so exporting a hit cannot restore incompatible parent mappings.
A fixed terminal modification fills only an end without
an RNase group, again as in enumeration. `routes:terminal_states=all` additionally tries all six
physical assignments for every oligonucleotide: OH/OH, OH/P, OH/cyclic-P, P/OH,
P/P and P/cyclic-P. Linkage modifications (codes ending in `*`) have
neighbor-dependent masses and are rejected. Explicit input/fixed termini
remain constrained; variable terminal/linkage modifications use enumeration/tag
modes.

## Scoring

`scoring:method=route` ranks coherent anchored ladders using
`n * (1 + matched_intensity / local_intensity)`. The longest route per ion
series wins, with Gaussian mass-weighted intensity breaking length ties. The
best prefix/suffix pair receives one node/intensity credit per original peak,
including alternative charge hypotheses. Models never mix within a ladder;
prefix and suffix may use distinct diagnostic-supported labile scaffolds.
The reported precursor retains its full modification composition.

Route matching retains alternative prefixes with distinct still-relevant used
peak identities. A peak identity expires only after its last possible rung;
the full route is reconstructed from backpointers for scoring and annotation.
`routes:max_states` also bounds retained scoring prefixes per model query.
Exhaustion fails the search instead of silently selecting an incomplete result.

Defaults are a 50-Da local radius, ten local peaks, two matched nodes per route,
and at most one missing intermediate rung. `scoring:route_mass_accuracy_tolerance`
can separate Gaussian weighting from the matching tolerance; zero inherits it.
`scoring:route_base_loss` adds coherent a-B prefix ladders if a-B ions are enabled.
These are experimental NASE implementation choices. Scores are rounded to eight
decimal places and all top score ties are retained.

For `graph_routes` with MVH/HyperScore, `scoring:labile_scaffolds` instead selects
the best complete intact or diagnostic-supported scaffold model using the
existing native score. It never pools incompatible models. The option is off
by default. `native_scoring_scaffold` records the selected fragment model.

Graph and fallback candidates pass through the same hit competition and native
FDR path. All candidate, diagnostic and scoring rules treat targets and decoys
identically. A route score or nominal q-value yield is not independent evidence
of identification accuracy.

## Validation and reproducibility

`out_graph` exports graph size, traversal counts and completion/fallback status.
Hit metadata records tag witnesses, diagnostic permissions, route node counts,
intensities and scaffold identities. Tests cover graph/placement oracle agreement,
resource exhaustion, cleavage-dependent modifications, exact score ties, peak
identity, charge hypotheses, scaffold consistency, native scoring and digest
cache provenance. The shared ordinal counter is checked exhaustively against
an independent small-space reference and is exposed in pyOpenMS.

`doc/benchmarking/nase/hap1_legacy_sensitivity_corrected_2026.json` retains the
exploratory legacy HAP1 recipe, input hashes and historical results. Its explicit
settings are experimental and do not change the defaults. A fresh replay
of clone1 replicates 1–3 after the reuse refactor preserved every exported winning
assignment, score and native q-value, and byte-identical graph/precursor diagnostics.
Combined yields were 5,454 / 7,587 / 8,422 PSMs at 1% / 5% / 10% native q-value.
This three-file replay does not constitute an independent accuracy benchmark or
a new nine-file comparison against another engine.

A focused follow-up with the revised route defaults is recorded in
`doc/benchmarking/nase/hap1_route_defaults_review_2026.json`. On one complete
HAP1 recording, enumeration/HyperScore gave 546 / 1,154 / 1,711 target-containing
spectra at 1% / 5% / 10% native q-value; default graph routes with HyperScore gave
546 / 1,177 / 1,708. All baseline accepted calls at 1% and 5% were retained.
The graph run completed without resource overrides but took 136 seconds versus
15 seconds for enumeration (16 threads, this machine). On its RT 35–39 minute
slice, default graph routes reproduced enumeration's 109 / 143 / 154 accepted
calls and assignments. The optional route score gave 169 / 181 / 191, but changed
27 of the 109 baseline assignments at 1%; increased nominal-q yield is not
independent evidence of accuracy. This single-recording comparison has no
entrapment control and does not validate general transferability. Historical
1% yields also varied strongly over decoy seeds; interpret the recorded 5%/10%
results and seed spread together. Enumeration and HyperScore remain the defaults.
