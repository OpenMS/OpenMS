# OpenDIA private records, unstable v0

These contracts implement roadmap steps 2, 3 and the step 4 private record shapes.
They are tool-local, not installed, not a public OpenMS API, and **not frozen**.
Scientific freeze is pending the WP-07 experiment review. Family inference,
scientific scoring and production adapters are deliberately not implemented here.

## Extraction-query identity

`ExtractionQuery` is a checked plain value key, not an extractor, storage layer or
cache integration. `specification()` returns a detached value. All nested equality
is standard defaulted equality: exact finite doubles, signed zero equal, adjacent
representable doubles unequal. Equality is conservative, not proof of scientific
interchangeability. No fuzzy bins, unit conversions or inferred calibration apply.
Caller calibration snapshot IDs must name immutable snapshots in their namespace.

The specification includes context, acquisition selection, nonzero MS level,
optional absolute isolation m/z interval, optional RT interval in seconds,
absolute m/z interval in Th, optional continuous mobility interval with explicit
existing unit, optional discrete FAIMS voltage, full-width nonnegative tolerance
in Th or ppm, aggregation and isotope policy. Absent isolation means no isolation
restriction, absent RT means all RT, absent mobility means no mobility restriction,
and absent FAIMS means no voltage restriction, not 0 V. Acquisition selection is
an optional set of opaque `(frame_id, scan_id)` pairs, not ordinal ranges. Absent
means all acquisitions, present empty means none. Only this set is sorted and
deduplicated. Neither caller inputs nor getter copies can mutate a retained key.

Each interval explicitly selects OPEN, CLOSED, LOWER_CLOSED or UPPER_CLOSED edges.
Bounds must be finite and ordered, with equal bounds allowed. LOWER_CLOSED includes
only the lower endpoint and UPPER_CLOSED only the upper endpoint. Existing tophat
extraction uses CLOSED RT and OPEN m/z and mobility edges. This unstable v0 supports
only the fixed tophat profile. Unsupported filters require a contract extension,
not a speculative enum or silent fallback. RAW_POINTS retains individual points
and SUM requests their intensity sum over the selected support. These are request
semantics only, with no algorithm supplied here.

Isotope NONE requires index and charge zero. SELECTED identifies one nonnegative
isotope index (zero is monoisotopic) and positive charge magnitude. The producer
supplies the corresponding canonical physical bounds. No isotope shift, widening,
calibration or tolerance arithmetic is inferred. Tolerance remains request policy
even when bounds coincide: consumers must honor every field or reject unsupported
requests, never silently ignore one. Invalid schema, enum, nonfinite values,
reversed bounds, negative tolerance and ambiguous isotope state are rejected.
There is deliberately no production hasher until a real WP-05 consumer exists.

## Coordinates and identities

`RecordContext` carries schema version (only 0 accepted), run, acquisition lane,
and calibration snapshot IDs. `RecordId` is an opaque caller-assigned uint64,
including zero. Producers must use a consistent namespace across pages/tiles:
run IDs identify runs in the caller's dataset, lane IDs within a run, frame IDs
within a lane, scan IDs within a frame. Calibration IDs identify the caller's
immutable calibration snapshot. Family, precursor and fragment IDs identify
entities within the caller's library/provenance namespace, not vector indices.
The caller must preserve that namespace when combining libraries. These records
neither generate IDs nor assert that different provenance namespaces are equal.

`Coordinates` preserves exact supplied double values: RT in **seconds**, m/z in
Th, frame/scan identities, exact source `native_id` string, MS level (nonzero), and isolation lower/upper absolute
m/z bounds (ordered, equal bounds allowed). No binning, calibration, aggregation,
or tolerance is applied. Non-isolated scans use equal isolation bounds supplied
by the producer (typically zero), not an invented precursor window.
Continuous mobility is optional and explicitly either milliseconds or inverse
reduced mobility in V s / cm². FAIMS is a separate optional **discrete voltage in
volts**, never a continuous mobility axis. Absent mobility means no axis, not
zero drift time. FAIMS 0 V is a valid present stratum. These are measurement
coordinates, not search indices or extraction query equivalence keys.

`RawPoint::intensity` is a finite supplied detector/baseline-corrected value in
producer-defined intensity units. Negative intensities and finite negative RT
are permitted. `Hypothesis` keeps caller family/precursor IDs, target/decoy,
precursor m/z, library provenance string and optional prediction provenance
string. Provenance strings are opaque source/revision or model identifiers,
not paths to be opened by this contract. Optional RT predictions use seconds,
optional mobility predictions use the same explicit units as coordinates.
Normalized/iRT values must not be put into the seconds field without conversion
to the run's reference. Unit-unspecified library IM must likewise remain absent
until the producer establishes units. Original library metadata stays with the
source/provenance store, this is not a lossless replacement for that store.
`Fragment` keeps a stable fragment ID, m/z and optional predicted intensity in
producer-defined units. No intensity normalization is inferred.

## Sparse shape and missingness

`HypothesisPage` uses flattened fragments with CSR offsets: exactly N+1 entries,
first 0, monotone, last equal to fragment count. Repeated offsets allow assays
with no fragments. `fragment(i,j)` checks both indices before offset addition.
There is no dense hypothesis-by-fragment-by-time tensor.

`EvidenceBlock` contains only supplied support rows with family/precursor/fragment
IDs and exact coordinates. `value` absent is an explicit missing observation,
0 is a measured zero, and an absent row makes no acquisition or coverage claim.
Mobility absence is independent of evidence missingness. Values use the producer's
intensity/evidence units. This v0 does not define an aggregation or score.

## Ownership and checks

`FamilyBatch` retains caller-declared `CandidateFamily` values. Each region has
exactly the batch context and uses the extraction-query units above. Family IDs
are unique within the batch. Member precursor IDs and alternative IDs are unique
within each family, not globally. They retain stable caller namespaces across
leases and must not be treated as positions. Zero and the full uint64 range are
valid IDs. Charge is an optional nonzero signed integer, not a magnitude or an
inferred charge. A present adduct is a nonempty opaque string. Target/decoy is
per member, never a derived family label. Provenance retains opaque source and
optional prediction identities.

An alternative's empty precursor list explicitly represents null. One precursor
is a singleton, more than one represents coexisting candidates. Several
alternatives preserve unresolved possibilities, not a selected winner. An empty
alternative vector declares no alternatives and is different from a null
alternative. Each alternative references only members of its own family, without
duplicates within that list. Sparse fragment incidence references family members,
with unique `(precursor_id, fragment_id)` pairs. The same fragment may support
several members. SHARED and DISCRIMINATING are caller declarations, not computed
roles. Fragment IDs use the caller's fragment namespace, no fragment catalogue is
inferred. UNKNOWN, PARTIAL and COMPLETE likewise record declared candidate
completeness, not scientific truth, confidence or exhaustive library coverage.

`ScoreBatch` owns family and evidence copy leases and requires equal complete
contexts. Every score names an existing family and one of its alternatives.
Duplicate `(family_id, alternative_id)` score rows are rejected. Evidence row
indices are relative to the retained evidence block. Referenced evidence must
have that family ID and a precursor belonging to that family. They need not be
restricted to the alternative's members: surrounding family evidence can inform
an alternative, including null. Empty evidence references are valid for every
alternative, including null. Repeated evidence indices are retained, not silently
deduplicated. Unreferenced evidence is not required to belong to this batch.

Identification, localization and quantity are independent optional observations.
Absent is not zero. Identification/localization are finite raw method-specific
scores, not probabilities and not bounded to [0,1]. Quantity is a finite value
with a nonempty explicit opaque unit, such as detector area. No normalization,
unit conversion or positivity is inferred. A nonempty `scoring_method` identifies
the caller's scoring procedure/revision. RESOLVED/UNRESOLVED is an explicit caller
ambiguity declaration, not a winner inference or calibrated confidence. These
records supply no scoring algorithm, wire format, cache or production adapter.

Family `records` counts families. Family `values` counts members + alternatives
+ all alternative precursor references + fragment incidence rows across the
batch. Score `records` counts score rows and score `values` counts flattened
evidence references only. Counts are checked by subtraction before accumulation
could overflow. Nested strings, query acquisition selections and retained input
leases are not a byte budget. Empty family/evidence/score batches are valid with
zero caps. Family and score accessors return detached values or owning leases,
including after their producers are destroyed. Malformed enums, duplicate identities,
dangling references, context mismatch, caps, bounds and moved-from access fail
with `InvalidValue` in production, including NDEBUG builds.

Constructors take vectors by value and move them into immutable shared storage.
Each concrete class is itself a cheap owning copy lease. All accessors return
values, including strings and optionals, never borrowed references/spans/pointers.
A copied lease survives producer destruction, assignment and moves. The final
lease releases storage. Access on a moved-from lease and out-of-range access
throws `OpenMS::Exception::InvalidValue`, as do malformed records. An empty valid
container is distinct from a moved-from lease.

The caller must supply `RecordLimits{records, values}` at construction. Raw and
evidence row counts must fit **both** caps. Hypothesis count fits `records` and
flattened fragment count fits `values`; offsets have the checked N+1 shape.
Zero caps permit empty containers. These are count limits, not byte accounting,
input decoding limits, or an admission scheduler. In particular, caller-owned
input allocations and provenance string bytes are not budgeted here.

Construction checks schema, finite coordinates/values, enum values, isolation
ordering, caps and offsets. It does not impose positivity, global uniqueness,
ordering of scans, scientific provenance validation or library completeness.
The family and score types additionally check the local references described above.
Consumers retain their own owning leases,
and adapters must establish scientific/unit compatibility before combining data.
Validation lives in `OpenDIA_impl`, shared by the tool and tests. The same test
body uses native ClassTest when available and active runtime checks against a
standalone SDK, including under NDEBUG.

## Real-loader boundary check

`OpenDIARecordsInterop_test` reads the existing workflow-1 non-IM, workflow-23 IM
and `IM_FAIMS_test.mzML` fixtures through the actual `MzMLFile` public interface.
The consumer checks every retained raw point against a separate load after the
producer experiment is destroyed, including exact source native IDs and all
coordinate fields. The mzML loader already converts source RT units to seconds.
An explicit test-owned identity registry maps workflow-17 transition source IDs
to caller IDs. `TransitionTSVFile` supplies the real library measurements and
target/decoy state. Its normalized RT and unit-unspecified IM are deliberately
not relabelled as raw seconds or known mobility units.

Query checks preserve selection, coordinates, mobility units and FAIMS presence
from these real acquisitions and distinguish a changed calibration snapshot.
A separate existing mzML fixture exercises the actual public
`ChromatogramExtractorAlgorithm`: the original request and a query-value
round trip must produce exactly equal time/intensity arrays with nonzero signal.
This checks the supported all-acquisition, all-RT, non-IM SUM request only, not
production cache integration or a universal query-to-extractor adapter.

The sparse evidence check transfers selected real observations and retains them
after its raw tile and hypothesis page are destroyed. It does not perform or
validate scientific matching. The fixture's full-file load, one-window spectra,
numeric scan mapping, count limits and by-value layout are not the future
WP-03/04 production adapters, a packed representation, byte-bounded streaming,
or a performance claim. Family and score handoff checks are test-owned fixture
adapters, not production adapters.

The complete handoff retains declared singleton, competing, coexisting and null
alternatives, partial completeness and unresolved scores after raw/page/local
producer destruction, including explicit missing and zero support in each
acquisition mode. These declarations do not assert scientific matching or null
calibration. Two native TOPP checks also consume genuine `results.tsv` exports
from the existing SQLite and Parquet CLI workflows, using the real `CsvFile`
interface. They preserve run/feature/precursor IDs, target/decoy, charge, RT/mz,
`d_score` and integrated intensity against an independent re-read. Localization
is left absent, not substituted with an identification confidence column.
Export rows do not establish a complete candidate population, so completeness
stays UNKNOWN. Unscored null alternatives are explicit test-owned additions,
not inferred scientific nulls. This test-only conversion is not a production
results importer or a lossless assay metadata store.
