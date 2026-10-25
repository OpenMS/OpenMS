# OpenDIA private records, unstable v0

These contracts implement roadmap step 2 only. They are tool-local, not installed,
not a public OpenMS API, and **not frozen**. Scientific freeze is pending the WP-07
experiment review. Extraction keys, family inference, score semantics and full
WP-01 completion are deliberately not implemented here.

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
ordering of scans, scientific provenance validation, library completeness or
cross-record referential integrity. Consumers retain their own owning leases,
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

The sparse evidence check transfers selected real observations and retains them
after its raw tile and hypothesis page are destroyed. It does not perform or
validate scientific matching. The fixture's full-file load, one-window spectra,
numeric scan mapping, count limits and by-value layout are not the future
WP-03/04 production adapters, a packed representation, byte-bounded streaming,
or a performance claim. Keys, alternatives and scores remain later WP-01 steps.
