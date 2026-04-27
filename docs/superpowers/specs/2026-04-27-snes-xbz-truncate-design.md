# SNES mother generation: truncate to unambiguous spans on X/B/Z

Tracks issue [#9192](https://github.com/OpenMS/OpenMS/issues/9192) item 2.

## Problem

`FragmentIndex::generateSNESMothers_` currently rejects an entire mother peptide
when its proposed window `[start, start+length)` overlaps any ambiguous residue
(X / B / Z). Valid shorter realizations anchored at the same terminus, ending
before (Single-N) or starting after (Single-C) the ambiguous residue, become
unsearchable.

Example: protein `ACDEFGHIXKLMNPQSTVWY`, `peptide_max_size=12`, Single-N anchor
at position 0. Today the length-12 mother is dropped because position 8 is `X`;
the still-valid length-8 prefix `ACDEFGHI` at the same anchor is lost. Symmetric
loss on the suffix side for Single-C.

## Approach

Approach C from brainstorming: keep the existing protein-level fast path for
the overwhelmingly common no-X/B/Z case; for proteins that contain at least one
ambiguous residue, split the protein into contiguous unambiguous spans up front
and run the existing Single-N / Single-C sweeps within each span.

`emitMother` is simplified: it loses its inner X/B/Z check (which is now
redundant — the fast path is gated at the protein level, the slow path is gated
structurally by span boundaries). It becomes a pure mass-compute + filter +
emit lambda.

## Architecture

`generateSNESMothers_` keeps two literal branches inside the per-protein OMP
loop body:

1. **Fast path — `!protein_has_ambiguous`**: existing sweeps over `[0, L)` run
   unchanged.
2. **Slow path — `protein_has_ambiguous`**: compute spans, then run the two
   sweeps per span with bounds `[span_start, span_end)` substituted for
   `[0, L)`.

The hot path (no ambiguous residues) sees zero change in instruction count —
the X/B/Z check moves from a per-mother cost to a single protein-level
`find_first_of` already present today.

## Span computation

```cpp
// Yields contiguous unambiguous spans of `seq`. Empty spans are not emitted.
auto spans = computeUnambiguousSpans(seq);  // local helper

// Equivalent to:
spans = []
p = 0
while (p < L) {
  bad = seq.find_first_of("XBZ", p)
  if (bad == npos) { spans.push_back([p, L));  break; }
  if (bad > p)     { spans.push_back([p, bad)); }
  p = bad + 1
}
```

The helper is a small file-local lambda or static function inside
`FragmentIndex.cpp`.

## Sweep semantics per span

For each span `[s, e)` with `e - s >= peptide_min_length_`:

- **Single-N**: `for (i = s; i + peptide_min_length_ <= e; ++i)`,
  `length = min(effective_max_length, e - i)`,
  `emitMother(i, length, /*single_c=*/false)`.
- **Single-C**: `for (j = s + snes_min_length - 1; j < e; ++j)`,
  `length = min(effective_max_length, j + 1 - s)`,
  `start = j + 1 - length`,
  `emitMother(start, length, /*single_c=*/true)`.

`effective_max_length` is computed per-protein (unchanged); the per-span min
with `e - i` / `j + 1 - s` enforces "stay inside the unambiguous span".

`snes_min_length = max(1, peptide_min_length_)` — preserves the current local
clamp guarding the `peptide_min_length_=0` corner case.

## Skipped-peptides counter

Today: `skipped_peptides.fetch_add(1)` is incremented once per mother whose
proposed window overlaps X/B/Z (typically thousands per ambiguous protein).

After: increment once per span that falls below `peptide_min_length_` (a much
smaller count). The log message ("skipped due to ambiguous residues or mass
filter") is informational; the count is used only for the user-facing summary
line. No callers depend on its magnitude.

## Edge cases

| Case | Behavior |
|------|----------|
| Whole protein is X/B/Z | Span list empty → zero mothers (matches today). |
| Protein shorter than `peptide_min_length_` | Existing `if (L < peptide_min_length_) continue` covers it. |
| Span shorter than `peptide_min_length_` | Skipped silently (counter incremented per skip). |
| `peptide_max_length_ == 0` ("no max") | `effective_max_length = L` already; per-span cap `e - s` (or `j + 1 - s`) holds. |
| `peptide_min_length_ == 0` | `snes_min_length = max(1, 0)` clamp inside the Single-C loop preserves the corner-case guard. |
| Multiple X/B/Z in protein | Each gap between consecutive ambiguous residues is its own span. |
| Adjacent X/B/Z (e.g. `...AXXB...`) | Empty interior spans are skipped. |

## Test coverage

Two new tests in `src/tests/class_tests/openms/source/FragmentIndex_test.cpp`,
beside the existing `(SNES mother generation rejects ambiguous residue spans
(X/B/Z))` block (which remains unchanged and still passes):

Both new tests share the parameter block from test 2419 (`peptide:enzyme_specificity=none`,
`peptide:min_mass=0`, `peptide:max_mass=50000`, empty fixed/variable mods,
`snes_enabled=true`) but with `peptide:min_size=8` and `peptide:max_size=12`
(currently both =8 in 2419). Mother-classification uses
`FragmentIndex_test::isSingleCMother(mother.mod_bitmask_)` — `is_single_c` is
not a direct field on `Peptide`.

1. **Single-N short realization in unambiguous prefix is kept.** Protein
   `ACDEFGHIXKLMNPQSTVWY` (X at index 8). Today: the start-0 Single-N mother
   of proposed length 12 spans X → dropped, with no shorter alternative kept.
   After: a Single-N mother at `start=0, length=8` must exist. Assertion: scan
   `fi.getPeptides()`, find at least one mother with
   `!isSingleCMother(mother.mod_bitmask_)`, `sequence_.first == 0`,
   `sequence_.second == 8`.
2. **Single-C short realization in unambiguous suffix is kept.** Same protein
   and parameters. The end-anchored mother at `j=19` (last residue) currently
   produces span `[8, 20)` length 12 — covers X → dropped. After: a Single-C
   mother at `start=9, length=11` must exist (length capped at
   `effective_max_length=12` but truncated to span `[9, 20)`). Assertion
   mirror: `isSingleCMother(mother.mod_bitmask_)`, `sequence_.first == 9`,
   `sequence_.second == 11`.

No regression-only assertion needed; the existing `(start > 8 || end <= 8)`
invariant in the prior test still holds (no kept mother spans the X) and is
exercised by the existing test.

## Out of scope

- Issue #9192 item 1 (`PROTEIN_N_TERM` / `PROTEIN_C_TERM` fixed-mod specificity).
  Tracked separately; will be implemented on its own branch after this lands.
- Non-SNES paths: `generatePeptides` already digests with `digestUnmodified`
  and rejects X/B/Z at the resulting peptide level (FragmentIndex.cpp:909-918).
  No change there.

## Files touched

- `src/openms/source/ANALYSIS/ID/FragmentIndex.cpp` — `generateSNESMothers_`
  body, `emitMother` lambda.
- `src/tests/class_tests/openms/source/FragmentIndex_test.cpp` — two new
  `START_SECTION` blocks adjacent to test 2419.

No header changes (`FragmentIndex.h` untouched).
