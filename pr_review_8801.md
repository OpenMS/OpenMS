## Review of PR #8801: Refactor isEnz_ to use OpenMS ProteaseDigestion (#5571)

### Summary

This PR replaces the hardcoded enzyme cleavage logic in `PercolatorInfile::isEnz_()` with OpenMS's `ProteaseDigestion::isValidProduct()`, reusing enzyme definitions from `ProteaseDB`. It also adds two new enzyme entries (Thermolysin and Proteinase K) to `BuiltInProteaseDataProvider.cpp` and removes the duplicate `isEnz_()` from `PercolatorAdapter.cpp`.

The refactoring approach is sound, but **there are three enzyme mappings that produce different cleavage behavior** compared to the original Percolator-compatible hardcoded logic. These must be fixed before merging to avoid changing Percolator feature computation.

---

### Critical: Enzyme Specificity Mismatches

The core of this PR maps Percolator enzyme names to OpenMS `ProteaseDB` names. The problem is that several OpenMS enzyme definitions **do not match** the cleavage rules that Percolator uses. Below is a detailed comparison:

#### 1. `elastase` → `"leukocyte elastase"` — **MISMATCH**

| Property | Percolator (old code) | OpenMS `leukocyte elastase` |
|---|---|---|
| Regex | Cuts after {L, V, A, G}, not before P | `(?<=[ALIJVX])(?!P)` — cuts after {A, L, I, J, V} |
| **Missing** | **G** (Glycine) — Percolator includes it | Not in OpenMS regex |
| **Extra** | — | **I** (Isoleucine) — OpenMS includes it, Percolator does not |

**Impact**: Peptides with G at P1 will no longer be recognized as enzymatic; peptides with I at P1 will now incorrectly be recognized as enzymatic. Both change the `enzN`/`enzC`/`enzInt` Percolator features.

**Fix**: Either (a) create a new enzyme entry like `"Percolator elastase"` with regex `(?<=[LVAGX])(?!P)`, or (b) map to a different existing enzyme, or (c) add a note explaining the intentional behavior change.

#### 2. `pepsin` → `"PepsinA"` — **MISMATCH**

| Property | Percolator (old code) | OpenMS `PepsinA` |
|---|---|---|
| Logic | Cuts when P1={F,L,W,Y} OR P1'={F,L,W,Y}, AND P1≠R | `(?<=[FLJX])` — cuts after {F, L, J} only |
| **Missing** | **W, Y** at P1 — Percolator includes them | Not in OpenMS PepsinA |
| **Missing** | P1' check entirely — Percolator checks C-terminal residue too | OpenMS only checks P1 (N-terminal of cut) |
| **Missing** | **R exclusion** at P1 — Percolator blocks R | Not in OpenMS PepsinA |

**Impact**: This is a significant semantic mismatch. Percolator's pepsin model considers **both** sides of the cleavage site (P1 and P1'), while OpenMS's PepsinA only considers P1. Additionally, W and Y at P1 are lost. This will substantially change enzymatic feature computation for pepsin digests.

**Fix**: Create a `"Percolator PepsinA"` enzyme entry with a regex that matches the original dual-sided logic: `(?<=[FLWYX])|(?=[FLWYX])` with an R exclusion, or handle pepsin as a special case.

#### 3. `glu-c` → `"glutamyl endopeptidase"` — **MISMATCH**

| Property | Percolator (old code) | OpenMS `glutamyl endopeptidase` |
|---|---|---|
| Regex | Cuts after E, not before P | `(?<=[DBEZX])` — cuts after {D, B, E, Z} |
| **Extra** | — | **D** (and B=ambiguous D/N) — OpenMS includes Asp, Percolator does not |
| **Missing** | P exclusion at P1' | OpenMS has **no** P1' restriction |

**Impact**: Peptides with D at P1 will now be recognized as enzymatic when they shouldn't be (per Percolator's model). The missing proline restriction also changes behavior.

**Fix**: Map to an E-only enzyme entry with proline restriction, or create `"Percolator Glu-C"` with regex `(?<=[EX])(?!P)`.

---

### Enzymes That Match Correctly ✓

These mappings are verified correct:

| Percolator name | OpenMS name | Status |
|---|---|---|
| `trypsin` | `Trypsin` | ✓ `(?<=[KRX])(?!P)` matches `(K\|R) && !P` |
| `trypsinp` | `Trypsin/P` | ✓ `(?<=[KRX])` matches `K\|R` (no P restriction) |
| `chymotrypsin` | `Chymotrypsin` | ✓ `(?<=[FYWLJX])(?!P)` matches `(F\|W\|Y\|L) && !P` |
| `lys-n` | `Lys-N` | ✓ `(?=[KX])` matches `c==K` |
| `lys-c` | `Lys-C` | ✓ `(?<=[KX])(?!P)` matches `K && !P` |
| `arg-c` | `Arg-C` | ✓ `(?<=[RX])(?!P)` matches `R && !P` |
| `asp-n` | `Asp-N` | ✓ `(?=[DBX])` matches `c==D` |

Note: Arg-C's `(?!P)` is a simplification (biochemically Arg-C *does* cleave before Pro), but it's **consistent** between both systems, so the mapping is correct for this refactor.

---

### New Enzyme Entries (Thermolysin & Proteinase K) ✓

The new entries in `BuiltInProteaseDataProvider.cpp` correctly capture the Percolator logic:

- **Thermolysin**: `(?<![DE])(?=[AFILMV])|(?<=R)(?=G)` — correctly models "before {A,F,I,L,M,V} unless preceded by D/E, plus R-G bonds"
- **Proteinase K**: `(?<=[AEFILTVWY])` — correctly models "after {A,E,F,I,L,T,V,W,Y}"

Minor note on Thermolysin: the synonyms set is empty (`set<String>()`). Consider adding `"thermolysin"` as a synonym for lookup convenience, similar to how other enzymes have lowercase/alternate names.

---

### Code Quality Observations

1. **`thread_local` caching** — Good optimization for the digest object cache, avoids repeated `setEnzyme()` calls. However, note that `thread_local` variables persist for the thread's lifetime; if this code is used in a long-running multi-threaded application, consider whether the cache invalidation logic (`enz != cached_enz`) is sufficient.

2. **`isValidProduct()` with mini protein** — The approach of constructing a 2-character `mini_protein` and calling `isValidProduct(mini_protein, 0, 1, true)` is clever but relies on the assumption that `isValidProduct` correctly handles 2-char strings. Worth verifying with a unit test for edge cases.

3. **Removed duplicate in PercolatorAdapter.cpp** — Good cleanup. The duplicate was a maintenance hazard.

4. **Terminal handling** — The addition of `n == '[' || c == ']'` for terminal positions is a nice improvement over just checking `'-'`.

---

### Recommendations

1. **Must fix** (before merge): Resolve the three enzyme mismatches (elastase, pepsin, glu-c). The safest approach is to add Percolator-specific enzyme entries to `BuiltInProteaseDataProvider.cpp` that exactly replicate the old hardcoded logic.

2. **Should add**: Unit tests that verify `isEnz_()` produces identical results to the old implementation for all enzyme types. A simple table-driven test with representative amino acid pairs would catch regressions.

3. **Nice to have**: Add `"thermolysin"` as a synonym in the Thermolysin enzyme entry.

---

### Verdict

The refactoring direction is correct and desirable (resolving the TODO from #5571). However, the PR introduces **behavioral changes** in 3 out of 12 enzyme mappings. These need to be fixed before merging. Once the enzyme mismatches are resolved and test coverage is added, this should be ready to merge.
