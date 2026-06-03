# Making the chemistry singleton databases immutable

OpenMS keeps several process-wide singleton "databases" under
`src/openms/include/OpenMS/CHEMISTRY/`. Most are already effectively read-only
after construction (`DigestionEnzymeDB`/`ProteaseDB`/`RNaseDB`,
`RibonucleotideDB`, `MonosaccharideDB` — their `getInstance()` returns a
`const` pointer and they expose no state-mutating public API).

Three were mutable:

| DB | Mutating API | Reason it mutates | Mutation class |
|----|--------------|-------------------|----------------|
| `ResidueDB` | `getModifiedResidue()` (×3) | Lazily creates + caches a modified `Residue` on first request | **Pure memoization** (observationally read-only) |
| `ModificationsDB` | `addModification()` (×2), `addNewModification_()` | Registers user-defined / unknown mods discovered while parsing files | **Genuine runtime registration** |
| `ElementDB` | `addElement(..., replace_existing)` | Adds, or replaces *in place*, an element's isotope distribution | **Genuine runtime registration** |

The distinction drives everything: a pure memoization cache can be made
immutable cheaply (logical-const), whereas genuine runtime registration cannot
be removed without either dropping a feature or relocating the mutable state
out of the global singleton.

---

## 1. ResidueDB — done (logical-const)

`getModifiedResidue()` is a transparent memoization: for a given
`(residue, modification)` it always yields the same modified residue, callers
only ever receive a `const Residue*`, and the cache exists purely for pointer
stability + de-duplication (already guarded by the `ResidueDB` OpenMP critical
section).

Change implemented:

* `getInstance()` now returns `const ResidueDB*`.
* All three `getModifiedResidue()` overloads are `const`.
* The cache members (`residue_mod_names_`, `const_modified_residues_`) are
  `mutable`; the lazy-insert helper was narrowed to
  `addModifiedResidue_(Residue*) const` (modified residues only) and
  `addModifiedResidueNames_(...) const`.

Result: no public method can change the observable state of the singleton — the
class is immutable to callers — while memoization stays as an internally
synchronized implementation detail. The handful of call sites that stored the
result in a non-`const ResidueDB*` were updated (`AASequence`,
`StaticModification`, the pyOpenMS binding, and three class tests). No behavior
change.

---

## 2. ElementDB — reclassified: NOT a low-risk removal

The original triage assumed `addElement()` was effectively dead (only a unit
test and the pyOpenMS binding call it). On closer inspection it is
**load-bearing**:

* With `replace_existing = true` it overwrites an existing element's isotope
  distribution **in place, preserving the pointer identity** — `ElementDB_test`
  asserts exactly this: *"ptr addresses cannot change, otherwise we are in
  trouble since `EmpiricalFormula` uses those"*.
* It is exposed in pyOpenMS specifically so users can supply custom isotope
  abundances (e.g. labelling / quant workflows) and have every existing
  `EmpiricalFormula` referencing that element pick up the change.
* The mutation therefore happens **after** `getInstance()`, at runtime — the
  same shape as `ModificationsDB`, not a pure cache.

Dropping `addElement` would be a real regression with no replacement, so it is
deliberately left unchanged here and folded into the Option C work below
(build-from-provider + a relocated override mechanism), rather than removed.

---

## 3. ModificationsDB — target design (Option C: reference DB + session registry)

`addModification()` registers modifications discovered *from data* at runtime —
unknown mass mods (`N[12345.6]`, `.n[...]`, `.c[...]` in `ResidueModification`),
mzIdentML mods, NuXL params, cross-links — that must later be looked up by
name/id and need singleton-lifetime pointers. These cannot be precomputed, so a
true fix must move the mutable state out of the global singleton rather than
hide it.

### Shape

1. **Immutable reference DB (global).** The singleton holds only the static
   UniMod / PSI-MOD / XLMOD data, loaded once through the *already-existing*
   `ModificationDataProvider` injection (`initializeModificationsDB`,
   `InMemoryDataProvider`). After load it is frozen and `getInstance()` returns
   `const ModificationsDB*`. The `addModification`/`addNewModification_` mutators
   and the `Residue`/`AASequence`/`CrossLinksDB` friendships are removed from the
   global type.

2. **Mutable session registry (non-global, caller-owned).** Runtime-discovered
   mods live in a small `ModificationRegistry` owned by the parsing/analysis
   context (natural home: `IdentificationData`, or a context threaded through the
   file handlers). Lookups consult the local registry first, then fall back to
   the immutable global reference DB. Interning stays idempotent (same `fullId`
   ⇒ same pointer).

3. **Migration.** Introduce the registry + a lookup facade, then convert call
   sites incrementally:
   - `ResidueModification::create*FromMassString` and friends take a registry.
   - `MzIdentMLDOMHandler`, `NuXLParameterParsing`, `CrossLinksDB`,
     `AASequence`/`Residue` thread the registry through instead of mutating the
     global.
   - `ElementDB`'s custom-isotope override is given the same provider-based
     build path (and, if post-construction override must remain, a per-context
     element overlay analogous to the mod registry).

### Why

* True immutability + thread-safety (no hidden global growth behind OpenMP
  critical sections).
* Removes a long-standing source of non-determinism: today a modification's
  index in `mods_` depends on which files a process happened to parse first.

### Cost / risk

High and cross-cutting — it touches identification-data flow broadly, so it is
intentionally **not** bundled with the ResidueDB change. The pre-existing
`ModificationDataProvider` dependency-injection scaffolding indicates this is the
already-intended direction.

---

## Status

* [x] ResidueDB — immutable (logical-const), implemented.
* [ ] ElementDB — Option C (provider-based build + relocated isotope override).
* [ ] ModificationsDB — Option C (reference DB + session registry).
