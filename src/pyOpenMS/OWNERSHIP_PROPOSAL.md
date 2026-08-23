# Ownership policy for pyOpenMS — proposal, **revised after adversarial audit**

Companion to `CYTHON_TO_NANOBIND_SEMANTICS_AUDIT.md`; addresses the open design
question in [#9792](https://github.com/OpenMS/OpenMS/issues/9792) Part VI/VII.

> **Revision note.** v1 of this document proposed *self-validating handles*
> (`{owner, index, data_snap, size_snap}`, re-resolved on every access) as the
> **default** for object-level access, and claimed *"All memory-safety failures
> are closed."* A 14-agent adversarial audit against every part of #9792
> refuted that claim in four independent ways and refuted the v1
> recommendation outright. **The mechanism survives; the recommendation does
> not.** v1's claims are corrected in §2 rather than quietly deleted, because
> several of them are the kind of plausible-but-wrong that this issue keeps
> producing.

---

## 1. Revised verdict

**The handle is a sound construction inside a narrow envelope, and the envelope
is far smaller than v1 implied. It should not be the default.**

The envelope: a handle rooted at a Python-**owning**, index-addressable
container, accessed through receiver-position get/set, at depth 1,
single-threaded. Inside it, the handle does what v1 said and is strictly better
than today's retained raw pointer. Outside it, it fails structurally.

The number that decides this — **verified directly, not taken from the audit**:

| | count |
|---|---|
| `rv_policy::reference_internal` sites in `bindings/*.cpp` | **356** |
| `__getitem__` bindings total / of those `reference_internal` | 21 / **16** |
| reference-returned containers (`-> std::vector<T>&`, list caster aliases per element) | **105** |

**The index locator natively expresses ~16–18 of 356 sites (~5%).** The rest are
member sub-object accessors (`getInstrumentSettings`, `getSample`, `getOptions`
— owners with no `data()`/`size()`, so v1's validator body is ill-formed for
them), reference-returned containers, and self-returns. `rv_policy` is a
one-token annotation; a handle is a rewritten return expression. Migration is
not a policy sweep — it is hundreds of hand-written lambdas.

**Recommendation: always-copy at element level plus the bulk table view
(#9792 VII.5 step 5) is now the better default.** Handles demote to a named,
distinct-type, opt-in accessor at a handful of sites — *contingent on a
spectrum-copy benchmark that #9792 VII.1 still lists as a TODO and that nobody
has run.*

---

## 2. Confirmed defects in v1

Each was verified against source or a compiled demo by an agent tasked with
refuting it; the ones marked **[verified here]** I re-checked directly.

**D1 — "All memory-safety failures are closed" is false in four ways.**
1. *Argument position is unguarded.* v1 scoped validation to "every attribute
   get and set". Argument passing is neither. For a `nb::class_`-bound type,
   argument conversion is `type_caster_base<T>::from_python` → `nb_type_get`:
   a typeid check, a state check, then `*out = inst_ptr(inst)`. **No hook
   exists.** `algo.process(exp[0])` hands the stale address to C++ with the
   fingerprint never consulted.
2. *The guard dereferences its own owner before it can validate.* Reading
   `c.data()`/`c.size()` requires `nb::cast<Container&>(owner)` first. If the
   owner is itself a non-owning wrapper — which all 356 `reference_internal`
   accessors hand out — the validator performs the use-after-free it exists to
   prevent. `keep_alive` protects the *root*, not a sub-object's validity, and
   nanobind sets `destruct=false` for `reference`/`reference_internal`, so
   ordinary refcounting keeps the *wrapper* alive, not the C++ storage.
3. *Buffer exports stay unguarded* — the tier v1 itself calls "the entire
   performance win". `clearMetaDataArrays()` dangles every `get_data_view()` in
   an experiment while a nested handle catches it: R-4.5 satisfied for handles
   and violated for views **in the same call**.
4. *Retention (R-4.6) and callbacks (R-4.7) are never addressed.* Grep over v1
   returns zero hits for `R-4.6`, `R-4.7`, `V.4`, `V.5`, `L6`, `III.8`.

**D2 — L3 is earned by `swap` alone.** The fingerprint tests container *shape*,
not *identity*. These preserve `data()` and `size()` bit-identically, with
probability 1 and no allocator luck: `setSpectra`/`setChromatograms`
(`MSExperiment.cpp` `spectra_ = spectra;`), `__setitem__` (13 sites),
`operator=`, `assign()`, `clear()`+refill (`[vector.capacity]` guarantees
`clear()` cannot change capacity), `resize`-down-then-up, `erase()+insert()`,
and all six `set_peaks` overloads. v1's §4 called this "one ABA case … allocator
reuses the block" — wrong twice: it is a family, and it is guaranteed by the
standard. Reallocating paths are not reliably caught either: two same-size
`select()` sorts return the original block from a LIFO free list.

**D3 — the single-type variant v1 recommended cannot be built.**
`nb::class_<T>` and `type_caster<T>` are mutually exclusive by `static_assert`.
Keeping `class_` means storing the element pointer in `inst->value` — i.e.
retaining a raw pointer, collapsing v1's own thesis — and registers the handle
at the element address in `inst_c2p`, making it a **participant** in #9792 III.8
rather than a bystander. The escape hatch `nb::implicitly_convertible<Handle,T>`
is worse: it builds a **temporary** destroyed at call return and hands it out as
a non-const lvalue, so every mutating algorithm writes into a doomed copy,
silently, and the handle then reports "valid".
**Rule to adopt: an implicit conversion may never target a type passed by
non-const reference.**

**D4 — "there is no check-then-act" is false, and the window is open today.**
The guard reads shared state, decides, then acts on that same state. That is
exactly R-9.1's shape. It does not need free-threading: **62 bindings release
the GIL** and ~11 of them structurally mutate a target container while released
(`FileHandler.loadExperiment` → `map.reset()`; Arrow feature import →
`push_back`). What v1 eliminated is a lease *table*, not the hazard.

**D5 — v1 violates the placement rule #9792 already decided.** R-4.1 and R-9.2
both require the counter to live *"in a stable heap control block — **not on the
Python wrapper**, which may be destroyed independently."* v1's struct lives
entirely on the wrapper.

**D6 — smaller corrections.** `keep_alive<0,1>` becomes unnecessary at **6**
sites, not 13 — 7 are outside the tiering, and `Param` cannot be index-converted
at all (its iterator is a forward tree walk; its `__getitem__` takes a string
key). The error type would be `ReferenceError`, not the `RuntimeError` #9792
VI.5 specifies for mutation-during-iteration — and nanobind has no builtin
`ReferenceError`. Chain depth is **6**, not "≤ 3", and the 5-hop idiom ships in
`pyopenms/addons/consensusmap.py`. v1's own ~455 argument-position figure did
not actually subtract `const` as its footnote claimed.

**D7 — the guard-signal inversion (found by the critic pass, not the audit).**
`MSExperiment::sortSpectra` is an unconditional in-place `std::sort`, bound with
the GIL released. The outer `data()`/`size()` are unchanged, so **spectrum
handles stay silent while writing to the wrong spectrum** — but `std::sort`
moves `MSSpectrum` objects, swapping their peak-buffer pointers, so **every
peak-level handle raises** even though its peak is alive, unmoved and
uncorrupted. One call: the guard fires exactly where there is no damage and
stays quiet exactly where there is.

**D8 — the alarm is data-dependent, which #9792 V.2 calls the worst possible
property for a safety rule.** Whether a sort reallocates depends on whether
metadata arrays are present. So `sortByPosition()` raises `ReferenceError` on a
spectrum carrying a float array and is silent on one that does not — same call,
same corruption, opposite signal, and the caller cannot predict which. #9792
established this about *guard rules*; nobody had applied it to the fingerprint.

---

## 3. The finding that changes the design: OpenMS already ships the identity key

**[verified here]** `FeatureMap` and `ConsensusMap` — the exact two containers
where #9792 *measured* L4 — already inherit `UniqueIdIndexer<T>`
(`FeatureMap.h:76`, `ConsensusMap.h:65`), and
`UniqueIdIndexer::uniqueIdToIndex()` (`CONCEPT/UniqueIdIndexer.h:61`) is
**self-healing by construction**:

> consult the hash map → *if the element at that index does not actually carry
> this unique id*, rebuild via `updateUniqueIdToIndex()` and retry → else return
> `Size(-1)` (invalid).

Keying a handle on `getUniqueId()` instead of an index therefore closes **L4
exactly**, on the containers where it was measured, with **zero core change,
zero annotation, and no ABA** — because it tracks identity rather than shape.
`Feature`/`ConsensusFeature` carry the id via `RichPeak2D : UniqueIdInterface`;
guard `UniqueIdInterface::INVALID == 0` by falling back to a copy. Spectra have
the analogous key in `getNativeID()`.

This is also #9792 R-8.2's own prescription ("stable keys must be explicit
immutable representations, never mutable OpenMS objects") and supplies the
replacement for the `exp[0] is exp[0]` identity/hash regression a handle
otherwise causes. **It is strictly better than the fingerprint on the residue
the audit calls "the mechanism", cheaper than a version counter, and it went
unmentioned in three documents.**

---

## 4. Live bug found on the way (unrelated to any of the above)

**[verified here]** `installIonMobilityArray` (`bind_spectrum.cpp:74`):

```cpp
self.getFloatDataArrays()[self.getIMData().first] = std::move(fda);
```

`FloatDataArray` publicly derives from `std::vector<float>`
(`METADATA/DataArrays.h:22-24`), so this move-assignment **deallocates the
destination buffer**. Any ndarray previously handed out by
`get_drift_time_array_view()` or `FloatDataArray.get_data_view()` — both capture
`data()` at call time with `self_obj` only as a keep-alive — then reads freed
memory, with the nanobind ownership graph fully intact. Two lines to reproduce:

```python
v = s.get_drift_time_array_view()
s.set_drift_time_array(new_vals)      # v now reads freed memory
```

Also reachable through `set_peaks(..., ion_mobility=...)` on the **default**
`metadata="error"` policy, since the alignment check exempts the array being
replaced. R-6.2's ordering rule does not cover it — that rule moved the
`"clear"` drop *after* the peak write; this free happens in the install step
afterwards. **Fix: assign in place when sizes match.** This is the most
actionable item in the whole audit and is independent of the ownership decision.

---

## 5. Revised recommendation

### Tier 1 — do now, independent of the ownership decision
1. **Fix `installIonMobilityArray`** (§4). ~10 lines. Live UAF in `develop`.
2. **Peak `__iter__` → copy**, for all three peak-like types (`Peak1D`,
   `ChromatogramPeak`, `MobilityPeak1D`). Restores 3.5, removes the
   indexing-vs-iteration split, justified independently by #9792 VII.1.
   *Still endorsed; ship regardless of everything else.*
3. **`__hash__ = None` sweep.** nanobind never installs `tp_hash`, and `.def` is
   a post-creation `PyObject_SetAttr` which — unlike a class body — does not
   null `__hash__`. ~66 value-comparable types are identity-hashable today;
   R-8.1 is violated now, with no handle involved.
4. **Strip the ~36 inert `reference_internal` policies** (primitives,
   `std::string`, `DPosition<1>` — value casters ignore `rv_policy`). Zero
   behaviour change; makes every VI.11 grep mean what it says.
5. **Fix exception translation before any `ReferenceError` contract is
   designed.** There are zero `register_exception_translator` hits outside the
   consumer header, so every `OpenMS::Exception` maps to Python `RuntimeError` —
   the slot VI.5 assigns to "modified during iteration". `pyopenms.set_option`
   does not exist in the tree either, though #9792 VI.5 and any deprecation plan
   both spend it.

> **Correction to the audit:** its Tier-1 also asked for R-9.5 GIL-release
> documentation. That is **already merged** and ships as
> `src/pyOpenMS/THREAD_SAFETY.md` **[verified here]**. #9792's status table is
> stale.

### Tier 2 — the measurement that actually decides this
6. **Benchmark object-level copy.** #9792 VII.1 measured 16-byte peak copies and
   left *"Real benchmark is a TODO"*. The entire case for handles rests on an
   unmeasured number. A free control group already exists and nobody used it:
   `OnDiscMSExperiment.getSpectrum()` returns **by value**, i.e. option A is
   already shipping on the largest files users have.
7. **Build the bulk table view** (#9792 VII.5 step 5). Largest win, no storage
   change, no semantic commitment — and it makes aligned-metadata operations
   correct by construction (R-6.1, R-6.3), a correctness argument independent of
   the copy/view default. This should be the next substantial work.

### Tier 3 — only if Tier 2 shows spectrum copies are unaffordable
8. Add handles as a **named, distinct-type, opt-in** accessor
   (`exp.spectrum_ref(i)`), leaving `__getitem__`/`__iter__` as owned values —
   satisfying R-3.1/R-3.2 instead of overriding them. Required amendments: owner
   must be Python-owning; **key on `getUniqueId()` where available (§3)** rather
   than index; explicit `.value()` for argument position; full dunder protocol;
   a new VI.2 category and VI.11 step, since a factory returning a handle
   matches none of the existing greps and would be invisible to the audit
   #9792 exists to make mechanical.

### Do NOT attempt
- **The single-type variant** (D3) — unbuildable, and unsound if forced.
- **`nb::implicitly_convertible<Handle, T>`** — silently mutates a temporary.
- **A generic handle sweep over the 356 sites** — research project, not rollout.
- **Removing any `keep_alive<0,1>`** — 7 of 13 sites are outside any tiering.
- **Adding `FREE_THREADED`** while any handle exists.

---

## 6. Open decisions for humans

1. **What is pyOpenMS for?** If the answer is read-and-extract, always-copy plus
   the bulk table wins outright and Tier 3 is unnecessary. @jcharkow already
   argued this position. **It must be answered first, not last** — v1 sequenced
   a prototype ahead of the question.
2. **Is `isinstance(s, MSSpectrum)` breakable in a major release?** No shape
   preserves it soundly: a nanobind subtype must physically carry an
   `MSSpectrum` at offset 0, the facade route costs the addon registry and the
   stub pipeline, and the single-type variant is unsound.
3. **R-6.4: should pyOpenMS be stricter than C++?** Making structured arrays
   read-only or detached closes the largest unguarded surface and costs the
   zero-copy write path. This is #9792 VII.4's conflict and is not resolvable
   technically.
4. **What deprecation budget exists?** `tests/unittests/test_mutable_references.py`
   **[verified here]** asserts today's aliasing as *correct* — `exp[0].setRT(99)`
   must mutate `exp`, with the failure message *"__getitem__ returned a copy
   instead of a mutable reference"*. Six such tests. **Always-copy fails them;
   handles pass them but invert their intent.** That file is the concrete answer
   to "what would we be breaking", and no prior document mentioned it.
5. **Should pyOpenMS promise a validated reference at all?** A guard that fires
   on shape but never on identity or invariants trains users to read silence as
   safety — the exact criterion v1 used to disqualify generation counters
   ("converts unsafe into believed safe"). D7 and D8 make this concrete.

---

## 7. Alternatives worth evaluating before committing

- **Copy-on-write handle** — re-resolve and read through with zero copy; on
  first *write*, detach permanently. Reads (the dominant path) pay nothing;
  writes never reach the container, so the L4 residue degrades from *silent
  wrong write* to *stale read*, and it restores the 3.5 baseline exactly.
  #9792 VII.4 already points here ("pandas is itself abandoning
  `SettingWithCopyWarning` for copy-on-write") and nobody followed. **Strongest
  unexplored option.**
- **Explicit transaction** — `with exp.edit(i) as spec:`, invalidated on exit.
  Makes retention lexically impossible, and since the handle owns nothing there
  is no `inst_c2p` entry, sidestepping III.8 — which is what defeated
  invalidate-on-return in V.5.
- **Immutable snapshot + explicit write-back** — soundest and cheapest; pays one
  copy per edit. Note R-5.2 already bans the obvious "you forgot to write back"
  `__del__` warning, so it must be taught, not detected.

---

## Appendix — verification status

**Verified directly in this checkout:** the 356/21/16/105 counts; `THREAD_SAFETY.md`
already merged; `test_mutable_references.py` contents; `UniqueIdIndexer`
inheritance and self-healing lookup; `FloatDataArray : std::vector<float>` and
the `installIonMobilityArray` move-assign.

**From the audit, not independently re-verified here:** the six-way
classification of all 356 sites and the resulting 5.1% figure (direction
confirmed by the 16/356 `__getitem__` count, which is a lower bound); the
compiled-mock results for the dangling-owner certification and the `-O2` codegen
substitution; the 233/204 argument-position recount; chain-depth 6; the ~66
identity-hashable types. Cost estimates ("~245 resolver expressions", "~21
receiver-position mutators") are order-of-magnitude, not derived.

**Known stale in #9792 itself:** the status table omits the merged R-9.5
documentation; V.5's cited line numbers for the `python_error` flattening
predate #9804; a later PR shifted `installIonMobilityArray` by one line.
