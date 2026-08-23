# Adversarial audit of the ownership proposal — full findings

Raw output of the 14-agent audit described in `OWNERSHIP_PROPOSAL.md`. 62 gaps survived
adversarial verification, 9 were refuted. Retained in full because several findings are
actionable independently of the ownership decision (see the proposal's Tier 1).

**Not all claims here were independently re-verified — see the proposal's Appendix for
which were.** Line numbers drift.

---

## 1. VERDICT

**Survives as a mechanism; fails as the ownership policy. It needs amendments large enough that the recommendation itself should change.**

The self-validating handle is a real and correct construction inside a narrow envelope: a handle rooted at a Python-*owning*, index-addressable container, accessed through receiver-position attribute get/set, at depth 1, single-threaded. Inside that envelope it does what §2 says, it is a strict improvement over today's retained raw pointer (it never converts a loud failure into a quiet one), and §3's central insight — read actual state instead of trusting annotations — is sound. Outside that envelope it fails in five structurally different ways, and the envelope is much smaller than the document implies: **18 of the 356 `rv_policy::reference_internal` sites are natively expressible by the `{owner, index, data_snap, size_snap}` struct (5.1%)**, argument-position conversion cannot be interposed at all for a `nb::class_`-bound type, the guard must dereference its owner before it can validate anything, and the fingerprint is a predicate over container *shape*, not container *identity* — so `exp.setSpectra(same_length_list)`, `exp[i] = spec`, `clear()+refill`, `assign()`, `resize`-down-then-up, `erase()+insert()` and `set_peaks()` all pass validation deterministically, with no allocator luck, while every element is different. The decisive finding is that §3 refutes §6.1 by its own argument: whatever makes a generation counter unannotatable at the argument-position sites makes a handle unguardable at those *same* sites, and §6.1 never notices. Combined with §6.2's admission that the tier holding "the entire performance win" stays unguarded, the sentence at `OWNERSHIP_PROPOSAL.md:63` — *"All memory-safety failures are closed"* — is not defensible in any reading. The honest position after this audit is that **always-copy plus the bulk table (VII.5 step 5) now looks better than handles as a default**, and handles should be demoted to a named, distinct-type, opt-in accessor at a handful of sites, contingent on a spectrum-copy benchmark that #9792 VII.1 still lists as a TODO and that nobody has run.

---

## 2. THE LOAD-BEARING GAPS

Ranked by how much each moves the decision.

### G1. Argument position is outside the mechanism — and §6.1's recommendation forecloses the only fix
**What breaks.** §2 scopes validation to *"Every attribute get **and** set"* (`:43`). Argument passing is neither. For a type bound with `nb::class_`, argument conversion is `type_caster_base<T>::from_python` → `nb_type_get`, whose entire success path is a typeid/subtype check, a 2-bit `inst->state` check, then `*out = inst_ptr(inst)` (nanobind 2.12.0 `nb_cast.h:464-469`, `nb_type.cpp:1478-1541`, `nb_internals.h:519-522`). No virtual, no per-instance hook, no user callback. `algo.process(exp[0])` hands the stale address straight to C++ with the fingerprint never consulted. Measured surface for the five target containers: **233 mutable-`T&` argument occurrences on 224 registration lines, plus 204 non-self `const T&`** — argument position is the *majority* of the call surface for four of the five (FeatureMap 40 receiver / 85 argument).

**Why it is fatal to the recommendation.** §6.1 recommends *"one type with two internal states"* (`:172-175`). The only inbound hook nanobind offers is a `type_caster<T>` specialisation, and `nb_class.h:549-557` static_asserts that `nb::class_<T>` and `type_caster<T>` cannot coexist. So the recommended variant must keep `class_<MSSpectrum>` and store the element pointer in `inst->value` — which retains a raw C++ pointer, collapsing §1's *"the axis that actually matters is retention"*, and additionally registers the handle at the element address in `inst_c2p`, making it a first-class participant in III.8 rather than a bystander. The escape hatch (`nb::implicitly_convertible<Handle, T>`) is worse: `nb_type_get_implicit` builds a **temporary** and appends it to the cleanup list (`nb_type.cpp:1444-1452`), `cleanup_list::release()` destroys it at call return (`nb_lib.h:59-60`), and `operator Type&()` hands it out as a non-const lvalue (`cast.h:501-504`) — so every mutating algorithm writes into a doomed copy, silently, and the handle's own fingerprint still reports "valid" afterwards.

**Smallest amendment.** Delete the §6.1 recommendation. Use a *distinct* handle type (which R-3.3 already prescribes) with an explicit `.value()`/`.copy()`; state as a rule that a handle may never be implicitly convertible to `T&` at argument position, because that conversion is unguardable. Add the corollary to §2: validation covers receiver position only.

### G2. The fingerprint tests shape, not identity — and it is defeated deterministically by ordinary bound APIs
**What breaks.** `(data(), size())` is preserved, bit-identically and with probability 1, by:

| operation | in-tree route | evidence |
|---|---|---|
| `setSpectra` / `setChromatograms` | `spectra_ = spectra;` — `MSExperiment.cpp:1087-1090`, bound twice at `bind_experiment.cpp:105-106` | vector copy-assign reuses the buffer at or below capacity (compiled) |
| `__setitem__` | `self[i] = val;` — `bind_experiment.cpp:147-150`, 13 sites tree-wide | `MSSpectrum::operator=` copy-assigns the peak vector (`MSSpectrum.cpp:503`) |
| `clear()` + refill | `spectra_.clear(); chromatograms_.clear();` — `MSExperiment.cpp:1276-1286`, no shrink; same for `FeatureMap.cpp:466-478`, `ConsensusMap.cpp:255-267`, `ExposedVector.h:199-203` | `[vector.capacity]`: `clear()` cannot change capacity |
| `assign()`, `resize(smaller); resize(larger)`, `erase()+insert()` | element-level | compiled, four variants |
| `set_peaks` (all 6 overloads) | `checkPeakMetadataAlignment(...); self.resize(n); fill(self, n);` — `binding_utils.h:351-358` | `resize(n)` with `n<=capacity()` never reallocates |
| in-place sorts | `std::sort(this->begin(), this->end(), Feature::RTLess());` — `FeatureMap.cpp:238-241`, unconditional | no allocation at all |

Two consequences. First, **the L3 row (`:59`) is earned only by `swap`.** #9792's L3 section names four members of that class in one sentence — *"`setSpectra`, `operator=` and `__setitem__` … belong to the same class"* — and three of them leave `data()` unchanged. Second, **§4's characterisation of the residue is wrong twice**: it is not *"one ABA case … allocator reuses the block"* (`:122`) but a family, and the reuse is guaranteed by the standard rather than granted by the allocator. Reallocating operations are not reliably caught either: `MSSpectrum::select()` allocates-fills-swaps-frees (`MSSpectrum.cpp:21-33`, and the same pattern on each data array at `:48-55`, `:69-76`), so two same-size passes return the original block from a LIFO free list and the fingerprint matches again.

**Smallest amendment.** Restate the residue honestly and move `setSpectra`/`setChromatograms`/`__setitem__`/`operator=`/`clear`/`resize`/`erase` out of the L3 row into it. Then note that the cheap sound close is *not* the ~455-site counter §3 argues against: these are all receiver-position bindings, roughly 15 of them plus ~6 clear/reset sites, enumerable by grep. A `uint64 version` bumped by exactly those, plus a `version_snap` in the struct, converts the worst silent-wrong-write case into a `ReferenceError` and has no ABA. That is not a "strict, additive refinement" (`:126`) — it is the difference between a shape check and an identity check, i.e. the mechanism.

### G3. The guard dereferences its own owner before it can validate anything
**What breaks.** The validator's first line is `nb::cast<Container&>(owner)`, and `nb_type_get` has no per-instance C++-liveness state to consult — `inst->state` is nanobind's own construction bookkeeping and stays `state_ready` forever on a `reference_internal` wrapper whose referent has been destroyed. So the guard must dereference in order to read `data()`/`size()`; it can validate the container's *contents*, never the container *pointer*. `keep_alive` (installed at `nb_type.cpp:1736`) protects the **root**, not the validity of a sub-object inside it. Reachable in four lines with no threads and no nesting:

```python
c = exp.getChromatogram(0)   # bind_experiment.cpp:115, reference_internal -> NON-OWNING wrapper
h = c[0]                     # handle: owner = c, fingerprint over c.data()/c.size()
exp.reset()                  # MSExperiment.cpp:858-864 destroys chromatograms_' elements
h.getIntensity()             # cast<MSChromatogram&>(c) -> dangling; read c.data() -> UAF INSIDE THE VALIDATOR
```

A compiled mock with MSExperiment's real member order (`offsetof(spectra_) = 184`) rooting a handle on a non-owning spectrum wrapper and then calling the equivalent of `exp.clear(False)` prints `match=YES` and reads through the freed peak block — because `~vector` deallocates but does not null the member bytes, so the fingerprint matches bit-identically while `exp` is alive and every refcount and keep-alive edge is satisfied. (Honest limit: when the *outer* vector reallocates, `std::vector`'s move constructor must null the source's `_M_start` and glibc then writes tcache bookkeeping over it, so that variant raises rather than certifying — verified in 4/4 configurations. The certification failure is specific to destruction-without-move, which `clear()`/`reset()` supply.)

**Smallest amendment.** Make it a stated precondition of §2 that a handle's owner must be a Python-owning instance — checkable as `nb::inst_state(owner).second == true`, the `destruct` bit nanobind sets false for `reference`/`reference_internal` (`nb_type.cpp:1732`) — with a fallback to a copy otherwise. That precondition excludes every one of the 356 `reference_internal` accessors as a handle root, which forces `getSpectrum`, `getChromatogram`, `getSpectra`, `getMSData` etc. to be migrated *in the same step* as `__getitem__`, not deferred to §7 step 4.

### G4. Coverage: the struct is a locator for one of six shapes
Mechanical classification of all 356 sites: **127 member sub-object accessors** (`FileHandler::getOptions`, `HPLC::getGradient`, `Instrument::getSoftware`, `MSSpectrum::getInstrumentSettings` — owners with no `data()`/`size()` at all, so §2's body is ill-formed for them), **93 containers-of-objects**, **40 self-returns**, **35 inert**, **25 keyed-or-index**, **12 iterators**, plus 24 with no explicit `-> T&`. The 93 container getters are the hard structural case: nanobind's list caster calls `Caster::from_cpp(value, policy, cleanup)` per element (`stl/detail/nb_list.h:60`) and the only owner information reachable is `cleanup_list::self()`, a bare `PyObject*` — no index, no member identity, and for `spec.getPrecursors()` the container that would have to be fingerprinted is `precursors_`, which the caster never sees. **`rv_policy` is a one-token annotation; a handle is a rewritten return expression.** Migration is 93 hand-written lambdas (or ~56 registered sequence types), not a policy sweep — and VII.2's stated case for option C over B was *"no per-element view class explosion"*, which a lazy container-getter class reintroduces one level up.

Related, and unpriced: re-resolution requires re-invoking the original accessor (`AASequence::getResidue` returns `*peptide_[index]` over `std::vector<const Residue*> peptide_`, `AASequence.cpp:130-136` / `AASequence.h:622` — the referent is not at `data() + i*sizeof(Residue)`), so the full sweep needs ~245 resolver expressions, against a §7 that estimates the work in units of 4 containers.

**Smallest amendment.** Replace `size_t index` with a tagged locator (`WholeObject | Index | Key | Pixel | Member`), re-derive §5's tiering over the six classes rather than three type names, and state explicitly which classes stay as today's unguarded aliases (VII.3's "then it must be documented as unsafe, not as guarded").

### G5. The headline sentence is contradicted by the document's own §6.2, and by 62 GIL-release sites
`:63` is unqualified. `:141` and `:176-178` keep `_view`/`_struct`/`get_data_view` zero-copy and unleased, and every one of those bindings captures a raw pointer at call time with `self_obj` only as a keep-alive, never re-deriving it (`bind_spectrum.cpp:300-317`, `:319-341`, `:450-459`; `bind_kernel.cpp:1378-1385`). `MSExperiment::clearMetaDataArrays()` — bound at `bind_experiment.cpp:100`, implemented at `MSExperiment.cpp:871-889` as `clear()` + `shrink_to_fit()` on all three array vectors of every spectrum — dangles every `get_data_view` in the experiment while a nested handle catches it. That is R-4.5 satisfied for handles and violated for views **in the same call**, on the tier §5 says holds "the entire performance win".

Separately: 62 bindings release the GIL (verified: `arrow_zerocopy.cpp` 16, `bind_experiment.cpp` 2, `bind_format.cpp` 13, `bind_misc.cpp` 31), and ~11 of them structurally mutate one of the four target containers with the GIL released — `FileHandler.loadExperiment` (`bind_format.cpp:477-484` → `MzMLFile.cpp:133` `map.reset()`), `featuremap_import_features_from_arrow` (`arrow_zerocopy.cpp:335-342` → `FeatureMapArrowIO.cpp:1471` `push_back`), `Biosaur2Algorithm.run` (`bind_misc.cpp:401` → `Biosaur2Algorithm.cpp:199` `feature_map.clear(true)`). The GIL serialises only threads that hold it. §2's bullet at `:69-72` — *"there is no check-then-act"* — is simply false: the guard reads shared state, decides, then acts on that same shared state; what was eliminated is a lease *table*, not the hazard R-9.1 names. (#9792's L5 scope decision does retire general thread-safety, so this is an unqualified-claim defect rather than an unmet obligation — but L5 also says *"the requirement that anything pyOpenMS invents be sound"* stands.)

**Smallest amendment.** Qualify the table rows as "single-threaded; element/object path only", and add the row §2 is missing: buffer export — **no**.

### G6. §6.4 is not expressible by §2's struct, and its depth bound is wrong
The struct holds one `owner`, one `index`, one snapshot pair. There is no chain to walk. Meanwhile the alias graph over all 356 sites admits an 8-hop path; discounting the two terminal hops that land in process-lifetime `ResidueDB` memory (Category S), the genuine Category-H depth is **6**, and 5 from `FeatureMap`. The 5-hop idiom is not hypothetical — `pyopenms/addons/consensusmap.py:98-107` and `:333-341` do `for f in cmap: f.getPeptideIdentifications()[0].getHits()[0].getSequence()` inside the shipped `to_df` path. Chain validation is also unsound in one direction and over-strict in the other: `exp[i] = other` with an equal peak count passes at *both* levels while every element differs (compiled: `outer same=1 | inner same=1 | inner[0].mz=99.0 (was 1.0)`), whereas an outer reallocation invalidates the chain even though #9792's own L2 measurement shows the inner peak buffer did not move.

### G7. Engineering surface never priced
- **Protocol.** One member is specified (`.copy()`, `:180`). A handle needs `__hash__ = None`, symmetric `__eq__`/`__ne__` across handle and value (nanobind's `nb::self == nb::self` is same-type-only, so `handle == owned_copy` falls back to identity `False`), a non-raising `__repr__` (or pdb/pytest rendering blows up), `__copy__`/`__deepcopy__` redefined (all 434 existing pairs return `T(self)`; an inherited `__deepcopy__` would deep-copy the whole experiment), no `__init__`, and a decision on `__bool__` (0 sites; truthiness falls to `__len__` at `bind_spectrum.cpp:286`, so `if spec:` on a stale handle now raises where today it cannot).
- **Identity regression.** `nb_type_put` returns the *existing* wrapper from `inst_c2p` for any non-copy rvp (`nb_type.cpp:1790-1806`), so `exp[0] is exp[0]` is `True` today and dict/set membership works off it. A handle has no resolved address to intern by; both silently become `False`, and R-8.1 forbids repairing it with a value hash.
- **Downstream machinery.** `pyopenms/addons/__init__.py:88-92` injects Python methods **by class name**, so a distinct type silently loses 11 MSSpectrum / 12 ConsensusMap / 8 MSExperiment / 7 FeatureMap / 6 MSChromatogram methods with no diagnostic. `fix_stubs.py:173-210` dedupes only *consecutive identical* overload blocks, so const/non-const pairs that gain differing return types stop deduplicating. Neither is mentioned.
- **VII.4's actual question — "where does the layer live?" — is never answered**, and §7 scopes the rollout to ~6 of 356 aliasing sites with no generator, no R-10.1 docstrings and no R-10.2 CI rule.

---

## 3. WHERE THE PROPOSAL OVERCLAIMED

Ordered by document position. Original quoted; correction follows.

**`:63-64`** — *"**All memory-safety failures are closed. The residue is purely semantic**"*
> False in four independent ways: (a) argument position is unguarded at ~233 mutable + 204 const sites (G1); (b) the validator dereferences an unvalidated owner whenever the root is a non-owning wrapper (G3); (c) the zero-copy view path the same document keeps at `:141`/`:176-178` is untouched, and `clearMetaDataArrays` dangles every `get_data_view` (G5); (d) R-4.6 retention-past-the-call is never addressed. Correct statement: *"All memory-safety failures in #9792's L1–L3 taxonomy are closed for receiver-position access to a handle rooted at a Python-owning, index-addressable container, single-threaded. Buffer exports, argument position, retained pointers and non-owning roots are not covered."*

**`:57`** — L1 row, *"strong Python ref; ordinary refcounting"*
> True only when `owner` is a Python-owning instance. `nb_type.cpp:1732` sets `destruct = false` for `reference`/`reference_internal`; ordinary refcounting keeps the *wrapper* alive, not the C++ storage. 356 bindings hand out such wrappers.

**`:59`** — L3 row, *"`data()` changed (measured in #9792 L3: `0x…470 → 0x…a70`)"*
> Earned by `swap` alone. `MSExperiment.cpp:1089` is `spectra_ = spectra;` — verified — and `setSpectra`, `operator=` and `__setitem__`, which #9792's L3 sentence names alongside `swap`, all leave `data()` and `size()` bit-identical at or below capacity.

**`:43`** — *"Every attribute get **and** set performs"*
> Correct and complete as written, which is the problem: attribute get/set is receiver position, and §3's own table (`:83-94`) exists to enumerate the argument-position sites the mechanism then never returns to.

**`:47`** — `if (c.data() != data_snap || ...)`
> `MSExperiment`, `FeatureMap` and `ConsensusMap` expose no `data()`. Routes exist (`getSpectra().data()`, `MSExperiment.h:1214`; `getData().data()`, `ExposedVector.h:333`) so §3's zero-core-change conclusion survives, but the FeatureMap/ConsensusMap route is never named.

**`:69-72`** — *"There is no lease table … so #9792's R-9.1 check-then-act race cannot arise, because there is no check-then-act."*
> The guard *is* check-then-act: read `c.data()`/`c.size()`, decide, act on `c[index]`. What was eliminated is a lease table. And under §6.1's recommended variant the "act" is an entire `MSSpectrum` method (e.g. `sortByPosition` over 100k peaks), i.e. a far longer validated interval than §2's `getRT()` example implies, which no single CAS could fuse.

**`:104-106`** — *"an unannotated third-party mutator is caught for free"*
> Only if the mutator has fully **completed** before the reader's check, and only if it changed `data()` or `size()`. Eleven GIL-releasing mutators can be in flight; seven classes of completed mutator preserve the fingerprint exactly (G2).

**`:115-117`** — *"one type check plus two integer compares"*
> Depth 1 only. At the 5-hop idiom shipped in `addons/consensusmap.py`, per-attribute cost is 5 casts + 10 compares + 5 bounds checks before any C++ work, because each level must resolve its parent to compute its own fingerprint.

**`:121-123`** — *"one ABA case (`clear()` + refill to the same size, allocator reuses the block)"*
> Not one case and not allocator-dependent. `[vector.capacity]` guarantees `clear()` does not change capacity, so the match is probability 1 with no allocator involved. The family additionally contains `setSpectra`, `operator=`, `assign()`, `resize`-down-then-up, `erase()+insert()`, `__setitem__` and `set_peaks` — and double-`select()` ABA even on reallocating paths.

**`:126`** — *"a strict, additive refinement of the same handle, not a redesign"*
> A version word is what converts the fingerprint from a *shape* predicate into an *identity* predicate. Under the corrected residue it is the mechanism, not an increment on it. It is also the same core change that free-threading readiness would require, so deferring §4 defers that too.

**`:152-155`** — *"**Object-level: source-compatible for correct code.** … Code that was *silently broken* now raises."*
> Baseline is unreleased `develop`. Against release/3.5.0 — the only baseline users have — `CYTHON_TO_NANOBIND_SEMANTICS_AUDIT.md:16-24, :162-168` records that object access returned an independent copy, so `s = exp[i]; s.setRT(t)` as a deliberate scratch edit now mutates the experiment, silently, with nothing structurally detectable and therefore nothing to raise. Additionally, a `size_snap` mismatch means one `addSpectrum` invalidates *every* live handle into the container even when no alias was actually broken — a false-positive class the document never mentions.

**`:160-162`** — *"`keep_alive<0,1>` at 13 sites becomes unnecessary"*
> 6, not 13. Verified: exactly 13 `make_iterator` sites, each with the annotation; §5's tiering plausibly covers MSSpectrum, MSChromatogram, Mobilogram, MSExperiment, ConsensusMap, FeatureMap. The other 7 — DeconvolvedSpectrum, PeakGroup, AASequence, IsotopeDistribution, NASequence, Param, PeptideIdentificationList — are outside the tiering. `Param` cannot convert at all: `ParamIterator` is a forward tree walker (`Param.h:171`, `operator++` only), `Param` has no `operator[](size_t)`, and its Python `__getitem__` takes a string key. Also note that holding `nb::object owner` **is** a keep-alive — internalised, not eliminated.

**`:162-164`** — *"mutation during iteration produces a clean `RuntimeError`, which is … what #9792 VI.5 already asks for"*
> Wrong on type and on trigger. §2's mechanism raises `ReferenceError` (`:48`), and it fires only on a `data()`/`size()` delta — so `FeatureMap::sortByRT` (`FeatureMap.cpp:238-241`, unconditional in-place `std::sort`) and `MSSpectrum::sortByPosition`'s fast path (`MSSpectrum.cpp:411-415`) produce **no signal at all**, which is precisely #9792's L4 iteration demo (`before: [30.0, 10.0, 20.0]` → `after: [10.0, 20.0, 30.0]`). nanobind has no builtin ReferenceError either (`nb_func.cpp:157-172`), so it needs a raw `PyErr_SetString`.

**`:172-175`** — *"**Recommend the single-type variant** — R-3.3 was written for buffer view-vs-value, where confusing the two is a data-corruption risk; here both states are safe."*
> Self-refuting. §4 concedes eleven lines earlier that L4 and ABA are *"memory-safe and semantically wrong"*, and #9792's L4 section measures exactly that as data corruption (`after f0.setMZ(999) : [(10.0, 999.0), …]`, RT=30 feature untouched). Both states are **not** safe, so R-3.3's rationale transfers exactly. Independently, the variant is the one shape that forecloses the inbound hook (G1) and makes the handle a first-class participant in III.8. Also unaddressed: R-3.1/R-3.2 have a *legibility* purpose independent of safety — R-3.5's representation extension and Status lesson 4 (from V.6/#9809, which had no aliasing in it) show the rulebook already treats naming as load-bearing on its own.

**`:181-183`** — *"must validate the whole chain, depth ≤ 3"*
> Measured Category-H depth is 6 (8 binding hops, 2 of them Category S), and the struct at `:34-41` cannot express a chain of any depth.

**`:83-94, :96-97`** — the ~455 argument-position table, footnote *"minus `const`"*
> `const` was not subtracted. The published column reproduces exactly as (all `OpenMS::T&`) − (`OpenMS::T& self`), which also counts `-> T&` return positions, `nb::init<>` and internal `nb::cast<>`. True mutable-argument count for the five containers is **233 on 224 lines**, with 204 const-argument occurrences alongside. The table is also *understated in coverage*: `PeptideIdentificationList` — a `reference_internal __getitem__` container (`bind_metadata.cpp:532-535`) — contributes 88 more mutable-argument lines than any listed container and is absent. §3's conclusion survives the correction; its arithmetic does not.

**`:32` and `:50`** — *"re-derive the address on every access"* / `// freshly resolved address`
> At `-O2` both GCC 13.3 and Clang 18 prove `c.data() == data_snap` on the taken branch and substitute the snapshot: `cmpq (%rdi),%rsi / jne / movq 8(%rdi),%rax / subq %rsi,%rax / … / movsd (%rsi,%rcx),%xmm0`. Single-threaded this is value-identical and harmless; the wording is nonetheless inaccurate, and it is fixable — an out-of-line call between check and act makes GCC emit a genuine reload.

**Absent entirely** (grep over the whole file returns zero hits): `R-4.6`, `R-4.7`, `V.4`, `V.5`, `callback`, `consumer`, `caster`, `#9803`, `#9804`, `L6`, `invariant`, `updateRanges`, `sorted`, `range`, `R-7`, `III.8`. That is the entire inbound half of the rulebook plus the invariant-coherence level. Two consequences worth naming: (i) re-resolution fires only when Python touches the object, so a pointer retained by C++ (R-4.6) or a value handed into a callback (R-4.7) is by construction outside the mechanism — and #9804 already shipped an *owning move-in/move-back* contract for what §6.1 would make the same Python type, with the opposite post-condition (a retained callback object silently stops writing through; a handle raises). (ii) The proposal blesses `s = exp[0]; s.setRT(5)` at `:153` as *"keeps behaving identically"* — a write that provably cannot move `spectra_.data()`/`size()`, so the guard cannot fire, while breaking both L6 failures one level up (`MSExperiment::RTBegin`/`RTEnd` are `lower_bound` over `spectra_` with the documented *"Make sure the spectra are sorted … Otherwise the result is undefined"*, `MSExperiment.cpp:598-640`, backing `areaBegin`/`areaEnd`; and `updateRanges` is manual, `MSExperiment.h:1025`). The handle adds no new corruption there — `exp[0]` already aliases — but it trains the user that silence means safety, which is the exact criterion §4 used at `:100-102` to disqualify generation counters: *"converts 'unsafe' into 'believed safe'"*.

---

## 4. HAZARDS IN NEITHER THE PROPOSAL NOR #9792

1. **A live use-after-free in shipped `develop`, in the code that fixed R-6.2.** `installIonMobilityArray`'s replace branch is `self.getFloatDataArrays()[self.getIMData().first] = std::move(fda);` (`bind_spectrum.cpp:73`) — a `std::vector` move-assignment over a live `FloatDataArray` (`class FloatDataArray : public MetaInfoDescription, public std::vector<float>`), which deallocates the destination buffer. Any ndarray previously returned by `get_drift_time_array_view()` (`:450-459`, captures `arr.data()` at call time, `self_obj` only as keep-alive) or `FloatDataArray.get_data_view()` (`bind_kernel.cpp:1378-1385`) then reads freed memory with the nanobind ownership graph fully intact. Two lines: `v = s.get_drift_time_array_view(); s.set_drift_time_array(new_vals)`. Also reachable via `set_peaks(..., ion_mobility=...)` on the **default** `metadata="error"` policy, since `checkPeakMetadataAlignment` exempts the replaced IM index (`binding_utils.h:288-292`). R-6.2's ordering rule does not cover it — that rule moved the `"clear"` drop *after* the peak write; this free happens in the install step afterwards. **Fix: assign in place when sizes match.**
2. **R-8.1 is violated today on ~66 types.** nanobind never installs `tp_hash` (zero `Py_tp_hash` hits in `nb_type.cpp`/`nb_class.h`), and `.def` is a post-creation `PyObject_SetAttr` (`nb_func.cpp:487`) which — unlike a class body — does **not** null `__hash__` (measured on CPython 3.11.15). 105 classes bind `nb::self == nb::self`; 39 bind `__hash__`. So ~66 value-comparable types including `MSSpectrum` (`bind_spectrum.cpp:204`) are identity-hashable. `Peak1D` is worse: a content-based `__hash__` (`bind_kernel.cpp:1951-1956`) sitting directly below writable `mz`/`intensity` setters (`:1941-1948`).
3. **#9792 V.5's line numbers for the `python_error` flattening are stale post-#9804** — `:22,:32,:42,:52` are now two doc comments, a `static_assert` and a brace; the real sites are `:63, :164, :178, :188, :198`. And the consequence is worse than V.5 records: the original exception is **destroyed**, not preserved as `__context__`. A compiled replay of `e.restore()` followed by nanobind's `default_exception_translator` (`nb_internals.cpp:158-159`) yields `RuntimeError(...)` with `__context__ = None`, and `except ReferenceError:` does not catch it. Any `ReferenceError` contract is unenforceable inside the only callback path in the tree — and VI.5 assigns `RuntimeError` a *different* meaning.
4. **36 of the 356 `reference_internal` sites are policy no-ops** — primitives, `std::string`, `std::streampos`, `std::vector<double>`, and `DPosition<1>` (custom value caster, `openms_dposition_caster.h:28-66`). nanobind's value casters ignore `rv_policy`; there is no address to retain. So I.3's headline overstates the alias surface by ~10%, VI.11 step 1's grep over-counts, and it is a live R-10.1 defect (a reader of `bind_kernel.cpp:1350` reasonably concludes `fda[0]` aliases).
5. **Handle-returning bindings would be invisible to VI.11.** Every VI.2 row is keyed to an `rv_policy`; VI.11 step 1 greps `rv_policy::reference[^_]` and step 3 greps `-> *OpenMS::[A-Za-z:]* *&`; R-10.1 is scoped to *"Every Category S, H, B binding"*. A factory returning a handle object matches none of those keys, so the audit the issue exists to make mechanical would report the migrated accessors as non-aliasing. A sixth VI.2 category (or a checked manifest) is required.
6. **Cycle-collector invisibility.** `inst_traverse` visits only the instance `__dict__` and the type (`nb_type.cpp:48-55`), and `tp_traverse`/`tp_clear` are installed only under `dynamic_attr`/`is_weak_referenceable` (`:1264-1298`). Any bound type storing an `nb::object` member leaks the cycle `obj -> __dict__ -> handle -> owner -> obj`. Shared with nanobind's existing keep-alive edges, but multiplied one-per-handle.
7. **ODR trap for any caster-based design.** `arrow_zerocopy.cpp` and `main_module.cpp` include no caster header while 13 TUs include `type_casters/all_casters.h`, and `arrow_zerocopy.cpp:177, 251, 304, 321, 338` do `nb::cast<const OpenMS::MSExperiment&>` / `nb::cast<OpenMS::FeatureMap&>`. A `type_caster<T>` specialisation dropped into `all_casters.h` would silently not apply there — an IFNDR whose failure mode is the unguarded path.
8. **Empty-container asymmetry is not uniform.** R-7.3 documents `get_peaks_struct()` returning a detached array (`bind_spectrum.cpp:334-336`, `base is None`), but `_get_peaks_view`, `get_drift_time_array_view` and `get_data_view` pass a `nullptr` data pointer *with* `self_obj`, so `.base` is not None. Two adjacent APIs, different answer to "is this a live view?". Low severity (all four return zero-length arrays of the right dtype), but R-7.3 should name all four.

---

## 5. REVISED RECOMMENDATION

### Tier 1 — do now, independent of the ownership decision (all cheap, all verified)
1. **Fix `installIonMobilityArray`** (`bind_spectrum.cpp:73`): assign into the existing array in place when sizes match; document the free when they differ. ~10 lines. This is a live UAF in `develop` and is the single most actionable item in the audit.
2. **Peak `__iter__` → copy** — §7 step 1, endorsed. Amend to name all three peak-like types (`Peak1D`, `ChromatogramPeak`, `MobilityPeak1D`), matching VI.1 step 3's Category V. Justified independently by VII.1. Ship regardless.
3. **`__hash__ = None` sweep** over the ~66 `nb::self == nb::self` types with no `__hash__`; decide `Peak1D`'s content hash separately. R-8.1, live today.
4. **Strip the 36 inert `reference_internal` policies.** Zero behaviour change; makes every later audit number and every VI.11 grep mean what it says.
5. **Documentation debt already unblocked**: R-9.5 (the 62 GIL-release sites), R-7.1 (write-through invariants on the paths that already write through), R-7.3 extended to all four view bindings. Correct #9792 V.5's stale line numbers.
6. **Preserve exception type at the callback boundary** before any `ReferenceError`-based contract is designed, or R-5.1 is unenforceable there.

### Tier 2 — the measurement that actually decides this
7. **Benchmark object-level copy.** #9792 VII.1 measured 16-byte peak copies and explicitly left *"Real benchmark is a TODO"*. Nobody has measured what copying an `MSSpectrum` costs in a real mzML workflow, and the entire case for handles at the object level rests on that unmeasured number. Measure it against handle per-access overhead on the same data. Everything below branches on the result.
8. **Build VII.5 step 5 — the bulk table view.** Largest win, no storage change, no semantic commitment, and it is what makes copies affordable at *both* the element and object level. It also makes aligned-metadata operations correct by construction (R-6.3, R-6.1) — a correctness argument independent of the copy/view default. This should be the next substantial piece of work, ahead of any ownership machinery.

### Tier 3 — only if Tier 2 shows spectrum copies are unaffordable
9. Then, and only then, add handles as a **named, distinct-type, opt-in accessor** (`exp.spectrum_ref(i)`, `exp.iter_refs()`), leaving `__getitem__`/`__iter__` as owned values. That satisfies R-3.1/R-3.2 instead of overriding them, keeps `isinstance` intact for the value type, and confines the mechanism to the handful of sites where in-place editing is genuinely wanted. Required amendments, all from §2: owner must be Python-owning; version word bumped by the ~21 receiver-position replacement mutators; chain locator with per-level snapshots rooted at the root; explicit `.value()` for argument position; the full dunder protocol; a VI.2 category and a VI.11 step so the audit can see it.
10. Deprecation: one release with `pyopenms.set_option("element_access", ...)` and `DeprecationWarning` through the existing helper (`pyopenms/addons/deprecated_aliases.py:9-20`, already working and unused here), plus a migration table stating for each of `exp[i]`, `exp.getSpectrum(i)`, `for s in exp`, `exp.getSpectra()`, `spec[i]`, `for p in spec` what 3.5 did, what develop does, and what ships. VII.4 requires this; §7 has no such step.

### Do NOT attempt
- **The §6.1 single-type variant.** `nb::class_<T>` and `type_caster<T>` are mutually exclusive by static_assert (`nb_class.h:549-557`); keeping `class_` means retaining a raw pointer (collapsing §1) and leaving argument position unguarded; and it registers the handle at the element address, making it a participant in III.8 rather than a bystander.
- **`nb::implicitly_convertible<Handle, T>`.** Silently mutates a temporary destroyed at call return. Strictly worse than the `TypeError` it replaces. Make it a rule: *implicit conversions may never target a type passed by non-const reference.*
- **A generic handle sweep over the 356 sites.** 93 list-caster sites a policy change structurally cannot touch, 127 member accessors with no fingerprint, ~245 resolver expressions. That is a research project, not a rollout.
- **Removing any `keep_alive<0,1>`.** 7 of 13 sites are outside the tiering and `Param` cannot be index-converted at all. Keep VI.11 step 2's CI check, retargeted to assert `make_iterator ⇒ keep_alive<0,1>` for every unconverted site.
- **Adding `FREE_THREADED`** while any handle exists. Verified absent today (`CMakeLists.txt:290-294, 299-303, 334-338`); the one `NB_FREE_THREADED` hit in the tree is a defensive guard, not a build flag.
- **Treating §4 as deferrable.** The version word is what makes the fingerprint an identity check; deferring it defers the correctness of the L3 row, not just L4.

---

## 6. OPEN DECISIONS FOR HUMANS

1. **What is pyOpenMS for?** §6.5 names this as a product call and is right — but the evidence says it must be answered **first**, not last. If the answer is read-and-extract, always-copy plus the bulk table wins outright and none of Tier 3 is needed. §7 currently sequences the prototype ahead of the question.
2. **Is `isinstance(s, MSSpectrum)` breakable in a major release?** The audit finds no shape that preserves it soundly: a nanobind subtype must physically carry an `MSSpectrum` at offset 0 (`nb_type.cpp:1502-1506`, `:1538` — no pointer adjustment), the facade route costs the addon registry, the stub pipeline and ~205 lambdas, and the single-type variant is unsound. Pick which one to pay for.
3. **R-6.4: should pyOpenMS be stricter than C++?** Making public structured arrays read-only or detached closes the largest remaining unguarded surface and costs the zero-copy write path that §5 calls "the entire performance win". This is the same conflict VII.4 flags between hendrikweisser's principle and safety, and it is not resolvable technically.
4. **Should pyOpenMS promise a validated reference at all?** A guard that fires on shape and never on identity or invariants trains users to read silence as safety. §4 used exactly that reasoning at `:100-102` to disqualify generation counters. Applying the same standard to the handle is a judgement about what a partial guarantee is worth versus documented unsafety (VII.3's fallback).
5. **What deprecation budget exists?** The baseline users have is 3.5, where object access copied. Any object-level decision — handle *or* always-copy — is a silent semantic change against that baseline for code that is correct under any reading.

---

# Completeness critic — what the audit itself missed

# COMPLETENESS AUDIT OF THE AUDIT

Everything below is grounded in files I read in this checkout or in text I read in #9792. Nothing is run against pyopenms.

---

## A. COVERAGE CHECKLIST

Legend: **✓** examined substantively · **~** mentioned in passing only · **✗** never touched by any dimension.

### Levels
| | | Note |
|---|---|---|
| L1 | ✓ | G3 is the strongest result in the whole audit. |
| L2 | ✓ | |
| L3 | ✓ | G2. |
| L4 | ✓ | |
| L5 | ✓ | GIL sites + free-threading. |
| L6 | ~ | Only in §3's closing paragraph, as a rhetorical point ("trains the user that silence means safety"). No dimension asked what the handle does to `findNearest`. |

### Parts
| | | Note |
|---|---|---|
| I.1 | ~ | The **dual-`__getitem__` dead-code trap is untouched** and still live: `bind_processing.cpp:50` (by-value) shadows `:58` (`-> const DataFilter&`), same C++ signature, registration order decides. Any handle rollout adds overloads to exactly these sites. |
| I.2 | ~ | |
| I.3 / I.4 | ✓ | G4. Verified: `grep -o 'rv_policy::reference_internal' *.cpp` = **356**; incl. `binding_utils.h` = 361. |
| II | ~ | Only hroest/hendrikweisser. **jpfeuffer's 2020-08-02 bioconda-build-timeout objection is untouched** — and a per-container handle class × 13 modules is exactly a compile-time cost. |
| III.1–III.6 | ✓ | |
| III.7 | ✓ | |
| III.7 free-function 1-arg fast path | ✗ | A handle *factory* is naturally `m.def`/`def_static`. Non-overloaded, one unnamed arg ⇒ nanobind installs a keep-alive on arg 0 whether or not it owns the referent (`nb_func.cpp:1035`). Never checked against any proposed factory shape. |
| III.8 | ✓ | |
| IV | ✓ | 13 `make_iterator`, 13 `keep_alive<0,1>` — verified. |
| V.1 / #9795 | ✗ | |
| **V.2** | ✗ | **Decisive and missed.** Whether a sort reallocates is *data-dependent* (metadata arrays present ⇒ `select()` ⇒ new buffer). So the handle's **own alarm is data-dependent**: `sortByPosition()` raises `ReferenceError` on a spectrum carrying a float array and is silent on one that doesn't — same call, same corruption, opposite signal, and the caller cannot predict which. #9792 calls that "the worst possible property for a safety rule". Nobody applied it to the fingerprint. |
| V.3 / #9815 | ✗ | The dtype/`offsetof` contract underwrites the bulk tier the proposal keeps. |
| V.4 / #9803 | ✓ | |
| V.5 / #9804 | ✓ | |
| V.6 / #9809 | ~ | |
| VI.1 decision procedure | ✗ | The handle has **no category** in it. Step 1 asks "does the result alias memory not owned by the returned object?" — a handle owns nothing and aliases nothing, so the procedure terminates at Category V, which is wrong. Neither document notices the procedure doesn't classify its own proposal. |
| VI.2–VI.11 | ~/✓ | VI.2 (hazard 5) and VI.11 ✓; the rest partial. |
| VII.1–VII.3 | ✓ | |
| VII.4 | ~ | 3 of 4 bullets. **The 4th — "Is warn-on-questionable-operation right at all? It is pandas' `SettingWithCopyWarning`, which pandas is itself abandoning for copy-on-write" — is untouched.** That bullet is #9792 literally pointing at the alternative design (D1 below) and nobody followed it. |
| VII.5 | ~ | Steps 4/5 endorsed; 6 and 7 not analysed. |
| VII.6 | ✗ | Extending snake_case scalar properties to `MSChromatogram`/`Mobilogram` collides with R-3.5 and with §1's "keep property syntax for reads *and* writes". 18 `def_prop_rw` sites today, **0 with `reference_internal`** — verified. |

### Numbered rules (32 total)
Touched: R-3.1–3.5, R-4.2, R-4.5, R-4.6, R-4.7, R-5.1, R-6.1–6.4, R-7.1, R-7.3, R-8.1, R-9.1, R-9.5, R-10.1, R-10.2.

**Never touched (12 live):**

| Rule | Why it matters |
|---|---|
| **R-2.2** | *"Category V is the default. A binding may only leave V with a written justification in its docstring."* The proposal moves hundreds of bindings out of V by design and §7 has no docstring step. The rulebook's default-setting rule is the one rule the proposal must argue against, and neither document cites it. |
| **R-4.1 + R-9.2** | **The single largest normative miss.** R-4.1 requires the counter to live *"in a stable heap control block — **not on the Python wrapper**, which may be destroyed independently"*; R-9.2 repeats it. The proposal's `{owner,index,data_snap,size_snap}` lives entirely on the Python wrapper, and the audit's own G2 amendment ("a `uint64 version` … plus a `version_snap` in the struct") puts the counter's mirror there too — without ever checking it against the placement rule the issue already decided. |
| **R-4.3** | *"Classification is per operation per class, determined by reading the implementation."* The fingerprint is a deliberate refusal to classify per operation. Combined with V.2 this is a policy conflict, not a detail. |
| **R-4.4** | *"Default to blocking."* The handle defaults to **permitting everything and reporting afterwards**. That is a straight inversion of the rulebook's default and neither document names it. |
| **R-5.2** | *"Do not infer user intent in `__del__`."* This forecloses the obvious mitigation for always-copy's worst property (a silent no-op write). The audit recommends always-copy without noting that the rulebook has already banned the cheap warning. |
| **R-7.2** | Eight search bindings with an unenforced sortedness precondition. A guard that fires on `data()`/`size()` cannot fire on sortedness, so every handle-mediated position write leaves `findNearest` returning a wrong index or a false "not found". |
| **R-8.2** | *"Stable keys must be explicit immutable representations."* This is the rule that supplies the **replacement** for the identity/hash regression the audit found — and it points straight at machinery OpenMS already ships (see D3). Cited by nobody. |
| **R-9.3** | *"Operations touching two ownership graphs (`swap`) must lock both."* `swap` is the one operation carrying the L3 row. |
| **R-9.4** | GIL-releasing bindings audited *against every lease type*. |
| R-2.1, R-2.3, R-3.6 | Minor / settled. |

### Merged & open PRs
| PR | |
|---|---|
| #9794, #9803, #9804, #9808 | ✓ |
| #9809 | ~ |
| #9795, #9815 | ✗ |
| #9807 (open) | ✗ |
| **#9857** | ✗ — **and the audit re-recommends it as new work.** Tier 1 item 5 asks for "R-9.5 (the 62 GIL-release sites)" documentation. It was merged in `33593da` and ships as `src/pyOpenMS/THREAD_SAFETY.md`, with the 62-site inventory table already broken down 16/2/13/31. The issue's Status table is stale; nobody checked `git log`. |
| **#9908** | ✗ — merged after #9792's last update. Touches `bind_format.cpp`, adds per-peak IM aux arrays through imzML and `OnDiscImzMLExperiment`. It **widened hazard 1's producer surface** and shifted the audit's own citation: the move-assign is now `bind_spectrum.cpp:74`, not `:73`. |

---

## B. Synthesis claims resting on assertion rather than source or demo

1. **"18 of the 356 … are natively expressible (5.1%)"** and the whole six-way classification (127/93/40/35/25/12/24) — no method stated, no sample shown, no script. It is the headline number of the entire verdict.
2. **"233 mutable-`T&` argument occurrences on 224 registration lines, plus 204 non-self `const T&`"** — the audit corrects #9792's arithmetic with arithmetic that is itself uncited. (Same objection it raises against `:83-94`.)
3. **"an 8-hop path … genuine Category-H depth is 6, and 5 from FeatureMap"** — the graph is asserted; only the 5-hop `addons/consensusmap.py:98-107` idiom is cited, and that one is real (verified: `for f in cmap: … pep[0].getHits()[0]…`).
4. **"~245 resolver expressions"**, **"~21 receiver-position replacement mutators"**, **"~56 registered sequence types"**, **"~205 lambdas"** — four cost estimates driving the Tier-3 decision, none derived.
5. **"~66 value-comparable types are identity-hashable"** — the inputs check out (105 `nb::self==nb::self`, 39 `"__hash__"` — verified), but 105−39=66 assumes disjointness, which is not shown.
6. **The G3 compiled mock** ("`offsetof(spectra_) = 184`", "4/4 configurations") and the **`-O2` codegen claim** (GCC 13.3 / Clang 18 substituting the snapshot) — the strongest and the most fragile results in the audit, and neither ships a reproducer.
7. **"nanobind 2.12.0"** line citations throughout — #9792 warns explicitly that three nanobind trees exist here at 2.4.0 / 2.10.0 / 2.12.0 and that line numbers shift. There is no `flags.make` check recorded.

---

## C. Whole categories nobody attacked

**C1. Error/exception translation — the infrastructure does not exist at all.**
`grep -rn "register_exception_translator\|nb::exception\|PyErr_SetString"` over `bindings/**` returns **zero hits outside `nanobind_ms_data_consumer.h`**. So:
- `OpenMS::Exception::BaseException : public std::runtime_error` (`Exception.h:62`) ⇒ nanobind's default translator maps **every core exception to Python `RuntimeError`**, `IndexOverflow` included. VI.5 assigns `RuntimeError` the meaning *"container structurally modified during iteration"* — a slot already fully occupied. `except RuntimeError` cannot distinguish a mutation-during-iteration from a missing file. Nobody noticed, in either document or the audit.
- `pyopenms.set_option(...)` appears **nowhere in the tree** — verified. #9792 VI.5 uses it, and the audit's Tier-3 deprecation plan (item 10) is built on it. It is unbuilt infrastructure being spent twice.
- No `OpenMSViewInvalidationWarning` / `ImplicitCopyWarning` / `AlignmentWarning` classes exist either.
- And there is a **13-module problem**: a shared `ReferenceError`-like exception object must be created once and seen by all 13 `NB_MODULE`s under `NB_DOMAIN pyopenms`. Nobody said where it lives or which module owns it.

**C2. Build/packaging, stubs, and the module split.**
- `fix_stubs.py` (303 lines) + `stubs.patch` (25 lines) run as a post-build custom target (`CMakeLists.txt:387-399`). A new handle type needs a `stubs.patch` entry, and `fix_singleton_types`/`fix_cpp_type_annotations` show that stubgen already fails to render pyOpenMS-invented Python types.
- `pyopenms/__init__.py:194-256` imports the 13 submodules by name and re-exports; addon injection is **by class name** (`addons/__init__.py:88-92`). Nobody asked which module registers the handle class or what happens when `_pyopenms_experiment` needs a type owned by `_pyopenms_spectrum`.
- ABI/wheel: cibuildwheel ships manylinux/macOS/Windows wheels (`.github/workflows/pyopenms-wheels-cibuildwheel.yml`); a handle changes the exported type set of every wheel. `NOMINSIZE` (`CMakeLists.txt:292`) means `-O3`, so the audit's `-O2` codegen claim is measured at the wrong optimisation level for the shipping build.

**C3. Test strategy — nobody asked how any of this gets regression-tested, and the tests already contradict the proposal.**
- **`src/pyOpenMS/tests/unittests/test_mutable_references.py` (108 lines) exists and asserts today's aliasing as correct**: `exp[0].setRT(99.0)` must mutate `exp` (:27), `exp.getSpectra()[0].setRT(42.0)` must mutate `exp` (:43), plus `getHits()` (:59, :75), `fm[0]` (:90), `cm[0]` (:105). Under always-copy these six tests fail; under handles they pass but their docstrings ("must return a reference, not a copy") become the opposite of policy. Neither document mentions the file.
- The repo's leak detection is `tests/memoryleaktests/testAll.py`, a **free-memory delta with a 10% threshold** — structurally incapable of catching the audit's hazard 6 (one uncollectable `obj→__dict__→handle→owner→obj` cycle per handle).
- `tests/extract_api.py` + `tests/compare_api.py` compare the Cython-3.5 API JSON to current. That harness is the natural CI gate for exactly the deprecation table Tier-3 item 10 asks for, and nobody connected them.
- #9792's own two methodology warnings (a `git stash push --` "revert" that reverts nothing; a failed compile leaving a stale `.so`) apply with full force to a change whose entire observable is *an exception that sometimes does not fire*. A negative control for "the guard fires" is much harder than for "the test crashes", and no dimension proposed one.

**C4. Doctests / documentation.** `fix_stubs.py` has a whole `.. code-block::` reindenting pass, i.e. docstrings carry executable examples into the `.pyi`. Grep for `doctest` in `src/pyOpenMS` returns one hit, in a test helper. R-10.1 requires per-binding ownership prose on 356 sites; nobody costed the docstring diff or noticed there is no doctest runner to keep the examples honest.

**C5. `arrow_zerocopy.cpp` and `_dataframes.py` as a second, unmodelled export surface.**
- `arrow_zerocopy.cpp` takes `nb::object` and does **manual `nb::cast<const MSExperiment&>` / `nb::cast<FeatureMap&>`** at `:177,251,304,321,338,352,370,387,404,418` — a *third* conversion path, neither receiver nor registered argument, and 16 of the 62 GIL releases live here. A distinct handle type fails these casts at runtime, not compile time.
- `addons/msexperiment.py:43-49` slices `mz_flat[offsets[i]:offsets[i+1]]` from `get2DPeakDataSemiWide` and stores the slices **inside a pandas DataFrame** — comment says "zero-copy views". Whatever `get2DPeakDataSemiWide` returns is retained by a DataFrame with a lifetime nobody models. `_dataframes.py` and `addons/*.py` are the largest *consumer* of element semantics in the tree and neither document opens them.

**C6. The shipped always-copy control group nobody measured.** `OnDiscMSExperiment.getSpectrum(i)` returns **by value** (`bind_kernel.cpp:1826`), as does `OnDiscImzMLExperiment.getSpectrum` (`:1888`), and `getSpectrumByNativeId` is a key-addressed value accessor (`:1837`, `:1860`). An entire parallel container family already implements option A, on the largest files users have, at a per-call cost far above a vector copy — and nobody used it as evidence in a debate whose deciding number (VII.1's spectrum-copy TODO) is still unmeasured.

**C7. A guard-signal inversion nobody derived.** `MSExperiment::sortSpectra` is an unconditional in-place `std::sort` over `spectra_` (`MSExperiment.cpp:792-804`), bound with the GIL released (`bind_experiment.cpp:96`). Outer `data()`/`size()` unchanged ⇒ **spectrum handles stay silent while writing to the wrong spectrum**. But `std::sort` moves `MSSpectrum` objects, swapping their peak-buffer pointers ⇒ **every peak-level handle's inner fingerprint mismatches and raises**, even though its peak is alive, unmoved and uncorrupted. One call, and the guard fires exactly where there is no damage and stays quiet exactly where there is. Note also `MSExperiment` is `final` with **two** independent vectors (`chromatograms_` at `MSExperiment.h:1289`, `spectra_` at `:1291`) and no `data()`, so `size_t index` is under-specified for it at the root.

---

## D. Designs nobody proposed

**D1. Copy-on-write / detach-on-first-write handle.** `exp[i]` returns an object that re-resolves and reads through with zero copy, and on the *first write* materialises a private copy and detaches permanently. Reads pay nothing (the dominant path: `to_df`, `get_peaks`, extraction), writes never reach the container, so **the L4 residue degrades from "silent wrong write" to "stale read"** — and it exactly restores the release/3.5.0 baseline the audit names as the deprecation problem. #9792 VII.4 already points at it ("pandas is itself abandoning [`SettingWithCopyWarning`] for copy-on-write") and no dimension followed. **Verdict: dominates the §2 handle on the safety axis at equal read cost; its cost is that `s = exp[0]; s.setRT(5)` becomes a no-op, i.e. it is always-copy semantics at handle performance. The strongest unexplored option.**

**D2. Explicit transaction / context manager** — `with exp.edit(i) as spec: …`, handle invalidated on `__exit__`. Retention (the axis §1 says is the only one that matters) becomes lexically impossible, and since the handle is a distinct type that owns nothing, invalidation is a flag with no `inst_c2p` entry — sidestepping III.8 entirely, which is what defeated invalidate-on-return in V.5. **Verdict: the only design that makes §1's own thesis enforceable rather than aspirational; costs an idiom change and still does nothing about mutation from inside the block.**

**D3. Locate by stable identity, not by index — using identity OpenMS already ships.** `FeatureMap` and `ConsensusMap` both inherit `UniqueIdInterface` **and `UniqueIdIndexer<T>`** (`FeatureMap.h:75-76`, `ConsensusMap.h:64-65`), and `UniqueIdIndexer::uniqueIdToIndex()` is documented as **self-healing** (`UniqueIdIndexer.h:53,61,76,92`: on a miss it calls `updateUniqueIdToIndex()` and retries). `Feature`/`ConsensusFeature` carry the id via `RichPeak2D : UniqueIdInterface`. So on exactly the two containers where #9792 *measured* L4 (`fmap.sortByRT()` → `f0.setMZ(999)` lands on the wrong feature), **L4 can be closed exactly, with zero core change, zero annotation, and no ABA**, by keying the handle on `getUniqueId()` instead of `index`. Spectra have the analogous key in `getNativeID()` (`SpectrumSettings.h:114`). Guard against `UniqueIdInterface::INVALID == 0` (`UniqueIdInterface.h:37,45`) by falling back to a copy. This is also the R-8.2 answer to the identity/hash regression. **Verdict: strictly better than the fingerprint on the residue the audit calls "the mechanism", cheaper than a version word, already in the core — and nobody, in three documents, mentioned it.**

**D4. Immutable snapshot + explicit write-back** (`s = exp[i]` frozen; `exp.set_spectrum(i, s)` / `exp.update(i, fn)`). **Verdict: soundest and cheapest to build; pays one copy per edit and one API idiom; R-5.2 confirms the obvious "you forgot to write back" warning is off the table, so it must be taught, not detected.**

**D5. Columnar table facade** (VII.5 step 5). Half-built already: `get2DPeakDataSemiWide` (`bind_experiment.cpp:279`) plus the slicing in `addons/msexperiment.py:43-49`. **Verdict: the audit is right to promote it, but it inherits R-6.3/R-6.4 unresolved and the V.3/#9815 dtype contract nobody re-checked.**

**D6. Arena / slot-map with generation-tagged ids.** **Verdict: strictly worse than D3 — it invents the identity `UniqueIdIndexer` already provides, and it is a core storage change, i.e. hendrikweisser's objection.**

**D7. Python object as sole owner (`std::vector<nb::object>` storage).** **Verdict: non-starter — it is #4800's rejected proposal in a worse form, breaks every C++ algorithm, and requires the GIL to touch storage, which contradicts all 62 `gil_scoped_release` sites.**

---

## E. Three things I would add to the audit before it is delivered

1. **Check the rulebook's own placement rule before proposing a version word.** R-4.1/R-9.2 already decided the counter must live in a stable heap control block, not on the wrapper. Both the proposal and the audit's amendment put it on the wrapper.
2. **Re-read `git log` before writing recommendations.** Tier 1 item 5 asks for work merged in `33593da` and shipping as `THREAD_SAFETY.md`; #9908 moved a cited line number under the audit's feet.
3. **Open `tests/`.** `test_mutable_references.py` asserts, as a maintained regression suite, precisely the semantics both options break — and it is the concrete answer to "what deprecation budget exists?" that Open Decision 5 asks for.