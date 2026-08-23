# Proposal: safe, consistent ownership via self-validating handles

Companion to `CYTHON_TO_NANOBIND_SEMANTICS_AUDIT.md`; addresses the open
design question in [#9792](https://github.com/OpenMS/OpenMS/issues/9792)
Part VI/VII (the "ownership policy" the issue exists to decide, and its
blocked step 4).

---

## 1. The question that started this

> *"Is there a way to still directly set members via a property without any of
> the known hazards? E.g. if a user can't acquire a reference / get the values
> but sets via a reference — would that work?"*

**Read/write direction is not the safety axis.** A write-only reference is
*more* dangerous than a read-write one, not less: a stale read may produce a
visibly absurd number, while a stale write silently corrupts a neighbouring
object and is never noticed. #9792's own measured L4 case is exactly a bad
*write* (`f.setMZ(999)` landing on the wrong feature).

The axis that actually matters is **retention**: is a raw C++ address stored in
a Python object across time? Every hazard in #9792 Parts III–V is downstream of
that one act.

So the answer is better than the question hoped for: once you stop retaining
addresses, **both** directions become safe, and you keep property syntax for
reads *and* writes.

## 2. The mechanism: re-resolve instead of retain

Store a *locator*, not a pointer, and re-derive the address on every access:

```cpp
struct Handle {
    nb::object owner;      // strong Python ref to the container  -> L1
    size_t     index;      // logical position
    const void* data_snap; // container.data() when handle was made -> L2/L3
    size_t      size_snap; // container.size() when handle was made
};
```

Every attribute get *and* set performs:

```cpp
auto& c = nb::cast<Container&>(owner);
if (c.data() != data_snap || c.size() != size_snap)
    throw ReferenceError("container was structurally modified");
if (index >= c.size()) throw nb::index_error();
return c[index].getRT();       // freshly resolved address
```

What this buys, against #9792's own level taxonomy:

| Level | Failure | Handled? | How |
|---|---|---|---|
| **L1** owner freed | UAF | **yes** | strong Python ref; ordinary refcounting |
| **L2** realloc moves storage | UAF | **yes** | `data()` changed → `ReferenceError` |
| **L3** `swap` reparents storage | UAF | **yes** | `data()` changed (measured in #9792 L3: `0x…470 → 0x…a70`) |
| **L4** in-place sort permutes | silent wrong write | **no** | ptr+size unchanged; needs §4 |
| bounds | OOB | **yes** | explicit check, not `OPENMS_PRECONDITION` |

**All memory-safety failures are closed. The residue is purely semantic** — a
handle after an in-place sort writes to a valid, live object, just not the one
you meant.

Two properties worth calling out:

- **It is self-validating.** Nothing has to be annotated, registered or
  released. There is no lease table, no refcount to leak, no admission
  protocol — so #9792's R-9.1 check-then-act race cannot arise, because
  there is no check-then-act.
- **The stale case raises instead of corrupting.** Today it segfaults (L2) or
  writes to the wrong object (L4).

## 3. Why a fingerprint, not a generation counter — the load-bearing finding

#9792 assumed the mechanism must be a generation counter in a heap control
block, and correctly judged that "may be too invasive" (contributor comment,
2026-08-22). It is worse than invasive — **for this codebase it cannot be made
sound by annotation.** Counted in the current bindings:

| Container | mutable-ref sites | as receiver (`self`) | **as argument** |
|---|---|---|---|
| MSExperiment | 192 | 78 | **114** |
| MSSpectrum | 219 | 92 | **127** |
| FeatureMap | 128 | 40 | **88** |
| ConsensusMap | 129 | 49 | **80** |
| MSChromatogram | 99 | 53 | **46** |
| | | | **≈455** |

*(`grep -o "OpenMS::T *& *ident"` minus `const`, minus `& self`; counts
overloads separately, so approximate — but the order of magnitude is the
point.)*

A generation counter needs a bump at **every mutation site**. The ~312
receiver-position sites are enumerable. The **~455 argument-position** sites
are not realistically annotatable: these are algorithms — `filterPeakMap(exp)`,
`pickExperiment(exp)` — that take your container and restructure it, with no
`self` to hang a hook on, spread across every `bind_*.cpp`. Miss one and the
counter is silently wrong, which is *worse than no counter*, because it
converts "unsafe" into "believed safe".

The fingerprint inverts the burden: it reads the container's **actual current
state** at use time, so an unannotated third-party mutator is caught for free.
Coverage does not depend on anyone remembering anything.

**It needs no core change.** `MSExperiment::getSpectra()` is public
(`MSExperiment.h:1214`) and `MSSpectrum` re-exports `data()` and `size()`
(`MSSpectrum.h:148`, `:137`). hendrikweisser's principle from #4800 — do not
degrade the C++ library for the bindings' benefit — is respected exactly, which
is the conflict #9792 VII.4 flagged as needing a decision. It does not need
one: the memory-safety half is achievable with zero core impact.

Cost per access: one type check plus two integer compares, against Python's own
attribute-dispatch overhead. Structurally negligible — but *should be measured*,
not asserted, before rollout.

## 4. L4 (index identity), if it is wanted later

The fingerprint deliberately misses in-place permutation, and one ABA case
(`clear()` + refill to the same size, allocator reuses the block). Both are
memory-safe and semantically wrong.

Closing them needs a real counter, and that decision can be **deferred** — it
is a strict, additive refinement of the same handle, not a redesign. Note the
argument-position problem above applies in full, so the honest options are
a conservative bump (any binding taking a mutable ref bumps; over-invalidates,
but sound) or accepting L4 as documented behaviour. Recommend shipping §2
first and deciding L4 on evidence.

## 5. Recommended tiering

Not "always copy". Copy is right at one level and wrong at another, which is
why a single global rule keeps failing:

| Level | Recommendation | Rationale |
|---|---|---|
| **Peaks** (`Peak1D`, 16 B) | **copy**, both `__getitem__` *and* `__iter__` | copying is cheaper than any guard; #9792 VII.1 |
| **Objects** (`MSSpectrum`, `Feature`, `ConsensusFeature`) | **handle** (§2) | copying a spectrum copies its whole peak buffer; this is where in-place editing is genuinely wanted |
| **Bulk numeric** | keep zero-copy `_view`/`_struct` | the entire performance win lives here |

Two consequences worth stating:

**Peak-level: this restores 3.5.** In 3.5 both peak paths copied. Making
`__iter__` copy is not a new invention — it un-breaks the split introduced by
the port, where `spec[0]` copies (`bind_spectrum.cpp:402`) but `for p in spec`
aliases (`:285`). It is small, self-contained, and independent of everything
else here. **This is the one change worth making regardless of how the rest is
decided.**

**Object-level: source-compatible for correct code.** Code that does
`s = exp[0]; s.setRT(5)` keeps behaving identically — it still mutates `exp`.
Code that was *silently broken* now raises. Compare "always copy", which
silently inverts that same line into a no-op — reintroducing precisely the
complaint hroest raised on #4800 in 2020, at the level where copies are most
expensive.

**Bonus: iteration gets simpler, not harder.** Index-based iteration yielding
handles retains no C++ iterators at all, so the iterator-invalidation hazard
class (#9792 Part IV / R-4.2 row 3) stops existing rather than being guarded —
`keep_alive<0,1>` at 13 sites becomes unnecessary. And mutation during
iteration produces a clean `RuntimeError`, which is Python's own idiom for
`dict`/`list` and what #9792 VI.5 already asks for.

## 6. Honest costs and open decisions

1. **Type identity.** Does `exp[0]` return `MSSpectrum` or `MSSpectrumHandle`?
   A distinct type is cleaner and satisfies R-3.3, but breaks
   `isinstance(s, MSSpectrum)` and every existing annotation. One type with two
   internal states (owning vs. bound) preserves compatibility at the cost of an
   indirection on every `MSSpectrum` method. **Recommend the single-type
   variant** — R-3.3 was written for buffer view-vs-value, where confusing the
   two is a data-corruption risk; here both states are safe, so the rule's
   rationale does not transfer. This is the main thing to decide before coding.
2. **numpy views are not covered.** A retained `get_peaks_struct()` ndarray
   still holds a raw pointer; the handle does not fix it. Needs buffer-export
   pinning (#9792 step 7) — genuinely separate work.
3. **Detach must be explicit and cheap**: `handle.copy()` → owned value
   (R-3.4); `copy.copy(handle)` likewise.
4. **Nested handles** (a `FloatDataArray` inside a spectrum inside an
   experiment) must validate the whole chain, depth ≤ 3 — this is R-4.5, and
   chain validation is the natural form of it.
5. **The counter-position deserves a real answer.** @jcharkow argued for
   copies-everywhere plus warnings, accepting that complex in-place workflows
   (PyTOPP) move off pyOpenMS. That is coherent and much cheaper to build. The
   case against is the spectrum-level copy cost, which is what #4800 was
   *originally* opened about. If the project's answer to "what is pyOpenMS for"
   is read-and-extract rather than workflow-building, always-copy wins and this
   proposal is over-engineering. **That is a product call, not a technical one,
   and it should be made explicitly rather than by default.**

## 7. Suggested sequencing

1. **Peak `__iter__` → copy.** Small, independent, restores 3.5, removes the
   worst inconsistency. Ship regardless of the rest.
2. **Prototype the handle on one container** (`MSExperiment` → `MSSpectrum`)
   behind the existing accessors; measure the per-access overhead for real.
3. **Decide §6.1** (type identity) on the prototype.
4. **Roll out** to `FeatureMap` / `ConsensusMap` / `MSChromatogram`.
5. **Then** revisit L4 (§4) and buffer pinning (§6.2) on evidence.
