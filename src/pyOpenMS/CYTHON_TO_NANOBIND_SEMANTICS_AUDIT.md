# Audit: did copy-vs-alias semantics change from Cython (3.5) to nanobind?

Scope: the specific question raised in
[#9792](https://github.com/OpenMS/OpenMS/issues/9792) — *"were spectra copied
before with the Cython bindings in 3.5 and are now not copied?"* This is the
**before/after** comparison that the RFC in #9792 references (Part II) but never
lays out surface-by-surface. It complements the RFC: the RFC audits the *current*
nanobind state and proposes a forward policy; this document establishes what the
port actually *changed* relative to `release/3.5.0`.

- **"3.5"** = `release/3.5.0`, the last Cython/[autowrap] release.
- **"nanobind"** = current `develop` (the hand-written nanobind bindings).

---

## Verdict

**Yes — the default access semantics were inverted for the container element
types.** In 3.5, obtaining a spectrum, chromatogram, feature or consensus feature
from its container gave you an **independent copy**; mutating it did nothing to the
container. On `develop` the same access gives you a **live alias** into the
container's storage; mutating it mutates the container (and can dangle — that is
the whole subject of #9792).

The user's intuition is correct, with three precise caveats that matter and are
easy to get wrong (detailed below):

1. The flip is **not universal**: `MSSpectrum.__getitem__` (peak-by-index) still
   returns a copy today, so peak *indexing* did **not** change — but peak
   *iteration* did.
2. The bulk numpy path (`get_peaks()`) was a copy in 3.5 and is **still** a copy;
   it never aliased.
3. Explicit zero-copy views (`_mv`/`_view`/`get_peaks_struct`) existed in **both**
   eras as an opt-in; they are not what changed.

So the headline is narrower and sharper than "3.5 copied, nanobind aliases": **the
implicit, unsuffixed object-access paths flipped from copy to alias, while the
explicit copy path and the explicit view path were both already there and stayed
put.**

---

## Root cause — why it flipped

The copy behaviour in 3.5 was not a designed policy; it was **forced by the
wrapper representation** and could not have been otherwise without hand-written
code. Every autowrap wrapper owns its C++ instance through a smart pointer:

```cython
cdef shared_ptr[_MSSpectrum] inst        # every wrapper: sole owner of a heap C++ object
```

(Visible throughout the 3.5 addons: `self.inst.get()`, and
`py_result.inst = shared_ptr[_MSSpectrum](_r)`.)

A wrapper in this model *cannot borrow* a subobject owned by another wrapper — it
must hold a `shared_ptr` to a heap object it owns. Therefore any method returning a
wrapped type has to **heap-allocate a fresh copy** and wrap that. This is directly
visible in the two hand-written accessors in 3.5 (`addons/MSExperiment.pyx`):

```cython
def getSpectrum(self, id_):
    cdef _MSSpectrum * _r = new _MSSpectrum(self.inst.get().getSpectrum(<size_t>id_))  # explicit COPY
    cdef MSSpectrum py_result = MSSpectrum.__new__(MSSpectrum)
    py_result.inst = shared_ptr[_MSSpectrum](_r)
    return py_result
```

The autowrap-generated `operator[]` and `__iter__` produce the *same*
`shared_ptr`-owned wrapper, so they copy identically. The **only** way 3.5 ever
handed out an alias was a hand-written method that returned a Python `memoryview`
or numpy array (not a wrapped C++ type) — always `_mv`/`_view`-suffixed and
documented "fast, unsafe" (see `test_MSSpectrum.py`: `get_data()` "returns copy -
safe" vs `get_data_mv()` "memory view - fast, unsafe").

nanobind removed the structural constraint: `rv_policy::reference_internal` lets a
returned wrapper *borrow* a subobject and merely keep the parent alive. The
hand-written nanobind bindings then declared element access as `-> T&` with
`reference_internal`, which produces aliases. **The C++ core returned references in
both eras** (`MSSpectrum& operator[]`, `MSSpectrum& getSpectrum(Size)`); what
changed is the binding layer's decision about whether to copy that reference or
expose it. Per #9792 Part II, that decision "is not a considered design decision;
it is a side effect of the nanobind port, noticed once, contested [by @jpfeuffer
on #4800], and left unresolved."

---

## The audit — surface by surface

Copy = mutating the returned object leaves the container unchanged.
Alias = mutating it mutates the container (and may dangle per #9792 L1–L4).

### Object containers (the meaningful cases)

| Access path | 3.5 (Cython) | nanobind (`develop`) | Changed? | nanobind site |
|---|---|---|---|---|
| `exp.getSpectrum(i)` | **copy** (`new _MSSpectrum(...)`) | **alias** (`-> MSSpectrum&`, `reference_internal`) | **YES** | `bind_experiment.cpp:151` |
| `exp[i]` | **copy** (autowrap `operator[]`) | **alias** (`-> MSSpectrum&`, `reference_internal`) | **YES** | `bind_experiment.cpp:143` |
| `for s in exp` | **copy** (autowrap iter) | **alias** (`make_iterator<reference_internal>`) | **YES** | `bind_experiment.cpp:140` |
| `exp.getSpectra()` | **list of copies** (returns `vector[MSSpectrum]` by value) | **list of aliases** (`-> vector<MSSpectrum>&`; caster propagates policy per element) | **YES** | `bind_experiment.cpp:109` (+ #9792 I.4) |
| `exp.getChromatogram(i)` / `getChromatograms()` | **copy** / list of copies | **alias** / list of aliases | **YES** | `bind_experiment.cpp:115` / `:114` |
| `fmap[i]`, `for f in fmap` | **copy** | **alias** (`-> Feature&`) | **YES** | `bind_kernel.cpp:3598` / `:3596` |
| `cmap[i]`, `for cf in cmap` | **copy** | **alias** (`-> ConsensusFeature&`) | **YES** | `bind_kernel.cpp:3392` / `:3390` |
| `chrom[i]`, `for p in chrom` | **copy** | **alias** (`-> ChromatogramPeak&`) | **YES** | `bind_chromatogram.cpp:137` / `:135` |

### Peak-level access on MSSpectrum (the inconsistent one)

| Access path | 3.5 | nanobind | Changed? | nanobind site |
|---|---|---|---|---|
| `spec[i]` (Peak1D by index) | **copy** | **copy** ("Return by value (copy)") | **NO** | `bind_spectrum.cpp:402` |
| `for p in spec` (Peak1D iter) | **copy** | **alias** (`make_iterator<reference_internal>`) | **YES** | `bind_spectrum.cpp:285` |
| `spec.get_peaks()` (numpy) | **copy** | **copy** | **NO** | `bind_spectrum.cpp:345` |
| `spec.get_peaks_struct()` / `_view` | view (opt-in, `_mv` era) | view (opt-in) | ~same (explicit) | `bind_spectrum.cpp:319` |

### What did **not** change (stated to avoid overclaiming)

- **Input/write semantics.** `addSpectrum(spec)`, `push_back(peak)`, `set_peaks(...)`
  copy *into* the container in both eras. `exp[i] = spec` writes a copy into the
  slot in both eras (`bind_experiment.cpp:147`). The port changed *reads*, not
  *writes*.
- **The bulk numpy copy** (`get_peaks`, `get_mz_array`, `get_intensity_array`,
  `FloatDataArray.get_data`) — copy then and now.
- **Explicit views existed in 3.5 already** (`get_data_mv`, `get_mz_array_view`,
  `get_drift_time_array_mv`; see `performance.rst`, `test_MSSpectrum.py`). nanobind
  renamed `_mv`→`_view` (#8857) and added `get_peaks_struct`, but the *category* —
  an opt-in, explicitly-named alias — is not new.
- **Scalar data-array indexing** (`fda[i]` → a number). nanobind returns `float&`
  but Python numbers are immutable values, so no observable aliasing either way.

---

## The internal inconsistency the port introduced

Within a single class, the two peak-access paths now **disagree**:

```python
p = spec[0]        # COPY  -> p.setMZ(9) does nothing to spec
for p in spec:     # ALIAS -> p.setMZ(9) mutates spec[0]
```

This is `MSSpectrum` only, and it is a direct artifact of the migration: peak
*indexing* was deliberately kept as a copy (restored via #9777's revert — #9792
I.1), while peak *iteration* took the `make_iterator<reference_internal>` default.
No other container has this split; the rest alias on both paths. In 3.5 both paths
were copies, so the class was self-consistent.

---

## Why it matters (link to #9792)

Every hazard catalogued in #9792 Parts I–V is downstream of exactly this flip. In
3.5, `spec = exp[0]` was a detached copy, so none of L1–L4 could arise: there was
no alias to dangle (L1/L2), no owner to reparent under it (L3), and no index whose
meaning could drift (L4). The port made the common, unsuffixed access path an
alias, which is what turned these into live defects — e.g. #9792's measured
`fmap.sortByRT()` silently-wrong-write (L4) and the `addSpectrum`-realloc
use-after-free (L2). The behaviour hroest reported on
[#4800](https://github.com/OpenMS/OpenMS/issues/4800) in 2020 (`s =
exp.getSpectrum(5); s.setRT(5.0)` "modifies a copy") is the *same call* that today
modifies the experiment — the exact inversion this audit documents.

There is also a backwards-compatibility consequence worth stating plainly: **code
written against 3.5 that relied on copy semantics changes meaning silently on
nanobind.** A 3.5 script that did `s = exp[i]; s.setRT(t)` as a deliberate no-op
scratch edit now mutates the experiment; conversely, a 3.5 write-back idiom `s =
exp[i]; s.setRT(t); exp[i] = s` still works but the middle object is now live. This
is the "all three options are breaking … needs a deprecation cycle" point in #9792
VII.4, seen from the 3.5 side.

---

## Evidence appendix (exact citations)

**3.5 / Cython — `release/3.5.0`**

- Wrapper owns via `shared_ptr` → copy is forced:
  `addons/MSExperiment.pyx` `getSpectrum`/`getChromatogram`
  (`new _MSSpectrum(...)`, `new _MSChromatogram(...)`).
- `pxds/MSExperiment.pxd`: `MSSpectrum& operator[](size_t)` (autowrap → copy),
  `MSSpectrum getSpectrum(Size) # wrap-ignore` (→ hand addon copy),
  `libcpp_vector[MSSpectrum] getSpectra()` (by value → list of copies),
  `begin()/end() # wrap-iter-begin:__iter__(MSSpectrum)` (→ copies).
- `pxds/MSSpectrum.pxd`: `Peak1D& operator[](size_t)` (→ copy),
  `begin()/end() # wrap-iter:__iter__(Peak1D)` (→ copies).
- `pxds/FeatureMap.pxd`: `Feature & operator[](size_t)`, `__iter__(Feature)`.
- Copy-vs-view philosophy in tests: `test_MSSpectrum.py` — `get_data()` "(returns
  copy - safe)" asserts `data_copy[0]=100; fda[0]==1.0`; `get_data_mv()` "(memory
  view - fast, unsafe)" asserts `data_view[0]=100; fda[0]==100.0`.
- Independent behavioural corroboration: hroest on #4800 (2020-07-31),
  `getSpectrum` "modifies a copy".

**nanobind — `develop`**

- `bind_experiment.cpp`: `__iter__` `:140`, `__getitem__` `:143` (`-> MSSpectrum&`),
  `__setitem__` `:147`, `getSpectrum` `:151`, `getSpectra` `:109`,
  `getChromatograms` `:114`, `getChromatogram` `:115` — all `reference_internal`.
- `bind_spectrum.cpp`: `__iter__` `:285` (alias), `__getitem__` `:402` (`return
  self[i]; // Return by value (copy)`), `get_peaks_struct` `:319`, `get_peaks` `:345`.
- `bind_chromatogram.cpp`: `__iter__` `:135`, `__getitem__` `:137` (`-> ChromatogramPeak&`).
- `bind_kernel.cpp`: FeatureMap `:3596`/`:3598`, ConsensusMap `:3390`/`:3392`,
  Mobilogram `:1618`/`:1620`.
- Runtime confirmation of the alias behaviour: #9792 (measured on `develop`
  throughout Parts I–V, and the comment thread's `fmap.sortByRT()` /
  `addSpectrum`-realloc demonstrations).

---

## Method and limitations

- The 3.5 side is established from checked-in source (pxd declarations + the
  hand-written addons) and the wrapper representation, corroborated by the 3.5
  test suite and hroest's 2020 first-hand report. I did **not** execute a 3.5
  build — no `pyopenms` is installed or built in this container — so the Cython
  behavioural claims are source/representation-level, not freshly re-run here. The
  two accessors that are hand-written (`getSpectrum`/`getChromatogram`) are direct
  witnesses; `operator[]`/`__iter__` copy is inferred from the identical
  `shared_ptr`-owned representation they all share.
- The nanobind side is established from checked-in `develop` source and is
  corroborated by the runtime measurements already recorded in #9792.
- Line numbers are against the checkout at audit time and will drift.
