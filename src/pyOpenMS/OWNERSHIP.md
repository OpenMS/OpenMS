# Copies and references in pyOpenMS

> **Every getter returns a copy. To change something, put it back.**

```python
s = exp[0]        # a copy
s.setRT(5.0)      # edits the copy
exp[0] = s        # now the experiment has it
```

Writing `exp[0].setRT(5.0)` on one line edits a temporary and the change is
discarded. So does `for s in exp: s.setRT(...)`. There is no accessor that
behaves differently, so there is nothing to memorise: if you did not put it
back, it did not happen.

## Nested edits need a write-back at each level

Each getter hands you a copy, so an edit two levels down has two levels to
climb:

```python
s = exp[0]
settings = s.getInstrumentSettings()
settings.setZoomScan(True)
s.setInstrumentSettings(settings)   # into the spectrum
exp[0] = s                          # into the experiment
```

## Reads are always safe

Reading through a temporary works — `exp[0].getInstrumentSettings().getZoomScan()`
is well-defined, because the parent stays alive for as long as the object you
took from it. Only writes are discarded.

The same holds while iterating. `for spec in exp:` walks the container by
index, taking an owned copy per step against the container's *current* size,
so mutating the container inside the loop is defined exactly like mutating a
Python list during iteration: appended elements are visited, shrinking ends
the loop early. (`Param` has no index access; its iterator snapshots owned
entry copies up front instead.)

Reads are not free, though. `exp[0].getRT()` copies a whole spectrum to read one
number, and `consumer.getData()` copies an entire experiment. Fetch once into a
variable rather than calling repeatedly in a loop, and for bulk numeric work use
the array API, which copies nothing per element:

```python
mz, intensity = spec.get_peaks()
spec.set_peaks(mz * 1.0001, intensity)
```

## Attributes are not getters

A handful of plain data-holder structs expose their members as Python
**attributes** rather than through `getX()`. Those behave the way Python
attributes normally do — they *are* the object, so edits land and you can assign
to them directly:

```python
pair.target.setRT(5.0)      # lands — `target` is an attribute, not a getter
pair.target = other         # also fine
```

So the full picture is two lines, and they match ordinary Python intuition:

> **A method hands you a copy — write it back with the matching setter.**
> **An attribute is the object itself — edit or assign it directly.**

This affects 31 members on structs such as `SiriusTargetDecoySpectra`,
`RangeSet`, `PreprocessedPairSpectra` and `AQS_featureConcentration`. The main
container classes have no such attributes.

## Three exceptions, all visible at the call site

| | Example | Why |
|---|---|---|
| **Fluent builders** | `ParquetFilter().eq("ms_level", 1).andNext()` | returns *itself* so the chain continues — that is the API |
| **Asking a database for its entry** | `ResidueDB().getResidue("A")`, `ModificationsDB().getModification(...)` | you asked the shared, process-lifetime database for *its* entry |
| **Shared pointers** | `spec.getDataProcessing()` | the list is a copy, but its entries are `shared_ptr`s that keep pointing at the same `DataProcessing` |

The first two are not a getter handing you part of an object it owns, which is
what the rule above is about. Note the second is narrow: it covers the database
singletons' own lookups. A database entry reached through *your* object —
`seq[0]`, `seq.getNTerminalModification()`, `residue.getModification()`,
`na_seq.getSequence()` — is a copy like everything else, so editing it cannot
corrupt the database for the rest of the process.

The third exception is a real one: `getDataProcessing()` copies the *list*, so
appending to or removing from the returned list does not touch the object — but
the `DataProcessing` entries inside it are shared by design (OpenMS deliberately
lets many spectra reference one provenance record instead of duplicating it), so
editing an entry in place *is* visible through the object:

```python
dp = spec.getDataProcessing()
dp[0].setMetaValue("note", "x")   # reaches spec, and every other spectrum sharing that record
```

Deep-copying on read would break that sharing, so the write-back rule still
applies here as the safe habit: build the list you want and call
`setDataProcessing()`. This affects the five `getDataProcessing()` getters
(`MSSpectrum`, `MSChromatogram`, `SpectrumSettings`, `ChromatogramSettings`,
`MetaInfoDescription`) — and, in the same category, the OpenSwath
interchange types, whose C++ data model is `shared_ptr` throughout:
`OSSpectrum`/`OSChromatogram` accessors (`getMZArray()`,
`getIntensityArray()`, `getTimeArray()`, `getDriftTimeArray()`,
`get_data_arrays()`) and the `ISpectrumAccess` implementations'
`getSpectrumById()`/`getChromatogramById()` hand back shared pointers
because sharing those arrays between pipeline stages without copying is
the point of that data model. Nothing else shares.

Outside these three, no *getter* hands out a live alias. `reference_internal`
does still appear in the bindings, in three places that are not getters: the
zero-copy numpy views (`data_view()`, `peaks_struct()` and friends), where the policy keeps
the owning object alive for as long as the view exists — which is what makes
those views safe; the element views of the `_view` family described in the
next section (`spectrum_view(i)`, `float_data_array_view(i)`, …), which opt
into aliasing by name; and the in-place operators and fluent builders above
(`AASequence.__iadd__`, `EmpiricalFormula.__iadd__` / `__isub__` /
`addChargeAdduct`, the `ParquetFilter` chain), which hand back the very object
they were called on. A future audit of `reference_internal` should expect all
three.

## Opting into views: the `_view` family

When copying is too expensive and you accept the hazards, ask for a **view** —
an accessor whose *name* announces that it aliases the object's storage. View
accessors never carry the `get_` prefix: `get_*` always hands you an owned
copy, `*_view` / `*_views` / `*_struct` always alias. The prefix *is* the
ownership contract.

| Form | Naming | Example |
|---|---|---|
| Single aliased element | `<singular>_view(i)` | `exp.spectrum_view(i)`, `fm.feature_view(i)`, `exp.chromatogram_view(i)` |
| List of aliased elements | `<singular>_views()` | `exp.spectrum_views()`, `fm.feature_views()`, `exp.chromatogram_views()` |
| Iterator of aliased elements | `iter_<singular>_views()` | `exp.iter_spectrum_views()`, `fm.iter_feature_views()`, `exp.iter_chromatogram_views()` |

The contract, stated once for the whole family:

* Edits through a view land **immediately** — no write-back needed.
* Every view keeps its parent object alive automatically, so a view never
  outlives the storage it points into.
* A view is only valid until the underlying list is **resized or reordered**:
  `addSpectrum`, `setSpectra`, `push_back`, `sort*` and friends invalidate every
  outstanding view into that list. Using a view after that is undefined
  behaviour — this is the #9792 hazard you are opting into, by name.

```python
for spec in exp.iter_spectrum_views():
    spec.setRT(spec.getRT() + shift)      # lands directly, nothing copied
```

The same convention already marks the zero-copy numpy views (`data_view()`,
`matrix_view()`, `peaks_struct()`); this extends it to object element access.
The family is
available on `MSExperiment` (spectra, chromatograms), `FeatureMap` (features),
`ConsensusMap` (consensus features), `PeptideIdentificationList`
(identifications), `MRMTransitionGroup` (features, chromatograms), and
`MSSpectrum`/`MSChromatogram` (float/integer/string data arrays) — where
`spec.float_data_array_view(i).data_view()` chains into a fully
zero-copy numpy view of spectrum-owned storage. Anything named `get_*`
returns a copy you own, as the rule above says (with only its three
documented exceptions).

pyOpenMS's own DataFrame exports (`to_df`, `get_ion_df`, ...) iterate views
internally, so the copy rule costs nothing on those paths.

## Why it works this way

Until release 3.5 every one of these returned a copy: the Cython wrapper owned
its C++ object outright and could not borrow. The nanobind port turned many of
them into live references — a side effect of the port rather than a decision —
and that made ordinary code unsafe. A spectrum taken from an experiment pointed
into the experiment's storage, so adding spectra moved that storage and left the
spectrum reading freed memory; sorting a feature map silently re-pointed a held
feature at a different feature; and a residue taken from a sequence aliased the
global `ResidueDB`, so editing it changed that residue for every sequence in the
process.

Returning owned values removes all of them, and restores what 3.5 did.
Aliasing was introduced and removed inside the same unreleased development
cycle, so no released version ever behaved differently and existing code
needs no changes.

See issue #9792 for the full analysis.
