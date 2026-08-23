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

## Three further exceptions, all visible at the call site

| | Example | Why |
|---|---|---|
| **Fluent builders** | `ParquetFilter().eq("ms_level", 1).andNext()` | returns *itself* so the chain continues — that is the API |
| **Database lookups** | `ResidueDB`, `ModificationsDB` entries | a shared, process-lifetime entry, not part of your object |
| **Shared pointers** | `spec.getDataProcessing()` | the list is a copy, but its entries are `shared_ptr`s that keep pointing at the same `DataProcessing` |

The first two are not a getter handing you part of an object it owns, which is
what the rule above is about.

The third is: `getDataProcessing()` copies the *list*, so appending to or
removing from the returned list does not touch the object — but the
`DataProcessing` entries inside it are shared by design (OpenMS deliberately
lets many spectra reference one provenance record instead of duplicating it),
so editing an entry in place *is* visible through the object:

```python
dp = spec.getDataProcessing()
dp[0].setMetaValue("note", "x")   # reaches spec, and every other spectrum sharing that record
```

Deep-copying on read would break that sharing, so the write-back rule still
applies here as the safe habit: build the list you want and call
`setDataProcessing()`. This affects the five `getDataProcessing()` getters
(`MSSpectrum`, `MSChromatogram`, `SpectrumSettings`, `ChromatogramSettings`,
`MetaInfoDescription`) and nothing else.

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

See issue #9792 for the full analysis.

## Migrating from 3.6.0

If your code assigned through a getter and relied on it landing, add the
write-back:

```python
# before
exp[0].setRT(5.0)
for s in exp: s.setRT(0.0)
spec.getPrecursors()[0].setMZ(500.0)

# after
s = exp[0]; s.setRT(5.0); exp[0] = s
for i in range(len(exp)):
    s = exp[i]; s.setRT(0.0); exp[i] = s
p = spec.getPrecursors(); p[0].setMZ(500.0); spec.setPrecursors(p)
```

Nothing raises when you get this wrong — the write is simply discarded — so it
is worth grepping for `get...().set...(` and for assignments to `[...]` results
when upgrading.
