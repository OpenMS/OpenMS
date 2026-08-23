# Copies and references in pyOpenMS

When you get something out of a pyOpenMS object, do you get your own copy, or a
live reference into the original? The answer decides whether your edits stick.

## The short version

> **Getting an element out of a container gives you a copy. To change the
> container, put it back.**

```python
s = exp[0]        # a copy
s.setRT(5.0)      # edits the copy
exp[0] = s        # now the experiment has it
```

Writing `exp[0].setRT(5.0)` in one line edits a temporary copy and the change is
discarded — silently. The same is true of `for s in exp: s.setRT(...)`.

## The two kinds of accessor

| Kind | Examples | You get | Writes reach the parent? |
|---|---|---|---|
| **Container access** | `exp[0]`, `for s in exp`, `at`/`front`/`back`, `getSpectrum(i)`, `fmap[0]`, `spec[0]` | an owned **copy** | no — write back |
| **Collection getter** | `getPrecursors()`, `getHits()`, `getFloatDataArrays()`, `getSpectra()` | a list of owned **copies** | no — write back with the matching setter |
| **Single sub-object getter** | `getInstrumentSettings()`, `getAcquisitionInfo()`, `getSourceFile()` | a live **reference** | **yes** |

So in a chain, **a write lands exactly when it does not cross a container access
or a collection getter**:

```python
spec.getInstrumentSettings().setZoomScan(True)     # lands   (spec is yours)
exp[0].getInstrumentSettings().setZoomScan(True)   # LOST    (exp[0] is a copy)
spec.getPrecursors()[0].setMZ(999.0)               # LOST    (getPrecursors copies)
```

**Reads are always safe.** Reading through a temporary — `exp[0].getInstrumentSettings().getZoomScan()`
— works: the parent stays alive for as long as the object you got from it.

## Doing it correctly

```python
# one element
s = exp[0]
s.setRT(77.0)
s.getInstrumentSettings().setZoomScan(True)   # nested edits ride along
exp[0] = s

# every element
for i in range(len(exp)):
    s = exp[i]
    s.setRT(s.getRT() / 60.0)
    exp[i] = s

# a collection
precs = spec.getPrecursors()
precs[0].setMZ(999.0)
spec.setPrecursors(precs)
```

For bulk numeric work, don't loop at all — use the array API, which is both
faster and unaffected by any of the above:

```python
mz, intensity = spec.get_peaks()
spec.set_peaks(mz * 1.0001, intensity)
```

## Why it works this way

Until release 3.5 every one of these returned a copy: the Cython wrapper owned
its C++ object outright and could not borrow. The nanobind port made many of
them live references, which was a side effect of the port rather than a
decision — and it made ordinary code unsafe. A spectrum obtained from an
experiment pointed into the experiment's storage, so growing the experiment
moved that storage and left the spectrum reading freed memory; sorting a feature
map silently re-pointed a held feature at a different feature. Returning owned
values removes both, and restores what 3.5 did.

See issue #9792 for the full analysis.

## Known rough edge

The table above has two rules, not one, and which one applies is not visible at
the call site: `getPrecursors()` copies while `getInstrumentSettings()` does not,
and nothing in the names says so.

The split is real but it is not a rule about *your* code. A collection getter
produced one alias per element, while a single sub-object sits at a fixed offset
inside its parent and cannot be invalidated on its own — so the two needed
different treatment for **memory safety**. What you experience is a different
question, "does my write stick?", and on that axis the split is arbitrary.

The single sub-object getters are the last group still returning references
(~166 bindings). Converting them too would reduce this document to its first
sentence. Two things worth knowing about what that would cost:

* Across this entire repository, exactly **two lines** of Python chain a write
  through such a getter (`chrom.getPrecursor().setMZ(...)` and its `getProduct()`
  sibling, in one test).
* The idiom usually cited as the reason not to — `algo.getParameters().setValue(...)`
  — **already does not work**, and no code here uses it.
  `DefaultParamHandler::setParameters()` applies defaults, validates and then
  calls `updateMembers_()`, which is where an algorithm copies parameter values
  into its own members. Editing the Param in place skips all of that, so the
  algorithm never sees the change. `setParameters()` was always required.

Until the decision is made, prefer the write-back form everywhere: it is correct
under either rule, and for parameters it is the only form that has ever worked.
