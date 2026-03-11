# Multimer Detection in MetaboliteAdductDecharger

## Problem

MetaboliteAdductDecharger cannot detect or annotate multimeric adducts like [2M+H]+.
The existing PR #8431 adds a `mol_multiplier_` field to Adduct but this is purely cosmetic
annotation — it does not change the mass model, labels all protonated features as dimers if
configured that way, and introduces hash/equality bugs. We need a correct, mass-aware
implementation.

Resolves: [#7151](https://github.com/OpenMS/OpenMS/issues/7151)
Supersedes: [PR #8431](https://github.com/OpenMS/OpenMS/pull/8431)

## Design Decisions

1. **Full mass-aware detection** — not just cosmetic annotation
2. **`max_multimer` parameter as sole control** (default 1 = disabled, preserves backward compat)
3. **No `mol_multiplier_` on Adduct** — multiplier is a property of the feature pairing, not the adduct
4. **Multiplier stored on ChargePair** (the edge between feature pairs)
5. **Linear scan for cross-multiplier queries** — simple, correct, optimize later if needed
6. **Multimer penalty parameter** — ILP solver can prefer monomer explanations

## Mass Equation

For two features of the same analyte M with charges q1, q2 and multipliers n1, n2:

```
mz1 × |q1| = n1×M + adducts1_mass
mz2 × |q2| = n2×M + adducts2_mass
```

Setting M equal from both sides and cross-multiplying:

```
n2 × (mz1×|q1|) - n1 × (mz2×|q2|) = n2 × adducts1_mass - n1 × adducts2_mass
```

M cancels exactly. No estimation needed. The left side is computed from observed data and
candidate (n1, n2). The right side is computed from the compomer's per-side adduct masses.

For n1 == n2, this reduces to `m2 - m1 = adducts2 - adducts1` — the existing algorithm unchanged.

## Components Changed

### 1. MetaboliteFeatureDeconvolution

New parameters:
- `max_multimer` (Int, default 1): Maximum molecular multiplier. Setting to 2 enables dimer
  detection, 3 enables trimers, etc. Default 1 preserves existing behavior.
- `multimer_log_penalty` (double, default -2.0): Log-probability penalty per multiplier step
  for cross-multiplier edges. Applied as `(max(n1,n2) - 1) × penalty` when scoring edges.
  Makes the ILP solver prefer monomer explanations unless mass evidence strongly supports a
  multimer.

Changes to `candidateEdges_()`:
```
for each feature pair (f1, f2), for each (q1, q2):
  m1 = mz1 × |q1|
  m2 = mz2 × |q2|

  for n1 = 1 to max_multimer:
    for n2 = n1 to max_multimer:
      if n1 == n2:
        // existing path: M cancels, binary search
        naive_mass_diff = m2 - m1
        query(net_charge, naive_mass_diff, tolerance, ...)
      else:
        // new path: exact cross-multiplier equation, linear scan
        queryMultimer(net_charge, m1, m2, n1, n2, tolerance, ...)

      process hits -> create ChargePairs with (n1, n2)
      apply multimer_log_penalty to edge score
```

Symmetry: n2 >= n1 avoids duplicate pairs. The LEFT/RIGHT compomer sides handle the
directionality (feature 1 as monomer + feature 2 as dimer, and vice versa).

Tolerance for cross-multiplier: `tolerance = n2 × tol1 + n1 × tol2` where tol1/tol2 are
per-feature tolerances already computed for q1/q2.

Changes to `annotate_feature_()`:
- Read multiplier from ChargePair
- Pass to `toAdductString(formula, charge, multiplier)`
- Set `mol_multiplier` metadata on feature when multiplier > 1

### 2. MassExplainer

New method:
```cpp
SignedSize queryMultimer(
    Int net_charge,
    double m1,              // mz1 × |q1|
    double m2,              // mz2 × |q2|
    Int n1, Int n2,         // multiplier combo (n1 < n2)
    double tolerance,       // absolute mass tolerance
    float thresh_log_p,
    std::vector<std::vector<Compomer>::const_iterator>& hits) const;
```

Algorithm:
1. Binary search to find compomers with matching net_charge
2. Linear scan within that group
3. For each compomer, compute left_mass and right_mass via `getSideMass()`
4. Check: `|n2×m1 - n1×m2 - (n2×left_mass - n1×right_mass)| < tolerance`
5. Collect matching compomer iterators into hits vector

Returns non-contiguous matches (unlike existing query which returns a contiguous range),
hence the vector-of-iterators return type.

### 3. ChargePair

Add two fields:
- `Int mol_multiplier_left_` (default 1)
- `Int mol_multiplier_right_` (default 1)

With getter/setter accessors. Constructor gets optional parameters with default 1 so all
existing call sites are unaffected.

### 4. Compomer

New method:
```cpp
double getSideMass(UInt side) const;
```

Computes `sum(amount × singleMass)` for all adducts on the given side. Purely derived from
existing data, no new stored fields.

### 5. Adduct::toAdductString

Add a 3-parameter overload (or static method):
```cpp
static String toAdductString(const String& ion_string, const Int& charge, UInt mol_multiplier);
```

When `mol_multiplier > 1`, prepends the number before M: `[2M+H]+`.
The existing 2-parameter method stays unchanged (backward compatible, assumes multiplier=1).

## Components NOT Changed

- **Adduct** — no `mol_multiplier_` field, no hash/equality changes
- **Compomer** — no structural changes (just the derived `getSideMass()` helper)
- **MassExplainer::compute()** — same compomer table generation
- **FeatureDeconvolution** (protein decharger) — multimers are a metabolomics concern
- **ILP solver** — cross-multiplier edges participate like any other edge
- **pyOpenMS bindings** — follow-up if needed

## Runtime Impact

With `max_multimer=1` (default): zero impact, identical code path.

With `max_multimer=2`: one extra (n1=1, n2=2) combo per feature pair per charge combo.
Cross-multiplier query does linear scan of net_charge group (~50-200 compomers) vs binary
search (~10 comparisons). Estimated ~16x slowdown on the inner loop, but the outer loop
(feature pair iteration) remains the bottleneck. Acceptable for an opt-in feature.

## Testing

1. **ChargePair_test** — multiplier fields: constructor, getters/setters, defaults, equality
2. **MassExplainer_test** — `queryMultimer()`: monomer-dimer match, different-adduct cross-match,
   no-match, filter by probability threshold
3. **Compomer_test** — `getSideMass()`: empty, single adduct, multiple adducts
4. **Adduct_test** — `toAdductString` with multiplier: `[M+H]+`, `[2M+H]+`, `[3M+Na]+`, `[2M-H]-`
5. **MetaboliteFeatureDeconvolution_test** — integration: feature map with monomer + dimer at
   correct m/z, verify grouping and annotation with `max_multimer=2`. Verify `max_multimer=1`
   produces identical results to current behavior (regression test).
