# GlyPy reference fixtures

These fixtures exercise `TheoreticalGlycanSpectrumGenerator_test` without Python
packages or network access at test time. Expected masses come from GlyPy, not
from OpenMS. The C++ test code and `regenerate.py` are BSD-3-Clause; the adapted
GlyPy data retain Apache-2.0. `$Maintainer: Timo Sachsenberg $`

## Provenance and license

Source: [GlyPy](https://github.com/mobiusklein/glypy), revision
`8d129a8c950e8635165cda9b4b1af392d4e7289e`.

- `stored_masses.tsv` is adapted from
  [`test_data/fragments-example.json`](https://github.com/mobiusklein/glypy/blob/8d129a8c950e8635165cda9b4b1af392d4e7289e/test_data/fragments-example.json).
  It retains the 96 supported glycosidic entries, their numerical masses and
  duplicate multiplicities. Unsupported kinds were removed and JSON was
  converted to TSV; masses were not regenerated.
- `structures.json` contains the unchanged GlycoCT strings for `common_glycan`,
  `branchy_glycan` and `G58143RL`, extracted from
  [`tests/common.py`](https://github.com/mobiusklein/glypy/blob/8d129a8c950e8635165cda9b4b1af392d4e7289e/tests/common.py).
  This selects three entries and changes the container from Python to JSON.
- `structural.tsv` is a derived snapshot generated with **GlyPy 1.0.17** from
  those structures. It maps nodes to parent-first zero-based indices, maps
  stereochemical residue names to OpenMS composition symbols, represents the
  sulfated residue by its formula, and records neutral fragment masses to
  11 decimal places. No OpenMS results are used to create this file.
- `LICENSE` and `NOTICE` are verbatim copies of GlyPy's
  [Apache-2.0 license](https://github.com/mobiusklein/glypy/blob/8d129a8c950e8635165cda9b4b1af392d4e7289e/LICENSE)
  and [upstream notice](https://github.com/mobiusklein/glypy/blob/8d129a8c950e8635165cda9b4b1af392d4e7289e/NOTICE).
  Modified data files also identify the changes in their headers or metadata.

The GlyPy authors and contributors retain their rights in the upstream data.
The upstream NOTICE credits the Pyteomics-derived composition implementation and
cites Goloborodko et al. (2013),
[doi:10.1007/s13361-012-0516-6](https://doi.org/10.1007/s13361-012-0516-6).

The charge conversion uses the proton mass `1.00727646677` Da from
[Pyteomics 5.0.1 `nist_mass`](https://pyteomics.readthedocs.io/en/latest/api/auxiliary.html#pyteomics.auxiliary.constants.nist_mass).
See also Levitsky et al. (2019), *Pyteomics 4.0*,
[doi:10.1021/acs.jproteome.8b00717](https://doi.org/10.1021/acs.jproteome.8b00717).
No Pyteomics source or glycopeptidepy fixtures are copied here.

## Coverage and conventions

The historical fixture checks 96 neutral masses, grouped by fragmentation kind
and compared as multisets. The structural snapshot has 524 unique cleavage
interpretations across the three glycans. Filtering it to cleavage limits 1, 2
and 3 produces 964 interpretations, tested at charges 1, 2 and 3: **2,892 charged
comparisons**, each checking both neutral mass and m/z. The historical comparison
overlaps the `common_glycan` structural case.

The structures cover internal fragments, distinct isobaric branches, fucose,
Neu5Ac and a sulfated HexNAc formula. Matching includes B/C ions with optional Y
branch cuts, and reducing-end Y-only or Z-only cuts. Cross-ring and mixed Y/Z
products are excluded. OpenMS's virtual attachment bond below node zero has no
GlyPy counterpart and is excluded from this comparison. Cleavage identity,
duplicate interpretations, missing fragments and extra fragments are checked.

The absolute tolerance is **0.0001 Da**, with relative tolerance disabled,
matching GlyPy's stored-fixture test tolerance. This accommodates small
differences in the libraries' atomic mass tables. These are theoretical
fragmentation references, not experimental spectra or intensity predictions.

## Regeneration

In a separate Python environment, install `glypy==1.0.17` and run:

```sh
python regenerate.py
```

The script uses only the three local GlycoCT inputs and GlyPy. It regenerates
`structural.tsv` deterministically and does not read OpenMS output. Keep the
historical `stored_masses.tsv` unchanged unless deliberately updating its pinned
upstream source and provenance.
