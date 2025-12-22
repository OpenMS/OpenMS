# Breaking API Changes in pyOpenMS 3.5

This document summarizes the breaking API changes identified in pyOpenMS 3.5 and the solutions implemented.

## 1. File I/O Classes: PeptideIdentificationList Requirement

### Breaking Change
In pyOpenMS 3.5, file I/O methods changed from accepting Python `list[PeptideIdentification]` to requiring `PeptideIdentificationList`.

### Affected Classes
- `IdXMLFile` - **FIXED**
- `PepXMLFile` - **FIXED**
- `MzIdentMLFile` - **FIXED**  
- `OMSSAXMLFile` - Not fixed (less common format)
- `OMSSACSVFile` - Not fixed (less common format)
- `XTandemXMLFile` - Not fixed (less common format)
- `MascotXMLFile` - Not fixed (less common format)
- `InspectOutfile` - Not fixed (less common format)
- `MSPFile` - Not fixed (less common format)
- `PercolatorOutfile` - Not fixed (less common format)

### Old API (pre-3.5)
```python
protein_ids = []
peptide_ids = []
oms.IdXMLFile().load(filename, protein_ids, peptide_ids)
```

### New API (3.5+)
```python
protein_ids = []
peptide_ids = oms.PeptideIdentificationList()
oms.IdXMLFile().load(filename, protein_ids, peptide_ids)
```

### Solution Implemented
Added backward-compatible wrapper methods in addon files for the three most common file formats:
1. `IdXMLFile` - Most common format for OpenMS peptide identifications
2. `PepXMLFile` - Common format from search engines like Mascot, X!Tandem
3. `MzIdentMLFile` - Standard mzIdentML format

The wrappers:
- Detect whether `peptide_ids` is a Python list or PeptideIdentificationList
- If it's a list, convert to PeptideIdentificationList, call C++ method, then convert back
- If it's PeptideIdentificationList, call C++ method directly

Both old and new styles now work for these three formats.

## 2. DataFrame Column Names Standardized

### Breaking Change
DataFrame column names were standardized to lowercase snake_case for PEP 8 compliance.

### Affected Methods
- `MSSpectrum.get_df()`
- `MSChromatogram.get_df()`
- `ConsensusMap.get_metadata_df()`
- `FeatureMap.get_df()`

### Changes
- `'RT'` → `'rt'`
- `'MZ'` → `'mz'`
- `'nativeID'` → `'native_id'`
- `'precursorMZ'` → `'precursor_mz'`
- And other similar changes

### Solution
**No code changes needed** - this is a documentation-only breaking change. Users need to update their code that references old column names.

### Migration Example
```python
# Old code
df = spectrum.get_df()
rt_values = df['RT']
mz_values = df['MZ']

# New code
df = spectrum.get_df()
rt_values = df['rt']
mz_values = df['mz']
```

## 3. Rank-Related Methods Removed

### Breaking Change
Removed `assignRanks`, `sortByRanks` methods in PeptideIdentifications, and `updateHitRanks` in IDFilter.

### Impact
Code that uses these methods will fail with AttributeError.

### Solution
**No fix implemented** - these methods were removed intentionally. Users should:
- Sort and filter by score instead
- Use rank as a meta value for backwards compatibility if needed

### Migration
```python
# Old code
peptide_ids.assignRanks()
peptide_ids.sortByRanks()

# New code
# Sort by score directly
hits = peptide_id.getHits()
hits.sort(key=lambda x: x.getScore(), reverse=True)
peptide_id.setHits(hits)
```

## 4. PeptideHit Rank Member Removed

### Breaking Change
The `rank` member in `PeptideHit` was removed and replaced with a meta value.

### Solution
**No fix needed** - rank is now stored as a meta value for backwards compatibility.

### Migration
```python
# Old code
rank = hit.getRank()
hit.setRank(1)

# New code
rank = hit.getMetaValue(b"rank") if hit.metaValueExists(b"rank") else None
hit.setMetaValue(b"rank", 1)
```

## Summary of Fixes

### Implemented
✅ `IdXMLFile.load()` and `.store()` - backward compatible with Python lists
✅ `PepXMLFile.load()` and `.store()` - backward compatible with Python lists  
✅ `MzIdentMLFile.load()` and `.store()` - backward compatible with Python lists
✅ Test cases added in `test_backward_compatibility.py`

### Not Implemented (Out of Scope)
- Other less common file I/O classes (can be added if users request)
- DataFrame column name changes (documentation-only, no code fix possible)
- Rank-related method removals (intentional API change)
- PeptideHit rank member change (already has backward compatibility via meta values)

## Testing

To test the changes, build pyOpenMS and run:
```bash
cd src/pyOpenMS
pytest tests/unittests/test_backward_compatibility.py
```

## Recommendations for Users

1. **Immediate action**: Update file I/O code to use the examples shown above
2. **DataFrame users**: Update column name references in your code
3. **Rank users**: Migrate away from rank methods, use score-based sorting
4. **Version check**: Use this pattern for conditional code:
   ```python
   if hasattr(oms, "PeptideIdentificationList"):
       # pyOpenMS 3.5+ code
       peptide_ids = oms.PeptideIdentificationList()
   else:
       # pre-3.5 code
       peptide_ids = []
   ```
