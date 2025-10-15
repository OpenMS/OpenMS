# QC Classes in pyOpenMS

This document describes the QC (Quality Control) metric classes that have been wrapped for use in pyOpenMS.

## Overview

The QC classes provide various quality control metrics for LC-MS data analysis. These classes can be used to:
- Summarize MSExperiments and FeatureMaps
- Calculate quality metrics for peptide identifications
- Annotate features with quality information
- Assess data quality for proteomics experiments

## Available QC Classes

### Base Class

#### QCBase
Abstract base class for all QC metric classes.

**Nested Classes:**
- `SpectraMap` - Map to find a spectrum via its NativeID
- `Requires` - Enum for input file requirements
- `ToleranceUnit` - Enum for tolerance units (AUTO, PPM, DA)

### QC Metric Classes

#### 1. TIC (Total Ion Count)
Calculates the Total Ion Count of an MSExperiment.

**Usage:**
```python
import pyopenms as oms

tic = oms.TIC()
result = tic.compute(exp, bin_size=0.0, ms_level=1)
print(f"Area under TIC: {result.area}")
print(f"Number of jumps: {result.jump}")
print(f"Number of falls: {result.fall}")
```

**Methods:**
- `compute(exp, bin_size, ms_level)` - Compute TIC with optional RT binning
- `getName()` - Returns "TIC"
- `getResults()` - Returns list of Result objects

**Result Structure:**
- `intensities` - TIC intensities
- `relative_intensities` - Relative TIC intensities  
- `retention_times` - RT values in seconds
- `area` - Area under TIC curve
- `fall` - MS1 signal fall (10x) count
- `jump` - MS1 signal jump (10x) count

---

#### 2. FWHM (Full Width at Half Maximum)
Transfers FWHM metavalues from features to their PeptideIdentifications.

**Usage:**
```python
fwhm = oms.FWHM()
fwhm.compute(feature_map)
```

**Methods:**
- `compute(features)` - Moves FWHM metavalues to PeptideIdentifications
- `getName()` - Returns "FWHM"

---

#### 3. PeptideMass
Annotates PeptideHits with their theoretical mass.

**Usage:**
```python
peptide_mass = oms.PeptideMass()
peptide_mass.compute(feature_map)
```

**Methods:**
- `compute(features)` - Sets 'mass' metavalue on all PeptideHits
- `getName()` - Returns "PeptideMass"

---

#### 4. Contaminants
Identifies contaminant peptides in a FeatureMap.

**Usage:**
```python
contaminants = oms.Contaminants()
fasta_entries = []  # Load from FASTAFile
contaminants.compute(feature_map, fasta_entries)
results = contaminants.getResults()
```

**Methods:**
- `compute(features, contaminants)` - Checks peptides against contaminant DB
- `getName()` - Returns "Contaminants"
- `getResults()` - Returns list of ContaminantsSummary objects

**ContaminantsSummary Structure:**
- `assigned_contaminants_ratio` - Ratio in assigned peptides
- `unassigned_contaminants_ratio` - Ratio in unassigned peptides
- `all_contaminants_ratio` - Overall ratio
- `assigned_contaminants_intensity_ratio` - Intensity-based ratio

---

#### 5. Ms2IdentificationRate
Computes the MS2 identification rate (identified PSMs / total MS2 scans).

**Usage:**
```python
ms2_rate = oms.Ms2IdentificationRate()
ms2_rate.compute(feature_map, exp, assume_all_target=False)
results = ms2_rate.getResults()
```

**Methods:**
- `compute(feature_map, exp, assume_all_target)` - Calculate ID rate from FeatureMap
- `compute(pep_ids, exp, assume_all_target)` - Calculate ID rate from PeptideIdentifications
- `getName()` - Returns "Ms2IdentificationRate"
- `getResults()` - Returns list of IdentificationRateData objects

**IdentificationRateData Structure:**
- `num_peptide_identification` - Number of peptide IDs
- `num_ms2_spectra` - Number of MS2 spectra
- `identification_rate` - Calculated rate

---

#### 6. MzCalibration
Annotates (un)calibrated m/z errors in PeptideIdentifications.

**Usage:**
```python
mz_cal = oms.MzCalibration()
spec_map = oms.SpectraMap(exp)
mz_cal.compute(feature_map, exp, spec_map)
```

**Methods:**
- `compute(features, exp, map_to_spectrum)` - Annotate m/z errors as metavalues
- `getName()` - Returns "MzCalibration"

---

#### 7. RTAlignment
Annotates raw and aligned retention times.

**Usage:**
```python
rt_align = oms.RTAlignment()
rt_align.compute(feature_map, transformation_description)
```

**Methods:**
- `compute(feature_map, trafo)` - Annotate RT values in FeatureMap
- `compute(pep_ids, trafo)` - Annotate RT values in PeptideIdentificationList
- `getName()` - Returns "RTAlignment"

---

#### 8. MissedCleavages
Counts missed cleavages per PeptideIdentification.

**Usage:**
```python
mc = oms.MissedCleavages()
search_params = protein_id.getSearchParameters()
mc.compute(feature_map, search_params, max_mc=2)
results = mc.getResults()
```

**Methods:**
- `compute(fmap, search_params, max_mc)` - Count missed cleavages
- `getName()` - Returns "MissedCleavages"
- `getResults()` - Returns list of maps (missed cleavages -> counts)

---

#### 9. FragmentMassError
Calculates fragment mass error for identified spectra.

**Usage:**
```python
fme = oms.FragmentMassError()
spec_map = oms.SpectraMap(exp)
search_params = protein_id.getSearchParameters()
fme.compute(pep_ids, search_params, exp, spec_map, 
            oms.ToleranceUnit.PPM, 20.0)
results = fme.getResults()
```

**Methods:**
- `compute(pep_ids, search_params, exp, map_to_spectrum, tolerance_unit, tolerance)`
- `getName()` - Returns "FragmentMassError"
- `getResults()` - Returns list of Statistics objects

**Statistics Structure:**
- `average_ppm` - Average fragment mass error in ppm
- `variance_ppm` - Variance of fragment mass error

---

#### 10. PSMExplainedIonCurrent
Computes the PSM explained ion current.

**Usage:**
```python
psm = oms.PSMExplainedIonCurrent()
spec_map = oms.SpectraMap(exp)
psm.compute(feature_map, exp, spec_map, 
            oms.ToleranceUnit.PPM, 20.0)
results = psm.getResults()
```

**Methods:**
- `compute(fmap, exp, map_to_spectrum, tolerance_unit, tolerance)` - From FeatureMap
- `compute(pep_ids, search_params, exp, map_to_spectrum, tolerance_unit, tolerance)` - From PeptideIDs
- `getName()` - Returns "PSMExplainedIonCurrent"
- `getResults()` - Returns list of Statistics objects

**Statistics Structure:**
- `average_correctness` - Average PSM explained ion current
- `variance_correctness` - Variance

---

#### 11. Ms2SpectrumStats
Collects MS2 spectrum statistics and scan event numbers.

**Usage:**
```python
ms2_stats = oms.Ms2SpectrumStats()
spec_map = oms.SpectraMap(exp)
unassigned_ids = ms2_stats.compute(exp, feature_map, spec_map)
```

**Methods:**
- `compute(exp, features, map_to_spectrum)` - Collect MS2 statistics
- `getName()` - Returns "Ms2SpectrumStats"

**ScanEvent Structure:**
- `scan_event_number` - Consecutive MS2 scan number
- `ms2_presence` - Boolean for MS2 presence

---

## Common Patterns

### Loading Data
```python
import pyopenms as oms

# Load mzML file
exp = oms.MSExperiment()
oms.MzMLFile().load("input.mzML", exp)

# Load featureXML file
feature_map = oms.FeatureMap()
oms.FeatureXMLFile().load("input.featureXML", feature_map)

# Create spectrum map for efficient lookup
spec_map = oms.SpectraMap(exp)
```

### Computing Multiple Metrics
```python
# TIC
tic = oms.TIC()
tic_result = tic.compute(exp, 0.0, 1)

# FWHM
fwhm = oms.FWHM()
fwhm.compute(feature_map)

# Peptide Mass
pm = oms.PeptideMass()
pm.compute(feature_map)

# MS2 ID Rate
ms2_rate = oms.Ms2IdentificationRate()
ms2_rate.compute(feature_map, exp, False)
ms2_results = ms2_rate.getResults()
```

## Notes

- Most QC classes require specific input data (e.g., MSExperiment, FeatureMap)
- Many metrics annotate data by adding metavalues to PeptideIdentifications or PeptideHits
- Results can be retrieved via `getResults()` methods where applicable
- The `SpectraMap` class provides efficient spectrum lookup by NativeID

## See Also

- OpenMS QualityControl TOPP tool documentation
- FeatureMap and MSExperiment documentation
- PeptideIdentification and ProteinIdentification documentation
