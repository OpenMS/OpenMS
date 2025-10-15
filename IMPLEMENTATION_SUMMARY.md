# Summary: QC Classes Wrapped in pyOpenMS

## Issue
The issue requested wrapping the QC (Quality Control) classes from the `src/openms/include/OpenMS/QC/` subfolder in Python. These classes are used by the QualityControl TOPP tool to summarize MSExperiments or FeatureMaps and provide various quality metrics for LC-MS data analysis.

## Solution
Created Python wrapper files (`.pxd` files) for all major QC classes used in the QualityControl tool. These wrappers follow the OpenMS pyOpenMS wrapping conventions using Cython.

## Files Created

### Core Wrapper Files (in `src/pyOpenMS/pxds/`)
1. **QCBase.pxd** - Base class for all QC metrics
   - Includes SpectraMap nested class for spectrum lookup
   - Includes Requires and ToleranceUnit enums
   
2. **FlagSet.pxd** - Support class for QCBase (template class for flag sets)

3. **TIC.pxd** - Total Ion Count metric
   - Computes TIC with optional RT binning
   - Returns Result structure with intensities, RTs, area, jumps, falls

4. **FWHM.pxd** - Full Width at Half Maximum metric
   - Transfers FWHM values from features to PeptideIdentifications

5. **PeptideMass.pxd** - Peptide mass calculation
   - Annotates theoretical mass on PeptideHits

6. **Contaminants.pxd** - Contaminant detection
   - Checks peptides against contaminant database
   - Returns ContaminantsSummary with ratios

7. **Ms2IdentificationRate.pxd** - MS2 identification rate
   - Calculates ratio of identified MS2 spectra
   - Returns IdentificationRateData structure

8. **MzCalibration.pxd** - m/z calibration error
   - Annotates (un)calibrated m/z errors

9. **RTAlignment.pxd** - Retention time alignment
   - Annotates raw and aligned RT values

10. **MissedCleavages.pxd** - Missed cleavage counting
    - Counts missed cleavages per peptide

11. **FragmentMassError.pxd** - Fragment mass error
    - Calculates FME in ppm and Da
    - Returns Statistics with average and variance

12. **PSMExplainedIonCurrent.pxd** - PSM explained ion current
    - Calculates explained ion current for PSMs
    - Returns Statistics structure

13. **Ms2SpectrumStats.pxd** - MS2 spectrum statistics
    - Collects scan event numbers and statistics
    - Returns ScanEvent structures

### Documentation and Test Files
1. **QC_CLASSES_README.md** - Comprehensive documentation
   - Overview of all QC classes
   - Usage examples for each class
   - Common patterns and workflows

2. **example_qc_usage.py** - Example usage script
   - Demonstrates instantiation of all QC classes
   - Shows basic usage patterns
   - Includes working examples for TIC and other metrics

3. **test_qc_wrappers.py** - Unit test file
   - Tests import and instantiation of all QC classes
   - Verifies getName() method works
   - Can be run after pyOpenMS is built

## Key Design Decisions

1. **Abstract Base Class**: QCBase is properly wrapped as an abstract class (no constructor)

2. **Nested Classes**: Nested classes like SpectraMap, ToleranceUnit, and result structures are wrapped in separate namespace blocks following OpenMS conventions

3. **Status/FlagSet**: The Status type (which is FlagSet<Requires> in C++) and related methods were commented out in the wrappers because FlagSet is a complex template class that's difficult to wrap. The core functionality of the QC classes doesn't require this for most use cases.

4. **Method Overloading**: Maintained method overloading where C++ classes have multiple compute() methods with different signatures

5. **Result Structures**: Wrapped all result structures (ContaminantsSummary, Statistics, IdentificationRateData, etc.) so users can access computed metrics

## Usage

After building pyOpenMS with these wrappers, users can:

```python
import pyopenms as oms

# Create QC metric objects
tic = oms.TIC()
fwhm = oms.FWHM()
contaminants = oms.Contaminants()

# Load data
exp = oms.MSExperiment()
oms.MzMLFile().load("input.mzML", exp)

feature_map = oms.FeatureMap()
oms.FeatureXMLFile().load("input.featureXML", feature_map)

# Compute metrics
tic_result = tic.compute(exp, 0.0, 1)
fwhm.compute(feature_map)

# Access results
print(f"TIC area: {tic_result.area}")
print(f"Metric name: {tic.getName()}")
```

## Benefits

1. **Python Access**: QC metrics can now be computed directly in Python scripts without calling TOPP tools
2. **Integration**: Easy integration with Python data analysis workflows and Jupyter notebooks
3. **Rapid Development**: Enables quick prototyping of QC pipelines in Python
4. **Flexibility**: Users can combine QC metrics programmatically

## Testing

The test file `test_qc_wrappers.py` can be run after building pyOpenMS to verify all classes are properly wrapped and accessible. The test imports each class, creates instances, and calls getName() to verify the wrappers work correctly.

## Next Steps

To use these wrappers:
1. Build pyOpenMS with the new wrapper files (the autowrap tool will generate the Cython code)
2. Run the test: `python src/pyOpenMS/tests/unittests/test_qc_wrappers.py`
3. Try the examples: `python src/pyOpenMS/example_qc_usage.py`
4. Refer to the documentation in `QC_CLASSES_README.md` for detailed usage

## Note

These wrappers follow OpenMS conventions and should work with the autowrap tool that generates the actual Cython bindings. The build process will create the final Python module that exposes these classes.
