# IsobaricWorkflow Test Data Requirements

This document describes the test data needed for comprehensive testing of the IsobaricWorkflow tool, specifically for validating the id_merge_index metadata tracking when processing multiple mzML/idXML file pairs.

## Background

When IsobaricWorkflow processes multiple mzML/idXML file pairs, it needs to:
1. Set `id_merge_index` on each PeptideIdentification to track the source file (0, 1, 2, ... for N files)
2. Ensure consensus map column indices are consecutive across all files
3. Enable proper mzTab export with correct spectrum references

## Required Test Files

### Input Files

1. **IsobaricWorkflow_input_1.mzML** (file index 0)
   - TMT or iTRAQ labeled MS2 (or MS3 for SPS-MS3) spectra
   - Should contain ~5-10 MS2/MS3 spectra with reporter ions
   - Peak picked (centroided) data
   - Example: TMT 10-plex with reporter ions at 126-131 m/z range

2. **IsobaricWorkflow_input_1.idXML** (file index 0)
   - Peptide identifications matching spectra in input_1.mzML
   - SpectrumReference must match the NativeID in the mzML file
   - Example format: `spectrum=<index>` or `controllerType=0 controllerNumber=1 scan=<scan>`

3. **IsobaricWorkflow_input_2.mzML** (file index 1)
   - Similar structure to input_1.mzML
   - Different spectra (different scan numbers/native IDs)
   - Same labeling method (e.g., TMT 10-plex)

4. **IsobaricWorkflow_input_2.idXML** (file index 1)
   - Peptide identifications for spectra in input_2.mzML
   - Must have matching spectrum references

5. **IsobaricWorkflow.ini** (optional)
   - Tool parameters for the test
   - Should specify the quantitation method (type parameter)

### Expected Output Files

1. **IsobaricWorkflow_output.consensusXML**
   - Expected consensus map output
   - Should contain:
     - Column headers for all channels from both files (consecutive map indices)
     - ConsensusFeatures with peptide identifications
     - Each PeptideIdentification with `id_merge_index` UserParam (0 for file 1, 1 for file 2)
     - Map indices in ConsensusFeature elements should be consecutive

2. **IsobaricWorkflow_output.mzTab**
   - Expected mzTab output
   - Should have correct ms_run references
   - PSM section should reference correct spectrum IDs

## Creating Test Data

### Option 1: Use Real Data (Recommended)
1. Start with a small TMT/iTRAQ dataset (e.g., from ProteomeXchange)
2. Extract a subset of scans (~10 MS2/MS3 spectra per file)
3. Run identification search to create idXML files
4. Manually verify the data is correct

### Option 2: Create Synthetic Data
1. Create minimal mzML files with:
   - MS1 spectra (precursors)
   - MS2/MS3 spectra with TMT/iTRAQ reporter ions at expected m/z values
   - Proper native IDs
2. Create matching idXML with:
   - ProteinIdentification section
   - PeptideIdentification entries matching spectrum native IDs
   - Reasonable scores

### Tools for Data Creation
- `FileFilter` to extract scan ranges from real data
- `MSExperiment` Python bindings (pyOpenMS) for programmatic creation
- `FileConverter` for format conversions

## Validation

Use the validation script to check output:
```bash
python3 IsobaricWorkflow_validate_merge_index.py IsobaricWorkflow_output.consensusXML 2
```

This checks:
- All PeptideIdentifications have id_merge_index
- Values are consecutive (0, 1, ..., N-1)
- Map indices in ConsensusFeatures are consecutive

## Test Execution

Once test data is available, add to CMakeLists.txt:

```cmake
add_test("TOPP_IsobaricWorkflow_1" ${TOPP_BIN_PATH}/IsobaricWorkflow -test 
  -in ${DATA_DIR_TOPP}/IsobaricWorkflow_input_1.mzML ${DATA_DIR_TOPP}/IsobaricWorkflow_input_2.mzML
  -in_id ${DATA_DIR_TOPP}/IsobaricWorkflow_input_1.idXML ${DATA_DIR_TOPP}/IsobaricWorkflow_input_2.idXML
  -type tmt10plex
  -out IsobaricWorkflow_output.tmp.consensusXML
  -out_mzTab IsobaricWorkflow_output.tmp.mzTab)

add_test("TOPP_IsobaricWorkflow_1_out1" ${DIFF} 
  -whitelist "id=" "<map" "?xml-stylesheet" "uniqueId" 
  -in1 IsobaricWorkflow_output.tmp.consensusXML 
  -in2 ${DATA_DIR_TOPP}/IsobaricWorkflow_output.consensusXML)
set_tests_properties("TOPP_IsobaricWorkflow_1_out1" PROPERTIES DEPENDS "TOPP_IsobaricWorkflow_1")

add_test("TOPP_IsobaricWorkflow_1_out2" ${DIFF} -whitelist "COM"
  -in1 IsobaricWorkflow_output.tmp.mzTab 
  -in2 ${DATA_DIR_TOPP}/IsobaricWorkflow_output.mzTab)
set_tests_properties("TOPP_IsobaricWorkflow_1_out2" PROPERTIES DEPENDS "TOPP_IsobaricWorkflow_1")
```

## Expected Behavior

For 2 input file pairs with TMT 10-plex:
- 20 columns in consensus map (10 channels × 2 files)
- Column map indices: 0-9 (file 0), 10-19 (file 1)
- All PeptideIdentifications from file 0 have id_merge_index=0
- All PeptideIdentifications from file 1 have id_merge_index=1
- mzTab export succeeds without errors
- PSM rows in mzTab correctly reference spectra from both files
