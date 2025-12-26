# OpenMS Test Coverage Analysis Report

## Executive Summary

This report analyzes test coverage across the OpenMS codebase and identifies areas requiring improved testing. The analysis covers C++ unit tests, Python bindings tests, TOPP integration tests, and GUI components.

### Key Metrics

| Category | Coverage | Status |
|----------|----------|--------|
| C++ Unit Tests (class_tests) | ~44% (599/1,372 files) | ⚠️ Needs Improvement |
| GUI Components | 3.4% (5/147 classes) | 🔴 Critical Gap |
| pyOpenMS Bindings | ~31% (145/476 classes documented) | 🔴 Critical Gap |
| TOPP Integration Tests | 70% (103/148 tools) | 🟡 Moderate |
| OpenSWATH Algorithms | 23% (3/13 classes) | 🔴 Critical Gap |
| NUXL Module | 15% (2/13 classes) | 🔴 Critical Gap |

---

## Priority 1: Critical Gaps (Immediate Attention Required)

### 1.1 GUI Module - 3.4% Coverage

**Location:** `src/openms_gui/include/OpenMS/VISUAL/`
**Impact:** High - User-facing components completely untested

The GUI module has **142 of 147 classes without unit tests**. This is the most significant coverage gap.

**Untested Categories:**

| Category | Count | Examples |
|----------|-------|----------|
| Dialog Classes | 20 | DataFilterDialog, FeatureEditDialog, SaveImageDialog |
| Visualizer Classes | 27 | AcquisitionVisualizer, PeptideHitVisualizer, SpectrumSettingsVisualizer |
| Plot/Canvas Classes | 12 | Plot1DCanvas, Plot2DCanvas, Plot3DCanvas, PlotWidget |
| Layer Classes | 11 | LayerData1DPeak, LayerDataFeature, LayerDataConsensus |
| Annotation Classes | 7 | Annotation1DPeakItem, Annotation1DTextItem |
| Tab/Tree Classes | 9 | SpectraTreeTab, DIATreeTab, LayerListView |
| TOPPAS Classes | 12 | TOPPASEdge, TOPPASToolVertex, TOPPASScene |

**Recommendation:**
- Implement Qt Test framework for GUI testing
- Start with critical user-facing dialogs (SaveImageDialog, DataFilterDialog)
- Add canvas/painter tests for rendering validation

---

### 1.2 NUXL Module - 15% Coverage

**Location:** `src/openms/include/OpenMS/ANALYSIS/NUXL/`
**Impact:** High - Nuclear cross-linking analysis functionality untested

**Only 2 of 13 classes have tests:**
- ✅ NuXLModificationsGenerator
- ✅ NuXLParameterParsing

**Missing Tests (11 classes):**
- NuXLAnnotateAndLocate
- NuXLAnnotatedHit
- NuXLConstants
- NuXLDeisotoper
- NuXLFDR
- NuXLFragmentAdductDefinition
- NuXLFragmentAnnotationHelper
- NuXLFragmentIonGenerator
- NuXLMarkerIonExtractor
- NuXLPresets
- NuXLReport

**Recommendation:** Create comprehensive test suite for NUXL module as it's a specialized analysis feature.

---

### 1.3 pyOpenMS Bindings - 31% Coverage

**Location:** `src/pyOpenMS/`
**Impact:** High - Python users are a significant portion of the user base

**Issues Identified:**

1. **Architecture Mismatch:** Two incompatible test systems
   - Legacy `test000.py` with `@tests` annotations (tracked)
   - Modern `test_*.py` files (26 files, NOT tracked by coverage checker)

2. **Coverage Checker Broken:** `check_test_coverage.py` uses Python 2 syntax

3. **Critical Untested Classes:**
   - File I/O: `MzMLFile`, `MzXMLFile`, `MzDataFile`
   - Algorithms: `ConsensusIDAlgorithm` (+ 8 subclasses)
   - Feature Finding: `FeatureFinderIdentificationAlgorithm`
   - Protein Inference: `BasicProteinInferenceAlgorithm`, `BayesianProteinInferenceAlgorithm`
   - Data Processing: `Deisotoper`
   - Databases: `ElementDB`, `ResidueDB`

**Recommendation:**
- Migrate `check_test_coverage.py` to Python 3
- Unify test tracking across both test styles
- Prioritize testing of file I/O classes (most commonly used)

---

### 1.4 OpenSWATH Algorithms - 23% Coverage

**Location:** `src/openswathalgo/include/OpenMS/OPENSWATHALGO/`
**Impact:** Medium-High - Core DIA analysis functionality

**Only 3 of 13 classes tested:**
- ✅ Scoring
- ✅ DataStructures
- ✅ SwathMap

**Missing Tests (10 classes):**
- DATAACCESS/DataFrameWriter
- DATAACCESS/ISpectrumAccess
- DATAACCESS/ITransition
- DATAACCESS/MockObjects
- DATAACCESS/SpectrumHelpers
- DATAACCESS/TransitionExperiment
- DATAACCESS/TransitionHelper
- DATAACCESS/Transitions
- ALGO/StatsHelpers
- Macros

**Recommendation:** Add tests for transition and spectrum access classes which are fundamental to OpenSWATH workflows.

---

## Priority 2: Moderate Gaps (Should Address Soon)

### 2.1 CONCEPT Module - 54% Coverage

**Location:** `src/openms/include/OpenMS/CONCEPT/`

**Missing 11 of 24 tests:**
- CommonEnums
- Constants
- EnumHelpers
- Helpers
- Init
- Macros
- MacrosTest
- PrecisionWrapper
- ProgressLogger
- Qt5Port
- RAIICleanup

**Impact:** These are foundational framework classes used throughout the codebase.

---

### 2.2 SYSTEM Module - 58% Coverage

**Location:** `src/openms/include/OpenMS/SYSTEM/`

**Missing 5 of 12 tests:**
- BuildInfo
- NetworkGetRequest
- RWrapper
- SIMDe
- UpdateCheck

---

### 2.3 ID Submodule (ANALYSIS) - 76% Coverage

**Location:** `src/openms/include/OpenMS/ANALYSIS/ID/`

**Missing 10 of 43 tests:**
- ConsensusIDAlgorithm (base class)
- ConsensusIDAlgorithmIdentity
- ConsensusIDAlgorithmSimilarity
- IDScoreGetterSetter
- IonIdentityMolecularNetworking
- MessagePasserFactory
- OpenSearchModificationAnalysis
- PeptideSearchEngineFIAlgorithm
- SiriusExportAlgorithm
- SiriusMSConverter

---

### 2.4 OPENSWATH Submodule (ANALYSIS) - 81% Coverage

**Location:** `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/`

**Missing 8 of 43 tests:**
- DATAACCESS/MRMFeatureAccessOpenMS
- DATAACCESS/SimpleOpenMSSpectraAccessFactory
- DATAACCESS/SpectrumAccessOpenMS
- DATAACCESS/SpectrumAccessOpenMSCached
- DATAACCESS/SpectrumAccessOpenMSInMemory
- DATAACCESS/SpectrumAccessTransforming
- OpenSwathOSWWriter
- OpenSwathWorkflow

---

### 2.5 TOPP Tools - 70% Coverage

**45 of 148 TOPP tools lack integration tests:**

**External Adapters (9):** CometAdapter, LuciphorAdapter, MSFraggerAdapter, MSGFPlusAdapter, MascotAdapterOnline, NovorAdapter, PercolatorAdapter, SageAdapter, SpectraSTSearchAdapter

**QC Tools (7):** QCEmbedder, QCExporter, QCExtractor, QCImporter, QCMerger, QCShrinker, QuantmsIOConverter

**Spectral Filters (4):** SpectraFilterNLargest, SpectraFilterNormalizer, SpectraFilterThresholdMower, SpectraMerger

**Other (25):** Various clustering, workflow, and analysis tools

---

## Priority 3: Easy Wins (Quick Improvements)

### 3.1 Single Missing Tests (96%+ Coverage Modules)

| Module | Missing Class | Location |
|--------|---------------|----------|
| ML/REGRESSION | LinearRegressionWithoutIntercept | `src/openms/include/OpenMS/ML/REGRESSION/` |
| MAPMATCHING | FeatureMapping | `src/openms/include/OpenMS/ANALYSIS/MAPMATCHING/` |

### 3.2 Degenerate Case Tests Underutilized

**Available edge case test files (15):**
- Empty files: `empty.fasta`, `empty.idXML`, `empty.mzML`
- Invalid values: `neg_int.mzML`, `neg_mz.mzML`, `neg_rt.mzML`
- Invalid parameters: `PPHiRes_invalidParamName.ini`, `PPHiRes_invalidValue.ini`
- Edge cases: `one_scan.mzML`, `empty_scan.mzML`

**Current usage:** Only 8 tests use these files.

**Recommendation:** Expand degenerate case testing to improve robustness.

---

## Priority 4: Well-Covered Modules (Maintain Quality)

The following modules have excellent coverage (90%+):

| Module | Coverage | Notes |
|--------|----------|-------|
| PROCESSING | 100% (31/31) | All signal processing classes tested |
| FEATUREFINDER | 93.5% (43/46) | Core feature detection well-tested |
| DATASTRUCTURES | 93% (42/45) | Core data structures tested |
| KERNEL | 88% (29/33) | Core spectrum/feature classes tested |
| CHEMISTRY | 82% (41/50) | Chemical utilities well-tested |
| FORMAT | 78% (108/139) | File I/O extensively tested |
| ML | 96% (25/26) | Machine learning algorithms tested |

---

## Recommendations Summary

### Immediate Actions (Priority 1)

1. **GUI Module:** Establish Qt Test framework and create tests for critical dialog classes
2. **NUXL Module:** Write comprehensive test suite for all 11 untested classes
3. **pyOpenMS:** Fix `check_test_coverage.py` for Python 3 and unify test tracking
4. **OpenSWATH Algo:** Add tests for transition and spectrum access classes

### Short-Term Actions (Priority 2)

5. **CONCEPT/SYSTEM:** Add tests for foundational framework classes
6. **ID Submodule:** Test ConsensusID algorithm variants and Sirius integration
7. **TOPP Tools:** Add basic integration tests for QC tools and spectral filters

### Quick Wins (Priority 3)

8. **ML/REGRESSION:** Add LinearRegressionWithoutIntercept test
9. **MAPMATCHING:** Add FeatureMapping test
10. **Edge Cases:** Expand use of degenerate case test files

### Process Improvements

11. **Test Coverage CI:** Add coverage reporting to CI pipeline with thresholds
12. **New Code Policy:** Require tests for all new classes (enforce in PR reviews)
13. **Documentation:** Update developer docs with test coverage expectations

---

## Test Infrastructure Notes

### Existing Frameworks
- **C++ Tests:** OpenMS ClassTest framework (custom macros)
- **Python Tests:** pytest/unittest
- **Integration Tests:** CTest with FuzzyDiff for output comparison
- **Coverage Tools:** gcov/lcov (enabled via `OPENMS_COVERAGE=ON`)

### Test Data Location
- Class tests: `src/tests/class_tests/openms/data/`
- TOPP tests: `src/tests/topp/`
- Degenerate cases: `src/tests/class_tests/openms/data/degenerate_cases/`

### Running Coverage Analysis
```bash
cmake -DOPENMS_COVERAGE=ON -DCMAKE_BUILD_TYPE=Debug ..
make
ctest
make OpenMS_coverage
```

---

## Appendix: Complete List of Untested Classes

### GUI Module (142 classes)
<details>
<summary>Click to expand</summary>

**ANNOTATION (7):** Annotation1DCaret, Annotation1DDistanceItem, Annotation1DItem, Annotation1DPeakItem, Annotation1DTextItem, Annotation1DVerticalLineItem, Annotations1DContainer

**APPLICATION (6):** FLASHDeconvWizardBase, INIFileEditorWindow, QApplicationTOPP, SwathWizardBase, TOPPASBase, TOPPViewBase

**DIALOG (20):** DataFilterDialog, FLASHDeconvTabWidget, FeatureEditDialog, HistogramDialog, LayerStatisticsDialog, ListFilterDialog, Plot1DGoToDialog, Plot1DPrefDialog, Plot2DGoToDialog, Plot2DPrefDialog, Plot3DPrefDialog, PythonModuleRequirement, PythonSelector, SaveImageDialog, SpectrumAlignmentDialog, SwathTabWidget, TOPPASIOMappingDialog, TOPPASInputFileDialog, TOPPASInputFilesDialog, TOPPASOutputFilesDialog, TOPPASToolConfigDialog, TOPPASVertexNameDialog, TOPPViewOpenDialog, TOPPViewPrefDialog, ToolsDialog, WizardHelper

**VISUALIZER (27):** AcquisitionInfoVisualizer, AcquisitionVisualizer, BaseVisualizer, BaseVisualizerGUI, ContactPersonVisualizer, DataProcessingVisualizer, DocumentIdentifierVisualizer, ExperimentalSettingsVisualizer, GradientVisualizer, HPLCVisualizer, InstrumentSettingsVisualizer, InstrumentVisualizer, IonDetectorVisualizer, IonSourceVisualizer, MassAnalyzerVisualizer, MetaInfoDescriptionVisualizer, MetaInfoVisualizer, PeptideHitVisualizer, PeptideIdentificationVisualizer, PrecursorVisualizer, ProductVisualizer, ProteinHitVisualizer, ProteinIdentificationVisualizer, SampleVisualizer, ScanWindowVisualizer, SoftwareVisualizer, SourceFileVisualizer, SpectrumSettingsVisualizer

**PLOT/CANVAS (12):** Plot1DCanvas, Plot1DWidget, Plot2DCanvas, Plot2DWidget, Plot3DCanvas, Plot3DOpenGLCanvas, Plot3DWidget, PlotCanvas, PlotWidget, Painter1DBase, Painter2DBase, PainterBase

**LAYER (11):** LayerData1DBase, LayerData1DChrom, LayerData1DIonMobility, LayerData1DPeak, LayerDataBase, LayerDataChrom, LayerDataConsensus, LayerDataFeature, LayerDataIdent, LayerDataIonMobility, LayerDataPeak

**TAB/TREE (9):** DIATreeTab, SpectraIDViewTab, SpectraTreeTab, EnhancedTabBar, EnhancedTabBarWidgetInterface, LayerListView, TOPPASTreeView, TableView, TreeView

**OTHER (50):** AxisPainter, AxisWidget, ColorSelector, EnhancedWorkspace, ExternalProcessMBox, FilterList, FilterableList, GUIProgressLoggerImpl, HistogramWidget, InputFile, InputFileList, ListEditor, LogWindow, MetaDataBrowser, MultiGradientSelector, OutputDirectory, ParamEditor, RecentFilesMenu, SequenceVisualizer, SwathLibraryStats, TOPPASEdge, TOPPASInputFileListVertex, TOPPASMergerVertex, TOPPASOutputFileListVertex, TOPPASOutputFolderVertex, TOPPASOutputVertex, TOPPASResource, TOPPASResources, TOPPASScene, TOPPASSplitterVertex, TOPPASToolVertex, TOPPASWidget, TOPPViewMenu, TVControllerBase, TVDIATreeTabController, TVIdentificationViewController, TVSpectraViewController, TVToolDiscovery, IPeptideIds, CommonDefs, and more...

</details>

### NUXL Module (11 classes)
NuXLAnnotateAndLocate, NuXLAnnotatedHit, NuXLConstants, NuXLDeisotoper, NuXLFDR, NuXLFragmentAdductDefinition, NuXLFragmentAnnotationHelper, NuXLFragmentIonGenerator, NuXLMarkerIonExtractor, NuXLPresets, NuXLReport

### CONCEPT Module (11 classes)
CommonEnums, Constants, EnumHelpers, Helpers, Init, Macros, MacrosTest, PrecisionWrapper, ProgressLogger, Qt5Port, RAIICleanup

### SYSTEM Module (5 classes)
BuildInfo, NetworkGetRequest, RWrapper, SIMDe, UpdateCheck

---

*Report generated: December 2024*
*OpenMS Version: Based on current develop branch*
