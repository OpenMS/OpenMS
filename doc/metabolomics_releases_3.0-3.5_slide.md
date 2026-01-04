# Metabolomics in OpenMS 3.0 - 3.5

## New Tools & Major Features

| Version | Highlights |
|---------|------------|
| **3.2** | SiriusExport, AssayGeneratorMetaboSirius, IonMobilityBinning |
| **3.3** | FeatureFinderMetabo: `report_smoothed_intensities` parameter |
| **3.5** | PeakPickerIM (TimsTOF), FAIMS support across FeatureFinders |

## Ion Mobility & FAIMS Support (3.5)

- **FeatureFinderMetabo / FeatureFinderMetaboIdent**: FAIMS + Bruker TimsTOF support
- **MassTraceExtractor**: Bruker Ion Mobility support
- **PeakPickerIM**: New tool with 3 methods (mobilogram, clustering, elution profiles)
- **IonMobilityBinning**: Auto-detect FAIMS, split by compensation voltage

## pyOpenMS Enhancements

- **Mobilogram & MobilityPeak1D** bindings (3.5)
- **FeatureFinderMetaboIdent** accessible from Python (3.3)
- **DataFrame API**: `get_df()` with snake_case columns (`'MZ'`→`'mz'`, `'RT'`→`'rt'`)
- **Parquet export** for high-performance data handling (3.5)
- **Drift time accessors** on MSSpectrum (3.5)

## Workflow Improvements

- **AssayGeneratorMetabo** split into heuristic + SIRIUS-based workflows (3.2)
- **MetaboliteAdductDecharger**: Fixed column headers for multi-map output (3.5)
- **mzTab-M**: Fixed validation errors for metabolomics metadata (3.5)
