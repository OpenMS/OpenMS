# Metabolomics in OpenMS 3.0 - 3.5

## New Tools & Major Features

| Version | Highlights |
|---------|------------|
| **3.0** | FLASHDeconv (ultra-fast deconvolution), TOPPView isotope patterns for metabolites |
| **3.2** | SiriusExport, AssayGeneratorMetaboSirius, IonMobilityBinning |
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
- **arm64 Linux wheels** + Python 3.14 support (3.5)

## Workflow Improvements

- **AssayGeneratorMetabo** split into heuristic + SIRIUS-based workflows (3.2)
- **OpenSwath**: IM scoring, mobilogram peak-picking, automated iRT calibration
- **FLASHDeconv**: FDR estimation, isobaric quantification (3.5)
