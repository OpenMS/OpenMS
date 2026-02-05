"""
pyOpenMS2: Python bindings for OpenMS using nanobind

This package provides Python bindings for the OpenMS C++ library,
enabling access to mass spectrometry data structures and algorithms.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

__version__ = "3.6.0"

# Set up OpenMS data path if not already set
if "OPENMS_DATA_PATH" not in os.environ:
    # Try to find the data path relative to this package
    _possible_paths = [
        Path(__file__).parent / "share" / "OpenMS",
        Path(__file__).parent.parent / "share" / "OpenMS",
        Path(sys.prefix) / "share" / "OpenMS",
    ]
    for _path in _possible_paths:
        if _path.exists():
            os.environ["OPENMS_DATA_PATH"] = str(_path)
            break


def _import_submodules():
    """Import all nanobind submodules and merge into this namespace."""
    import importlib
    import pkgutil

    _imported_modules = []
    # Discover all _pyopenms2_* modules (domain-based naming)
    for finder, name, ispkg in pkgutil.iter_modules(__path__):
        if name.startswith("_pyopenms2_"):
            try:
                mod = importlib.import_module(f".{name}", package=__name__)
                _imported_modules.append(mod)
                for attr in dir(mod):
                    if not attr.startswith("_"):
                        globals()[attr] = getattr(mod, attr)
            except Exception as e:
                import warnings, traceback
                warnings.warn(f"Failed to import {name}: {type(e).__name__}: {e}")
                traceback.print_exc()

    if not _imported_modules:
        raise ImportError(
            "Failed to import any pyOpenMS2 bindings.\n"
            "Make sure the package was built correctly with nanobind."
        )


# Import bindings
_import_submodules()

# Import pure Python modules
from ._dataframes import DataFrameMixin

# Apply Python addons to classes
from .addons import apply_addons

apply_addons(globals())

# Add nested enum aliases for backwards compatibility with pyOpenMS (autowrap)
# These enums are defined inside C++ classes but exposed at module level in nanobind
# For compatibility, also add them as class attributes
_NESTED_ENUM_ALIASES = {
    "SpectrumSettings": ["SpectrumType"],
    "ChromatogramSettings": ["ChromatogramType"],
    "IonSource": ["Polarity", "InletType", "IonizationMethod"],
    "SourceFile": ["ChecksumType"],
    "ProteinIdentification": ["PeakMassType"],
    "Instrument": ["IonOpticsType"],
    "InstrumentSettings": ["ScanMode"],
    "IonDetector": ["Type", "AcquisitionMode"],
    "PercolatorOutfile": ["ScoreType"],
    "MassAnalyzer": ["AnalyzerType", "ResolutionMethod", "ResolutionType", "ScanDirection", "ScanLaw", "ReflectronState"],
    "Sample": ["SampleState"],
    "Precursor": ["ActivationMethod"],
    "MZTrafoModel": ["MODELTYPE"],
    "MultipleTesting": ["Pi0Method", "LfdrTransform"],
    "RankData": ["Method", "NaNPolicy"],
    "MSNumpressCoder": ["NumpressCompression"],
    "DataProcessing": ["ProcessingAction"],
    "FileTypes": ["Type"],
    "PeakGroup": ["TargetDecoyType"],
}

for class_name, enum_names in _NESTED_ENUM_ALIASES.items():
    if class_name in globals():
        cls = globals()[class_name]
        for enum_name in enum_names:
            if enum_name in globals() and not hasattr(cls, enum_name):
                try:
                    setattr(cls, enum_name, globals()[enum_name])
                except (TypeError, AttributeError):
                    pass  # Some classes may not allow attribute assignment

del _NESTED_ENUM_ALIASES

# Export nested enums at module level for backward compatibility
# These are enums defined inside classes that old pyOpenMS exposed at module level
_MODULE_LEVEL_NESTED_ENUMS = {
    "Scores": ["IDType"],
    "ProgressLogger": ["LogType"],
}

for class_name, enum_names in _MODULE_LEVEL_NESTED_ENUMS.items():
    if class_name in globals():
        cls = globals()[class_name]
        for enum_name in enum_names:
            if hasattr(cls, enum_name) and enum_name not in globals():
                globals()[enum_name] = getattr(cls, enum_name)

del _MODULE_LEVEL_NESTED_ENUMS

# Add backward-compatible aliases
if "MSExperiment" in globals():
    globals()["PeakMap"] = globals()["MSExperiment"]
if "MSSpectrum" in globals():
    globals()["PeakSpectrum"] = globals()["MSSpectrum"]
# FileTypes.FileType nested enum is exported as "FileType" in old pyOpenMS for convenience
if "FileTypes" in globals() and hasattr(globals()["FileTypes"], "FileType"):
    globals()["FileType"] = globals()["FileTypes"].FileType

# Create Interfaces namespace for pyopenms.Interfaces.Spectrum / Chromatogram
class _InterfacesNamespace:
    """Namespace for OpenMS::Interfaces data structures."""
    pass

_interfaces = _InterfacesNamespace()
if "Spectrum" in globals():
    _interfaces.Spectrum = globals()["Spectrum"]
if "Chromatogram" in globals():
    _interfaces.Chromatogram = globals()["Chromatogram"]
globals()["Interfaces"] = _interfaces
del _interfaces, _InterfacesNamespace

# Add module-level DataFrame utility functions (backward compatibility)
from ._dataframes_compat import peptide_identifications_to_df, update_scores_from_df

# Clean up namespace
del _import_submodules, apply_addons, DataFrameMixin
