"""Python bindings to the OpenMS C++ library.

The pyOpenMS package contains Python bindings for a large part of the OpenMS
library (https://openms.de) for mass spectrometry based proteomics. It thus
provides providing facile access to a feature-rich, open-source algorithm
library for mass-spectrometry based proteomics analysis. These Python bindings
allow raw access to the data-structures and algorithms implemented in OpenMS,
specifically those for file access (mzXML, mzML, TraML, mzIdentML among
others), basic signal processing (smoothing, filtering, de-isotoping and
peak-picking) and complex data analysis (including label-free, SILAC, iTRAQ and
SWATH analysis tools).

For further documentation, please see https://pyopenms.readthedocs.io

Please cite:

    Röst HL, Schmitt U, Aebersold R, Malmström L.
    pyOpenMS: a Python-based interface to the OpenMS mass-spectrometry algorithm library.
    Proteomics. 2014 Jan;14(1):74-7. doi: 10.1002/pmic.201300246.

"""
from __future__ import annotations

import os
import sys
import warnings

from ._sysinfo import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)

try:
    import importlib.metadata
    __version__ = importlib.metadata.version("pyopenms")
except Exception:
    __version__ = "0+unknown"

here = os.path.abspath(os.path.dirname(__file__))

# Set up OpenMS data path if not already set
default_openms_data_path = os.path.join(here, "share", "OpenMS")
env_openms_data_path = os.environ.get("OPENMS_DATA_PATH")

if os.path.exists(default_openms_data_path):
    if not env_openms_data_path:
        os.environ["OPENMS_DATA_PATH"] = default_openms_data_path
    elif os.path.abspath(env_openms_data_path) != os.path.abspath(default_openms_data_path):
        warnings.warn(
            "Warning: OPENMS_DATA_PATH environment exists and points to a different location "
            "than the default share directory. "
            "pyOpenMS will use it ({env}) to locate data in the OpenMS share folder "
            "(e.g., the unimod database), instead of the default ({default})."
            .format(env=env_openms_data_path, default=default_openms_data_path)
        )
else:
    if not env_openms_data_path:
        # Try to find the data path relative to this package
        from pathlib import Path
        _possible_paths = [
            Path(__file__).parent.parent / "share" / "OpenMS",
            Path(sys.prefix) / "share" / "OpenMS",
        ]
        for _path in _possible_paths:
            if _path.exists():
                os.environ["OPENMS_DATA_PATH"] = str(_path)
                break
        else:
            warnings.warn(
                "Warning: OPENMS_DATA_PATH environment variable not found and no share directory was installed. "
                "Some functionality might not work as expected."
            )
        del _possible_paths, _path, Path


# Patch pyarrow to handle filesystem registration conflicts with OpenMS Arrow C++.
# Both pyarrow and C++ Arrow try to register
# filesystem factories for the 'file' scheme, causing ArrowKeyError.
# This patch caches LocalFileSystem instances created before C++ loads and reuses them.
# See: https://github.com/apache/arrow/issues/44696#issuecomment-3532192904
def _patch_pyarrow_filesystem():
    try:
        import pyarrow.fs as _pa_fs
        from pyarrow._fs import FileSystem, LocalFileSystem
        from pyarrow.lib import ArrowKeyError
    except ImportError:
        return  # pyarrow not installed, no conflict possible

    # Verify required internal APIs exist (may change between pyarrow versions)
    if not all(hasattr(_pa_fs, attr) for attr in
               ('_resolve_filesystem_and_path', '_is_path_like', '_stringify_path', '_ensure_filesystem')):
        return  # Incompatible pyarrow version, skip patching

    # Create cached LocalFileSystem instances BEFORE C++ Arrow loads
    # This claims the 'file://' scheme registration for pyarrow
    try:
        _cached_local_fs = {
            False: LocalFileSystem(use_mmap=False),
            True: LocalFileSystem(use_mmap=True),
        }
    except ArrowKeyError:
        # C++ Arrow already registered - too late to cache, patching won't help
        return

    def _patched_resolve(path, filesystem=None, *, memory_map=False):
        """Patched version that reuses cached LocalFileSystem instances."""
        # Handle file-like objects
        if not _pa_fs._is_path_like(path):
            if filesystem is not None:
                raise ValueError(
                    "'filesystem' passed but the specified path is file-like, so"
                    " there is nothing to open with 'filesystem'."
                )
            return filesystem, path

        # Handle explicit filesystem
        if filesystem is not None:
            filesystem = _pa_fs._ensure_filesystem(filesystem, use_mmap=memory_map)
            if isinstance(filesystem, LocalFileSystem):
                path = _pa_fs._stringify_path(path)
            elif not isinstance(path, str):
                raise TypeError(
                    "Expected string path; path-like objects are only allowed "
                    "with a local filesystem"
                )
            path = filesystem.normalize_path(path)
            return filesystem, path

        path = _pa_fs._stringify_path(path)

        # PATCH: Reuse cached LocalFileSystem instead of creating new one
        # Cache is populated before C++ Arrow loads, so this should always hit
        filesystem = _cached_local_fs.get(memory_map)
        if filesystem is None:
            # Cache miss - create and cache (may fail if C++ already registered)
            filesystem = LocalFileSystem(use_mmap=memory_map)
            _cached_local_fs[memory_map] = filesystem

        try:
            file_info = filesystem.get_file_info(path)
        except ValueError:
            file_info = None
            exists_locally = False
        else:
            exists_locally = file_info.type != _pa_fs.FileType.NotFound

        # If file doesn't exist locally, try parsing as URI
        if not exists_locally:
            try:
                filesystem, path = FileSystem.from_uri(path)
            except ValueError as e:
                if "empty scheme" in str(e) or "Cannot parse URI" in str(e):
                    pass  # Neither URI nor local path - will propagate file not found
                else:
                    raise
        else:
            path = filesystem.normalize_path(path)

        return filesystem, path

    _pa_fs._resolve_filesystem_and_path = _patched_resolve

_patch_pyarrow_filesystem()
del _patch_pyarrow_filesystem

# On conda the libs will be installed to the general conda lib path which is available during load.
# Try to skip this loading if we do not ship the libraries in the package (e.g. as wheel via pip).
if sys.platform.startswith("linux") and os.path.exists(os.path.join(here, "libOpenMS.so")):
    # load local shared libraries before we import pyopenms*.so, else
    # those are not found. setting LD_LIBRARY_PATH does not work,
    # see: http://stackoverflow.com/questions/1178094
    import ctypes
    ctypes.cdll.LoadLibrary(os.path.join(here, "libOpenSwathAlgo.so"))
    ctypes.cdll.LoadLibrary(os.path.join(here, "libOpenMS.so"))


def _import_submodules():
    """Import all nanobind submodules and merge into this namespace."""
    import importlib
    import pkgutil

    _imported_modules = []

    # Priority loading order for nanobind shared-domain type resolution.
    #
    # Hard constraint (crash on violation):
    #   nb::class_<Derived, Base> requires Base to be registered first.
    #   Currently the ONLY cross-module inheritance is:
    #     MSExperiment (bind_experiment.cpp) -> ExperimentalSettings (bind_kernel.cpp)
    #   so _pyopenms_kernel MUST load before _pyopenms_experiment.
    #   All other nb::class_ inheritance pairs live in the same binding file.
    #
    # Soft constraint (no crash, but silent None returns):
    #   Functions accepting/returning types from other modules resolve lazily
    #   at call time via nb_type_c2p(). Since import pyopenms loads ALL modules
    #   before any user code runs, order among the remaining modules doesn't
    #   matter in practice.
    #
    # The full list below is kept as defensive documentation of the logical
    # dependency chain. Only the kernel-before-experiment pair is strictly
    # required; the rest could be removed without breakage.
    _priority_modules = [
        "_pyopenms_datastructures",
        "_pyopenms_kernel",           # MUST load before experiment (ExperimentalSettings base)
        "_pyopenms_spectrum",
        "_pyopenms_chromatogram",
        "_pyopenms_experiment",       # requires kernel (ExperimentalSettings)
        "_pyopenms_metadata",
        "_pyopenms_chemistry",
    ]

    all_names = sorted(
        name for _, name, _ in pkgutil.iter_modules(__path__)
        if name.startswith("_pyopenms_")
    )
    ordered = [n for n in _priority_modules if n in all_names]
    ordered += [n for n in all_names if n not in _priority_modules]

    for name in ordered:
        try:
            mod = importlib.import_module(f".{name}", package=__name__)
            _imported_modules.append(mod)
            for attr in dir(mod):
                if not attr.startswith("_") or attr.startswith("__static_"):
                    globals()[attr] = getattr(mod, attr)
        except Exception as e:
            import traceback
            warnings.warn(f"Failed to import {name}: {type(e).__name__}: {e}")
            traceback.print_exc()

    if not _imported_modules:
        raise ImportError(
            "Failed to import any pyOpenMS bindings.\n"
            "Make sure the package was built correctly with nanobind."
        )


# Import bindings
_import_submodules()

# Import pure Python modules
from ._dataframes import DataFrameMixin
from ._python_extras import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)

# Apply Python addons to classes
from .addons import apply_addons

apply_addons(globals())

# Add nested enum aliases for backwards compatibility
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
    "BaseFeature": ["AnnotationState"],
    "BSpline2d": ["BoundaryCondition"],
    "ReactionMonitoringTransition": ["DecoyTransitionType"],
    "Peak2D": ["DimensionDescription"],
    "IDMapper": ["Measure"],
    "ConsensusMapNormalizerAlgorithmMedian": ["NormalizationMethod"],
    "IDRipper": ["OriginAnnotationFormat"],
    "SourceFile": ["ChecksumType"],
}

for class_name, enum_names in _MODULE_LEVEL_NESTED_ENUMS.items():
    if class_name in globals():
        cls = globals()[class_name]
        for enum_name in enum_names:
            if hasattr(cls, enum_name) and enum_name not in globals():
                globals()[enum_name] = getattr(cls, enum_name)

del _MODULE_LEVEL_NESTED_ENUMS
# Clean up loop variables to prevent leaking into stubs
for _cleanup_var in ('class_name', 'enum_names', 'cls', 'enum_name'):
    globals().pop(_cleanup_var, None)
del _cleanup_var

# Add backward-compatible aliases
if "MSExperiment" in globals():
    globals()["PeakMap"] = globals()["MSExperiment"]
if "MSSpectrum" in globals():
    globals()["PeakSpectrum"] = globals()["MSSpectrum"]
if "MassTrace" in globals():
    globals()["Kernel_MassTrace"] = globals()["MassTrace"]
# Date: Now a proper OpenMS::Date class with day(), month(), year(), today()
# In old Cython pyOpenMS, Date was an alias for DateTime; now it's its own class.
# MZTrafoModel_MODELTYPE alias (Cython used ClassName_EnumName convention for some nested enums)
if "MZTrafoModel" in globals() and hasattr(globals()["MZTrafoModel"], "MODELTYPE"):
    globals()["MZTrafoModel_MODELTYPE"] = globals()["MZTrafoModel"].MODELTYPE
# FileTypes.FileType nested enum is exported as "FileType" in old pyOpenMS for convenience
if "FileTypes" in globals() and hasattr(globals()["FileTypes"], "FileType"):
    globals()["FileType"] = globals()["FileTypes"].FileType
# TargetedExperiment_Instrument: Now a separate nanobind class (no longer an alias for Instrument)
# Both "Instrument" (OpenMS::Instrument from METADATA) and "TargetedExperiment_Instrument"
# (TargetedExperimentHelper::Instrument) are exported natively by the kernel module.

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
del os, here, sys
