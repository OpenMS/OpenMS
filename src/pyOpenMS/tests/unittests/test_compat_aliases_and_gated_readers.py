"""Regression tests for the nanobind rewrite's backward-compatibility aliases
and the graceful absence of build-gated vendor readers (OpenMS issue #9460, §2).

The Cython->nanobind rewrite ships compatibility aliases (``PeakMap``,
``PeakSpectrum``, ``Kernel_MassTrace``, ``FileType``, the ``Interfaces``
namespace) and turned ``Date`` into its own class (it used to alias
``DateTime``). It also exposes vendor readers only behind build flags:
``BrukerTimsFile`` is present only with ``WITH_OPENTIMS``; ``ThermoRawFile``
is not exposed to Python at all.

These tests pin two contracts the general suite only exercises incidentally:
1. the compat aliases still resolve to the right objects, and
2. ``import pyopenms`` works regardless of build flags and probing a gated
   reader returns a clean ``bool`` (never raises ``AttributeError``).

The deprecation-*warning* contract for the DataFrame shims
(``peptide_identifications_to_df`` / ``update_scores_from_df``) is covered
separately in test_deprecated_dataframe_shims.py; the aliases below are
permanent compatibility aliases and intentionally do not warn.
"""
import pyopenms


# -----------------------------------------------------------------------
# Backward-compatibility alias identities
# -----------------------------------------------------------------------


def test_kernel_type_aliases_resolve():
    """PeakMap/PeakSpectrum/Kernel_MassTrace must alias their new classes."""
    assert pyopenms.PeakMap is pyopenms.MSExperiment
    assert pyopenms.PeakSpectrum is pyopenms.MSSpectrum
    assert pyopenms.Kernel_MassTrace is pyopenms.MassTrace


def test_filetype_alias_resolves():
    """FileType is the convenience alias for the nested FileTypes.Type enum."""
    assert pyopenms.FileType is pyopenms.FileTypes.FileType
    # a couple of representative enum members must be reachable through it
    assert hasattr(pyopenms.FileType, "MZML")
    assert hasattr(pyopenms.FileType, "FEATUREXML")


def test_interfaces_namespace():
    """pyopenms.Interfaces.{Spectrum,Chromatogram} must mirror the top-level types."""
    assert hasattr(pyopenms, "Interfaces")
    assert hasattr(pyopenms.Interfaces, "Spectrum")
    assert hasattr(pyopenms.Interfaces, "Chromatogram")
    assert pyopenms.Interfaces.Spectrum is pyopenms.Spectrum
    assert pyopenms.Interfaces.Chromatogram is pyopenms.Chromatogram


# -----------------------------------------------------------------------
# Date is now its own class (used to alias DateTime in Cython pyOpenMS)
# -----------------------------------------------------------------------


def test_date_is_distinct_class():
    """Date is a proper OpenMS::Date class, not an alias for DateTime."""
    assert pyopenms.Date is not pyopenms.DateTime
    today = pyopenms.Date.today()
    assert type(today) is pyopenms.Date
    assert today.year() >= 2020
    assert 1 <= today.month() <= 12
    assert 1 <= today.day() <= 31


# -----------------------------------------------------------------------
# Graceful absence of build-gated vendor readers
# -----------------------------------------------------------------------


def test_import_pyopenms_succeeds():
    """import pyopenms must succeed regardless of build flags."""
    assert pyopenms is not None
    # a core, always-present class is importable
    assert hasattr(pyopenms, "MSExperiment")


def test_bruker_tims_file_presence_is_clean_bool():
    """hasattr(pyopenms, 'BrukerTimsFile') must be a clean bool, never raise.

    The class is only exposed with WITH_OPENTIMS; probing for it must not
    raise AttributeError on a reader-absent build.
    """
    present = hasattr(pyopenms, "BrukerTimsFile")
    assert isinstance(present, bool)
    if present:
        # when built in, it must be usable (mirrors test_BrukerTimsFile.py)
        cfg = pyopenms.BrukerTimsFile.Config()
        assert cfg is not None
    else:
        # when absent, attribute access raises a clean AttributeError
        import pytest
        with pytest.raises(AttributeError):
            _ = pyopenms.BrukerTimsFile


def test_thermo_raw_file_not_exposed_to_python():
    """ThermoRawFile is intentionally not exposed to Python (issue #9460, §2/§6).

    Probing must still be a clean bool (no AttributeError).
    """
    present = hasattr(pyopenms, "ThermoRawFile")
    assert isinstance(present, bool)
    assert present is False


if __name__ == "__main__":
    test_kernel_type_aliases_resolve()
    test_filetype_alias_resolves()
    test_interfaces_namespace()
    test_date_is_distinct_class()
    test_import_pyopenms_succeeds()
    test_bruker_tims_file_presence_is_clean_bool()
    test_thermo_raw_file_not_exposed_to_python()
    print("All compat-alias / gated-reader tests passed!")
