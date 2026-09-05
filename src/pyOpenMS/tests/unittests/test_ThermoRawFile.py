"""Tests for native Thermo RAW reading through FileHandler.

Reading a .raw file needs three things at run time: OpenMS built with
``WITH_THERMO_RAW=ON``, the openms-thermo-bridge managed assemblies (shipped
in ``pyopenms/share/OpenMS/openms_thermo_bridge/managed``), and a .NET 8 or
newer runtime on the machine. The real-data tests therefore only run when
``PYOPENMS_THERMO_RAW_TEST_FILE`` names a Thermo .raw file; the wheel CI sets
it to PNNL's Angiotensin_AllScans.raw and provides the runtime.

When the variable is set the tests must pass. A wheel that was built without
Thermo support, or that lacks the managed assemblies, is exactly what they are
meant to catch, so there is no feature-detection fallback here.
"""
import os

import pytest
import pyopenms

RAW_FILE = os.environ.get("PYOPENMS_THERMO_RAW_TEST_FILE", "")

needs_raw_file = pytest.mark.skipif(
    not RAW_FILE,
    reason="PYOPENMS_THERMO_RAW_TEST_FILE not set",
)


def test_raw_extension_is_recognized():
    """.raw is typed as RAW regardless of WITH_THERMO_RAW, so tools can hand it
    to an external converter even without the in-process reader."""
    assert pyopenms.FileHandler.getTypeByFileName("run.raw") == pyopenms.FileTypes.RAW
    assert pyopenms.FileHandler.getTypeByFileName("run.RAW") == pyopenms.FileTypes.RAW


def test_managed_bridge_is_bundled_with_share():
    """The managed assemblies travel with share/OpenMS; if the share folder is
    present, the bridge files must be too, otherwise every .raw load fails at
    run time with 'Managed bridge runtime files are missing'."""
    share = os.environ.get("OPENMS_DATA_PATH", "")
    if not share or not os.path.isdir(share):
        pytest.skip("no OpenMS share directory available")
    managed = os.path.join(share, "openms_thermo_bridge", "managed")
    if not os.path.isdir(managed):
        if RAW_FILE:
            pytest.fail("PYOPENMS_THERMO_RAW_TEST_FILE is set but %s is missing" % managed)
        pytest.skip("no bundled Thermo bridge (WITH_THERMO_RAW=OFF or in-tree build)")
    for name in ("ThermoWrapperManaged.dll",
                 "ThermoWrapperManaged.runtimeconfig.json",
                 "ThermoFisher.CommonCore.RawFileReader.dll"):
        assert os.path.isfile(os.path.join(managed, name)), name
    # pyopenms/__init__.py hands this directory to the reader explicitly, so the
    # lookup does not depend on which share/ folder libOpenMS resolves on its own.
    assert os.path.samefile(os.environ.get("OPENMS_THERMO_MANAGED_DIR", ""), managed)


@needs_raw_file
def test_load_experiment_from_raw():
    """FileHandler dispatches .raw to ThermoRawFile: spectra, RTs, precursors
    and instrument metadata are populated from Angiotensin_AllScans.raw."""
    assert os.path.isfile(RAW_FILE), RAW_FILE

    exp = pyopenms.MSExperiment()
    pyopenms.FileHandler().loadExperiment(RAW_FILE, exp)

    # Angiotensin_AllScans.raw: 87 MS1 scans plus hundreds of MS2 scans.
    assert exp.getNrSpectra() > 1000
    assert len(exp.getSourceFiles()) == 1
    assert exp.getSourceFiles()[0].getNameOfFile() == os.path.basename(RAW_FILE)

    levels = [s.getMSLevel() for s in exp]
    assert levels.count(1) > 0
    assert sum(1 for level in levels if level >= 2) > 0

    rts = [s.getRT() for s in exp]
    assert rts == sorted(rts)

    ms1 = next(s for s in exp if s.getMSLevel() == 1)
    mz, intensity = ms1.get_peaks()
    assert len(mz) > 0
    assert len(mz) == len(intensity)

    ms2_with_precursor = [s for s in exp if s.getMSLevel() >= 2 and len(s.getPrecursors()) > 0]
    assert ms2_with_precursor
    assert ms2_with_precursor[0].getPrecursors()[0].getMZ() > 0.0

    # Instrument metadata is filled from the RAW header (Orbitrap, ESI).
    instrument = exp.getInstrument()
    assert len(instrument.getMassAnalyzers()) > 0
    assert len(instrument.getIonSources()) > 0


@needs_raw_file
def test_load_experiment_from_raw_roundtrips_through_mzml(tmp_path):
    """An experiment read from RAW survives a store/load round trip through mzML."""
    exp = pyopenms.MSExperiment()
    pyopenms.FileHandler().loadExperiment(RAW_FILE, exp)

    out = str(tmp_path / "roundtrip.mzML")
    pyopenms.MzMLFile().store(out, exp)

    reloaded = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(out, reloaded)
    assert reloaded.getNrSpectra() == exp.getNrSpectra()
    assert reloaded[0].getMSLevel() == exp[0].getMSLevel()
    assert reloaded[0].getRT() == pytest.approx(exp[0].getRT())
