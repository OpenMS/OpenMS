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
def test_load_experiment_from_raw_applies_vendor_peak_picking():
    """FileHandler applies Thermo's peak picking, like FileConverter; ThermoRawFile
    keeps the scans as acquired unless its centroid option is set."""
    exp = pyopenms.MSExperiment()
    pyopenms.FileHandler().loadExperiment(RAW_FILE, exp)
    centroid = pyopenms.SpectrumSettings.SpectrumType.CENTROID
    assert all(s.getType() == centroid for s in exp)

    # Angiotensin_AllScans.raw has profile MS1 scans and centroided MS2 scans.
    acquired = pyopenms.MSExperiment()
    pyopenms.ThermoRawFile().load(RAW_FILE, acquired)
    assert [s.getNativeID() for s in acquired] == [s.getNativeID() for s in exp]
    profile = pyopenms.SpectrumSettings.SpectrumType.PROFILE
    ms1 = [(a, c) for a, c in zip(acquired, exp) if a.getMSLevel() == 1]
    assert ms1 and all(a.getType() == profile for a, _ in ms1)
    assert sum(c.size() for _, c in ms1) < sum(a.size() for a, _ in ms1)


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


# ---------------------------------------------------------------------------
# Public options and instrument-configuration API; independent of RAW fixtures.
# ---------------------------------------------------------------------------

def test_thermo_options():
    if not hasattr(pyopenms, "ThermoRawFile"):
        pytest.skip("OpenMS built without Thermo RAW support")
    reader = pyopenms.ThermoRawFile()
    options = reader.getOptions()
    assert options.preserve_trailers and options.instrument_methods and options.checksum
    assert not options.centroid
    options.centroid = options.charge_data = options.noise_data = options.all_detectors = True
    reader.setOptions(options)
    assert reader.getOptions().centroid
    assert reader.getOptions().all_detectors
    # ProgressLogger methods are bound explicitly
    reader.setLogType(pyopenms.LogType.NONE)
    assert reader.getLogType() == pyopenms.LogType.NONE


def test_thermo_selection_options():
    if not hasattr(pyopenms, "ThermoRawFile"):
        pytest.skip("OpenMS built without Thermo RAW support")
    options = pyopenms.ThermoRawFileOptions()
    assert options.first_scan == 0 and options.last_scan == 0
    assert options.rt_min == float("-inf") and options.rt_max == float("inf")
    assert options.precursor_mz_min == float("-inf") and options.precursor_mz_max == float("inf")
    assert options.ms_levels == []
    assert not options.metadata_only and not options.skip_chromatograms
    options.ms_levels = [2, 3]
    options.first_scan, options.last_scan = 10, 20
    options.rt_min, options.rt_max = 1.0, 2.0
    options.metadata_only = options.skip_chromatograms = True
    copy = pyopenms.ThermoRawFileOptions(options)
    assert copy.ms_levels == [2, 3]
    assert (copy.first_scan, copy.last_scan) == (10, 20)
    assert (copy.rt_min, copy.rt_max) == (1.0, 2.0)
    assert copy.metadata_only and copy.skip_chromatograms


@needs_raw_file
def test_raw_selection_matches_full_load():
    """Selected scans are the same as in a full load, precursors included."""
    options = pyopenms.ThermoRawFileOptions()
    options.centroid = True
    options.checksum = False
    reader = pyopenms.ThermoRawFile()
    reader.setOptions(options)
    full = pyopenms.MSExperiment()
    reader.load(RAW_FILE, full)

    options.ms_levels = [2]
    options.skip_chromatograms = True
    reader.setOptions(options)
    ms2 = pyopenms.MSExperiment()
    reader.load(RAW_FILE, ms2)
    expected = [s for s in full if s.getMSLevel() == 2]
    assert ms2.getNrSpectra() == len(expected) > 0
    assert ms2.getNrChromatograms() == 0
    for got, want in zip(ms2, expected):
        assert got.getNativeID() == want.getNativeID()
        assert got.getPrecursors()[0].getMZ() == want.getPrecursors()[0].getMZ()
        assert got.getPrecursors()[0].getIntensity() == want.getPrecursors()[0].getIntensity()
        assert got.size() == want.size()

    options.ms_levels = []
    options.metadata_only = True
    reader.setOptions(options)
    headers = pyopenms.MSExperiment()
    reader.load(RAW_FILE, headers)
    assert headers.getNrSpectra() == full.getNrSpectra()
    assert all(s.size() == 0 for s in headers)


def test_instrument_configurations():
    experiment = pyopenms.MSExperiment()
    instrument = pyopenms.Instrument()
    instrument.setName("Orbitrap Astral")
    experiment.setInstrumentConfigurations({"astral": instrument})
    configurations = experiment.getInstrumentConfigurations()
    assert configurations["astral"].getName() == "Orbitrap Astral"

    instrument.setName("changed input")
    assert experiment.getInstrumentConfigurations()["astral"].getName() == "Orbitrap Astral"
    configurations["astral"].setName("changed copy")
    assert experiment.getInstrumentConfigurations()["astral"].getName() == "Orbitrap Astral"
    with pytest.raises(TypeError):
        experiment.setInstrumentConfigurations({"invalid": None})
    experiment.setInstrumentConfigurations({})
    assert experiment.getInstrumentConfigurations() == {}
