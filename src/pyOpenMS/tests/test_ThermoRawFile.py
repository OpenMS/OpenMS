"""Public options and instrument-configuration API, independent of RAW fixtures."""
import pytest
import pyopenms as oms


def test_thermo_options():
    if not hasattr(oms, "ThermoRawFile"):
        pytest.skip("OpenMS built without Thermo RAW support")
    reader = oms.ThermoRawFile()
    options = reader.getOptions()
    assert options.preserve_trailers and options.instrument_methods and options.checksum
    assert not options.centroid
    options.centroid = options.charge_data = options.noise_data = options.all_detectors = True
    reader.setOptions(options)
    assert reader.getOptions().centroid
    assert reader.getOptions().all_detectors


def test_instrument_configurations():
    experiment = oms.MSExperiment()
    instrument = oms.Instrument()
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
