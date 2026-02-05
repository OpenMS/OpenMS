#pragma once

#include <nanobind/nanobind.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>

namespace nb = nanobind;

/// Duck-typing wrapper: bridges a Python consumer object (with consumeSpectrum,
/// consumeChromatogram, setExpectedSize, setExperimentalSettings methods) to the
/// C++ IMSDataConsumer interface. This enables MzMLFile.transform() / MzXMLFile.transform()
/// to call back into Python.
class NanobindMSDataConsumer : public OpenMS::Interfaces::IMSDataConsumer {
    nb::object py_consumer_;
public:
    explicit NanobindMSDataConsumer(nb::object consumer) : py_consumer_(std::move(consumer)) {}

    void consumeSpectrum(SpectrumType& s) override {
        nb::gil_scoped_acquire gil;
        py_consumer_.attr("consumeSpectrum")(nb::cast(s, nb::rv_policy::reference));
    }

    void consumeChromatogram(ChromatogramType& c) override {
        nb::gil_scoped_acquire gil;
        py_consumer_.attr("consumeChromatogram")(nb::cast(c, nb::rv_policy::reference));
    }

    void setExpectedSize(size_t expectedSpectra, size_t expectedChromatograms) override {
        nb::gil_scoped_acquire gil;
        py_consumer_.attr("setExpectedSize")(expectedSpectra, expectedChromatograms);
    }

    void setExperimentalSettings(const OpenMS::ExperimentalSettings& exp) override {
        nb::gil_scoped_acquire gil;
        py_consumer_.attr("setExperimentalSettings")(nb::cast(exp, nb::rv_policy::copy));
    }
};
