from MSSpectrum cimport *
from MSChromatogram cimport *
from IMSDataConsumer cimport *
from ExperimentalSettings cimport *

cdef extern from "python_ms_data_consumer.hpp":

    cdef cppclass PythonMSDataConsumer(IMSDataConsumer[Peak1D, ChromatogramPeak]):
        # wrap-ignore

        PythonMSDataConsumer(object py_consumer,
                             object (*spectrum_wrapper)(const MSSpectrum[Peak1D] &),
                             object (*chromatogram_wrapper)(const MSChromatogram[ChromatogramPeak] &),
                             object (*experimental_settings_wrapper)(const ExperimentalSettings &)
                            )

        void consumeSpectrum(MSSpectrum[Peak1D] &) except +

        void consumeChromatogram(MSChromatogram[ChromatogramPeak] &) except +

        void setExpectedSize(Size ns, Size, nc) except +

        void setExperimentalSettings(ExperimentalSettings &) except +












