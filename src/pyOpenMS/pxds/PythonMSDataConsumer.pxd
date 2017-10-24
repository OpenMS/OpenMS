from MSSpectrum cimport *
from MSChromatogram cimport *
from IMSDataConsumer cimport *
from ExperimentalSettings cimport *

# see ../extra_includes/python_ms_data_consumer.hpp for actual wrapped C++ code
cdef extern from "python_ms_data_consumer.hpp":

    cdef cppclass PythonMSDataConsumer(IMSDataConsumer[Peak1D, ChromatogramPeak]):
        # wrap-ignore
        # no-pxd-import

        PythonMSDataConsumer(object py_consumer,
                             object (*spectrum_wrapper)(const MSSpectrum &),
                             object (*chromatogram_wrapper)(const MSChromatogram &),
                             object (*experimental_settings_wrapper)(const ExperimentalSettings &)
                            )

        void consumeSpectrum(MSSpectrum &) except +

        void consumeChromatogram(MSChromatogram &) except +

        void setExpectedSize(Size ns, Size, nc) except +

        void setExperimentalSettings(ExperimentalSettings &) except +












