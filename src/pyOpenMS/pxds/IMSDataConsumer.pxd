from ExperimentalSettings cimport *
from Types cimport *

cdef extern from "<OpenMS/INTERFACES/IMSDataConsumer.h>" namespace "OpenMS::Interfaces":

    cdef cppclass IMSDataConsumer[SpectrumType, ChromatogramType]:
        # wrap-ignore
        # ABSTRACT class

        void consumeSpectrum(SpectrumType &)

        void consumeChromatogram(ChromatogramType &)

        void setExpectedSize(Size, Size)

        void setExperimentalSettings(ExperimentalSettings &)

