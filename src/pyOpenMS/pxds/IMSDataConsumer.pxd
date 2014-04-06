from ExperimentalSettings cimport *
from Types cimport *

cdef extern from "<OpenMS/FORMAT/MzMLFile.h>" namespace "OpenMS::Interfaces":

    cdef cppclass IMSDataConsumer[SpectrumType, ChromatogramType]:
        # wrap-ignore

        void consumeSpectrum(SpectrumType &)

        void consumeChromatogram(ChromatogramType &)

        void setExpectedSize(Size, Size)

        void setExperimentalSettings(ExperimentalSettings &)

