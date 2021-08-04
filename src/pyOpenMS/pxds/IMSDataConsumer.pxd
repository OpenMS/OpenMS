from ExperimentalSettings cimport *
from Types cimport *

cdef extern from "<OpenMS/INTERFACES/IMSDataConsumer.h>" namespace "OpenMS::Interfaces":

    cdef cppclass IMSDataConsumer[SpectrumType, ChromatogramType]:
        # wrap-ignore
        # no-pxd-import
        # ABSTRACT class

        void consumeSpectrum(SpectrumType &) nogil except + # wrap-doc:Consume a spectrum. The spectrum will be consumed by the implementation and possibly modified

        void consumeChromatogram(ChromatogramType &) nogil except + # wrap-doc:Consume a chromatogram. The chromatogram will be consumed by the implementation and possibly modified

        void setExpectedSize(Size, Size) nogil except + 
            # wrap-doc:
                #   Set expected size of spectra and chromatograms to be consumed
                #   -----
                #   Some implementations might care about the number of spectra and
                #   chromatograms to be consumed and need to be informed about this
                #   (usually before consuming starts)
                #   -----
                #   :param expectedSpectra: Number of spectra expected
                #   :param expectedChromatograms: Number of chromatograms expected

        void setExperimentalSettings(ExperimentalSettings &) nogil except +
            # wrap-doc:
                #   Set experimental settings (meta-data) of the data to be consumed
                #   -----
                #   Some implementations might need to know about the meta-data (or the
                #   context) of the spectra and chromatograms to be consumed. This method
                #   allows them learn this
                #   -----
                #   :param exp: Experimental settings meta data for the data to be consumed

