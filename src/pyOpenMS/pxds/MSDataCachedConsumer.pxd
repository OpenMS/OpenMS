from MSSpectrum cimport *
from MSChromatogram cimport *
from ExperimentalSettings cimport *

cdef extern from "<OpenMS/FORMAT/DATAACCESS/MSDataCachedConsumer.h>" namespace "OpenMS":

    cdef cppclass MSDataCachedConsumer:
        # wrap-doc:
            #   Transforming and cached writing consumer of MS data
            #   -----
            #   Is able to transform a spectrum on the fly while it is read using a
            #   function pointer that can be set on the object. The spectra is then
            #   cached to disk using the functions provided in CachedMzMLHandler.


        MSDataCachedConsumer(String filename) nogil except +
        MSDataCachedConsumer(String filename, bool clear) nogil except +
        # copy constructor of 'MSDataCachedConsumer' is implicitly deleted because field 'ofs_' has a deleted copy constructor
        MSDataCachedConsumer(MSDataCachedConsumer &) nogil except + # wrap-ignore

        void consumeSpectrum(MSSpectrum & s) nogil except +
            # wrap-doc:
                #   Write a spectrum to the output file
                #   -----
                #   May delete data from spectrum (if clearData is set)

        void consumeChromatogram(MSChromatogram & c) nogil except +
            # wrap-doc:
                #   Write a chromatogram to the output file
                #   -----
                #   May delete data from chromatogram (if clearData is set)

        void setExperimentalSettings(ExperimentalSettings& exp) nogil except +
        void setExpectedSize(Size expectedSpectra, Size expectedChromatograms) nogil except +

