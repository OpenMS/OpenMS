from MSSpectrum cimport *
from MSChromatogram cimport *
from ExperimentalSettings cimport *

cdef extern from "<OpenMS/FORMAT/DATAACCESS/MSDataCachedConsumer.h>" namespace "OpenMS":

    cdef cppclass MSDataCachedConsumer:
        # wrap-doc:
            #  Transforming and cached writing consumer of MS data
            #  
            #  Is able to transform a spectrum on the fly while it is read using a
            #  function pointer that can be set on the object. The spectra is then
            #  cached to disk using the functions provided in CachedMzMLHandler.


        MSDataCachedConsumer(String filename) except + nogil 
        MSDataCachedConsumer(String filename, bool clear) except + nogil 
        # copy constructor of 'MSDataCachedConsumer' is implicitly deleted because field 'ofs_' has a deleted copy constructor
        MSDataCachedConsumer(MSDataCachedConsumer &) except + nogil  # wrap-ignore

        void consumeSpectrum(MSSpectrum & s) except + nogil 
            # wrap-doc:
                #  Write a spectrum to the output file
                #  
                #  May delete data from spectrum (if clearData is set)

        void consumeChromatogram(MSChromatogram & c) except + nogil 
            # wrap-doc:
                #  Write a chromatogram to the output file
                #  
                #  May delete data from chromatogram (if clearData is set)

        void setExperimentalSettings(ExperimentalSettings& exp) except + nogil 
        void setExpectedSize(Size expectedSpectra, Size expectedChromatograms) except + nogil 

