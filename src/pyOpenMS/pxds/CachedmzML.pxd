from MSExperiment  cimport *
from MSSpectrum  cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from String cimport *
from ProgressLogger cimport *
from streampos cimport *

cdef extern from "<OpenMS/FORMAT/CachedMzML.h>" namespace "OpenMS":

    cdef cppclass CachedmzML:

        CachedmzML() except + nogil  # wrap-doc:A class that uses on-disk caching to read and write spectra and chromatograms
        CachedmzML(CachedmzML &) except + nogil 

        CachedmzML(String filename) except + nogil 

        Size getNrSpectra() except + nogil 
        Size getNrChromatograms() except + nogil 

        MSSpectrum getSpectrum(Size idx) except + nogil 
        MSChromatogram getChromatogram(Size idx) except + nogil 

        # COMMENT: only retrieves experiment meta data (no actual data in spectra/chromatograms)
        # COMMENT: useful for filtering by attributes to then retrieve data
        MSExperiment getMetaData() except + nogil 

# COMMENT: wrap static methods
cdef extern from "<OpenMS/FORMAT/CachedMzML.h>" namespace "OpenMS::CachedmzML":
    
    void store(const String& filename, MSExperiment exp) except + nogil  # wrap-attach:CachedmzML
    void load(const String& filename, CachedmzML& exp) except + nogil  # wrap-attach:CachedmzML
