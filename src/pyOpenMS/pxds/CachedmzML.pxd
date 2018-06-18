from MSExperiment  cimport *
from MSSpectrum  cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from String cimport *
from ProgressLogger cimport *
from streampos cimport *

cdef extern from "<OpenMS/FORMAT/CachedMzML.h>" namespace "OpenMS":

    cdef cppclass CachedmzML:

        CachedmzML() nogil except +
        CachedmzML(CachedmzML) nogil except +
        CachedmzML(String filename) nogil except +

        Size getNrSpectra() nogil except +
        Size getNrChromatograms() nogil except +

        MSSpectrum getSpectrum(Size idx) nogil except +
        MSChromatogram getChromatogram(Size idx) nogil except +

        # COMMENT: only retrieves experiment meta data (no actual data in spectra/chromatograms)
        # COMMENT: useful for filtering by attributes to then retrieve data
        MSExperiment getMetaData() nogil except +

# COMMENT: wrap static methods
cdef extern from "<OpenMS/FORMAT/CachedMzML.h>" namespace "OpenMS::CachedmzML":
    
    void store(const String& filename, MSExperiment exp) nogil except + # wrap-attach:CachedmzML
    void load(const String& filename, CachedmzML& exp) nogil except + # wrap-attach:CachedmzML

