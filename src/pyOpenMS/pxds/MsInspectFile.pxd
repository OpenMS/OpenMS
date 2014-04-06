from Types cimport *
from String cimport *
from Feature cimport *
from FeatureMap cimport *
from MSSpectrum cimport *

cdef extern from "<OpenMS/FORMAT/MsInspectFile.h>" namespace "OpenMS":
    
    cdef cppclass MsInspectFile "OpenMS::MsInspectFile":
        MsInspectFile() nogil except +
        MsInspectFile(MsInspectFile) nogil except + #wrap-ignore
        void load(String & filename, FeatureMap[Feature] & feature_map) nogil except +
        void store(String & filename, MSSpectrum[Peak1D] & spectrum) nogil except +

