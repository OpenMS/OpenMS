from String cimport *

from MSSpectrum cimport *
from Peak1D cimport *
from FeatureMap cimport *
from Feature cimport *

cdef extern from "<OpenMS/FORMAT/KroenikFile.h>" namespace "OpenMS":

    cdef cppclass KroenikFile:

        KroenikFile() nogil except +

        void store(String filename, MSSpectrum[Peak1D] & spectrum) 
        void load(String filename, FeatureMap[Feature] & feature_map)
