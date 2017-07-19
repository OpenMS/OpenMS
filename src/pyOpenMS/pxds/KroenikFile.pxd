from String cimport *

from MSSpectrum cimport *
from Peak1D cimport *
from FeatureMap cimport *
from Feature cimport *

cdef extern from "<OpenMS/FORMAT/KroenikFile.h>" namespace "OpenMS":

    cdef cppclass KroenikFile:

        KroenikFile() nogil except +

        void store(String filename, MSSpectrum & spectrum)  nogil except +
        void load(String filename, FeatureMap & feature_map) nogil except +
