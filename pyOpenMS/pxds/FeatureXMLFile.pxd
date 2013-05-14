from Feature cimport *
from FeatureMap cimport *
from String cimport *
from FeatureFileOptions cimport *

cdef extern from "<OpenMS/FORMAT/FeatureXMLFile.h>" namespace "OpenMS":

    cdef cppclass FeatureXMLFile:
        FeatureXMLFile() nogil except +
        # cython does not support free template args, so Peak1D has
        # to be used as a fixed argument
        void load(String, FeatureMap[Feature] &) nogil except+
        void store(String, FeatureMap[Feature] &) nogil except+

        FeatureFileOptions getOptions() nogil except +
        void setOptions(FeatureFileOptions) nogil except +

        Size loadSize(String path) nogil except +

