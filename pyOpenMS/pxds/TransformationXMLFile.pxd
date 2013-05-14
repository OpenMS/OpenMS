from Feature cimport *
from FeatureMap cimport *
from String cimport *
from TransformationDescription cimport *

cdef extern from "<OpenMS/FORMAT/TransformationXMLFile.h>" namespace "OpenMS":

    cdef cppclass TransformationXMLFile:
        TransformationXMLFile() nogil except +
        # cython does not support free template args, so Peak1D has
        # to be used as a fixed argument
        void load(String, TransformationDescription &) nogil except+
        void store(String, TransformationDescription) nogil except+
