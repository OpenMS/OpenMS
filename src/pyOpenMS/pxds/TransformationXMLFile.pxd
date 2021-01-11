from Types cimport *
from Feature cimport *
from FeatureMap cimport *
from String cimport *
from TransformationDescription cimport *

cdef extern from "<OpenMS/FORMAT/TransformationXMLFile.h>" namespace "OpenMS":

    cdef cppclass TransformationXMLFile:
        TransformationXMLFile() nogil except +

        void load(String, TransformationDescription &, bool fit_model) nogil except+
        void store(String, TransformationDescription) nogil except+

