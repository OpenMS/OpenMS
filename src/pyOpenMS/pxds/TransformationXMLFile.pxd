from Types cimport *
from Feature cimport *
from FeatureMap cimport *
from String cimport *
from TransformationDescription cimport *

cdef extern from "<OpenMS/FORMAT/TransformationXMLFile.h>" namespace "OpenMS":

    cdef cppclass TransformationXMLFile:
        TransformationXMLFile() except + nogil 
        # copy constructor of 'TransformationXMLFile' is implicitly deleted because base class 'Internal::XMLHandler' has a deleted copy constructor protected Internal::XMLHandler
        TransformationXMLFile(TransformationXMLFile &) except + nogil  # wrap-ignore

        void load(String, TransformationDescription &, bool fit_model) except + nogil 
        void store(String, TransformationDescription) except + nogil 
