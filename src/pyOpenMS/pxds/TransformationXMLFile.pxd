from Types cimport *
from Feature cimport *
from FeatureMap cimport *
from String cimport *
from TransformationDescription cimport *

cdef extern from "<OpenMS/FORMAT/TransformationXMLFile.h>" namespace "OpenMS":

    cdef cppclass TransformationXMLFile:
        TransformationXMLFile() nogil except +
        # copy constructor of 'TransformationXMLFile' is implicitly deleted because base class 'Internal::XMLHandler' has a deleted copy constructor protected Internal::XMLHandler
        TransformationXMLFile(TransformationXMLFile &) nogil except + # wrap-ignore

        void load(String, TransformationDescription &, bool fit_model) nogil except +
        void store(String, TransformationDescription) nogil except +
