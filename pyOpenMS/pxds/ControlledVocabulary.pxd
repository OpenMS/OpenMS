from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool

from String cimport *
from StringList cimport *
from CVMappings cimport *

cdef extern from "<OpenMS/FORMAT/ControlledVocabulary.h>" namespace "OpenMS":

    cdef cppclass ControlledVocabulary:

        ControlledVocabulary() nogil except +
        ControlledVocabulary(ControlledVocabulary & voc) nogil except +

        void loadFromOBO(String& name, String& filename) nogil except +

