from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool

from String cimport *
from StringList cimport *
from CVMappings cimport *
from ControlledVocabulary cimport *

cdef extern from "<OpenMS/FORMAT/VALIDATORS/SemanticValidator.h>" namespace "OpenMS::Internal":

    cdef cppclass SemanticValidator:

        SemanticValidator(CVMappings & mapping, ControlledVocabulary & cv) nogil except +

        void setCheckTermValueTypes(bool check) nogil except +
        void setCheckUnits(bool check) nogil except +

        bool validate(String & filename, StringList & errors, StringList & warnings) nogil except +

