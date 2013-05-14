from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from SemanticValidator cimport *
from ControlledVocabulary cimport *

cdef extern from "<OpenMS/FORMAT/VALIDATORS/MzMLValidator.h>" namespace "OpenMS::Internal":
    
    cdef cppclass Internal_MzMLValidator "OpenMS::Internal::MzMLValidator":
        # Internal_MzMLValidator() nogil except +
        Internal_MzMLValidator(Internal_MzMLValidator) nogil except + #wrap-ignore
        Internal_MzMLValidator(CVMappings &mapping, ControlledVocabulary &cv) nogil except +
        # ~MzMLValidator() nogil except +
