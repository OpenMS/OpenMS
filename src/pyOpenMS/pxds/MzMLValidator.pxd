from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from SemanticValidator cimport *
from ControlledVocabulary cimport *

cdef extern from "<OpenMS/FORMAT/VALIDATORS/MzMLValidator.h>" namespace "OpenMS::Internal":
    
    cdef cppclass Internal_MzMLValidator "OpenMS::Internal::MzMLValidator":
        # private
        Internal_MzMLValidator() except + nogil  # wrap-ignore
        # private
        Internal_MzMLValidator(Internal_MzMLValidator &) except + nogil  # wrap-ignore

        Internal_MzMLValidator(CVMappings &mapping, ControlledVocabulary &cv) except + nogil 
        # ~MzMLValidator() except + nogil 
