from Types cimport *
from ConsensusMap cimport *
from TransformationDescription cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/BaseSuperimposer.h>" namespace "OpenMS":
    
    cdef cppclass BaseSuperimposer(DefaultParamHandler) :
        # wrap-ignore
        # ABSTRACT class
        # no-pxd-import
        # wrap-inherits:
        #  DefaultParamHandler
        BaseSuperimposer() except + nogil 
        BaseSuperimposer(BaseSuperimposer &) except + nogil  # compiler
 
