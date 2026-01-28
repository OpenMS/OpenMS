from Types cimport *
from ConsensusMap cimport *
from TransformationDescription cimport *
from ProgressLogger cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/BaseSuperimposer.h>" namespace "OpenMS":
    
    cdef cppclass BaseSuperimposer(DefaultParamHandler,ProgressLogger) :
        # wrap-ignore
        # ABSTRACT class
        # no-pxd-import
        # wrap-inherits:
        #  DefaultParamHandler
        #  ProgressLogger
        BaseSuperimposer() except + nogil 
        BaseSuperimposer(BaseSuperimposer &) except + nogil  # compiler
 
