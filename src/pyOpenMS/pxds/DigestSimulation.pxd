from DefaultParamHandler cimport *
from SimTypes cimport *
from FeatureMap cimport *
from Feature cimport *

cdef extern from "<OpenMS/SIMULATION/DigestSimulation.h>" namespace "OpenMS":
    
    cdef cppclass DigestSimulation "OpenMS::DigestSimulation":
        DigestSimulation() except + nogil 
        DigestSimulation(DigestSimulation &) except + nogil 
        void digest(FeatureMap & feature_map) except + nogil 

