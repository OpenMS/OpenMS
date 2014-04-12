from DefaultParamHandler cimport *
from SimTypes cimport *
from FeatureMap cimport *
from Feature cimport *

cdef extern from "<OpenMS/SIMULATION/DigestSimulation.h>" namespace "OpenMS":
    
    cdef cppclass DigestSimulation "OpenMS::DigestSimulation":
        DigestSimulation() nogil except +
        DigestSimulation(DigestSimulation) nogil except +
        void digest(FeatureMap[Feature] & feature_map) nogil except +

