from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from ChargePair cimport *
from FeatureMap cimport *
from Feature cimport *

# ctypedef std::vector<ChargePair> PairsType;
cdef extern from "<OpenMS/ANALYSIS/DECHARGING/ILPDCWrapper.h>" namespace "OpenMS":
    
    cdef cppclass ILPDCWrapper "OpenMS::ILPDCWrapper":
        ILPDCWrapper() except + nogil 
        ILPDCWrapper(ILPDCWrapper &) except + nogil  # compiler
        double compute(FeatureMap fm, libcpp_vector[ChargePair] & pairs, Size verbose_level) except + nogil  # wrap-doc:Compute optimal solution and return value of objective function. If the input feature map is empty, a warning is issued and -1 is returned

