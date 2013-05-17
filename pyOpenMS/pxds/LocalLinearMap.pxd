from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from Matrix cimport *

cdef extern from "<OpenMS/ANALYSIS/PIP/LocalLinearMap.h>" namespace "OpenMS":
    
    cdef cppclass LocalLinearMap "OpenMS::LocalLinearMap":
        LocalLinearMap() nogil except +
        LocalLinearMap(LocalLinearMap) nogil except + #wrap-ignore
        # TODO doesnt work
        # LLMParam getLLMParam() nogil except +
        # TODO wrap Matrix
        # Matrix[ double ]  getCodebooks() nogil except +
        # Matrix[ double ]  getMatrixA() nogil except +
        libcpp_vector[ double ]  getVectorWout() nogil except +
        # Matrix[ UInt ]  getCord() nogil except +
        void normalizeVector(libcpp_vector[ double ] & aaIndexVariables) nogil except +
        # libcpp_vector[ double ] neigh(Matrix[ UInt ] & cord, Size win, DoubleReal radius) nogil except +


cdef extern from "<OpenMS/ANALYSIS/PIP/LocalLinearMap.h>" namespace "OpenMS::LocalLinearMap":
    
    cdef cppclass LLMParam "OpenMS::LocalLinearMap::LLMParam":
        LLMParam() #wrap-ignore
        LLMParam(LLMParam) nogil except + #wrap-ignore
        UInt xdim
        UInt ydim
        DoubleReal radius

