from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from Matrix cimport *

cdef extern from "<OpenMS/ANALYSIS/PIP/LocalLinearMap.h>" namespace "OpenMS":
    
    cdef cppclass LocalLinearMap "OpenMS::LocalLinearMap":
        LocalLinearMap() nogil except +
        LocalLinearMap(LocalLinearMap) nogil except + #wrap-ignore
        LLMParam getLLMParam() nogil except +
        Matrix[ double ]  getCodebooks() nogil except +
        Matrix[ double ]  getMatrixA() nogil except +
        libcpp_vector[ double ]  getVectorWout() nogil except +
        # TODO STL attributes unsigned int 
        # Matrix[ UInt ]  getCord() nogil except +
        void normalizeVector(libcpp_vector[ double ] & aaIndexVariables) nogil except +
        # libcpp_vector[ double ] neigh(Matrix[ unsigned int ] & cord, Size win, DoubleReal radius) nogil except +


cdef extern from "<OpenMS/ANALYSIS/PIP/LocalLinearMap.h>" namespace "OpenMS::LocalLinearMap":
    
    cdef cppclass LLMParam "OpenMS::LocalLinearMap::LLMParam":
        LLMParam() nogil except +
        LLMParam(LLMParam) nogil except + #wrap-ignore
        UInt xdim
        UInt ydim
        DoubleReal radius

