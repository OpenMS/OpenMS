from libcpp.vector cimport vector as libcpp_vector
from Param cimport *
from String cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/DefaultParamHandler.h>" namespace "OpenMS":

    cdef cppclass DefaultParamHandler:
        #wrap-ignore

        # DefaultParamHandler(String & name) nogil except +
        # DefaultParamHandler(DefaultParamHandler & rhs) nogil except +
        # libcpp_vector[ String ]  getSubsections() nogil except +

        void setParameters(Param &param)  nogil except +
        Param getParameters()  nogil except +
        Param getDefaults()  nogil except +
        String getName()  nogil except +
        void setName(String)  nogil except +
