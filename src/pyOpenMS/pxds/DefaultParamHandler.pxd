from libcpp.vector cimport vector as libcpp_vector
from Param cimport *
from String cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/DefaultParamHandler.h>" namespace "OpenMS":

    cdef cppclass DefaultParamHandler:
        #wrap-ignore
        #no-pxd-import

        DefaultParamHandler(String name) nogil except +
        DefaultParamHandler(DefaultParamHandler &) nogil except +
        libcpp_vector[ String ] getSubsections() nogil except +

        void setParameters(Param &param)  nogil except + # wrap-doc:Sets the parameters
        Param getParameters()  nogil except + # wrap-doc:Returns the parameters
        Param getDefaults()  nogil except + # wrap-doc:Returns the default parameters
        String getName()  nogil except + # wrap-doc:Returns the name
        void setName(const String&)  nogil except + # wrap-doc:Sets the name
