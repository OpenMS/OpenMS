from Param cimport *
from String cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/DefaultParamHandler.h>" namespace "OpenMS":

    cdef cppclass DefaultParamHandler:
        #wrap-ignore

        void setParameters(Param &param)  nogil except +
        Param getParameters()  nogil except +
        Param getDefaults()  nogil except +
        String getName()  nogil except +
        void setName(String)  nogil except +
