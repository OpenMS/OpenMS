from libcpp.vector cimport vector as libcpp_vector
from Param cimport *
from String cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/DefaultParamHandler.h>" namespace "OpenMS":

    cdef cppclass DefaultParamHandler:
        #wrap-ignore
        #no-pxd-import

        DefaultParamHandler(String name) except + nogil 
        DefaultParamHandler(DefaultParamHandler &) except + nogil 
        libcpp_vector[ String ] getSubsections() except + nogil 

        void setParameters(Param &param)  except + nogil  # wrap-doc:Sets the parameters
        Param getParameters()  except + nogil  # wrap-doc:Returns the parameters
        Param getDefaults()  except + nogil  # wrap-doc:Returns the default parameters
        String getName()  except + nogil  # wrap-doc:Returns the name
        void setName(const String&)  except + nogil  # wrap-doc:Sets the name
