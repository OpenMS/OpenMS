from libcpp cimport set as libcpp_set
from MSExperiment cimport *

cdef extern from "<OpenMS/IONMOBILITY/FAIMSHelper.h>" namespace "OpenMS":

    cdef cppclass FAIMSHelper:

        FAIMSHelper() nogil except +
        FAIMSHelper(FAIMSHelper) nogil except + # wrap-ignore


# COMMENT: wrap static methods
cdef extern from "<OpenMS/IONMOBILITY/FAIMSHelper.h>" namespace "OpenMS::FAIMSHelper":
        
        # static members
        libcpp_set[ double ] getCompensationVoltages(const MSExperiment & exp) nogil except +  # wrap-attach:FAIMSHelper
