from libcpp cimport bool
from Types cimport *
from String cimport *

cdef extern from "<OpenMS/SYSTEM/BuildInfo.h>" namespace "OpenMS::Internal":

    cdef cppclass OpenMSOSInfo:

        OpenMSOSInfo() except + nogil 
        OpenMSOSInfo(OpenMSOSInfo &) except + nogil  # compiler
        String getOSAsString() except + nogil 
        String getArchAsString() except + nogil 
        String getOSVersionAsString() except + nogil 

    cdef cppclass OpenMSBuildInfo:

        OpenMSBuildInfo() except + nogil 
        OpenMSBuildInfo(OpenMSBuildInfo &) except + nogil  # compiler


cdef extern from "<OpenMS/SYSTEM/BuildInfo.h>" namespace "OpenMS::Internal::OpenMSOSInfo":

    OpenMSOSInfo getOSInfo() except + nogil  # wrap-attach:OpenMSOSInfo

    String getBinaryArchitecture() except + nogil  # wrap-attach:OpenMSOSInfo


cdef extern from "<OpenMS/SYSTEM/BuildInfo.h>" namespace "OpenMS::Internal::OpenMSBuildInfo":

    bool isOpenMPEnabled() except + nogil  # wrap-attach:OpenMSBuildInfo

    String getBuildType() except + nogil  # wrap-attach:OpenMSBuildInfo

    Size getOpenMPMaxNumThreads() except + nogil  # wrap-attach:OpenMSBuildInfo

    void setOpenMPNumThreads(Int num_threads) except + nogil  # wrap-attach:OpenMSBuildInfo
