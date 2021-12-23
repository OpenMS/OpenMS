from libcpp cimport bool
from Types cimport *
from String cimport *

cdef extern from "<OpenMS/SYSTEM/BuildInfo.h>" namespace "OpenMS::Internal":

    cdef cppclass OpenMSOSInfo:

        OpenMSOSInfo() nogil except +
        OpenMSOSInfo(OpenMSOSInfo &) nogil except + # compiler
        String getOSAsString() nogil except +
        String getArchAsString() nogil except +
        String getOSVersionAsString() nogil except +

    cdef cppclass OpenMSBuildInfo:

        OpenMSBuildInfo() nogil except +
        OpenMSBuildInfo(OpenMSBuildInfo &) nogil except + # compiler


cdef extern from "<OpenMS/SYSTEM/BuildInfo.h>" namespace "OpenMS::Internal::OpenMSOSInfo":

    OpenMSOSInfo getOSInfo() nogil except + # wrap-attach:OpenMSOSInfo

    String getBinaryArchitecture() nogil except + # wrap-attach:OpenMSOSInfo


cdef extern from "<OpenMS/SYSTEM/BuildInfo.h>" namespace "OpenMS::Internal::OpenMSBuildInfo":

    bool isOpenMPEnabled() nogil except + # wrap-attach:OpenMSBuildInfo

    String getBuildType() nogil except + # wrap-attach:OpenMSBuildInfo

    Size getOpenMPMaxNumThreads() nogil except + # wrap-attach:OpenMSBuildInfo
