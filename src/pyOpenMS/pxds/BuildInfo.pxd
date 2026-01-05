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

        @staticmethod
        OpenMSOSInfo getOSInfo() except + nogil

        @staticmethod
        String getBinaryArchitecture() except + nogil

    cdef cppclass OpenMSBuildInfo:

        OpenMSBuildInfo() except + nogil
        OpenMSBuildInfo(OpenMSBuildInfo &) except + nogil  # compiler

        @staticmethod
        bool isOpenMPEnabled() except + nogil

        @staticmethod
        String getBuildType() except + nogil

        @staticmethod
        Size getOpenMPMaxNumThreads() except + nogil

        @staticmethod
        void setOpenMPNumThreads(Int num_threads) except + nogil
