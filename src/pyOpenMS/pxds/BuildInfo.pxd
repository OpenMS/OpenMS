from libcpp cimport bool
from Types cimport *
from String cimport *

cdef extern from "<OpenMS/SYSTEM/BuildInfo.h>" namespace "OpenMS::Internal":

    cdef enum class OpenMS_OS "OpenMS::Internal::OpenMS_OS":
        # wrap-attach:
        #    OpenMSOSInfo
        OS_UNKNOWN, OS_MACOS, OS_WINDOWS, OS_LINUX, SIZE_OF_OPENMS_OS

    cdef enum class OpenMS_Architecture "OpenMS::Internal::OpenMS_Architecture":
        # wrap-attach:
        #    OpenMSOSInfo
        ARCH_UNKNOWN, ARCH_32BIT, ARCH_64BIT, SIZE_OF_OPENMS_ARCHITECTURE

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
