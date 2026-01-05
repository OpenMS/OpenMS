from libcpp cimport bool
from String cimport *
from Types cimport *

cdef extern from "<OpenMS/CONCEPT/VersionInfo.h>" namespace "OpenMS":

    cdef cppclass VersionInfo:

        @staticmethod
        VersionDetails getVersionStruct() except + nogil

        @staticmethod
        String getVersion() except + nogil

        @staticmethod
        String getTime() except + nogil

        @staticmethod
        String getRevision() except + nogil

        @staticmethod
        String getBranch() except + nogil

cdef extern from "<OpenMS/CONCEPT/VersionInfo.h>" namespace "OpenMS::VersionInfo":

    cdef cppclass VersionDetails:
        Int version_major
        Int version_minor
        Int version_patch
        String pre_release_identifier

        VersionDetails() except + nogil
        VersionDetails(VersionDetails &) except + nogil
        bool operator<(VersionDetails) except + nogil
        bool operator==(VersionDetails) except + nogil
        bool operator>(VersionDetails) except + nogil

        @staticmethod
        VersionDetails create(String) except + nogil
