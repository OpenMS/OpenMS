from libcpp cimport bool
from String cimport *
from Types cimport *

cdef extern from "<OpenMS/CONCEPT/VersionInfo.h>" namespace "OpenMS":

    cdef cppclass VersionInfo:
        pass

cdef extern from "<OpenMS/CONCEPT/VersionInfo.h>" namespace "OpenMS::VersionInfo":

    cdef cppclass VersionDetails:
        Int version_major
        Int version_minor
        Int version_patch
        # VersionDetails EMPTY

        VersionDetails() nogil except +
        VersionDetails(VersionDetails) nogil except +   #wrap-ignore
        bool operator<(VersionDetails) nogil except +
        bool operator==(VersionDetails) nogil except +
        bool operator>(VersionDetails) nogil except +

    VersionDetails getVersionStruct() nogil except +  #wrap-attach:VersionInfo
    String getVersion()  nogil except +  #wrap-attach:VersionInfo
    String getTime()     nogil except +  #wrap-attach:VersionInfo
    String getRevision() nogil except +  #wrap-attach:VersionInfo

cdef extern from "<OpenMS/CONCEPT/VersionInfo.h>" namespace "OpenMS::VersionInfo::VersionDetails":

    VersionDetails create(String) #wrap-attach:VersionDetails


