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
        String pre_release_identifier

        VersionDetails() except + nogil 
        VersionDetails(VersionDetails &) except + nogil 
        bool operator<(VersionDetails) except + nogil 
        bool operator==(VersionDetails) except + nogil 
        bool operator>(VersionDetails) except + nogil 

    VersionDetails getVersionStruct() except + nogil   #wrap-attach:VersionInfo
    String getVersion()  except + nogil   #wrap-attach:VersionInfo
    String getTime()     except + nogil   #wrap-attach:VersionInfo
    String getRevision() except + nogil   #wrap-attach:VersionInfo
    String getBranch()   except + nogil   #wrap-attach:VersionInfo

cdef extern from "<OpenMS/CONCEPT/VersionInfo.h>" namespace "OpenMS::VersionInfo::VersionDetails":

    VersionDetails create(String) #wrap-attach:VersionDetails

