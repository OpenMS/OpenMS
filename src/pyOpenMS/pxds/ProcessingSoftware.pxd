from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from Software cimport *
from String cimport *
from ScoreType cimport *
from MetaData cimport *


cdef extern from "<OpenMS/METADATA/ID/ProcessingSoftware.h>" namespace "OpenMS::IdentificationDataInternal":
  cdef cppclass ProcessingSoftware(Software):

    ProcessingSoftware() nogil except +

    ProcessingSoftware(ProcessingSoftware other) nogil except +

    libcpp_vector[ ScoreTypeRef] assigned_scores 

    ProcessingSoftware(String name, String version, libcpp_vector[ ScoreTypeRef ] assigned_scores) nogil except +

  ctypedef libcpp_set[ ProcessingSoftware ] ProcessingSoftwares
  #ctypedef libcpp_set[ ProcessingSoftware ].iterator 

  cdef cppclass ProcessingSoftwareRef:
    ProcessingSoftwareRef() nogil except +
    ProcessingSoftwareRef(const ProcessingSoftwareRef & other) nogil except +
    bool operator!=(const ProcessingSoftwareRef & other) nogil except +
    bool operator<(const ProcessingSoftwareRef & other) nogil except +
    ProcessingSoftware deref() nogil except +


