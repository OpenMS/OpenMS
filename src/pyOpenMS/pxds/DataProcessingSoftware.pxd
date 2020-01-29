from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from Software cimport *
from ScoreType cimport *

cdef extern from "<OpenMS/METADATA/ID/DataProcessingSoftware.h>" namespace "OpenMS::IdentificationDataInternal":

    cdef cppclass DataProcessingSoftware(Software) :
        # wrap-inherits:
        #   Software
        DataProcessingSoftware() nogil except +
        DataProcessingSoftware(DataProcessingSoftware) nogil except + # auto gen code has this as ignore
        libcpp_vector[ ScoreTypeRef ] assigned_scores

        DataProcessingSoftware(const String & name, const String & version, libcpp_vector[ ScoreTypeRef ] assigned_scores) nogil except +

    ctypedef libcpp_set[ DataProcessingSoftware ] DataProcessingSoftwares
    ctypedef IteratorWrapper[ DataProcessingSoftwares].Iterator ProcessingSoftwareRef