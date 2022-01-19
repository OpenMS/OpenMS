from Types cimport *
from libcpp cimport bool
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from MetaInfoInterface cimport *
from ScoreType cimport *

cdef extern from "<OpenMS/METADATA/ID/IdentificationData.h>" namespace "OpenMS":
    
    cdef cppclass IdentificationData(MetaInfoInterface):
        # wrap-inherits:
        #    MetaInfoInterface

        IdentificationData() nogil except + # wrap-doc:Representation of spectrum identification results and associated data

        IdentificationData(IdentificationData &) nogil except + # wrap-doc:Copy constructor

cdef extern from "<OpenMS/METADATA/ID/MetaData.h>" namespace "OpenMS::IdentificationDataInternal":
  cdef cppclass IteratorWrapper[Iterator]:
    # wrap-instances:
    #  ScoreTypeRef := IteratorWrapper[setSTit]

    # wrap-doc:
    #   Class for ScoreTypeRef
    IteratorWrapper() nogil except +
    IteratorWrapper(Iterator it) nogil except +