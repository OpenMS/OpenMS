from String cimport *
from MetaData cimport *
from libcpp.set cimport set as libcpp_set
#from libcpp.optional cimport optional as libcpp_optional
from ObservationMatch cimport *


cdef extern from "<OpenMS/METADATA/ID/ObservationMatchGroup.h>" namespace "OpenMS::IdentificationDataInternal":
    
  cdef cppclass ObservationMatchGroup(ScoredProcessingResult):
    ObservationMatchGroup() nogil except +
    ObservationMatchGroup(const ObservationMatchGroup & other) nogil except +
    libcpp_set[ObservationMatchRef] observation_match_refs
    bool allSameMolecule() nogil except +
    bool operator==(const ObservationMatchGroup) nogil except +
    bool operator!=(const ObservationMatchGroup) nogil except +


  cdef cppclass ObservationMatchGroups:
    ObservationMatchGroups() nogil except +
    ObservationMatchGroups(ObservationMatchGroups & other) nogil except +
    #ObservationMatchGroupes operator[](size_t index) #wrap-upper-limit:size() #TODO: Add some sort of access to get the ObservationMatchGroups back out


    
  cdef cppclass MatchGroupRef:
      MatchGroupRef() nogil except +
      MatchGroupRef(const MatchGroupRef & other) nogil except +
      bool operator!=(const MatchGroupRef & other) nogil except +
      bool operator<(const MatchGroupRef & other) nogil except +
      ObservationMatchGroup deref() nogil except +
      