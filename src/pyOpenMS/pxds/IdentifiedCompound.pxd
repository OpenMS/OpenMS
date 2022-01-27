from String cimport *
from MetaData cimport *
from libcpp.set cimport set as libcpp_set
from EmpiricalFormula cimport *
from AppliedProcessingStep cimport *
from ScoredProcessingResult cimport *


cdef extern from "<OpenMS/METADATA/ID/IdentifiedCompound.h>" namespace "OpenMS::IdentificationDataInternal":
  cdef cppclass IdentifiedCompound(ScoredProcessingResult):

    IdentifiedCompound() nogil except +

    IdentifiedCompound(IdentifiedCompound & other) nogil except +

    String identifier

    EmpiricalFormula formula

    String name

    String smile

    String inchi

    IdentifiedCompound(const String & identifier, const EmpiricalFormula & formula, const String & name, const String & smile, const String & inchi, const AppliedProcessingSteps & steps_and_scores) # Note that default args don't work

  cdef cppclass IdentifiedCompounds:
    IdentifiedCompounds() nogil except +
    IdentifiedCompounds(IdentifiedCompounds & other) nogil except +
    #IdentifiedCompound operator[](size_t index) #wrap-upper-limit:size() #TODO: Add some sort of access to get the IdentifiedCompounds back out


    
  cdef cppclass IdentifiedCompoundRef:
      IdentifiedCompoundRef() nogil except +
      IdentifiedCompoundRef(const IdentifiedCompoundRef & other) nogil except +
      bool operator!=(const IdentifiedCompoundRef & other) nogil except +
      bool operator<(const IdentifiedCompoundRef & other) nogil except +
      IdentifiedCompound deref()