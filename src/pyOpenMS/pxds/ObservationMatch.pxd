from String cimport *
from MetaData cimport *
from libcpp.set cimport set as libcpp_set
#from libcpp.optional cimport optional as libcpp_optional
from InputFile cimport *
from Observation cimport *
from IdentifiedMolecule cimport *
from PeptideHit cimport *
from AdductInfo cimport *

#need to add optional to get PeakAnnotationSteps

cdef extern from "<OpenMS/METADATA/ID/ObservationMatch.h>" namespace "OpenMS::IdentificationDataInternal":


  cdef cppclass PeakAnnotations:
    PeakAnnotations() nogil except +
    PeakAnnotations(const PeakAnnotations & other) nogil except +
    bool operator!=(const PeakAnnotations & other) nogil except +
    bool operator<(const PeakAnnotations & other) nogil except +

  cdef cppclass PeakAnnotationSteps:
    PeakAnnotationSteps() nogil except +
    PeakAnnotationSteps(const PeakAnnotationSteps & other) nogil except +
  
    #TODO add deref here?

  #cdef AdductCompare:
  #  bool operator()(AMSE_AdductInfo & left, AMSE_AdductInfo & right) #wrap-cast:toBool

  #cdef Adducts:
  #  Adducts() nogil except + 
  #  Adducts(const Adducts & other) nogil except +

  cdef cppclass AdductRef:
    AdductRef() nogil except +
    AdductRef(const AdductRef & other) nogil except +
    bool operator!=(const AdductRef & other) nogil except +
    bool operator<(const AdductRef & other) nogil except +
    #bool operator=(const AdductRef & other) nogil except +
    AMSE_AdductInfo deref()
    
  #ctypedef licpp_optional[AdductRef] AdductOpt FIXME add optioanl
    
  cdef cppclass ObservationMatch(ScoredProcessingResult):

    IdentifiedMolecule identified_molecule_var
    ObservationRef observation_ref
    int charge
    #AdductOpt adduct_opt
    PeakAnnotationSteps peak_annotations

    ObservationMatch() nogil except +

    ObservationMatch(ObservationMatch & other) nogil except +

    #ObservationMatch(IdentifiedMolecule identified_molecule_var, ObservationRef observation_ref, int charge, const libcpp_optional[AdductRef] & adduct_opt, const AppliedProcessingSteps & steps_and_scores, const PeakAnnotationSteps & peak_annotations) nogil except +

    ObservationMatch & merge(const ObservationMatch & other) nogil except +


  cdef cppclass ObservationMatches:
    ObservationMatches() nogil except +
    ObservationMatches(ObservationMatches & other) nogil except +
    #ObservationMatches operator[](size_t index) #wrap-upper-limit:size() #TODO: Add some sort of access to get the ObservationMatchess back out


    
  cdef cppclass ObservationMatchRef:
      ObservationMatchRef() nogil except +
      ObservationMatchRef(const ObservationMatchRef & other) nogil except +
      bool operator!=(const ObservationMatchRef & other) nogil except +
      bool operator<(const ObservationMatchRef & other) nogil except +
      ObservationMatch deref()
      