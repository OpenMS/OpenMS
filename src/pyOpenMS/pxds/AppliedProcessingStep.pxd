from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from libcpp.pair cimport pair as libcpp_pair
from libcpp.map cimport map as libcpp_map
from DataProcessing cimport *
from MetaInfoInterface cimport *
from MetaData cimport *
from ProcessingSoftware cimport *
from ProcessingStep cimport *
from ScoreType cimport *
from InputFile cimport *
from DateTime cimport *

cdef extern from "<OpenMS/METADATA/ID/AppliedProcessingStep.h>" namespace "OpenMS::IdentificationDataInternal":
  cdef cppclass AppliedProcessingStep:

    libcpp_map[IteratorWrapper[setSTit, ScoreType], double] scores

    # boost_optional[IteratorWrapper[setPSoftSit, ProcessingStep]] processing_step_opt #FIXME

    AppliedProcessingStep() nogil except +
    
    AppliedProcessingStep(AppliedProcessingStep & other)

    # AppliedProcessingStep(boost_optional[IteratorWrapper[setPSoftSit, ProcessingStep]]& processing_step_opt, libcpp_map[IteratorWrapper[setSTit, ScoreType], double] scores)

    bool operator==(AppliedProcessingStep & other) nogil except +

  #  libcpp_vector[libcpp_pair[IteratorWrapper[setSTit, ScoreType],double]] getScoresInOrder(bool primary_only) nogil except + #FIXME

  cdef cppclass AppliedProcessingSteps:
    AppliedProcessingSteps() nogil except +
    AppliedProcessingSteps(AppliedProcessingSteps & other) nogil except +

