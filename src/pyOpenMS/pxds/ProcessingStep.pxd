from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from DataProcessing cimport *
from MetaInfoInterface cimport *
from MetaData cimport *
from ProcessingSoftware cimport *
from InputFile cimport *
from DateTime cimport *

cdef extern from "<OpenMS/METADATA/ID/ProcessingStep.h>" namespace "OpenMS::IdentificationDataInternal":
  cdef cppclass ProcessingStep(MetaInfoInterface):

    ProcessingStep() nogil except +
    
    IteratorWrapper[setPSit,ProcessingSoftware] software_ref
    
    libcpp_vector[InputFileRef] input_file_refs

    DateTime date_time
    libcpp_set[ProcessingAction] actions

    ProcessingStep(IteratorWrapper[setPSit,ProcessingSoftware] software_ref, \
      libcpp_vector[InputFileRef] input_file_refs, \
      DateTime & date_time, \
      libcpp_set[ProcessingAction] actions) #Note that we don't have default args

    ProcessingStep(ProcessingStep & other)

    bool operator<(ProcessingStep & other) nogil except +
    bool operator==(ProcessingStep & other) nogil except +
  ctypedef libcpp_set[ProcessingStep] ProcessingSteps
  ctypedef libcpp_set[ProcessingStep].iterator setPSoftSit