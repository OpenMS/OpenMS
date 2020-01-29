from Types cimport *
from libcpp cimport bool
from libcpp.set cimport set as libcpp_set
from libcpp.vector cimport vector as libcpp_vector
from DataProcessing cimport *
from DataProcessingSoftware cimport *
from MetaInfoInterface cimport *

cdef extern from "<Openms/METADATA/ID/DataProcessingStep.h>" namespace "OpenMS::IdentificationDataInternal":

    cdef cppclass DataProcessingStep(MetaInfoInterface) :
        # wrap-inherits:
        #   MetaInfoInterface
        DataProcessingStep() nogil except +
        DataProcessingStep(DataProcessingStep) nogil except +

        ProcessingSoftwareRef software_ref
        libcpp_vector[ InputFileRef ] input_file_refs
        libcpp_vector[ String ] primary_files
        DateTime date_time
        libcpp_set[ ProcessingAction ] actions
        DataProcessingStep(ProcessingSoftwareRef software_ref, libcpp_vector[ InputFileRef ] & input_file_refs, libcpp_vector[ String ] & primary_files, DateTime & date_time, libcpp_set[ ProcessingAction ] actions) nogil except +
        bool operator <(DataProcessingStep & other) nogil except +
        bool operator ==(DataProcessingStep & other) nogil except +