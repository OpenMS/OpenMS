from String cimport *
from MetaInfoInterface cimport *
from MetaData cimport *
from libcpp.set cimport set as libcpp_set
from InputFile cimport *


cdef extern from "<OpenMS/METADATA/ID/Observation.h>" namespace "OpenMS::IdentificationDataInternal":
  cdef cppclass Observation(MetaInfoInterface):

    Observation() nogil except +

    Observation(Observation & other) nogil except +

    String data_id

    InputFileRef input_file

    double rt

    double mz

    Observation(const String & data_id, const InputFileRef & input_file, double rt, double mz) # Note that default args don't work

  cdef cppclass Observations:
    Observations() nogil except +
    Observations(Observations & other) nogil except +
    #Observation operator[](size_t index) #wrap-upper-limit:size() #TODO: Add some sort of access to get the Observations back out


    
  cdef cppclass ObservationRef:
      ObservationRef() nogil except +
      ObservationRef(const ObservationRef & other) nogil except +
      bool operator!=(const ObservationRef & other) nogil except +
      bool operator<(const ObservationRef & other) nogil except +
      #bool operator=(const ObservationRef & other) nogil except +
      Observation deref()
      