from libcpp cimport bool
from String cimport *
from Types cimport *
# from SampleTreatment cimport *

cdef extern from "<OpenMS/METADATA/Digestion.h>" namespace "OpenMS":
    
    cdef cppclass Digestion "OpenMS::Digestion":
        Digestion() nogil except +
        Digestion(Digestion) nogil except +
        # bool operator==(SampleTreatment &rhs) nogil except +
        # SampleTreatment * clone() nogil except +
        String  getEnzyme() nogil except +
        void setEnzyme(String &enzyme) nogil except +
        DoubleReal getDigestionTime() nogil except +
        void setDigestionTime(DoubleReal digestion_time) nogil except +
        DoubleReal getTemperature() nogil except +
        void setTemperature(DoubleReal temperature) nogil except +
        DoubleReal getPh() nogil except +
        void setPh(DoubleReal ph) nogil except +

