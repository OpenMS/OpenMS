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
        double getDigestionTime() nogil except +
        void setDigestionTime(double digestion_time) nogil except +
        double getTemperature() nogil except +
        void setTemperature(double temperature) nogil except +
        double getPh() nogil except +
        void setPh(double ph) nogil except +

