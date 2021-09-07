from libcpp cimport bool
from String cimport *
from Types cimport *
# from SampleTreatment cimport *

cdef extern from "<OpenMS/METADATA/Digestion.h>" namespace "OpenMS":
    
    cdef cppclass Digestion "OpenMS::Digestion":
        Digestion() nogil except +
        Digestion(Digestion &) nogil except +
        # bool operator==(SampleTreatment &rhs) nogil except +
        # SampleTreatment * clone() nogil except +
        String  getEnzyme() nogil except + # wrap-doc:Returns the enzyme name (default is "")
        void setEnzyme(const String& enzyme) nogil except + # wrap-doc:Sets the enzyme name
        double getDigestionTime() nogil except + # wrap-doc:Returns the digestion time in minutes (default is 0.0)
        void setDigestionTime(double digestion_time) nogil except + # wrap-doc:Sets the digestion time in minutes
        double getTemperature() nogil except + # wrap-doc:Returns the temperature during digestion in degree C (default is 0.0)
        void setTemperature(double temperature) nogil except + # wrap-doc:Sets the temperature during digestion in degree C
        double getPh() nogil except + # wrap-doc:Returns the pH value (default is 0.0)
        void setPh(double ph) nogil except + # wrap-doc:Sets the pH value

