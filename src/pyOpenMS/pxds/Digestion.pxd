from libcpp cimport bool
from String cimport *
from Types cimport *
# from SampleTreatment cimport *

cdef extern from "<OpenMS/METADATA/Digestion.h>" namespace "OpenMS":
    
    cdef cppclass Digestion "OpenMS::Digestion":
        Digestion() except + nogil 
        Digestion(Digestion &) except + nogil 
        # bool operator==(SampleTreatment &rhs) except + nogil 
        # SampleTreatment * clone() except + nogil 
        String  getEnzyme() except + nogil  # wrap-doc:Returns the enzyme name (default is "")
        void setEnzyme(const String& enzyme) except + nogil  # wrap-doc:Sets the enzyme name
        double getDigestionTime() except + nogil  # wrap-doc:Returns the digestion time in minutes (default is 0.0)
        void setDigestionTime(double digestion_time) except + nogil  # wrap-doc:Sets the digestion time in minutes
        double getTemperature() except + nogil  # wrap-doc:Returns the temperature during digestion in degree C (default is 0.0)
        void setTemperature(double temperature) except + nogil  # wrap-doc:Sets the temperature during digestion in degree C
        double getPh() except + nogil  # wrap-doc:Returns the pH value (default is 0.0)
        void setPh(double ph) except + nogil  # wrap-doc:Sets the pH value

