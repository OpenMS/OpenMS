from libcpp.string cimport string as libcpp_string
#from InstrumentSettings cimport *

cdef extern from "<OpenMS/METADATA/Precursor.h>" namespace "OpenMS":

    cdef cppclass Precursor:
        Precursor()           nogil except +
        Precursor(Precursor)           nogil except +
        double getMZ()           nogil except +
        double getIntensity()           nogil except +
        void setMZ(double ) nogil except +
        void setIntensity(double ) nogil except +
