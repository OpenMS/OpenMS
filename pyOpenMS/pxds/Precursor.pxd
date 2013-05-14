from libcpp.string cimport string as libcpp_string
#from InstrumentSettings cimport *
from CVTermList cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/METADATA/Precursor.h>" namespace "OpenMS":

    cdef cppclass Precursor(CVTermList, Peak1D):
        # wrap-inherits:
        #    CVTermList
        #    Peak1D
        Precursor()           nogil except +
        Precursor(Precursor)           nogil except +

        double getActivationEnergy()
        void setActivationEnergy(double activation_energy)

        double getIsolationWindowLowerOffset()
        void setIsolationWindowLowerOffset(double bound)

        double getIsolationWindowUpperOffset()
        void setIsolationWindowUpperOffset(double bound)

        int getCharge()
        void setCharge(int charge)

        # Inherited from MetaInfoInterface - copyNpaste
        void getKeys(libcpp_vector[String] & keys)
        void getKeys(libcpp_vector[unsigned int] & keys)
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +
