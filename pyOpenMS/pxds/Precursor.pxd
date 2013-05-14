from libcpp.string cimport string as libcpp_string
#from InstrumentSettings cimport *
from CVTermList cimport *
from Peak1D cimport *
from Map cimport *

cdef extern from "<OpenMS/METADATA/Precursor.h>" namespace "OpenMS":

    cdef cppclass Precursor(Peak1D):
        # wrap-inherits:
        #    Peak1D
        Precursor()           nogil except +
        Precursor(Precursor)           nogil except +

        double getActivationEnergy() nogil except +
        void setActivationEnergy(double activation_energy) nogil except +

        double getIsolationWindowLowerOffset() nogil except +
        void setIsolationWindowLowerOffset(double bound) nogil except +

        double getIsolationWindowUpperOffset() nogil except +
        void setIsolationWindowUpperOffset(double bound) nogil except +

        int getCharge() nogil except +
        void setCharge(int charge) nogil except +

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

        # as cython has problems with inheriting overriden methods
        # we can not inherit from CVTerm but have to repeat the
        # methods here:
        void setCVTerms(libcpp_vector[CVTerm] & terms)  nogil except +
        void replaceCVTerm(CVTerm & term)               nogil except +

        void replaceCVTerms(libcpp_vector[CVTerm] cv_terms,
                             String accession
                            ) nogil except +

        void replaceCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map
                            ) nogil except +

        Map[String, libcpp_vector[CVTerm] ] getCVTerms()
        void addCVTerm(CVTerm & term)                   nogil except +

        bool operator==(Precursor)  nogil except +
        bool operator!=(Precursor)  nogil except +

        bool hasCVTerm(String accessoin)  nogil except +
        bool empty()                      nogil except +
