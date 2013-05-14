from IonSource cimport *
from MetaInfoInterface cimport *

cdef extern from "<OpenMS/METADATA/InstrumentSettings.h>" namespace "OpenMS":

    cdef cppclass InstrumentSettings(MetaInfoInterface):
        # wrap-inherits:
        #    MetaInfoInterface

        InstrumentSettings()     nogil except +
        InstrumentSettings(InstrumentSettings)     nogil except +

        Polarity getPolarity()     nogil except +
        void setPolarity(Polarity)  nogil except +

        # declare again: cython complains for overloaded methods in base
        # classes
        void getKeys(libcpp_vector[String] & keys) nogil except +
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except +
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +
