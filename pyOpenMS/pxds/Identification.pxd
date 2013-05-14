from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool

from MetaInfoInterface cimport *

from DateTime cimport *
from SpectrumIdentification cimport *

cdef extern from "<OpenMS/METADATA/Identification.h>" namespace "OpenMS":

    cdef cppclass Identification(MetaInfoInterface):
        # wrap-inherits:
        #   MetaInfoInterface

        Identification()   nogil except +
        Identification(Identification) nogil except + # wrap-ignore


        void setCreationDate(DateTime date)
        DateTime getCreationDate()

        # /// sets the spectrum identifications
        void setSpectrumIdentifications(libcpp_vector[SpectrumIdentification] & ids)

        # /// adds a spectrum identification
        void addSpectrumIdentification(SpectrumIdentification & id)

        # /// returns the spectrum identifications stored
        libcpp_vector[SpectrumIdentification] getSpectrumIdentifications() 

        # COPY-PASTE from MetaInfoInterface
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
