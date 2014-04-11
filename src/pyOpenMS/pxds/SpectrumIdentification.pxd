from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool

from MetaInfoInterface cimport *
from IdentificationHit cimport *

cdef extern from "<OpenMS/METADATA/SpectrumIdentification.h>" namespace "OpenMS":

    cdef cppclass SpectrumIdentification(MetaInfoInterface):
        # wrap-inherits:
        #   MetaInfoInterface

        SpectrumIdentification()   nogil except +
        SpectrumIdentification(SpectrumIdentification) nogil except + # wrap-ignore

        # /// sets the identification hits of this spectrum identification (corresponds to single peptide hit in the list)
        void setHits(libcpp_vector[IdentificationHit] & hits) nogil except +

        # /// adds a single identification hit to the hits
        void addHit(IdentificationHit & hit) nogil except +

        # /// returns the identificatio hits of this spectrum identification
        libcpp_vector[IdentificationHit] getHits() nogil except +

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
