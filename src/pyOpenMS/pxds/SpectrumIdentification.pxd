from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool

from MetaInfoInterface cimport *
from IdentificationHit cimport *

cdef extern from "<OpenMS/METADATA/SpectrumIdentification.h>" namespace "OpenMS":

    cdef cppclass SpectrumIdentification(MetaInfoInterface):
        # wrap-inherits:
        #   MetaInfoInterface

        SpectrumIdentification() nogil except +
        SpectrumIdentification(SpectrumIdentification &) nogil except +

        void setHits(libcpp_vector[IdentificationHit] & hits) nogil except + # wrap-doc:Sets the identification hits of this spectrum identification (corresponds to single peptide hit in the list)

        void addHit(IdentificationHit & hit) nogil except + # wrap-doc:Adds a single identification hit to the hits

        libcpp_vector[IdentificationHit] getHits() nogil except + # wrap-doc:Returns the identification hits of this spectrum identification
