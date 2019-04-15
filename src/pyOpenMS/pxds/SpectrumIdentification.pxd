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

        # /// returns the identification hits of this spectrum identification
        libcpp_vector[IdentificationHit] getHits() nogil except +

