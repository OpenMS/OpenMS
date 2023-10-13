from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool

from MetaInfoInterface cimport *
from IdentificationHit cimport *

cdef extern from "<OpenMS/METADATA/SpectrumIdentification.h>" namespace "OpenMS":

    cdef cppclass SpectrumIdentification(MetaInfoInterface):
        # wrap-inherits:
        #  MetaInfoInterface

        SpectrumIdentification() except + nogil 
        SpectrumIdentification(SpectrumIdentification &) except + nogil 

        void setHits(libcpp_vector[IdentificationHit] & hits) except + nogil  # wrap-doc:Sets the identification hits of this spectrum identification (corresponds to single peptide hit in the list)

        void addHit(IdentificationHit & hit) except + nogil  # wrap-doc:Adds a single identification hit to the hits

        libcpp_vector[IdentificationHit] getHits() except + nogil  # wrap-doc:Returns the identification hits of this spectrum identification
