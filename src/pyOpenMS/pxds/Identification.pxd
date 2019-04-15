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


        void setCreationDate(DateTime date) nogil except +
        DateTime getCreationDate() nogil except +

        # /// sets the spectrum identifications
        void setSpectrumIdentifications(libcpp_vector[SpectrumIdentification] & ids) nogil except +

        # /// adds a spectrum identification
        void addSpectrumIdentification(SpectrumIdentification & id) nogil except +

        # /// returns the spectrum identifications stored
        libcpp_vector[SpectrumIdentification] getSpectrumIdentifications()  nogil except +

