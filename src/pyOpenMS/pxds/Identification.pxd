from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool

from MetaInfoInterface cimport *

from DateTime cimport *
from SpectrumIdentification cimport *

cdef extern from "<OpenMS/METADATA/Identification.h>" namespace "OpenMS":

    cdef cppclass Identification(MetaInfoInterface):
        # wrap-inherits:
        #   MetaInfoInterface

        Identification() nogil except + # wrap-doc:Represents a object which can store the information of an analysisXML instance
        Identification(Identification &) nogil except +

        void setCreationDate(DateTime date) nogil except + # wrap-doc:Sets the date and time the file was written
        DateTime getCreationDate() nogil except + # wrap-doc:Returns the date and time the file was created

        void setSpectrumIdentifications(libcpp_vector[SpectrumIdentification] & ids) nogil except + # wrap-doc:Sets the spectrum identifications

        void addSpectrumIdentification(SpectrumIdentification & id) nogil except + # wrap-doc:Adds a spectrum identification

        libcpp_vector[SpectrumIdentification] getSpectrumIdentifications()  nogil except + # wrap-doc:Returns the spectrum identifications stored
