from Types cimport *
from libcpp.map cimport map as libcpp_map
from libcpp.vector cimport vector as libcpp_vector
from ProteinIdentification cimport *
from PeptideIdentification cimport *
from XMLHandler cimport *
from XMLFile cimport *
from AASequence cimport *

cdef extern from "<OpenMS/FORMAT/PepXMLFileMascot.h>" namespace "OpenMS":
    
    cdef cppclass PepXMLFileMascot :

        PepXMLFileMascot() nogil except +
        PepXMLFileMascot(PepXMLFileMascot) nogil except + #wrap-ignore

        # TODO map
        # void load(const String & filename, libcpp_map[ String, libcpp_vector[ AASequence ] ] & peptides) nogil except +

