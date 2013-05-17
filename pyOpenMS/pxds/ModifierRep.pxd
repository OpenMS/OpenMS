from Types cimport *
from libcpp.map cimport map as libcpp_map
from libcpp.vector cimport vector as libcpp_vector
from String cimport *

cdef extern from "<OpenMS/CHEMISTRY/ModifierRep.h>" namespace "OpenMS":
    
    cdef cppclass ModifierRep "OpenMS::ModifierRep":
        ModifierRep() nogil except +
        ModifierRep(ModifierRep) nogil except +
        void setNumberOfModifications(Size i) nogil except +
        Size getNumberOfModifications() nogil except +
        libcpp_vector[ libcpp_vector[ double ] ]  getModificationTable() nogil except +
        void refreshModificationList(libcpp_map[ double, ptrdiff_t ] & mod_map, char & c) nogil except +
        Size getMaxModificationMasses() nogil except +
        # TODO immutable types by reference
        # libcpp_vector[ String ] getModificationsForMass(double & m) nogil except +
        # libcpp_vector[ String ] getModificationsForMass(double & m, String & seq) nogil except +

