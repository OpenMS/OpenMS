from Types cimport *
from String cimport *
from DigestionEnzymeRNA cimport *

cdef extern from "<OpenMS/CHEMISTRY/RNaseDB.h>" namespace "OpenMS":

    cdef cppclass RNaseDB "OpenMS::RNaseDB":
        # wrap-manual-memory:
        #    cdef AutowrapPtrHolder[_RNaseDB] inst

        # protected
        RNaseDB() except + nogil  # wrap-ignore
        # due to wrap-manual-memory
        RNaseDB(RNaseDB &) except + nogil  # wrap-ignore

        const DigestionEnzymeRNA* getEnzyme(const String& name) except + nogil 
        const DigestionEnzymeRNA* getEnzymeByRegEx(const String& cleavage_regex) except + nogil 
        void getAllNames(libcpp_vector[ String ]& all_names) except + nogil 
        bool hasEnzyme(const String& name) except + nogil 
        # bool hasRegEx(const String& cleavage_regex) except + nogil # We don't use regexes for RNA
        # bool hasEnzyme(DigestionEnzymeRNA* enzyme) except + nogil  # does not make sense as the ptr wont match

        # ConstEnzymeIterator beginEnzyme() except + nogil 
        # ConstEnzymeIterator endEnzyme() except + nogil 

## wrap static methods
cdef extern from "<OpenMS/CHEMISTRY/RNaseDB.h>" namespace "OpenMS::RNaseDB":

    RNaseDB* getInstance() except + nogil  # wrap-ignore
