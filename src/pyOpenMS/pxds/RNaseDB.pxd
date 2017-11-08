from Types cimport *
from String cimport *
from DigestionEnzymeRNA cimport *

cdef extern from "<OpenMS/CHEMISTRY/RNaseDB.h>" namespace "OpenMS":

    cdef cppclass RNaseDB "OpenMS::RNaseDB":
        # wrap-manual-memory:
        #     cdef AutowrapPtrHolder[_RNaseDB] inst

        RNaseDB() nogil except + #wrap-ignore
        RNaseDB(RNaseDB) nogil except + #wrap-ignore

        const DigestionEnzymeRNA* getEnzyme(String& name) nogil except +
        const DigestionEnzymeRNA* getEnzymeByRegEx(String& cleavage_regex) nogil except +
        void getAllNames(libcpp_vector[ String ]& all_names) nogil except +
        bool hasEnzyme(String& name) nogil except +
        bool hasRegEx(String& cleavage_regex) nogil except +
        # bool hasEnzyme(DigestionEnzymeRNA* enzyme) nogil except + # does not make sense as the ptr wont match

        # ConstEnzymeIterator beginEnzyme() nogil except +
        # ConstEnzymeIterator endEnzyme() nogil except +

## wrap static methods
cdef extern from "<OpenMS/CHEMISTRY/RNaseDB.h>" namespace "OpenMS::RNaseDB":

    RNaseDB* getInstance() nogil except + # wrap-ignore
