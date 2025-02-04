from Types cimport *
from String cimport *
from DigestionEnzymeProtein cimport *

cdef extern from "<OpenMS/CHEMISTRY/ProteaseDB.h>" namespace "OpenMS":

    cdef cppclass ProteaseDB "OpenMS::ProteaseDB":
        # wrap-manual-memory:
        #    cdef AutowrapPtrHolder[_ProteaseDB] inst

        # protected
        ProteaseDB() except + nogil  #wrap-ignore

        const DigestionEnzymeProtein* getEnzyme(const String& name) except + nogil 
        const DigestionEnzymeProtein* getEnzymeByRegEx(const String& cleavage_regex) except + nogil 
        void getAllNames(libcpp_vector[ String ]& all_names) except + nogil 
        void getAllXTandemNames(libcpp_vector[ String ]& all_names) except + nogil  # wrap-doc:Returns all the enzyme names available for XTandem
        void getAllOMSSANames(libcpp_vector[ String ]& all_names) except + nogil  # wrap-doc:Returns all the enzyme names available for OMSSA
        void getAllCometNames(libcpp_vector[ String ]& all_names) except + nogil  # wrap-doc:Returns all the enzyme names available for Comet
        void getAllMSGFNames(libcpp_vector[ String ] & all_names) except + nogil  # wrap-doc:Returns all the enzyme names available for MSGFPlus
        bool hasEnzyme(const String& name) except + nogil 
        bool hasRegEx(const String& cleavage_regex) except + nogil 
        # bool hasEnzyme(DigestionEnzymeProtein* enzyme) except + nogil  # does not make sense as the ptr wont match

        # ConstEnzymeIterator beginEnzyme() except + nogil 
        # ConstEnzymeIterator endEnzyme() except + nogil 

## wrap static methods
cdef extern from "<OpenMS/CHEMISTRY/ProteaseDB.h>" namespace "OpenMS::ProteaseDB":

    ProteaseDB* getInstance() except + nogil  # wrap-ignore
