from Types cimport *
from String cimport *
from DigestionEnzymeProtein cimport *

cdef extern from "<OpenMS/CHEMISTRY/ProteaseDB.h>" namespace "OpenMS":

    cdef cppclass ProteaseDB "OpenMS::ProteaseDB":
        # wrap-manual-memory:
        #     cdef AutowrapPtrHolder[_ProteaseDB] inst

        # protected
        ProteaseDB() nogil except + #wrap-ignore

        const DigestionEnzymeProtein* getEnzyme(const String& name) nogil except +
        const DigestionEnzymeProtein* getEnzymeByRegEx(const String& cleavage_regex) nogil except +
        void getAllNames(libcpp_vector[ String ]& all_names) nogil except +
        void getAllXTandemNames(libcpp_vector[ String ]& all_names) nogil except + # wrap-doc:Returns all the enzyme names available for XTandem
        void getAllOMSSANames(libcpp_vector[ String ]& all_names) nogil except + # wrap-doc:Returns all the enzyme names available for OMSSA
        void getAllCometNames(libcpp_vector[ String ]& all_names) nogil except + # wrap-doc:Returns all the enzyme names available for Comet
        void getAllMSGFNames(libcpp_vector[ String ] & all_names) nogil except + # wrap-doc:Returns all the enzyme names available for MSGFPlus
        bool hasEnzyme(const String& name) nogil except +
        bool hasRegEx(const String& cleavage_regex) nogil except +
        # bool hasEnzyme(DigestionEnzymeProtein* enzyme) nogil except + # does not make sense as the ptr wont match

        # ConstEnzymeIterator beginEnzyme() nogil except +
        # ConstEnzymeIterator endEnzyme() nogil except +

## wrap static methods
cdef extern from "<OpenMS/CHEMISTRY/ProteaseDB.h>" namespace "OpenMS::ProteaseDB":

    ProteaseDB* getInstance() nogil except + # wrap-ignore
