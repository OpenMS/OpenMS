from Types cimport *
from String cimport *
from DigestionEnzymeProtein cimport *

cdef extern from "<OpenMS/CHEMISTRY/ProteaseDB.h>" namespace "OpenMS":

    cdef cppclass ProteaseDB "OpenMS::ProteaseDB":
        # wrap-manual-memory:
        #    cdef AutowrapPtrHolder[_ProteaseDB] inst
        # wrap-doc:
        #  Database of enzymes that digest proteins (proteases). This is a singleton class.
        #  The enzymes are read from share/CHEMISTRY/Enzymes.xml.

        # protected
        ProteaseDB() except + nogil  #wrap-ignore

        const DigestionEnzymeProtein* getEnzyme(const String& name) except + nogil  # wrap-doc:Returns the enzyme with the given name (supports synonym names)
        const DigestionEnzymeProtein* getEnzymeByRegEx(const String& cleavage_regex) except + nogil  # wrap-doc:Returns the enzyme with the given cleavage regex
        void getAllNames(libcpp_vector[ String ]& all_names) except + nogil  # wrap-doc:Returns all enzyme names (does not include synonym names)
        void getAllXTandemNames(libcpp_vector[ String ]& all_names) except + nogil  # wrap-doc:Returns all the enzyme names available for XTandem
        void getAllOMSSANames(libcpp_vector[ String ]& all_names) except + nogil  # wrap-doc:Returns all the enzyme names available for OMSSA
        void getAllCometNames(libcpp_vector[ String ]& all_names) except + nogil  # wrap-doc:Returns all the enzyme names available for Comet
        void getAllMSGFNames(libcpp_vector[ String ] & all_names) except + nogil  # wrap-doc:Returns all the enzyme names available for MSGFPlus
        bool hasEnzyme(const String& name) except + nogil  # wrap-doc:Returns True if the database contains an enzyme with the given name
        bool hasRegEx(const String& cleavage_regex) except + nogil  # wrap-doc:Returns True if the database contains an enzyme with the given cleavage regex
        # bool hasEnzyme(DigestionEnzymeProtein* enzyme) except + nogil  # does not make sense as the ptr wont match

        # ConstEnzymeIterator beginEnzyme() except + nogil
        # ConstEnzymeIterator endEnzyme() except + nogil

## wrap static methods
cdef extern from "<OpenMS/CHEMISTRY/ProteaseDB.h>" namespace "OpenMS::ProteaseDB":

    ProteaseDB* getInstance() except + nogil  # wrap-ignore
