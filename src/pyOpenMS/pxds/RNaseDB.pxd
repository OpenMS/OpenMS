from Types cimport *
from String cimport *
from DigestionEnzymeRNA cimport *

cdef extern from "<OpenMS/CHEMISTRY/RNaseDB.h>" namespace "OpenMS":

    cdef cppclass RNaseDB "OpenMS::RNaseDB":
        # wrap-manual-memory:
        #    cdef AutowrapPtrHolder[_RNaseDB] inst
        # wrap-doc:
        #  Database of enzymes that digest RNA (RNases). This is a singleton class.
        #  The enzymes are read from share/CHEMISTRY/Enzymes_RNA.xml.

        # protected
        RNaseDB() except + nogil  # wrap-ignore
        # due to wrap-manual-memory
        RNaseDB(RNaseDB &) except + nogil  # wrap-ignore

        const DigestionEnzymeRNA* getEnzyme(const String& name) except + nogil  # wrap-doc:Returns the enzyme with the given name (supports synonym names)
        const DigestionEnzymeRNA* getEnzymeByRegEx(const String& cleavage_regex) except + nogil  # wrap-doc:Returns the enzyme with the given cleavage regex
        void getAllNames(libcpp_vector[ String ]& all_names) except + nogil  # wrap-doc:Returns all enzyme names (does not include synonym names)
        bool hasEnzyme(const String& name) except + nogil  # wrap-doc:Returns True if the database contains an enzyme with the given name
        # bool hasRegEx(const String& cleavage_regex) except + nogil # We don't use regexes for RNA
        # bool hasEnzyme(DigestionEnzymeRNA* enzyme) except + nogil  # does not make sense as the ptr wont match

        # ConstEnzymeIterator beginEnzyme() except + nogil
        # ConstEnzymeIterator endEnzyme() except + nogil

## wrap static methods
cdef extern from "<OpenMS/CHEMISTRY/RNaseDB.h>" namespace "OpenMS::RNaseDB":

    RNaseDB* getInstance() except + nogil  # wrap-ignore
