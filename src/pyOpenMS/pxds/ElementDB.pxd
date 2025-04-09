from Types cimport *
from libcpp cimport bool
from String cimport *
from Element cimport *
from IsotopeDistribution cimport *

cdef extern from "<OpenMS/CHEMISTRY/ElementDB.h>" namespace "OpenMS":
    
    cdef cppclass ElementDB "OpenMS::ElementDB":
        # wrap-manual-memory:
        #   cdef AutowrapPtrHolder[_ElementDB] inst

        # private
        ElementDB() except + nogil  # wrap-ignore
        # private 
        ElementDB(ElementDB) except + nogil  #wrap-ignore

        # No wrapping of const ref
        # const Map[ String, Element * ]  getNames() except + nogil 
        # const Map[ String, Element * ] getSymbols() except + nogil 
        # const Map[unsigned int, Element * ] getAtomicNumbers() except + nogil 
        const Element * getElement(const String & name) except + nogil 
        const Element * getElement(UInt atomic_number) except + nogil 
        void addElement(libcpp_string name, libcpp_string symbol,
                        unsigned int an,
                        libcpp_map[unsigned int, double] abundance,
                        libcpp_map[unsigned int, double] mass,
                        bool replace_existing) except + nogil 
        bool hasElement(const String & name) except + nogil  # wrap-doc:Returns true if the db contains an element with the given name, else false
        bool hasElement(UInt atomic_number) except + nogil  # wrap-doc:Returns true if the db contains an element with the given atomic_number, else false

## wrap static methods
cdef extern from "<OpenMS/CHEMISTRY/ElementDB.h>" namespace "OpenMS::ElementDB":
    
    ElementDB* getInstance() except + nogil  # wrap-ignore

