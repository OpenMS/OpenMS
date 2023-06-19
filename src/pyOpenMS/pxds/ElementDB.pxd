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
        ElementDB() nogil except + # wrap-ignore
        # private 
        ElementDB(ElementDB) nogil except + #wrap-ignore

        # No wrapping of const ref
        # const Map[ String, Element * ]  getNames() nogil except +
        # const Map[ String, Element * ] getSymbols() nogil except +
        # const Map[unsigned int, Element * ] getAtomicNumbers() nogil except +
        const Element * getElement(const String & name) nogil except +
        const Element * getElement(UInt atomic_number) nogil except +
        void addElement(libcpp_string name, libcpp_string symbol,
                        unsigned int an,
                        libcpp_map[unsigned int, double] abundance,
                        libcpp_map[unsigned int, double] mass,
                        bool replace_existing) nogil except +
        bool hasElement(const String & name) nogil except + # wrap-doc:Returns true if the db contains an element with the given name, else false
        bool hasElement(UInt atomic_number) nogil except + # wrap-doc:Returns true if the db contains an element with the given atomic_number, else false

## wrap static methods
cdef extern from "<OpenMS/CHEMISTRY/ElementDB.h>" namespace "OpenMS::ElementDB":
    
    ElementDB* getInstance() nogil except + # wrap-ignore

