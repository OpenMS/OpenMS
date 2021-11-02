from Types cimport *
from libcpp cimport bool
from String cimport *
from Map cimport *
from Element cimport *
from IsotopeDistribution cimport *

cdef extern from "<OpenMS/CHEMISTRY/ElementDB.h>" namespace "OpenMS":
    
    cdef cppclass ElementDB "OpenMS::ElementDB":
        # wrap-manual-memory:
        #    cdef AutowrapConstPtrHolder[_ElementDB] inst

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
        bool hasElement(const String & name) nogil except + # wrap-doc:Returns true if the db contains an element with the given name, else false
        bool hasElement(UInt atomic_number) nogil except + # wrap-doc:Returns true if the db contains an element with the given atomic_number, else false

## wrap static methods
cdef extern from "<OpenMS/CHEMISTRY/ElementDB.h>" namespace "OpenMS::ElementDB":
    
    const ElementDB* getInstance() nogil except + # wrap-ignore

