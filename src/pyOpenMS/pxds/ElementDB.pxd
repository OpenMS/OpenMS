from Types cimport *
from libcpp cimport bool
from String cimport *
from Map cimport *
from Element cimport *
from IsotopeDistribution cimport *

cdef extern from "<OpenMS/CHEMISTRY/ElementDB.h>" namespace "OpenMS":
    
    cdef cppclass ElementDB "OpenMS::ElementDB":
        ElementDB(ElementDB) nogil except + #wrap-ignore
        # Cannot get ElementDB since its destructor is private ... 
        # POINTER # ElementDB * getInstance() nogil except +
        # POINTER # Map[ String, Element * ]  getNames() nogil except +
        # POINTER # Map[ String, Element * ]  getSymbols() nogil except +
        # POINTER # Map[ UInt, Element * ]  getAtomicNumbers() nogil except +
        # CONST POINTER # Element * getElement(String & name) nogil except +
        # CONST POINTER # Element * getElement(UInt atomic_number) nogil except +
        bool hasElement(String & name) nogil except +
        bool hasElement(UInt atomic_number) nogil except +

