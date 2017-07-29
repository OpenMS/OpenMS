from Types cimport *
from String cimport *
from Enzyme cimport *

cdef extern from "<OpenMS/CHEMISTRY/EnzymesDB.h>" namespace "OpenMS":
    
    cdef cppclass EnzymesDB "OpenMS::EnzymesDB":
        # wrap-manual-memory

        EnzymesDB() nogil except + #wrap-ignore
        EnzymesDB(EnzymesDB) nogil except + #wrap-ignore

        const Enzyme * getEnzyme(String & name) nogil except +
        const Enzyme * getEnzymeByRegEx(String & cleavage_regex) nogil except +
        void setEnzymes(String & filename) nogil except +
        void addEnzyme(Enzyme & enzyme) nogil except +
        void clear() nogil except +
        void getAllNames(libcpp_vector[ String ] & all_names) nogil except + 
        void getAllXTandemNames(libcpp_vector[ String ] & all_names) nogil except +
        void getAllOMSSANames(libcpp_vector[ String ] & all_names) nogil except +
        void getAllCometNames(libcpp_vector[ String ] & all_names) nogil except +
        bool hasEnzyme(String & name) nogil except +
        bool hasRegEx(String & cleavage_regex) nogil except +
        # bool hasEnzyme(Enzyme * enzyme) nogil except + # does not make sense as the ptr wont match

        # ConstEnzymeIterator beginEnzyme() nogil except +
        # ConstEnzymeIterator endEnzyme() nogil except +

## wrap static methods
cdef extern from "<OpenMS/CHEMISTRY/EnzymesDB.h>" namespace "OpenMS::EnzymesDB":
    
    EnzymesDB* getInstance() nogil except + # wrap-ignore

