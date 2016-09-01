from Types cimport *
from String cimport *
from Enzyme cimport *

cdef extern from "<OpenMS/CHEMISTRY/EnzymesDB.h>" namespace "OpenMS":
    
    cdef cppclass EnzymesDB "OpenMS::EnzymesDB":
        EnzymesDB(EnzymesDB) nogil except + #wrap-ignore

        # const Enzyme * getEnzyme(String & name) nogil except +
        # const Enzyme * getEnzymeByRegEx(String & cleavage_regex) nogil except +
        void setEnzymes(String & filename) nogil except +
        void addEnzyme(Enzyme & enzyme) nogil except +
        void clear() nogil except +
        void getAllNames(libcpp_vector[ String ] & all_names) nogil except +
        void getAllXTandemNames(libcpp_vector[ String ] & all_names) nogil except +
        void getAllOMSSANames(libcpp_vector[ String ] & all_names) nogil except +
        bool hasEnzyme(String & name) nogil except +
        bool hasRegEx(String & cleavage_regex) nogil except +
        bool hasEnzyme(Enzyme * enzyme) nogil except +
        # ConstEnzymeIterator beginEnzyme() nogil except +
        # ConstEnzymeIterator endEnzyme() nogil except +
        # POINTER # EnzymesDB * getInstance() nogil except +

