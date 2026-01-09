from Types cimport *
from String cimport *
from Ribonucleotide cimport *
from libcpp.string cimport string as libcpp_utf8_string

cdef extern from "<OpenMS/CHEMISTRY/RibonucleotideDB.h>" namespace "OpenMS":
    
    cdef cppclass RibonucleotideDB "OpenMS::RibonucleotideDB":
        # wrap-manual-memory:
        #  cdef AutowrapPtrHolder[_RibonucleotideDB] inst


        # protected
        RibonucleotideDB() except + nogil  #wrap-ignore
        # deleted
        RibonucleotideDB(RibonucleotideDB) except + nogil  #wrap-ignore

        const Ribonucleotide * getRibonucleotide(const libcpp_utf8_string& code) except + nogil
        const Ribonucleotide * getRibonucleotidePrefix(const libcpp_utf8_string& code) except + nogil
        libcpp_pair[const Ribonucleotide *, const Ribonucleotide *] getRibonucleotideAlternatives(const libcpp_utf8_string& code) except + nogil  # wrap-ignore

# COMMENT: wrap static methods
cdef extern from "<OpenMS/CHEMISTRY/RibonucleotideDB.h>" namespace "OpenMS::RibonucleotideDB":
    
    RibonucleotideDB* getInstance() except + nogil  # wrap-ignore
