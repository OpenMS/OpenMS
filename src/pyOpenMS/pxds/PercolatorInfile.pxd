from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from libcpp.string cimport string as libcpp_utf8_string
from PeptideIdentification cimport *
from ProteinIdentification cimport *
from StringList cimport *

cdef extern from "<OpenMS/FORMAT/PercolatorInfile.h>" namespace "OpenMS":
    
    cdef cppclass PercolatorInfile "OpenMS::PercolatorInfile":
        # wrap-doc:
            #   Class for storing Percolator tab-delimited input files

        PercolatorInfile() nogil except +
        PercolatorInfile(PercolatorInfile &) nogil except + 


# COMMENT: wrap static methods
cdef extern from "<OpenMS/FORMAT/PercolatorInfile.h>" namespace "OpenMS::PercolatorInfile":
        
        # static members
        void store(String pin_file, libcpp_vector[PeptideIdentification] peptide_ids, StringList feature_set, libcpp_string, int min_charge, int max_charge) nogil except +  # wrap-attach:PercolatorInfile
        
