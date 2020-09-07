from Types cimport *
from libcpp.map cimport map as libcpp_map
from libcpp.string cimport string as libcpp_string
from libcpp.vector cimport vector as libcpp_vector
from String cimport *

cdef extern from "<OpenMS/FORMAT/OSWFile.h>" namespace "OpenMS":
    
    cdef cppclass OSWFile "OpenMS::OSWFile":
        OSWFile() nogil except +
        OSWFile(OSWFile) nogil except + #wrap-ignore

        # Cannot wrap string stream
        # void read(const libcpp_string & in_osw,
        #          const libcpp_string & osw_level,
        #          libcpp_stringstream & pin_output,
        #          const double & ipf_max_peakgroup_pep,
        #          const double & ipf_max_transition_isotope_overlap,
        #          const double & ipf_min_transition_sn) nogil except +

        # NESTED STL # void write(const libcpp_string & in_osw, const libcpp_string & osw_level, const libcpp_map[ libcpp_string, libcpp_vector[ double ] ] & features) nogil except +

