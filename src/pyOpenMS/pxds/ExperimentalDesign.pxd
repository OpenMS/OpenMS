from Types cimport *
from String cimport *
from MSRun cimport *

cdef extern from "<OpenMS/METADATA/ExperimentalDesign.h>" namespace "OpenMS":
    
    cdef cppclass ExperimentalDesign "OpenMS::ExperimentalDesign":
        ExperimentalDesign(ExperimentalDesign) nogil except + #wrap-ignore
        libcpp_vector[ MSRun ] runs
        # libcpp_map[ unsigned, libcpp_set[ unsigned ] ] getFractionToRunsMapping() nogil except +
        bool sameNrOfRunsPerFraction() nogil except +
        void load(String & tsv_file, ExperimentalDesign & design) nogil except +

cdef extern from "<OpenMS/METADATA/ExperimentalDesign.h>" namespace "OpenMS::ExperimentalDesign":
    
    cdef cppclass MSRun "OpenMS::ExperimentalDesign::MSRun":
        # wrap-attach:
        #    ExperimentalDesign
        MSRun() nogil except +
        MSRun(MSRun) nogil except + #wrap-ignore
        libcpp_string file
        unsigned fraction
        unsigned technical_replicate

