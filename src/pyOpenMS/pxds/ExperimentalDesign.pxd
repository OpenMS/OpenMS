from Types cimport *
from String cimport *

cdef extern from "<OpenMS/METADATA/ExperimentalDesign.h>" namespace "OpenMS":
    
    cdef cppclass ExperimentalDesign "OpenMS::ExperimentalDesign":
        ExperimentalDesign(ExperimentalDesign) nogil except + #wrap-ignore
        libcpp_vector[ ExperimentalDesign_MSRun] runs
        # libcpp_map[ unsigned, libcpp_set[unsigned] ] getFractionToRunsMapping() nogil except +
        bool sameNrOfRunsPerFraction() nogil except +
        void load(const String & tsv_file, ExperimentalDesign & design) nogil except +

cdef extern from "<OpenMS/METADATA/ExperimentalDesign.h>" namespace "OpenMS::ExperimentalDesign":
    
    cdef cppclass ExperimentalDesign_MSRun "OpenMS::ExperimentalDesign::MSRun":
        ExperimentalDesign_MSRun() nogil except +
        ExperimentalDesign_MSRun(ExperimentalDesign_MSRun) nogil except + #wrap-ignore
        libcpp_string file
        unsigned fraction
        unsigned technical_replicate

