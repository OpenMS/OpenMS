from smart_ptr cimport shared_ptr
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from libcpp.pair cimport pair as libcpp_pair
from TheoreticalSpectrumGenerator cimport *
from SimTypes cimport *
from SVMWrapper cimport *
from IonType cimport *

cdef extern from "<OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGenerator.h>" namespace "OpenMS::SvmTheoreticalSpectrumGenerator":
    
    cdef cppclass SvmModelParameterSet "OpenMS::SvmTheoreticalSpectrumGenerator::SvmModelParameterSet":
        SvmModelParameterSet(SvmModelParameterSet) nogil except + #wrap-ignore
        # libcpp_vector[ shared_ptr[ SVMWrapper ] ] class_models
        # libcpp_vector[ shared_ptr[ SVMWrapper ] ] reg_models
        # TODO STL map with wrapped key
        # libcpp_map[ ResidueType, double ] _intensities
        # libcpp_vector[ IonType ] ion_types
        # TODO STL map with wrapped key
        # libcpp_map[ IonType, libcpp_vector[ IonType ] ] secondary_types
        Size number_intensity_levels
        Size number_regions
        libcpp_vector[ double ] feature_max
        libcpp_vector[ double ] feature_min
        double scaling_lower
        double scaling_upper
        libcpp_vector[ double ] intensity_bin_boarders
        libcpp_vector[ double ] intensity_bin_values
        # TODO nested STL
        # libcpp_map[ libcpp_pair[ IonType, Size ], libcpp_vector[ libcpp_vector[ double ] ] ] conditional_prob

