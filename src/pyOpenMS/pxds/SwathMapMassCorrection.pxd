from Types cimport *
from libcpp cimport bool
from libcpp.string cimport string as libcpp_string
from libcpp.vector cimport vector as libcpp_vector
from SwathMap cimport *
from MRMFeatureFinderScoring cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/SwathMapMassCorrection.h>" namespace "OpenMS":
    
    cdef cppclass SwathMapMassCorrection "OpenMS::SwathMapMassCorrection":
        SwathMapMassCorrection(SwathMapMassCorrection) nogil except + #wrap-ignore
        # NAMESPACE # void correctMZ(const OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType & transition_group_map, libcpp_vector[ OpenSwath::SwathMap ] & swath_maps, const libcpp_string & corr_type, const double mz_extr_window, const bool ppm) nogil except +

