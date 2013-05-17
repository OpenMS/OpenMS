from Types cimport *
from Types cimport *
from String cimport *
from StringList cimport *
from ConvexHull2D cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>" namespace "OpenMS::FeatureFinderAlgorithmPickedHelperStructs":
    
    cdef cppclass TheoreticalIsotopePattern "OpenMS::FeatureFinderAlgorithmPickedHelperStructs::TheoreticalIsotopePattern":
        TheoreticalIsotopePattern() nogil except + # wrap-ignore
        TheoreticalIsotopePattern(TheoreticalIsotopePattern) nogil except + #wrap-ignore

        # libcpp_vector[ double ] intensity
        Size optional_begin
        Size optional_end
        DoubleReal max
        Size trimmed_left

        Size size() nogil except +

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>" namespace "OpenMS::FeatureFinderAlgorithmPickedHelperStructs":
    
    cdef cppclass MassTrace[PeakType]:
        MassTrace(MassTrace) nogil except + #wrap-ignore
        # POINTER # PeakType * max_peak
        DoubleReal max_rt
        DoubleReal theoretical_int
        # POINTER # libcpp_vector[ libcpp_pair[ DoubleReal, PeakType * ] ] peaks
        ConvexHull2D getConvexhull() nogil except +
        void updateMaximum() nogil except +
        DoubleReal getAvgMZ() nogil except +
        bool isValid() nogil except +

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>" namespace "OpenMS::FeatureFinderAlgorithmPickedHelperStructs":
    
    cdef cppclass Seed "OpenMS::FeatureFinderAlgorithmPickedHelperStructs::Seed":
        Seed(Seed) nogil except + #wrap-ignore
        Size spectrum
        Size peak
        Real intensity
        bool operator<(Seed & rhs) nogil except +

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>" namespace "OpenMS::FeatureFinderAlgorithmPickedHelperStructs":
    
    cdef cppclass IsotopePattern "OpenMS::FeatureFinderAlgorithmPickedHelperStructs::IsotopePattern":
        IsotopePattern(IsotopePattern) nogil except + #wrap-ignore
        # TODO 
        libcpp_vector[ ptrdiff_t ] peak
        # libcpp_vector[ size_t ] spectrum
        # libcpp_vector[ double ] intensity
        # libcpp_vector[ double ] mz_score
        # libcpp_vector[ double ] theoretical_mz
        TheoreticalIsotopePattern theoretical_pattern
        IsotopePattern(Size size) nogil except +

