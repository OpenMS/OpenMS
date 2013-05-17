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

