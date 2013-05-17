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

