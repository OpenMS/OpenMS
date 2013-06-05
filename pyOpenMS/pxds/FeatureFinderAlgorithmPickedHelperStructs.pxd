from libcpp.pair cimport pair as libcpp_pair
from Types cimport *
from String cimport *
from StringList cimport *
from ConvexHull2D cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>" namespace "OpenMS::FeatureFinderAlgorithmPickedHelperStructs":
    
    cdef cppclass TheoreticalIsotopePattern "OpenMS::FeatureFinderAlgorithmPickedHelperStructs::TheoreticalIsotopePattern":
        TheoreticalIsotopePattern() nogil except + # wrap-ignore
        TheoreticalIsotopePattern(TheoreticalIsotopePattern) nogil except + #wrap-ignore

        libcpp_vector[ double ] intensity
        Size optional_begin
        Size optional_end
        DoubleReal max
        Size trimmed_left

        Size size() nogil except +

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>" namespace "OpenMS::FeatureFinderAlgorithmPickedHelperStructs":
    
    # Since this is a templated class, we cannot tell Cython what the C++
    # equivalent would be and we need to name it MassTrace
    cdef cppclass MassTrace[PeakType]:
        # wrap-instances:
        #   MassTrace := MassTrace[Peak1D]
        MassTrace(MassTrace) nogil except + #wrap-ignore
        # POINTER # PeakType * max_peak
        DoubleReal max_rt
        DoubleReal theoretical_int
        # POINTER # libcpp_vector[ libcpp_pair[ DoubleReal, PeakType * ] ] peaks
        ConvexHull2D getConvexhull() nogil except +
        void updateMaximum() nogil except +
        DoubleReal getAvgMZ() nogil except +
        bool isValid() nogil except +

    # Since this is a templated class, we cannot tell Cython what the C++
    # equivalent would be and we need to name it MassTraces
    cdef cppclass MassTraces[PeakType]:
        # wrap-instances:
        #   MassTraces := MassTraces[Peak1D]
        MassTraces() nogil except +
        MassTraces(MassTraces) nogil except + #wrap-ignore
        Size max_trace
        DoubleReal baseline
        Size getPeakCount() nogil except +
        bool isValid(DoubleReal seed_mz, DoubleReal trace_tolerance) nogil except +
        Size getTheoreticalmaxPosition() nogil except +
        void updateBaseline() nogil except +
        libcpp_pair[ DoubleReal, DoubleReal ] getRTBounds() nogil except +

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
        # TODO STL attributes -- Signed size does not work either!
        # vector.from_py:33:13: 'ptrdiff_t' is not a type identifier
        # libcpp_vector[ SignedSize ] peak
        # libcpp_vector[ size_t ] spectrum
        # libcpp_vector[ double ] intensity
        # libcpp_vector[ double ] mz_score
        # libcpp_vector[ double ] theoretical_mz
        TheoreticalIsotopePattern theoretical_pattern
        IsotopePattern(Size size) nogil except +


