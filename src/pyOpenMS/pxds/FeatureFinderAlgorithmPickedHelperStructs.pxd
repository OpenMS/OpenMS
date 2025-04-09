from libcpp.pair cimport pair as libcpp_pair
from Types cimport *
from String cimport *
from StringList cimport *
from ConvexHull2D cimport *

cdef extern from "<OpenMS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>" namespace "OpenMS::FeatureFinderAlgorithmPickedHelperStructs":

    cdef cppclass TheoreticalIsotopePattern "OpenMS::FeatureFinderAlgorithmPickedHelperStructs::TheoreticalIsotopePattern":
        TheoreticalIsotopePattern() except + nogil  # compiler
        TheoreticalIsotopePattern(TheoreticalIsotopePattern &) except + nogil  # compiler

        libcpp_vector[ double ] intensity
        Size optional_begin
        Size optional_end
        double max
        Size trimmed_left

        Size size() except + nogil 

cdef extern from "<OpenMS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>" namespace "OpenMS::FeatureFinderAlgorithmPickedHelperStructs":

    # Since this is a templated class, we cannot tell Cython what the C++
    # equivalent would be and we need to name it MassTrace
    cdef cppclass MassTrace:
        MassTrace() except + nogil  # compiler
        MassTrace(MassTrace &) except + nogil  # compiler

        # POINTER # PeakType * max_peak
        double max_rt
        double theoretical_int
        # POINTER # libcpp_vector[ libcpp_pair[ double, PeakType * ] ] peaks
        ConvexHull2D getConvexhull() except + nogil 
        void updateMaximum() except + nogil 
        double getAvgMZ() except + nogil 
        bool isValid() except + nogil 

    # Since this is a templated class, we cannot tell Cython what the C++
    # equivalent would be and we need to name it MassTraces
    cdef cppclass MassTraces:
        MassTraces() except + nogil  # compiler
        MassTraces(MassTraces &) except + nogil  # compiler

        Size max_trace
        double baseline
        Size getPeakCount() except + nogil 
        bool isValid(double seed_mz, double trace_tolerance) except + nogil 
        Size getTheoreticalmaxPosition() except + nogil 
        void updateBaseline() except + nogil 
        libcpp_pair[ double, double ] getRTBounds() except + nogil 

cdef extern from "<OpenMS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>" namespace "OpenMS::FeatureFinderAlgorithmPickedHelperStructs":

    cdef cppclass Seed "OpenMS::FeatureFinderAlgorithmPickedHelperStructs::Seed":
        Seed() except + nogil  # compiler
        Seed(Seed &) except + nogil  # compiler
        Size spectrum
        Size peak
        float intensity
        bool operator<(Seed & rhs) except + nogil 

cdef extern from "<OpenMS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>" namespace "OpenMS::FeatureFinderAlgorithmPickedHelperStructs":

    cdef cppclass IsotopePattern "OpenMS::FeatureFinderAlgorithmPickedHelperStructs::IsotopePattern":
        IsotopePattern(Size size) except + nogil 
        IsotopePattern(IsotopePattern &) except + nogil  # compiler

        # TODO STL attributes -- Signed size does not work either!
        # vector.from_py:33:13: 'ptrdiff_t' is not a type identifier
        # libcpp_vector[ SignedSize ] peak
        libcpp_vector[ size_t ] spectrum
        libcpp_vector[ double ] intensity
        libcpp_vector[ double ] mz_score
        libcpp_vector[ double ] theoretical_mz
        TheoreticalIsotopePattern theoretical_pattern



