from Types cimport *
from libcpp cimport bool
from libcpp.pair cimport pair as libcpp_pair
# from TraceFitter cimport *
from Peak1D cimport *
from FeatureFinderAlgorithmPickedHelperStructs cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussTraceFitter.h>" namespace "OpenMS":

    cdef cppclass GaussTraceFitter:
        GaussTraceFitter() nogil except +
        GaussTraceFitter(GaussTraceFitter) nogil except +
        void fit(MassTraces& traces) nogil except +
        double getLowerRTBound() nogil except +
        double getUpperRTBound() nogil except +
        double getHeight() nogil except +
        double getCenter() nogil except +
        double getFWHM() nogil except +
        double getSigma() nogil except +
        bool checkMaximalRTSpan(double max_rt_span) nogil except +
        bool checkMinimalRTSpan(libcpp_pair[ double, double ] & rt_bounds, double min_rt_span) nogil except +
        double computeTheoretical(MassTrace & trace, Size k) nogil except +
        double getArea() nogil except +
        String getGnuplotFormula(MassTrace & trace, char function_name, double baseline, double rt_shift) nogil except +
        double getValue(double rt) nogil except +

