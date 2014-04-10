from Types cimport *
from libcpp cimport bool
from libcpp.pair cimport pair as libcpp_pair
# from TraceFitter cimport *
from Peak1D cimport *
from FeatureFinderAlgorithmPickedHelperStructs cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussTraceFitter.h>" namespace "OpenMS":
    
    cdef cppclass GaussTraceFitter[PeakType]:
        # wrap-instances:
        #   GaussTraceFitter := GaussTraceFitter[Peak1D]

        # possible inheritance ..
        #  OpenMS::TraceFitter< PeakType >
        GaussTraceFitter() nogil except +
        GaussTraceFitter(GaussTraceFitter) nogil except +
        void fit(MassTraces[ PeakType ] & traces) nogil except +
        double getLowerRTBound() nogil except +
        double getUpperRTBound() nogil except +
        double getHeight() nogil except +
        double getCenter() nogil except +
        double getFWHM() nogil except +
        double getSigma() nogil except +
        bool checkMaximalRTSpan(double max_rt_span) nogil except +
        bool checkMinimalRTSpan(libcpp_pair[ double, double ] & rt_bounds, double min_rt_span) nogil except +
        double computeTheoretical(MassTrace[ Peak1D ] & trace, Size k) nogil except +
        double getArea() nogil except +
        String getGnuplotFormula(MassTrace[ Peak1D ] & trace, char function_name, double baseline, double rt_shift) nogil except +

