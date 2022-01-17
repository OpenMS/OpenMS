from Types cimport *
from libcpp cimport bool
from libcpp.pair cimport pair as libcpp_pair
# from TraceFitter cimport *
from Peak1D cimport *
from FeatureFinderAlgorithmPickedHelperStructs cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussTraceFitter.h>" namespace "OpenMS":

    cdef cppclass GaussTraceFitter:

        GaussTraceFitter() nogil except + # wrap-doc:Fitter for RT profiles using a Gaussian background model
        GaussTraceFitter(GaussTraceFitter &) nogil except +
        void fit(MassTraces& traces) nogil except + # wrap-doc:Override important methods
        double getLowerRTBound() nogil except + # wrap-doc:Returns the lower RT bound
        double getUpperRTBound() nogil except + # wrap-doc:Returns the upper RT bound
        double getHeight() nogil except + # wrap-doc:Returns height of the fitted gaussian model
        double getCenter() nogil except + # wrap-doc:Returns center of the fitted gaussian model
        double getFWHM() nogil except + # wrap-doc:Returns FWHM of the fitted gaussian model
        double getSigma() nogil except + # wrap-doc:Returns Sigma of the fitted gaussian model
        bool checkMaximalRTSpan(double max_rt_span) nogil except + # TODO
        bool checkMinimalRTSpan(libcpp_pair[ double, double ] & rt_bounds, double min_rt_span) nogil except + # TODO
        double computeTheoretical(MassTrace & trace, Size k) nogil except + # TODO
        double getArea() nogil except + # wrap-doc:Returns area of the fitted gaussian model
        String getGnuplotFormula(MassTrace & trace, char function_name, double baseline, double rt_shift) nogil except + # TODO
        double getValue(double rt) nogil except + # wrap-doc:Returns value of the fitted gaussian model
