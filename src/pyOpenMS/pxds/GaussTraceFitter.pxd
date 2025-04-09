from Types cimport *
from libcpp cimport bool
from libcpp.pair cimport pair as libcpp_pair
# from TraceFitter cimport *
from Peak1D cimport *
from FeatureFinderAlgorithmPickedHelperStructs cimport *

cdef extern from "<OpenMS/FEATUREFINDER/GaussTraceFitter.h>" namespace "OpenMS":

    cdef cppclass GaussTraceFitter:

        GaussTraceFitter() except + nogil  # wrap-doc:Fitter for RT profiles using a Gaussian background model
        GaussTraceFitter(GaussTraceFitter &) except + nogil 
        void fit(MassTraces& traces) except + nogil  # wrap-doc:Override important methods
        double getLowerRTBound() except + nogil  # wrap-doc:Returns the lower RT bound
        double getUpperRTBound() except + nogil  # wrap-doc:Returns the upper RT bound
        double getHeight() except + nogil  # wrap-doc:Returns height of the fitted gaussian model
        double getCenter() except + nogil  # wrap-doc:Returns center of the fitted gaussian model
        double getFWHM() except + nogil  # wrap-doc:Returns FWHM of the fitted gaussian model
        double getSigma() except + nogil  # wrap-doc:Returns Sigma of the fitted gaussian model
        bool checkMaximalRTSpan(double max_rt_span) except + nogil  # TODO
        bool checkMinimalRTSpan(libcpp_pair[ double, double ] & rt_bounds, double min_rt_span) except + nogil  # TODO
        double computeTheoretical(MassTrace & trace, Size k) except + nogil  # TODO
        double getArea() except + nogil  # wrap-doc:Returns area of the fitted gaussian model
        String getGnuplotFormula(MassTrace & trace, char function_name, double baseline, double rt_shift) except + nogil  # TODO
        double getValue(double rt) except + nogil  # wrap-doc:Returns value of the fitted gaussian model
