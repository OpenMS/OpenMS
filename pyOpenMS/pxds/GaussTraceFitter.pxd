from Types cimport *
from libcpp cimport bool
from libcpp.pair cimport pair as libcpp_pair
# from TraceFitter cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussTraceFitter.h>" namespace "OpenMS":
    
    cdef cppclass GaussTraceFitter[PeakType]:
        # wrap-instances:
        #   GaussTraceFitter := GaussTraceFitter[Peak1D]

        # possible inheritance ..
        #  OpenMS::TraceFitter< PeakType >
        GaussTraceFitter() nogil except +
        GaussTraceFitter(GaussTraceFitter) nogil except +
        # TODO
        # NAMESPACE # void fit(FeatureFinderAlgorithmPickedHelperStructs::MassTraces[ PeakType ] & traces) nogil except +
        DoubleReal getLowerRTBound() nogil except +
        DoubleReal getUpperRTBound() nogil except +
        DoubleReal getHeight() nogil except +
        DoubleReal getCenter() nogil except +
        DoubleReal getFWHM() nogil except +
        DoubleReal getSigma() nogil except +
        bool checkMaximalRTSpan(DoubleReal max_rt_span) nogil except +
        bool checkMinimalRTSpan(libcpp_pair[ DoubleReal, DoubleReal ] & rt_bounds, DoubleReal min_rt_span) nogil except +
        # NAMESPACE # DoubleReal computeTheoretical(FeatureFinderAlgorithmPickedHelperStructs::MassTrace[ PeakType ] & trace, Size k) nogil except +
        DoubleReal getFeatureIntensityContribution() nogil except +
        # NAMESPACE # String getGnuplotFormula(FeatureFinderAlgorithmPickedHelperStructs::MassTrace[ PeakType ] & trace, char function_name, DoubleReal baseline, DoubleReal rt_shift) nogil except +

