from Types cimport *
from libcpp cimport bool
from PeakShape cimport *
# from OptimizePeakDeconvolution cimport *
from MSExperiment cimport *
from MSSpectrum cimport *
# from PeakIndex cimport *
from IsotopeCluster cimport *
from DefaultParamHandler cimport *
# from OptimizePick cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/TwoDOptimization.h>" namespace "OpenMS":
    
    cdef cppclass TwoDOptimization(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        TwoDOptimization() nogil except +
        TwoDOptimization(TwoDOptimization &) nogil except +

        double getMZTolerance() nogil except + # wrap-doc:Returns the matching epsilon
        void setMZTolerance(double tolerance_mz) nogil except + # wrap-doc:Sets the matching epsilon
        double getMaxPeakDistance() nogil except + # wrap-doc:Returns the maximal peak distance in a cluster
        void setMaxPeakDistance(double max_peak_distance) nogil except + # wrap-doc:Sets the maximal peak distance in a cluster
        UInt getMaxIterations() nogil except + # wrap-doc:Returns the maximal number of iterations
        void setMaxIterations(UInt max_iteration) nogil except + # wrap-doc:Sets the maximal number of iterations
        # NAMESPACE # OptimizationFunctions::PenaltyFactorsIntensity  getPenalties() nogil except +
        # NAMESPACE # void setPenalties(OptimizationFunctions::PenaltyFactorsIntensity & penalties) nogil except +
        # TEMPLATE # void optimize(InputSpectrumIterator first, InputSpectrumIterator last, MSExperiment & ms_exp, bool real2D) nogil except +
