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
        TwoDOptimization() except + nogil 
        TwoDOptimization(TwoDOptimization &) except + nogil 

        double getMZTolerance() except + nogil  # wrap-doc:Returns the matching epsilon
        void setMZTolerance(double tolerance_mz) except + nogil  # wrap-doc:Sets the matching epsilon
        double getMaxPeakDistance() except + nogil  # wrap-doc:Returns the maximal peak distance in a cluster
        void setMaxPeakDistance(double max_peak_distance) except + nogil  # wrap-doc:Sets the maximal peak distance in a cluster
        UInt getMaxIterations() except + nogil  # wrap-doc:Returns the maximal number of iterations
        void setMaxIterations(UInt max_iteration) except + nogil  # wrap-doc:Sets the maximal number of iterations
        # NAMESPACE # OptimizationFunctions::PenaltyFactorsIntensity  getPenalties() except + nogil 
        # NAMESPACE # void setPenalties(OptimizationFunctions::PenaltyFactorsIntensity & penalties) except + nogil 
        # TEMPLATE # void optimize(InputSpectrumIterator first, InputSpectrumIterator last, MSExperiment & ms_exp, bool real2D) except + nogil 
