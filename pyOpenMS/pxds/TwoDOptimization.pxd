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
        TwoDOptimization(TwoDOptimization) nogil except +
        double getMZTolerance() nogil except +
        void setMZTolerance(double tolerance_mz) nogil except +
        double getMaxPeakDistance() nogil except +
        void setMaxPeakDistance(double max_peak_distance) nogil except +
        UInt getMaxIterations() nogil except +
        void setMaxIterations(UInt max_iteration) nogil except +
        # NAMESPACE # OptimizationFunctions::PenaltyFactorsIntensity  getPenalties() nogil except +
        # NAMESPACE # void setPenalties(OptimizationFunctions::PenaltyFactorsIntensity & penalties) nogil except +
        # TEMPLATE # void optimize(InputSpectrumIterator first, InputSpectrumIterator last, MSExperiment[ OutputPeakType ] & ms_exp, bool real2D) nogil except +

