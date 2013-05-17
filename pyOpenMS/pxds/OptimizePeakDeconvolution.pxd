from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from PeakShape cimport *
from OptimizePick cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePeakDeconvolution.h>" namespace "OpenMS":
    
    cdef cppclass OptimizePeakDeconvolution(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        OptimizePeakDeconvolution() nogil except +
        OptimizePeakDeconvolution(OptimizePeakDeconvolution) nogil except +
        # NAMESPACE # OptimizationFunctions::PenaltyFactorsIntensity  getPenalties() nogil except +
        # NAMESPACE # void setPenalties(OptimizationFunctions::PenaltyFactorsIntensity & penalties) nogil except +
        Int getCharge() nogil except +
        void setCharge(Int charge) nogil except +
        bool optimize(libcpp_vector[ PeakShape ] & peaks, OptimizePeakDeconvolution_Data & data) nogil except +


cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePeakDeconvolution.h>" namespace "OpenMS::OptimizationFunctions":
    
    cdef cppclass PenaltyFactorsIntensity(PenaltyFactors) :
        # wrap-inherits:
        #  PenaltyFactors
        PenaltyFactorsIntensity() nogil except +
        PenaltyFactorsIntensity(PenaltyFactorsIntensity) nogil except +
        DoubleReal height


cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePeakDeconvolution.h>" namespace "OpenMS::OptimizePeakDeconvolution":

    cdef cppclass OptimizePeakDeconvolution_Data "OpenMS::OptimizePeakDeconvolution::Data"::
      libcpp_vector[PeakShape] peaks
      libcpp_vector[double] positions
      libcpp_vector[double] signal
      PenaltyFactorsIntensity penalties 
      Int charge 

