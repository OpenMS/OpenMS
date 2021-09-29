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
        OptimizePeakDeconvolution(OptimizePeakDeconvolution &) nogil except +
        PenaltyFactorsIntensity  getPenalties() nogil except + # wrap-doc:Return the penalty parameter
        void setPenalties(PenaltyFactorsIntensity & penalties) nogil except + # wrap-doc:Set the penalty parameter
        Int getCharge() nogil except + # wrap-doc:Return the charge
        void setCharge(Int charge) nogil except + # wrap-doc:Set the charge
        bool optimize(libcpp_vector[ PeakShape ] & peaks, OptimizePeakDeconvolution_Data & data) nogil except + # wrap-doc:Performs a nonlinear optimization of the peaks that belong to the current isotope pattern
        Size getNumberOfPeaks_(Int charge, libcpp_vector[ PeakShape ] & temp_shapes, OptimizePeakDeconvolution_Data & data) nogil except +

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePeakDeconvolution.h>" namespace "OpenMS::OptimizationFunctions":
    
    cdef cppclass PenaltyFactorsIntensity(OptimizationFunctions_PenaltyFactors):
        # wrap-inherits:
        #  OptimizationFunctions_PenaltyFactors
        # wrap-doc:
        #   Class for the penalty factors used during the optimization
        #   -----
        #   A great deviation (squared deviation) of a peak shape's position or its left or right width parameter can be penalised
        #   During the optimization negative heights may occur, they are penalised as well

        PenaltyFactorsIntensity() nogil except +
        PenaltyFactorsIntensity(PenaltyFactorsIntensity) nogil except +
        double height


cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePeakDeconvolution.h>" namespace "OpenMS::OptimizePeakDeconvolution":

    cdef cppclass OptimizePeakDeconvolution_Data "OpenMS::OptimizePeakDeconvolution::Data":
      OptimizePeakDeconvolution_Data() nogil except +
      OptimizePeakDeconvolution_Data(OptimizePeakDeconvolution_Data) nogil except + # no-wrap

      libcpp_vector[PeakShape] peaks
      libcpp_vector[double] positions
      libcpp_vector[double] signal
      PenaltyFactorsIntensity penalties 
      Int charge 

