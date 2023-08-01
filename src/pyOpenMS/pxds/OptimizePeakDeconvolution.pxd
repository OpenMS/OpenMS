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
        # wrap-doc:
        #  This class provides the deconvolution of peak regions using non-linear optimization
        #  
        #  Given a vector of peak shapes, this class optimizes all peak shapes parameters using a non-linear optimization.
        #  For the non-linear optimization we use the Levenberg-Marquardt algorithm.
        #  There are a few constraints for the parameters: the positions are equidistant according to the peptide
        #  mass rule, e.g. two consecutive isotopic peaks are 1.003/charge away from each other. Besides the
        #  peaks have all the same left and right width, respectively

        OptimizePeakDeconvolution() except + nogil 
        OptimizePeakDeconvolution(OptimizePeakDeconvolution &) except + nogil 
        PenaltyFactorsIntensity  getPenalties() except + nogil  # wrap-doc:Returns the penalty parameter
        void setPenalties(PenaltyFactorsIntensity & penalties) except + nogil  # wrap-doc:Sets the penalty parameter
        Int getCharge() except + nogil  # wrap-doc:Returns the charge
        void setCharge(Int charge) except + nogil  # wrap-doc:Sets the charge
        bool optimize(libcpp_vector[ PeakShape ] & peaks, OptimizePeakDeconvolution_Data & data) except + nogil  # wrap-doc:Performs a nonlinear optimization of the peaks that belong to the current isotope pattern
        Size getNumberOfPeaks_(Int charge, libcpp_vector[ PeakShape ] & temp_shapes, OptimizePeakDeconvolution_Data & data) except + nogil 

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePeakDeconvolution.h>" namespace "OpenMS::OptimizationFunctions":
    
    cdef cppclass PenaltyFactorsIntensity(OptimizationFunctions_PenaltyFactors):
        # wrap-inherits:
        #  OptimizationFunctions_PenaltyFactors
        # wrap-doc:
        #  Class for the penalty factors used during the optimization
        #  
        #  A great deviation (squared deviation) of a peak shape's position or its left or right width parameter can be penalised
        #  During the optimization negative heights may occur, they are penalised as well

        PenaltyFactorsIntensity() except + nogil 
        PenaltyFactorsIntensity(PenaltyFactorsIntensity) except + nogil 
        double height


cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePeakDeconvolution.h>" namespace "OpenMS::OptimizePeakDeconvolution":

    cdef cppclass OptimizePeakDeconvolution_Data "OpenMS::OptimizePeakDeconvolution::Data":
      OptimizePeakDeconvolution_Data() except + nogil 
      OptimizePeakDeconvolution_Data(OptimizePeakDeconvolution_Data) except + nogil  # no-wrap

      libcpp_vector[PeakShape] peaks
      libcpp_vector[double] positions
      libcpp_vector[double] signal
      PenaltyFactorsIntensity penalties 
      Int charge 

