from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from PeakShape cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePick.h>" namespace "OpenMS":
    
    cdef cppclass OptimizePick "OpenMS::OptimizePick":
        # wrap-doc:
            #  This class provides the non-linear optimization of the peak parameters
            #  
            #  Given a vector of peak shapes, this class optimizes all peak shapes parameters using a non-linear optimization
            #  For the non-linear optimization we use the Levenberg-Marquardt algorithm provided by the Eigen

        OptimizePick() except + nogil 
        OptimizePick(OptimizePick &) except + nogil  # compiler
        OptimizePick(OptimizationFunctions_PenaltyFactors penalties_, int max_iteration_) except + nogil 
        OptimizationFunctions_PenaltyFactors getPenalties() except + nogil  # wrap-doc:Returns the penalty factors
        void setPenalties(OptimizationFunctions_PenaltyFactors penalties) except + nogil  # wrap-doc:Sets the penalty factors
        unsigned int  getNumberIterations() except + nogil  # wrap-doc:Returns the number of iterations
        void setNumberIterations(int max_iteration) except + nogil  # wrap-doc:Sets the number of iterations
        # void optimize(libcpp_vector[ PeakShape ] & peaks, OptimizePick_Data & data) except + nogil 

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePick.h>" namespace "OpenMS::OptimizationFunctions":
    
    cdef cppclass OptimizationFunctions_PenaltyFactors "OpenMS::OptimizationFunctions::PenaltyFactors":
        OptimizationFunctions_PenaltyFactors() except + nogil 
        OptimizationFunctions_PenaltyFactors(OptimizationFunctions_PenaltyFactors) except + nogil  # wrap-ignore
        double pos
        double lWidth
        double rWidth

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePick.h>" namespace "OpenMS::OptimizePick":
    
    cdef cppclass OptimizePick_Data "OpenMS::OptimizePick::Data":
        OptimizePick_Data(OptimizePick_Data) except + nogil  #wrap-ignore
        libcpp_vector[ double ] positions
        libcpp_vector[ double ] signal

        # TODO STL attribute
        # OptimizationFunctions_PenaltyFactors penalties
        # libcpp_vector[ PeakShape ] peaks

