from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from PeakShape cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePick.h>" namespace "OpenMS":
    
    cdef cppclass OptimizePick "OpenMS::OptimizePick":
        OptimizePick() nogil except +
        OptimizePick(OptimizePick) nogil except + #wrap-ignore
        OptimizePick(OptimizationFunctions_PenaltyFactors penalties_, int max_iteration_) nogil except +
        OptimizationFunctions_PenaltyFactors getPenalties() nogil except +
        void setPenalties(OptimizationFunctions_PenaltyFactors penalties) nogil except +
        unsigned int  getNumberIterations() nogil except +
        void setNumberIterations(int max_iteration) nogil except +
        # void optimize(libcpp_vector[ PeakShape ] & peaks, OptimizePick_Data & data) nogil except +

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePick.h>" namespace "OpenMS::OptimizationFunctions":
    
    cdef cppclass OptimizationFunctions_PenaltyFactors "OpenMS::OptimizationFunctions::PenaltyFactors":
        OptimizationFunctions_PenaltyFactors() nogil except +
        OptimizationFunctions_PenaltyFactors(OptimizationFunctions_PenaltyFactors) nogil except +
        double pos
        double lWidth
        double rWidth

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePick.h>" namespace "OpenMS::OptimizePick":
    
    cdef cppclass OptimizePick_Data "OpenMS::OptimizePick::Data":
        OptimizePick_Data(OptimizePick_Data) nogil except + #wrap-ignore
        libcpp_vector[ double ] positions
        libcpp_vector[ double ] signal

        # TODO STL attribute
        # OptimizationFunctions_PenaltyFactors penalties
        # libcpp_vector[ PeakShape ] peaks

