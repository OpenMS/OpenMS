from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from PeakShape cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePick.h>" namespace "OpenMS":
    
    cdef cppclass OptimizePick "OpenMS::OptimizePick":
        OptimizePick() nogil except +
        OptimizePick(OptimizePick) nogil except + #wrap-ignore
        # NAMESPACE #  OptimizePick(struct OptimizationFunctions::PenaltyFactors & penalties_, int max_iteration_, double eps_abs_, double eps_rel_) nogil except +
        # NAMESPACE # struct OptimizationFunctions::PenaltyFactors  getPenalties() nogil except +
        # NAMESPACE # struct OptimizationFunctions::PenaltyFactors  getPenalties() nogil except +
        # NAMESPACE # void setPenalties(struct OptimizationFunctions::PenaltyFactors & penalties) nogil except +
        unsigned int  getNumberIterations() nogil except +
        void setNumberIterations(int max_iteration) nogil except +
        # void optimize(libcpp_vector[ PeakShape ] & peaks, Data & data) nogil except +

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePick.h>" namespace "OpenMS::OptimizationFunctions":
    
    cdef cppclass OptimizationFunctions_PenaltyFactors "OpenMS::OptimizationFunctions::PenaltyFactors":
        OptimizationFunctions_PenaltyFactors() nogil except +
        OptimizationFunctions_PenaltyFactors(OptimizationFunctions_PenaltyFactors) nogil except +
        double pos
        double lWidth
        double rWidth

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePick.h>" namespace "OpenMS::OptimizePick":
    
    cdef cppclass Data "OpenMS::OptimizePick::Data":
        Data(Data) nogil except + #wrap-ignore
        libcpp_vector[ double ] positions
        libcpp_vector[ double ] signal
        # TODO STL attribute
        # libcpp_vector[ PeakShape ] peaks
        # NAMESPACE # OptimizationFunctions::PenaltyFactors penalties

