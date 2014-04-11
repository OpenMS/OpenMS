from Types cimport *
from IntegerMassDecomposer cimport *

# Cannot use unsigned int as template argument ?
# ctypedef libcpp_map[unsigned int, libcpp_pair[unsigned int, unsigned int] ] constraints_type
ctypedef unsigned long long number_of_decompositions_type

cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/RealMassDecomposer.h>" namespace "OpenMS::ims":
    
    cdef cppclass RealMassDecomposer "OpenMS::ims::RealMassDecomposer":
        RealMassDecomposer() nogil except + #wrap-ignore
        RealMassDecomposer(RealMassDecomposer) nogil except + #wrap-ignore
        RealMassDecomposer(IMSWeights & weights) nogil except +
        # libcpp_vector[ libcpp_vector[unsigned int] ] getDecompositions(double mass, double error) nogil except +
        # libcpp_vector[int] getDecompositions(double mass, double error, constraints_type & constraints) nogil except +
        unsigned long long getNumberOfDecompositions(double mass, double error) nogil except +

