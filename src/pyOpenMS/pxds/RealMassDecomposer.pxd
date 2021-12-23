from Types cimport *
from IntegerMassDecomposer cimport *

# Cannot use unsigned int as template argument ?
# ctypedef libcpp_map[unsigned int, libcpp_pair[unsigned int, unsigned int] ] constraints_type
ctypedef UInt64 number_of_decompositions_type

cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/RealMassDecomposer.h>" namespace "OpenMS::ims":
    
    cdef cppclass RealMassDecomposer "OpenMS::ims::RealMassDecomposer":
        RealMassDecomposer() nogil except + #wrap-ignore
        RealMassDecomposer(RealMassDecomposer) nogil except + # compiler
        RealMassDecomposer(IMSWeights & weights) nogil except +
        # libcpp_vector[ libcpp_vector[unsigned int] ] getDecompositions(double mass, double error) nogil except +
        # libcpp_vector[int] getDecompositions(double mass, double error, constraints_type & constraints) nogil except +
        UInt64 getNumberOfDecompositions(double mass, double error) nogil except +
            # wrap-doc:
                #   Gets a number of all decompositions for amass with an error
                #   allowed. It's similar to thegetDecompositions(double,double) function
                #   but less space consuming, since doesn't use container to store decompositions
                #   -----
                #   :param mass: Mass to be decomposed
                #   :param error: Error allowed between given and result decomposition
                #   :returns: Number of all decompositions for a given mass and error
