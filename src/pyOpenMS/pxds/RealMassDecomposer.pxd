from Types cimport *
from IntegerMassDecomposer cimport *

# Cannot use unsigned int as template argument ?
# ctypedef libcpp_map[unsigned int, libcpp_pair[unsigned int, unsigned int] ] constraints_type
ctypedef UInt64 number_of_decompositions_type

cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/RealMassDecomposer.h>" namespace "OpenMS::ims":
    
    cdef cppclass RealMassDecomposer "OpenMS::ims::RealMassDecomposer":
        RealMassDecomposer() except + nogil  #wrap-ignore
        RealMassDecomposer(RealMassDecomposer) except + nogil  # compiler
        RealMassDecomposer(IMSWeights & weights) except + nogil 
        # libcpp_vector[ libcpp_vector[unsigned int] ] getDecompositions(double mass, double error) except + nogil 
        # libcpp_vector[int] getDecompositions(double mass, double error, constraints_type & constraints) except + nogil 
        UInt64 getNumberOfDecompositions(double mass, double error) except + nogil 
            # wrap-doc:
                #  Gets a number of all decompositions for amass with an error
                #  allowed. It's similar to thegetDecompositions(double,double) function
                #  but less space consuming, since doesn't use container to store decompositions
                #  
                #  
                #  :param mass: Mass to be decomposed
                #  :param error: Error allowed between given and result decomposition
                #  :return: Number of all decompositions for a given mass and error
