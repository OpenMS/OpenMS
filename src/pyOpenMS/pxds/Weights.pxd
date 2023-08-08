from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from IMSIsotopeDistribution cimport *

cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/Weights.h>" namespace "OpenMS::ims::Weights":

    ctypedef double alphabet_mass_type
    ctypedef long unsigned int weight_type
    ctypedef libcpp_vector[weight_type] weights_type

cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/Weights.h>" namespace "OpenMS::ims":
    
    cdef cppclass IMSWeights "OpenMS::ims::Weights":
        IMSWeights() except + nogil 
        IMSWeights(IMSWeights) except + nogil 

        # IMSWeights(alphabet_masses_type & masses, alphabet_mass_type precision) except + nogil 
        size_type size() except + nogil  # wrap-doc:Gets size of a set of weights
        weight_type getWeight(size_type i) except + nogil  # wrap-doc:Gets a scaled integer weight by index
        void setPrecision(alphabet_mass_type precision) except + nogil  # wrap-doc:Sets a new precision to scale double values to integer
        alphabet_mass_type getPrecision() except + nogil  # wrap-doc:Gets precision.
        # weight_type operator[](size_type i) except + nogil 
        weight_type back() except + nogil  # wrap-doc:Gets a last weight
        alphabet_mass_type getAlphabetMass(size_type i) except + nogil  # wrap-doc:Gets an original (double) alphabet mass by index
        alphabet_mass_type getParentMass(libcpp_vector[ unsigned int ] & decomposition) except + nogil  # wrap-doc:Returns a parent mass for a given `decomposition`
        void swap(size_type index1, size_type index2) except + nogil  # wrap-doc:Exchanges weight and mass at index1 with weight and mass at index2
        bool divideByGCD() except + nogil  # wrap-doc:Divides the integer weights by their gcd. The precision is also adjusted
        alphabet_mass_type getMinRoundingError() except + nogil 
        alphabet_mass_type getMaxRoundingError() except + nogil 
