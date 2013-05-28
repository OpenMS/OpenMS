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
        IMSWeights() nogil except +
        IMSWeights(IMSWeights) nogil except +
        # IMSWeights(alphabet_masses_type & masses, alphabet_mass_type precision) nogil except +
        size_type size() nogil except +
        weight_type getWeight(size_type i) nogil except +
        void setPrecision(alphabet_mass_type precision) nogil except +
        alphabet_mass_type getPrecision() nogil except +
        # weight_type operator[](size_type i) nogil except +
        weight_type back() nogil except +
        alphabet_mass_type getAlphabetMass(size_type i) nogil except +
        alphabet_mass_type getParentMass(libcpp_vector[ unsigned int ] & decomposition) nogil except +
        void swap(size_type index1, size_type index2) nogil except +
        bool divideByGCD() nogil except +
        alphabet_mass_type getMinRoundingError() nogil except +
        alphabet_mass_type getMaxRoundingError() nogil except +

