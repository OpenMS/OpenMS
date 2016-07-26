from Types cimport *
from AASequence cimport *
from String cimport *
from DefaultParamHandler cimport *
from MultiplexDeltaMasses cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMassesGenerator.h>" namespace "OpenMS":
    
    cdef cppclass MultiplexDeltaMassesGenerator(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        MultiplexDeltaMassesGenerator() nogil except +
        MultiplexDeltaMassesGenerator(MultiplexDeltaMassesGenerator) nogil except + #wrap-ignore
        MultiplexDeltaMassesGenerator(String labels, int missed_cleavages, libcpp_map[ String, double ] label_mass_shift) nogil except +

        void generateKnockoutDeltaMasses() nogil except +
        void printSamplesLabelsList() nogil except +
        void printDeltaMassesList() nogil except +
        libcpp_vector[ MultiplexDeltaMasses ] getDeltaMassesList() nogil except +
        libcpp_vector[ MultiplexDeltaMasses ]  getDeltaMassesList() nogil except +
        libcpp_vector[ libcpp_vector[ String ] ] getSamplesLabelsList() nogil except +
        libcpp_vector[ libcpp_vector[ String ] ]  getSamplesLabelsList() nogil except +
        String getLabelShort(String label) nogil except +
        String getLabelLong(String label) nogil except +

        # NAMESPACE # MultiplexDeltaMasses::LabelSet extractLabelSet(AASequence sequence) nogil except +

