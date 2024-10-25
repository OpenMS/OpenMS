from Types cimport *
from String cimport *
from DefaultParamHandler cimport *
from MultiplexDeltaMasses cimport *

cdef extern from "<OpenMS/FEATUREFINDER/MultiplexDeltaMassesGenerator.h>" namespace "OpenMS":
    
    cdef cppclass MultiplexDeltaMassesGenerator(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        MultiplexDeltaMassesGenerator() except + nogil 
        MultiplexDeltaMassesGenerator(MultiplexDeltaMassesGenerator &) except + nogil  # compiler
        MultiplexDeltaMassesGenerator(String labels, int missed_cleavages, libcpp_map[ String, double ] label_mass_shift) except + nogil 

        void generateKnockoutDeltaMasses() except + nogil 
        libcpp_vector[ MultiplexDeltaMasses ] getDeltaMassesList() except + nogil 
        # libcpp_vector[ libcpp_vector[ String ] ] getSamplesLabelsList() except + nogil 
        String getLabelShort(String label) except + nogil 
        String getLabelLong(String label) except + nogil 

        # missing multiset support
        # NAMESPACE # MultiplexDeltaMasses::LabelSet extractLabelSet(AASequence sequence) except + nogil 

cdef extern from "<OpenMS/FEATUREFINDER/MultiplexDeltaMassesGenerator.h>" namespace "OpenMS::MultiplexDeltaMassesGenerator":
    
    cdef cppclass MultiplexDeltaMassesGenerator_Label "OpenMS::MultiplexDeltaMassesGenerator::Label":

        MultiplexDeltaMassesGenerator_Label(MultiplexDeltaMassesGenerator_Label) except + nogil  #wrap-ignore
        MultiplexDeltaMassesGenerator_Label(String sn, String ln, String d, double dm) except + nogil 

        String short_name
        String long_name
        String description
        double delta_mass


