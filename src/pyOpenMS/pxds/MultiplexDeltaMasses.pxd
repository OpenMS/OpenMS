from Types cimport *
from AASequence cimport *
from String cimport *
from DefaultParamHandler cimport *
from MultiplexDeltaMasses cimport *

# typedef std::multiset<String> LabelSet;

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMasses.h>" namespace "OpenMS":
    
    cdef cppclass MultiplexDeltaMasses "OpenMS::MultiplexDeltaMasses":
        MultiplexDeltaMasses() nogil except +
        MultiplexDeltaMasses(MultiplexDeltaMasses) nogil except + #wrap-ignore
        MultiplexDeltaMasses(libcpp_vector[ DeltaMass ] & dm) nogil except +
        libcpp_vector[ DeltaMass ]  getDeltaMasses() nogil except +
        libcpp_vector[ DeltaMass ]  getDeltaMasses() nogil except +
        # String labelSetToString(LabelSet ls) nogil except +

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMassesGenerator.h>" namespace "OpenMS::MultiplexDeltaMassesGenerator":
    
    cdef cppclass MultiplexDeltaMassesGenerator_Label "OpenMS::MultiplexDeltaMassesGenerator::Label":
        MultiplexDeltaMassesGenerator_Label(MultiplexDeltaMassesGenerator_Label) nogil except + #wrap-ignore
        MultiplexDeltaMassesGenerator_Label(String sn, String ln, String d, double dm) nogil except +

        String short_name
        String long_name
        String description
        double delta_mass

