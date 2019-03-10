from Types cimport *
from String cimport *

# TODO: support multiset
# typedef std::multiset<String> LabelSet;

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMasses.h>" namespace "OpenMS":
    
    cdef cppclass MultiplexDeltaMasses "OpenMS::MultiplexDeltaMasses":
        MultiplexDeltaMasses() nogil except +
        MultiplexDeltaMasses(MultiplexDeltaMasses) nogil except + #wrap-ignore

        MultiplexDeltaMasses(libcpp_vector[ MultiplexDeltaMasses_DeltaMass ] & dm) nogil except +
        libcpp_vector[ MultiplexDeltaMasses_DeltaMass ]  getDeltaMasses() nogil except +

        # String labelSetToString(LabelSet ls) nogil except +

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMasses.h>" namespace "OpenMS::MultiplexDeltaMasses":
    
    cdef cppclass MultiplexDeltaMasses_DeltaMass "OpenMS::MultiplexDeltaMasses::DeltaMass":

        MultiplexDeltaMasses_DeltaMass(MultiplexDeltaMasses_DeltaMass) nogil except + #wrap-ignore
        MultiplexDeltaMasses_DeltaMass(double dm, String l) nogil except +
        # DeltaMass(double dm, LabelSet ls);

        double delta_mass
        # LabelSet label_set
      
      
