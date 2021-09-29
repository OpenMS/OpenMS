from Types cimport *
from String cimport *

# TODO: support multiset
# typedef std::multiset<String> LabelSet;

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMasses.h>" namespace "OpenMS":
    
    cdef cppclass MultiplexDeltaMasses "OpenMS::MultiplexDeltaMasses":
        # wrap-doc:
            #   Data structure for mass shift pattern
            #   -----
            #   Groups of labelled peptides appear with characteristic mass shifts
            #   -----
            #   For example, for an Arg6 labeled SILAC peptide pair we expect to see
            #   mass shifts of 0 and 6 Da. Or as second example, for a 
            #   peptide pair of a dimethyl labelled sample with a single lysine
            #   we will see mass shifts of 56 Da and 64 Da.
            #   28 Da (N-term) + 28 Da (K) and 34 Da (N-term) + 34 Da (K)
            #   for light and heavy partners respectively
            #   -----
            #   The data structure stores the mass shifts and corresponding labels
            #   for a group of matching peptide features

        MultiplexDeltaMasses() nogil except +
        MultiplexDeltaMasses(MultiplexDeltaMasses &) nogil except +

        MultiplexDeltaMasses(libcpp_vector[ MultiplexDeltaMasses_DeltaMass ] & dm) nogil except +
        libcpp_vector[ MultiplexDeltaMasses_DeltaMass ]  getDeltaMasses() nogil except +

        # String labelSetToString(LabelSet ls) nogil except +

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMasses.h>" namespace "OpenMS::MultiplexDeltaMasses":
    
    cdef cppclass MultiplexDeltaMasses_DeltaMass "OpenMS::MultiplexDeltaMasses::DeltaMass":

        MultiplexDeltaMasses_DeltaMass(double dm, String l) nogil except +
        MultiplexDeltaMasses_DeltaMass(MultiplexDeltaMasses_DeltaMass &) nogil except +

        # DeltaMass(double dm, LabelSet ls);

        double delta_mass
        # LabelSet label_set
