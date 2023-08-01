from libcpp.map cimport map as libcpp_map
from libcpp.pair cimport pair as libcpp_pair
from libcpp cimport bool
from Types cimport *
from ProgressLogger cimport *
from DefaultParamHandler cimport *
from FASTAFile cimport *
from MetaInfoInterface cimport *

cdef extern from "<OpenMS/SIMULATION/SimTypes.h>" namespace "OpenMS::SimTypes":

    cdef cppclass SimRandomNumberGenerator:

        SimRandomNumberGenerator() except + nogil  # compiler
        SimRandomNumberGenerator(SimRandomNumberGenerator &) except + nogil  # compiler
        void initialize(bool biological_random, bool technical_random) except + nogil 

    cdef cppclass SimProtein:
        SimProtein(FASTAEntry entry, MetaInfoInterface meta)
        # TODO does this work? 
        # FASTAFile::FASTAEntry entry
        # MetaInfoInterface meta

    ctypedef libcpp_vector[SimProtein] SampleProteins
    ctypedef libcpp_vector[libcpp_vector[SimProtein]] SampleChannels
