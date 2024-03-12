from Residue cimport *
from ResidueDB cimport *
from DefaultParamHandler cimport *
from FASTAFile cimport *
from libcpp.vector cimport vector as libcpp_vector
from DeconvolvedSpectrum cimport *

cdef extern from "<OpenMS/ANALYSIS/TOPDOWN/FLASHTagger.h>" namespace "OpenMS":

    cdef cppclass FLASHTagger(DefaultParamHandler):
        # wrap-inherits:

        # default constructor
        FLASHTagger() except + nogil
        # copy constructor
        FLASHTagger(FLASHTagger &) except + nogil
        
        void run(DeconvolvedSpectrum & dspec, double ppm) except + nogil
        void run(std::vector<DeconvolvedSpectrum> & dspecs, double ppm) except + nogil