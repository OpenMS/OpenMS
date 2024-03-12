from Residue cimport *
from ResidueDB cimport *
from DefaultParamHandler cimport *
from FASTAFile cimport *
from libcpp.vector cimport vector as libcpp_vector
from DeconvolvedSpectrum cimport *

cdef extern from "<OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>" namespace "OpenMS":

    cdef cppclass FLASHTaggerAlgorithm(DefaultParamHandler):
        # wrap-inherits:

        # default constructor
        FLASHTaggerAlgorithm() except + nogil
        # copy constructor
        FLASHTaggerAlgorithm(FLASHTagger &) except + nogil
        
        void run(DeconvolvedSpectrum & dspec, double ppm) except + nogil
        void run(libcpp_vector[DeconvolvedSpectrum] & dspecs, double ppm) except + nogil
