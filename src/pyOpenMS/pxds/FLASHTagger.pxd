from Residue cimport *
from ResidueDB cimport *
from DefaultParamHandler cimport *
from FASTAFile cimport *
from libcpp.vector cimport vector as libcpp_vector
from DeconvolvedSpectrum cimport *
from FLASHHelperClasses cimport Tag
from ProteinHit cimport ProteinHit

cdef extern from "<OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>" namespace "OpenMS":

    cdef cppclass FLASHTaggerAlgorithm(DefaultParamHandler):
        # wrap-inherits:
        #   DefaultParamHandler

        # default constructor
        FLASHTaggerAlgorithm() except + nogil
        # copy constructor
        FLASHTaggerAlgorithm(FLASHTaggerAlgorithm &) except + nogil
        
        void run(DeconvolvedSpectrum & deconvolved_spectrum, double ppm) except + nogil


