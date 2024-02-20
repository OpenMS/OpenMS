from Residue cimport *
from ResidueDB cimport *
from DefaultParamHandler cimport *
from FASTAFile cimport *
from libcpp.vector cimport vector as libcpp_vector
from DeconvolvedSpectrum cimport *

cdef extern from "<OpenMS/ANALYSIS/TOPDOWN/TopDownTagger.h>" namespace "OpenMS":

    cdef cppclass TopDownTagger(DefaultParamHandler):
        # wrap-inherits:

        # default constructor
        TopDownTagger() except + nogil
        # copy constructor
        TopDownTagger(TopDownTagger &) except + nogil  
        
        void run(DeconvolvedSpectrum & dspec, double ppm , libcpp_vector[Tag_FDHS] & Tag) except + nogil