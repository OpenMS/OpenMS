from Residue cimport *
from ResidueDB cimport *
from DefaultParamHandler cimport *
from FASTAFile cimport *
from libcpp.vector cimport vector as libcpp_vector
from DeconvolvedSpectrum cimport *
from FLASHDeconvHelperStructs cimport Tag
from ProteinHit cimport ProteinHit

cdef extern from "<OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>" namespace "OpenMS":

    cdef cppclass FLASHTaggerAlgorithm(DefaultParamHandler):
        # wrap-inherits:
        #   DefaultParamHandler

        # default constructor
        FLASHTaggerAlgorithm() except + nogil
        # copy constructor
        FLASHTaggerAlgorithm(FLASHTaggerAlgorithm &) except + nogil
        
        void run(DeconvolvedSpectrum & dspec, double ppm) except + nogil
        void run(libcpp_vector[DeconvolvedSpectrum] & dspecs, double ppm) except + nogil
        void runMatching(String fasta) except + nogil

        libcpp_vector[Tag] getTags() except + nogil
        libcpp_vector[Tag] getTags(ProteinHit & hit) except + nogil
        libcpp_vector[ProteinHit] getProteinHits() except + nogil
        libcpp_vector[ProteinHit] getProteinHits(Tag & tag) except + nogil
        int getProteinIndex(ProteinHit & hit) except + nogil
        int getTagIndex(Tag & tag) except + nogil

