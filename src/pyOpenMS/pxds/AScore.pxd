from libcpp.vector cimport vector as libcpp_vector
from ConsensusMap cimport *
from DefaultParamHandler cimport *
from Feature cimport *
from FeatureMap cimport *
from String cimport *

from ProteinIdentification cimport *
from PeptideIdentification cimport *

from PeptideHit cimport *

from MSExperiment cimport *
from MSSpectrum cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/AScore.h>" namespace "OpenMS":

    cdef cppclass AScore(DefaultParamHandler):
        # wrap-inherits:
        #  DefaultParamHandler
        AScore() except + nogil 
        AScore(AScore &) except + nogil  # compiler

        PeptideHit compute(PeptideHit & hit,
                           MSSpectrum & real_spectrum) except + nogil 


    cdef cppclass ProbablePhosphoSites:

        ProbablePhosphoSites() except + nogil  # compiler
        ProbablePhosphoSites(ProbablePhosphoSites &) except + nogil  # compiler

        Size first
        Size second
        Size seq_1
        Size seq_2
        Size peak_depth
        Size AScore

