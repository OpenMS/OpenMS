from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool

from ConsensusMap cimport *
from DefaultParamHandler cimport *
from Feature cimport *
from FeatureMap cimport *
from String cimport *

from ProteinIdentification cimport *
from PeptideIdentification cimport *

from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *

cdef extern from "<OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h>" namespace "OpenMS::Math":

    cdef cppclass PosteriorErrorProbabilityModel(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        PosteriorErrorProbabilityModel() nogil except +
        PosteriorErrorProbabilityModel(PosteriorErrorProbabilityModel) nogil except +   # wrap-ignore

        bool fit(libcpp_vector[double] & search_engine_scores) nogil except +
        bool fit(libcpp_vector[double] & search_engine_scores, libcpp_vector[double] & probabilities) nogil except +
        DoubleReal computeProbability(DoubleReal score) nogil except +

