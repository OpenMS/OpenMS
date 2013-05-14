from MSSpectrum cimport *
from FeatureMap cimport *
from ConsensusMap cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from Param cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ItraqQuantifier.h>" namespace "OpenMS":

    cdef cppclass ItraqQuantifier(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        ItraqQuantifier()                  nogil except +
        ItraqQuantifier(ItraqQuantifier)   nogil except + #wrap-ignore

        void run(ConsensusMap & map_in, ConsensusMap & map_out) nogil except +

