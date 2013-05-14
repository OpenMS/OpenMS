from MSSpectrum cimport *
from FeatureMap cimport *
from ConsensusMap cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from Param cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *
from ProteinIdentification cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/PeptideAndProteinQuant.h>" namespace "OpenMS":

    cdef cppclass PeptideAndProteinQuant(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        PeptideAndProteinQuant()                  nogil except +
        PeptideAndProteinQuant(PeptideAndProteinQuant)   nogil except + #wrap-ignore

        void quantifyPeptides(FeatureMap[Feature] & map_in) nogil except +
        void quantifyPeptides(ConsensusMap & map_in) nogil except +
        void quantifyProteins(ProteinIdentification & proteins) nogil except +

