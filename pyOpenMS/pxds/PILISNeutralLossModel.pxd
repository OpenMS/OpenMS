from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from Map cimport *
from String cimport *
from Types cimport *
from HiddenMarkovModel cimport *
from DefaultParamHandler cimport *
from MSSpectrum cimport *
from RichPeak1D cimport *
from AASequence cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/PILISNeutralLossModel.h>" namespace "OpenMS":
    
    cdef cppclass PILISNeutralLossModel(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        PILISNeutralLossModel() nogil except +
        PILISNeutralLossModel(PILISNeutralLossModel) nogil except +
        DoubleReal train(MSSpectrum[RichPeak1D] & spec, AASequence & peptide, DoubleReal ion_weight, UInt charge, DoubleReal peptide_weight) nogil except +
        void getIons(libcpp_vector[ RichPeak1D ] & peaks, AASequence & peptide, DoubleReal initial_prob) nogil except +
        void setHMM(HiddenMarkovModel & model) nogil except +
        HiddenMarkovModel  getHMM() nogil except +
        void generateModel() nogil except +
        void evaluate() nogil except +

