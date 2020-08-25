from Types cimport *
from libcpp cimport bool
from DefaultParamHandler cimport *
from MSSpectrum cimport *
from AASequence cimport *
from ProteinProteinCrossLink cimport *

cdef extern from "<OpenMS/CHEMISTRY/TheoreticalSpectrumGeneratorXLMS.h>" namespace "OpenMS":

    cdef cppclass TheoreticalSpectrumGeneratorXLMS(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler

        TheoreticalSpectrumGeneratorXLMS() nogil except +
        TheoreticalSpectrumGeneratorXLMS(TheoreticalSpectrumGeneratorXLMS) nogil except +

        void getLinearIonSpectrum(MSSpectrum& spectrum, AASequence peptide,
                Size link_pos, bool frag_alpha, int charge, Size link_pos_2) nogil except +

        void getXLinkIonSpectrum(MSSpectrum& spectrum, AASequence peptide,
                Size link_pos, double precursor_mass, bool frag_alpha,
                int mincharge, int maxcharge, Size link_pos_2) nogil except +

        void getXLinkIonSpectrum(MSSpectrum& spectrum, ProteinProteinCrossLink crosslink,
                bool frag_alpha, int mincharge, int maxcharge) nogil except +
