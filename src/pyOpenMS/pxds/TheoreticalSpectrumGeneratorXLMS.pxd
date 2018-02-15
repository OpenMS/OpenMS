from Types cimport *
from libcpp cimport bool
from Residue cimport *
from DefaultParamHandler cimport *
from MSSpectrum cimport *
from AASequence cimport *

cdef extern from "<OpenMS/CHEMISTRY/TheoreticalSpectrumGeneratorXLMS.h>" namespace "OpenMS":
    
    cdef cppclass TheoreticalSpectrumGeneratorXLMS(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler

        TheoreticalSpectrumGeneratorXLMS() nogil except +
        TheoreticalSpectrumGeneratorXLMS(TheoreticalSpectrumGeneratorXLMS) nogil except +

        void getCommonIonSpectrum(MSSpectrum& spectrum, AASequence peptide,
                Size link_pos, bool frag_alpha, int charge, Size link_pos_2) nogil except +

        void getXLinkIonSpectrum(MSSpectrum& spectrum, AASequence peptide,
                Size link_pos, double precursor_mass, bool frag_alpha,
                int mincharge, int maxcharge, Size link_pos_2) nogil except +

