from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from MSSpectrum cimport *
from AASequence cimport *
from ProteinProteinCrossLink cimport *

cdef extern from "<OpenMS/CHEMISTRY/SimpleTSGXLMS.h>" namespace "OpenMS":

    cdef cppclass SimpleTSGXLMS(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler

        SimpleTSGXLMS() nogil except +
        SimpleTSGXLMS(SimpleTSGXLMS) nogil except +

        void getLinearIonSpectrum(libcpp_vector[ SimplePeak ]& spectrum, AASequence peptide,
                Size link_pos, int charge, Size link_pos_2) nogil except +

        void getXLinkIonSpectrum(libcpp_vector[ SimplePeak ]& spectrum, AASequence peptide,
                Size link_pos, double precursor_mass,
                int mincharge, int maxcharge, Size link_pos_2) nogil except +

        void getXLinkIonSpectrum(libcpp_vector[ SimplePeak ]& spectrum, ProteinProteinCrossLink crosslink,
                bool frag_alpha, int mincharge, int maxcharge) nogil except +


cdef extern from "<OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>" namespace "OpenMS::SimpleTSGXLMS":

    cdef cppclass SimplePeak "OpenMS::SimpleTSGXLMS::SimplePeak":

        SimplePeak(SimplePeak) nogil except +
        SimplePeak() nogil except +
        SimplePeak(double mz, int charge) nogil except +

        double mz
        int charge
