from ProgressLogger cimport *
from DefaultParamHandler cimport *
from MSSpectrum cimport *
from Peak1D cimport *
from MSExperiment cimport *
from MzTab cimport *
from PeptideHit cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/ANALYSIS/ID/MetaboliteSpectralMatching.h>" namespace "OpenMS":

    cdef cppclass MetaboliteSpectralMatching(ProgressLogger, DefaultParamHandler):
        # wrap-inherits:
        #    ProgressLogger
        #    DefaultParamHandler

        MetaboliteSpectralMatching() nogil except +
        MetaboliteSpectralMatching(MetaboliteSpectralMatching) nogil except + 

        void run(MSExperiment & exp, MSExperiment & speclib, MzTab & mz_tab) nogil except +
        
    cdef cppclass SpectralMatch:

        SpectralMatch() nogil except +
        SpectralMatch(SpectralMatch) nogil except + 

        double getObservedPrecursorMass() nogil except +
        void setObservedPrecursorMass(double)nogil except +

        double getObservedPrecursorRT() nogil except +
        void setObservedPrecursorRT(double)nogil except +

        double getFoundPrecursorMass() nogil except +
        void setFoundPrecursorMass(double)nogil except +

        Int getFoundPrecursorCharge() nogil except +
        void setFoundPrecursorCharge(Int)nogil except +

        double getMatchingScore() nogil except +
        void setMatchingScore(double)nogil except +

        Size getObservedSpectrumIndex() nogil except +
        void setObservedSpectrumIndex(Size)nogil except +

        Size getMatchingSpectrumIndex() nogil except +
        void setMatchingSpectrumIndex(Size)nogil except +

        String getPrimaryIdentifier() nogil except +
        void setPrimaryIdentifier(String)nogil except +

        String getSecondaryIdentifier() nogil except +
        void setSecondaryIdentifier(String)nogil except +

        String getCommonName() nogil except +
        void setCommonName(String)nogil except +

        String getSumFormula() nogil except +
        void setSumFormula(String)nogil except +

        String getInchiString() nogil except +
        void setInchiString(String)nogil except +

        String getSMILESString() nogil except +
        void setSMILESString(String)nogil except +

        String getPrecursorAdduct() nogil except +
        void setPrecursorAdduct(String)nogil except +

# COMMENT: wrap static methods
cdef extern from "<OpenMS/ANALYSIS/ID/MetaboliteSpectralMatching.h>" namespace"OpenMS::MetaboliteSpectralMatching":
        
        # static members        
        # double computeHyperScore(double fragment_mass_error,
        #                          bool fragment_mass_tolerance_unit_ppm,
        #                          MSSpectrum exp_spectrum,
        #                          MSSpectrum db_spectrum) nogil except +  # wrap-attach:MetaboliteSpectralMatching
        # 
        # double computeHyperScore(double fragment_mass_error,
        #                          bool fragment_mass_tolerance_unit_ppm,
        #                          MSSpectrum exp_spectrum,
        #                          MSSpectrum db_spectrum,
        #                          double mz_lower_bound) nogil except + # wrap-attach:MetaboliteSpectralMatching
        # 
        # double computeHyperScore(double fragment_mass_error,
        #                          bool fragment_mass_tolerance_unit_ppm,
        #                          MSSpectrum exp_spectrum,
        #                          MSSpectrum db_spectrum,
        #                          libcpp_vector[PeptideHit_PeakAnnotation]& annotations) nogil except + 
        #     
        double computeHyperScore(double fragment_mass_error,
                                 bool fragment_mass_tolerance_unit_ppm,
                                 MSSpectrum exp_spectrum,
                                 MSSpectrum db_spectrum,
                                 libcpp_vector[PeptideHit_PeakAnnotation]& annotations,
                                 double mz_lower_bound) nogil except + # wrap-attach:MetaboliteSpectralMatching
        
