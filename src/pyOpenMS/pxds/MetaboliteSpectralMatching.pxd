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
        #   ProgressLogger
        #   DefaultParamHandler

        MetaboliteSpectralMatching() except + nogil 
        MetaboliteSpectralMatching(MetaboliteSpectralMatching &) except + nogil  # compiler

        void run(MSExperiment & exp, MSExperiment & speclib, MzTab & mz_tab, String & out_spectra) except + nogil 
        
    cdef cppclass SpectralMatch:

        SpectralMatch() except + nogil # TODO(Whole file)
        SpectralMatch(SpectralMatch &) except + nogil 

        double getObservedPrecursorMass() except + nogil 
        void setObservedPrecursorMass(double) except + nogil 

        double getObservedPrecursorRT() except + nogil 
        void setObservedPrecursorRT(double) except + nogil 

        double getFoundPrecursorMass() except + nogil 
        void setFoundPrecursorMass(double) except + nogil 

        Int getFoundPrecursorCharge() except + nogil 
        void setFoundPrecursorCharge(Int) except + nogil 

        double getMatchingScore() except + nogil 
        void setMatchingScore(double) except + nogil 

        Size getObservedSpectrumIndex() except + nogil 
        void setObservedSpectrumIndex(Size) except + nogil 

        Size getMatchingSpectrumIndex() except + nogil 
        void setMatchingSpectrumIndex(Size) except + nogil 

        String getPrimaryIdentifier() except + nogil 
        void setPrimaryIdentifier(String) except + nogil 

        String getSecondaryIdentifier() except + nogil 
        void setSecondaryIdentifier(String) except + nogil 

        String getCommonName() except + nogil 
        void setCommonName(String) except + nogil 

        String getSumFormula() except + nogil 
        void setSumFormula(String) except + nogil 

        String getInchiString() except + nogil 
        void setInchiString(String) except + nogil 

        String getSMILESString() except + nogil 
        void setSMILESString(String) except + nogil 

        String getPrecursorAdduct() except + nogil 
        void setPrecursorAdduct(String) except + nogil 

# COMMENT: wrap static methods
cdef extern from "<OpenMS/ANALYSIS/ID/MetaboliteSpectralMatching.h>" namespace"OpenMS::MetaboliteSpectralMatching":
        
        # static members        
        # double computeHyperScore(double fragment_mass_error,
        #                         bool fragment_mass_tolerance_unit_ppm,
        #                         MSSpectrum exp_spectrum,
        #                         MSSpectrum db_spectrum) except + nogil   # wrap-attach:MetaboliteSpectralMatching
        # 
        # double computeHyperScore(double fragment_mass_error,
        #                         bool fragment_mass_tolerance_unit_ppm,
        #                         MSSpectrum exp_spectrum,
        #                         MSSpectrum db_spectrum,
        #                         double mz_lower_bound) except + nogil  # wrap-attach:MetaboliteSpectralMatching
        # 
        # double computeHyperScore(double fragment_mass_error,
        #                         bool fragment_mass_tolerance_unit_ppm,
        #                         MSSpectrum exp_spectrum,
        #                         MSSpectrum db_spectrum,
        #                         libcpp_vector[PeptideHit_PeakAnnotation]& annotations) except + nogil  
        #    
        double computeHyperScore(double fragment_mass_error,
                                 bool fragment_mass_tolerance_unit_ppm,
                                 MSSpectrum exp_spectrum,
                                 MSSpectrum db_spectrum,
                                 libcpp_vector[PeptideHit_PeakAnnotation]& annotations,
                                 double mz_lower_bound) except + nogil  # wrap-attach:MetaboliteSpectralMatching
        
