from Types cimport *
from libcpp cimport bool
from String cimport *
from ModificationDefinitionsSet cimport *
# from XMLFile cimport *

cdef extern from "<OpenMS/FORMAT/XTandemInfile.h>" namespace "OpenMS":
    
    cdef cppclass XTandemInfile: # (XMLFile):
        XTandemInfile() except + nogil 
        # protected
        XTandemInfile(XTandemInfile &) except + nogil  # wrap-ignore

        void setFragmentMassTolerance(double tolerance) except + nogil 
        double getFragmentMassTolerance() except + nogil 

        void setPrecursorMassTolerancePlus(double tol) except + nogil 
        double getPrecursorMassTolerancePlus() except + nogil 
        void setPrecursorMassToleranceMinus(double tol) except + nogil 
        double getPrecursorMassToleranceMinus() except + nogil 

        void setPrecursorErrorType(MassType mono_isotopic) except + nogil 
        MassType getPrecursorErrorType() except + nogil 

        void setFragmentMassErrorUnit(ErrorUnit unit) except + nogil 
        ErrorUnit getFragmentMassErrorUnit() except + nogil 
        void setPrecursorMassErrorUnit(ErrorUnit unit) except + nogil 
        ErrorUnit getPrecursorMassErrorUnit() except + nogil 

        void setNumberOfThreads(UInt threads) except + nogil 
        UInt getNumberOfThreads() except + nogil 

        void setModifications(ModificationDefinitionsSet & mods) except + nogil 
        ModificationDefinitionsSet  getModifications() except + nogil 

        void setOutputFilename(const String & output) except + nogil 
        String getOutputFilename() except + nogil 
        void setInputFilename(const String & input_file) except + nogil 
        String getInputFilename() except + nogil 
        void setTaxonomyFilename(const String & filename) except + nogil 
        String getTaxonomyFilename() except + nogil 
        void setDefaultParametersFilename(const String & filename) except + nogil 
        String getDefaultParametersFilename() except + nogil 

        void setTaxon(const String & taxon) except + nogil 
        String getTaxon() except + nogil 

        void setMaxPrecursorCharge(Int max_charge) except + nogil 
        Int getMaxPrecursorCharge() except + nogil 

        void setNumberOfMissedCleavages(UInt missed_cleavages) except + nogil 
        UInt getNumberOfMissedCleavages() except + nogil 

        void setOutputResults(String result) except + nogil 
        String getOutputResults() except + nogil 

        void setMaxValidEValue(double value) except + nogil 
        double getMaxValidEValue() except + nogil 

        void setSemiCleavage(bool semi_cleavage) except + nogil 
        void setAllowIsotopeError(bool allow_isotope_error) except + nogil 

        # bool getNoiseSuppression() except + nogil 
        # void setNoiseSuppression(bool noise_suppression) except + nogil 

        void write(String filename, bool ignore_member_parameters, bool force_default_mods) except + nogil 

        void setCleavageSite(String cleavage_site) except + nogil 
        String getCleavageSite() except + nogil 

cdef extern from "<OpenMS/FORMAT/XTandemInfile.h>" namespace "OpenMS::XTandemInfile":
    cdef enum ErrorUnit "OpenMS::XTandemInfile::ErrorUnit":
        #wrap-attach:
        #   XTandemInfile
        DALTONS
        PPM

cdef extern from "<OpenMS/FORMAT/XTandemInfile.h>" namespace "OpenMS::XTandemInfile":
    cdef enum MassType "OpenMS::XTandemInfile::MassType":
        #wrap-attach:
        #   XTandemInfile
        MONOISOTOPIC
        AVERAGE

