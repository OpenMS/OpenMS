from Types cimport *
from libcpp cimport bool
from String cimport *
from ModificationDefinitionsSet cimport *
# from XMLFile cimport *

cdef extern from "<OpenMS/FORMAT/XTandemInfile.h>" namespace "OpenMS":
    
    cdef cppclass XTandemInfile: # (XMLFile):
        XTandemInfile() nogil except +
        XTandemInfile(XTandemInfile) nogil except + #wrap-ignore

        void setFragmentMassTolerance(double tolerance) nogil except +
        double getFragmentMassTolerance() nogil except +

        void setPrecursorMassTolerancePlus(double tol) nogil except +
        double getPrecursorMassTolerancePlus() nogil except +
        void setPrecursorMassToleranceMinus(double tol) nogil except +
        double getPrecursorMassToleranceMinus() nogil except +

        void setPrecursorErrorType(MassType mono_isotopic) nogil except +
        MassType getPrecursorErrorType() nogil except +

        void setFragmentMassErrorUnit(ErrorUnit unit) nogil except +
        ErrorUnit getFragmentMassErrorUnit() nogil except +
        void setPrecursorMassErrorUnit(ErrorUnit unit) nogil except +
        ErrorUnit getPrecursorMassErrorUnit() nogil except +

        void setNumberOfThreads(UInt threads) nogil except +
        UInt getNumberOfThreads() nogil except +

        void setModifications(ModificationDefinitionsSet & mods) nogil except +
        ModificationDefinitionsSet  getModifications() nogil except +

        void setOutputFilename(const String & output) nogil except +
        String getOutputFilename() nogil except +
        void setInputFilename(const String & input_file) nogil except +
        String getInputFilename() nogil except +
        void setTaxonomyFilename(const String & filename) nogil except +
        String getTaxonomyFilename() nogil except +
        void setDefaultParametersFilename(const String & filename) nogil except +
        String getDefaultParametersFilename() nogil except +

        void setTaxon(const String & taxon) nogil except +
        String getTaxon() nogil except +

        void setMaxPrecursorCharge(Int max_charge) nogil except +
        Int getMaxPrecursorCharge() nogil except +

        void setNumberOfMissedCleavages(UInt missed_cleavages) nogil except +
        UInt getNumberOfMissedCleavages() nogil except +

        void setOutputResults(String result) nogil except +
        String getOutputResults() nogil except +

        void setMaxValidEValue(double value) nogil except +
        double getMaxValidEValue() nogil except +

        void setSemiCleavage(bool semi_cleavage) nogil except +
        void setAllowIsotopeError(bool allow_isotope_error) nogil except +

        # bool getNoiseSuppression() nogil except +
        # void setNoiseSuppression(bool noise_suppression) nogil except +

        void write(String filename, bool ignore_member_parameters, bool force_default_mods) nogil except +

        void setCleavageSite(String cleavage_site) nogil except +
        String getCleavageSite() nogil except +

cdef extern from "<OpenMS/FORMAT/XTandemInfile.h>" namespace "OpenMS::XTandemInfile":
    cdef enum ErrorUnit "OpenMS::XTandemInfile::ErrorUnit":
        #wrap-attach:
        #    XTandemInfile
        DALTONS
        PPM

cdef extern from "<OpenMS/FORMAT/XTandemInfile.h>" namespace "OpenMS::XTandemInfile":
    cdef enum MassType "OpenMS::XTandemInfile::MassType":
        #wrap-attach:
        #    XTandemInfile
        MONOISOTOPIC
        AVERAGE

