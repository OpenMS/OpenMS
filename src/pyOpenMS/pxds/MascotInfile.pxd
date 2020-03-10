from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *
from String cimport *
from File cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/FORMAT/MascotInfile.h>" namespace "OpenMS":
    
    cdef cppclass MascotInfile(ProgressLogger) :
        # wrap-inherits:
        #  ProgressLogger
        MascotInfile() nogil except +
        MascotInfile(MascotInfile) nogil except + #wrap-ignore
        void store(const String & filename, MSSpectrum & spec, double mz, double retention_time, String search_title) nogil except +
        void store(const String & filename, MSExperiment & experiment, String search_title) nogil except +
        void load(const String & filename, MSExperiment & exp) nogil except +
        String  getBoundary() nogil except +
        void setBoundary(const String & boundary) nogil except +
        String  getDB() nogil except +
        void setDB(const String & db) nogil except +
        String  getSearchType() nogil except +
        void setSearchType(const String & search_type) nogil except +
        String  getHits() nogil except +
        void setHits(const String & hits) nogil except +
        String  getCleavage() nogil except +
        void setCleavage(const String & cleavage) nogil except +
        String  getMassType() nogil except +
        void setMassType(const String & mass_type) nogil except +
        libcpp_vector[ String ]  getModifications() nogil except +
        void setModifications(libcpp_vector[ String ] & mods) nogil except +
        libcpp_vector[ String ]  getVariableModifications() nogil except +
        void setVariableModifications(libcpp_vector[ String ] & mods) nogil except +
        String  getInstrument() nogil except +
        void setInstrument(const String & instrument) nogil except +
        UInt getMissedCleavages() nogil except +
        void setMissedCleavages(UInt missed_cleavages) nogil except +
        float getPrecursorMassTolerance() nogil except +
        void setPrecursorMassTolerance(float precursor_mass_tolerance) nogil except +
        float getPeakMassTolerance() nogil except +
        void setPeakMassTolerance(float ion_mass_tolerance) nogil except +
        String  getTaxonomy() nogil except +
        void setTaxonomy(const String & taxonomy) nogil except +
        String  getFormVersion() nogil except +
        void setFormVersion(const String & form_version) nogil except +
        String  getCharges() nogil except +
        void setCharges(libcpp_vector[ int ] & charges) nogil except +

