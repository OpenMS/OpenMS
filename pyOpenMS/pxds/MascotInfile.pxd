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
        void store(String & filename, MSSpectrum[Peak1D] & spec, DoubleReal mz, DoubleReal retention_time, String search_title) nogil except +
        void store(String & filename, MSExperiment[Peak1D, ChromatogramPeak] & experiment, String search_title) nogil except +
        void load(String & filename, MSExperiment[Peak1D, ChromatogramPeak] & exp) nogil except +
        String  getBoundary() nogil except +
        void setBoundary(String & boundary) nogil except +
        String  getDB() nogil except +
        void setDB(String & db) nogil except +
        String  getSearchType() nogil except +
        void setSearchType(String & search_type) nogil except +
        String  getHits() nogil except +
        void setHits(String & hits) nogil except +
        String  getCleavage() nogil except +
        void setCleavage(String & cleavage) nogil except +
        String  getMassType() nogil except +
        void setMassType(String & mass_type) nogil except +
        libcpp_vector[ String ]  getModifications() nogil except +
        void setModifications(libcpp_vector[ String ] & mods) nogil except +
        libcpp_vector[ String ]  getVariableModifications() nogil except +
        void setVariableModifications(libcpp_vector[ String ] & mods) nogil except +
        String  getInstrument() nogil except +
        void setInstrument(String & instrument) nogil except +
        UInt getMissedCleavages() nogil except +
        void setMissedCleavages(UInt missed_cleavages) nogil except +
        Real getPrecursorMassTolerance() nogil except +
        void setPrecursorMassTolerance(Real precursor_mass_tolerance) nogil except +
        Real getPeakMassTolerance() nogil except +
        void setPeakMassTolerance(Real ion_mass_tolerance) nogil except +
        String  getTaxonomy() nogil except +
        void setTaxonomy(String & taxonomy) nogil except +
        String  getFormVersion() nogil except +
        void setFormVersion(String & form_version) nogil except +
        String  getCharges() nogil except +
        void setCharges(libcpp_vector[ int ] & charges) nogil except +

