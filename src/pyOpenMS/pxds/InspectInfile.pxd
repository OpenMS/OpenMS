from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from String cimport *

cdef extern from "<OpenMS/FORMAT/InspectInfile.h>" namespace "OpenMS":

    cdef cppclass InspectInfile "OpenMS::InspectInfile":

        InspectInfile() nogil except +
        InspectInfile(InspectInfile) nogil except +

        bool operator==(InspectInfile & inspect_infile) nogil except +
        void store(const String & filename) nogil except +
        void handlePTMs(const String & modification_line, const String & modifications_filename, bool monoisotopic) nogil except +
        String  getSpectra() nogil except +
        void setSpectra(const String & spectra) nogil except +
        String  getDb() nogil except +
        void setDb(const String & db) nogil except +
        String  getEnzyme() nogil except +
        void setEnzyme(const String & enzyme) nogil except +
        Int getModificationsPerPeptide() nogil except +
        void setModificationsPerPeptide(Int modifications_per_peptide) nogil except +
        UInt getBlind() nogil except +
        void setBlind(UInt blind) nogil except +
        float getMaxPTMsize() nogil except +
        void setMaxPTMsize(float maxptmsize) nogil except +
        float getPrecursorMassTolerance() nogil except +
        void setPrecursorMassTolerance(float precursor_mass_tolerance) nogil except +
        float getPeakMassTolerance() nogil except +
        void setPeakMassTolerance(float peak_mass_tolerance) nogil except +
        UInt getMulticharge() nogil except +
        void setMulticharge(UInt multicharge) nogil except +
        String  getInstrument() nogil except +
        void setInstrument(const String & instrument) nogil except +
        Int getTagCount() nogil except +
        void setTagCount(Int TagCount) nogil except +

        libcpp_map[ String, libcpp_vector[ String ] ]  getModifications() nogil except + # wrap-ignore

