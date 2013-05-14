from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from String cimport *
from Map cimport *

cdef extern from "<OpenMS/FORMAT/InspectInfile.h>" namespace "OpenMS":
    
    cdef cppclass InspectInfile "OpenMS::InspectInfile":
        InspectInfile() nogil except +
        InspectInfile(InspectInfile) nogil except +
        bool operator==(InspectInfile & inspect_infile) nogil except +
        void store(String & filename) nogil except +
        void handlePTMs(String & modification_line, String & modifications_filename, bool monoisotopic) nogil except +
        String  getSpectra() nogil except +
        void setSpectra(String & spectra) nogil except +
        String  getDb() nogil except +
        void setDb(String & db) nogil except +
        String  getEnzyme() nogil except +
        void setEnzyme(String & enzyme) nogil except +
        Int getModificationsPerPeptide() nogil except +
        void setModificationsPerPeptide(Int modifications_per_peptide) nogil except +
        UInt getBlind() nogil except +
        void setBlind(UInt blind) nogil except +
        Real getMaxPTMsize() nogil except +
        void setMaxPTMsize(Real maxptmsize) nogil except +
        Real getPrecursorMassTolerance() nogil except +
        void setPrecursorMassTolerance(Real precursor_mass_tolerance) nogil except +
        Real getPeakMassTolerance() nogil except +
        void setPeakMassTolerance(Real peak_mass_tolerance) nogil except +
        UInt getMulticharge() nogil except +
        void setMulticharge(UInt multicharge) nogil except +
        String  getInstrument() nogil except +
        void setInstrument(String & instrument) nogil except +
        Int getTagCount() nogil except +
        void setTagCount(Int TagCount) nogil except +
        # TODO nested map/STL
        # Map[ String, libcpp_vector[ String ] ]  getModifications() nogil except +

