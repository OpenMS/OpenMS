from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from String cimport *

cdef extern from "<OpenMS/FORMAT/InspectInfile.h>" namespace "OpenMS":

    cdef cppclass InspectInfile "OpenMS::InspectInfile":

        InspectInfile() except + nogil  # wrap-doc:Inspect input file adapter
        InspectInfile(InspectInfile &) except + nogil 

        bool operator==(InspectInfile & inspect_infile) except + nogil 
        void store(const String & filename) except + nogil  # wrap-doc:Stores the experiment data in an Inspect input file that can be used as input for Inspect shell execution
        void handlePTMs(const String & modification_line, const String & modifications_filename, bool monoisotopic) except + nogil 
            # wrap-doc:
                #  Retrieves the name, mass change, affected residues, type and position for all modifications from a string
                #  
                #  
                #  :param modification_line:
                #  :param modifications_filename:
                #  :param monoisotopic: if true, masses are considered to be monoisotopic
                #  :raises:
                #    Exception: FileNotReadable if the modifications_filename could not be read
                #  :raises:
                #    Exception: FileNotFound if modifications_filename could not be found
                #  :raises:
                #    Exception: ParseError if modifications_filename could not be parsed

        String  getSpectra() except + nogil  # wrap-doc:Specifies a spectrum file to search
        void setSpectra(const String & spectra) except + nogil  # wrap-doc:Specifies a spectrum file to search
        String  getDb() except + nogil  # wrap-doc:Specifies the name of a database (.trie file) to search
        void setDb(const String & db) except + nogil  # wrap-doc:Specifies the name of a database (.trie file) to search
        String  getEnzyme() except + nogil  # wrap-doc:Specifies the name of a enzyme. "Trypsin", "None", and "Chymotrypsin" are the available values
        void setEnzyme(const String & enzyme) except + nogil  # wrap-doc:Specifies the name of a enzyme. "Trypsin", "None", and "Chymotrypsin" are the available values
        Int getModificationsPerPeptide() except + nogil  # wrap-doc:Number of PTMs permitted in a single peptide
        void setModificationsPerPeptide(Int modifications_per_peptide) except + nogil  # wrap-doc:Number of PTMs permitted in a single peptide
        UInt getBlind() except + nogil  # wrap-doc:Run inspect in a blind mode
        void setBlind(UInt blind) except + nogil  # wrap-doc:Run inspect in a blind mode
        float getMaxPTMsize() except + nogil  # wrap-doc:The maximum modification size (in Da) to consider in a blind search
        void setMaxPTMsize(float maxptmsize) except + nogil  # wrap-doc:The maximum modification size (in Da) to consider in a blind search
        float getPrecursorMassTolerance() except + nogil  # wrap-doc:Specifies the parent mass tolerance, in Daltons
        void setPrecursorMassTolerance(float precursor_mass_tolerance) except + nogil  # wrap-doc:Specifies the parent mass tolerance, in Daltons
        float getPeakMassTolerance() except + nogil  # wrap-doc:How far b and y peaks can be shifted from their expected masses.
        void setPeakMassTolerance(float peak_mass_tolerance) except + nogil  # wrap-doc:How far b and y peaks can be shifted from their expected masses
        UInt getMulticharge() except + nogil  # wrap-doc:If set to true, attempt to guess the precursor charge and mass, and consider multiple charge states if feasible
        void setMulticharge(UInt multicharge) except + nogil  # wrap-doc:If set to true, attempt to guess the precursor charge and mass, and consider multiple charge states if feasible
        String  getInstrument() except + nogil  # wrap-doc:If set to QTOF, uses a QTOF-derived fragmentation model, and does not attempt to correct the parent mass
        void setInstrument(const String & instrument) except + nogil  # wrap-doc:If set to QTOF, uses a QTOF-derived fragmentation model, and does not attempt to correct the parent mass
        Int getTagCount() except + nogil  # wrap-doc:Number of tags to generate
        void setTagCount(Int TagCount) except + nogil  # wrap-doc:Number of tags to generate

        libcpp_map[ String, libcpp_vector[ String ] ]  getModifications() except + nogil  # wrap-ignore

