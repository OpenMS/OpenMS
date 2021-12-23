from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from String cimport *

cdef extern from "<OpenMS/FORMAT/InspectInfile.h>" namespace "OpenMS":

    cdef cppclass InspectInfile "OpenMS::InspectInfile":

        InspectInfile() nogil except + # wrap-doc:Inspect input file adapter
        InspectInfile(InspectInfile &) nogil except +

        bool operator==(InspectInfile & inspect_infile) nogil except +
        void store(const String & filename) nogil except + # wrap-doc:Stores the experiment data in an Inspect input file that can be used as input for Inspect shell execution
        void handlePTMs(const String & modification_line, const String & modifications_filename, bool monoisotopic) nogil except +
            # wrap-doc:
                #   Retrieves the name, mass change, affected residues, type and position for all modifications from a string
                #   -----
                #   :param modification_line
                #   :param modifications_filename
                #   :param monoisotopic: if true, masses are considered to be monoisotopic
                #   :raises:
                #     Exception: FileNotReadable if the modifications_filename could not be read
                #   :raises:
                #     Exception: FileNotFound if modifications_filename could not be found
                #   :raises:
                #     Exception: ParseError if modifications_filename could not be parsed

        String  getSpectra() nogil except + # wrap-doc:Specifies a spectrum file to search
        void setSpectra(const String & spectra) nogil except + # wrap-doc:Specifies a spectrum file to search
        String  getDb() nogil except + # wrap-doc:Specifies the name of a database (.trie file) to search
        void setDb(const String & db) nogil except + # wrap-doc:Specifies the name of a database (.trie file) to search
        String  getEnzyme() nogil except + # wrap-doc:Specifies the name of a enzyme. "Trypsin", "None", and "Chymotrypsin" are the available values
        void setEnzyme(const String & enzyme) nogil except + # wrap-doc:Specifies the name of a enzyme. "Trypsin", "None", and "Chymotrypsin" are the available values
        Int getModificationsPerPeptide() nogil except + # wrap-doc:Number of PTMs permitted in a single peptide
        void setModificationsPerPeptide(Int modifications_per_peptide) nogil except + # wrap-doc:Number of PTMs permitted in a single peptide
        UInt getBlind() nogil except + # wrap-doc:Run inspect in a blind mode
        void setBlind(UInt blind) nogil except + # wrap-doc:Run inspect in a blind mode
        float getMaxPTMsize() nogil except + # wrap-doc:The maximum modification size (in Da) to consider in a blind search
        void setMaxPTMsize(float maxptmsize) nogil except + # wrap-doc:The maximum modification size (in Da) to consider in a blind search
        float getPrecursorMassTolerance() nogil except + # wrap-doc:Specifies the parent mass tolerance, in Daltons
        void setPrecursorMassTolerance(float precursor_mass_tolerance) nogil except + # wrap-doc:Specifies the parent mass tolerance, in Daltons
        float getPeakMassTolerance() nogil except + # wrap-doc:How far b and y peaks can be shifted from their expected masses.
        void setPeakMassTolerance(float peak_mass_tolerance) nogil except + # wrap-doc:How far b and y peaks can be shifted from their expected masses
        UInt getMulticharge() nogil except + # wrap-doc:If set to true, attempt to guess the precursor charge and mass, and consider multiple charge states if feasible
        void setMulticharge(UInt multicharge) nogil except + # wrap-doc:If set to true, attempt to guess the precursor charge and mass, and consider multiple charge states if feasible
        String  getInstrument() nogil except + # wrap-doc:If set to QTOF, uses a QTOF-derived fragmentation model, and does not attempt to correct the parent mass
        void setInstrument(const String & instrument) nogil except + # wrap-doc:If set to QTOF, uses a QTOF-derived fragmentation model, and does not attempt to correct the parent mass
        Int getTagCount() nogil except + # wrap-doc:Number of tags to generate
        void setTagCount(Int TagCount) nogil except + # wrap-doc:Number of tags to generate

        libcpp_map[ String, libcpp_vector[ String ] ]  getModifications() nogil except + # wrap-ignore

