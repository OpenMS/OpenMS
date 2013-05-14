from libcpp cimport bool
from libcpp.map cimport map as libcpp_map
from libcpp.vector cimport vector as libcpp_vector
from String cimport *
from Types cimport *

cdef extern from "<OpenMS/FORMAT/SequestInfile.h>" namespace "OpenMS":
    
    cdef cppclass SequestInfile "OpenMS::SequestInfile":
        SequestInfile() nogil except +
        SequestInfile(SequestInfile) nogil except +
        bool operator==(SequestInfile &sequest_infile) nogil except +
        void store(String &filename) nogil except +
        String getEnzymeInfoAsString() nogil except +
        String  getDatabase() nogil except +
        void setDatabase(String &database) nogil except +
        String  getNeutralLossesForIons() nogil except +
        void setNeutralLossesForIons(String &neutral_losses_for_ions) nogil except +
        String  getIonSeriesWeights() nogil except +
        void setIonSeriesWeights(String &ion_series_weights) nogil except +
        String  getPartialSequence() nogil except +
        void setPartialSequence(String &partial_sequence) nogil except +
        String  getSequenceHeaderFilter() nogil except +
        void setSequenceHeaderFilter(String &sequence_header_filter) nogil except +
        String  getProteinMassFilter() nogil except +
        void setProteinMassFilter(String &protein_mass_filter) nogil except +
        Real getPeakMassTolerance() nogil except +
        void setPeakMassTolerance(Real peak_mass_tolerance) nogil except +
        Real getPrecursorMassTolerance() nogil except +
        void setPrecursorMassTolerance(Real precursor_mass_tolerance) nogil except +
        Real getMatchPeakTolerance() nogil except +
        void setMatchPeakTolerance(Real match_peak_tolerance) nogil except +
        Real getIonCutoffPercentage() nogil except +
        void setIonCutoffPercentage(Real ion_cutoff_percentage) nogil except +
        Size getPeptideMassUnit() nogil except +
        void setPeptideMassUnit(Size peptide_mass_unit) nogil except +
        Size getOutputLines() nogil except +
        void setOutputLines(Size output_lines) nogil except +
        Size getEnzymeNumber() nogil except +
        String getEnzymeName() nogil except +
        Size setEnzyme(String enzyme_name) nogil except +
        Size getMaxAAPerModPerPeptide() nogil except +
        void setMaxAAPerModPerPeptide(Size max_aa_per_mod_per_peptide) nogil except +
        Size getMaxModsPerPeptide() nogil except +
        void setMaxModsPerPeptide(Size max_mods_per_peptide) nogil except +
        Size getNucleotideReadingFrame() nogil except +
        void setNucleotideReadingFrame(Size nucleotide_reading_frame) nogil except +
        Size getMaxInternalCleavageSites() nogil except +
        void setMaxInternalCleavageSites(Size max_internal_cleavage_sites) nogil except +
        Size getMatchPeakCount() nogil except +
        void setMatchPeakCount(Size match_peak_count) nogil except +
        Size getMatchPeakAllowedError() nogil except +
        void setMatchPeakAllowedError(Size match_peak_allowed_error) nogil except +
        bool getShowFragmentIons() nogil except +
        void setShowFragmentIons(bool show_fragments) nogil except +
        bool getPrintDuplicateReferences() nogil except +
        void setPrintDuplicateReferences(bool print_duplicate_references) nogil except +
        bool getRemovePrecursorNearPeaks() nogil except +
        void setRemovePrecursorNearPeaks(bool remove_precursor_near_peaks) nogil except +
        bool getMassTypeParent() nogil except +
        void setMassTypeParent(bool mass_type_parent) nogil except +
        bool getMassTypeFragment() nogil except +
        void setMassTypeFragment(bool mass_type_fragment) nogil except +
        bool getNormalizeXcorr() nogil except +
        void setNormalizeXcorr(bool normalize_xcorr) nogil except +
        bool getResiduesInUpperCase() nogil except +
        void setResiduesInUpperCase(bool residues_in_upper_case) nogil except +
        void addEnzymeInfo(libcpp_vector[ String ] &enzyme_info) nogil except +
        # libcpp_map[ String, libcpp_vector[ String ] ]  getModifications() nogil except +
        void handlePTMs(String &modification_line, String &modifications_filename, bool monoisotopic) nogil except +

