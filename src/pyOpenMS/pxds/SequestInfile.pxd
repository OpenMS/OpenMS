from libcpp cimport bool
from libcpp.map cimport map as libcpp_map
from libcpp.vector cimport vector as libcpp_vector
from String cimport *
from Types cimport *

cdef extern from "<OpenMS/FORMAT/SequestInfile.h>" namespace "OpenMS":
    
    cdef cppclass SequestInfile "OpenMS::SequestInfile":
        SequestInfile() except + nogil  # wrap-doc:Sequest input file adapter
        SequestInfile(SequestInfile &) except + nogil 
        bool operator==(SequestInfile &sequest_infile) except + nogil 
        void store(const String &filename) except + nogil 
            # wrap-doc:
                #  Stores the experiment data in a Sequest input file that can be used as input for Sequest shell execution
                #  
                #  :param filename: the name of the file in which the infile is stored into

        String getEnzymeInfoAsString() except + nogil  # wrap-doc:Returns the enzyme list as a string
        String  getDatabase() except + nogil  # wrap-doc:Returns the used database
        void setDatabase(const String &database) except + nogil  # wrap-doc:Sets the used database
        String  getNeutralLossesForIons() except + nogil  # wrap-doc:Returns whether neutral losses are considered for the a-, b- and y-ions
        void setNeutralLossesForIons(const String &neutral_losses_for_ions) except + nogil  # wrap-doc:Sets whether neutral losses are considered for the a-, b- and y-ions
        String  getIonSeriesWeights() except + nogil  # wrap-doc:Returns the weights for the a-, b-, c-, d-, v-, w-, x-, y- and z-ion series
        void setIonSeriesWeights(const String &ion_series_weights) except + nogil  # wrap-doc:Sets the weights for the a-, b-, c-, d-, v-, w-, x-, y- and z-ion series
        String  getPartialSequence() except + nogil  # wrap-doc:Returns the partial sequences (space delimited) that have to occur in the theoretical spectra
        void setPartialSequence(const String &partial_sequence) except + nogil  # wrap-doc:Sets the partial sequences (space delimited) that have to occur in the theoretical spectra
        String  getSequenceHeaderFilter() except + nogil  # wrap-doc:Returns the sequences (space delimited) that have to occur, or be absent (preceded by a tilde) in the header of a protein to be considered
        void setSequenceHeaderFilter(const String &sequence_header_filter) except + nogil  # wrap-doc:Sets the sequences (space delimited) that have to occur, or be absent (preceded by a tilde) in the header of a protein to be considered
        String  getProteinMassFilter() except + nogil  # wrap-doc:Returns the protein mass filter (either min and max mass, or mass and tolerance value in percent)
        void setProteinMassFilter(const String &protein_mass_filter) except + nogil  # wrap-doc:Sets the protein mass filter (either min and max mass, or mass and tolerance value in percent)
        float getPeakMassTolerance() except + nogil  # wrap-doc:Returns the peak mass tolerance
        void setPeakMassTolerance(float peak_mass_tolerance) except + nogil  # wrap-doc:Sets the peak mass tolerance
        float getPrecursorMassTolerance() except + nogil  # wrap-doc:Returns the precursor mass tolerance
        void setPrecursorMassTolerance(float precursor_mass_tolerance) except + nogil  # wrap-doc:Sets the precursor mass tolerance
        float getMatchPeakTolerance() except + nogil  # wrap-doc:Returns the match peak tolerance
        void setMatchPeakTolerance(float match_peak_tolerance) except + nogil  # wrap-doc:Sets the match peak tolerance
        float getIonCutoffPercentage() except + nogil  # wrap-doc:Returns the the cutoff of the ratio matching theoretical peaks/theoretical peaks
        void setIonCutoffPercentage(float ion_cutoff_percentage) except + nogil  # wrap-doc:Sets the ion cutoff of the ratio matching theoretical peaks/theoretical peaks
        Size getPeptideMassUnit() except + nogil  # wrap-doc:Returns the peptide mass unit
        void setPeptideMassUnit(Size peptide_mass_unit) except + nogil  # wrap-doc:Sets the peptide mass unit
        Size getOutputLines() except + nogil  # wrap-doc:Returns the number of peptides to be displayed
        void setOutputLines(Size output_lines) except + nogil  # wrap-doc:Sets the number of peptides to be displayed
        Size getEnzymeNumber() except + nogil  # wrap-doc:Returns the enzyme used for cleavage (by means of the number from a list of enzymes)
        String getEnzymeName() except + nogil  # wrap-doc:Returns the enzyme used for cleavage
        Size setEnzyme(String enzyme_name) except + nogil  # wrap-doc:Sets the enzyme used for cleavage (by means of the number from a list of enzymes)
        Size getMaxAAPerModPerPeptide() except + nogil  # wrap-doc:Returns the maximum number of amino acids containing the same modification in a peptide
        void setMaxAAPerModPerPeptide(Size max_aa_per_mod_per_peptide) except + nogil  # wrap-doc:Sets the maximum number of amino acids containing the same modification in a peptide
        Size getMaxModsPerPeptide() except + nogil  # wrap-doc:Returns the maximum number of modifications that are allowed in a peptide
        void setMaxModsPerPeptide(Size max_mods_per_peptide) except + nogil  # wrap-doc:Sets the maximum number of modifications that are allowed in a peptide
        Size getNucleotideReadingFrame() except + nogil  # wrap-doc:Returns the nucleotide reading frame
        void setNucleotideReadingFrame(Size nucleotide_reading_frame) except + nogil  # wrap-doc:Sets the nucleotide reading frame
        Size getMaxInternalCleavageSites() except + nogil  # wrap-doc:Returns the maximum number of internal cleavage sites
        void setMaxInternalCleavageSites(Size max_internal_cleavage_sites) except + nogil  # wrap-doc:Sets the maximum number of internal cleavage sites
        Size getMatchPeakCount() except + nogil  # wrap-doc:Returns the number of top abundant peaks to match with theoretical ones
        void setMatchPeakCount(Size match_peak_count) except + nogil  # wrap-doc:Sets the number of top abundant peaks to with theoretical ones
        Size getMatchPeakAllowedError() except + nogil  # wrap-doc:Returns the number of top abundant peaks that are allowed not to match with a theoretical peak
        void setMatchPeakAllowedError(Size match_peak_allowed_error) except + nogil  # wrap-doc:Sets the number of top abundant peaks that are allowed not to match with a theoretical peak
        bool getShowFragmentIons() except + nogil  # wrap-doc:Returns whether fragment ions shall be displayed
        void setShowFragmentIons(bool show_fragments) except + nogil  # wrap-doc:Sets whether fragment ions shall be displayed
        bool getPrintDuplicateReferences() except + nogil  # wrap-doc:Returns whether all proteins containing a found peptide should be displayed
        void setPrintDuplicateReferences(bool print_duplicate_references) except + nogil  # wrap-doc:Sets whether all proteins containing a found peptide should be displayed
        bool getRemovePrecursorNearPeaks() except + nogil  # wrap-doc:Returns whether peaks near (15 amu) the precursor peak are removed
        void setRemovePrecursorNearPeaks(bool remove_precursor_near_peaks) except + nogil  # wrap-doc:Sets whether peaks near (15 amu) the precursor peak are removed
        bool getMassTypeParent() except + nogil  # wrap-doc:Returns the mass type of the parent (0 - monoisotopic, 1 - average mass)
        void setMassTypeParent(bool mass_type_parent) except + nogil  # wrap-doc:Sets the mass type of the parent (0 - monoisotopic, 1 - average mass)
        bool getMassTypeFragment() except + nogil  # wrap-doc:Returns the mass type of the fragments (0 - monoisotopic, 1 - average mass)
        void setMassTypeFragment(bool mass_type_fragment) except + nogil  # wrap-doc:Sets the mass type of the fragments (0 - monoisotopic, 1 - average mass)
        bool getNormalizeXcorr() except + nogil  # wrap-doc:Returns whether normalized xcorr values are displayed
        void setNormalizeXcorr(bool normalize_xcorr) except + nogil  # wrap-doc:Sets whether normalized xcorr values are displayed
        bool getResiduesInUpperCase() except + nogil  # wrap-doc:Returns whether residues are in upper case
        void setResiduesInUpperCase(bool residues_in_upper_case) except + nogil  # wrap-doc:Sets whether residues are in upper case
        void addEnzymeInfo(libcpp_vector[ String ] &enzyme_info) except + nogil  # wrap-doc:Adds an enzyme to the list and sets is as used

        libcpp_map[ String, libcpp_vector[ String ] ]  getModifications() except + nogil  # wrap-ignore
        void handlePTMs(const String &modification_line, const String &modifications_filename, bool monoisotopic) except + nogil 

