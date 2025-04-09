// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow, Hendrik Weisser $
// $Authors: Chris Bielow, Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/SpectrumMetaDataLookup.h>

#include <vector>
#include <map>
#include <set>


namespace OpenMS
{
  /**
    @brief Used to load and store PepXML files

    This class is used to load and store documents that implement the schema of PepXML files.

    A documented schema for this format comes with the TPP and can also be
    found at https://github.com/OpenMS/OpenMS/tree/develop/share/OpenMS/SCHEMAS

    @ingroup FileIO
  */
  class OPENMS_DLLAPI PepXMLFile :
    protected Internal::XMLHandler,
    public Internal::XMLFile
  {
public:

    /// Constructor
    PepXMLFile();

    /// Destructor
    ~PepXMLFile() override;

    /**
        @brief Loads peptide sequences with modifications out of a PepXML file

        @param filename PepXML file to load
        @param proteins Protein identification output
        @param peptides Peptide identification output
        @param experiment_name Experiment file name, which is used to extract the corresponding search results from the PepXML file.
        @param lookup Helper for looking up retention times (PepXML may contain only scan numbers).

        @exception Exception::FileNotFound is thrown if the file could not be opened
        @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String& filename,
              std::vector<ProteinIdentification>& proteins,
              std::vector<PeptideIdentification>& peptides,
              const String& experiment_name,
              const SpectrumMetaDataLookup& lookup);

    /**
        @brief @a load function with empty defaults for some parameters (see above)

        @exception Exception::FileNotFound is thrown if the file could not be opened
        @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String& filename,
              std::vector<ProteinIdentification>& proteins,
              std::vector<PeptideIdentification>& peptides,
              const String& experiment_name = "");

    /**
        @brief Stores idXML as PepXML file

        @exception Exception::UnableToCreateFile is thrown if the file could not be opened for writing
    */
    void store(const String& filename, std::vector<ProteinIdentification>& protein_ids,
               std::vector<PeptideIdentification>& peptide_ids, const String& mz_file = "",
               const String& mz_name = "", bool peptideprophet_analyzed = false, double rt_tolerance = 0.01);

    /**
        @brief Whether we should keep the native spectrum name of the pepXML

        @note This will lead to a "pepxml_spectrum_name" meta value being added
        to each PeptideIdentification containing the original name of the
        spectrum in TPP format.
    */
    void keepNativeSpectrumName(bool keep)
    {
      keep_native_name_ = keep;
    }

    /// sets the preferred fixed modifications
    void setPreferredFixedModifications(const std::vector<const ResidueModification*>& mods);

    /// sets the preferred variable modifications
    void setPreferredVariableModifications(const std::vector<const ResidueModification*>& mods);

    /// sets if during load, unknown scores should be parsed
    void setParseUnknownScores(bool parse_unknown_scores);

protected:

    /// Docu in base class
    void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname) override;

    /// Docu in base class
    void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes) override;

private:

    /// Fill @p scan_map_
    void makeScanMap_();

    /// Read RT, m/z, charge information from attributes of "spectrum_query"
    void readRTMZCharge_(const xercesc::Attributes& attributes);

    struct AminoAcidModification
    {
      private:

      String aminoacid_;
      double massdiff_;
      double mass_;
      bool is_variable_;
      String description_;
      String terminus_;
      bool is_protein_terminus_; // "true" if protein terminus, "false" if peptide terminus
      ResidueModification::TermSpecificity term_spec_;
      std::vector<String> errors_;
      const ResidueModification* registered_mod_;

      const ResidueModification* lookupModInPreferredMods_(const std::vector<const ResidueModification*>& preferred_fixed_mods,
                                                           const String& aminoacid,
                                                           double massdiff,
                                                           const String& description,
                                                           const ResidueModification::TermSpecificity term_spec,
                                                           double tolerance);

      public:
      AminoAcidModification() = delete;

      /// Creates an AminoAcidModification object from the pepXML attributes in
      /// EITHER aminoacid_modification elements
      /// OR terminal_modification elements
      /// since we use them ambiguously
      AminoAcidModification(
          const String& aminoacid, const String& massdiff, const String& mass,
          String variable, const String& description, String terminus, const String& protein_terminus,
          const std::vector<const ResidueModification*>& preferred_fixed_mods,
          const std::vector<const ResidueModification*>& preferred_var_mods,
          double tolerance);

      AminoAcidModification(const AminoAcidModification& rhs) = default;

      virtual ~AminoAcidModification() = default;

      AminoAcidModification& operator=(const AminoAcidModification& rhs) = default;

      String toUnimodLikeString() const;

      const String& getDescription() const;

      bool isVariable() const;

      const ResidueModification* getRegisteredMod() const;

      double getMassDiff() const;

      double getMass() const;

      const String& getTerminus() const;

      const String& getAminoAcid() const;

      const std::vector<String>& getErrors() const;
    };

    /// Pointer to the list of identified proteins
    std::vector<ProteinIdentification>* proteins_;

    /// Pointer to the list of identified peptides
    std::vector<PeptideIdentification>* peptides_;

    /// Pointer to wrapper for looking up spectrum meta data
    const SpectrumMetaDataLookup* lookup_;

    /// Name of the associated experiment (filename of the data file, extension will be removed)
    String exp_name_;

    /// Set name of search engine
    String search_engine_;

    /// Several optional attributes of spectrum_query
    String native_spectrum_name_;
    String experiment_label_;
    String swath_assay_;
    String status_;

    /// Get RT and m/z for peptide ID from precursor scan (should only matter for RT)?
    bool use_precursor_data_{};

    /// Mapping between scan number in the pepXML file and index in the corresponding MSExperiment
    std::map<Size, Size> scan_map_;

    /// Hydrogen data (for mass types)
    Element hydrogen_;

    /// Are we currently in an "analysis_summary" element (should be skipped)?
    bool analysis_summary_;

    /// Whether we should keep the native spectrum name of the pepXML
    bool keep_native_name_;

    /// Are we currently in an "search_score_summary" element (should be skipped)?
    bool search_score_summary_;

    /// Are we currently in an "search_summary" element (should be skipped)?
    bool search_summary_{};

    /// Do current entries belong to the experiment of interest (for pepXML files that bundle results from different experiments)?
    bool wrong_experiment_{};

    /// Have we seen the experiment of interest at all?
    bool seen_experiment_{};

    /// Have we checked the "base_name" attribute in the "msms_run_summary" element?
    bool checked_base_name_{};

    /// Does the file have decoys (e.g. from Comet's internal decoy search)
    bool has_decoys_{};

    /// Also parse unknown scores as metavalues?
    bool parse_unknown_scores_{};

    /// In case it has decoys, what is the prefix?
    String decoy_prefix_;

    /// current base name
    String current_base_name_;

    /// References to currently active ProteinIdentifications
    std::vector<std::vector<ProteinIdentification>::iterator> current_proteins_;

    /// Search parameters of the current identification run
    ProteinIdentification::SearchParameters params_;

    /// Enzyme name associated with the current identification run
    String enzyme_;
    String enzyme_cuttingsite_;

    /// PeptideIdentification instance currently being processed
    PeptideIdentification current_peptide_;

    /// Analysis result instance currently being processed
    PeptideHit::PepXMLAnalysisResult current_analysis_result_;

    /// PeptideHit instance currently being processed
    PeptideHit peptide_hit_;

    /// Sequence of the current peptide hit
    String current_sequence_;

    /// RT and m/z of current PeptideIdentification (=spectrum)  
    double rt_{}, mz_{};

    /// 1-based scan nr. of current PeptideIdentification (=spectrum). Scannr is usually from the start_scan attribute
    Size scannr_{};

    /// Precursor ion charge
    Int charge_{};

    /// ID of current search result
    UInt search_id_{};

    /// Identifier linking PeptideIdentifications and ProteinIdentifications
    String prot_id_;

    /// Date the pepXML file was generated
    DateTime date_;

    /// Mass of a hydrogen atom (monoisotopic/average depending on case)
    double hydrogen_mass_{};

    /// The modifications of the current peptide hit (position is 1-based)
    std::vector<std::pair<const ResidueModification*, Size> > current_modifications_;

    /// Fixed aminoacid modifications as parsed from the header
    std::vector<AminoAcidModification> fixed_modifications_;

    /// Variable aminoacid modifications as parsed from the header
    std::vector<AminoAcidModification> variable_modifications_;

    /// Fixed modifications that should be preferred when parsing the header
    /// (e.g. when pepXML was produced through an adapter)
    std::vector<const ResidueModification*> preferred_fixed_modifications_;

    /// Variable modifications that should be preferred when parsing the header
    /// (e.g. when pepXML was produced through an adapter)
    std::vector<const ResidueModification*> preferred_variable_modifications_;

    //@}

    static const double mod_tol_;
    static const double xtandem_artificial_mod_tol_;

    /// looks up modification by @p modification_mass and aminoacid of current_sequence_[ @p modification_position ]
    /// and adds it to the current_modifications_
    bool lookupAddFromHeader_(double modification_mass,
                              Size modification_position,
                              std::vector<AminoAcidModification> const& header_mods);

    //static std::vector<int> getIsotopeErrorsFromIntSetting_(int intSetting);
  };
} // namespace OpenMS
