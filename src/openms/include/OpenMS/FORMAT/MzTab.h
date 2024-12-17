// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/SVOutStream.h>
#include <OpenMS/FORMAT/MzTabBase.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/METADATA/PeptideEvidence.h>

#include <optional>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wnon-virtual-dtor"

namespace OpenMS
{
/**
      @brief Data model of MzTab files.

      Please see the official MzTab specification at https://code.google.com/p/mztab/

      @ingroup FileIO
  */

  class OPENMS_DLLAPI MzTabModification
  {
public:
    MzTabModification();

    bool isNull() const;

    void setNull(bool b);

    /// set (potentially ambiguous) position(s) with associated parameter (might be null if not set)
    void setPositionsAndParameters(const std::vector<std::pair<Size, MzTabParameter> >& ppp);

    std::vector<std::pair<Size, MzTabParameter> > getPositionsAndParameters() const;

    void setModificationIdentifier(const MzTabString& mod_id);

    MzTabString getModOrSubstIdentifier() const;

    String toCellString() const;

    void fromCellString(const String& s);

    ~MzTabModification() = default;
protected:
    std::vector<std::pair<Size, MzTabParameter> > pos_param_pairs_;
    MzTabString mod_identifier_;
  };

  class OPENMS_DLLAPI MzTabModificationList
  {
public:
    bool isNull() const;

    void setNull(bool b);

    String toCellString() const;

    void fromCellString(const String& s);

    std::vector<MzTabModification> get() const;

    void set(const std::vector<MzTabModification>& entries);

    ~MzTabModificationList() = default;
protected:
    std::vector<MzTabModification> entries_;
  };

// MTD
  struct OPENMS_DLLAPI MzTabModificationMetaData
  {
    MzTabParameter modification;
    MzTabString site;
    MzTabString position;
  };

  struct OPENMS_DLLAPI MzTabAssayMetaData
  {
    MzTabParameter quantification_reagent;
    std::map<Size, MzTabModificationMetaData> quantification_mod;
    MzTabString sample_ref;
    std::vector<int> ms_run_ref; // adapted to address https://github.com/HUPO-PSI/mzTab/issues/26
  };

  struct OPENMS_DLLAPI MzTabMSRunMetaData
  {
    MzTabParameter format;
    MzTabString location;
    MzTabParameter id_format;
    MzTabParameterList fragmentation_method;
  };

  struct OPENMS_DLLAPI MzTabStudyVariableMetaData
  {
    std::vector<int> assay_refs;
    std::vector<int> sample_refs;
    MzTabString description;
  };

  /// all meta data of a mzTab file. Please refer to specification for documentation.
  class OPENMS_DLLAPI MzTabMetaData
  {
public:
    MzTabMetaData();

    MzTabString mz_tab_version;
    MzTabString mz_tab_mode;
    MzTabString mz_tab_type;
    MzTabString mz_tab_id;
    MzTabString title;
    MzTabString description;

    std::map<Size, MzTabParameter> protein_search_engine_score;
    std::map<Size, MzTabParameter> peptide_search_engine_score;
    std::map<Size, MzTabParameter> psm_search_engine_score;
    std::map<Size, MzTabParameter> smallmolecule_search_engine_score;
    std::map<Size, MzTabParameter> nucleic_acid_search_engine_score;
    std::map<Size, MzTabParameter> oligonucleotide_search_engine_score;
    std::map<Size, MzTabParameter> osm_search_engine_score;

    std::map<Size, MzTabParameterList> sample_processing;

    std::map<Size, MzTabInstrumentMetaData> instrument;

    std::map<Size, MzTabSoftwareMetaData> software;

    MzTabParameterList false_discovery_rate;

    std::map<Size, MzTabString> publication;

    std::map<Size, MzTabContactMetaData> contact;

    std::map<Size, MzTabString> uri;

    std::map<Size, MzTabModificationMetaData> fixed_mod;

    std::map<Size, MzTabModificationMetaData> variable_mod;

    MzTabParameter quantification_method;

    MzTabParameter protein_quantification_unit;
    MzTabParameter peptide_quantification_unit;
    MzTabParameter small_molecule_quantification_unit;

    std::map<Size, MzTabMSRunMetaData> ms_run;

    std::map<Size, MzTabParameter> custom;

    std::map<Size, MzTabSampleMetaData> sample;

    std::map<Size, MzTabAssayMetaData> assay;

    std::map<Size, MzTabStudyVariableMetaData> study_variable;

    std::map<Size, MzTabCVMetaData> cv;

    std::vector<String> colunit_protein;
    std::vector<String> colunit_peptide;
    std::vector<String> colunit_psm;
    std::vector<String> colunit_small_molecule;
  };

  /// PRT - Protein section (Table based)
  struct OPENMS_DLLAPI MzTabProteinSectionRow
  {
    MzTabProteinSectionRow();
    MzTabString accession; ///< The protein’s accession.
    MzTabString description; ///< Human readable description (i.e. the name)
    MzTabInteger taxid; ///< NEWT taxonomy for the species.
    MzTabString species; ///< Human readable name of the species
    MzTabString database; ///< Name of the protein database.
    MzTabString database_version; ///< String Version of the protein database.
    MzTabParameterList search_engine; ///< Search engine(s) identifying the protein.
    std::map<Size, MzTabDouble>  best_search_engine_score; ///< best_search_engine_score[1-n]
    std::map<Size, std::map<Size, MzTabDouble> > search_engine_score_ms_run; ///< search_engine_score[index1]_ms_run[index2]
    MzTabInteger reliability;
    std::map<Size, MzTabInteger> num_psms_ms_run;
    std::map<Size, MzTabInteger> num_peptides_distinct_ms_run;
    std::map<Size, MzTabInteger> num_peptides_unique_ms_run;
    MzTabStringList ambiguity_members; ///< Alternative protein identifications.
    MzTabModificationList modifications; ///< Modifications identified in the protein.
    MzTabString uri; ///< Location of the protein’s source entry.
    MzTabStringList go_terms; ///< List of GO terms for the protein.
    MzTabDouble coverage; ///< (0-1) Amount of protein sequence identified.
    std::map<Size, MzTabDouble> protein_abundance_assay;
    std::map<Size, MzTabDouble> protein_abundance_study_variable;
    std::map<Size, MzTabDouble> protein_abundance_stdev_study_variable;
    std::map<Size, MzTabDouble> protein_abundance_std_error_study_variable;
    std::vector<MzTabOptionalColumnEntry> opt_; ///< Optional Columns must start with “opt_”

    /// Comparison operator for sorting rows
    struct RowCompare
    {
      bool operator()(const MzTabProteinSectionRow& row1,
                      const MzTabProteinSectionRow& row2) const
      {
        return row1.accession.get() < row2.accession.get();
      }
    };
  };

  /// PEP - Peptide section (Table based)
  struct OPENMS_DLLAPI MzTabPeptideSectionRow
  {
    MzTabString sequence; ///< The peptide’s sequence.
    MzTabString accession; ///< The protein’s accession.
    MzTabBoolean unique; ///< 0=false, 1=true, null else: Peptide is unique for the protein.
    MzTabString database; ///< Name of the sequence database.
    MzTabString database_version; ///< Version (and optionally # of entries).
    MzTabParameterList search_engine; ///< Search engine(s) that identified the peptide.
    std::map<Size, MzTabDouble> best_search_engine_score; ///< Search engine(s) score(s) for the peptide.
    std::map<Size, std::map<Size, MzTabDouble> > search_engine_score_ms_run;
    MzTabInteger reliability; ///< (1-3) 0=null Identification reliability for the peptide.
    MzTabModificationList modifications; ///< Modifications identified in the peptide.
    MzTabDoubleList retention_time; ///< Time points in seconds. Semantics may vary.
    MzTabDoubleList retention_time_window;
    MzTabInteger charge; ///< Precursor ion’s charge.
    MzTabDouble mass_to_charge; ///< Precursor ion’s m/z.
    MzTabString uri; ///< Location of the PSMs source entry.
    MzTabSpectraRef spectra_ref; ///< Spectra identifying the peptide.
    std::map<Size, MzTabDouble> peptide_abundance_assay;
    std::map<Size, MzTabDouble> peptide_abundance_study_variable;
    std::map<Size, MzTabDouble> peptide_abundance_stdev_study_variable;
    std::map<Size, MzTabDouble> peptide_abundance_std_error_study_variable;
    std::vector<MzTabOptionalColumnEntry> opt_; ///< Optional columns must start with “opt_”.


    /// Comparison operator for sorting rows
    struct RowCompare
    {
      bool operator()(const MzTabPeptideSectionRow& row1,
                      const MzTabPeptideSectionRow& row2) const
      {
        return (std::make_pair(row1.sequence.get(), row1.accession.get()) <
                std::make_pair(row2.sequence.get(), row2.accession.get()));
      }
    };
  };

  /// PSM - PSM section (Table based)
  struct OPENMS_DLLAPI MzTabPSMSectionRow
  {
    MzTabString sequence; ///< The peptide’s sequence.
    MzTabInteger PSM_ID; ///< A unique ID of a PSM line
    MzTabString accession; ///< List of potential parent protein accessions as in the fasta DB.
    MzTabBoolean unique; ///< 0=false, 1=true, null else: Peptide is unique for the protein.
    MzTabString database; ///< Name of the sequence database.
    MzTabString database_version; ///< Version (and optionally # of entries).
    MzTabParameterList search_engine; ///< Search engine(s) that identified the peptide.
    std::map<Size, MzTabDouble> search_engine_score; ///< Search engine(s) score(s) for the peptide.
    MzTabInteger reliability; ///< (1-3) 0=null Identification reliability for the peptide.
    MzTabModificationList modifications; ///< Modifications identified in the peptide.
    MzTabDoubleList retention_time; ///< Time points in seconds. Semantics may vary.
    MzTabInteger charge; ///< The charge of the experimental precursor ion.
    MzTabDouble exp_mass_to_charge; ///< The observed m/z ratio of the experimental precursor ion (either directly from the raw data or corrected).
    MzTabDouble calc_mass_to_charge; ///< The calculated m/z ratio of the experimental precursor ion.
    MzTabString uri; ///< Location of the PSM’s source entry.
    MzTabSpectraRef spectra_ref; ///< Spectrum for this PSM
    MzTabString pre; ///< (List of) Amino acid in parent protein(s) before the start of the current PSM
    MzTabString post; ///< (List of) Amino acid in parent protein(s) after the start of the current PSM
    MzTabString start; ///< (List of) Start positions in parent protein(s)
    MzTabString end; ///< (List of) Start positions in parent protein(s)
    std::vector<MzTabOptionalColumnEntry> opt_; ///< Optional columns must start with “opt_”.

    /**
      @brief Gets peptide_evidences with data from internal structures adds their info to an MzTabPSMSectionRow (pre- or unfilled)

      @param peptide_evidences Vector of PeptideEvidence holding internal data.
    */

    void addPepEvidenceToRows(const std::vector<PeptideEvidence>& peptide_evidences);
    /// Comparison operator for sorting rows
    struct RowCompare
    {
      bool operator()(const MzTabPSMSectionRow& row1,
                      const MzTabPSMSectionRow& row2) const
      {
        // @TODO: sort by "PSM_ID"? what's the point of that field?
        return (std::make_tuple(row1.sequence.get(),
                                row1.spectra_ref.getMSFile(),
                                row1.spectra_ref.getSpecRef(),
                                row1.accession.get()) <
                std::make_tuple(row2.sequence.get(),
                                row2.spectra_ref.getMSFile(),
                                row2.spectra_ref.getSpecRef(),
                                row2.accession.get()));
      }
    };
  };

  /// SML Small molecule section (table based)
  struct OPENMS_DLLAPI MzTabSmallMoleculeSectionRow
  {
    MzTabStringList identifier; ///< The small molecule’s identifier.
    MzTabString chemical_formula; ///< Chemical formula of the identified compound.
    MzTabString smiles; ///< Molecular structure in SMILES format.
    MzTabString inchi_key; ///< InChi Key of the identified compound.
    MzTabString description; ///< Human readable description (i.e. the name)
    MzTabDouble exp_mass_to_charge; ///< Precursor ion’s m/z.
    MzTabDouble calc_mass_to_charge; ///< Precursor ion’s m/z.
    MzTabInteger charge; ///< Precursor ion’s charge.
    MzTabDoubleList retention_time; ///< Time points in seconds. Semantics may vary.
    MzTabInteger taxid; ///< NEWT taxonomy for the species.
    MzTabString species; ///< Human readable name of the species
    MzTabString database; ///< Name of the used database.
    MzTabString database_version; ///< String Version of the database (and optionally # of compounds).
    MzTabInteger reliability; ///< (1-3) The identification reliability.
    MzTabString uri; ///< The source entry’s location.
    MzTabSpectraRef spectra_ref; ///< Spectra identifying the small molecule.
    MzTabParameterList search_engine; ///< Search engine(s) identifying the small molecule.
    std::map<Size, MzTabDouble> best_search_engine_score; ///< Search engine(s) identifications score(s).
    std::map<Size, std::map<Size, MzTabDouble> > search_engine_score_ms_run;
    MzTabString modifications; ///< Modifications identified on the small molecule.
    std::map<Size, MzTabDouble> smallmolecule_abundance_assay;
    std::map<Size, MzTabDouble> smallmolecule_abundance_study_variable;
    std::map<Size, MzTabDouble> smallmolecule_abundance_stdev_study_variable;
    std::map<Size, MzTabDouble> smallmolecule_abundance_std_error_study_variable;
    std::vector<MzTabOptionalColumnEntry> opt_; ///< Optional columns must start with “opt_”.
  };

  /// NUC - Nucleic acid section (table-based)
  struct OPENMS_DLLAPI MzTabNucleicAcidSectionRow
  {
    MzTabString accession; ///< The nucleic acid’s accession.
    MzTabString description; ///< Human readable description (i.e. the name)
    MzTabInteger taxid; ///< NEWT taxonomy for the species.
    MzTabString species; ///< Human readable name of the species
    MzTabString database; ///< Name of the sequence database.
    MzTabString database_version; ///< Version of the sequence database.
    MzTabParameterList search_engine; ///< Search engine(s) that identified the nucleic acid.
    std::map<Size, MzTabDouble>  best_search_engine_score; ///< Best search engine(s) score(s) (over all MS runs)
    std::map<Size, std::map<Size, MzTabDouble> > search_engine_score_ms_run;
    MzTabInteger reliability;
    std::map<Size, MzTabInteger> num_osms_ms_run;
    std::map<Size, MzTabInteger> num_oligos_distinct_ms_run;
    std::map<Size, MzTabInteger> num_oligos_unique_ms_run;
    MzTabStringList ambiguity_members; ///< Alternative nucleic acid identifications.
    MzTabModificationList modifications; ///< Modifications identified in the nucleic acid.
    MzTabString uri; ///< Location of the nucleic acid’s source entry.
    // do GO terms make sense for nucleic acid sequences?
    MzTabStringList go_terms; ///< List of GO terms for the nucleic acid.
    MzTabDouble coverage; ///< (0-1) Fraction of nucleic acid sequence identified.
    std::vector<MzTabOptionalColumnEntry> opt_; ///< Optional Columns must start with “opt_”

    /// Comparison operator for sorting rows
    struct RowCompare
    {
      bool operator()(const MzTabNucleicAcidSectionRow& row1,
                      const MzTabNucleicAcidSectionRow& row2) const
      {
        return row1.accession.get() < row2.accession.get();
      }
    };
  };

  /// OLI - Oligonucleotide section (table-based)
  struct OPENMS_DLLAPI MzTabOligonucleotideSectionRow
  {
    MzTabString sequence; ///< The oligonucleotide’s sequence.
    MzTabString accession; ///< The nucleic acid’s accession.
    MzTabBoolean unique; ///< 0=false, 1=true, null else: Oligonucleotide maps uniquely to the nucleic acid sequence.
    MzTabParameterList search_engine; ///< Search engine(s) that identified the match.
    std::map<Size, MzTabDouble> best_search_engine_score; ///< Search engine(s) score(s) for the match.
    std::map<Size, std::map<Size, MzTabDouble>> search_engine_score_ms_run; ///< Search engine(s) score(s) per individual MS run
    MzTabInteger reliability; ///< (1-3) 0=null Identification reliability for the match.
    MzTabModificationList modifications; ///< Modifications identified in the oligonucleotide.
    MzTabDoubleList retention_time; ///< Time points in seconds. Semantics may vary.
    MzTabDoubleList retention_time_window;
    MzTabString uri; ///< Location of the oligonucleotide's source entry.
    MzTabString pre;
    MzTabString post;
    MzTabInteger start;
    MzTabInteger end;
    std::vector<MzTabOptionalColumnEntry> opt_; ///< Optional columns must start with “opt_”.

    /// Comparison operator for sorting rows
    struct RowCompare
    {
      bool operator()(const MzTabOligonucleotideSectionRow& row1,
                      const MzTabOligonucleotideSectionRow& row2) const
        {
          return (std::make_tuple(row1.sequence.get(), row1.accession.get(),
                                  row1.start.get(), row1.end.get()) <
                  std::make_tuple(row2.sequence.get(), row2.accession.get(),
                                  row2.start.get(), row2.end.get()));
        }
    };

  };

  /// OSM - OSM (oligonucleotide-spectrum match) section (table-based)
  struct OPENMS_DLLAPI MzTabOSMSectionRow
  {
    MzTabString sequence; ///< The oligonucleotide’s sequence.
    MzTabParameterList search_engine; ///< Search engine(s) that identified the match.
    std::map<Size, MzTabDouble> search_engine_score; ///< Search engine(s) score(s) for the match.
    MzTabInteger reliability; ///< (1-3) 0=null Identification reliability for the match.
    MzTabModificationList modifications; ///< Modifications identified in the oligonucleotide.
    MzTabDoubleList retention_time; ///< Time points in seconds. Semantics may vary.
    MzTabInteger charge; ///< The charge of the experimental precursor ion.
    MzTabDouble exp_mass_to_charge; ///< The m/z ratio of the experimental precursor ion.
    MzTabDouble calc_mass_to_charge; ///< The theoretical m/z ratio of the oligonucleotide.
    MzTabString uri; ///< Location of the OSM’s source entry.
    MzTabSpectraRef spectra_ref; ///< Reference to the spectrum underlying the match.
    std::vector<MzTabOptionalColumnEntry> opt_; ///< Optional columns must start with “opt_”.

    /// Comparison operator for sorting rows
    struct RowCompare
    {
      bool operator()(const MzTabOSMSectionRow& row1,
                      const MzTabOSMSectionRow& row2) const
      {
        return (std::make_tuple(row1.sequence.get(),
                                row1.spectra_ref.getMSFile(),
                                row1.spectra_ref.getSpecRef()) <
                std::make_tuple(row2.sequence.get(),
                                row2.spectra_ref.getMSFile(),
                                row2.spectra_ref.getSpecRef()));
      }
    };
  };

  typedef std::vector<MzTabProteinSectionRow> MzTabProteinSectionRows;
  typedef std::vector<MzTabPeptideSectionRow> MzTabPeptideSectionRows;
  typedef std::vector<MzTabPSMSectionRow> MzTabPSMSectionRows;
  typedef std::vector<MzTabSmallMoleculeSectionRow> MzTabSmallMoleculeSectionRows;
  typedef std::vector<MzTabNucleicAcidSectionRow> MzTabNucleicAcidSectionRows;
  typedef std::vector<MzTabOligonucleotideSectionRow> MzTabOligonucleotideSectionRows;
  typedef std::vector<MzTabOSMSectionRow> MzTabOSMSectionRows;


  /**
      @brief Data model of MzTab files.
      Please see the official MzTab specification at https://code.google.com/p/mztab/

      @ingroup FileIO
 */
  class OPENMS_DLLAPI MzTab : public MzTabBase
  {
  public:
    /// Default constructor
    MzTab() = default;
    ~MzTab() = default;

    const MzTabMetaData& getMetaData() const;

    void setMetaData(const MzTabMetaData& md);

    MzTabProteinSectionRows& getProteinSectionRows();

    const MzTabProteinSectionRows& getProteinSectionRows() const;

    void setProteinSectionRows(const MzTabProteinSectionRows& psd);

    MzTabPeptideSectionRows& getPeptideSectionRows();

    const MzTabPeptideSectionRows& getPeptideSectionRows() const;

    void setPeptideSectionRows(const MzTabPeptideSectionRows& psd);

    MzTabPSMSectionRows& getPSMSectionRows();

    const MzTabPSMSectionRows& getPSMSectionRows() const;

    /// Returns the number of PSMs in the PSM section (which is not necessarily the number of rows in the section, due to duplication of rows for each protein)
    /// @note Relies on the PSM_ID to be set correctly for each PSM row
    size_t getNumberOfPSMs() const;

    void setPSMSectionRows(const MzTabPSMSectionRows& psd);

    const MzTabSmallMoleculeSectionRows& getSmallMoleculeSectionRows() const;

    void setSmallMoleculeSectionRows(const MzTabSmallMoleculeSectionRows& smsd);

    const MzTabNucleicAcidSectionRows& getNucleicAcidSectionRows() const;

    void setNucleicAcidSectionRows(const MzTabNucleicAcidSectionRows& nasd);

    const MzTabOligonucleotideSectionRows& getOligonucleotideSectionRows() const;

    void setOligonucleotideSectionRows(const MzTabOligonucleotideSectionRows& onsd);

    const MzTabOSMSectionRows& getOSMSectionRows() const;

    void setOSMSectionRows(const MzTabOSMSectionRows& osd);

    void setCommentRows(const std::map<Size, String>& com);

    void setEmptyRows(const std::vector<Size>& empty);

    const std::vector<Size>& getEmptyRows() const;

    const std::map<Size, String>& getCommentRows() const;

    /// Extract opt_ (custom, optional column names)
    std::vector<String> getProteinOptionalColumnNames() const;

    /// Extract opt_ (custom, optional column names)
    std::vector<String> getPeptideOptionalColumnNames() const;

    /// Extract opt_ (custom, optional column names)
    std::vector<String> getPSMOptionalColumnNames() const;

    /// Extract opt_ (custom, optional column names)
    std::vector<String> getSmallMoleculeOptionalColumnNames() const;

    /// Extract opt_ (custom, optional column names)
    std::vector<String> getNucleicAcidOptionalColumnNames() const;

    /// Extract opt_ (custom, optional column names)
    std::vector<String> getOligonucleotideOptionalColumnNames() const;

    static void addMetaInfoToOptionalColumns(const std::set<String>& keys, std::vector<MzTabOptionalColumnEntry>& opt, const String& id, const MetaInfoInterface& meta);

    /// Extract opt_ (custom, optional column names)
    std::vector<String> getOSMOptionalColumnNames() const;

    static std::map<Size, MzTabModificationMetaData> generateMzTabStringFromModifications(const std::vector<String>& mods);

    static std::map<Size, MzTabModificationMetaData> generateMzTabStringFromVariableModifications(const std::vector<String>& mods);

    static std::map<Size, MzTabModificationMetaData> generateMzTabStringFromFixedModifications(const std::vector<String>& mods);

    static MzTab exportFeatureMapToMzTab(const FeatureMap& feature_map, const String& filename);

    /**
      * @brief Export peptide and protein identifications to mzTab
      *
      * Additionally this function fills two std::maps with mappings for external usage.
      *
      * @param[in] prot_ids Data structure containing protein identifications
      * @param[in] peptide_ids Data structure containing peptide identifications
      * @param[in] filename Input idXML file name
      * @param[in] first_run_inference_only Is all protein inference information (groups and scores) stored in the first run?
      * @param[in] export_empty_pep_ids		Export spectra without PSMs as well?
      * @param[in] export_all_psms		Instead of just the best PSM per spectrum, should other PSMs be exported as well?
      * @param[in] title The title for the metadata section
      * @return mzTab object
    */
    static MzTab exportIdentificationsToMzTab(
        const std::vector<ProteinIdentification>& prot_ids,
        const std::vector<PeptideIdentification>& peptide_ids,
        const String& filename,
        bool first_run_inference_only,
        bool export_empty_pep_ids = false,
        bool export_all_psms = false,
        const String& title = "ID export from OpenMS");


    /// Generate MzTab style list of PTMs from PeptideHit (PSM) object.
    /// All passed fixed modifications are not reported (as suggested by the standard for the PRT and PEP section).
    /// In contrast, all modifications are reported in the PSM section (see standard document for details).
    /// If meta values for modification localization are found, this information is added.
    static MzTabModificationList extractModificationList(const PeptideHit& pep_hit, const std::vector<String>& fixed_mods, const std::vector<String>& localization_mods);

	/**
	 * @brief export linked peptide features aka consensus map
	 *
	 * @param consensus_map		data structure of the linked peptide features
	 * @param filename		input consensusXML file name
   * @param first_run_inference_only Is all protein inference information (groups and scores) stored in the first run?
	 * @param export_unidentified_features		Should not identified peptide features be exported?
	 * @param export_unassigned_ids		Should unassigned identifications be exported?
	 * @param export_subfeatures		The position of the consensus feature will always be exported. Should the individual subfeatures be exported as well?
   * @param export_empty_pep_ids		Export spectra without PSMs as well?
	 * @param export_all_psms		Instead of just the best PSM per spectrum, should other PSMs be exported as well?
   * @param title The title for the metadata section
	 *
	 * @return mzTab object
	 */
    static MzTab exportConsensusMapToMzTab(
      const ConsensusMap& consensus_map,
      const String& filename,
      const bool first_run_inference_only,
      const bool export_unidentified_features,
      const bool export_unassigned_ids,
      const bool export_subfeatures,
      const bool export_empty_pep_ids = false,
      const bool export_all_psms = false,
      const String& title = "ConsensusMap export from OpenMS");

    class IDMzTabStream
    {
       public:
        IDMzTabStream(
          const std::vector<const ProteinIdentification*>& prot_ids,
          const std::vector<const PeptideIdentification*>& peptide_ids,
          const String& filename,
          bool first_run_inference_only,
          bool export_empty_pep_ids = false,
          bool export_all_psms = false,
          const String& title = "ID export from OpenMS");

         const MzTabMetaData& getMetaData() const;

         const std::vector<String>& getProteinOptionalColumnNames() const; 
         const std::vector<String>& getPeptideOptionalColumnNames() const;
         const std::vector<String>& getPSMOptionalColumnNames() const;

         bool nextPRTRow(MzTabProteinSectionRow& row);
         bool nextPEPRow(MzTabPeptideSectionRow& row);
         bool nextPSMRow(MzTabPSMSectionRow& row);
       private:
         std::set<String> protein_hit_user_value_keys_;
         std::set<String> peptide_id_user_value_keys_;
         std::set<String> peptide_hit_user_value_keys_;

         // beautiful mapping structs
         std::map<Size, std::set<Size>> ind2prot_;
         std::map<Size, std::set<Size>> pg2prot_;
         std::map<String, size_t> idrunid_2_idrunindex_;
         std::map<Size, std::vector<std::pair<String, String>>> run_to_search_engines_;
         std::map<Size, std::vector<std::vector<std::pair<String, String>>>> run_to_search_engines_settings_;
         std::map<std::pair<size_t,size_t>,size_t> map_id_run_fileidx_2_msfileidx_;
         std::map<std::pair< String, unsigned >, unsigned> path_label_to_assay_;

         std::vector<const ProteinIdentification*> prot_ids_;
         std::vector<const PeptideIdentification*> peptide_ids_;

         StringList ms_runs_;
         bool first_run_inference_;
         String filename_;
         StringList fixed_mods_;
         /* currently unused
         bool export_unidentified_features_; 
         bool export_subfeatures_; */
         bool export_empty_pep_ids_;
         bool export_all_psms_;
         size_t quant_study_variables_ = 0;
         // size_t n_study_variables_ = 0; //currently unused
         size_t PRT_STATE_ = 0;
         size_t prt_run_id_ = 0; // current (protein) identification run
         size_t prt_hit_id_ = 0; // current protein in (protein) identification run
         size_t prt_group_id_ = 0;
         size_t prt_indistgroup_id_ = 0;
         size_t pep_id_ = 0;
         size_t psm_id_ = 0;
         size_t current_psm_idx_ = 0;
         MzTabString db_, db_version_;

         std::vector<String> prt_optional_column_names_;
         std::vector<String> pep_optional_column_names_;
         std::vector<String> psm_optional_column_names_;

         MzTabMetaData meta_data_;
    };

    class CMMzTabStream
    {
       public:
        CMMzTabStream(
          const ConsensusMap& consensus_map,
          const String& filename,
          const bool first_run_inference_only,
          const bool export_unidentified_features,
          const bool export_unassigned_ids,
          const bool export_subfeatures,
          const bool export_empty_pep_ids = false,
          const bool export_all_psms = false,
          const String& title = "ConsensusMap export from OpenMS");

         const MzTabMetaData& getMetaData() const;

         const std::vector<String>& getProteinOptionalColumnNames() const; 
         const std::vector<String>& getPeptideOptionalColumnNames() const;
         const std::vector<String>& getPSMOptionalColumnNames() const;

         bool nextPRTRow(MzTabProteinSectionRow& row);
         bool nextPEPRow(MzTabPeptideSectionRow& row);
         bool nextPSMRow(MzTabPSMSectionRow& row);

       private:
         const ConsensusMap& consensus_map_;
         std::set<String> protein_hit_user_value_keys_;
         std::set<String> consensus_feature_user_value_keys_;
         std::set<String> consensus_feature_peptide_identification_user_value_keys_;
         std::set<String> consensus_feature_peptide_hit_user_value_keys_;

         // beautiful mapping structs
         std::map<Size, std::set<Size>> ind2prot_;
         std::map<Size, std::set<Size>> pg2prot_;
         std::map<String, size_t> idrunid_2_idrunindex_;
         std::map<Size, std::vector<std::pair<String, String>>> run_to_search_engines_;
         std::map<Size, std::vector<std::vector<std::pair<String, String>>>> run_to_search_engines_settings_;
         std::map<std::pair<size_t,size_t>,size_t> map_id_run_fileidx_2_msfileidx_;
         std::map<std::pair< String, unsigned >, unsigned> path_label_to_assay_;

         std::vector<const ProteinIdentification*> prot_ids_;
         std::vector<const PeptideIdentification*> peptide_ids_;

         StringList ms_runs_;
         bool first_run_inference_;
         String filename_;
         StringList fixed_mods_;
         bool export_unidentified_features_; 
         bool export_subfeatures_;
         bool export_empty_pep_ids_;
         bool export_all_psms_;
         size_t quant_study_variables_ = 0;
         size_t n_study_variables_ = 0;
         size_t PRT_STATE_ = 0;
         size_t prt_run_id_ = 0; // current (protein) identification run
         size_t prt_hit_id_ = 0; // current protein in (protein) identification run
         size_t prt_group_id_ = 0;
         size_t prt_indistgroup_id_ = 0;
         size_t pep_id_ = 0;
         size_t pep_counter_ = 0;
         size_t psm_id_ = 0;
         size_t current_psm_idx_ = 0;
         MzTabString db_, db_version_;

         std::vector<String> prt_optional_column_names_;
         std::vector<String> pep_optional_column_names_;
         std::vector<String> psm_optional_column_names_;

         MzTabMetaData meta_data_;
    };

     
  protected:
    // extract basic mappings

    static std::map<String, Size> mapIDRunIdentifier2IDRunIndex_(const std::vector<const ProteinIdentification*>& prot_ids);

    static std::optional<MzTabPSMSectionRow> PSMSectionRowFromPeptideID_(
      PeptideIdentification const& pid,
      std::vector<ProteinIdentification const*> const& prot_id,
      std::map<String, size_t>& idrun_2_run_index,
      std::map<std::pair<size_t, size_t>, size_t>& map_run_fileidx_2_msfileidx,
      std::map<Size, std::vector<std::pair<String, String>>>& run_to_search_engines,
      Size const current_psm_idx,
      Size const psm_id,
      MzTabString const& db,
      MzTabString const& db_version,
      bool const export_empty_pep_ids,
      bool const export_all_psms);

    static MzTabPeptideSectionRow peptideSectionRowFromConsensusFeature_(
      const ConsensusFeature& c, 
      const ConsensusMap& consensus_map,
      const StringList& ms_runs,
      const Size n_study_variables,
      const std::set<String>& consensus_feature_user_value_keys,
      const std::set<String>& peptide_identifications_user_value_keys,
      const std::set<String>& peptide_hit_user_value_keys,
      const std::map<String, size_t>& idrun_2_run_index,
      const std::map<std::pair<size_t,size_t>,size_t>& map_run_fileidx_2_msfileidx,
      const std::map< std::pair< String, unsigned >, unsigned>& path_label_to_assay,
      const std::vector<String>& fixed_mods,
      bool export_subfeatures);

    static MzTabPeptideSectionRow peptideSectionRowFromFeature_(
      const Feature& c, 
      const std::set<String>& feature_user_value_keys,
      const std::set<String>& peptide_identifications_user_value_keys,
      const std::set<String>& peptide_hit_user_value_keys,
      const std::vector<String>& fixed_mods);

    static MzTabProteinSectionRow proteinSectionRowFromProteinHit_(
      const ProteinHit& hit,
      const MzTabString& db,
      const MzTabString& db_version,
      const std::set<String>& protein_hit_user_value_keys);

    static MzTabProteinSectionRow nextProteinSectionRowFromProteinGroup_(
      const ProteinIdentification::ProteinGroup& group,
      const MzTabString& db,
      const MzTabString& db_version);

    static MzTabProteinSectionRow nextProteinSectionRowFromIndistinguishableGroup_(
      const std::vector<ProteinHit>& protein_hits,
      const ProteinIdentification::ProteinGroup& group,
      const size_t g,
      const std::map<Size, std::set<Size>>& ind2prot,
      const MzTabString& db,
      const MzTabString& db_version);

    static void addMSRunMetaData_(
      const std::map<size_t, String>& msrunindex_2_msfilename,
      MzTabMetaData& meta_data);

    static void mapBetweenMSFileNameAndMSRunIndex_(
      const std::vector<const ProteinIdentification*>& prot_ids, 
      bool skip_first, 
      std::map<String, size_t>& msfilename_2_msrunindex,
      std::map<size_t, String>& msrunindex_2_msfilename);

    static size_t getQuantStudyVariables_(const ProteinIdentification& pid);

    static MzTabParameter getProteinScoreType_(const ProteinIdentification& prot_id);

    // TODO: move to core classes?
    static void getConsensusMapMetaValues_(const ConsensusMap& consensus_map, 
      std::set<String>& consensus_feature_user_value_keys,
      std::set<String>& peptide_identification_user_value_keys,
      std::set<String>& peptide_hit_user_value_keys);

    static void getFeatureMapMetaValues_(const FeatureMap& feature_map,
      std::set<String>& feature_user_value_keys,
      std::set<String>& peptide_identification_user_value_keys,
      std::set<String>& peptide_hit_user_value_keys);

    static void getIdentificationMetaValues_(
      const std::vector<const ProteinIdentification*>& prot_ids, 
      std::vector<const PeptideIdentification*>& peptide_ids_,
      std::set<String>& protein_hit_user_value_keys,
      std::set<String>& peptide_id_user_value_keys,
      std::set<String>& peptide_hit_user_value_keys);

    // determine spectrum reference identifier type (e.g., Thermo nativeID) from spectrum references
    static MzTabParameter getMSRunSpectrumIdentifierType_(const std::vector<const PeptideIdentification*>& peptide_ids_);

    static void mapBetweenRunAndSearchEngines_(
      const std::vector<const ProteinIdentification*>& prot_ids,
      const std::vector<const PeptideIdentification*>& pep_ids,
      bool skip_first_run,
      std::map<std::tuple<String, String, String>, std::set<Size>>& search_engine_to_runs,
      std::map<Size, std::vector<std::pair<String, String>>>& run_to_search_engines,
      std::map<Size, std::vector<std::vector<std::pair<String, String>>>>& run_to_search_engines_settings,
      std::map<String, std::vector<std::pair<String, String>>>& search_engine_to_settings);

    static std::map<Size, std::set<Size>> mapGroupsToProteins_(
      const std::vector<ProteinIdentification::ProteinGroup>& groups, 
      const std::vector<ProteinHit>& proteins);

    static void addSearchMetaData_(
        const std::vector<const ProteinIdentification*>& prot_ids,
        const std::map<std::tuple<String, String, String>, std::set<Size>>& search_engine_to_runs,
        const std::map<String, std::vector<std::pair<String,String>>>& search_engine_to_settings,
        MzTabMetaData& meta_data,
        bool first_run_inference_only);

    static void mapIDRunFileIndex2MSFileIndex_(
      const std::vector<const ProteinIdentification*>& prot_ids,
      const std::map<String, size_t>& msfilename_2_msrunindex,
      bool skip_first_run, 
      std::map<std::pair<size_t, size_t>, size_t>& map_run_fileidx_2_msfileidx);

    static void getSearchModifications_(
      const std::vector<const ProteinIdentification*>& prot_ids,
      StringList& var_mods, 
      StringList& fixed_mods);

    // create MzTab compatible modification identifier from ResidueModification
    // If the Modification has a unimod identifier it will be prefixed as UNIMOD
    // otherwise as CHEMMOD (see MzTab specification for details)
    static MzTabString getModificationIdentifier_(const ResidueModification& r);

    static void checkSequenceUniqueness_(const std::vector<PeptideIdentification>& curr_pep_ids);

    MzTabMetaData meta_data_;
    MzTabProteinSectionRows protein_data_;
    MzTabPeptideSectionRows peptide_data_;
    MzTabPSMSectionRows psm_data_;
    MzTabSmallMoleculeSectionRows small_molecule_data_;
    MzTabNucleicAcidSectionRows nucleic_acid_data_;
    MzTabOligonucleotideSectionRows oligonucleotide_data_;
    MzTabOSMSectionRows osm_data_; ///</ oligonucleotide-spectrum matches
    std::vector<Size> empty_rows_; ///< index of empty rows
    std::map<Size, String> comment_rows_; ///< comments
  };

} // namespace OpenMS

#pragma clang diagnostic pop
