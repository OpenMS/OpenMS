// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/MzTabBase.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <set>

namespace OpenMS
{
  /**
      @brief Data model of MzTabM files.
      Please see the official MzTabM specification at
      https://github.com/HUPO-PSI/mzTab/tree/master/specification_document-releases/2_0-Metabolomics-Release
      @ingroup FileIO
  */

  struct CompareMzTabMMatchRef
  {
    bool operator() (const IdentificationDataInternal::ObservationMatchRef& lhs, const IdentificationDataInternal::ObservationMatchRef& rhs)  const
    {
      return lhs->identified_molecule_var.getIdentifiedCompoundRef()->identifier < rhs->identified_molecule_var.getIdentifiedCompoundRef()->identifier;
    }
  };

  /**
    @brief MztabM Assay Metadata
  */
  class OPENMS_DLLAPI MzTabMAssayMetaData
  {
  public:
    MzTabString name; ///< Name of the assay
    std::map<Size,MzTabParameter> custom; ///< Additional parameters or values for a given assay
    MzTabString external_uri; ///< A reference to further information about the assay
    MzTabInteger sample_ref; ///< An association from a given assay to the sample analysed
    MzTabInteger ms_run_ref; ///< An association from a given assay to the source MS run
  };

  /**
    @brief MztabM MSRun Metadata
  */
  class OPENMS_DLLAPI MzTabMMSRunMetaData
  {
  public:
    MzTabString location; ///< Location of the external data file
    MzTabInteger instrument_ref; ///< Link to a specific instrument
    MzTabParameter format; ///< Parameter specifying the data format of the external MS data file
    MzTabParameter id_format; ///< Parameter specifying the id format used in the external data file
    std::map<Size, MzTabParameter> fragmentation_method; ///< The type of fragmentation used in a given ms run
    std::map<Size, MzTabParameter> scan_polarity; ///< The polarity mode of a given run
    MzTabString hash; ///< Hash value of the corresponding external MS data file
    MzTabParameter hash_method; ///< Parameter specifying the hash methods
  };

  /**
    @brief MztabM StudyVariable Metadata
  */
  class OPENMS_DLLAPI MzTabMStudyVariableMetaData
  {
  public:
    MzTabString name; ///< Name of the study variable
    std::vector<int> assay_refs; ///< References to the IDs of assays grouped in the study variable
    MzTabParameter average_function; ///< The function used to calculate the study variable quantification value
    MzTabParameter variation_function; ///< The function used to calculate the study variable quantification variation value
    MzTabString description; ///< A textual description of the study variable
    MzTabParameterList factors; ///< Additional parameters or factors
  };

  /**
    @brief MztabM Database Metadata
  */
  class OPENMS_DLLAPI MzTabMDatabaseMetaData // mztab-m
  {
  public:
    MzTabParameter database; ///< The description of databases used
    MzTabString prefix; ///< The prefix used in the “identifier” column of data tables
    MzTabString version; ///< The database version
    MzTabString uri; ///< The URI to the database
  };

  /**
    @brief MztabM Metadata
  */
  class OPENMS_DLLAPI MzTabMMetaData
  {
  public:
    MzTabMMetaData();

    MzTabString mz_tab_version; ///< MzTab-M Version
    MzTabString mz_tab_id; ///<  MzTab-M file id (e.g. repository-, local identifier)
    MzTabString title; ///< Title
    MzTabString description; ///< Description
    std::map<Size, MzTabParameterList> sample_processing; ///< List of parameters describing the sample processing/preparation/handling
    std::map<Size, MzTabInstrumentMetaData> instrument; ///< List of parameters describing the instrument
    std::map<Size, MzTabSoftwareMetaData> software; ///< Software used to analyze the data
    std::map<Size, MzTabString> publication; ///< Associated publication(s)
    std::map<Size, MzTabContactMetaData> contact; ///< Contact name
    std::map<Size, MzTabString> uri; ///< Pointing to file source (e.g. MetaboLights)
    std::map<Size, MzTabString> external_study_uri; ///< Pointing to an external file with more details about the study (e.g. ISA-TAB file)
    MzTabParameter quantification_method; ///< Quantification method used in the experiment
    std::map<Size, MzTabSampleMetaData> sample; ///< Sample details
    std::map<Size, MzTabMMSRunMetaData> ms_run; ///< MS run details
    std::map<Size, MzTabMAssayMetaData> assay; ///< Assay details
    std::map<Size, MzTabMStudyVariableMetaData> study_variable; ///< Study Variable details
    std::map<Size, MzTabParameter> custom; ///< Custom parameters
    std::map<Size, MzTabCVMetaData> cv; ///< Controlled Vocabulary details
    std::map<Size, MzTabMDatabaseMetaData> database; ///< Database details
    std::map<Size, MzTabParameter> derivatization_agent; ///< A description of derivatization agents applied to small molecules
    MzTabParameter small_molecule_quantification_unit; ///< Description of the unit type used
    MzTabParameter small_molecule_feature_quantification_unit; ///< Description of the unit type used
    MzTabParameter small_molecule_identification_reliability; ///< Reliability of identification (4-level schema)
    std::map<Size, MzTabParameter> id_confidence_measure; ///< Confidence measures / scores
    std::vector<MzTabString> colunit_small_molecule; ///< Defines the unit used for a specific column
    std::vector<MzTabString> colunit_small_molecule_feature; ///< Defines the unit used for a specific column
    std::vector<MzTabString> colunit_small_molecule_evidence; ///< Defines the unit used for a specific column
  };

  /**
    @brief SML Small molecule section (mztab-m)
  */
  class OPENMS_DLLAPI MzTabMSmallMoleculeSectionRow
  {
  public:
    MzTabString sml_identifier; ///< The small molecule’s identifier.
    MzTabStringList smf_id_refs; ///< References to all the features on which quantification has been based.
    MzTabStringList database_identifier; ///< Names of the used databases.
    MzTabStringList chemical_formula; ///< Potential chemical formula of the reported compound.
    MzTabStringList smiles; ///< Molecular structure in SMILES format.
    MzTabStringList inchi; ///< InChi of the potential compound identifications.
    MzTabStringList chemical_name; ///< Possible chemical/common names or general description
    MzTabStringList uri; ///< The source entry’s location.
    MzTabDoubleList theoretical_neutral_mass; ///< Precursor theoretical neutral mass
    MzTabStringList adducts; ///< Adducts
    // Reliability information of the used indentificavtion method has to be stored in the ID data structure
    MzTabString reliability; ///< Reliability of the given small molecule identification
    MzTabParameter best_id_confidence_measure; ///< The identification approach with the highest confidence
    MzTabDouble best_id_confidence_value; ///< The best confidence measure
    std::map<Size, MzTabDouble> small_molecule_abundance_assay; ///< The small molecule’s abundance in every assay described in the metadata section
    std::map<Size, MzTabDouble> small_molecule_abundance_study_variable; ///< The small molecule’s abundance in all the study variables described in the metadata section
    std::map<Size, MzTabDouble> small_molecule_abundance_variation_study_variable; ///< A measure of the variability of the study variable abundance measurement
    std::vector<MzTabOptionalColumnEntry> opt_; ///< Optional columns must start with “opt_”.
  };

  /**
    @brief SMF Small molecule feature section (mztab-m)
  */
  class OPENMS_DLLAPI MzTabMSmallMoleculeFeatureSectionRow
  {
  public:
    MzTabString smf_identifier; ///< Within file unique identifier for the small molecule feature.
    MzTabStringList sme_id_refs; ///< Reference to the identification evidence.
    MzTabInteger sme_id_ref_ambiguity_code; ///< Ambiguity in identifications.
    MzTabString adduct; ///< Adduct
    MzTabParameter isotopomer; ///< If de-isotoping has not been performed, then the isotopomer quantified MUST be reported here.
    MzTabDouble exp_mass_to_charge; ///< Precursor ion’s m/z.
    MzTabInteger charge; ///< Precursor ion’s charge.
    MzTabDouble retention_time; ///< Time point in seconds.
    MzTabDouble rt_start; ///< The start time of the feature on the retention time axis.
    MzTabDouble rt_end; ///< The end time of the feature on the retention time axis
    std::map<Size, MzTabDouble> small_molecule_feature_abundance_assay; ///< Feature abundance in every assay
    std::vector<MzTabOptionalColumnEntry> opt_; ///< Optional columns must start with “opt_”.
  };

  /**
    @brief SME Small molecule evidence section (mztab-m)
  */
  class OPENMS_DLLAPI MzTabMSmallMoleculeEvidenceSectionRow
  {
  public:
    MzTabString sme_identifier; ///< Within file unique identifier for the small molecule evidence result.
    MzTabString evidence_input_id; ///< Within file unique identifier for the input data used to support this identification e.g. fragment spectrum, RT and m/z pair.
    MzTabString database_identifier; ///< The putative identification for the small molecule sourced from an external database.
    MzTabString chemical_formula; ///< The putative molecular formula.
    MzTabString smiles; ///< Potential molecular structure as SMILES.
    MzTabString inchi; ///< InChi of the potential compound identifications.
    MzTabString chemical_name; ///< Possible chemical/common names or general description
    MzTabString uri; ///< The source entry’s location.
    MzTabParameter derivatized_form; ///< derivatized form.
    MzTabString adduct; ///< Adduct
    MzTabDouble exp_mass_to_charge; ///< Precursor ion’s m/z.
    MzTabInteger charge; ///< Precursor ion’s charge.
    MzTabDouble calc_mass_to_charge; ///< Precursor ion’s m/z.
    MzTabSpectraRef spectra_ref; ///< Reference to a spectrum
    MzTabParameter identification_method; ///< Database search, search engine or process that was used to identify this small molecule
    MzTabParameter ms_level; ///< The highest MS level used to inform identification
    std::map<Size, MzTabDouble> id_confidence_measure; ///< Statistical value or score for the identification
    MzTabInteger rank; ///< Rank of the identification (1 = best)
    std::vector<MzTabOptionalColumnEntry> opt_; ///< Optional columns must start with “opt_”.
  };

  typedef std::vector<MzTabMSmallMoleculeSectionRow> MzTabMSmallMoleculeSectionRows;
  typedef std::vector<MzTabMSmallMoleculeFeatureSectionRow> MzTabMSmallMoleculeFeatureSectionRows;
  typedef std::vector<MzTabMSmallMoleculeEvidenceSectionRow> MzTabMSmallMoleculeEvidenceSectionRows;

  /**
  @brief Data model of MzTab-M files
  Please see the MzTab-M specification at https://github.com/HUPO-PSI/mzTab/blob/master/specification_document-releases/2_0-Metabolomics-Release/mzTab_format_specification_2_0-M_release.adoc#use-cases-for-mztab
  */
  class OPENMS_DLLAPI MzTabM : public MzTabBase
  {
  public:
    /// Default constructor
    MzTabM() = default;

    /// Destructor
    ~MzTabM() = default;

    /// Extract MzTabMMetaData
    const MzTabMMetaData& getMetaData() const;

    /// Set MzTabMMetaData
    void setMetaData(const MzTabMMetaData& m_md);

    /// Extract MzTabMSmallMoleculeSectionRows
    const MzTabMSmallMoleculeSectionRows& getMSmallMoleculeSectionRows() const;

    /// Set MzTabMSmallMoleculeSectionRows
    void setMSmallMoleculeSectionRows(const MzTabMSmallMoleculeSectionRows& m_smsd);

    /// Extract MzTabMSmallMoleculeFeatureSectionRows
    const MzTabMSmallMoleculeFeatureSectionRows& getMSmallMoleculeFeatureSectionRows() const;

    /// Set MzTabMSmallMoleculeFeatureSectionRows
    void setMSmallMoleculeFeatureSectionRows(const MzTabMSmallMoleculeFeatureSectionRows& m_smfsd);

    /// Extract MzTabMSmallMoleculeEvidenceSectionRows
    const MzTabMSmallMoleculeEvidenceSectionRows& getMSmallMoleculeEvidenceSectionRows() const;

    /// Set MzTabMSmallMoleculeEvidenceSectionRows
    void setMSmallMoleculeEvidenceSectionRows(const MzTabMSmallMoleculeEvidenceSectionRows& m_smesd);

    /// Set comment rows
    void setCommentRows(const std::map<Size, String>& com);

    /// Set empty rows
    void setEmptyRows(const std::vector<Size>& empty);

    /// Get empty rows
    const std::vector<Size>& getEmptyRows() const;

    /// Get comment rows
    const std::map<Size, String>& getCommentRows() const;

    /// Extract opt_ (custom, optional column names)
    std::vector<String> getMSmallMoleculeOptionalColumnNames() const;

    /// Extract opt_ (custom, optional column names)
    std::vector<String> getMSmallMoleculeFeatureOptionalColumnNames() const;

    /// Extract opt_ (custom, optional column names)
    std::vector<String> getMSmallMoleculeEvidenceOptionalColumnNames() const;

    static void addMetaInfoToOptionalColumns(const std::set<String>& keys,
                                             std::vector<MzTabOptionalColumnEntry>& opt,
                                             const String& id, const MetaInfoInterface& meta);

  /**
  * @brief Export FeatureMap with Identifications to MzTabM
  *
  * @return MzTabM object
  */
  static MzTabM exportFeatureMapToMzTabM(const FeatureMap& feature_map);

  protected:
    MzTabMMetaData m_meta_data_;
    MzTabMSmallMoleculeSectionRows m_small_molecule_data_;
    MzTabMSmallMoleculeFeatureSectionRows m_small_molecule_feature_data_;
    MzTabMSmallMoleculeEvidenceSectionRows m_small_molecule_evidence_data_;
    std::vector<Size> empty_rows_; ///< index of empty rows
    std::map<Size, String> comment_rows_; ///< comments
    std::vector<String> sml_optional_column_names_;
    std::vector<String> smf_optional_column_names_;
    std::vector<String> sme_optional_column_names_;

    static String getAdductString_(const IdentificationDataInternal::ObservationMatchRef& match_ref);

    static void getFeatureMapMetaValues_(const FeatureMap& feature_map,
                                         std::set<String>& feature_user_value_keys,
                                         std::set<String>& observationmatch_user_value_keys,
                                         std::set<String>& compound_user_value_keys);

  };
} // namespace OpenMS
