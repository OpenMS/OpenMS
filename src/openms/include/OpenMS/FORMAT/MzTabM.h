// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

  class OPENMS_DLLAPI MzTabMAssayMetaData
  {
  public:
    MzTabString name;
    MzTabParameter custom; // mztab-m
    MzTabString external_uri; // mztab-m
    std::vector<int> sample_ref;
    std::vector<int> ms_run_ref; // adapted to address https://github.com/HUPO-PSI/mzTab/issues/26
  };

  class OPENMS_DLLAPI MzTabMMSRunMetaData
  {
  public:
    MzTabString location;
    MzTabInteger instrument_ref; // mztab-m
    MzTabParameter format;
    MzTabParameter id_format;
    MzTabParameterList fragmentation_method;
    MzTabParameterList scan_polarity; // mztab-m
    MzTabString hash; // mztab-m
    MzTabParameter hash_method; // mztab-m
  };

  class OPENMS_DLLAPI MzTabMStudyVariableMetaData
  {
  public:
    MzTabString name;
    std::vector<int> assay_refs;
    MzTabParameter average_function; // mztab-m
    MzTabParameter variation_function; // mztab-m
    MzTabString description;
    MzTabParameterList factors; // mztab-m
  };

  class OPENMS_DLLAPI MzTabMDatabaseMetaData // mztab-m
  {
  public:
    MzTabParameter database;
    MzTabString prefix;
    MzTabString version;
    MzTabString uri;
  };

  /// Metadata for MzTab-M
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
    std::map<Size, MzTabParameter> derivatization_agent; ///<
    MzTabParameter small_molecule_quantification_unit; ///< Description of the unit type used
    MzTabParameter small_molecule_feature_quantification_unit; ///< Description of the unit type used
    MzTabParameter small_molecule_identification_reliability; ///< Reliability of identification (4-level schema)
    std::map<Size, MzTabParameter> id_confidence_measure; ///< Confidence measures / scores
    // https://github.com/HUPO-PSI/mzTab/blob/master/specification_document-releases/2_0-Metabolomics-Release/mzTab_format_specification_2_0-M_release.adoc#6260-colunit-small_molecule_feature
    // TODO: Not sure how that would be best encoded?
    // for a specific optional column?
    // Should be in optional column metadata?
    // e.g. opt_global_mass_error=[UO, UO:0000169, parts per million, ]
    std::vector<MzTabString> colunit_small_molecule; ///< Defines the unit used for a specific column
    std::vector<MzTabString> colunit_small_molecule_feature; ///< Defines the unit used for a specific column
    std::vector<MzTabString> colunit_small_molecule_evidence; ///< Defines the unit used for a specific column
  };

  /// SML Small molecule section (mztab-m)
  class OPENMS_DLLAPI MzTabMSmallMoleculeSectionRow
  {
  public:
    MzTabInteger identifier; ///< The small molecule’s identifier.
    MzTabStringList smf_id_refs; ///< References to all the features on which quantification has been based.
    MzTabStringList database_identifier; ///< Names of the used databases.
    MzTabStringList chemical_formula; ///< Potential chemical formula of the reported compound.
    MzTabStringList smiles; ///< Molecular structure in SMILES format.
    MzTabStringList inchi; ///< InChi of the potential compound identifications.
    MzTabStringList chemical_name; ///< Possible chemical/common names or general description
    MzTabStringList uri; ///< The source entry’s location. // TODO: URI List ?

    MzTabDoubleList theoretical_neutral_mass; ///< Precursor theoretical neutral mass
    MzTabStringList adducts; ///< Adducts
    // TODO: https://github.com/HUPO-PSI/mzTab/blob/master/specification_document-releases/2_0-Metabolomics-Release/mzTab_format_specification_2_0-M_release.adoc#6311-reliability
    MzTabString reliability; ///< Reliability of the given small molecule identification
    // TODO: e.g. use best search_engine score
    MzTabParameter best_id_confidence_measure; ///< The identification approach with the highest confidence
    MzTabDouble best_id_confidence_value; ///< The best confidence measure

    std::map<Size, MzTabDouble> small_molecule_abundance_assay; ///<
    std::map<Size, MzTabDouble> small_molecule_abundance_study_variable; ///<
    std::map<Size, MzTabDouble> small_molecule_abundance_stdev_study_variable; ///<
    std::map<Size, MzTabDouble> small_molecule_abundance_std_error_study_variable; ///<
    std::vector<MzTabOptionalColumnEntry> opt_; ///< Optional columns must start with “opt_”.
  };

  /// SMF Small molecule feature section (mztab-m)
  class OPENMS_DLLAPI MzTabMSmallMoleculeFeatureSectionRow
  {
  public:
    MzTabInteger smf_identifier; ///< Within file unique identifier for the small molecule feature.
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

  /// SME Small molecule evidence section (mztab-m)
  class OPENMS_DLLAPI MzTabMSmallMoleculeEvidenceSectionRow
  {
  public:
    MzTabInteger sme_identifier; ///< Within file unique identifier for the small molecule evidence result.
    MzTabString evidence_input_id; ///< Within file unique identifier for the input data used to support this identification e.g. fragment spectrum, RT and m/z pair.
    MzTabString database_identifier; ///< The putative identification for the small molecule sourced from an external database.
    MzTabString chemical_formula; ///< The putative molecular formula.
    MzTabString smiles; ///< Potential molecular structure as SMILES.
    MzTabString inchi; ///< InChi of the potential compound identifications.
    MzTabString chemical_name; ///< Possible chemical/common names or general description
    MzTabString uri; ///< The source entry’s location.
    MzTabParameter derivatized_form; ///< //TODO: What has to be added here?
    MzTabString adduct; ///< Adduct
    MzTabDouble exp_mass_to_charge; ///< Precursor ion’s m/z.
    MzTabInteger charge; ///< Precursor ion’s charge.
    MzTabDouble calc_mass_to_charge; ///< Precursor ion’s m/z.
    MzTabStringList spectra_ref; ///< Reference to a spectrum
    MzTabParameter identification_method; ///< Database search, search engine or process that was used to identify this small molecule
    MzTabParameter ms_level; ///< The highest MS level used to inform identification
    MzTabDouble id_confidence_measure; ///< Statistical value or score for the identification
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

    const MzTabMMetaData& getMetaData() const;

    void setMetaData(const MzTabMMetaData& m_md); // TODO: check if the metadata section is the same or if additional / other stuff is needed as well

    const MzTabMSmallMoleculeSectionRows& getMSmallMoleculeSectionRows() const;

    void setMSmallMoleculeSectionRows(const MzTabMSmallMoleculeSectionRows& m_smsd);

    const MzTabMSmallMoleculeFeatureSectionRows& getMSmallMoleculeFeatureSectionRows() const;

    void setMSmallMoleculeFeatureSectionRows(const MzTabMSmallMoleculeFeatureSectionRows& m_smfsd);

    const MzTabMSmallMoleculeEvidenceSectionRows& getMSmallMoleculeEvidenceSectionRows() const;

    void setMSmallMoleculeEvidenceSectionRows(const MzTabMSmallMoleculeEvidenceSectionRows& m_smesd);

    void setCommentRows(const std::map<Size, String>& com);

    void setEmptyRows(const std::vector<Size>& empty);

    const std::vector<Size>& getEmptyRows() const;

    const std::map<Size, String>& getCommentRows() const;

    /// Extract opt_ (custom, optional column names)
    std::vector<String> getMSmallMoleculeOptionalColumnNames() const;

    /// Extract opt_ (custom, optional column names)
    std::vector<String> getMSmallMoleculeFeatureOptionalColumnNames() const;

    /// Extract opt_ (custom, optional column names)
    std::vector<String> getMSmallMoleculeEvidenceOptionalColumnNames() const;

    static void addMetaInfoToOptionalColumns(const std::set<String>& keys, std::vector<MzTabOptionalColumnEntry>& opt, const String& id, const MetaInfoInterface& meta);

  /**
  * @brief Export metabolite identifications to MzTabM
  *
  * @return MzTabM object
  */
  static MzTabM exportFeatureMapToMzTabM(const FeatureMap& feature_map, const String& filename);

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
  };

} // namespace OpenMS
