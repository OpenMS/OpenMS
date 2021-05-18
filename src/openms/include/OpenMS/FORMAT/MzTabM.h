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

#include <OpenMS/FORMAT/MzTab.h>

namespace OpenMS
{

  /**
      @brief Data model of MzTabM files.
      Please see the official MzTabM specification at
      https://github.com/HUPO-PSI/mzTab/tree/master/specification_document-releases/2_0-Metabolomics-Release
      @ingroup FileIO
  */

  struct OPENMS_DLLAPI MzTabMAssayMetaData
  {
    MzTabParameter quantification_reagent;
    std::map<Size, MzTabModificationMetaData> quantification_mod;
    MzTabString sample_ref;
    std::vector<int> ms_run_ref; // adapted to address https://github.com/HUPO-PSI/mzTab/issues/26
    MzTabParameter custom; // mztab-m
    MzTabString external_uri; // mztab-m
  };

  struct OPENMS_DLLAPI MzTabMMSRunMetaData
  {
    MzTabParameter format;
    MzTabString location;
    MzTabParameter id_format;
    MzTabParameterList fragmentation_method;
    MzTabInteger instrument_ref; // mztab-m
    MzTabParameter scan_polarity; // mztab-m
    MzTabString hash; // mztab-m
    MzTabParameter hash_method; // mztab-m
  };

    struct OPENMS_DLLAPI MzTabMStudyVariableMetaData
  {
    std::vector<int> assay_refs;
    std::vector<int> sample_refs;
    MzTabString description;
    MzTabParameter average_function; // mztab-m
    MzTabParameter variation_function; // mztab-m
    MzTabParameterList factors; // mztab-m
  };

  struct OPENMS_DLLAPI MzTabMDatabaseMetaData // mztab-m
  {
    MzTabParameter database;
    MzTabString prefix;
    MzTabString version;
    MzTabString uri;
  };

  /// all meta data of a mzTab file. Please refer to specification for documentation.
  // TODO: Check https://github.com/HUPO-PSI/mzTab/blob/master/specification_document-releases/2_0-Metabolomics-Release/mzTab_format_specification_2_0-M_release.adoc#62-metadata-section
  // TODO: Check if the Type used here and the Type in the specification are the same
  class OPENMS_DLLAPI MzTabMMetaData
  {
  public:
    MzTabMMetaData();

    MzTabString mz_tab_version;
    MzTabString mz_tab_id;
    MzTabString title;
    MzTabString description;
    std::map<Size, MzTabParameterList> sample_processing;
    std::map<Size, MzTabInstrumentMetaData> instrument;
    std::map<Size, MzTabSoftwareMetaData> software;
    std::map<Size, MzTabString> publication;
    std::map<Size, MzTabContactMetaData> contact;
    std::map<Size, MzTabString> uri;
    std::map<Size, MzTabString> external_study_uri;
    MzTabParameter quantification_method;
    std::map<Size, MzTabSampleMetaData> sample;
    std::map<Size, MzTabMMSRunMetaData> ms_run;
    std::map<Size, MzTabMAssayMetaData> assay;
    std::map<Size, MzTabMStudyVariableMetaData> study_variable;
    std::map<Size, MzTabParameter> custom;
    std::map<Size, MzTabCVMetaData> cv;
    std::map<Size, MzTabMDatabaseMetaData> database;
    std::map<Size, MzTabParameter> derivatization_agent;
    MzTabParameter small_molecule_quantification_unit;
    MzTabParameter small_molecule_feature_quantification_unit;
    MzTabParameter small_molecule_identification_reliability;
    std::map<Size, MzTabParameter> id_confidence_measure; // TODO: (ADD)
    // https://github.com/HUPO-PSI/mzTab/blob/master/specification_document-releases/2_0-Metabolomics-Release/mzTab_format_specification_2_0-M_release.adoc#6260-colunit-small_molecule_feature
    std::vector<MzTabString> colunit_small_molecule; // TODO: ?
    std::vector<MzTabString> colunit_small_molecule_feature; // TODO: ?
    std::vector<MzTabString> colunit_small_molecule_evidence; // TODO: ?
  };

  /// SML Small molecule section (mztab-m)
  struct OPENMS_DLLAPI MzTabMSmallMoleculeSectionRow
  {
    MzTabInteger identifier; ///< The small molecule’s identifier.
    MzTabStringList smf_id_refs; ///< References to all the features on which quantification has been based.
    MzTabStringList database_identifier; ///< Names of the used databases.
    MzTabStringList chemical_formula; ///< Potential chemical formula of the reported compound.
    MzTabStringList smiles; ///< Molecular structure in SMILES format.
    MzTabStringList inchi; ///< InChi of the potential compound identifications.
    MzTabStringList chemical_name; ///< Possible chemical/common names or general description
    MzTabStringList uri; ///< The source entry’s location. // TODO: URI List ?

    MzTabDoubleList theoretical_neutral_mass; ///< Precursor theoretical neutral mass
    MzTabStringList adducts; ///> Adducts
    MzTabString reliability; ///> Reliability of the given small molecule identification
    // TODO: e.g. use best search_engine score
    MzTabParameter best_id_confidence_measure; ///> The identification approach with the highest confidence
    MzTabDouble best_id_confidence_value; ///> The best confidence measure

    std::map<Size, MzTabDouble> small_molecule_abundance_assay;
    std::map<Size, MzTabDouble> small_molecule_abundance_study_variable;
    std::map<Size, MzTabDouble> small_molecule_abundance_stdev_study_variable;
    std::map<Size, MzTabDouble> small_molecule_abundance_std_error_study_variable;
    std::vector<MzTabOptionalColumnEntry> opt_; ///< Optional columns must start with “opt_”.
  };

  /// SMF Small molecule feature section (mztab-m)
  struct OPENMS_DLLAPI MzTabMSmallMoleculeFeatureSectionRow
  {
    MzTabInteger smf_identifier; ///< Within file unique identifier for the small molecule feature.
    MzTabStringList sme_id_refs; ///< Reference to the identification evidence.
    // 1=Ambiguous identification; 2=Only different evidence streams for the same molecule with no ambiguity; 3=Both ambiguous identification and multiple evidence streams.
    MzTabInteger sme_id_ref_ambiguity_code; ///< Ambiguity in identifications.
    MzTabString adduct; ///< Adduct
    MzTabParameter isotopomer; ///< //TODO? - usually used monoisotopic trace for quantification - always de-isotoped?
    MzTabDouble exp_mass_to_charge; ///< Precursor ion’s m/z.
    MzTabInteger charge; ///< Precursor ion’s charge.
    MzTabDouble retention_time; ///< Time point in seconds.
    MzTabDouble rt_start; ///< The start time of the feature on the retention time axis.
    MzTabDouble rt_end; ///< The end time of the feature on the retention time axis
    std::map<Size, MzTabDouble> small_molecule_feature_abundance_assay; // Feature abundance in every assay
    std::vector<MzTabOptionalColumnEntry> opt_; ///< Optional columns must start with “opt_”.
  };

  /// SME Small molecule evidence section (mztab-m)
  struct OPENMS_DLLAPI MzTabMSmallMoleculeEvidenceSectionRow
  {
    MzTabInteger sme_identifier; ///< Within file unique identifier for the small molecule evidence result.
    MzTabString evidence_input_id; ///< Within file unique identifier for the input data used to support this identification e.g. fragment spectrum, RT and m/z pair.
    MzTabString database_identifier; ///< The putative identification for the small molecule sourced from an external database.
    MzTabString chemical_formula; ///< The putative molecular formula.
    MzTabString smiles; ///< Potential molecular structure as SMILES.
    MzTabString inchi; ///< InChi of the potential compound identifications.
    MzTabString chemical_name; ///< Possible chemical/common names or general description
    MzTabString uri; ///< The source entry’s location.
    MzTabParameter derivatized_form; ///< //TODO?
    MzTabString adduct; ///< Adduct
    MzTabDouble exp_mass_to_charge; ///< Precursor ion’s m/z.
    MzTabInteger charge; ///< Precursor ion’s charge.
    MzTabDouble calc_mass_to_charge; ///< Precursor ion’s m/z.
    MzTabStringList spectra_ref; ///< Reference to a spectrum
    MzTabParameter identification_method; ///<
    MzTabParameter ms_level; ///<
    MzTabDouble id_confidence_measure; ///<
    MzTabInteger rank; ///< Rank of the identification (1 = best)
    std::vector<MzTabOptionalColumnEntry> opt_; ///< Optional columns must start with “opt_”.
  };

  typedef std::vector<MzTabMSmallMoleculeSectionRow> MzTabMSmallMoleculeSectionRows;
  typedef std::vector<MzTabMSmallMoleculeFeatureSectionRow> MzTabMSmallMoleculeFeatureSectionRows;
  typedef std::vector<MzTabMSmallMoleculeEvidenceSectionRow> MzTabMSmallMoleculeEvidenceSectionRows;

  /**
  @brief Data model of MzTab-M files
  Please see the offical MzTab-M specification at https://github.com/HUPO-PSI/mzTab/blob/master/specification_document-releases/2_0-Metabolomics-Release/mzTab_format_specification_2_0-M_release.adoc#use-cases-for-mztab
  */
  class OPENMS_DLLAPI MzTabM
  {
  public:
    /// Default constructor
    MzTabM();

    /// Destructor
    virtual ~MzTabM();

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

    // TODO: check if all levels (feature, evidence) can have optional columns

    /// Extract opt_ (custom, optional column names)
    std::vector<String> getMSmallMoleculeOptionalColumnNames() const;

    /// Extract opt_ (custom, optional column names)
    std::vector<String> getMSmallMoleculeFeatureOptionalColumnNames() const;

    /// Extract opt_ (custom, optional column names)
    std::vector<String> getMSmallMoleculeEvidenceOptionalColumnNames() const;

    static void addMetaInfoToOptionalColumns(const std::set<String>& keys, std::vector<MzTabOptionalColumnEntry>& opt, const String& id, const MetaInfoInterface& meta);

    // TODO: check if Modification functions are needed for Metabolomics (I guess not)
    // static std::map<Size, MzTabModificationMetaData> generateMzTabStringFromModifications(const std::vector<String>& mods);
    // static std::map<Size, MzTabModificationMetaData> generateMzTabStringFromVariableModifications(const std::vector<String>& mods);
    //static std::map<Size, MzTabModificationMetaData> generateMzTabStringFromFixedModifications(const std::vector<String>& mods);


    // TODO: see what has to be changed for metebolomics?
    // TODO: is all the info/metadata available at that point?
    static MzTab exportFeatureMapToMzTabM(const FeatureMap& feature_map, const String& filename);

    /**
      * @brief Export metabolite identifications to mzTab
      *
      * @return mzTabM object
    */
    static MzTabM exportIdentificationsToMzTabM(
        const std::vector<ProteinIdentification>& prot_ids,
        const std::vector<PeptideIdentification>& peptide_ids,
        const String& filename,
        bool first_run_inference_only,
        std::map<std::pair<size_t,size_t>,size_t>& map_run_fileidx_2_msfileidx,
        std::map<String, size_t>& idrun_2_run_index,
        bool export_empty_pep_ids = false);

    /// Generate MzTab style list of PTMs from AASequence object.
    /// All passed fixed modifications are not reported (as suggested by the standard for the PRT and PEP section).
    /// In contrast, all modifications are reported in the PSM section (see standard document for details).
    static MzTabModificationList extractModificationListFromAASequence(const AASequence& aas, const std::vector<String>& fixed_mods = std::vector<String>());

    /**
     * @brief export linked peptide features aka consensus map
     *
     * @param consensus_map		data structure of the linked peptide features
     * @param filename		input consensusXML file name
     * @param export_unidentified_features		Should not identified peptide features be exported?
     * @param export_unassigned_ids		Should unassigned identifications be exported?
     * @param export_subfeatures		The position of the consensus feature will always be exported. Should the individual subfeatures be exported as well?
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
        const String& title = "ConsensusMap export from OpenMS");


  protected:
    /// Helper function for "get...OptionalColumnNames" functions
    template <typename SectionRows>
    std::vector<String> getOptionalColumnNames_(const SectionRows& rows) const
    {
      // vector is used to preserve the column order
      std::vector<String> names;
      if (!rows.empty())
      {
        for (typename SectionRows::const_iterator it = rows.begin(); it != rows.end(); ++it)
        {
          for (std::vector<MzTabOptionalColumnEntry>::const_iterator it_opt = it->opt_.begin(); it_opt != it->opt_.end(); ++it_opt)
          {
            if (std::find(names.begin(), names.end(), it_opt->first) == names.end())
            {
              names.push_back(it_opt->first);
            }
          }
        }
      }
      return names;
    }

    static void checkSequenceUniqueness_(const std::vector<PeptideIdentification>& curr_pep_ids);

    MzTabMetaData m_meta_data_;
    MzTabMSmallMoleculeSectionRows m_small_molecule_data_;
    MzTabMSmallMoleculeFeatureSectionRows m_small_molecule_feature_data_;
    MzTabMSmallMoleculeEvidenceSectionRows m_small_molecule_evidence_data_;
    std::vector<Size> empty_rows_; ///< index of empty rows
    std::map<Size, String> comment_rows_; ///< comments

  };

} // namespace OpenMS
