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

#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/FORMAT/MzTabM.h>
#include <OpenMS/SYSTEM/File.h>

namespace OpenMS
{

  MzTabMMetaData::MzTabMMetaData()
  {
    mz_tab_version.fromCellString(String("2.0.0-M"));
  }

  const MzTabMMetaData& MzTabM::getMetaData() const
  {
    return m_meta_data_;
  }

  void MzTabM::setMetaData(const MzTabMMetaData& md)
  {
    m_meta_data_ = md;
  }

  const MzTabMSmallMoleculeSectionRows& MzTabM::getMSmallMoleculeSectionRows() const
  {
    return m_small_molecule_data_;
  }

  void MzTabM::setMSmallMoleculeSectionRows(const MzTabMSmallMoleculeSectionRows &m_smsd)
  {
    m_small_molecule_data_ = m_smsd;
  }

  const MzTabMSmallMoleculeFeatureSectionRows& MzTabM::getMSmallMoleculeFeatureSectionRows() const
  {
    return m_small_molecule_feature_data_;
  }

  void MzTabM::setMSmallMoleculeFeatureSectionRows(const MzTabMSmallMoleculeFeatureSectionRows &m_smfsd)
  {
    m_small_molecule_feature_data_ = m_smfsd;
  }

  const MzTabMSmallMoleculeEvidenceSectionRows& MzTabM::getMSmallMoleculeEvidenceSectionRows() const
  {
    return m_small_molecule_evidence_data_;
  }

  void MzTabM::setMSmallMoleculeEvidenceSectionRows(const MzTabMSmallMoleculeEvidenceSectionRows &m_smesd)
  {
    m_small_molecule_evidence_data_ = m_smesd;
  }

  std::vector<String> MzTabM::getMSmallMoleculeOptionalColumnNames() const
  {

    return getOptionalColumnNames_(m_small_molecule_data_);
  }

  std::vector<String> MzTabM::getMSmallMoleculeFeatureOptionalColumnNames() const
  {
    return getOptionalColumnNames_(m_small_molecule_feature_data_);
  }

  std::vector<String> MzTabM::getMSmallMoleculeEvidenceOptionalColumnNames() const
  {
    return getOptionalColumnNames_(m_small_molecule_evidence_data_);
  }

  void MzTabM::addMetaInfoToOptionalColumns(const std::set<String>& keys,
                                            std::vector<MzTabOptionalColumnEntry>& opt,
                                            const String& id,
                                            const MetaInfoInterface& meta)
  {
    for (String const & key : keys)
    {
      MzTabOptionalColumnEntry opt_entry;
      // column names must not contain spaces
      opt_entry.first = "opt_" + id + "_" + String(key).substitute(' ','_');
      if (meta.metaValueExists(key))
      {
        opt_entry.second = MzTabString(meta.getMetaValue(key).toString());
      } // otherwise it is default ("null")
      opt.push_back(opt_entry);
    }
  }

  // FeatureMap with associated identification data
   MzTabM MzTabM::exportFeatureMapToMzTabM(const FeatureMap& feature_map,
                                           const String& filename)
   {
     OPENMS_LOG_INFO << "exporting identification data: \"" << filename << "\" to MzTab-M: " << std::endl;

     MzTabM mztabm;
     MzTabMMetaData m_meta_data;

     // extract identification data from FeatureMap
     IdentificationData id_data = feature_map.getIdentificationData();

     UInt64 local_id = feature_map.getUniqueId();
     // mz_tab_id (mandatory)
     m_meta_data.mz_tab_id.set("local_id: " + String(local_id));

     // title (not mandatory)
     // description (not mandatory)
     // sample_processing (not mandatory)
     // instrument-name (not mandatory)
     // instrument-source (not mandatory)
     // instrument-analyzer (not mandatory)
     // instrument-detector (not mandatory)
     // meta_software.setting[0] (not mandatory)

     MzTabSoftwareMetaData meta_software;
     ControlledVocabulary cv;
     cv.loadFromOBO("PSI-MS", File::find("/CV/psi-ms.obo"));
     for (const auto& software : id_data.getProcessingSoftwares())
     {
       MzTabParameter p_software;
       ControlledVocabulary::CVTerm cvterm;
       // add TOPP - all OpenMS Tools have TOPP attached in the PSI-OBO
       std::string topp_tool = "TOPP " + software.getName();
       if (cv.hasTermWithName(topp_tool)) // asses CV-term based on tool name
       {
         cvterm = cv.getTermByName(topp_tool);
       }
       else
       {
         // use "analysis software" instead
         OPENMS_LOG_WARN << "The tool: " << topp_tool << " is currently not registered in the PSI-OBO.\n";
         OPENMS_LOG_WARN << "'The general term 'analysis software' will be used instead.\n";
         OPENMS_LOG_WARN << "Please register the tool as soon as possible in the psi-ms.obo (https://github.com/HUPO-PSI/psi-ms-CV)" << std::endl;
         cvterm = cv.getTermByName("analysis software");
       }
       p_software.fromCellString("[MS, " + cvterm.id + ", " + cvterm.name + ", " + software.getVersion() + "]");
       meta_software.software = p_software;
       m_meta_data.software[m_meta_data.software.size() + 1] = meta_software; // starts at 1
     }

     // publication (not mandatory)
     // contact name (not mandatory)
     // contact aff (not mandatory)
     // contact mail (not mandatory)
     // uri (not mandatory)
     // ext. study uri (not mandatory)
     // quantification_method (mandatory)

     MzTabParameter quantification_method;
     quantification_method.setNull(true);

     // TODO: Check processing action
     // TODO: add map using ProcessingAction and Tool!
     std::map<String, std::vector<String>> action_softwaren_ame;
     for (const auto& step : id_data.getProcessingSteps())
     {
       IdentificationDataInternal::ProcessingSoftwareRef s_ref = step.software_ref;
       std::cout << "Software: " << s_ref->getName() << std::endl;

       for (const auto& action : step.actions)
       {
         std::cout << "ProcessingAction: " << DataProcessing::NamesOfProcessingAction[action] << std::endl;
         action_softwaren_ame[action].emplace_back(s_ref->getName());
       }
     };

     // TODO: use map instead -> check for Quantification Tools
     // set quantification method based on OpenMS Tool
     // current only FeatureFinderMetabo is used
     for (const auto& software : id_data.getProcessingSoftwares())
     {
       if (software.getName() == "FeatureFinderMetabo")
       {
         ControlledVocabulary::CVTerm cvterm;
         cvterm = cv.getTermByName("LC-MS label-free quantitation analysis");
         quantification_method.fromCellString("[MS, " + cvterm.id + ", " + cvterm.name + ", ]");
       }
     }
     if (quantification_method.isNull())
     {
       OPENMS_LOG_WARN << "If the quantification of your computational analysis is not 'LC-MS label-free quantitation analysis'.\n"
                       << "Please contact a OpenMS Developer to add the appropriate tool and description to MzTab-M." << std::endl;
     }
     m_meta_data.quantification_method = quantification_method;

     // sample[1-n] (not mandatory)
     // sample[1-n]-species[1-n] (not mandatory)
     // sample[1-n]-tissue[1-n] (not mandatory)
     // sample[1-n]-cell_type[1-n] (not mandatory)
     // sample[1-n]-disease[1-n] (not mandatory)
     // sample[1-n]-description (not mandatory)
     // sample[1-n]-custom[1-n] (not mandatory)

     MzTabMMSRunMetaData meta_ms_run;
     String input_file_name;
     auto input_files = id_data.getInputFiles();
     for (const auto& input_file : input_files) // should only be one in featureXML
     {
       input_file_name = input_file.name;
       meta_ms_run.location.set(input_file_name);
     }
     // meta_ms_run.location.set(input_files[0].name);

     // ms_run[1-n]-instrument_ref (not mandatory)
     // ms_run[1-n]-format (not mandatory)
     // ms_run[1-n]-id_format (not mandatory)
     // ms_run[1-n]-fragmentation_method[1-n] (not mandatory)

     // ms_run[1-n]-scan_polarity[1-n] (mandatory)
     // assess scan polarity based on the adducts
     std::vector<MzTabParameter> polartiy;
     auto adducts = id_data.getAdducts();
     if (!adducts.empty())
     {
       for (const auto& adduct : adducts)
       {
         String adduct_suffix = adduct.getName().suffix(';').trim();
         if (adduct_suffix == "+")
         {
           ControlledVocabulary::CVTerm cvterm;
           cvterm = cv.getTermByName("positive scan");
           MzTabParameter spol;
           spol.fromCellString("[MS, " + cvterm.id + ", " + cvterm.name + ", ]");
           polartiy.emplace_back(spol);
         }
         else
         {
           ControlledVocabulary::CVTerm cvterm;
           cvterm = cv.getTermByName("negative scan");
           MzTabParameter spol;
           spol.fromCellString("[MS, " + cvterm.id + ", " + cvterm.name + ", ]");
           polartiy.emplace_back(spol);
         }
       }
     }
     else
     {
         // if no adduct information is available warn, but assume positive mode.
         OPENMS_LOG_WARN << "No adduct information available: scan polarity will be assumed to be positive." << std::endl;
         ControlledVocabulary::CVTerm cvterm;
         cvterm = cv.getTermByName("positive scan");
         MzTabParameter spol;
         spol.fromCellString("[MS, " + cvterm.id + ", " + cvterm.name + ", ]");
         polartiy.emplace_back(spol);
     }
     meta_ms_run.scan_polarity.set(polartiy);

     // ms_run[1-n]-hash (not mandatory)
     // ms_run[1-n]-hash_method (not mandatory)

     MzTabMAssayMetaData meta_ms_assay;
     // assay[1-n] (mandatory)
     meta_ms_assay.name = MzTabString("assay_" + File::basename(input_file_name).prefix('.').trim());
     // assay[1-n]-custom[1-n] (not mandatory)
     // assay[1-n]-external_uri (not mandatory)
     // assay[1-n]-sample_ref (not mandatory)

     // assay[1-n]-ms_run_ref (mandatory)
     std::vector<int> ms_run_ref(1);
     meta_ms_assay.ms_run_ref = ms_run_ref;

     MzTabMStudyVariableMetaData meta_ms_study_variable;
     // study_variable[1-n] (mandatory)
     meta_ms_study_variable.name = MzTabString("study_variable_" + File::basename(input_file_name).prefix('.').trim());

     // study_variable[1-n]-assay_refs (mandatory)
     std::vector<int> assay_refs(1);
     meta_ms_study_variable.assay_refs = assay_refs;

     // study_variable[1-n]-average_function (not mandatory)
     // study_variable[1-n]-variation_function (not mandatory)

     // study_variable[1-n]-description (mandatory)
     meta_ms_study_variable.description = MzTabString("study_variable_" + File::basename(input_file_name).prefix('.').trim());

     // study_variable[1-n]-factors (not mandatory)
     // custom[1-n] (not mandatory)

     MzTabCVMetaData meta_cv;
     // cv[1-n]-label (mandatory)
     meta_cv.label = MzTabString(cv.name());
     // cv[1-n]-full_name (mandatory)
     meta_cv.full_name = MzTabString(cv.label());
     // cv[1-n]-version (mandatory)
     meta_cv.version = MzTabString(cv.version());
     // cv[1-n]-uri (mandatory)
     meta_cv.url = MzTabString(cv.url());

     // these have to be added to the identification data
     // in the actual tool writes the mztam-m
     MzTabMDatabaseMetaData meta_db;
     for (const auto& db : id_data.getDBSearchParams())
     {
       if (db.database == "") // no database
       {
         meta_db.prefix = MzTabString("null");
         meta_db.version = MzTabString("Unknown");
         meta_db.database.fromCellString("[,, no database , null]");
       }
       else if (db.database.find("custom") != std::string::npos) // custom database
       {
         meta_db.prefix = MzTabString("null");
         meta_db.version = MzTabString(db.database_version);
         meta_db.database.fromCellString("[,, " + db.database + ", ]");
       }
       else // assumption that prefix is the same as database name
       {
         meta_db.prefix = MzTabString(db.database);
         meta_db.version = MzTabString(db.database_version);
         meta_db.database.fromCellString("[,," + db.database + ", ]");
       }
       meta_db.uri = MzTabString("null"); // URL is not available at this point
       m_meta_data.database[m_meta_data.database.size() + 1] = meta_db; // starts at 1
     }

     // derivatization_agent[1-n] (not mandatory)

     // TODO: use quantification tools map!
     // small_molecule-quantification_unit (mandatory)
     // small_molecule_feature-quantification_unit (mandatory)
     MzTabParameter quantification_unit;
     quantification_unit.setNull(true);
     for (const auto& software : id_data.getProcessingSoftwares())
     {
       if (software.getName() == "FeatureFinderMetabo")
       {
         if (software.metaValueExists("parameter: algorithm:mtd:quant_method"))
         {
           String quant_method = software.getMetaValue("parameter: algorithm:mtd:quant_method");
           if (quant_method == "area")
           {
             ControlledVocabulary::CVTerm cvterm;
             cvterm = cv.getTermByName("MS1 feature area");
             quantification_unit.fromCellString("[MS, " + cvterm.id + ", " + cvterm.name + ", ]");
           }
           else if (quant_method == "median")
           {
             ControlledVocabulary::CVTerm cvterm;
             cvterm = cv.getTermByName("median");
             quantification_unit.fromCellString("[MS, " + cvterm.id + ", " + cvterm.name + ", ]");
           }
           else // max_height
           {
             ControlledVocabulary::CVTerm cvterm;
             cvterm = cv.getTermByName("MS1 feature maximum intensity");
             quantification_unit.fromCellString("[MS, " + cvterm.id + ", " + cvterm.name + ", ]");
           }
         }
       }
     }
     if (quantification_unit.isNull())
     {
       OPENMS_LOG_WARN << "It was not possible to assess the quantification_unit - MS1 feature area - will be used as default." << std::endl;
       ControlledVocabulary::CVTerm cvterm;
       cvterm = cv.getTermByName("MS1 feature area");
       quantification_unit.fromCellString("[MS, " + cvterm.id + ", " + cvterm.name + ", ]");
     }
     m_meta_data.small_molecule_quantification_unit = quantification_unit;
     m_meta_data.small_molecule_feature_quantification_unit = quantification_unit;

     // small_molecule-identification_reliability (not mandatory)

     // TODO: How to use actual registered CV-terms?
     int score_counter = 1;
     for (const auto& software : id_data.getProcessingSoftwares())
     {
       std::vector<IdentificationDataInternal::ScoreTypeRef> score_type_refs = software.assigned_scores;
       for (const auto& score_type_ref : score_type_refs)
       {
         m_meta_data.id_confidence_measure[score_counter].fromCellString("[,, " + score_type_ref->cv_term.getName()+ ", ]");
         score_counter += 1;
       }
     }
     // colunit-small_molecule (not mandatory)
     // colunit-small_molecule_feature (not mandatory)
     // colunit-small_molecule_evidence (not mandatory)

     m_meta_data.ms_run[1] = meta_ms_run;
     m_meta_data.assay[1] = meta_ms_assay;
     m_meta_data.study_variable[1] = meta_ms_study_variable;

     mztabm.setMetaData(m_meta_data);

     // iterate over features and construct the feature, summary and evidence section
     MzTabMSmallMoleculeSectionRows smss;
     MzTabMSmallMoleculeFeatureSectionRows smfs;
     MzTabMSmallMoleculeEvidenceSectionRows smes;

     int summary_section_counter = 1;
     int feature_section_entry_counter = 1;
     int evidence_section_entry_counter = 1;

     // iterate over features and fill all sections
     for (auto& f : feature_map)
     {
       // iterate over the identification of the ObservationMatches
       auto match_refs = f.getIDMatches();

       // evidence section:
       // MzTabInteger sme_identifier; ///< Within file unique identifier for the small molecule evidence result.
       // MzTabString evidence_input_id; ///< Within file unique identifier for the input data used to support this identification e.g. fragment spectrum, RT and m/z pair.
       // MzTabString database_identifier; ///< The putative identification for the small molecule sourced from an external database.
       // MzTabString chemical_formula; ///< The putative molecular formula.
       // MzTabString smiles; ///< Potential molecular structure as SMILES.
       // MzTabString inchi; ///< InChi of the potential compound identifications.
       // MzTabString chemical_name; ///< Possible chemical/common names or general description
       // MzTabString uri; ///< The source entry’s location.
       // MzTabParameter derivatized_form; ///<
       // MzTabString adduct; ///< Adduct //<
       // MzTabDouble exp_mass_to_charge; ///< Precursor ion’s m/z.
       // MzTabInteger charge; ///< Precursor ion’s charge.
       // MzTabDouble calc_mass_to_charge; ///< Precursor ion’s m/z.
       // MzTabStringList spectra_ref; ///< Reference to a spectrum
       // MzTabParameter identification_method; ///< Database search, search engine or process that was used to identify this small molecule
       // MzTabParameter ms_level; ///< The highest MS level used to inform identification
       // MzTabDouble id_confidence_measure; ///< Statistical value or score for the identification
       // MzTabInteger rank; ///< Rank of the identification (1 = best)

       // feature row based on number of individual adducts!
       std::map<String, std::vector<int>> evidence_id_ref_per_adduct;
       for (const IdentificationDataInternal::ObservationMatchRef& match_ref : match_refs) // iterate over all identifications of a feature
       {
         // evidence section
         MzTabMSmallMoleculeEvidenceSectionRow sme;

         // IdentifiedCompound
         IdentificationData::IdentifiedMolecule molecule = match_ref->identified_molecule_var;
         IdentificationData::IdentifiedCompoundRef compound_ref = molecule.getIdentifiedCompoundRef();

         sme.sme_identifier = MzTabInteger(evidence_section_entry_counter);
         sme.evidence_input_id = MzTabString("mass=" + String(f.getMZ()) + ",rt=" + String(f.getRT()));
         sme.database_identifier = MzTabString(compound_ref->identifier);
         sme.chemical_formula = MzTabString(compound_ref->formula.toString());
         sme.smiles = MzTabString(compound_ref->smile);
         sme.inchi = MzTabString(compound_ref->inchi);
         sme.chemical_name = MzTabString(compound_ref->name);
         sme.uri.setNull(true);
         sme.derivatized_form.setNull(true);
         String adduct = getAdductString_(match_ref);
         sme.adduct = MzTabString(adduct);
         sme.exp_mass_to_charge = MzTabDouble(f.getMZ());
         sme.charge = MzTabInteger(f.getCharge());
         sme.calc_mass_to_charge = MzTabDouble(compound_ref->formula.getMonoWeight());
         // TODO: ISSUE: IdentificationData only one ref per identifiedmolecule?
         // TODO: ISSUE: What about e.g. SIRIUS using multiple MS2 spectra for one identification?
         sme.spectra_ref.fromCellString(match_ref->observation_ref->data_id); // MzTabStringList
         // based on tool used for identification (CV-Term)
         // TODO: Would actually have to be set per ID
         // TODO: Since there is a possibility to more that on identification
         // TODO: Get that from the tools used?
         // TODO: How to add that information here?
         // TODO: Should that be added upstream of the pipeline?
         // use processing action identification
         // to get the correct tool? And then use the Tool Name?
         // MzTabParameter identification_method;
         // TODO: not sure what to do here!
         // TODO: depends on the methods seen above!
         // MzTabParameter ms_level; ///< The highest MS level used to inform identification
         // TODO: Score of the identification?!
         // MzTabDouble id_confidence_measure; ///< Statistical value or score for the identification
         sme.rank = MzTabInteger(1); // defaults to 1 if no rank system is used.
         // TODO: How to add opt_ columns

         evidence_id_ref_per_adduct[adduct].emplace_back(evidence_section_entry_counter);
         evidence_section_entry_counter += 1;
         smes.emplace_back(sme);
       }

       // feature section

       // MzTabInteger smf_identifier; ///< Within file unique identifier for the small molecule feature.
       // MzTabStringList sme_id_refs; ///< Reference to the identification evidence.
       // MzTabInteger sme_id_ref_ambiguity_code; ///< Ambiguity in identifications.
       // MzTabString adduct; ///< Adduct
       // MzTabParameter isotopomer; ///< If de-isotoping has not been performed, then the isotopomer quantified MUST be reported here.
       // MzTabDouble exp_mass_to_charge; ///< Precursor ion’s m/z.
       // MzTabInteger charge; ///< Precursor ion’s charge.
       // MzTabDouble retention_time; ///< Time point in seconds.
       // MzTabDouble rt_start; ///< The start time of the feature on the retention time axis.
       // MzTabDouble rt_end; ///< The end time of the feature on the retention time axis
       // std::map<Size, MzTabDouble> small_molecule_feature_abundance_assay; ///< Feature abundance in every assay
       // std::vector<MzTabOptionalColumnEntry> opt_; ///< Optional columns must start with “opt_”

       // one feature entry per adduct - iterate evidences_per_adduct
       for (const auto& epa : evidence_id_ref_per_adduct)
       {
         MzTabMSmallMoleculeFeatureSectionRow smf;
         smf.smf_identifier = MzTabInteger(feature_section_entry_counter);
         std::vector<MzTabString> corresponding_evidences;
         for (const auto& evidence : epa.second)
         {
           corresponding_evidences.emplace_back(String(evidence));
         }
         smf.sme_id_refs.set(corresponding_evidences);
         smf.adduct = MzTabString(epa.first);
         if (epa.second.size() <= 1)
         {
           smf.sme_id_ref_ambiguity_code.setNull(true);
         }
         else
         {
           smf.sme_id_ref_ambiguity_code = MzTabInteger(1);
         }
         smf.isotopomer.setNull(true);
         smf.exp_mass_to_charge = MzTabDouble(f.getMZ());
         smf.charge = MzTabInteger(f.getCharge());
         smf.retention_time = MzTabDouble(f.getRT());
         smf.rt_start.setNull(true); // TODO: how to get that information in the future
         smf.rt_end.setNull(true); // TODO: haw to get that information in the future
         smf.small_molecule_feature_abundance_assay[1] = MzTabDouble(f.getIntensity()); // only one map in featureXML
         // TODO: how do we add opt_ columns check MzTab

         smfs.emplace_back(smf);
       }

       // small mol section (summary)
       MzTabMSmallMoleculeSectionRows sms;

        //  MzTabInteger identifier; ///< The small molecule’s identifier.
        //  MzTabStringList smf_id_refs; ///< References to all the features on which quantification has been based.
        //  MzTabStringList database_identifier; ///< Names of the used databases.
        //  MzTabStringList chemical_formula; ///< Potential chemical formula of the reported compound.
        //  MzTabStringList smiles; ///< Molecular structure in SMILES format.
        //  MzTabStringList inchi; ///< InChi of the potential compound identifications.
        //  MzTabStringList chemical_name; ///< Possible chemical/common names or general description
        //  MzTabStringList uri; ///< The source entry’s location. // TODO: URI List ?
        //  MzTabDoubleList theoretical_neutral_mass; ///< Precursor theoretical neutral mass
        //  MzTabStringList adducts; ///< Adducts
        //  // TODO: https://github.com/HUPO-PSI/mzTab/blob/master/specification_document-releases/2_0-Metabolomics-Release/mzTab_format_specification_2_0-M_release.adoc#6311-reliability
        //  MzTabString reliability; ///< Reliability of the given small molecule identification
        //  // TODO: e.g. use best search_engine score
        //  MzTabParameter best_id_confidence_measure; ///< The identification approach with the highest confidence
        //  MzTabDouble best_id_confidence_value; ///< The best confidence measure
        //  std::map<Size, MzTabDouble> small_molecule_abundance_assay; ///<
        //  std::map<Size, MzTabDouble> small_molecule_abundance_study_variable; ///<
        //  std::map<Size, MzTabDouble> small_molecule_abundance_stdev_study_variable; ///<
        //  std::map<Size, MzTabDouble> small_molecule_abundance_std_error_study_variable; ///<

     }
     mztabm.setMSmallMoleculeFeatureSectionRows(smfs);
     mztabm.setMSmallMoleculeEvidenceSectionRows(smes);
     // mztabm.setMSmallMoleculeSectionRows();
     return mztabm;
   }

  String MzTabM::getAdductString_(const IdentificationDataInternal::ObservationMatchRef& match_ref)
  {
    String adduct_name;
    if (match_ref->adduct_opt)
    {
      adduct_name = (*match_ref->adduct_opt)->getName();
      // M+H;1+ -> [M+H]1+
      if (adduct_name.find(';') != std::string::npos) // wrong format -> reformat
      {
        String prefix = adduct_name.substr(0, adduct_name.find(';'));
        String suffix = adduct_name.substr(adduct_name.find(';') + 1, adduct_name.size());
        adduct_name = "[" + prefix + "]" + suffix;
      }
    }
    else
    {
      adduct_name = "null";
    }
    return adduct_name;
  }
}