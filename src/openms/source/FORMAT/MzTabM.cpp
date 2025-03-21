// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/FORMAT/MzTabM.h>
#include <OpenMS/SYSTEM/File.h>
#include <regex>

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

  void MzTabM::setMSmallMoleculeSectionRows(const MzTabMSmallMoleculeSectionRows &m_smlsd)
  {
    m_small_molecule_data_ = m_smlsd;
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

  void MzTabM::getFeatureMapMetaValues_(const FeatureMap& feature_map,
                                        std::set<String>& feature_user_value_keys,
                                        std::set<String>& observationmatch_user_value_keys,
                                        std::set<String>& compound_user_value_keys)
  {
    for (Size i = 0; i < feature_map.size(); ++i)
    {
      // feature section optional columns
      const Feature& f = feature_map[i];
      std::vector<String> keys;
      f.getKeys(keys);
      // replace whitespaces with underscore
      std::transform(keys.begin(), keys.end(), keys.begin(), [&](String& s) { return s.substitute(' ', '_'); });
      feature_user_value_keys.insert(keys.begin(), keys.end());

      auto match_refs = f.getIDMatches();
      for (const IdentificationDataInternal::ObservationMatchRef& match_ref : match_refs)
      {
        // feature section optional columns
        std::vector<String> obsm_keys;
        match_ref->getKeys(obsm_keys);
        // replace whitespaces with underscore
        std::transform(obsm_keys.begin(), obsm_keys.end(), obsm_keys.begin(), [&](String& s) { return s.substitute(' ', '_'); });

        // remove "IDConverter_trace" metadata from the ObservationMatch
        // introduced by the IdentificationDataConverter
        // since it leads to convolution of IDConverter_trace_* optional columns
        for (const auto& key : obsm_keys)
        {
          if (!key.hasSubstring("IDConverter_trace"))
          {
            observationmatch_user_value_keys.insert(key);
          }
        }

        // evidence section optional columns
        IdentificationData::IdentifiedMolecule molecule = match_ref->identified_molecule_var;
        IdentificationData::IdentifiedCompoundRef compound_ref = molecule.getIdentifiedCompoundRef();
        std::vector<String> compound_keys;
        compound_ref->getKeys(compound_keys);
        // replace whitespaces with underscore
        std::transform(compound_keys.begin(), compound_keys.end(), compound_keys.begin(), [&](String& s) { return s.substitute(' ', '_'); });
        compound_user_value_keys.insert(compound_keys.begin(), compound_keys.end());
      }
    }
  }

  // FeatureMap with associated identification data
  MzTabM MzTabM::exportFeatureMapToMzTabM(const FeatureMap& feature_map)
  {
    MzTabM mztabm;
    MzTabMMetaData m_meta_data;

    // extract identification data from FeatureMap
    const IdentificationData& id_data = feature_map.getIdentificationData();

    OPENMS_PRECONDITION(!id_data.empty(),
                        "The FeatureMap has to have a non empty IdentificationData object attached!")

    // extract MetaValues from FeatureMap
    std::set<String> feature_user_value_keys;
    std::set<String> observationmatch_user_value_keys;
    std::set<String> compound_user_value_keys;
    MzTabM::getFeatureMapMetaValues_(feature_map, feature_user_value_keys, observationmatch_user_value_keys, compound_user_value_keys);

    // ####################################################
    // MzTabMetaData
    // ####################################################

    std::regex reg_backslash{R"(\\)"};
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
    MzTabString reliability = MzTabString("2"); // initialize at 2 (should be valid for all tools - putatively annotated compound)
    cv.loadFromOBO("PSI-MS", File::find("/CV/psi-ms.obo"));
    for (const auto& software : id_data.getProcessingSoftwares())
    {
      if (software.metaValueExists("reliability"))
      {
        reliability = MzTabString(std::string(software.getMetaValue("reliability")));
      }
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
    std::map<String, std::vector<String>> action_software_name;
    for (const auto& step : id_data.getProcessingSteps())
    {
      IdentificationDataInternal::ProcessingSoftwareRef s_ref = step.software_ref;
      for (const auto& action : step.actions)
      {
        action_software_name[action].emplace_back(s_ref->getName());
      }
    };

    // set quantification method based on OpenMS Tool(s)
    // current only FeatureFinderMetabo is used
    for (const auto& quantification_software : action_software_name[DataProcessing::QUANTITATION])
    {
      if (quantification_software == "FeatureFinderMetabo")
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
    std::string input_file_name;
    auto input_files = id_data.getInputFiles();
    for (const auto& input_file : input_files) // should only be one in featureXML
    {
      input_file_name = input_file.name;
      input_file_name = String(std::regex_replace(input_file_name, reg_backslash, "/"));
      if (!String(input_file_name).hasPrefix("file://")) input_file_name = "file://" + input_file_name;
      meta_ms_run.location.set(input_file_name);
    }
    // meta_ms_run.location.set(input_files[0].name);

    // ms_run[1-n]-instrument_ref (not mandatory)
    // ms_run[1-n]-format (not mandatory)
    // ms_run[1-n]-id_format (not mandatory)
    // ms_run[1-n]-fragmentation_method[1-n] (not mandatory)

    // ms_run[1-n]-scan_polarity[1-n] (mandatory)
    // assess scan polarity based on the first adduct
    auto adducts = id_data.getAdducts();
    if (!adducts.empty())
    {
      std::string_view first_adduct;
      for (const auto& adduct : adducts)
      {
        first_adduct = adduct.getName();
        break;
      }
      if (first_adduct.at(first_adduct.size() - 1) == '+')
      {
        ControlledVocabulary::CVTerm cvterm;
        cvterm = cv.getTermByName("positive scan");
        MzTabParameter spol;
        spol.fromCellString("[MS, " + cvterm.id + ", " + cvterm.name + ", ]");
        meta_ms_run.scan_polarity[1] = spol;
      }
      else
      {
        ControlledVocabulary::CVTerm cvterm;
        cvterm = cv.getTermByName("negative scan");
        MzTabParameter spol;
        spol.fromCellString("[MS, " + cvterm.id + ", " + cvterm.name + ", ]");
        meta_ms_run.scan_polarity[1] = spol;
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
      meta_ms_run.scan_polarity[1] = spol;
    }

    // ms_run[1-n]-hash (not mandatory)
    // ms_run[1-n]-hash_method (not mandatory)

    MzTabMAssayMetaData meta_ms_assay;
    // assay[1-n] (mandatory)
    meta_ms_assay.name = MzTabString("assay_" + File::basename(input_file_name).prefix('.').trim());
    // assay[1-n]-custom[1-n] (not mandatory)
    // assay[1-n]-external_uri (not mandatory)
    // assay[1-n]-sample_ref (not mandatory)

    // assay[1-n]-ms_run_ref (mandatory)
    MzTabInteger ms_run_ref(1);
    meta_ms_assay.ms_run_ref = ms_run_ref;

    MzTabMStudyVariableMetaData meta_ms_study_variable;
    // study_variable[1-n] (mandatory)
    meta_ms_study_variable.name = MzTabString("study_variable_" + File::basename(input_file_name).prefix('.').trim());

    // study_variable[1-n]-assay_refs (mandatory)
    std::vector<int> assay_refs;
    assay_refs.emplace_back(1);
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
    m_meta_data.cv[1] = meta_cv;

    // these have to be added to the identification data
    // in the actual tool writes the mztam-m
    MzTabMDatabaseMetaData meta_db;
    for (const auto& db : id_data.getDBSearchParams())
    {
      if (db.database == "") // no database
      {
        meta_db.prefix.setNull(true);
        meta_db.version = MzTabString("Unknown");
        meta_db.database.fromCellString("[,, no database , null]");
      }
      else if (db.database.find("custom") != std::string::npos) // custom database
      {
        meta_db.prefix.setNull(true);
        meta_db.version = MzTabString(db.database_version);
        meta_db.database.fromCellString("[,, " + db.database + ", ]");
      }
      else // assumption that prefix is the same as database name
      {
        meta_db.prefix = MzTabString(db.database);
        meta_db.version = MzTabString(db.database_version);
        meta_db.database.fromCellString("[,," + db.database + ", ]");
      }
      if (db.metaValueExists("database_location"))
      {
        std::vector<std::string> db_loc = ListUtils::create<std::string>(db.getMetaValue("database_location"), '|');
        for (auto& loc : db_loc)
        {
          loc = String(std::regex_replace(loc, reg_backslash, "/"));
          if (!String(loc).hasPrefix("file://")) loc = "file://" + loc;
        }
        String db_location_uri = ListUtils::concatenate(db_loc, '|');
        meta_db.uri = MzTabString(db_location_uri);
      }
      else
      {
        meta_db.uri.setNull(true);
      }
      m_meta_data.database[m_meta_data.database.size() + 1] = meta_db; // starts at 1
    }

    // derivatization_agent[1-n] (not mandatory)
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
    // small_molecule-quantification_unit (mandatory)
    // small_molecule_feature-quantification_unit (mandatory)
    m_meta_data.small_molecule_quantification_unit = quantification_unit;
    m_meta_data.small_molecule_feature_quantification_unit = quantification_unit;

    // small_molecule-identification_reliability (mandatory)
    MzTabParameter rel;
    ControlledVocabulary::CVTerm cvterm_rel = cv.getTermByName("compound identification confidence level");
    rel.fromCellString("[MS, " + cvterm_rel.id + ", " + cvterm_rel.name + ", ]");
    m_meta_data.small_molecule_identification_reliability = rel;

    int software_score_counter = 0;
    std::vector<String> identification_tools = action_software_name[DataProcessing::IDENTIFICATION];
    std::vector<IdentificationDataInternal::ScoreTypeRef> id_score_refs;
    for (const IdentificationDataInternal::ProcessingSoftware& software : id_data.getProcessingSoftwares())
    {
      // check if in "Identification Vector"
      if (std::find(identification_tools.begin(), identification_tools.end(), software.getName()) != identification_tools.end())
      {
        for (const IdentificationDataInternal::ScoreTypeRef& score_type_ref : software.assigned_scores)
        {
          ++software_score_counter;
          m_meta_data.id_confidence_measure[software_score_counter].fromCellString("[,, " + score_type_ref->cv_term.getName() + ", ]");
          id_score_refs.emplace_back(score_type_ref);
        }
      }
    }
    // colunit-small_molecule (not mandatory)
    // colunit-small_molecule_feature (not mandatory)
    // colunit-small_molecule_evidence (not mandatory)

    m_meta_data.ms_run[1] = meta_ms_run;
    m_meta_data.assay[1] = meta_ms_assay;
    m_meta_data.study_variable[1] = meta_ms_study_variable;

    // iterate over features and construct the feature, summary and evidence section
    MzTabMSmallMoleculeSectionRows smls;
    MzTabMSmallMoleculeFeatureSectionRows smfs;
    MzTabMSmallMoleculeEvidenceSectionRows smes;

    // set identification method based on OpenMS Tool(s)
    // usually only one identification_method in one featureXML
    MzTabParameter identification_method;
    identification_method.setNull(true);
    MzTabParameter ms_level;
    ms_level.setNull(true);
    for (const auto& identification_software : action_software_name[DataProcessing::IDENTIFICATION])
    {
      int id_mslevel = 0;
      if (identification_software == "AccurateMassSearch")
      {
        ControlledVocabulary::CVTerm cvterm;
        cvterm = cv.getTermByName("accurate mass");
        identification_method.fromCellString("[MS, " + cvterm.id + ", " + cvterm.name + ", ]");
        id_mslevel = 1;
      }
      if (identification_software == "SiriusAdapter")
      {
        ControlledVocabulary::CVTerm cvterm;
        cvterm = cv.getTermByName("de novo search");
        identification_method.fromCellString("[MS, " + cvterm.id + ", " + cvterm.name + ", ]");
        id_mslevel = 2;
      }
      if (identification_software == "MetaboliteSpectralMatcher")
      {
        ControlledVocabulary::CVTerm cvterm;
        cvterm = cv.getTermByName("TOPP SpecLibSearcher");
        identification_method.fromCellString("[MS, " + cvterm.id + ", " + cvterm.name + ", ]");
        id_mslevel = 2;
      }
      ControlledVocabulary::CVTerm cvterm_level = cv.getTermByName("ms level");
      ms_level.fromCellString("[MS, " + cvterm_level.id + ", " + cvterm_level.name + ", " + String(id_mslevel) + "]");
    }
    if (identification_method.isNull())
    {
      OPENMS_LOG_WARN << "The identification method of your computational analysis can not be assessed'.\n"
                      << "Please check if the ProcessingActions are set correctly!" << std::endl;
    }

    // ####################################################
    //
    // ####################################################

    int feature_section_entry_counter = 1;
    int evidence_section_entry_counter = 1;
    for (auto& f : feature_map) // iterate over features and fill all sections
    {
      auto match_refs = f.getIDMatches();
      if (match_refs.empty()) // features without identification
      {
        MzTabMSmallMoleculeFeatureSectionRow smf;
        smf.smf_identifier = MzTabString(feature_section_entry_counter);
        std::vector<MzTabString> corresponding_evidences;
        smf.sme_id_refs.setNull(true);
        if (f.metaValueExists("adducts"))
        {
          StringList adducts = f.getMetaValue("adducts");
          smf.adduct = MzTabString(ListUtils::concatenate(adducts,'|'));
        }
        else
        {
          smf.adduct.setNull(true);
        }
        smf.sme_id_ref_ambiguity_code.setNull(true);
        smf.isotopomer.setNull(true);
        smf.exp_mass_to_charge = MzTabDouble(f.getMZ());
        smf.charge = MzTabInteger(f.getCharge());
        smf.retention_time = MzTabDouble(f.getRT());
        smf.rt_start.setNull(true);
        smf.rt_end.setNull(true);
        smf.small_molecule_feature_abundance_assay[1] = MzTabDouble(f.getIntensity()); // only one map in featureXML

        addMetaInfoToOptionalColumns(feature_user_value_keys, smf.opt_, String("global"), f);

        smfs.emplace_back(smf);
        ++feature_section_entry_counter;
      }
      else
      {
        // feature row based on number of individual identifications and adducts!
        std::map<String, std::vector<int>> evidence_id_ref_per_adduct;

        // TODO: Remove copy operation (operator< IDData Ref)
        std::set<IdentificationDataInternal::ObservationMatchRef, CompareMzTabMMatchRef> sorted_match_refs(match_refs.begin(), match_refs.end());

        for (const auto& ref : sorted_match_refs) // iterate over all identifications of a feature
        {
          // evidence section
          MzTabMSmallMoleculeEvidenceSectionRow sme;

          // IdentifiedCompound
          IdentificationData::IdentifiedMolecule molecule = ref->identified_molecule_var;
          IdentificationData::IdentifiedCompoundRef compound_ref = molecule.getIdentifiedCompoundRef();

          sme.sme_identifier = MzTabString(evidence_section_entry_counter);
          sme.evidence_input_id = MzTabString("mass=" + String(f.getMZ()) + ",rt=" + String(f.getRT()));
          sme.database_identifier = MzTabString(compound_ref->identifier);
          sme.chemical_formula = MzTabString(compound_ref->formula.toString());
          sme.smiles = MzTabString(compound_ref->smile);
          sme.inchi = MzTabString(compound_ref->inchi);
          sme.chemical_name = MzTabString(compound_ref->name);
          sme.uri.setNull(true);
          sme.derivatized_form.setNull(true);
          String adduct = getAdductString_(ref);
          sme.adduct = MzTabString(adduct);
          sme.exp_mass_to_charge = MzTabDouble(f.getMZ());
          sme.charge = MzTabInteger(f.getCharge());
          sme.calc_mass_to_charge = MzTabDouble(compound_ref->formula.getMonoWeight());
          // For e.g. SIRIUS using multiple MS2 spectra for one identification
          // use the with pipe concatenated native_ids as spectra ref
          // this should  also be available match_ref
          MzTabSpectraRef sp_ref;
          sp_ref.setMSFile(1);
          sp_ref.setSpecRef(ref->observation_ref->data_id);
          sme.spectra_ref = sp_ref;
          sme.identification_method = identification_method; // based on tool used for identification (CV-Term)
          sme.ms_level = ms_level;
          int score_counter = 0;
          for (const auto& id_score_ref : id_score_refs) // vector of references based on the ProcessingStep
          {
            ++score_counter; //starts at 1 anyway
            sme.id_confidence_measure[score_counter] = MzTabDouble(ref->getScore(id_score_ref).first);
          }
          sme.rank = MzTabInteger(1); // defaults to 1 if no rank system is used

          addMetaInfoToOptionalColumns(observationmatch_user_value_keys, sme.opt_, String("global"), *ref);
          addMetaInfoToOptionalColumns(compound_user_value_keys, sme.opt_, String("global"), *compound_ref);

          evidence_id_ref_per_adduct[adduct].emplace_back(evidence_section_entry_counter);
          evidence_section_entry_counter += 1;
          smes.emplace_back(sme);
        }

        // feature section
        // one feature entry per adduct - iterate evidences_per_adduct
        for (const auto& epa : evidence_id_ref_per_adduct)
        {
          MzTabMSmallMoleculeFeatureSectionRow smf;
          smf.smf_identifier = MzTabString(feature_section_entry_counter);
          std::vector<MzTabString> corresponding_evidences;
          for (const auto& evidence : epa.second)
          {
            corresponding_evidences.emplace_back(evidence);
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
          smf.rt_start.setNull(true);
          smf.rt_end.setNull(true);
          smf.small_molecule_feature_abundance_assay[1] = MzTabDouble(f.getIntensity()); // only one map in featureXML

          addMetaInfoToOptionalColumns(feature_user_value_keys, smf.opt_, String("global"), f);

          smfs.emplace_back(smf);
          ++feature_section_entry_counter;
        }
      }
    }

    // based summary on available features and evidences
    // OpenMS does currently not aggregate the information of two features with corresponding adducts
    // e.g. F1 mz = 181.0712 rt = 10, [M+H]1+; F2 mz = 203.0532, rt = 10 [M+Na]1+ (neutral mass for both: 180.0634)
    // features will be represented individually here.

    for (const auto& smf : smfs)
    {
      MzTabMSmallMoleculeSectionRow sml;

      sml.sml_identifier = smf.smf_identifier;
      sml.smf_id_refs.set({smf.smf_identifier});
      std::vector<MzTabString> database_identifier;
      std::vector<MzTabString> chemical_formula;
      std::vector<MzTabString> smiles;
      std::vector<MzTabString> inchi;
      std::vector<MzTabString> chemical_name;
      std::vector<MzTabString> uri;
      std::vector<MzTabDouble> theoretical_neutral_mass;
      std::vector<MzTabString> adducts;
      for (const MzTabString & evidence : smf.sme_id_refs.get())
      {
        const auto& current_row_it = std::find_if(smes.begin(), smes.end(), [&evidence] (const MzTabMSmallMoleculeEvidenceSectionRow& sme) { return sme.sme_identifier.get() == evidence.get(); });
        database_identifier.emplace_back(current_row_it->database_identifier);
        chemical_formula.emplace_back(current_row_it->chemical_formula);
        smiles.emplace_back(current_row_it->smiles);
        inchi.emplace_back(current_row_it->inchi);
        chemical_name.emplace_back(current_row_it->chemical_name);
        uri.emplace_back(current_row_it->uri);
        MzTabString cm = current_row_it->chemical_formula;
        if (cm.toCellString() != "" && cm.toCellString() != "null" )
        {
          theoretical_neutral_mass.emplace_back(EmpiricalFormula(cm.toCellString()).getMonoWeight());
        }
        else
        {
          MzTabDouble dnull;
          dnull.setNull(true);
          theoretical_neutral_mass.emplace_back(dnull);
        }
        adducts.emplace_back(current_row_it->adduct);
      }
      sml.database_identifier.set(database_identifier);
      sml.chemical_formula.set(chemical_formula);
      sml.smiles.set(smiles);
      sml.inchi.set(inchi);
      sml.chemical_name.set(chemical_name);
      sml.uri.set(uri);
      sml.theoretical_neutral_mass.set(theoretical_neutral_mass);
      sml.adducts.set(adducts);
      sml.reliability = reliability;
      sml.best_id_confidence_measure.setNull(true);
      sml.best_id_confidence_value.setNull(true);
      sml.small_molecule_abundance_assay = smf.small_molecule_feature_abundance_assay;
      sml.small_molecule_abundance_study_variable[1].setNull(true);
      sml.small_molecule_abundance_variation_study_variable[1].setNull(true);

      smls.emplace_back(sml);
    }

    mztabm.setMetaData(m_meta_data);
    mztabm.setMSmallMoleculeEvidenceSectionRows(smes);
    mztabm.setMSmallMoleculeFeatureSectionRows(smfs);
    mztabm.setMSmallMoleculeSectionRows(smls);
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
