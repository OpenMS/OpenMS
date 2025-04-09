// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTabMFile.h>
#include <OpenMS/FORMAT/TextFile.h>

namespace OpenMS
{

  MzTabMFile::MzTabMFile()= default;

  MzTabMFile::~MzTabMFile()= default;

  void MzTabMFile::generateMzTabMMetaDataSection_(const MzTabMMetaData& md, StringList& sl) const
  {
    sl.push_back(String("MTD\tmzTab-version\t") + md.mz_tab_version.toCellString()); // mandatory
    sl.push_back(String("MTD\tmzTab-ID\t") + md.mz_tab_id.toCellString()); // mandatory

    if (!md.title.isNull())
    {
      String s = String("MTD\ttitle\t") + md.title.toCellString();
      sl.push_back(s);
    }

    if(!md.description.isNull())
    {
      String s = String("MTD\tdescription\t") + md.description.toCellString();
      sl.push_back(s);
    }

    for (const auto& sp : md.sample_processing)
    {
      String s = "MTD\tsample_processing[" + String(sp.first) + "]\t" + sp.second.toCellString();
      sl.push_back(s);
    }

    for (const auto& inst : md.instrument)
    {
      const MzTabInstrumentMetaData& imd = inst.second;

      if (!imd.name.isNull())
      {
        String s = "MTD\tinstrument[" + String(inst.first) + "]-name\t" +  imd.name.toCellString();
        sl.push_back(s);
      }

      if (!imd.source.isNull())
      {
        String s = "MTD\tinstrument[" + String(inst.first) + "]-source\t" + imd.source.toCellString();
        sl.push_back(s);
      }

      for (const auto& mit : imd.analyzer)
      {
        if (!mit.second.isNull())
        {
          String s = "MTD\tinstrument[" + String(inst.first) + "]-analyzer[" + String(mit.first) + "]\t" + mit.second.toCellString();
          sl.push_back(s);
        }
      }

      if (!imd.detector.isNull())
      {
        String s = "MTD\tinstrument[" + String(inst.first) + "]-detector\t" + imd.detector.toCellString();
        sl.push_back(s);
      }
    }

    for (const auto & sw : md.software)
    {
      MzTabSoftwareMetaData msmd = sw.second;

      String s = "MTD\tsoftware[" + String(sw.first) + "]\t" + msmd.software.toCellString(); // mandatory
      sl.push_back(s);

      for (const auto& setting : msmd.setting)
      {
        String s = "MTD\tsoftware[" + String(sw.first) + "]-setting[" + String(setting.first) + String("]\t") + setting.second.toCellString();
        sl.push_back(s);
      }
    }

    for (auto const& pub : md.publication)
    {
      String s = "MTD\tpublication[" + String(pub.first) + "]\t" + pub.second.toCellString();
      sl.push_back(s);
    }

    for (const auto& contact : md.contact)
    {
      const MzTabContactMetaData& md = contact.second;
      if (!md.name.isNull())
      {
        String s = "MTD\tcontact[" + String(contact.first) + "]-name\t" + md.name.toCellString();
        sl.push_back(s);
      }

      if (!md.affiliation.isNull())
      {
        String s = "MTD\tcontact[" + String(contact.first) + "]-affiliation\t" + md.affiliation.toCellString();
        sl.push_back(s);
      }

      if (!md.email.isNull())
      {
        String s = "MTD\tcontact[" + String(contact.first) + "]-email\t" + md.email.toCellString();
        sl.push_back(s);
      }
    }

    for (const auto& uri : md.uri)
    {
      String s = "MTD\turi[" + String(uri.first) + "]\t" + uri.second.toCellString();
      sl.push_back(s);
    }

    for (const auto& ext_study : md.external_study_uri)
    {
      String s = "MTD\texternal_study_uri[" + String(ext_study.first) + "]\t" + ext_study.second.toCellString();
      sl.push_back(s);
    }

    String s = String("MTD\tquantification_method\t") + md.quantification_method.toCellString(); // mandatory
    sl.push_back(s);

    for (const auto& sample : md.sample)
    {
      const MzTabSampleMetaData& msmd = sample.second;

      if (!msmd.description.isNull())
      {
        String s = "MTD\tsample[" + String(sample.first) + "]-description\t" + msmd.description.toCellString();
        sl.push_back(s);
      }

      for (const auto& species : msmd.species)
      {
        String s = "MTD\tsample[" + String(sample.first) + "]-species[" + String(species.first) + "]\t" + species.second.toCellString();
        sl.push_back(s);
      }

      for (const auto& tissue : msmd.tissue)
      {
        String s = "MTD\tsample[" + String(sample.first) + "]-tissue[" + String(tissue.first) + "]\t" + tissue.second.toCellString();
        sl.push_back(s);
      }

      for (const auto& cell_type : msmd.cell_type)
      {
        String s = "MTD\tsample[" + String(sample.first) + "]-cell_type[" + String(cell_type.first) + "]\t" + cell_type.second.toCellString();
        sl.push_back(s);
      }

      for (auto const& disease : msmd.disease)
      {
        String s = "MTD\tsample[" + String(sample.first) + "]-disease[" + String(disease.first) + "]\t" + disease.second.toCellString();
        sl.push_back(s);
      }

      for (const auto& custom : msmd.custom)
      {
        String s = "MTD\tsample[" + String(sample.first) + "]-custom[" + String(custom.first) + "]\t" + custom.second.toCellString();
        sl.push_back(s);
      }
    }

    for (const auto& ms_run : md.ms_run)
    {
      const MzTabMMSRunMetaData& rmmd = ms_run.second;

      String s = "MTD\tms_run[" + String(ms_run.first) + "]-location\t" + rmmd.location.toCellString(); // mandatory
      sl.push_back(s);

      if (!rmmd.instrument_ref.isNull())
      {
        String s = "MTD\tms_run[" + String(ms_run.first) + "]-instrument_ref\t" + rmmd.instrument_ref.toCellString();
        sl.push_back(s);
      }

      if (!rmmd.format.isNull())
      {
        String s = "MTD\tms_run[" + String(ms_run.first) + "]-format\t" + rmmd.format.toCellString();
        sl.push_back(s);
      }

      for (const auto& fragmentation_method : rmmd.fragmentation_method)
      {
        String s = "MTD\tms_run[" + String(ms_run.first) + "]-fragmentation_method[" + String(fragmentation_method.first) + "]\t" + fragmentation_method.second.toCellString();
        sl.push_back(s);
      }

      for (const auto& scan_polarity : rmmd.scan_polarity)
      {
        String s = "MTD\tms_run[" + String(ms_run.first) + "]-scan_polarity[" + String(scan_polarity.first) + "]\t" + scan_polarity.second.toCellString(); // mandatory
        sl.push_back(s);
      }

      if (!rmmd.hash.isNull())
      {
        String s = "MTD\tms_run[" + String(ms_run.first) + "]-hash\t" + rmmd.hash.toCellString();
        sl.push_back(s);
      }

      if (!rmmd.hash_method.isNull())
      {
        String s = "MTD\tms_run[" + String(ms_run.first) + "]-hash_method\t" + rmmd.hash_method.toCellString();
        sl.push_back(s);
      }
    }

    for (const auto& assay : md.assay)
    {
      const MzTabMAssayMetaData& amd = assay.second;

      String name = "MTD\tassay[" + String(assay.first) + "]\t" + amd.name.toCellString(); // mandatory
      sl.push_back(name);

      for (const auto& custom : amd.custom)
      {
        String s = "MTD\tms_run[" + String(assay.first) + "]-custom[" + String(custom.first) + "]\t" + custom.second.toCellString();
        sl.push_back(s);
      }

      if (!amd.external_uri.isNull())
      {
        String s = "MTD\tassay[" + String(assay.first) + "]-external_uri\t" + amd.external_uri.toCellString();
        sl.push_back(s);
      }

      if (!amd.sample_ref.isNull())
      {
        String s = "MTD\tassay[" + String(assay.first) + "]-sample_ref\tsample[" + amd.sample_ref.toCellString() + "]";
        sl.push_back(s);
      }

      String ms_run_ref = "MTD\tassay[" + String(assay.first) + "]-ms_run_ref\tms_run[" + amd.ms_run_ref.toCellString() + "]"; // mandatory
      sl.push_back(ms_run_ref);
    }

    for (const auto& sv : md.study_variable)
    {
      const MzTabMStudyVariableMetaData& svmd = sv.second;

      String name = "MTD\tstudy_variable[" + String(sv.first) + "]\t" + svmd.name.toCellString(); // mandatory
      sl.push_back(name);

      std::vector<MzTabString> strings;
      MzTabStringList refs_string;
      for (const auto& ref : svmd.assay_refs)
      {
        strings.emplace_back(MzTabString("assay[" + std::to_string(ref) + ']'));
      }
      refs_string.set(strings);
      String refs = "MTD\tstudy_variable[" + String(sv.first) + "]-assay_refs\t" + refs_string.toCellString(); // mandatory
      sl.push_back(refs);

      if (!svmd.average_function.isNull())
      {
        String s = "MTD\tstudy_variable[" + String(sv.first) + "]-average_function\t" + svmd.average_function.toCellString();
        sl.push_back(s);
      }

      if (!svmd.variation_function.isNull())
      {
        String s = "MTD\tstudy_variable[" + String(sv.first) + "]-variation_function\t" + svmd.variation_function.toCellString();
        sl.push_back(s);
      }

      String description = "MTD\tstudy_variable[" + String(sv.first) + "]-description\t" + svmd.description.toCellString(); // mandatory
      sl.push_back(description);

      for (const auto& factor : svmd.factors.get())
      {
        String s = "MTD\tstudy_variable[" + String(sv.first) + "]-factors\t" + factor.toCellString();
        sl.push_back(s);
      }
    }

    for (const auto& custom : md.custom)
    {
      String s = "MTD\tcustom[" + String(custom.first) + "]\t" + custom.second.toCellString();
      sl.push_back(s);
    }

    for (const auto& cv : md.cv)
    {
      const MzTabCVMetaData& cvmd = cv.second;

      String label = "MTD\tcv[" + String(cv.first) + "]-label\t" + cvmd.label.toCellString(); // mandatory
      sl.push_back(label);

      String full_name = "MTD\tcv[" + String(cv.first) + "]-full_name\t" + cvmd.full_name.toCellString();  // mandatory
      sl.push_back(full_name);

      String version = "MTD\tcv[" + String(cv.first) + "]-version\t" + cvmd.version.toCellString();  // mandatory
      sl.push_back(version);

      String url = "MTD\tcv[" + String(cv.first) + "]-uri\t" + cvmd.url.toCellString();  // mandatory
      sl.push_back(url);
    }

    for (const auto& db : md.database)
    {
      MzTabMDatabaseMetaData dbmd = db.second;

      String database = "MTD\tdatabase[" + String(db.first) + "]\t" + dbmd.database.toCellString();  // mandatory
      sl.push_back(database);

      String prefix = "MTD\tdatabase[" + String(db.first) + "]-prefix\t" + dbmd.prefix.toCellString();  // mandatory
      sl.push_back(prefix);

      String version = "MTD\tdatabase[" + String(db.first) + "]-version\t" + dbmd.version.toCellString();  // mandatory
      sl.push_back(version);

      String uri = "MTD\tdatabase[" + String(db.first) + "]-uri\t" + dbmd.uri.toCellString();  // mandatory
      sl.push_back(uri);
    }

    for (const auto& agent : md.derivatization_agent)
    {
      String s = "MTD\tderivatization_agent[" + String(agent.first) + "]-uri\t" + agent.second.toCellString();
      sl.push_back(s);
    }

    sl.push_back(String("MTD\tsmall_molecule-quantification_unit\t") + md.small_molecule_quantification_unit.toCellString()); // mandatory

    sl.push_back(String("MTD\tsmall_molecule_feature-quantification_unit\t") + md.small_molecule_feature_quantification_unit.toCellString()); // mandatory (feature section)

    sl.push_back(String("MTD\tsmall_molecule-identification_reliability\t") + md.small_molecule_identification_reliability.toCellString()); // mandatory

    for (const auto& id_conf : md.id_confidence_measure)
    {
      String s = "MTD\tid_confidence_measure[" + String(id_conf.first) + "]\t" + id_conf.second.toCellString(); // mandatory
      sl.push_back(s);
    }

    for (const auto& csm : md.colunit_small_molecule)
    {
      String s = "MTD\tcolunit_small_molecule\t" + csm.toCellString();
      sl.push_back(s);
    }

    for (const auto& csmf : md.colunit_small_molecule_feature)
    {
      String s = "MTD\tcolunit_small_molecule\t" + csmf.toCellString();
      sl.push_back(s);
    }

    for (const auto& csme : md.colunit_small_molecule_evidence)
    {
      String s = "MTD\tcolunit_small_molecule\t" + csme.toCellString();
      sl.push_back(s);
    }
  }

  String MzTabMFile::generateMzTabMSmallMoleculeHeader_(const MzTabMMetaData& meta, const std::vector<String>& optional_columns, size_t& n_columns) const
  {
    StringList header;
    header.emplace_back("SMH");
    header.emplace_back("SML_ID");
    header.emplace_back("SMF_ID_REFS");
    header.emplace_back("database_identifier");
    header.emplace_back("chemical_formula");
    header.emplace_back("smiles");
    header.emplace_back("inchi");
    header.emplace_back("chemical_name");
    header.emplace_back("uri");
    header.emplace_back("theoretical_neutral_mass");
    header.emplace_back("adduct_ions");
    header.emplace_back("reliability");
    header.emplace_back("best_id_confidence_measure");
    header.emplace_back("best_id_confidence_value");

    for (const auto& a : meta.assay)
    {
      header.emplace_back(String("abundance_assay[") + String(a.first) + String("]"));
    }

    for (const auto& a : meta.study_variable)
    {
      header.emplace_back(String("abundance_study_variable[") + String(a.first) + String("]"));
    }

    for (const auto& a : meta.study_variable)
    {
      header.emplace_back(String("abundance_variation_study_variable[") + String(a.first) + String("]"));
    }

    std::copy(optional_columns.begin(), optional_columns.end(), std::back_inserter(header));
    n_columns = header.size();
    return ListUtils::concatenate(header, "\t");
  }

  String MzTabMFile::generateMzTabMSmallMoleculeSectionRow_(const MzTabMSmallMoleculeSectionRow& row, const std::vector<String>& optional_columns, size_t& n_columns) const
  {
    StringList s;
    s.emplace_back("SML");
    s.emplace_back(row.sml_identifier.toCellString());
    s.emplace_back(row.smf_id_refs.toCellString());
    s.emplace_back(row.database_identifier.toCellString());
    s.emplace_back(row.chemical_formula.toCellString());
    s.emplace_back(row.smiles.toCellString());
    s.emplace_back(row.inchi.toCellString());
    s.emplace_back(row.chemical_name.toCellString());
    s.emplace_back(row.uri.toCellString());
    s.emplace_back(row.theoretical_neutral_mass.toCellString());
    s.emplace_back(row.adducts.toCellString());
    s.emplace_back(row.reliability.toCellString());
    s.emplace_back(row.best_id_confidence_measure.toCellString());
    s.emplace_back(row.best_id_confidence_value.toCellString());

    for (const auto& abundance_assay : row.small_molecule_abundance_assay)
    {
      s.emplace_back(abundance_assay.second.toCellString());
    }

    for (const auto& abundance_study_variable : row.small_molecule_abundance_study_variable)
    {
      s.emplace_back(abundance_study_variable.second.toCellString());
    }

    for (const auto& variation_study_variable : row.small_molecule_abundance_variation_study_variable)
    {
      s.emplace_back(variation_study_variable.second.toCellString());
    }

    MzTabFile::addOptionalColumnsToSectionRow_(optional_columns, row.opt_, s);
    n_columns = s.size();
    return ListUtils::concatenate(s, "\t");
  }

  String MzTabMFile::generateMzTabMSmallMoleculeFeatureHeader_(const MzTabMMetaData& meta, const std::vector<String>& optional_columns, size_t& n_columns) const
  {
   StringList header;
   header.emplace_back("SFH");
   header.emplace_back("SMF_ID");
   header.emplace_back("SME_ID_REFS");
   header.emplace_back("SME_ID_REF_ambiguity_code");
   header.emplace_back("adduct_ion");
   header.emplace_back("isotopomer");
   header.emplace_back("exp_mass_to_charge");
   header.emplace_back("charge");
   header.emplace_back("retention_time_in_seconds");
   header.emplace_back("retention_time_in_seconds_start");
   header.emplace_back("retention_time_in_seconds_end");

   for (const auto& a : meta.assay)
   {
     header.emplace_back(String("abundance_assay[") + String(a.first) + String("]"));
   }

   std::copy(optional_columns.begin(), optional_columns.end(), std::back_inserter(header));
   n_columns = header.size();
   return ListUtils::concatenate(header, "\t");
  }

  String MzTabMFile::generateMzTabMSmallMoleculeFeatureSectionRow_(const MzTabMSmallMoleculeFeatureSectionRow& row, const std::vector<String>& optional_columns, size_t& n_columns) const
  {
    StringList s;
    s.emplace_back("SMF");
    s.emplace_back(row.smf_identifier.toCellString());
    s.emplace_back(row.sme_id_refs.toCellString());
    s.emplace_back(row.sme_id_ref_ambiguity_code.toCellString());
    s.emplace_back(row.adduct.toCellString());
    s.emplace_back(row.isotopomer.toCellString());
    s.emplace_back(row.exp_mass_to_charge.toCellString());
    s.emplace_back(row.charge.toCellString());
    s.emplace_back(row.retention_time.toCellString());
    s.emplace_back(row.rt_start.toCellString());
    s.emplace_back(row.rt_end.toCellString());

    for (const auto& feature_abundance : row.small_molecule_feature_abundance_assay)
    {
      s.emplace_back(feature_abundance.second.toCellString());
    }

    MzTabFile::addOptionalColumnsToSectionRow_(optional_columns, row.opt_, s);
    n_columns = s.size();
    return ListUtils::concatenate(s, "\t");
  }

  String MzTabMFile::generateMzTabMSmallMoleculeEvidenceHeader_(const MzTabMMetaData& meta, const std::vector<String>& optional_columns, size_t& n_columns) const
  {
    StringList header;
    header.emplace_back("SEH");
    header.emplace_back("SME_ID");
    header.emplace_back("evidence_input_id");
    header.emplace_back("database_identifier");
    header.emplace_back("chemical_formula");
    header.emplace_back("smiles");
    header.emplace_back("inchi");
    header.emplace_back("chemical_name");
    header.emplace_back("uri");
    header.emplace_back("derivatized_form");
    header.emplace_back("adduct_ion");
    header.emplace_back("exp_mass_to_charge");
    header.emplace_back("charge");
    header.emplace_back("theoretical_mass_to_charge");
    header.emplace_back("spectra_ref");
    header.emplace_back("identification_method");
    header.emplace_back("ms_level");

    for (const auto& id_conf : meta.id_confidence_measure)
    {
      header.emplace_back(String("id_confidence_measure[") + String(id_conf.first) + String("]"));
    }

    header.emplace_back("rank");

    std::copy(optional_columns.begin(), optional_columns.end(), std::back_inserter(header));
    n_columns = header.size();
    return ListUtils::concatenate(header, "\t");
  }

  String MzTabMFile::generateMzTabMSmallMoleculeEvidenceSectionRow_(const MzTabMSmallMoleculeEvidenceSectionRow& row, const std::vector<String>& optional_columns, size_t& n_columns) const
  {
    StringList s;
    s.emplace_back("SME");
    s.emplace_back(row.sme_identifier.toCellString());
    s.emplace_back(row.evidence_input_id.toCellString());
    s.emplace_back(row.database_identifier.toCellString());
    s.emplace_back(row.chemical_formula.toCellString());
    s.emplace_back(row.smiles.toCellString());
    s.emplace_back(row.inchi.toCellString());
    s.emplace_back(row.chemical_name.toCellString());
    s.emplace_back(row.uri.toCellString());
    s.emplace_back(row.derivatized_form.toCellString());
    s.emplace_back(row.adduct.toCellString());
    s.emplace_back(row.exp_mass_to_charge.toCellString());
    s.emplace_back(row.charge.toCellString());
    s.emplace_back(row.calc_mass_to_charge.toCellString());
    s.emplace_back(row.spectra_ref.toCellString());
    s.emplace_back(row.identification_method.toCellString());
    s.emplace_back(row.ms_level.toCellString());

    for (const auto& id_conf : row.id_confidence_measure)
    {
      s.emplace_back(id_conf.second.toCellString());
    }

    s.emplace_back(row.rank.toCellString());

    MzTabFile::addOptionalColumnsToSectionRow_(optional_columns, row.opt_, s);
    n_columns = s.size();
    return ListUtils::concatenate(s, "\t");
  }

  void MzTabMFile::store(const String& filename, const MzTabM& mztab_m) const
  {
    OPENMS_LOG_INFO << "exporting identification data: \"" << filename << "\" to MzTab-M: " << std::endl;

    if (!(FileHandler::hasValidExtension(filename, FileTypes::MZTAB) || FileHandler::hasValidExtension(filename, FileTypes::TSV)))
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "invalid file extension, expected '"
                                                                                                    + FileTypes::typeToName(FileTypes::MZTAB) + "' or '" + FileTypes::typeToName(FileTypes::TSV) + "'");
    }

    StringList out;
    generateMzTabMMetaDataSection_(mztab_m.getMetaData(), out);

    size_t n_sml_header_columns = 0;
    out.emplace_back("");
    out.emplace_back(generateMzTabMSmallMoleculeHeader_(mztab_m.getMetaData(),
                                                        mztab_m.getMSmallMoleculeOptionalColumnNames(),
                                                        n_sml_header_columns));

    size_t n_sml_section_columns = 0;
    const MzTabMSmallMoleculeSectionRows& sm_section = mztab_m.getMSmallMoleculeSectionRows();
    for (const auto& sms_row: sm_section)
    {
      out.emplace_back(generateMzTabMSmallMoleculeSectionRow_(sms_row,
                                                              mztab_m.getMSmallMoleculeOptionalColumnNames(),
                                                              n_sml_section_columns));

      OPENMS_POSTCONDITION(n_sml_header_columns == n_sml_section_columns,
                           "The number of columns of the small molecule header do not assort to the number of columns of the small molecule section row.")
    }

    size_t n_smf_header_columns = 0;

    out.emplace_back("");
    out.emplace_back(generateMzTabMSmallMoleculeFeatureHeader_(mztab_m.getMetaData(),
                                                               mztab_m.getMSmallMoleculeFeatureOptionalColumnNames(),
                                                               n_smf_header_columns));

    size_t n_smf_section_columns = 0;
    const MzTabMSmallMoleculeFeatureSectionRows& feature_section = mztab_m.getMSmallMoleculeFeatureSectionRows();
    for (const auto& smf_row : feature_section)
    {
      out.emplace_back(generateMzTabMSmallMoleculeFeatureSectionRow_(smf_row,
                                                                     mztab_m.getMSmallMoleculeFeatureOptionalColumnNames(),
                                                                     n_smf_section_columns));

      OPENMS_POSTCONDITION(n_smf_header_columns == n_smf_section_columns,
                           "The number of columns of the small molecule feature header do not assort to the number of columns of the small molecule feature section row.")
    }

    size_t n_sme_header_columns = 0;
    out.emplace_back("");
    out.emplace_back(generateMzTabMSmallMoleculeEvidenceHeader_(mztab_m.getMetaData(),
                                                                mztab_m.getMSmallMoleculeEvidenceOptionalColumnNames(),
                                                                n_sme_header_columns));

    size_t n_sme_section_columns = 0;
    const MzTabMSmallMoleculeEvidenceSectionRows& evidence_section = mztab_m.getMSmallMoleculeEvidenceSectionRows();
    for (const auto& sme_row : evidence_section)
    {
      out.emplace_back(generateMzTabMSmallMoleculeEvidenceSectionRow_(sme_row,
                                                                      mztab_m.getMSmallMoleculeEvidenceOptionalColumnNames(),
                                                                      n_sme_section_columns));

      OPENMS_POSTCONDITION(n_sme_header_columns == n_sme_section_columns,
                           "The number of columns in the small molecule evidence header do not assort to the number of columns of the small molecule evidence section row.")
    }

    TextFile tmp_out;
    for (TextFile::ConstIterator it = out.begin(); it != out.end(); ++it)
    {
      tmp_out.addLine(*it);
    }
    tmp_out.store(filename);
  }

}
