// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow$
// $Authors: Valentin Noske, Vincent Musch$
// --------------------------------------------------------------------------

#include <OpenMS/QC/MQEvidenceExporter.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <QtCore/QDir>
#include <cmath> // isnan
#include <fstream>
#include <vector>


using namespace OpenMS;

MQEvidence::MQEvidence(const String& path)
{
  if (path.empty())
  {
    return;
  }
  filename_ = path + "/evidence.txt";
  try
  {
    QString evi_path = QString::fromStdString(path);
    QDir().mkpath(evi_path);
    file_ = std::fstream(filename_, std::fstream::out);
  }
  catch (...)
  {
    OPENMS_LOG_FATAL_ERROR << filename_ << " wasn’t created" << std::endl;
    throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "out_evd");
  }
  exportHeader_();
}

MQEvidence::~MQEvidence()
{
  file_.close();
}

void MQEvidence::exportHeader_()
{
  file_ << "Sequence" << "\t";
  file_ << "Length" << "\t";
  file_ << "Modifications" << "\t";
  file_ << "Modified sequence" << "\t";
  // file_ << "Oxidation (M) Probabilities" << "\t"; not supported by OpenMS
  // file_ << "Oxidation (M) Score Diffs" << "\t"; not supported by OpenMS
  file_ << "Acetyl (Protein N-term)" << "\t";
  file_ << "Oxidation (M)" << "\t";
  file_ << "Missed cleavages" << "\t";
  file_ << "Proteins" << "\t";
  file_ << "Leading Proteins" << "\t";
  file_ << "Leading Razor Protein" << "\t";
  file_ << "Gene Names" << "\t";
  file_ << "Protein Names" << "\t";
  file_ << "Type" << "\t";
  file_ << "Raw file" << "\t";
  // file_ << "Fraction" << "\t"; not in this workflow
  file_ << "MS/MS m/z" << "\t";
  file_ << "Charge" << "\t";
  file_ << "m/z" << "\t";
  file_ << "Mass" << "\t";
  file_ << "Resolution" << "\t";
  file_ << "Uncalibrated - Calibrated m/z [ppm]" << "\t";
  file_ << "Uncalibrated - Calibrated m/z [Da]" << "\t";
  file_ << "Mass Error [ppm]" << "\t";
  file_ << "Mass Error [Da]" << "\t";
  file_ << "Uncalibrated Mass Error [ppm]" << "\t";
  file_ << "Uncalibrated Mass Error [Da]" << "\t";
  // file_ << "Max intensity m/z 0" << "\t";
  file_ << "Retention time" << "\t";
  file_ << "Retention length" << "\t";
  file_ << "Calibrated retention time" << "\t";
  file_ << "Calibrated retention time start" << "\t";
  file_ << "Calibrated retention time finish" << "\t";
  file_ << "Retention time calibration" << "\t";
  file_ << "Match time difference" << "\t";
  file_ << "Match m/z difference" << "\t";
  file_ << "Match q-value" << "\t";
  file_ << "Match score" << "\t";
  file_ << "Number of data points" << "\t";
  // file_ << "Number of scans" << "\t"; not practical to implement
  file_ << "Number of isotopic peaks" << "\t";
  // file_ << "PIF" << "\t"; not practical to implement
  file_ << "Fraction of total spectrum" << "\t";
  file_ << "Base peak fraction" << "\t";
  file_ << "PEP" << "\t";
  file_ << "MS/MS Count" << "\t";
  file_ << "MS/MS Scan Number" << "\t";
  file_ << "Score" << "\t";
  file_ << "Delta score" << "\t";
  // file_ << "Combinatorics" << "\t"; not supported by OpenMS
  file_ << "Intensity" << "\t";
  /*file_ << "Reporter intensity 0" << "\t"; not supported by the given data
  file_ << "Reporter intensity 1" << "\t";
  file_ << "Reporter intensity 2" << "\t";
  file_ << "Reporter intensity 3" << "\t";
  file_ << "Reporter intensity 4" << "\t";
  file_ << "Reporter intensity 5" << "\t";
  file_ << "Reporter intensity not corrected 0" << "\t";
  file_ << "Reporter intensity not corrected 1" << "\t";
  file_ << "Reporter intensity not corrected 2" << "\t";
  file_ << "Reporter intensity not corrected 3" << "\t";
  file_ << "Reporter intensity not corrected 4" << "\t";
  file_ << "Reporter intensity not corrected 5" << "\t";
  file_ << "Reporter PIF" << "\t";
  file_ << "Reporter fraction" << "\t"; */
  file_ << "Reverse" << "\t";
  file_ << "Potential contaminant" << "\t";
  file_ << "id" << "\t";
  file_ << "Protein group IDs" << "\n";
  /*file_ << "Peptide ID" << "\t"; not useful without the other MQ files
  file_ << "Mod. peptide ID" << "\t";
  file_ << "MS/MS IDs" << "\t";
  file_ << "Best MS/MS" << "\t";
  file_ << "AIF MS/MS IDs" << "\t";
  file_ << "Oxidation (M) site IDs" << "\n"; */
}

void MQEvidence::exportRowFromFeature_(
  const Feature& f,
  const ConsensusMap& cmap,
  const Size c_feature_number,
  const String& raw_file,
  const std::multimap<String, std::pair<Size, Size>>& UIDs,
  const ProteinIdentification::Mapping& mp_f,
  const MSExperiment& exp,
  const std::map<String,String>& prot_mapper)
{

  MQExporterHelper::MQCommonOutputs common_outputs{f, cmap, c_feature_number, UIDs, mp_f, exp, prot_mapper};

  const PeptideHit* ptr_best_hit; // the best hit referring to score
  const ConsensusFeature& cf = cmap[c_feature_number];
  Size pep_ids_size = 0;
  String type;
  if (MQExporterHelper::hasValidPepID_(f, c_feature_number, UIDs, mp_f))
  {
    for (Size i = 1; i < f.getPeptideIdentifications().size(); ++i) // for msms-count
    {
      if (!f.getPeptideIdentifications()[i].getHits().empty())
      {
        if (f.getPeptideIdentifications()[i].getHits()[0].getSequence() == f.getPeptideIdentifications()[0].getHits()[0].getSequence())
        {
          ++pep_ids_size;
        }
        else
          break;
      }
    }
    type = "MULTI-MSMS";
    ptr_best_hit = &f.getPeptideIdentifications()[0].getHits()[0];
  }
  else if (MQExporterHelper::hasPeptideIdentifications_(cf))
  {
    type = "MULTI-MATCH";
    ptr_best_hit = &cf.getPeptideIdentifications()[0].getHits()[0];
  }
  else
  {
    return; // no valid PepID; nothing to export
  }
  const AASequence& pep_seq = ptr_best_hit->getSequence();

  if (pep_seq.empty())
  {
    return;
  }

  file_ << pep_seq.toUnmodifiedString() << "\t"; // Sequence
  file_ << pep_seq.size() << "\t";               // Length


  file_ << common_outputs.modifications.str() << "\t"; // Modifications
  file_ << "_" << pep_seq << "_" << "\t"; // Modified Sequence
  file_ << common_outputs.acetyl << "\t", // Acetyl (Protein N-term)
  file_ << common_outputs.oxidation.str() << "\t"; // Oxidation (M)

  file_ << ptr_best_hit->getMetaValue("missed_cleavages", "NA") << "\t"; // missed cleavages
  const std::set<String>& accessions = ptr_best_hit->extractProteinAccessionsSet();
  file_ << ListUtils::concatenate(accessions, ";") << "\t";  // Proteins
  file_ << ptr_best_hit->getPeptideEvidences()[0].getProteinAccession() << "\t"; // Leading Proteins
  file_ << ptr_best_hit->getPeptideEvidences()[0].getProteinAccession() << "\t"; // Leading Razor Proteins

  file_ << common_outputs.gene_names.str() << "\t"; // Gene Names
  file_ << common_outputs.protein_names.str() << "\t"; // Protein Names
  file_ << type << "\t"; //type

  file_ << raw_file << "\t"; // Raw File
  file_ << common_outputs.msms_mz.str() << "\t"; // MS/MS m/z
  file_ << f.getCharge() << "\t";           // Charge
  file_ << f.getMZ() << "\t";               // MZ
  file_ << pep_seq.getMonoWeight() << "\t"; // Mass
  file_ << f.getMZ() / f.getWidth() << "\t"; // Resolution

  file_ << common_outputs.uncalibrated_calibrated_mz_ppm.str() << "\t"; // Uncalibrated - Calibrated m/z [ppm]
  file_ << common_outputs.uncalibrated_calibrated_mz_mda.str() << "\t"; // Uncalibrated - Calibrated m/z [Da]
  file_ << common_outputs.mass_error_ppm.str() << "\t"; // Mass error [ppm]
  file_ << common_outputs.mass_error_da.str() << "\t"; // Mass error [Da]
  file_ << common_outputs.uncalibrated_mass_error_ppm.str() << "\t"; // Uncalibrated Mass error [ppm]
  file_ << common_outputs.uncalibrated_mass_error_da .str() << "\t"; // Uncalibrated Mass error [Da]

  file_ << f.getRT() / 60 << "\t"; // Retention time in min.

  // RET LENGTH
  f.metaValueExists("rt_raw_end") && f.metaValueExists("rt_raw_start") ?
    file_ << (double(f.getMetaValue("rt_raw_end")) - double(f.getMetaValue("rt_raw_start"))) / 60 << "\t" : file_
      << "NA" << "\t";

  if (f.metaValueExists("rt_align"))
  {
    file_ << double(f.getMetaValue("rt_align")) / 60 << "\t";   // Calibrated Retention Time
  }
  f.metaValueExists("rt_align_start") ? file_ << double(f.getMetaValue("rt_align_start")) / 60 << "\t" :
                                        file_ << "NA" << "\t"; //  Calibrated retention time start
  f.metaValueExists("rt_align_end") ? file_ << double(f.getMetaValue("rt_align_end")) / 60 << "\t" :
                                      file_ << "NA" << "\t"; // Calibrated retention time end
  if (f.metaValueExists("rt_align"))
  {
    file_ << ((double(f.getMetaValue("rt_align")) - f.getRT() ) / 60) << "\t"; // Retention time calibration
  }
  else
  {
    file_ << "NA" << "\t"; // calibrated retention time
    file_ << "NA" << "\t"; // Retention time calibration
  }

  if (type == "MULTI-MSMS")
  {
    file_ << "NA" << "\t"; // Match time diff
    file_ << "NA" << "\t"; // Match mz diff
  }
  else
  {
    f.metaValueExists("rt_align") ? file_ << double(f.getMetaValue("rt_align")) - cmap[c_feature_number].getRT() << "\t" :
                                    file_ << "NA"
                                          << "\t";               // Match time diff
    file_ << f.getMZ() - cmap[c_feature_number].getMZ() << "\t"; // Match mz diff
  }

  const PeptideHit* cf_ptr_best_hit = &cf.getPeptideIdentifications()[0].getHits()[0];
  cf_ptr_best_hit->metaValueExists("qvalue")? file_ << cf_ptr_best_hit->getMetaValue("qvalue") << "\t" : file_ << "\t"; // Match q-value

  file_ << cf_ptr_best_hit->getScore() << "\t"; // Match score

  f.metaValueExists(Constants::UserParam::NUM_OF_DATAPOINTS) ? file_ << (f.getMetaValue(Constants::UserParam::NUM_OF_DATAPOINTS)) << "\t": file_ << "\t"; // Number of data points
  file_ << f.getConvexHulls().size() << "\t"; // Number of isotopic peaks

  f.metaValueExists(Constants::UserParam::PSM_EXPLAINED_ION_CURRENT_USERPARAM) ? file_ << (f.getMetaValue(Constants::UserParam::PSM_EXPLAINED_ION_CURRENT_USERPARAM)) << "\t": file_ << "\t"; // Fraction of total spectrum

  file_ << common_outputs.base_peak_fraction.str() << "\t"; // Base peak fraction
  ptr_best_hit->metaValueExists("PEP")? file_ << ptr_best_hit->getMetaValue("PEP") << "\t" : file_ << "\t"; // PEP

  file_ << pep_ids_size << "\t"; // MS/MS count
  f.metaValueExists("spectrum_index") ? file_ << (f.getMetaValue("spectrum_index")) << "\t" : file_ << "\t";// MS/MS Scan Number
  file_ << ptr_best_hit->getScore() << "\t"; // Score

  f.metaValueExists("delta") ? file_ << (f.getMetaValue("delta")) << "\t" : file_ << "\t"; // Delta score

  file_ << f.getIntensity() << "\t"; // Intensity

  ptr_best_hit->getMetaValue("target_decoy") == "decoy" ? file_ << "1\t" : file_ << "\t"; // reverse

  String pot_containment = ptr_best_hit->getMetaValue("is_contaminant", "NA");
  if (pot_containment == "1")
  {
    file_ << "+"
          << "\t"; // Potential contaminant
  }
  else
  {
    file_ << "\t";
  }

  file_ << id_ << "\t"; // ID
  ++id_;

  file_ << ListUtils::concatenate(accessions, ";")  << "\n"; // Protein group IDs

}

void MQEvidence::exportFeatureMap(const FeatureMap& feature_map, const ConsensusMap& cmap, const MSExperiment& exp, const std::map<String,String>& prot_mapper)
{
if (!MQExporterHelper::isValid(filename_))
  {
    OPENMS_LOG_ERROR << "MqEvidence object is not valid." << std::endl;
    throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename_);
  }
  const std::map<Size, Size>& fTc = MQExporterHelper::makeFeatureUIDtoConsensusMapIndex_(cmap);
  StringList spectra_data;
  feature_map.getPrimaryMSRunPath(spectra_data);
  String raw_file = File::basename(spectra_data.empty() ? feature_map.getLoadedFilePath() : spectra_data[0]);

  ProteinIdentification::Mapping mp_f;
  mp_f.create(feature_map.getProteinIdentifications());

  std::multimap<String, std::pair<Size, Size>> UIDs = PeptideIdentification::buildUIDsFromAllPepIDs(cmap);

  for (const Feature& f : feature_map)
  {
    const Size& f_id = f.getUniqueId();
    const auto& c_id = fTc.find(f_id);
    if (c_id != fTc.end())
    {
      exportRowFromFeature_(f, cmap, c_id->second, raw_file, UIDs, mp_f, exp, prot_mapper);
    }
    else
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Feature in FeatureMap has no associated ConsensusFeature.");
    }
  }
  file_.flush();
}
