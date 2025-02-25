// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow$
// $Authors: Hendrik Beschorner, Lenny Kovac, Virginia Rossow$
// --------------------------------------------------------------------------

#include <OpenMS/QC/MQMsmsExporter.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <QtCore/QDir>
#include <cmath> // isnan
#include <fstream>
//#include <vector>

using namespace OpenMS;


MQMsms::MQMsms(const String& path)
{
  if (path.empty())
  {
    return;
  }
  filename_ = path + "/msms.txt";
  try
  {
    QString msms_path = QString::fromStdString(path);
    QDir().mkpath(msms_path);
    file_ = std::fstream(filename_, std::fstream::out);
  }
  catch (...)
  {
    OPENMS_LOG_FATAL_ERROR << filename_ << " wasn’t created" << std::endl;
    throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "out_msms");
  }
  exportHeader_();
}


MQMsms::~MQMsms()
{
  file_.close();
}


void MQMsms::exportHeader_()
{

  file_ << "Raw file" << "\t";
  file_ << "Scan number" << "\t"; // NA
  file_ << "Scan index" << "\t"; // NA
  file_ << "Sequence" << "\t";
  file_ << "Length" << "\t";
  file_ << "Missed cleavages" << "\t";
  file_ << "Modifications" << "\t";
  file_ << "Modified sequence" << "\t";
  //file_ << "Oxidation (M) Probabilities" << "\t"; --> not supported by OpenMS
  //file_ << "Oxidation (M) Score diffs" << "\t"; --> not supported by OpenMS
  file_ << "Acetyl (Protein N-term)" << "\t";
  file_ << "Oxidation (M)" << "\t";
  file_ << "Proteins" << "\t";
  file_ << "Charge" << "\t";
  file_ << "Fragmentation" << "\t"; // NA
  file_ << "Mass analyzer" << "\t"; // NA
  file_ << "Type" << "\t";
  file_ << "Scan event number" << "\t"; // NA
  file_ << "Isotope index" << "\t"; // NA
  file_ << "m/z" << "\t";
  file_ << "Mass" << "\t";
  file_ << "Mass error [ppm]" << "\t";
  file_ << "Mass error [Da]" << "\t";
  file_ << "Simple mass error [ppm]" << "\t"; // NA
  file_ << "Retention time" << "\t";
  file_ << "PEP" << "\t";
  file_ << "Score" << "\t";
  file_ << "Delta score" << "\t";
  file_ << "Score diff" << "\t";
  file_ << "Localization prob" << "\t";
  //file_ << "Combinatorics" << "\t"; --> not supported by OpenMS
  //file_ << "PIF" << "\t"; --> not practical to implement
  file_ << "Fraction of total spectrum" << "\t";
  file_ << "Base peak fraction" << "\t";
  file_ << "Precursor full scan number" << "\t"; // NA
  file_ << "Precursor Intensity" << "\t"; // NA
  file_ << "Precursor apex fraction" << "\t"; // NA
  file_ << "Precursor apex offset" << "\t"; // NA
  file_ << "Precursor apex offset time" << "\t"; // NA
  file_ << "Matches Intensities" << "\t"; // NA
  file_ << "Mass deviations [Da]" << "\t"; // NA
  file_ << "Mass deviations [ppm]" << "\t"; // NA
  file_ << "Masses" << "\t"; // NA
  file_ << "Number of matches" << "\t"; // NA
  file_ << "Intensity coverage" << "\t"; // NA
  file_ << "Peak coverage" << "\t"; // NA
  file_ << "Neutral loss level" << "\t"; // NA
  file_ << "ETD identification type" << "\t"; // NA
  file_ << "Reverse" << "\t"; // NA
  file_ << "All scores" << "\t"; // NA
  file_ << "All sequences" << "\t"; // NA
  file_ << "All modified sequences" << "\t"; // NA
  //file_ << "Reporter PIF" << "\t"; --> not supported by OpenMS
  //file_ << "Reporter fraction" << "\t"; --> not supported by OpenMS
  file_ << "id" << "\t";
  file_ << "Protein group IDs" << "\n";
  //file_ << "Peptide ID" << "\t"; --> not useful without the other MQ files
  //file_ << "Mod. peptide ID" << "\t"; --> not useful without the other MQ files
  //file_ << "Evidence ID" << "\t"; // --> in evidence.txt
  //file_ << "Oxidation (M) site IDs" << "\n"; --> not useful without the other MQ files

}

void MQMsms::exportRowFromFeature_(
  const Feature& f,
  const ConsensusMap& cmap,
  const Size c_feature_number,
  const String& raw_file,
  const std::multimap<String, std::pair<Size, Size>>& UIDs,
  const ProteinIdentification::Mapping& mp_f,
  const MSExperiment& exp,
  const std::map<String,String>& prot_mapper)
{

  // use struct common_outpots from the ExporterHelper
  MQExporterHelper::MQCommonOutputs common_outputs{f, cmap, c_feature_number, UIDs, mp_f, exp, prot_mapper};

  const PeptideHit* ptr_best_hit = nullptr; // the best hit referring to score
  const PeptideIdentification* ptr_best_id = nullptr;
  const ConsensusFeature& cf = cmap[c_feature_number];

  String type;
  if (MQExporterHelper::hasValidPepID_(f, c_feature_number, UIDs, mp_f))
  {
    type = "MULTI-MSMS";
    ptr_best_hit = &f.getPeptideIdentifications()[0].getHits()[0];
    ptr_best_id = &f.getPeptideIdentifications()[0];
  }
  else if (MQExporterHelper::hasPeptideIdentifications_(cf))
  {
    type = "MULTI-MATCH";
    ptr_best_hit = &cf.getPeptideIdentifications()[0].getHits()[0];
    ptr_best_id = &cf.getPeptideIdentifications()[0];
  }
  else
  {
    return; // no valid PepID; nothing to export
  }

  //const double& max_score = ptr_best_hit->getScore();
  const AASequence& pep_seq = ptr_best_hit->getSequence();

  if (pep_seq.empty())
  {
    return;
  }


  // what is written in the file in this exact order

  file_ << raw_file << "\t"; // raw file
  file_ << "NA" << "\t"; // Scan number
  file_ << "NA" << "\t"; // Scan index
  file_ << pep_seq.toUnmodifiedString() << "\t"; // Sequence
  file_ << pep_seq.size() << "\t";               // Length
  file_ << ptr_best_hit->getMetaValue("missed_cleavages", "NA") << "\t"; // missed cleavages


  file_ << common_outputs.modifications.str() << "\t"; // Modifications
  file_ << "_" << pep_seq << "_" << "\t"; // Modified Sequence
  file_ << common_outputs.acetyl << "\t"; // Acetyl (Protein N-term)
  file_ << common_outputs.oxidation.str() << "\t"; // Oxidation (M)

  const std::set<String>& accessions = ptr_best_hit->extractProteinAccessionsSet();
  file_ << ListUtils::concatenate(accessions, ";") << "\t";  // Proteins

  file_ << f.getCharge() << "\t"; // Charge

  file_ << ptr_best_id->getMetaValue("activation_method", "NA") << "\t"; // Fragmentation

  file_ << "NA" << "\t"; // Mass analyzer
  file_ << type << "\t"; // type

  file_ << ptr_best_id->getMetaValue("ScanEventNumber", "NA") << "\t"; // Scan event number

  file_ << "NA" << "\t"; // Isotope index
  file_ << f.getMZ() << "\t"; // M/Z
  file_ << pep_seq.getMonoWeight() << "\t"; // Mass
  file_ << common_outputs.mass_error_ppm.str() << "\t"; // Mass Error [ppm]
  file_ << common_outputs.mass_error_da.str() << "\t"; // Mass error [Da]

  file_ << "NA" << "\t"; // Simple mass error [ppm]

  f.metaValueExists("rt_raw_end") && f.metaValueExists("rt_raw_start") ?
    file_ << (double(f.getMetaValue("rt_raw_end")) - double(f.getMetaValue("rt_raw_start"))) / 60 << "\t" : file_
      << "NA" << "\t"; // Retention time

  ptr_best_hit->metaValueExists("PEP")? file_ << ptr_best_hit->getMetaValue("PEP") << "\t" : file_ << "\t"; // PEP
  file_ << ptr_best_hit->getScore() << "\t"; // Score
  f.metaValueExists("delta") ? file_ << (f.getMetaValue("delta")) << "\t" : file_ << "\t"; // Delta score

  file_ << "NA" << "\t"; // Score diff
  file_ << "NA" << "\t"; // Localization prob

  f.metaValueExists(Constants::UserParam::PSM_EXPLAINED_ION_CURRENT_USERPARAM) ? file_ << (f.getMetaValue(Constants::UserParam::PSM_EXPLAINED_ION_CURRENT_USERPARAM)) << "\t": file_ << "\t"; // Fraction of total spectrum

  file_ << common_outputs.base_peak_fraction.str() << "\t"; // Base peak fraction

  file_ << "NA" << "\t"; // Precursor full scan number
  file_ << "NA" << "\t"; // Precursor Intensity
  file_ << "NA" << "\t"; // Precursor apex fraction
  file_ << "NA" << "\t"; // Precursor apex offset
  file_ << "NA" << "\t"; // Precursor apex offset time
  file_ << "NA" << "\t"; // Matches Intensities
  file_ << "NA" << "\t"; // Mass deviations [Da]
  file_ << "NA" << "\t"; // Mass deviations [ppm]
  file_ << "NA" << "\t"; // Masses
  file_ << "NA" << "\t"; // Number of matches
  file_ << "NA" << "\t"; // Intensity coverage
  file_ << "NA" << "\t"; // Peak coverage
  file_ << "NA" << "\t"; // Neutral loss level
  file_ << "NA" << "\t"; // ETD identification type

  ptr_best_hit->getMetaValue("target_decoy") == "decoy" ? file_ << "1\t" : file_ << "\t"; // reverse

  file_ << "NA" << "\t"; // All scores
  file_ << "NA" << "\t"; // All sequences
  file_ << "NA" << "\t"; // All modified sequences

  file_ << id_ << "\t"; // ID
  ++id_;

  file_ << ListUtils::concatenate(accessions, ";")  << "\n"; // Protein group IDs


}

void MQMsms::exportFeatureMap(const FeatureMap& feature_map, const ConsensusMap& cmap, const MSExperiment& exp, const std::map<String,String>& prot_mapper)
{
  if (!MQExporterHelper::isValid(filename_))
  {
    OPENMS_LOG_ERROR << "MqMsms object is not valid." << std::endl;
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



