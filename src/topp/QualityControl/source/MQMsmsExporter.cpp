// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow$
// $Authors: Hendrik Beschorner, Lenny Kovac, Virginia Rossow$
// --------------------------------------------------------------------------

#include "MQMsmsExporter.h"

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>

#include <OpenMS/SYSTEM/PathUtils.h>
#include <filesystem>
#include <cmath> // isnan
#include <fstream>
//#include <vector>

using namespace OpenMS;


MQMsms::MQMsms(const std::string& path)
{
  if (path.empty())
  {
    return;
  }
  filename_ = path + "/msms.txt";
  try
  {
    std::filesystem::create_directories(to_path(path));
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
  const std::vector<MQExporterHelper::LinkedIdentification>& f_ids,
  const std::vector<MQExporterHelper::LinkedIdentification>& cf_ids,
  const std::string& raw_file,
  const MSExperiment& exp,
  const std::map<std::string, std::string>& prot_mapper)
{
  std::string type;
  const auto* best_id = MQExporterHelper::bestIdentification(f_ids, cf_ids, type);
  if (best_id == nullptr)
  {
    return; // no valid identification; nothing to export
  }
  const PeptideHit* ptr_best_hit = &*best_id->best; // the best hit referring to score
  const IdentificationData::Observation* ptr_best_id = &best_id->observation;

  // use struct common_outpots from the ExporterHelper
  MQExporterHelper::MQCommonOutputs common_outputs{f, *ptr_best_hit, exp, prot_mapper};

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

  const std::set<std::string>& accessions = ptr_best_hit->extractProteinAccessionsSet();
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

  ptr_best_hit->isDecoy() ? file_ << "1\t" : file_ << "\t"; // reverse

  file_ << "NA" << "\t"; // All scores
  file_ << "NA" << "\t"; // All sequences
  file_ << "NA" << "\t"; // All modified sequences

  file_ << id_ << "\t"; // ID
  ++id_;

  file_ << ListUtils::concatenate(accessions, ";")  << "\n"; // Protein group IDs


}

void MQMsms::exportFeatureMap(const FeatureMap& feature_map, const ConsensusMap& cmap, const MSExperiment& exp, const std::map<std::string, std::string>& prot_mapper)
{
  if (!MQExporterHelper::isValid(filename_))
  {
    OPENMS_LOG_ERROR << "MqMsms object is not valid." << std::endl;
    throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename_);
  }
  const std::map<Size, Size>& fTc = MQExporterHelper::makeFeatureUIDtoConsensusMapIndex_(cmap);
  StringList spectra_data;
  feature_map.getPrimaryMSRunPath(spectra_data);
  std::string raw_file = File::basename(spectra_data.empty() ? feature_map.getLoadedFilePath() : spectra_data[0]);

  // the identifications of the maps (maps with peptide identifications are read through a converted copy)
  std::optional<FeatureMap> converted_features;
  const FeatureMap& features = IdentificationDataConverter::withIdentificationData(feature_map, converted_features);
  std::optional<ConsensusMap> converted_consensus;
  const ConsensusMap& consensus = IdentificationDataConverter::withIdentificationData(cmap, converted_consensus);

  for (const Feature& f : features)
  {
    const Size& f_id = f.getUniqueId();
    const auto& c_id = fTc.find(f_id);
    if (c_id != fTc.end())
    {
      exportRowFromFeature_(f, MQExporterHelper::linkedIdentifications(f, features.getIdentificationData()),
                            MQExporterHelper::linkedIdentifications(consensus[c_id->second], consensus.getIdentificationData()),
                            raw_file, exp, prot_mapper);
    }
    else
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Feature in FeatureMap has no associated ConsensusFeature.");
    }
  }
  file_.flush();
}



