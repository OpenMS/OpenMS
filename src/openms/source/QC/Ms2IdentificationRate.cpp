// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow$
// $Authors: Patricia Scheil, Swenja Wagner$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/QC/Ms2IdentificationRate.h>
#include <algorithm>


namespace OpenMS
{
  // computes number of peptide identifications, number of ms2 spectra and ratio
  // data is stored in vector of structs
  void Ms2IdentificationRate::compute(const FeatureMap& feature_map, const MSExperiment& exp, bool assume_all_target)
  {
    // count ms2 spectra
    Size ms2_level_counter = getMS2Count_(exp);

    // count peptideIdentifications
    Size peptide_identification_counter {};

    auto f = [assume_all_target, &peptide_identification_counter](const PeptideIdentification& id) { peptide_identification_counter += isTargetPeptide_(id, assume_all_target); };

    // iterates through all PeptideIdentifications in FeatureMap, applies function f to all of them
    feature_map.applyFunctionOnPeptideIDs(f, true);

    writeResults_(peptide_identification_counter, ms2_level_counter);
  }

  void Ms2IdentificationRate::compute(const std::vector<PeptideIdentification>& pep_ids, const MSExperiment& exp, bool assume_all_target)
  {
    // count ms2 spectra
    Size ms2_level_counter = getMS2Count_(exp);

    // count peptideIdentifications
    Size peptide_identification_counter {};

    for (const auto& id : pep_ids)
    {
      peptide_identification_counter += isTargetPeptide_(id, assume_all_target);
    }

    writeResults_(peptide_identification_counter, ms2_level_counter);
  }


  const String& Ms2IdentificationRate::getName() const
  {
    return name_;
  }


  const std::vector<OpenMS::Ms2IdentificationRate::IdentificationRateData>& Ms2IdentificationRate::getResults() const
  {
    return rate_result_;
  }


  QCBase::Status Ms2IdentificationRate::requirements() const
  {
    return QCBase::Status() | QCBase::Requires::RAWMZML | QCBase::Requires::POSTFDRFEAT;
  }

  void Ms2IdentificationRate::addMetaDataMetricsToMzTab(MzTabMetaData& meta) const
  {
    // Adding MS2_ID_Rate to meta data
    const auto& ms2_irs = this->getResults();
    for (Size i = 0; i < ms2_irs.size(); ++i)
    {
      MzTabParameter ms2_ir {};
      ms2_ir.setCVLabel("MS2 identification rate");
      ms2_ir.setAccession("null");
      ms2_ir.setName("MS2_ID_Rate_" + String(i + 1));
      ms2_ir.setValue(String(100 * ms2_irs[i].identification_rate));
      meta.custom[meta.custom.size()] = ms2_ir;
    }
  }
  Size Ms2IdentificationRate::getMS2Count_(const MSExperiment& exp)
  {
    if (exp.empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MSExperiment is empty");
    }

    Size ms2_counter {};
    for (auto const& spec : exp.getSpectra())
    {
      if (spec.getMSLevel() == 2)
      {
        ++ms2_counter;
      }
    }

    if (ms2_counter == 0)
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No MS2 spectra found");
    }

    return ms2_counter;
  }

  bool Ms2IdentificationRate::isTargetPeptide_(const PeptideIdentification& id, bool all_targets)
  {
    if (id.getHits().empty())
    {
      return false;
    }
    if (all_targets)
    {
      return true;
    }
    if (!(id.getHits()[0].metaValueExists("target_decoy")))
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No target/decoy annotation found. If you want to continue regardless use -MS2_id_rate:assume_all_target");
    }
    // check for 'target' information, also allow "target+decoy" value
    String td_info(id.getHits()[0].getMetaValue("target_decoy"));
    return (td_info.find("target") == 0);
  }

  void Ms2IdentificationRate::writeResults_(Size pep_ids_count, Size ms2_spectra_count)
  {
    if (ms2_spectra_count < pep_ids_count)
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "There are more Identifications than MS2 spectra. Please check your data.");
    }

    // compute ratio
    double ratio = (double)pep_ids_count / ms2_spectra_count;

    // struct to store results
    IdentificationRateData id_rate_data {};

    // store results
    id_rate_data.num_peptide_identification = pep_ids_count;
    id_rate_data.num_ms2_spectra = ms2_spectra_count;
    id_rate_data.identification_rate = ratio;

    rate_result_.push_back(id_rate_data);
  }
} // namespace OpenMS
