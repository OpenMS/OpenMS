// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Swenja Wagner, Patricia Scheil $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/QC/MissedCleavages.h>
#include <iostream>

namespace OpenMS
{
  typedef std::map<UInt32, UInt32> MapU32;
  // digests the Sequence in PeptideHit and counts the number of missed cleavages
  void MissedCleavages::get_missed_cleavages_from_peptide_identification_(const ProteaseDigestion& digestor, MapU32& result, const UInt32& max_mc, PeptideIdentification& pep_id)
  {
    if (pep_id.getHits().empty())
    {
      OPENMS_LOG_WARN << "There is a Peptideidentification(RT: " << pep_id.getRT() << ", MZ: " << pep_id.getMZ() << ") without PeptideHits.\n";
      return;
    }
    std::vector<AASequence> digest_output;
    digestor.digest(pep_id.getHits()[0].getSequence(), digest_output);
    UInt32 num_mc = UInt32(digest_output.size() - 1);

    // warn if number of missed cleavages is greater than allowed maximum number of missed cleavages
    if (num_mc > max_mc)
    {
      OPENMS_LOG_WARN << "Observed number of missed cleavages: " << num_mc << " is greater than: " << max_mc
                      << " the allowed maximum number of missed cleavages during MS2-Search in: " << pep_id.getHits()[0].getSequence() << "\n";
    }

    ++result[num_mc];

    pep_id.getHits()[0].setMetaValue("missed_cleavages", num_mc);
  };

  void MissedCleavages::compute(std::vector<ProteinIdentification>& prot_ids, std::vector<PeptideIdentification>& pep_ids)
  {
    MapU32 result {};

    // Exception if ProteinIdentification is empty
    if (prot_ids.empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Missing information in ProteinIdentifications.");
    }

    String enzyme = prot_ids[0].getSearchParameters().digestion_enzyme.getName();
    auto max_mc = prot_ids[0].getSearchParameters().missed_cleavages;

    // Exception if digestion enzyme is not given
    if (enzyme == "unknown_enzyme")
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No digestion enzyme in ID data detected. No computation possible.");
    }

    // create a digestor, which doesn't allow any missed cleavages
    ProteaseDigestion digestor;
    digestor.setEnzyme(enzyme);
    digestor.setMissedCleavages(0);

    for (PeptideIdentification& pep_id : pep_ids)
    {
      get_missed_cleavages_from_peptide_identification_(digestor, result, max_mc, pep_id);
    }

    mc_result_.push_back(result);
  }


  void MissedCleavages::compute(FeatureMap& fmap)
  {
    MapU32 result {};

    bool has_pepIDs = QCBase::hasPepID(fmap);
    if (!has_pepIDs)
    {
      mc_result_.push_back(result);
      return;
    }

    // if the FeatureMap is empty, result is 0
    if (fmap.empty())
    {
      OPENMS_LOG_WARN << "FeatureXML is empty.\n";
      mc_result_.push_back(result);
      return;
    }

    // Exception if ProteinIdentification is empty
    if (fmap.getProteinIdentifications().empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Missing information in ProteinIdentifications.");
    }

    String enzyme = fmap.getProteinIdentifications()[0].getSearchParameters().digestion_enzyme.getName();
    auto max_mc = fmap.getProteinIdentifications()[0].getSearchParameters().missed_cleavages;

    // Exception if digestion enzyme is not given
    if (enzyme == "unknown_enzyme")
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No digestion enzyme in FeatureMap detected. No computation possible.");
    }

    // create a digestor, which doesn't allow any missed cleavages
    ProteaseDigestion digestor;
    digestor.setEnzyme(enzyme);
    digestor.setMissedCleavages(0);

    // small lambda function to apply get_missed_cleavages_from_peptide_identification on pep_ids
    auto l = [&](PeptideIdentification& pep_id) {
      get_missed_cleavages_from_peptide_identification_(digestor, result, max_mc, pep_id);
      return;
    };
    // iterate through all PeptideIdentifications of a given FeatureMap and applies the given lambda function
    fmap.applyFunctionOnPeptideIDs(l);

    mc_result_.push_back(result);
  }


  const String& MissedCleavages::getName() const
  {
    static const String& name = "MissedCleavages";
    return name;
  }


  const std::vector<MapU32>& MissedCleavages::getResults() const
  {
    return mc_result_;
  }


  QCBase::Status MissedCleavages::requirements() const
  {
    return QCBase::Status() | QCBase::Requires::POSTFDRFEAT;
  }
} // namespace OpenMS
