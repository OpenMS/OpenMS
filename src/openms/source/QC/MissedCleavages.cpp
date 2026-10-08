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
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include <OpenMS/QC/MissedCleavages.h>
#include <iostream>

namespace OpenMS
{
  typedef std::map<UInt32, UInt32> MapU32;
  // digests the Sequence in PeptideHit and counts the number of missed cleavages
  void MissedCleavages::get_missed_cleavages_from_peptide_identification_(const ProteaseDigestion& digestor, MapU32& result, const UInt32& max_mc, QCBase::AnnotatedIdentification& id)
  {
    if (id.top == nullptr)
    {
      OPENMS_LOG_WARN << "There is a Peptideidentification(RT: " << id.rt << ", MZ: " << id.mz << ") without PeptideHits.\n";
      return;
    }
    std::vector<AASequence> digest_output;
    digestor.digest(id.top->getSequence(), digest_output);
    UInt32 num_mc = UInt32(digest_output.size() - 1);

    // warn if number of missed cleavages is greater than allowed maximum number of missed cleavages
    if (num_mc > max_mc)
    {
      OPENMS_LOG_WARN << "Observed number of missed cleavages: " << num_mc << " is greater than: " << max_mc
                      << " the allowed maximum number of missed cleavages during MS2-Search in: " << id.top->getSequence() << "\n";
    }

    ++result[num_mc];

    id.top->setMetaValue("missed_cleavages", num_mc);
  };

  void MissedCleavages::compute(std::vector<ProteinIdentification>& prot_ids, PeptideIdentificationList& pep_ids)
  {
    MapU32 result {};

    // Exception if ProteinIdentification is empty
    if (prot_ids.empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Missing information in ProteinIdentifications.");
    }

    std::string enzyme = prot_ids[0].getSearchParameters().digestion_enzyme.getName();
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
      auto id = QCBase::annotated(pep_id);
      get_missed_cleavages_from_peptide_identification_(digestor, result, max_mc, id);
    }

    mc_result_.push_back(result);
  }


  void MissedCleavages::compute(FeatureMap& fmap)
  {
    IdentificationDataConverter::editAsIdentificationData(fmap, [&](FeatureMap& map) {
      MapU32 result {};

      bool has_pepIDs = QCBase::hasPepID(map);
      if (!has_pepIDs)
      {
        mc_result_.push_back(result);
        return;
      }

      // if the FeatureMap is empty, result is 0
      if (map.empty())
      {
        OPENMS_LOG_WARN << "FeatureXML is empty.\n";
        mc_result_.push_back(result);
        return;
      }

      // Exception if there is no identification run (ProteinIdentification)
      const auto* search = QCBase::searchParameters(map.getIdentificationData());
      if (search == nullptr)
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Missing information in ProteinIdentifications.");
      }

      std::string enzyme = search->digestion_enzyme.getName();
      auto max_mc = search->missed_cleavages;

      // Exception if digestion enzyme is not given
      if (enzyme == "unknown_enzyme")
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No digestion enzyme in FeatureMap detected. No computation possible.");
      }

      // create a digestor, which doesn't allow any missed cleavages
      ProteaseDigestion digestor;
      digestor.setEnzyme(enzyme);
      digestor.setMissedCleavages(0);

      // the identifications of the features and the unassigned ones
      QCBase::annotateIdentifications(map, [&](Feature*, std::vector<QCBase::AnnotatedIdentification>& identifications) {
        for (auto& id : identifications)
        {
          get_missed_cleavages_from_peptide_identification_(digestor, result, max_mc, id);
        }
      });

      mc_result_.push_back(result);
    });
  }


  const std::string& MissedCleavages::getName() const
  {
    static const std::string& name = "MissedCleavages";
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
