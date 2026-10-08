// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Dominik Schmitz, Chris Bielow$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include "Contaminants.h"
#include <algorithm>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

using namespace std;

namespace OpenMS
{

  void Contaminants::compute(FeatureMap& features, const std::vector<FASTAFile::FASTAEntry>& contaminants)
  {
    IdentificationDataConverter::editAsIdentificationData(features, [&](FeatureMap& map) { computeNative_(map, contaminants); });
  }

  void Contaminants::computeNative_(FeatureMap& features, const std::vector<FASTAFile::FASTAEntry>& contaminants)
  {
    // empty FeatureMap
    if (features.empty())
    {
      OPENMS_LOG_WARN << "FeatureMap is empty"
                      << "\n";
    }
    // empty contaminants database
    if (contaminants.empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No contaminants provided.");
    }
    // fill the unordered set once with the digested contaminants database
    if (digested_db_.empty())
    {
      const auto* search = QCBase::searchParameters(features.getIdentificationData());
      if (search == nullptr)
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No proteinidentifications in FeatureMap.");
      }
      ProteaseDigestion digestor;
      std::string enzyme = search->digestion_enzyme.getName();

      // no enzyme is given
      if (enzyme == "unknown_enzyme")
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No digestion enzyme in FeatureMap detected. No computation possible.");
      }

      digestor.setEnzyme(enzyme);

      // get the missed cleavages for the digestor. If none are given, its default is 0.
      UInt missed_cleavages(search->missed_cleavages);
      digestor.setMissedCleavages(missed_cleavages);

      // digest the contaminants database and add the peptides into the unordered set
      for (const FASTAFile::FASTAEntry& fe : contaminants)
      {
        vector<AASequence> current_digest;
        digestor.digest(AASequence::fromString(fe.sequence), current_digest);

        // fill unordered set digested_db_ with digested sequences
        for (auto const& s : current_digest)
        {
          digested_db_.insert(s.toUnmodifiedString());
        }
      }
    }
    Int64 total = 0;
    Int64 cont = 0;
    double sum_total = 0.0;
    double sum_cont = 0.0;
    Int64 feature_has_no_sequence = 0;

    // Check if peptides of featureMap are contaminants or not and add is_contaminant = 0/1 to the top hit of the identifications.
    // If so, raise contaminants ratio.
    ContaminantsSummary final;
    UInt64 utotal = 0;
    UInt64 ucont = 0;
    QCBase::annotateIdentifications(features, [&](Feature* f, std::vector<QCBase::AnnotatedIdentification>& identifications) {
      if (f != nullptr)
      {
        if (identifications.empty())
        {
          ++feature_has_no_sequence;
          return;
        }
        for (auto& id : identifications)
        {
          // the identification of feature f has no hit
          if (id.top == nullptr)
          {
            ++feature_has_no_sequence;
            continue;
          }
          std::string key = (id.top->getSequence().toUnmodifiedString());
          this->compare_(key, *id.top, total, cont, sum_total, sum_cont, f->getIntensity());
        }
        return;
      }

      // save the contaminants ratio in object before searching through the unassigned identifications
      final.assigned_contaminants_ratio = (cont / double(total));

      final.empty_features.first = feature_has_no_sequence;
      final.empty_features.second = features.size();

      // Change the assigned contaminants ratio to total contaminants ratio by adding the unassigned.
      // Additionally save the unassigned contaminants ratio and add the is_contaminant = 0/1 to the top hit of the unassigned identifications.
      for (auto& id : identifications)
      {
        if (id.top == nullptr)
        {
          continue;
        }
        std::string key = (id.top->getSequence().toUnmodifiedString());
        ++utotal;

        // peptide is not in contaminant database
        if (!digested_db_.contains(key))
        {
          id.top->setMetaValue("is_contaminant", 0);
          continue;
        }

        // peptide is contaminant
        ++ucont;
        id.top->setMetaValue("is_contaminant", 1);
      }
    });

    total += utotal;
    cont += ucont;

    // save all ratios and the intensity to the object
    final.all_contaminants_ratio = (cont / double(total));
    final.unassigned_contaminants_ratio = (ucont / double(utotal));
    final.assigned_contaminants_intensity_ratio = (sum_cont / sum_total);


    // add the object to the results vector
    results_.push_back(final);
  }

  const std::string& Contaminants::getName() const
  {
    return name_;
  }

  const std::vector<Contaminants::ContaminantsSummary>& Contaminants::getResults()
  {
    return results_;
  }


  // Check if peptide is in contaminants database or not and add the is_contaminant = 0/1.
  // If so, raise the contaminant ratio.
  void Contaminants::compare_(const std::string& key, PeptideHit& pep_hit, Int64& total, Int64& cont, double& sum_total, double& sum_cont, double intensity)
  {
    ++total;
    sum_total += intensity;
    // peptide is not in contaminant database
    if (!digested_db_.contains(key))
    {
      pep_hit.setMetaValue("is_contaminant", 0);
      return;
    }
    // peptide is contaminant
    ++cont;
    sum_cont += intensity;
    pep_hit.setMetaValue("is_contaminant", 1);
  }

  QCBase::Status Contaminants::requirements() const
  {
    return (QCBase::Status(QCBase::Requires::POSTFDRFEAT) | QCBase::Requires::CONTAMINANTS);
  }


} // namespace OpenMS
