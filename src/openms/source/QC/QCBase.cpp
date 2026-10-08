// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow, Tom Waschischeck $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/QC/QCBase.h>
#include <OpenMS/METADATA/DataProcessingUtils.h>
#include <OpenMS/METADATA/ID/IdentificationDataAdapter.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <optional>

namespace OpenMS
{
  const std::string QCBase::names_of_requires[] = {"fail", "raw.mzML", "postFDR.featureXML", "preFDR.featureXML", "contaminants.fasta", "trafoAlign.trafoXML", "id.idXML"};
  static_assert(sizeof(QCBase::names_of_requires) / sizeof(QCBase::names_of_requires[0])
                  == static_cast<Size>(QCBase::Requires::SIZE_OF_REQUIRES),
                "names_of_requires must have one entry per QCBase::Requires value");

  const std::string QCBase::names_of_toleranceUnit[] = {"auto", "ppm", "da"};

  QCBase::SpectraMap::SpectraMap(const MSExperiment& exp)
  {
    calculateMap(exp);
  }

  void QCBase::SpectraMap::calculateMap(const MSExperiment& exp)
  {
    nativeid_to_index_.clear();
    for (Size i = 0; i < exp.size(); ++i)
    {
      nativeid_to_index_[exp[i].getNativeID()] = i;
    }
  }

  UInt64 QCBase::SpectraMap::at(const std::string& identifier) const
  {
    if (const auto& it = nativeid_to_index_.find(identifier); it == nativeid_to_index_.end())
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,std::string("No spectrum with identifier '") + identifier + "' in MSExperiment!");
    }
    else
    {
      return it->second;
    }
  }

  void QCBase::SpectraMap::clear()
  {
    nativeid_to_index_.clear();
  }

  bool QCBase::SpectraMap::empty() const
  {
    return nativeid_to_index_.empty();
  }

  Size QCBase::SpectraMap::size() const
  {
    return nativeid_to_index_.size();
  }

  namespace
  {
    using ID = IdentificationData;

    bool hasIdentifications(const ID& data)
    {
      for (const auto& run : data.getRuns())
        for (const auto& source : run.getSources())
          if (!source.identifications.empty()) return true;
      return false;
    }

    /// The identifications of a feature map as QC metrics see them, with editable copies of their observations and top matches
    class Annotations
    {
    public:
      std::vector<QCBase::AnnotatedIdentification> identifications(const std::vector<ID::QueryMatches>& entries)
      {
        std::vector<QCBase::AnnotatedIdentification> result;
        result.reserve(entries.size());
        for (const auto& entry : entries)
        {
          QCBase::AnnotatedIdentification item;
          auto observation = observations_.try_emplace(ID::QueryReference {entry.run->getUuid(), entry.query->getId()});
          if (observation.second) observation.first->second = {ID::Observation(*entry.query), ID::Observation(*entry.query)};
          auto& annotated = observation.first->second.second;
          if (annotated.rt) item.rt = *annotated.rt;
          if (annotated.mz) item.mz = *annotated.mz;
          item.spectrum_reference = annotated.data_id;
          item.meta = &annotated;
          if (const auto* best = entry.getBestMatch())
          {
            auto hit = hits_.try_emplace(ID::MatchReference {entry.run->getUuid(), best->getId()});
            if (hit.second)
            {
              const auto peptide = IdentificationDataAdapter::materializePeptide(*entry.run, *best, *entry.run->getPrimaryScore());
              hit.first->second = {peptide, peptide};
            }
            item.top = &hit.first->second.second;
          }
          result.push_back(item);
        }
        return result;
      }

      /// Store the meta values set on the observations and top matches in @p data
      void store(ID& data) const
      {
        for (const auto& [reference, observation] : observations_)
        {
          if (observation.first == observation.second) continue;
          data.getRun(data.findRunByUuid(reference.run_uuid)->getIdentifier()).replaceObservation(reference.query, observation.second);
        }
        for (const auto& [reference, hit] : hits_)
        {
          std::vector<std::string> before, after;
          hit.first.getKeys(before);
          hit.second.getKeys(after);
          auto& run = data.getRun(data.findRunByUuid(reference.run_uuid)->getIdentifier());
          ID::MatchData match = run.getMatch(reference.match);
          bool changed = false;
          for (const auto& key : after)
          {
            if (hit.first.metaValueExists(key) && hit.first.getMetaValue(key) == hit.second.getMetaValue(key)) continue;
            match.setMetaValue(key, hit.second.getMetaValue(key));
            changed = true;
          }
          for (const auto& key : before)
          {
            if (hit.second.metaValueExists(key) || !match.metaValueExists(key)) continue;
            match.removeMetaValue(key);
            changed = true;
          }
          if (changed) run.replaceMatch(reference.match, match);
        }
      }

    private:
      std::map<ID::QueryReference, std::pair<ID::Observation, ID::Observation>> observations_; ///< original, annotated
      std::map<ID::MatchReference, std::pair<PeptideHit, PeptideHit>> hits_;                     ///< original, annotated
    };
  } // namespace

  QCBase::AnnotatedIdentification QCBase::annotated(PeptideIdentification& id)
  {
    AnnotatedIdentification item;
    item.rt = id.getRT();
    item.mz = id.getMZ();
    if (id.metaValueExists("spectrum_reference")) item.spectrum_reference = id.getSpectrumReference();
    item.meta = &id;
    if (!id.getHits().empty()) item.top = &id.getHits()[0];
    return item;
  }

  bool QCBase::hasPepID(const FeatureMap& fmap)
  {
    std::optional<FeatureMap> converted;
    return hasIdentifications(IdentificationDataConverter::withIdentificationData(fmap, converted).getIdentificationData());
  }

  bool QCBase::hasPepID(const ConsensusMap& cmap)
  {
    std::optional<ConsensusMap> converted;
    return hasIdentifications(IdentificationDataConverter::withIdentificationData(cmap, converted).getIdentificationData());
  }

  void QCBase::annotateIdentifications(FeatureMap& map, const std::function<void(Feature*, std::vector<AnnotatedIdentification>&)>& visit, bool include_unassigned)
  {
    IdentificationDataConverter::editAsIdentificationData(map, [&](FeatureMap& native) {
      Annotations annotations;
      const auto& data = native.getIdentificationData();
      for (auto& feature : native)
      {
        auto identifications = annotations.identifications(feature.getLinkedIdentifications(data));
        visit(&feature, identifications);
      }
      if (include_unassigned)
      {
        auto identifications = annotations.identifications(native.getUnassignedIdentifications());
        visit(nullptr, identifications);
      }
      annotations.store(native.getIdentificationData());
    });
  }

  void QCBase::visitIdentifications(const FeatureMap& map, const std::function<void(const Feature*, const std::vector<AnnotatedIdentification>&)>& visit, bool include_unassigned)
  {
    std::optional<FeatureMap> converted;
    const FeatureMap& native = IdentificationDataConverter::withIdentificationData(map, converted);
    Annotations annotations;
    const auto& data = native.getIdentificationData();
    for (const auto& feature : native)
    {
      visit(&feature, annotations.identifications(feature.getLinkedIdentifications(data)));
    }
    if (include_unassigned)
    {
      visit(nullptr, annotations.identifications(native.getUnassignedIdentifications()));
    }
  }

  const SearchParameters* QCBase::searchParameters(const IdentificationData& data)
  {
    return data.getRuns().empty() ? nullptr : &data.getRuns().front().getSettings().search;
  }

  // function tests if a metric has the required input files
  // gives a warning with the name of the metric that can not be performed
  bool QCBase::isRunnable(const Status& s) const
  {
    if (s.isSuperSetOf(this->requirements()))
    {
      return true;
    }
    for (Size i = 0; i < (UInt64)QCBase::Requires::SIZE_OF_REQUIRES; ++i)
    {
      if (this->requirements().isSuperSetOf(QCBase::Requires(i)) && !s.isSuperSetOf(QCBase::Requires(i)))
      {
        OPENMS_LOG_WARN << "Note: Metric '" << this->getName() << "' cannot run because input data '" << QCBase::names_of_requires[i] << "' is missing!\n";
      }
    }
    return false;
  }

  bool QCBase::isLabeledExperiment(const ConsensusMap& cm)
  {
    return DataProcessingUtils::hasIsobaricAnalyzer(cm.getDataProcessing());
  }

} // namespace OpenMS
