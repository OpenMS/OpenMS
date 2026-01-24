// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/USI.h>
#include <OpenMS/CHEMISTRY/ProForma.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/CONCEPT/Constants.h>



using namespace std;

namespace OpenMS
{

  PeptideIdentification::PeptideIdentification() :
    MetaInfoInterface(),
    id_(),
    hits_(),
    score_type_(),
    higher_score_better_(true),
    mz_(std::numeric_limits<double>::quiet_NaN()),
    rt_(std::numeric_limits<double>::quiet_NaN())
  {
  }

  PeptideIdentification::~PeptideIdentification() noexcept = default;

  // Equality operator
  bool PeptideIdentification::operator==(const PeptideIdentification& rhs) const
  {
    return MetaInfoInterface::operator==(rhs)
           && id_ == rhs.id_
           && hits_ == rhs.hits_
           && getSignificanceThreshold() == rhs.getSignificanceThreshold()
           && score_type_ == rhs.score_type_
           && higher_score_better_ == rhs.higher_score_better_
           && getExperimentLabel() == rhs.getExperimentLabel()
           && getBaseName() == rhs.getBaseName()
           && (mz_ == rhs.mz_ || (!this->hasMZ() && !rhs.hasMZ())) // might be NaN, so comparing == will always be false
           && (rt_ == rhs.rt_ || (!this->hasRT() && !rhs.hasRT()));// might be NaN, so comparing == will always be false
  }

  // Inequality operator
  bool PeptideIdentification::operator!=(const PeptideIdentification& rhs) const
  {
    return !(*this == rhs);
  }

  double PeptideIdentification::getRT() const
  {
    return rt_;
  }

  void PeptideIdentification::setRT(double rt)
  {
    rt_ = rt;
  }

  bool PeptideIdentification::hasRT() const
  {
    return !std::isnan(rt_);
  }

  double PeptideIdentification::getMZ() const
  {
    return mz_;
  }

  void PeptideIdentification::setMZ(double mz)
  {
    mz_ = mz;
  }

  bool PeptideIdentification::hasMZ() const
  {
    return !std::isnan(mz_);
  }

  const std::vector<PeptideHit>& PeptideIdentification::getHits() const
  {
    return hits_;
  }

  std::vector<PeptideHit>& PeptideIdentification::getHits()
  {
    return hits_;
  }

  void PeptideIdentification::insertHit(const PeptideHit& hit)
  {
    hits_.push_back(hit);
  }

  void PeptideIdentification::insertHit(PeptideHit&& hit)
  {
    hits_.push_back(std::move(hit));
  }

  void PeptideIdentification::setHits(const std::vector<PeptideHit>& hits)
  {
    hits_ = hits;
  }

  void PeptideIdentification::setHits(std::vector<PeptideHit>&& hits)
  {
    hits_ = std::move(hits);
  }

  double PeptideIdentification::getSignificanceThreshold() const
  {
    return getMetaValue(Constants::UserParam::SIGNIFICANCE_THRESHOLD, 0.0);
  }

  void PeptideIdentification::setSignificanceThreshold(double value)
  {
    if (value != 0.0) 
    {
      setMetaValue(Constants::UserParam::SIGNIFICANCE_THRESHOLD, value);
    }
    else
    {
      // Remove the meta value if the value is 0.0
      removeMetaValue(Constants::UserParam::SIGNIFICANCE_THRESHOLD);
    }
  }

  const String& PeptideIdentification::getScoreType() const
  {
    return score_type_;
  }

  void PeptideIdentification::setScoreType(const String& type)
  {
    score_type_ = type;
  }

  bool PeptideIdentification::isHigherScoreBetter() const
  {
    return higher_score_better_;
  }

  void PeptideIdentification::setHigherScoreBetter(bool value)
  {
    higher_score_better_ = value;
  }

  const String& PeptideIdentification::getIdentifier() const
  {
    return id_;
  }

  void PeptideIdentification::setIdentifier(const String& id)
  {
    id_ = id;
  }

  String PeptideIdentification::getSpectrumReference() const
  {
    return getMetaValue(Constants::UserParam::SPECTRUM_REFERENCE, "");
  }

  void PeptideIdentification::setSpectrumReference(const String& id)
  {
    setMetaValue(Constants::UserParam::SPECTRUM_REFERENCE, id);
  }

  String PeptideIdentification::getBaseName() const
  {
    return getMetaValue(Constants::UserParam::BASE_NAME, "");
  }

  void PeptideIdentification::setBaseName(const String& base_name)
  {
    // do not store empty base_name (default value)
    if (!base_name.empty())
    {
      setMetaValue(Constants::UserParam::BASE_NAME, base_name);
    }
    else
    {
      removeMetaValue(Constants::UserParam::BASE_NAME);
    }
  }

  const String PeptideIdentification::getExperimentLabel() const
  {
    // implement as meta value in order to reduce bloat of PeptideIdentification object
    //  -> this is mostly used for pepxml at the moment which allows each peptide id to belong to a different experiment
    return this->getMetaValue("experiment_label", "");
  }

  void PeptideIdentification::setExperimentLabel(const String& label)
  {
    // do not store empty label (default value)
    if (!label.empty())
    {
      setMetaValue("experiment_label", label);
    }
  }

  std::function<bool(const PeptideHit&, const PeptideHit&)>
  PeptideIdentification::getScoreComparator(bool higher_score_better)
  {
    if (higher_score_better)
    {
      return PeptideHit::ScoreMore();
    }
    else
    {
      return PeptideHit::ScoreLess();
    }
  }

  void PeptideIdentification::sort()
  {
    std::stable_sort(hits_.begin(), hits_.end(), getScoreComparator(higher_score_better_));
  }

  bool PeptideIdentification::empty() const
  {
    return id_.empty()
           && hits_.empty()
           && getSignificanceThreshold() == 0.0
           && score_type_.empty()
           && higher_score_better_ == true
           && getBaseName().empty();
  }

  std::vector<PeptideHit> PeptideIdentification::getReferencingHits(const std::vector<PeptideHit>& hits, const std::set<String>& accession)
  {
    std::vector<PeptideHit> filtered;
    for (const PeptideHit& h_it : hits)
    {
      set<String> hit_accessions = h_it.extractProteinAccessionsSet();
      set<String> intersect;
      set_intersection(hit_accessions.begin(), hit_accessions.end(), accession.begin(), accession.end(), std::inserter(intersect, intersect.begin()));
      if (!intersect.empty())
      {
        filtered.push_back(h_it);
      }
    }
    return filtered;
  }

  std::multimap<String, std::pair<Size, Size>>
  PeptideIdentification::buildUIDsFromAllPepIDs(const ConsensusMap &cmap)
  {
    multimap<String, std::pair<Size, Size>> customID_to_cpepID{};

    IdentifierMSRunMapper mapping(cmap.getProteinIdentifications());
    //Iterates of the vector of PeptideIdentification to build the UID
    //and the pep_index
    auto lambda = [](const PeptideIdentificationList &cpep_ids,
                     const IdentifierMSRunMapper& id_mapping,
                     multimap<String, std::pair<Size, Size>> &customID_to_cpepID,
                     const Size &cf_index) {
        Size pep_index = 0;
        for (const PeptideIdentification &cpep_id : cpep_ids)
        {
          std::pair<Size, Size> index = {cf_index, pep_index};
          auto uid = buildUIDFromPepID(cpep_id, id_mapping);
          customID_to_cpepID.insert(make_pair(uid, index));
          ++pep_index;
        }
    };
    //Build the multimap of all UIDs from the ConsensusMap
    //The multimap maps the UID of the PeptideIdentification to all occurrences in the ConsensusMap
    //An occurrence is described as a pair of an index of the ConsensusFeature and the PeptideIdentification
    for (Size i = 0; i < cmap.size(); ++i)
    {
      lambda(cmap[i].getPeptideIdentifications(), mapping, customID_to_cpepID, i);
    }
    //all unassigned PeptideIdentifications get -1 for the ConsensusFeature index
    lambda(cmap.getUnassignedPeptideIdentifications(), mapping, customID_to_cpepID, -1);

    return customID_to_cpepID;
  }

  String PeptideIdentification::buildUIDFromPepID(const PeptideIdentification &pep_id,
                                                  const IdentifierMSRunMapper& mapping)
  {
    if (!pep_id.metaValueExists("spectrum_reference"))
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Spectrum reference missing at PeptideIdentification.");
    }
    String UID;
    const StringList& ms_run_paths = mapping.getMSRunPaths(pep_id.getIdentifier());
    if (ms_run_paths.size() == 1)
    {
      UID = ms_run_paths[0] + '|' + pep_id.getSpectrumReference();
    }
    else if (pep_id.metaValueExists("map_index"))
    {
      UID = pep_id.getMetaValue("map_index").toString() + '|' +
            pep_id.getSpectrumReference();
    }
    else
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Multiple files in a run, but no map_index in PeptideIdentification found.");
    }
    return UID;
  }

  USI PeptideIdentification::buildUSI(const String& ms_run_name,
                                      const String& dataset_id,
                                      bool include_interpretation) const
  {
    // Get the spectrum reference (native ID)
    String spec_ref = getSpectrumReference();

    // Try to extract scan number from native ID
    auto scan_num = USI::extractScanNumberFromNativeID(spec_ref);

    // Build interpretation string if requested and hits are available
    // Note: The first hit is used as the interpretation. For best results,
    // call sort() before this method to ensure the best-scoring hit is first.
    String interpretation;
    if (include_interpretation && !hits_.empty())
    {
      const PeptideHit& first_hit = hits_.front();
      ProForma::PeptidoformIon ion;
      ion.chains.emplace_back(ProForma::fromAASequence(first_hit.getSequence()));
      if (first_hit.getCharge() != 0)
      {
        ion.charge = first_hit.getCharge();
      }
      interpretation = ProForma::toString(ion);
    }

    if (scan_num.has_value())
    {
      // Use scan number if extraction succeeded
      return USI(dataset_id, ms_run_name, USI::IndexType::SCAN, String(scan_num.value()), interpretation);
    }
    else if (!spec_ref.empty())
    {
      // Fall back to nativeId if scan number extraction failed
      return USI(dataset_id, ms_run_name, USI::IndexType::NATIVEID, spec_ref, interpretation);
    }
    else
    {
      // Return invalid/empty USI if no spectrum reference is available
      return USI();
    }
  }

  USI PeptideIdentification::buildUSI(const IdentifierMSRunMapper& mapping,
                                      const String& dataset_id,
                                      bool include_interpretation) const
  {
    // Get the primary MS run path for this identification using the mapping
    String ms_run_path = mapping.getPrimaryMSRunPath(*this);
    if (ms_run_path.empty())
    {
      return USI(); // Invalid USI if mapping not found
    }

    // Extract basename from the resolved path
    String ms_run_name = USI::extractBasename(ms_run_path);

    return buildUSI(ms_run_name, dataset_id, include_interpretation);
  }

} // namespace OpenMS
